/*
 * ft8d_airspyhf - Airspy HF+ front end feeding ft8d directly, no external
 * csdr/airspyhf_rx processes involved.
 *
 * Usage: ft8d_airspyhf -f <freq_kHz>[,<freq_kHz>...] [-h <home_locator>]
 *                       [-sf|-sd|-sn] [-l <log_dir>] [-q] [-d <minutes>]
 *   e.g. ft8d_airspyhf -f 14074,18100,21074 -h JO70 -sd -l /home/odroid/ft8logs -d 30
 *
 * -f takes one frequency (no rotation) or a comma-separated list. With a
 * list, bands rotate one minute each, retuning the live device handle --
 * no reopening. Each band gets its own log/diff file if -l is given.
 *
 * -h is optional. When given, a genuine locator found in a decoded
 * message gets a distance-in-km column appended (e.g. "JN53  1560km").
 * RR73/RRR/73 and signal reports are never treated as locators, even
 * though some are syntactically indistinguishable from a real grid.
 *
 * -sf/-sd/-sn pick how each 15s cycle's decodes are sorted before being
 * printed (all at once, once ft8d for that cycle exits):
 *   -sf  frequency ascending (default)
 *   -sd  distance descending (farthest first); rows without a distance
 *        (no -h given, or no locator in the message) sort last
 *   -sn  SNR descending (strongest first)
 *
 * -l <log_dir>  also write decodes (no tracing info) to a file named
 *               <freq_kHz>_<startdate>_<starttime>.ft8 in that
 *               directory, e.g. 18100_20260824_130923.ft8, one per band,
 *               kept open until the program exits. Each 15s interval,
 *               decoded or not, ends with a blank line so the log shows
 *               a continuous timeline.
 * -q            suppress screen output. Only useful together with -l --
 *               refused if given without it (nothing would go anywhere).
 * -d <minutes>  differential mode, window 5-180 minutes, reset on every
 *               restart. A callsign already flagged as new on its band
 *               within the window is suppressed from the screen and the
 *               per-band "<...>_d.ft8" file (only created if -l is also
 *               given) -- but always still goes to the main archive
 *               file, which stays a complete, unfiltered log regardless.
 *
 * Default (neither -l nor -q given): screen only, everything. -l alone:
 * both screen and file, everything. -l with -q: file only. Adding -d
 * filters screen and adds the "_d" file; the main archive is unaffected.
 *
 * -k <api_key>  persisted to ~/.ft8d_airspyhf.conf (mode 600) so you don't
 *               need to retype it -- but never turns uploading on by
 *               itself.
 * -c <cloudlog_url>  the actual per-run toggle: upload happens only when
 *               -c is given THIS run (never saved/remembered), using
 *               whatever key is available (this run's -k, or the saved
 *               one). Requires -d -- Cloudlog upload always reuses the
 *               same per-band dedup state as the diff filter. Uses curl
 *               as a fire-and-forget subprocess (no libcurl dependency,
 *               no blocking the audio thread). Fixed two apparent bugs
 *               from the original 2-year-old script: RST_SENT no longer
 *               gets a stray "-8" appended, and COMMENT no longer gets a
 *               trailing ":". FREQ/FREQ_RX are sent in ADIF-standard MHz
 *               (the old script sent raw Hz).
 */

#define _POSIX_C_SOURCE 200809L /* gmtime_r, strdup -- explicit, don't rely on GNU-dialect defaults */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <libairspyhf/airspyhf.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define SAMPLE_RATE_HZ   192000u   /* native Airspy HF+ rate we ask for */
#define DECIM            48        /* 192000 / 48 = 4000, exactly what ft8d wants */
#define NMAX             60000     /* 15 s at 4000 Hz -- matches ft8_params.f90 */

/* Single-stage FIR anti-alias/decimation filter. Generous tap count since
 * compute headroom is not a concern here (a few million MACs/sec, trivial
 * for the C5). Good enough to get real decodes flowing; a multistage
 * decimator (e.g. x4, x4, x3) would be a cleaner design to revisit later
 * if filter quality ever turns out to matter in practice. */
#define FIR_NTAPS        1025
#define FIR_CUTOFF_HZ    1900.0

#define FT8D_PATH        "./ft8d"    /* adjust if your layout differs */
#define MAX_CHILDREN     8
#define CHILD_TIMEOUT_S  45          /* decode should finish well inside 60 s */

/* ft8d.f90's search window is +/-1600 Hz, symmetric around whatever
 * "dialfreq" ends up embedded in the .c2 filename. Real FT8 activity on a
 * given band conventionally spans roughly dial+200 to dial+2900 Hz (the
 * usual USB-audio convention), which is wider than +/-1600 Hz centered on
 * the plain dial. Shifting the actual RF center (and the dialfreq we
 * embed -- both together, consistently) by this much centers the window
 * on that real activity band instead, with margin on both sides. 1500 Hz
 * matches the value from Jan's original working shell script. */
#define FREQ_OFFSET_HZ   1500.0

#define MAX_RECS_PER_CYCLE 128
#define MAX_BANDS 16
#define MAX_TRACKED_CALLS 512
#define MAX_UPLOADS 8

/* Fixed station identity for ADIF records, matching Jan's original
 * script. MY_GRIDSQUARE deliberately reuses g_home_grid (from -h)
 * instead of a separate hardcoded constant, so the two always agree. */
#define CLOUDLOG_OPERATOR  "SWL/OK/THEP111"
#define CLOUDLOG_CITY      "Hracholusky"
#define CLOUDLOG_COUNTRY   "CZECH REPUBLIC"
#define CLOUDLOG_DXCC      "503"
#define CLOUDLOG_CQ_ZONE   "15"
#define CLOUDLOG_ITU_ZONE  "28"

typedef struct { float re, im; } cplx32_t;

typedef enum { SORT_FREQ, SORT_DIST, SORT_SNR } sort_mode_t;

typedef struct {
    char   prefix[40]; /* dtime..freq, ~32 chars */
    char   msg[64];
    char   call[16];   /* from ft8d's msgcall -- parsed for future ADIF use, not displayed */
    char   grid[5];
    int    has_grid;
    double dist_km;
    int    has_dist;
    double freq_hz;
    double snr;
} decode_rec_t;

typedef struct {
    pid_t   pid;
    char    path[160];
    time_t  started;
    int     used;
    int     readfd;       /* read end of ft8d's redirected stdout */
    char    linebuf[256]; /* partial (not yet newline-terminated) output */
    int     linelen;
    int     band_idx;     /* which band this decode belongs to, tagged at launch */
    char    date_str[9];  /* YYYYMMDD at capture time -- ft8d's own dtime is HHMMSS only */
    decode_rec_t recs[MAX_RECS_PER_CYCLE];
    int          nrecs;
} child_t;

typedef struct { char call[16]; time_t last_seen; } call_seen_t;
typedef struct { call_seen_t seen[MAX_TRACKED_CALLS]; int count; } band_dedup_t;

typedef struct { pid_t pid; char path[160]; time_t started; int used; } upload_t;

static airspyhf_device_t *g_dev = NULL;

static double g_band_freq_khz[MAX_BANDS];
static int    g_nbands = 0;
static int    g_band_idx = 0;
static time_t g_band_start = 0;

static double        g_dial_freq_hz = 14074000.0;  /* nominal, for display only */
static double        g_rf_center_hz = 0.0;          /* dial + FREQ_OFFSET_HZ -- used for BOTH tuning and the .c2 filename, must always match */
static float          g_fir_coef[FIR_NTAPS];
static cplx32_t        g_fir_hist[FIR_NTAPS];
static int              g_fir_pos = 0;
static int              g_decim_ctr = 0;
static cplx32_t          g_outbuf[NMAX];
static int                g_outcount = 0;
static child_t             g_children[MAX_CHILDREN];
static upload_t             g_uploads[MAX_UPLOADS];
static volatile sig_atomic_t g_stop = 0;
static sort_mode_t            g_sort_mode = SORT_FREQ;
static int                     g_armed = 0; /* have we hit a real 15s UTC boundary yet? */
static int                     g_last_checked_sec = -1;

static int    g_have_home = 0;
static double g_home_lat = 0.0, g_home_lon = 0.0;
static char   g_home_grid[8] = {0};

static int    g_want_screen = 1; /* -q turns this off */
static FILE  *g_logfiles[MAX_BANDS];  /* main per-band archive, set when -l is given */
static FILE  *g_difffiles[MAX_BANDS]; /* per-band "_d" novelty-only file, needs -l and -d together */

static int          g_diff_mode = 0;
static int           g_diff_window_min = 0;
static band_dedup_t   g_dedup[MAX_BANDS];

static int   g_have_cloudlog = 0;
static char  g_cloudlog_key[128] = {0};
static char  g_cloudlog_url[256] = {0};

static void handle_sigint(int sig) { (void)sig; g_stop = 1; }

/* Hamming-windowed sinc lowpass, DC gain normalised to 1.0. Because it's
 * symmetric it doesn't matter which end of the circular history buffer
 * we start multiplying from in feed_sample() below. */
static void design_lowpass(float *coef, int ntaps, double fs, double fc)
{
    int M = ntaps - 1;
    double sum = 0.0;
    for (int n = 0; n <= M; n++) {
        double x = n - M / 2.0;
        double s = (x == 0.0) ? 2.0 * fc / fs
                               : sin(2.0 * M_PI * fc * x / fs) / (M_PI * x);
        double w = 0.54 - 0.46 * cos(2.0 * M_PI * n / M);
        coef[n] = (float)(s * w);
        sum += coef[n];
    }
    for (int n = 0; n < ntaps; n++)
        coef[n] = (float)(coef[n] / sum);
}

/* True if g[0..3] is a syntactically valid 4-char Maidenhead locator
 * (2 letters A-R, 2 digits 0-9). Does NOT by itself rule out RR73 --
 * that exclusion happens by exact-string-match in extract_locator(),
 * since RR73 is syntactically indistinguishable from a real locator. */
static int is_valid_grid4(const char *g)
{
    return g[0] >= 'A' && g[0] <= 'R' &&
           g[1] >= 'A' && g[1] <= 'R' &&
           g[2] >= '0' && g[2] <= '9' &&
           g[3] >= '0' && g[3] <= '9';
}

static void grid_to_latlon(const char *grid, double *lat, double *lon)
{
    char c1 = (char)toupper((unsigned char)grid[0]);
    char c2 = (char)toupper((unsigned char)grid[1]);
    *lon = (c1 - 'A') * 20.0 - 180.0 + (grid[2] - '0') * 2.0 + 1.0;
    *lat = (c2 - 'A') * 10.0 - 90.0 + (grid[3] - '0') * 1.0 + 0.5;
}

static double haversine_km(double lat1, double lon1, double lat2, double lon2)
{
    const double R = 6371.0; /* mean Earth radius, km */
    double dlat = (lat2 - lat1) * M_PI / 180.0;
    double dlon = (lon2 - lon1) * M_PI / 180.0;
    double a = sin(dlat / 2) * sin(dlat / 2) +
               cos(lat1 * M_PI / 180.0) * cos(lat2 * M_PI / 180.0) *
               sin(dlon / 2) * sin(dlon / 2);
    double c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));
    return R * c;
}

/* True if `call` on this band hasn't been flagged as new in the last
 * g_diff_window_min minutes (or has never been seen at all). Per-band,
 * matching the confirmed dedup scope. Mirrors Jan's old diff.sh
 * semantics: the window is measured from when a call was LAST flagged
 * as new, not refreshed on every intervening sighting -- so once it's
 * quiet for the window, it can be flagged again. */
static int dedup_is_new(int band, const char *call, time_t now)
{
    if (call[0] == '\0') return 1; /* nothing to key on, don't gate it */
    band_dedup_t *d = &g_dedup[band];
    for (int i = 0; i < d->count; i++) {
        if (strcmp(d->seen[i].call, call) == 0) {
            if (now - d->seen[i].last_seen < g_diff_window_min * 60)
                return 0; /* still within the window, not new */
            d->seen[i].last_seen = now; /* re-arm */
            return 1;
        }
    }
    int slot = d->count;
    if (slot >= MAX_TRACKED_CALLS) {
        int oldest = 0;
        for (int i = 1; i < d->count; i++)
            if (d->seen[i].last_seen < d->seen[oldest].last_seen) oldest = i;
        slot = oldest;
    } else {
        d->count++;
    }
    strncpy(d->seen[slot].call, call, sizeof(d->seen[slot].call) - 1);
    d->seen[slot].call[sizeof(d->seen[slot].call) - 1] = '\0';
    d->seen[slot].last_seen = now;
    return 1;
}

static void cloudlog_config_path(char *buf, size_t n)
{
    const char *home = getenv("HOME");
    snprintf(buf, n, "%s/.ft8d_airspyhf.conf", home ? home : ".");
}

/* Only fills g_cloudlog_key if it's still empty, so a value already
 * given on the command line always wins over the saved file. The URL
 * is deliberately never persisted here -- -c is a per-run toggle for
 * whether upload happens at all, not a remembered setting, so a saved
 * key alone can never silently turn uploading back on. */
static void load_cloudlog_key(void)
{
    char path[512];
    cloudlog_config_path(path, sizeof(path));
    FILE *f = fopen(path, "r");
    if (!f) return;
    char line[400];
    while (fgets(line, sizeof(line), f)) {
        char *nl = strchr(line, '\n');
        if (nl) *nl = '\0';
        if (strncmp(line, "api_key=", 8) == 0 && !g_cloudlog_key[0])
            strncpy(g_cloudlog_key, line + 8, sizeof(g_cloudlog_key) - 1);
    }
    fclose(f);
}

static void save_cloudlog_key(void)
{
    char path[512];
    cloudlog_config_path(path, sizeof(path));
    FILE *f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "warning: could not save %s: %s\n", path, strerror(errno));
        return;
    }
    fprintf(f, "api_key=%s\n", g_cloudlog_key);
    fclose(f);
    chmod(path, S_IRUSR | S_IWUSR); /* 600 -- this file holds a secret */
}

/* Minimal JSON string escaping for embedding the ADIF blob as a JSON
 * string value. Caller frees the result. */
static char *json_escape(const char *s)
{
    size_t len = strlen(s);
    char *out = malloc(len * 2 + 1); /* worst case: every char escaped */
    if (!out) return NULL;
    size_t j = 0;
    for (size_t i = 0; i < len; i++) {
        unsigned char c = (unsigned char)s[i];
        if (c == '"' || c == '\\') { out[j++] = '\\'; out[j++] = (char)c; }
        else if (c == '\n') { out[j++] = '\\'; out[j++] = 'n'; }
        else if (c == '\r') { /* drop */ }
        else if (c < 0x20) { /* drop other control chars */ }
        else out[j++] = (char)c;
    }
    out[j] = '\0';
    return out;
}

static void reap_uploads(int force_all)
{
    time_t now = time(NULL);
    for (int i = 0; i < MAX_UPLOADS; i++) {
        if (!g_uploads[i].used) continue;
        int status;
        pid_t r = waitpid(g_uploads[i].pid, &status, WNOHANG);
        if (r == g_uploads[i].pid) {
            unlink(g_uploads[i].path);
            g_uploads[i].used = 0;
            continue;
        }
        if (force_all || (now - g_uploads[i].started) > 15) {
            kill(g_uploads[i].pid, SIGKILL);
            waitpid(g_uploads[i].pid, &status, 0);
            unlink(g_uploads[i].path);
            g_uploads[i].used = 0;
        }
    }
}

/* Writes the JSON payload to /dev/shm and fork/execs curl to POST it,
 * fire-and-forget -- never blocks the audio thread on a slow/unreachable
 * Cloudlog server. adif_blob may contain multiple <eor>-terminated
 * records (one per newly-flagged decode this cycle). */
static void post_to_cloudlog(const char *adif_blob)
{
    if (!g_have_cloudlog) return;
    reap_uploads(0);

    int slot = -1;
    for (int i = 0; i < MAX_UPLOADS; i++)
        if (!g_uploads[i].used) { slot = i; break; }
    if (slot < 0) {
        reap_uploads(1);
        slot = 0;
    }

    char *escaped = json_escape(adif_blob);
    if (!escaped) return;

    char path[160];
    snprintf(path, sizeof(path), "/dev/shm/cloudlog_%ld_%d.json",
             (long)time(NULL), slot);
    FILE *f = fopen(path, "w");
    if (!f) {
        perror("fopen cloudlog payload");
        free(escaped);
        return;
    }
    fprintf(f, "{\"key\":\"%s\",\"station_profile_id\":\"1\",\"type\":\"adif\",\"string\":\"%s\"}",
            g_cloudlog_key, escaped);
    fclose(f);
    free(escaped);

    char data_arg[180];
    snprintf(data_arg, sizeof(data_arg), "@%s", path);

    pid_t pid = fork();
    if (pid == 0) {
        execlp("curl", "curl", "-s", "-o", "/dev/null",
               "-w", "cloudlog upload: HTTP %{http_code}\n",
               "-X", "POST", "-H", "Content-Type: application/json",
               "-d", data_arg, g_cloudlog_url, (char *)NULL);
        perror("execlp curl");
        _exit(127);
    } else if (pid > 0) {
        strncpy(g_uploads[slot].path, path, sizeof(g_uploads[slot].path) - 1);
        g_uploads[slot].path[sizeof(g_uploads[slot].path) - 1] = '\0';
        g_uploads[slot].pid = pid;
        g_uploads[slot].started = time(NULL);
        g_uploads[slot].used = 1;
    } else {
        perror("fork");
        unlink(path);
    }
}

/* Builds one ADIF record for a diff-filtered decode. Appends to out
 * (caller-managed buffer/offset) rather than returning a new string, so
 * a whole cycle's worth of records can be accumulated into one blob
 * before a single upload. */
static void append_adif_record(const decode_rec_t *r, const char *date_str,
                                char *out, size_t outsz, int *off)
{
    char snr_buf[8];
    snprintf(snr_buf, sizeof(snr_buf), "%d", (int)r->snr);

    char freq_mhz[16];
    snprintf(freq_mhz, sizeof(freq_mhz), "%.6f", r->freq_hz / 1000000.0);

    char time_str[7];
    memcpy(time_str, r->prefix, 6); /* dtime is HHMMSS, always the first 6 chars */
    time_str[6] = '\0';

    int n = *off;
    n += snprintf(out + n, outsz - (size_t)n, "<CALL:%zu>%s", strlen(r->call), r->call);
    n += snprintf(out + n, outsz - (size_t)n, "<MODE:3>FT8");
    if (r->has_grid)
        n += snprintf(out + n, outsz - (size_t)n, "<GRIDSQUARE:%zu>%s", strlen(r->grid), r->grid);
    n += snprintf(out + n, outsz - (size_t)n, "<OPERATOR:%zu>%s",
                  strlen(CLOUDLOG_OPERATOR), CLOUDLOG_OPERATOR);
    n += snprintf(out + n, outsz - (size_t)n, "<FREQ:%zu>%s", strlen(freq_mhz), freq_mhz);
    n += snprintf(out + n, outsz - (size_t)n, "<FREQ_RX:%zu>%s", strlen(freq_mhz), freq_mhz);
    n += snprintf(out + n, outsz - (size_t)n, "<RST_SENT:%zu>%s", strlen(snr_buf), snr_buf);
    n += snprintf(out + n, outsz - (size_t)n, "<RST_RCVD:3>---");
    n += snprintf(out + n, outsz - (size_t)n, "<QSO_DATE:8>%s", date_str);
    n += snprintf(out + n, outsz - (size_t)n, "<TIME_ON:6>%s", time_str);
    n += snprintf(out + n, outsz - (size_t)n, "<QSO_DATE_OFF:8>%s", date_str);
    n += snprintf(out + n, outsz - (size_t)n, "<TIME_OFF:6>%s", time_str);
    n += snprintf(out + n, outsz - (size_t)n, "<STATION_CALLSIGN:%zu>%s",
                  strlen(CLOUDLOG_OPERATOR), CLOUDLOG_OPERATOR);
    n += snprintf(out + n, outsz - (size_t)n, "<MY_CITY:%zu>%s",
                  strlen(CLOUDLOG_CITY), CLOUDLOG_CITY);
    n += snprintf(out + n, outsz - (size_t)n, "<MY_COUNTRY:%zu>%s",
                  strlen(CLOUDLOG_COUNTRY), CLOUDLOG_COUNTRY);
    n += snprintf(out + n, outsz - (size_t)n, "<MY_DXCC:%zu>%s",
                  strlen(CLOUDLOG_DXCC), CLOUDLOG_DXCC);
    if (g_have_home)
        n += snprintf(out + n, outsz - (size_t)n, "<MY_GRIDSQUARE:%zu>%s",
                      strlen(g_home_grid), g_home_grid);
    n += snprintf(out + n, outsz - (size_t)n, "<MY_CQ_ZONE:%zu>%s",
                  strlen(CLOUDLOG_CQ_ZONE), CLOUDLOG_CQ_ZONE);
    n += snprintf(out + n, outsz - (size_t)n, "<MY_ITU_ZONE:%zu>%s",
                  strlen(CLOUDLOG_ITU_ZONE), CLOUDLOG_ITU_ZONE);
    n += snprintf(out + n, outsz - (size_t)n, "<COMMENT:%zu>%s", strlen(r->msg), r->msg);
    n += snprintf(out + n, outsz - (size_t)n, "<eor>\n");
    *off = n;
}

/* Looks at the last whitespace-separated token of the (already-trimmed)
 * message text. Returns 1 and fills grid_out (5 bytes incl. NUL) if it's
 * a genuine locator -- explicitly excludes RR73/RRR/73 and signal
 * reports, which can be syntactically indistinguishable from (or
 * adjacent to) a real grid. */
static int extract_locator(const char *msg, char *grid_out)
{
    const char *end = msg + strlen(msg);
    while (end > msg && isspace((unsigned char)*(end - 1))) end--;
    const char *start = end;
    while (start > msg && !isspace((unsigned char)*(start - 1))) start--;
    size_t len = (size_t)(end - start);

    if (len != 4) return 0;

    char tok[5];
    memcpy(tok, start, 4);
    tok[4] = '\0';

    if (strcmp(tok, "RR73") == 0) return 0; /* roger/sign-off, not a grid */
    if (!is_valid_grid4(tok)) return 0;      /* also rejects things like R-09 */

    memcpy(grid_out, tok, 5);
    return 1;
}

/* Reads up to `len` chars starting at `start` from s (bounds-checked
 * against s's actual length) and parses them as a number. Returns 0.0
 * for anything out of range instead of reading past the string. */
static double parse_field(const char *s, int start, int len)
{
    int slen = (int)strlen(s);
    if (start >= slen) return 0.0;
    int avail = slen - start;
    int n = len < avail ? len : avail;
    if (n <= 0) return 0.0;
    char buf[16];
    int m = n < (int)sizeof(buf) - 1 ? n : (int)sizeof(buf) - 1;
    memcpy(buf, s + start, m);
    buf[m] = '\0';
    return atof(buf);
}

static int cmp_freq(const void *a, const void *b)
{
    const decode_rec_t *ra = a, *rb = b;
    if (ra->freq_hz < rb->freq_hz) return -1;
    if (ra->freq_hz > rb->freq_hz) return 1;
    return 0;
}

static int cmp_snr(const void *a, const void *b) /* highest SNR first */
{
    const decode_rec_t *ra = a, *rb = b;
    if (ra->snr > rb->snr) return -1;
    if (ra->snr < rb->snr) return 1;
    return 0;
}

static int cmp_dist(const void *a, const void *b) /* farthest first, no-distance rows last */
{
    const decode_rec_t *ra = a, *rb = b;
    if (ra->has_dist != rb->has_dist)
        return ra->has_dist ? -1 : 1;
    if (!ra->has_dist) return 0;
    if (ra->dist_km > rb->dist_km) return -1;
    if (ra->dist_km < rb->dist_km) return 1;
    return 0;
}

/* Parse a completed line from ft8d's stdout into a record and buffer it
 * on the child (rather than printing immediately), so a full 15s cycle's
 * worth of decodes can be sorted together once ft8d exits. Splits the
 * fixed-width Fortran prefix precisely rather than relying on its a20
 * message padding, which breaks down for messages over 20 chars.
 * msgcall (the trailing field) is peeled from the END of the line
 * instead, since its Fortran type (character*13) makes that width
 * genuinely guaranteed, unlike msg37's a20 which can overflow. */
static void process_decode_line(child_t *c, char *line)
{
    size_t len = strlen(line);
    while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r'))
        line[--len] = '\0';
    if (len == 0) return;

    const size_t PREFIX_LEN = 32;
    char prefix[40];
    size_t body_start, body_len;
    if (len > PREFIX_LEN + 1) {
        size_t n = PREFIX_LEN < sizeof(prefix) - 1 ? PREFIX_LEN : sizeof(prefix) - 1;
        memcpy(prefix, line, n);
        prefix[n] = '\0';
        body_start = PREFIX_LEN + 1; /* skip the 1x separator */
        body_len = len - body_start;
    } else {
        prefix[0] = '\0';
        body_start = 0;
        body_len = len;
    }

    /* peel msgcall off the end of the body, if there's room for it */
    char call[16] = {0};
    const size_t CALL_LEN = 13;
    if (body_len > CALL_LEN + 1) {
        size_t call_off = body_start + body_len - CALL_LEN;
        memcpy(call, line + call_off, CALL_LEN);
        call[CALL_LEN] = '\0';
        size_t clen = strlen(call);
        while (clen > 0 && isspace((unsigned char)call[clen - 1])) call[--clen] = '\0';
        body_len -= CALL_LEN + 1; /* also drop the 1x separator before it */
    }

    char msg[64];
    size_t n = body_len < sizeof(msg) - 1 ? body_len : sizeof(msg) - 1;
    memcpy(msg, line + body_start, n);
    msg[n] = '\0';
    size_t mlen = strlen(msg);
    while (mlen > 0 && isspace((unsigned char)msg[mlen - 1]))
        msg[--mlen] = '\0';
    if (mlen == 0) return;

    if (c->nrecs >= MAX_RECS_PER_CYCLE) return; /* silently drop past the cap */
    decode_rec_t *r = &c->recs[c->nrecs++];
    strncpy(r->prefix, prefix, sizeof(r->prefix) - 1);
    r->prefix[sizeof(r->prefix) - 1] = '\0';
    strncpy(r->msg, msg, sizeof(r->msg) - 1);
    r->msg[sizeof(r->msg) - 1] = '\0';
    strncpy(r->call, call, sizeof(r->call) - 1);
    r->call[sizeof(r->call) - 1] = '\0';
    r->freq_hz = parse_field(prefix, 23, 9);
    r->snr = parse_field(prefix, 13, 4);

    if (extract_locator(msg, r->grid)) {
        r->has_grid = 1;
        if (g_have_home) {
            double lat, lon;
            grid_to_latlon(r->grid, &lat, &lon);
            r->dist_km = haversine_km(g_home_lat, g_home_lon, lat, lon);
            r->has_dist = 1;
        } else {
            r->has_dist = 0;
        }
    } else {
        r->has_grid = 0;
        r->has_dist = 0;
    }
}

/* Sort this child's buffered decodes per g_sort_mode and write them out.
 * The per-band archive file (if -l given) always gets everything. Screen
 * and the per-band "_d" file (if -d given) only get lines whose callsign
 * passes the dedup filter -- i.e. wasn't already flagged new within the
 * window on this band. Every interval, including a quiet one, ends with
 * a blank-line separator on every active destination. */
static void flush_child_output(child_t *c)
{
    int band = c->band_idx;
    FILE *archive = g_logfiles[band];
    FILE *diff_f = g_difffiles[band];
    time_t now = time(NULL);
    char adif_blob[MAX_RECS_PER_CYCLE * 220];
    int  adif_off = 0;

    if (c->nrecs > 0) {
        int (*cmp)(const void *, const void *) = cmp_freq;
        if (g_sort_mode == SORT_DIST) cmp = cmp_dist;
        else if (g_sort_mode == SORT_SNR) cmp = cmp_snr;
        qsort(c->recs, c->nrecs, sizeof(c->recs[0]), cmp);

        for (int i = 0; i < c->nrecs; i++) {
            decode_rec_t *r = &c->recs[i];

            /* r->prefix is "dtime sync snr dt freq" (a6,1x,f6.1,i4,f6.2,i9).
             * Drop the sync field, positions 7 through 12, for display --
             * it's an internal sync8 metric, not useful once a line has
             * already passed CRC and made it this far. */
            char display_prefix[32];
            snprintf(display_prefix, sizeof(display_prefix), "%.7s%s",
                     r->prefix, r->prefix + 13);

            char outline[160];
            int off = snprintf(outline, sizeof(outline), "%s %-20s %-13s",
                                display_prefix, r->msg, r->call);
            if (r->has_grid && r->has_dist)
                snprintf(outline + off, sizeof(outline) - off, " %-6s %6.0fkm",
                         r->grid, r->dist_km);
            else if (r->has_grid)
                snprintf(outline + off, sizeof(outline) - off, " %-6s", r->grid);

            if (archive) fprintf(archive, "%s\n", outline);

            int is_new = !g_diff_mode || dedup_is_new(band, r->call, now);
            if (is_new) {
                if (g_want_screen) printf("%s\n", outline);
                if (diff_f) fprintf(diff_f, "%s\n", outline);
                if (g_have_cloudlog)
                    append_adif_record(r, c->date_str, adif_blob, sizeof(adif_blob), &adif_off);
            }
        }
    }
    if (g_want_screen) { printf("\n"); fflush(stdout); }
    if (archive) { fprintf(archive, "\n"); fflush(archive); }
    if (diff_f) { fprintf(diff_f, "\n"); fflush(diff_f); }
    if (g_have_cloudlog && adif_off > 0)
        post_to_cloudlog(adif_blob);
    c->nrecs = 0;
}

/* Non-blocking drain of whatever ft8d has written so far, processing
 * every complete (newline-terminated) line and keeping any trailing
 * partial line buffered for next time. */
static void drain_child_output(child_t *c)
{
    char buf[512];
    for (;;) {
        ssize_t n = read(c->readfd, buf, sizeof(buf));
        if (n <= 0) break; /* EAGAIN (nothing available) or EOF */
        for (ssize_t i = 0; i < n; i++) {
            if (c->linelen < (int)sizeof(c->linebuf) - 1)
                c->linebuf[c->linelen++] = buf[i];
            if (buf[i] == '\n') {
                c->linebuf[c->linelen] = '\0';
                process_decode_line(c, c->linebuf);
                c->linelen = 0;
            }
        }
    }
}

/* Reap finished/overdue ft8d children and unlink their .c2 file every
 * time -- whether or not anything decoded, per the "always clean up"
 * rule, not just on success. */
static void reap_children(int force_all)
{
    time_t now = time(NULL);
    for (int i = 0; i < MAX_CHILDREN; i++) {
        if (!g_children[i].used) continue;
        drain_child_output(&g_children[i]);
        int status;
        pid_t r = waitpid(g_children[i].pid, &status, WNOHANG);
        if (r == g_children[i].pid) {
            drain_child_output(&g_children[i]); /* catch final buffered output */
            flush_child_output(&g_children[i]);
            close(g_children[i].readfd);
            unlink(g_children[i].path);
            g_children[i].used = 0;
            continue;
        }
        if (force_all || (now - g_children[i].started) > CHILD_TIMEOUT_S) {
            fprintf(stderr, "ft8d pid %d over %ds, killing\n",
                    (int)g_children[i].pid, CHILD_TIMEOUT_S);
            kill(g_children[i].pid, SIGKILL);
            waitpid(g_children[i].pid, &status, 0);
            drain_child_output(&g_children[i]);
            flush_child_output(&g_children[i]);
            close(g_children[i].readfd);
            unlink(g_children[i].path);
            g_children[i].used = 0;
        }
    }
}

static void launch_decode(const char *path, const char *date_str)
{
    reap_children(0);

    int slot = -1;
    for (int i = 0; i < MAX_CHILDREN; i++)
        if (!g_children[i].used) { slot = i; break; }
    if (slot < 0) {
        fprintf(stderr, "no free child slot, forcing cleanup\n");
        reap_children(1);
        slot = 0;
    }

    int pfd[2];
    if (pipe(pfd) != 0) {
        perror("pipe");
        unlink(path);
        return;
    }

    pid_t pid = fork();
    if (pid == 0) {
        close(pfd[0]);
        dup2(pfd[1], STDOUT_FILENO);
        close(pfd[1]);
        execlp(FT8D_PATH, "ft8d", path, (char *)NULL);
        perror("execlp ft8d");
        _exit(127);
    } else if (pid > 0) {
        close(pfd[1]);
        fcntl(pfd[0], F_SETFL, O_NONBLOCK);
        strncpy(g_children[slot].path, path, sizeof(g_children[slot].path) - 1);
        g_children[slot].path[sizeof(g_children[slot].path) - 1] = '\0';
        g_children[slot].pid = pid;
        g_children[slot].started = time(NULL);
        g_children[slot].readfd = pfd[0];
        g_children[slot].linelen = 0;
        g_children[slot].nrecs = 0;
        g_children[slot].band_idx = g_band_idx;
        strncpy(g_children[slot].date_str, date_str, sizeof(g_children[slot].date_str) - 1);
        g_children[slot].date_str[sizeof(g_children[slot].date_str) - 1] = '\0';
        g_children[slot].used = 1;
    } else {
        perror("fork");
        close(pfd[0]);
        close(pfd[1]);
        unlink(path);
    }
}

/* If rotating (more than one band) and we've held the current band for
 * at least 60s, retune to the next one on the live device handle -- no
 * reopening. Discards any in-flight buffer/FIR history and forces a
 * fresh wait for the next 15s UTC boundary, since samples right after a
 * retune are transitional and shouldn't be mixed with the new band. */
static void maybe_rotate_band(void)
{
    if (g_nbands <= 1) return;
    time_t now = time(NULL);
    if (now - g_band_start < 60) return;

    g_band_idx = (g_band_idx + 1) % g_nbands;
    g_dial_freq_hz = g_band_freq_khz[g_band_idx] * 1000.0;
    g_rf_center_hz = g_dial_freq_hz + FREQ_OFFSET_HZ;
    airspyhf_set_freq(g_dev, (uint32_t)g_rf_center_hz);

    memset(g_fir_hist, 0, sizeof(g_fir_hist));
    g_fir_pos = 0;
    g_decim_ctr = 0;
    g_outcount = 0;
    g_armed = 0;
    fprintf(stderr, "rotating to %.3f kHz\n", g_band_freq_khz[g_band_idx]);
}

/* Filename must end in exactly HHMMSS_FREQ8.c2 -- ft8d.f90 parses those
 * last 15 characters before ".c2" directly (dtime = last-15..-10,
 * frequency = last-8..end). Anything before that is free-form. */
static void write_c2_and_launch(void)
{
    time_t now = time(NULL);
    struct tm tmv;
    gmtime_r(&now, &tmv);

    char path[160];
    snprintf(path, sizeof(path), "/dev/shm/ft8_%02d%02d%02d_%08.0f.c2",
             tmv.tm_hour, tmv.tm_min, tmv.tm_sec, g_rf_center_hz);

    FILE *f = fopen(path, "wb");
    if (!f) {
        perror("fopen c2");
        g_outcount = 0;
        return;
    }
    size_t n = fwrite(g_outbuf, sizeof(cplx32_t), NMAX, f);
    fclose(f);
    if (n != NMAX) {
        fprintf(stderr, "short write on %s (%zu/%d samples), skipping\n",
                path, n, NMAX);
        unlink(path);
        g_outcount = 0;
        return;
    }
    fprintf(stderr, "captured 15s -> %s, launching ft8d\n", path);

    char date_str[16]; /* worst-case int width for %04d, though the real
                           value is always exactly 8 digits (YYYYMMDD) */
    snprintf(date_str, sizeof(date_str), "%04d%02d%02d",
             tmv.tm_year + 1900, tmv.tm_mon + 1, tmv.tm_mday);

    launch_decode(path, date_str);
    g_outcount = 0;
    maybe_rotate_band();
}

static void append_output_sample(float re, float im)
{
    if (!g_armed) {
        time_t now = time(NULL);
        struct tm tmv;
        gmtime_r(&now, &tmv);
        if (tmv.tm_sec != g_last_checked_sec) {
            g_last_checked_sec = tmv.tm_sec;
            if (tmv.tm_sec % 15 == 0) {
                g_armed = 1;
                g_outcount = 0;
                g_band_start = now; /* dwell budget starts NOW, not at retune time --
                                       the alignment wait itself must not eat into it */
                fprintf(stderr, "aligned to 15s UTC boundary, starting capture\n");
            }
        }
        if (!g_armed)
            return; /* still waiting for alignment -- discard this sample */
    }

    if (g_outcount < NMAX) {
        g_outbuf[g_outcount].re = re;
        g_outbuf[g_outcount].im = im;
        g_outcount++;
    }
    if (g_outcount >= NMAX)
        write_c2_and_launch();
}

static void feed_sample(airspyhf_complex_float_t s)
{
    g_fir_hist[g_fir_pos].re = s.re;
    g_fir_hist[g_fir_pos].im = s.im;
    g_fir_pos = (g_fir_pos + 1) % FIR_NTAPS;

    if (++g_decim_ctr < DECIM)
        return;
    g_decim_ctr = 0;

    float acc_re = 0.0f, acc_im = 0.0f;
    int idx = g_fir_pos; /* oldest sample in the history */
    for (int i = 0; i < FIR_NTAPS; i++) {
        acc_re += g_fir_coef[i] * g_fir_hist[idx].re;
        acc_im += g_fir_coef[i] * g_fir_hist[idx].im;
        idx = (idx + 1) % FIR_NTAPS;
    }
    append_output_sample(acc_re, acc_im);
}

static int g_first_seen = 0;

static int rx_callback(airspyhf_transfer_t *transfer)
{
    if (!g_first_seen) {
        g_first_seen = 1;
        fprintf(stderr, "first callback: %d samples, dropped=%llu\n",
                transfer->sample_count,
                (unsigned long long)transfer->dropped_samples);
    }

    if (transfer->dropped_samples > 0)
        fprintf(stderr, "warning: %llu dropped samples\n",
                (unsigned long long)transfer->dropped_samples);

    for (int i = 0; i < transfer->sample_count; i++)
        feed_sample(transfer->samples[i]);

    reap_children(0);
    return 0;
}

int main(int argc, char **argv)
{
    const char *freq_arg = NULL;
    const char *home_arg = NULL;
    const char *log_dir = NULL;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-f") == 0 && i + 1 < argc)
            freq_arg = argv[++i];
        else if (strcmp(argv[i], "-h") == 0 && i + 1 < argc)
            home_arg = argv[++i];
        else if (strcmp(argv[i], "-sf") == 0)
            g_sort_mode = SORT_FREQ;
        else if (strcmp(argv[i], "-sd") == 0)
            g_sort_mode = SORT_DIST;
        else if (strcmp(argv[i], "-sn") == 0)
            g_sort_mode = SORT_SNR;
        else if (strcmp(argv[i], "-l") == 0) {
            /* -l with no path (end of args, or another flag right after)
             * means "current directory", not "silently ignore -l". */
            if (i + 1 < argc && argv[i + 1][0] != '-')
                log_dir = argv[++i];
            else
                log_dir = ".";
        }
        else if (strcmp(argv[i], "-q") == 0)
            g_want_screen = 0;
        else if (strcmp(argv[i], "-d") == 0 && i + 1 < argc) {
            g_diff_window_min = atoi(argv[++i]);
            g_diff_mode = 1;
        }
        else if (strcmp(argv[i], "-k") == 0 && i + 1 < argc)
            strncpy(g_cloudlog_key, argv[++i], sizeof(g_cloudlog_key) - 1);
        else if (strcmp(argv[i], "-c") == 0 && i + 1 < argc)
            strncpy(g_cloudlog_url, argv[++i], sizeof(g_cloudlog_url) - 1);
    }
    if (!freq_arg) {
        fprintf(stderr, "Usage: %s -f <freq_kHz>[,<freq_kHz>...] [-h <home_locator>] "
                "[-sf|-sd|-sn] [-l <log_dir>] [-q] [-d <minutes>] [-k <api_key>] "
                "[-c <cloudlog_url>]\n", argv[0]);
        return 1;
    }
    {
        char *fcopy = strdup(freq_arg);
        char *tok = strtok(fcopy, ",");
        while (tok && g_nbands < MAX_BANDS) {
            g_band_freq_khz[g_nbands++] = atof(tok);
            tok = strtok(NULL, ",");
        }
        free(fcopy);
    }
    if (g_nbands == 0) {
        fprintf(stderr, "no valid frequency in -f '%s'\n", freq_arg);
        return 1;
    }
    if (!g_want_screen && !log_dir) {
        fprintf(stderr, "-q with no -l means no output would go anywhere\n");
        return 1;
    }
    if (g_diff_mode && (g_diff_window_min < 5 || g_diff_window_min > 180)) {
        fprintf(stderr, "-d window must be between 5 and 180 minutes\n");
        return 1;
    }
    if (g_cloudlog_key[0]) {
        save_cloudlog_key(); /* -k given this run -- remember it for next time */
    } else {
        load_cloudlog_key(); /* -k not given -- try the saved one */
    }
    g_have_cloudlog = g_cloudlog_url[0] != '\0'; /* -c given THIS run turns upload on */
    if (g_have_cloudlog && !g_cloudlog_key[0]) {
        fprintf(stderr, "-c given but no API key available -- pass -k (once, to save it, "
                "or every run)\n");
        return 1;
    }
    if (g_have_cloudlog && !g_diff_mode) {
        fprintf(stderr, "-c (Cloudlog upload) requires -d -- only diff-filtered "
                "decodes get uploaded\n");
        return 1;
    }
    if (home_arg) {
        char g[5];
        size_t hl = strlen(home_arg);
        if (hl != 4) {
            fprintf(stderr, "invalid -h locator '%s' (need exactly 4 characters)\n",
                    home_arg);
            return 1;
        }
        g[0] = (char)toupper((unsigned char)home_arg[0]);
        g[1] = (char)toupper((unsigned char)home_arg[1]);
        g[2] = home_arg[2];
        g[3] = home_arg[3];
        g[4] = '\0';
        if (!is_valid_grid4(g)) {
            fprintf(stderr, "invalid -h locator '%s' (expected e.g. JO70)\n", home_arg);
            return 1;
        }
        grid_to_latlon(g, &g_home_lat, &g_home_lon);
        strncpy(g_home_grid, g, sizeof(g_home_grid) - 1);
        g_have_home = 1;
    }
    g_band_idx = 0;
    g_dial_freq_hz = g_band_freq_khz[0] * 1000.0;
    g_rf_center_hz = g_dial_freq_hz + FREQ_OFFSET_HZ;

    if (log_dir) {
        time_t start = time(NULL);
        struct tm tmv;
        gmtime_r(&start, &tmv);
        size_t dl = strlen(log_dir);
        const char *sep = (dl > 0 && log_dir[dl - 1] == '/') ? "" : "/";
        for (int b = 0; b < g_nbands; b++) {
            char logpath[512];
            snprintf(logpath, sizeof(logpath),
                     "%s%s%.0f_%04d%02d%02d_%02d%02d%02d.ft8",
                     log_dir, sep, g_band_freq_khz[b],
                     tmv.tm_year + 1900, tmv.tm_mon + 1, tmv.tm_mday,
                     tmv.tm_hour, tmv.tm_min, tmv.tm_sec);
            g_logfiles[b] = fopen(logpath, "a");
            if (!g_logfiles[b]) {
                fprintf(stderr, "could not open log file %s: %s\n",
                        logpath, strerror(errno));
                return 1;
            }
            fprintf(stderr, "logging decodes to %s\n", logpath);

            if (g_diff_mode) {
                char diffpath[520];
                snprintf(diffpath, sizeof(diffpath), "%.*s_d.ft8",
                         (int)(strlen(logpath) - 4), logpath);
                g_difffiles[b] = fopen(diffpath, "a");
                if (!g_difffiles[b]) {
                    fprintf(stderr, "could not open diff file %s: %s\n",
                            diffpath, strerror(errno));
                    return 1;
                }
                fprintf(stderr, "logging new-only decodes to %s\n", diffpath);
            }
        }
    }

    design_lowpass(g_fir_coef, FIR_NTAPS, (double)SAMPLE_RATE_HZ, FIR_CUTOFF_HZ);
    memset(g_children, 0, sizeof(g_children));
    signal(SIGINT, handle_sigint);

    if (airspyhf_open(&g_dev) != AIRSPYHF_SUCCESS) {
        fprintf(stderr, "airspyhf_open failed (device plugged in? plugdev group?)\n");
        return 1;
    }

    /* Confirm 192 kS/s is actually offered before committing to it --
     * the list varies by firmware, per the earlier discussion. */
    uint32_t rate_count = 0;
    airspyhf_get_samplerates(g_dev, &rate_count, 0);
    uint32_t *rates = malloc(rate_count * sizeof(uint32_t));
    airspyhf_get_samplerates(g_dev, rates, rate_count);
    int have_rate = 0;
    for (uint32_t i = 0; i < rate_count; i++)
        if (rates[i] == SAMPLE_RATE_HZ) have_rate = 1;
    free(rates);
    if (!have_rate) {
        fprintf(stderr, "device does not offer %u sps -- check firmware's rate list\n",
                SAMPLE_RATE_HZ);
        airspyhf_close(g_dev);
        return 1;
    }

    airspyhf_set_samplerate(g_dev, SAMPLE_RATE_HZ);
    /* Tune to g_rf_center_hz (dial + FREQ_OFFSET_HZ), the SAME value used
     * in the .c2 filename above -- hardware tuning and the embedded
     * dialfreq must always match, or the frequency column comes out
     * wrong even though decoding itself still works fine. */
    airspyhf_set_freq(g_dev, (uint32_t)g_rf_center_hz);
    airspyhf_set_hf_agc(g_dev, 1);   /* matches -g on */
    airspyhf_set_hf_lna(g_dev, 1);   /* matches -m on */

    if (airspyhf_start(g_dev, rx_callback, NULL) != AIRSPYHF_SUCCESS) {
        fprintf(stderr, "airspyhf_start failed\n");
        airspyhf_close(g_dev);
        return 1;
    }
    g_band_start = time(NULL);

    const char *sort_name = g_sort_mode == SORT_DIST ? "distance" :
                             g_sort_mode == SORT_SNR  ? "SNR" : "frequency";
    if (g_nbands == 1) {
        fprintf(stderr, "receiving dial %.3f kHz (RF center %.3f kHz, "
                "covering %.0f - %.0f Hz)%s%s, sort by %s",
                g_band_freq_khz[0], g_rf_center_hz / 1000.0,
                g_rf_center_hz - 1600.0, g_rf_center_hz + 1600.0,
                g_have_home ? ", home " : "",
                g_have_home ? g_home_grid : "",
                sort_name);
        if (g_diff_mode) fprintf(stderr, ", diff window %d min", g_diff_window_min);
        if (g_have_cloudlog) fprintf(stderr, ", uploading to Cloudlog");
        fprintf(stderr, ", Ctrl+C to stop\n");
    } else {
        fprintf(stderr, "rotating %d bands (1 min each):", g_nbands);
        for (int b = 0; b < g_nbands; b++)
            fprintf(stderr, " %.3fkHz", g_band_freq_khz[b]);
        fprintf(stderr, "%s%s, sort by %s", g_have_home ? ", home " : "",
                g_have_home ? g_home_grid : "", sort_name);
        if (g_diff_mode) fprintf(stderr, ", diff window %d min", g_diff_window_min);
        if (g_have_cloudlog) fprintf(stderr, ", uploading to Cloudlog");
        fprintf(stderr, ", Ctrl+C to stop\n");
    }
    while (!g_stop && airspyhf_is_streaming(g_dev))
        sleep(1);

    airspyhf_stop(g_dev);
    airspyhf_close(g_dev);
    reap_children(1);
    reap_uploads(1);
    for (int b = 0; b < g_nbands; b++) {
        if (g_logfiles[b]) fclose(g_logfiles[b]);
        if (g_difffiles[b]) fclose(g_difffiles[b]);
    }
    return 0;
}
