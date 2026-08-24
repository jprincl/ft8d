/*
 * ft8d_airspyhf - Airspy HF+ front end feeding ft8d directly, no external
 * csdr/airspyhf_rx processes involved.
 *
 * v1: single fixed frequency, no band rotation yet (that comes later).
 *
 * Usage: ft8d_airspyhf -f <freq_kHz> [-h <home_locator>] [-sf|-sd|-sn]
 *                       [-l <log_dir>] [-q]
 *   e.g. ft8d_airspyhf -f 14074 -h JO70 -sd -l /home/odroid/ft8logs
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
 *               ft8_<startdate><starttime>_<freq_Hz>.txt in that
 *               directory, e.g. ft8_20260824080614_14075500.txt.
 *               Each 15s interval, decoded or not, ends with a blank
 *               line so the log shows a continuous timeline.
 * -q            suppress screen output. Only useful together with -l --
 *               refused if given without it (nothing would go anywhere).
 * Default (neither given): screen only, as before. -l alone: both.
 * -l with -q: file only.
 */

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
#include <sys/wait.h>
#include <libairspyhf/airspyhf.h>

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
    decode_rec_t recs[MAX_RECS_PER_CYCLE];
    int          nrecs;
} child_t;

static double        g_dial_freq_hz = 14074000.0;  /* nominal, for display only */
static double        g_rf_center_hz = 0.0;          /* dial + FREQ_OFFSET_HZ -- used for BOTH tuning and the .c2 filename, must always match */
static float          g_fir_coef[FIR_NTAPS];
static cplx32_t        g_fir_hist[FIR_NTAPS];
static int              g_fir_pos = 0;
static int              g_decim_ctr = 0;
static cplx32_t          g_outbuf[NMAX];
static int                g_outcount = 0;
static child_t             g_children[MAX_CHILDREN];
static volatile sig_atomic_t g_stop = 0;
static sort_mode_t            g_sort_mode = SORT_FREQ;

static int    g_have_home = 0;
static double g_home_lat = 0.0, g_home_lon = 0.0;
static char   g_home_grid[8] = {0};

static int    g_want_screen = 1; /* -q turns this off */
static FILE  *g_logfile = NULL;  /* set when -l is given */

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

/* Sort this child's buffered decodes per g_sort_mode and write them to
 * whichever destination(s) are active (screen and/or log file -- the
 * log file only ever gets decode lines, never the stderr tracing).
 * Every interval, including a quiet one with zero decodes, ends with a
 * blank-line separator so the log shows a continuous timeline. */
static void flush_child_output(child_t *c)
{
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

            if (g_want_screen) printf("%s\n", outline);
            if (g_logfile) fprintf(g_logfile, "%s\n", outline);
        }
    }
    if (g_want_screen) { printf("\n"); fflush(stdout); }
    if (g_logfile) { fprintf(g_logfile, "\n"); fflush(g_logfile); }
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

static void launch_decode(const char *path)
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
        g_children[slot].used = 1;
    } else {
        perror("fork");
        close(pfd[0]);
        close(pfd[1]);
        unlink(path);
    }
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

    launch_decode(path);
    g_outcount = 0;
}

static int g_armed = 0;          /* have we hit a real 15s UTC boundary yet? */
static int g_last_checked_sec = -1;

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
    double freq_khz = -1.0;
    const char *home_arg = NULL;
    const char *log_dir = NULL;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-f") == 0 && i + 1 < argc)
            freq_khz = atof(argv[++i]);
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
    }
    if (freq_khz <= 0.0) {
        fprintf(stderr, "Usage: %s -f <freq_kHz> [-h <home_locator>] "
                "[-sf|-sd|-sn] [-l <log_dir>] [-q]\n", argv[0]);
        return 1;
    }
    if (!g_want_screen && !log_dir) {
        fprintf(stderr, "-q with no -l means no output would go anywhere\n");
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
    g_dial_freq_hz = freq_khz * 1000.0;
    g_rf_center_hz = g_dial_freq_hz + FREQ_OFFSET_HZ;

    if (log_dir) {
        time_t start = time(NULL);
        struct tm tmv;
        gmtime_r(&start, &tmv);
        size_t dl = strlen(log_dir);
        const char *sep = (dl > 0 && log_dir[dl - 1] == '/') ? "" : "/";
        char logpath[512];
        snprintf(logpath, sizeof(logpath),
                 "%s%sft8_%04d%02d%02d%02d%02d%02d_%.0f.txt",
                 log_dir, sep,
                 tmv.tm_year + 1900, tmv.tm_mon + 1, tmv.tm_mday,
                 tmv.tm_hour, tmv.tm_min, tmv.tm_sec, g_dial_freq_hz);
        g_logfile = fopen(logpath, "a");
        if (!g_logfile) {
            fprintf(stderr, "could not open log file %s: %s\n",
                    logpath, strerror(errno));
            return 1;
        }
        fprintf(stderr, "logging decodes to %s\n", logpath);
    }

    design_lowpass(g_fir_coef, FIR_NTAPS, (double)SAMPLE_RATE_HZ, FIR_CUTOFF_HZ);
    memset(g_children, 0, sizeof(g_children));
    signal(SIGINT, handle_sigint);

    airspyhf_device_t *dev = NULL;
    if (airspyhf_open(&dev) != AIRSPYHF_SUCCESS) {
        fprintf(stderr, "airspyhf_open failed (device plugged in? plugdev group?)\n");
        return 1;
    }

    /* Confirm 192 kS/s is actually offered before committing to it --
     * the list varies by firmware, per the earlier discussion. */
    uint32_t rate_count = 0;
    airspyhf_get_samplerates(dev, &rate_count, 0);
    uint32_t *rates = malloc(rate_count * sizeof(uint32_t));
    airspyhf_get_samplerates(dev, rates, rate_count);
    int have_rate = 0;
    for (uint32_t i = 0; i < rate_count; i++)
        if (rates[i] == SAMPLE_RATE_HZ) have_rate = 1;
    free(rates);
    if (!have_rate) {
        fprintf(stderr, "device does not offer %u sps -- check firmware's rate list\n",
                SAMPLE_RATE_HZ);
        airspyhf_close(dev);
        return 1;
    }

    airspyhf_set_samplerate(dev, SAMPLE_RATE_HZ);
    /* Tune to g_rf_center_hz (dial + FREQ_OFFSET_HZ), the SAME value used
     * in the .c2 filename above -- hardware tuning and the embedded
     * dialfreq must always match, or the frequency column comes out
     * wrong even though decoding itself still works fine. */
    airspyhf_set_freq(dev, (uint32_t)g_rf_center_hz);
    airspyhf_set_hf_agc(dev, 1);   /* matches -g on */
    airspyhf_set_hf_lna(dev, 1);   /* matches -m on */

    if (airspyhf_start(dev, rx_callback, NULL) != AIRSPYHF_SUCCESS) {
        fprintf(stderr, "airspyhf_start failed\n");
        airspyhf_close(dev);
        return 1;
    }

    const char *sort_name = g_sort_mode == SORT_DIST ? "distance" :
                             g_sort_mode == SORT_SNR  ? "SNR" : "frequency";
    fprintf(stderr, "receiving dial %.3f kHz (RF center %.3f kHz, "
            "covering %.0f - %.0f Hz)%s%s, sort by %s, Ctrl+C to stop\n",
            freq_khz, g_rf_center_hz / 1000.0,
            g_rf_center_hz - 1600.0, g_rf_center_hz + 1600.0,
            g_have_home ? ", home " : "",
            g_have_home ? g_home_grid : "",
            sort_name);
    while (!g_stop && airspyhf_is_streaming(dev))
        sleep(1);

    airspyhf_stop(dev);
    airspyhf_close(dev);
    reap_children(1);
    if (g_logfile) fclose(g_logfile);
    return 0;
}
