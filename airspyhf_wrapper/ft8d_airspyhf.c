/*
 * ft8d_airspyhf - Airspy HF+ front end feeding ft8d directly, no external
 * csdr/airspyhf_rx processes involved.
 *
 * v1: single fixed frequency, no band rotation yet (that comes later).
 *
 * Usage: ft8d_airspyhf -f <freq_kHz>
 *   e.g. ft8d_airspyhf -f 14074
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>
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

#define FT8D_PATH        "../ft8d"   /* adjust if your layout differs */
#define MAX_CHILDREN     8
#define CHILD_TIMEOUT_S  45          /* decode should finish well inside 60 s */

typedef struct { float re, im; } cplx32_t;

typedef struct {
    pid_t   pid;
    char    path[160];
    time_t  started;
    int     used;
} child_t;

static double        g_dial_freq_hz = 14074000.0;
static float          g_fir_coef[FIR_NTAPS];
static cplx32_t        g_fir_hist[FIR_NTAPS];
static int              g_fir_pos = 0;
static int              g_decim_ctr = 0;
static cplx32_t          g_outbuf[NMAX];
static int                g_outcount = 0;
static child_t             g_children[MAX_CHILDREN];
static volatile sig_atomic_t g_stop = 0;

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

/* Reap finished/overdue ft8d children and unlink their .c2 file every
 * time -- whether or not anything decoded, per the "always clean up"
 * rule, not just on success. */
static void reap_children(int force_all)
{
    time_t now = time(NULL);
    for (int i = 0; i < MAX_CHILDREN; i++) {
        if (!g_children[i].used) continue;
        int status;
        pid_t r = waitpid(g_children[i].pid, &status, WNOHANG);
        if (r == g_children[i].pid) {
            unlink(g_children[i].path);
            g_children[i].used = 0;
            continue;
        }
        if (force_all || (now - g_children[i].started) > CHILD_TIMEOUT_S) {
            fprintf(stderr, "ft8d pid %d over %ds, killing\n",
                    (int)g_children[i].pid, CHILD_TIMEOUT_S);
            kill(g_children[i].pid, SIGKILL);
            waitpid(g_children[i].pid, &status, 0);
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

    pid_t pid = fork();
    if (pid == 0) {
        execlp(FT8D_PATH, "ft8d", path, (char *)NULL);
        perror("execlp ft8d");
        _exit(127);
    } else if (pid > 0) {
        strncpy(g_children[slot].path, path, sizeof(g_children[slot].path) - 1);
        g_children[slot].path[sizeof(g_children[slot].path) - 1] = '\0';
        g_children[slot].pid = pid;
        g_children[slot].started = time(NULL);
        g_children[slot].used = 1;
    } else {
        perror("fork");
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
             tmv.tm_hour, tmv.tm_min, tmv.tm_sec, g_dial_freq_hz);

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

static void append_output_sample(float re, float im)
{
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

static int rx_callback(airspyhf_transfer_t *transfer)
{
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
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-f") == 0 && i + 1 < argc)
            freq_khz = atof(argv[++i]);
    }
    if (freq_khz <= 0.0) {
        fprintf(stderr, "Usage: %s -f <freq_kHz>\n", argv[0]);
        return 1;
    }
    g_dial_freq_hz = freq_khz * 1000.0;

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
    /* ft8d.f90 expects "audio zero" 2000 Hz below dialfreq (same convention
     * as a real SSB receiver -- see nfa+2000/nfb+2000/f1-2000+dialfreq in
     * ft8d.f90). Our complex receiver instead centers its +/-2000 Hz
     * passband on whatever it's tuned to, so the actual RF tuning needs to
     * sit 2000 Hz above the nominal dial frequency for the FT8 sub-band
     * (dial+400..dial+3600 Hz) to land inside that passband. dialfreq
     * itself (filename, decode math) stays the plain nominal value. */
    airspyhf_set_freq(dev, (uint32_t)(g_dial_freq_hz + 2000.0));
    airspyhf_set_hf_agc(dev, 1);   /* matches -g on */
    airspyhf_set_hf_lna(dev, 1);   /* matches -m on */

    if (airspyhf_start(dev, rx_callback, NULL) != AIRSPYHF_SUCCESS) {
        fprintf(stderr, "airspyhf_start failed\n");
        airspyhf_close(dev);
        return 1;
    }

    fprintf(stderr, "receiving on %.3f kHz, Ctrl+C to stop\n", freq_khz);
    while (!g_stop && airspyhf_is_streaming(dev))
        sleep(1);

    airspyhf_stop(dev);
    airspyhf_close(dev);
    reap_children(1);
    return 0;
}
