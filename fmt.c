/*
  Frequency estimator for ARRL Frequency Measurement Test (FMT)
  18 Apr 2026 - KA9Q (with help from ChatGPT)
  Requires FFTW3:
   gcc -O2 -Wall -o fmt fmt.c -lfftw3 -lm

   Sample use:
   fmt -v -w 3.0 -i 1.5 -S '2026-04-17 02:49:06.587' -d 60.078 -l 90. -h 110 7065k2026-04-17T02:38:34.8Z.wav

   args: -v: turn on verbose output
         -w [seconds] FFT window duration
	 -i [seconds] Amount to slide each FFT
	 -S [date]    UTC of keydown time (requires unixstarttime attribute)
	 or
	 -s [seconds] keydown time offset from file start
	 -d [seconds] duration of keydown
	 -l [Hz]      lower limit of search (may be negative)
	 -h [Hz]      upper limit of search (may be negative)
	 -t [0.0 to .99] Trim fraction (default 0.10)
	 -O [oversample] FFT zero-padding ratio, default 4 (pad time domain data to 4x length)
	 -o [Hz]      threshold for discarding frequency outliers from median (1 hz default)

    The file argument should be produced by pcmrecord (part of ka9q-radio)
    as a 16-bit PCM .wav file with 2 channels (IQ).
    fmt uses two attributes created by pcmrecord:
     "unixstarttime" (UTC start of file)
     "frequency" (radio frequency corresponding to 0 Hz in the IQ data)

    sample rate is extracted from the .wav header

    The -v option enables dumping of the individual FFT results: the time within the analysis
    interval, the frequency estimate, and the energy relative to the average.

    The algorithm first runs a series of FFTs across the analysis interval with specified width
    and spacing (windows can and probably should overlap, eg, by 50%). The estimated frequency
    is extracted using quadratic interpolation from the two bins surrounding the peak bin
    within the allowed range.

 */

#define _GNU_SOURCE 1
#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <stdbool.h>
#include <getopt.h>
#include <fcntl.h>
#include <errno.h>
#include <assert.h>
#include <locale.h>
#include <sys/xattr.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
  double time_sec;    // midpoint time of window */
  double freq_hz;     // estimated frequency */
  double power;       // peak power */
  int sequence;       // index in sort by frequency
  bool valid;         // true if accepted so far
  bool weak;          // True if rejected due to low energy
  bool trimmed;       // true if in trimmed tails of distribution
  bool outlier;       // true if discarded as outlier
} estimate_t;

typedef struct {
  double median_hz;
  double trimmed_weighted_mean_hz;
  double drift_ref_time_sec;
  double drift_ref_freq_hz;
  double drift_slope_hz_per_sec;
  int    n_total;
  int    n_good;
  double residual_rms_hz;
} summary_t;

/* ---------- Utilities ---------- */

static inline float normf(float complex x){
  return crealf(x) * crealf(x) + cimagf(x) * cimagf(x);
}
#if 0
static inline double norm(double complex x){
  return creal(x) * creal(x) + cimag(x) * cimag(x);
}
#endif

int load_wav(float **signal, int* num_samples, int *num_channels, int* sample_rate, const char* path,int fd);

int Verbose = 0;

static int cmp_double(const void *a, const void *b){
  double da = *(const double *)a;
  double db = *(const double *)b;
  if (da < db) return -1;
  if (da > db) return 1;
  return 0;
}

static int cmp_estimate_freq(const void *a, const void *b){
  const estimate_t *ea = *(const estimate_t **)a;
  const estimate_t *eb = *(const estimate_t **)b;
  if (ea->freq_hz < eb->freq_hz) return -1;
  if (ea->freq_hz > eb->freq_hz) return 1;
  return 0;
}

static double median_double(double *x, int n){
  assert(n > 0);
  qsort(x, n, sizeof(double), cmp_double);
  if (n & 1) {
    return x[n / 2];
  } else {
    return 0.5 * (x[n / 2 - 1] + x[n / 2]);
  }
}

static void make_hann_window(float *w, int n){
  for (int i = 0; i < n; i++) {
    w[i] = 0.5 - 0.5 * cos(2.0 * M_PI * i / (n - 1));
  }
}

static double quadratic_peak_offset_logpower(double ym1, double y0, double yp1){
  assert(!isnan(ym1) && !isnan(y0) && !isnan(yp1));
  /* 3-point parabolic interpolation in log-power */
  double denom = ym1 - 2.0 * y0 + yp1;
  if (fabs(denom) < 1e-30)
    return 0.0;
  return 0.5 * (ym1 - yp1) / denom;
}

// Estimate frequency within a single FFT window, using quadratic interpolation
static int estimate_window_frequency(const float complex *x,  // input window
				     int n,                    // window length
				     double fs,                // sample rate
				     int nfft,                 // FFT size, >= n
				     double search_lo_hz,
				     double search_hi_hz,
				     fftwf_plan plan,
				     float complex *fft_in,
				     float complex *fft_out,
				     const float *window,
				     double *freq_hz_out,
				     double *power_out){
  // Remove mean and apply window
  double complex mean = 0.0;
  for (int i = 0; i < n; i++){
    assert(!isnan(crealf(x[i])) && !isnan(cimag(x[i])));
    mean += x[i];
  }
  mean /= (double)n;

  for (int i = 0; i < n; i++) {
    fft_in[i] = (x[i] - mean) * window[i];
  }
  for (int i = n; i < nfft; i++) {
    fft_in[i] = 0;
  }

  fftwf_execute(plan);

  // Search only requested frequency span
  int best_k = -1;
  double best_p = -1.0;

  for (int k = 0; k < nfft; k++) {
    // fftfreq-like mapping
    double f = (k <= nfft / 2) ? (fs * k / nfft) : (fs * (k - nfft) / nfft);

    if (f < search_lo_hz || f > search_hi_hz)
      continue;

    double p = normf(fft_out[k]);
    if (p > best_p) {
      best_p = p;
      best_k = k;
    }
  }
  if (best_k < 0)
    return -1;

  // Need neighbors for interpolation
  if (best_k == 0 || best_k == nfft - 1) {
    double f = (best_k <= nfft / 2) ? (fs * best_k / nfft)
      : (fs * (best_k - nfft) / nfft);
    assert(!isnan(f) && !isnan(best_p));
    *freq_hz_out = f;
    *power_out = best_p;
    return 0;
  }

  double pm1 = normf(fft_out[best_k - 1]);
  double p0  = normf(fft_out[best_k]);
  double pp1 = normf(fft_out[best_k + 1]);

  assert(pm1 > 0);
  double ym1 = log(pm1 + 1e-300);
  assert(p0 > 0);
  double y0  = log(p0  + 1e-300);
  assert(pp1 > 0);
  double yp1 = log(pp1 + 1e-300);

  double delta = quadratic_peak_offset_logpower(ym1, y0, yp1);
  assert(!isnan(delta));
  double fbin  = (best_k <= nfft / 2) ? (fs * best_k / nfft)
    : (fs * (best_k - nfft) / nfft);
  double bin_hz = fs / nfft;

  assert(!isnan(fs) && nfft != 0);
  *freq_hz_out = fbin + delta * bin_hz;
  *power_out   = best_p;
  assert(!isnan(*freq_hz_out) && !isnan(*power_out));
  return 0;
}

/* Process list of estimates:
   1. Find median frequency
   2. Discard frequency outliers
   3. Discard weak energies (so many dB below mean)
   4. Sort by frequency, trim tails (currently 10%)
   5. Compute mean of survivors weighted by energies
   6. Compute RMS of weighted survivor distribution
*/
static int summarize_estimates(estimate_t *est,
			       int n_est,
			       double min_rel_power_db,
			       double outlier_hz,
			       double trim_fraction,
			       summary_t *out){
  assert(n_est > 0);
  if (n_est <= 0)
    return -1;

  out->n_total = n_est;
  out->n_good  = 0;

  // Find average power
  double pavg = 0;
  for (int i = 1; i < n_est; i++) {
    pavg += est[i].power;
  }
  pavg = 10 * log10(pavg / n_est); // convert to dB

  /* Power gate */
  int n_power = 0;
  for (int i = 0; i < n_est; i++) {
    double rel_db = 10.0 * log10(est[i].power + 1e-300) - pavg;
    //	fprintf(stdout,"%d %.1lf dB\n",i,rel_db);
    est[i].valid = (rel_db >= min_rel_power_db);
    est[i].weak = !est[i].valid;
    n_power += est[i].valid;
  }
  assert(n_power > 0);
  if (n_power <= 0)
    return -1;

  // Median of valid freqs
  double *tmp = malloc((size_t)n_power * sizeof(*tmp));
  assert(tmp != NULL);
  if (tmp == NULL)
    return -1;

  int j = 0;
  for (int i = 0; i < n_est; i++) {
    if (est[i].valid){
      assert(j < n_power);
      tmp[j++] = est[i].freq_hz;
    }
  }
  double med = median_double(tmp, n_power);
  out->median_hz = med;
  free(tmp);

  // Outlier gate around median
  int n_good = 0;
  for (int i = 0; i < n_est; i++) {
    // Count outliers for statistics even if already invalid due to weakness
    if(fabs(est[i].freq_hz - med) > outlier_hz){
      est[i].outlier = true;
      est[i].valid = false;
    } else
      n_good += est[i].valid;
  }
  assert(n_good > 0);
  if (n_good <= 0)
    return -1;

  out->n_good = n_good;

  // Collect valid estimates and sort by frequency for trimming
  estimate_t **good = malloc((size_t)n_good * sizeof(**good));
  assert(good != NULL);
  if (good == NULL)
    return -1;

  j = 0;
  for (int i = 0; i < n_est; i++) {
    if (est[i].valid){
      assert(j < n_good);
      good[j++] = &est[i];
    }
  }
  qsort(good, n_good, sizeof(estimate_t **), cmp_estimate_freq);

  int ktrim = (int)floor(trim_fraction * n_good);
  int lo = ktrim;
  int hi = n_good - ktrim;
  if (hi <= lo) {
    lo = 0;
   hi = n_good;
  }

  // Weighted mean after trimming
  for(int i=0; i < lo; i++){
    good[i]->sequence = i;
    good[i]->trimmed = true;
  }
  double sw = 0.0, sf = 0.0;
  for (int i = lo; i < hi; i++) {
    good[i]->sequence = i;
    double w = good[i]->power;
    sw += w;
    sf += w * good[i]->freq_hz;
  }
  for(int i=hi; i < n_good; i++){
    good[i]->sequence = i;
    good[i]->trimmed = true;
  }
  out->trimmed_weighted_mean_hz = (sw > 0.0) ? (sf / sw) : good[(lo + hi) / 2]->freq_hz;

  // Weighted linear fit: f(t) = a + b*(t - t0)
  double t0 = 0.0;
  double swt = 0.0;
  for (int i = lo; i < hi; i++) {
    double w = good[i]->power;
    sw += 0.0; /* keep compiler quiet if reused mentally; harmless */
    swt += w;
    t0 += w * good[i]->time_sec;
  }
  if (swt > 0.0)
    t0 /= swt;
  else
    t0 = good[(lo + hi) / 2]->time_sec;

  double S0 = 0.0, S1 = 0.0, S2 = 0.0;
  double T0 = 0.0, T1 = 0.0;

  for (int i = lo; i < hi; i++) {
    double w = good[i]->power;
    double tt = good[i]->time_sec - t0;
    double ff = good[i]->freq_hz;

    S0 += w;
    S1 += w * tt;
    S2 += w * tt * tt;
    T0 += w * ff;
    T1 += w * tt * ff;
  }

  double det = S0 * S2 - S1 * S1;
  double a, b;
  if (fabs(det) < 1e-30) {
    a = out->trimmed_weighted_mean_hz;
    b = 0.0;
  } else {
    a = (T0 * S2 - T1 * S1) / det;
    b = (S0 * T1 - S1 * T0) / det;
  }
  assert(!isnan(a) && !isnan(b) && !isnan(t0));
  out->drift_ref_time_sec   = t0;
  out->drift_ref_freq_hz    = a;
  out->drift_slope_hz_per_sec = b;
  double sse = 0.0;
  double sw_rms = 0.0;

  for (int i = lo; i < hi; i++) {
    double w  = good[i]->power;
    double tt = good[i]->time_sec - t0;
    double f_fit = a + b * tt;
    double err = good[i]->freq_hz - f_fit;

    sse += w * err * err;
    sw_rms += w;
  }
  double rms = 0.0;
  if (sw_rms > 0.0)
    rms = sqrt(sse / sw_rms);

  out->residual_rms_hz = rms;
  free(good);
  return 0;
}

// Top level function called by main()
static int estimate_track(const float complex *iq,
			  int64_t nsamp,
			  double fs,
			  double win_sec,
			  double hop_sec,
			  int nfft,
			  double search_lo_hz,
			  double search_hi_hz,
			  double min_rel_power_db,
			  double outlier_hz,
			  double trim_fraction,
			  estimate_t **est_out,
			  int *n_est_out,
			  summary_t *summary_out){
  int win_n = (int)llround(win_sec * fs);
  int hop_n = (int)llround(hop_sec * fs);

  assert (win_n > 8 && hop_n > 0 && nfft >= win_n);

  if (win_n <= 8 || hop_n <= 0 || nfft < win_n)
    return -1;

  int n_est = 0;
  if (nsamp >= win_n)
    n_est = 1 + (int)((nsamp - win_n) / hop_n);

  assert(n_est > 0);
  if (n_est <= 0)
    return -1;

  estimate_t *est = calloc((size_t)n_est, sizeof(*est));
  assert(est != NULL);
  if (est == NULL)
    return -1;

  float *window = malloc((size_t)win_n * sizeof(*window));
  assert(window != NULL);
  if (window == NULL) {
    free(est);
    return -1;
  }
  make_hann_window(window, win_n);

  float complex *fft_in  = fftwf_malloc((size_t)nfft * sizeof(*fft_in));
  float complex *fft_out = fftwf_malloc((size_t)nfft * sizeof(*fft_out));
  assert(fft_in != NULL && fft_out != NULL);
  if (fft_in == NULL || fft_out == NULL) {
    free(window);
    free(est);
    if (fft_in)  fftwf_free(fft_in);
    if (fft_out) fftwf_free(fft_out);
    return -1;
  }

  fftwf_plan plan = fftwf_plan_dft_1d(nfft, fft_in, fft_out, FFTW_FORWARD, FFTW_MEASURE);
  assert(plan != NULL);
  if (plan == NULL) {
    free(window);
    free(est);
    fftwf_free(fft_in);
    fftwf_free(fft_out);
    return -1;
  }

  int idx = 0;
  for (int64_t start = 0; start + win_n <= nsamp; start += hop_n) {
    double f_est = NAN, p_est = NAN;
    if (estimate_window_frequency(iq + start,
				  win_n,
				  fs,
				  nfft,
				  search_lo_hz,
				  search_hi_hz,
				  plan,
				  fft_in,
				  fft_out,
				  window,
				  &f_est,
				  &p_est) == 0) {
      assert(!isnan(f_est) && !isnan(p_est));
      est[idx].time_sec = ((double)start + 0.5 * win_n) / fs;
      est[idx].freq_hz  = f_est;
      est[idx].power    = p_est;
      est[idx].valid    = true;
      est[idx].sequence = -1; // overwritten before sort of valid entries
      idx++;
    }
  }

  fftwf_destroy_plan(plan);
  fftwf_free(fft_in);
  fftwf_free(fft_out);
  free(window);

  if (idx <= 0) {
    free(est);
    return -1;
  }

  if (summarize_estimates(est, idx, min_rel_power_db, outlier_hz, trim_fraction, summary_out) != 0) {
    free(est); // fails here
    return -1;
  }

  *est_out   = est;
  *n_est_out = idx;
  return 0;
}
const char *Command;

void usage(){
  fprintf(stdout,"Usage: %s -l low_freq -h high_freq [-s start_sec] [-d duration_sec] [-o outlier Hz] [-O oversample_n] [-m minpowerdB] [file]\n",Command);
}


char const *Locale = NULL;
int main(int argc,char *argv[]){
  Command = argv[0];
  {
    char const * const cp = getenv("LANG");
    if(cp != NULL)
      Locale = cp;
  }
#if __linux__
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

  double freq_lo = -INFINITY;
  double freq_hi = +INFINITY;
  double win_sec = 1.0;
  double hop_sec = 0;
  double start = -INFINITY;
  double duration_sec = INFINITY;
  double outlier_hz = 1.0; // Hz; estimates greater than this are discarded
  int oversample = 4;
  double min_rel_power_db = -15.0; // dB threshold for accepting a window estimate
  struct timespec start_time = {0};
  struct timespec file_time = {0};
  double trim_fraction = 0.1;

  int c;
  while((c = getopt(argc,argv,"l:h:w:i:s:S:d:o:O:t:m:v?")) != -1){
    switch(c){
    case 't':
      trim_fraction = strtod(optarg,NULL);
      if(trim_fraction < 0.0 || trim_fraction >= 1.0){
	fprintf(stdout,"Trim %lf out of range: should be 0-1\n",trim_fraction);
	exit(1);
      }
      break;
    case 'v':
      Verbose++;
      break;
    case 'l':
      freq_lo = strtod(optarg,NULL);
      break;
    case 'h':
      freq_hi = strtod(optarg,NULL);
      break;
    case 'w':
      win_sec = strtod(optarg,NULL);
      break;
    case 'i':
      hop_sec = strtod(optarg,NULL);
      break;
    case 's':
      start = strtod(optarg,NULL);
      break;
    case 'S':
      {
	struct tm tm = {0};
	char *ptr;
	if((ptr = strptime(optarg, "%Y-%m-%d %H:%M:%S", &tm)) != NULL){
	  start_time.tv_sec = timegm(&tm);  // or mktime()
	  start_time.tv_nsec = 0; // handle later?
	  if(*ptr == '.'){
	    // Fraction specified
	    double fract = strtod(ptr,NULL);
	    start_time.tv_nsec = llrint(fract * 1e9);
	  }

	} else {
	  fprintf(stdout,"Invalid -S time format, use YYYY-mm-dd HH:MM:SS. Must be UTC\n");
	}
      }
      break;
    case 'd':
      duration_sec = strtod(optarg,NULL);
      break;
    case 'o':
      outlier_hz = strtod(optarg,NULL);
      break;
    case 'O':
      oversample = strtol(optarg,NULL,0);
      break;
    case 'm':
      min_rel_power_db = strtod(optarg,NULL);
      break;
    case '?':
    default:;
      usage();
      exit(1);
    }
  }
  setlocale(LC_ALL,Locale); // Set either the hardwired default or the value of $LANG if it exists
  if(optind == argc){
    usage();
    exit(1);
  }
  if(!isfinite(freq_lo) || !isfinite(freq_hi)){
    usage();
    exit(1);
  }
  // ensure freq_hi > freq_lo
  if(freq_lo > freq_hi){
    double tmp = freq_lo;
    freq_lo = freq_hi;
    freq_hi = tmp;
  }
  char const * const file = argv[optind];
  int fd = open(file,O_RDONLY);
  if(fd == -1){
    fprintf(stdout,"can't open %s: %s\n",file,strerror(errno));
    exit(1);
  }
  fprintf(stdout,"%s\n",file);
  double base_frequency = 0;
  {
    // Look for extended file attribute "user.frequency" (linux) or "frequency" (macos)
    char att_buffer[1024] = {0}; // Shouldn't be anywhere near this long

#ifdef __linux__
    ssize_t as = getxattr(file,"user.frequency",att_buffer,sizeof(att_buffer) - 1);
#else
    ssize_t as = getxattr(file,"frequency",att_buffer,sizeof(att_buffer) - 1,0,0);
#endif
    if(as > 0){
      // Extract from attribute
      base_frequency = strtod(att_buffer,NULL);
      fprintf(stdout,"Base frequency %'lf Hz\n",base_frequency);
    }
  }
  {
    // User specified an absolute starting time in UTC, get file start and calculate
    char att_buffer[1024] = {0}; // Shouldn't be anywhere near this long
#ifdef __linux__
    int as = getxattr(file,"user.unixstarttime",att_buffer,sizeof(att_buffer) - 1);
#else
    int as = getxattr(file,"unixstarttime",att_buffer,sizeof(att_buffer) - 1,0,0);
#endif
    if(as > 0){
      char *ptr = NULL;
      file_time.tv_sec = strtol(att_buffer,&ptr,0); // integer part
      if(ptr != NULL && *ptr == '.')
	file_time.tv_nsec = strtol(ptr+1,&ptr,0);

      struct tm tm;
      struct tm *tmp = gmtime_r(&file_time.tv_sec,&tm);
      if(tmp != NULL)
	fprintf(stdout,"File start time %ld.%ld = %04d-%02d-%02d %02d:%02d:%02d.%09ld UTC\n",
		file_time.tv_sec,file_time.tv_nsec,
		tm.tm_year + 1900,
		tm.tm_mon + 1,
		tm.tm_mday,
		tm.tm_hour,
		tm.tm_min,
		tm.tm_sec,
		file_time.tv_nsec);
    }
  }
  if(!isfinite(start) && start_time.tv_sec != 0 && file_time.tv_sec != 0){
    // Calculate starting offset from absolute time and file timestamp
    struct timespec offset = {0};
    offset.tv_sec = start_time.tv_sec - file_time.tv_sec;
    offset.tv_nsec = start_time.tv_nsec - file_time.tv_nsec;
    if(offset.tv_nsec < 0){
      offset.tv_sec--;
      offset.tv_nsec += 1000000000LL; // billion ns = 1 sec
    }
    start = offset.tv_sec + offset.tv_nsec * 1e-9;
    if(start < 0){
      fprintf(stdout,"Specified start time is before file start\n");
      exit(0);
    }
    fprintf(stdout,"calculated key-down offset at %'.3lf sec\n",start);
  } else if(isfinite(start) && file_time.tv_sec != 0){
    struct tm tm = {0};
    struct tm *tmp;
    time_t t = file_time.tv_sec + (long)floor(start);
    long ns = file_time.tv_nsec + llrint(1e9*fmod(start,1.0));
    if(ns > 1000000000){
      ns -= 1000000000;
      t++;
    }
    tmp = gmtime_r(&t,&tm);
    (void)tmp;
    assert(tmp != NULL);

    fprintf(stdout,"calculated key-down time %04d-%02d-%02d %02d:%02d:%02d.%09ld UTC\n",
	    tm.tm_year + 1900,
	    tm.tm_mon + 1,
	    tm.tm_mday,
	    tm.tm_hour,
	    tm.tm_min,
	    tm.tm_sec,
	    ns);

  }
  if(!isfinite(start))
    start = 0; // Default to start of file if nothing else sets it

  float complex *iq = NULL;

  int samprate = 0; // frames/sec
  int num_samp = 0; // stereo frames (IQ pairs)
  int num_chan = 0; // must be 2; could implement real support later
  if(load_wav((float **)&iq, &num_samp, &num_chan, &samprate, file, fd) == -1){
    fprintf(stdout,"can't load %s: %s\n",file,strerror(errno));
    exit(1);
  }
  if (iq == NULL) {
    fprintf(stdout, "load_wav returned no data\n");
    return 1;
  }
  if(num_chan != 2){
    // Could eventually support 1 channel (real signal) with real FFT
    fprintf(stdout,"%s: %s has %d channels; must be 2\n",argv[0],file,num_chan);
    exit(1);
  }

  double const fs = (double)samprate;
  if(freq_lo < -0.5 * fs || freq_lo > 0.5 * fs){
    freq_lo = -0.5 * fs;
    fprintf(stdout,"freq_lo forced to -Fs/2 = %.1lf\n",freq_lo);
  }
  if(freq_hi < freq_lo || freq_hi > 0.5 * fs){
    freq_hi = 0.5 * fs;
    fprintf(stdout,"freq_hi forced to +Fs/2 = %.1lf\n",freq_hi);
  }

  long start_sample = llround(start * fs);
  if(start_sample >= num_samp){
    fprintf(stdout,"start after end of file\n");
    exit(1);
  }
  if(duration_sec > start + (double)num_samp / fs)
    duration_sec = start + (double)num_samp / fs;
  if(hop_sec == 0)
    hop_sec = win_sec * 0.5; // default to 50% overlap

  fprintf(stdout,"samples %'d; time %'.1lf sec; samprate %'.1lf Hz\n",num_samp,(double)num_samp/fs,fs);
  fprintf(stdout,"analysis start @ %'.3lf sec, duration %'.3lf sec, freq_lo %'.3lf Hz, freq_hi %'.3lf Hz, win %'.3lf sec, hop %'.3lf sec\n",
	  start,duration_sec,freq_lo,freq_hi, win_sec, hop_sec);
  fprintf(stdout,"outlier thresh %'.3lf Hz, FFT oversample %d, min power %.1lf dB trim fraction %.2lf\n",
	  outlier_hz,oversample,min_rel_power_db,trim_fraction);

  if(Verbose)
    fprintf(stdout,"Buffer at %p + offset %'ld = %p\n",iq,start_sample,iq + start_sample);

  estimate_t *est = NULL;
  int n_est = 0;
  summary_t s = {0};

  /* Typical starting values:
     win_sec = 1.0
     hop_sec = 0.5
     nfft    = 4 * win_n  (zero padding helps interpolation behave more smoothly)
     search range narrow around expected residual frequency
     */
  int win_n = (int)llround(win_sec * fs);
  int nfft = oversample * win_n;

  if (estimate_track(iq + start_sample,(int64_t)llrint(duration_sec * fs),
		     fs,
		     win_sec,
		     hop_sec,
		     nfft,
		     freq_lo,
		     freq_hi,
		     min_rel_power_db,
		     outlier_hz,
		     trim_fraction,
		     &est,
		     &n_est,
		     &s) != 0) {

        fprintf(stdout, "estimate_track failed\n");
        free(iq);
        return 1;
  }
  printf("Good windows:             %'d / %'d\n", s.n_good,s.n_total);
  printf("Median freq:              %'.6lf Hz\n", s.median_hz);
  printf("Trimmed weighted mean:    %'.6lf Hz\n", s.trimmed_weighted_mean_hz);

  printf("Drift fit ref time:       %'.3lf s\n", s.drift_ref_time_sec);
  printf("Drift fit ref freq:       %'.6lf Hz\n", s.drift_ref_freq_hz);
  printf("Drift slope:              %'.6lf Hz/s\n", s.drift_slope_hz_per_sec);
  printf("Residual RMS:           %'.6lf Hz\n", s.residual_rms_hz);
  printf("Estimates:                %'d\n",n_est);
  printf("Best estimate: %'.6lf +/- %'.6lf Hz\n",
	 base_frequency + s.trimmed_weighted_mean_hz,
	 s.residual_rms_hz);


  if(Verbose){
    // Average power of good windows
    double pavg = 0;
    int n_valid = 0;
    for (int i = 0; i < n_est; i++) {
      if(est[i].valid){
	n_valid++;
	pavg += est[i].power;
      }
    }
    if(n_valid > 0){
      pavg = 10 * log10(pavg / n_valid);
      printf("avg power = %.1lf dB\n",pavg);
    }

    printf("  sec         Hz    dB   seq flags\n");

    for (int i = 0; i < n_est; i++) {
      double db = 10*log10(est[i].power);
      double rel_db = db - pavg;
      printf("%'5.1lf %'10.3lf %'5.1lf",
	     (double)est[i].time_sec,
	     (double)est[i].freq_hz,(double)rel_db);
      if(est[i].sequence != -1)
	printf(" %5d",est[i].sequence);
      else
	printf("      ");
      if(est[i].weak)
	printf(" weak");
      if(est[i].outlier)
	printf(" outlier");
      if(est[i].trimmed)
	printf(" trimmed");
      printf("\n");
    }
  }
  free(est);
  free(iq);
  return 0;
}
