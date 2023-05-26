#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "constants.h"
#include "util.h"
#include "cache.h"
#include "sieve.h"
#include "twin_primes.h"
#include "inverse_interpolate.h"
#include "real.h"
#include "mathl.h"

/******************************************************************************/
/*                                TWIN PRIMES                                 */
/******************************************************************************/

/* Twin prime counts (X * 10^Y to (X+1) * 10^Y). */
#if BITS_PER_WORD < 64
static const UV twin_steps[] =
  {58980,48427,45485,43861,42348,41457,40908,39984,39640,39222,
   373059,353109,341253,332437,326131,320567,315883,312511,309244,
   2963535,2822103,2734294,2673728,
  };
static const unsigned int twin_num_exponents = 3;
static const unsigned int twin_last_mult = 4;      /* 4000M */
#else
static const UV twin_steps[] =
  {58980,48427,45485,43861,42348,41457,40908,39984,39640,39222,
   373059,353109,341253,332437,326131,320567,315883,312511,309244,
   2963535,2822103,2734294,2673728,2626243,2585752,2554015,2527034,2501469,
   /* pi2(1e10,2e10) = 24096420;  pi2(2e10,3e10) = 23046519;  ...  */
   24096420,23046519,22401089,21946975,21590715,21300632,21060884,20854501,20665634,
   199708605,191801047,186932018,183404596,180694619,178477447,176604059,174989299,173597482,
   1682185723,1620989842,1583071291,1555660927,1534349481,1517031854,1502382532,1489745250, 1478662752,
   14364197903,13879821868,13578563641,13361034187,13191416949,13053013447,12936030624,12835090276, 12746487898,
   124078078589,120182602778,117753842540,115995331742,114622738809,113499818125,112551549250,111732637241,111012321565,
   1082549061370,1050759497170,1030883829367,1016473645857,1005206830409,995980796683,988183329733,981441437376,975508027029,
   9527651328494, 9264843314051, 9100153493509, 8980561036751, 8886953365929, 8810223086411, 8745329823109, 8689179566509, 8639748641098,
   84499489470819, 82302056642520, 80922166953330, 79918799449753, 79132610984280, 78487688897426, 77941865286827, 77469296499217, 77053075040105,
   754527610498466, 735967887462370, 724291736697048,
  };
static const unsigned int twin_num_exponents = 12;
static const unsigned int twin_last_mult = 4;      /* 4e18 */
#endif

UV twin_prime_count(UV n)
{
  return (n < 3) ? 0 : twin_prime_count_range(0,n);
}
UV twin_prime_count_range(UV beg, UV end)
{
  unsigned char* segment;
  UV sum = 0;

  /* First use the tables of #e# from 1e7 to 4e18. */
  if (beg <= 3 && end >= 10000000) {
    UV mult, exp, step = 0, base = 10000000;
    for (exp = 0; exp < twin_num_exponents && end >= base; exp++) {
      for (mult = 1; mult < 10 && end >= mult*base; mult++) {
        sum += twin_steps[step++];
        beg = mult*base;
        if (exp == twin_num_exponents-1 && mult >= twin_last_mult) break;
      }
      base *= 10;
    }
  }
  if (beg <= 3 && end >= 3) sum++;
  if (beg <= 5 && end >= 5) sum++;
  if (beg < 11) beg = 7;
  if (beg <= end) {
    /* Make end points odd */
    beg |= 1;
    end = (end-1) | 1;
    /* Cheesy way of counting the partial-byte edges */
    while ((beg % 30) != 1) {
      if (is_prime(beg) && is_prime(beg+2) && beg <= end) sum++;
      beg += 2;
    }
    while ((end % 30) != 29) {
      if (is_prime(end) && is_prime(end+2) && beg <= end) sum++;
      end -= 2;  if (beg > end) break;
    }
  }
  if (beg <= end) {
    UV seg_base, seg_low, seg_high;
    void* ctx = start_segment_primes(beg, end, &segment);
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      UV bytes = seg_high/30 - seg_low/30 + 1;
      unsigned char s, x;
      const unsigned char* sp = segment;
      const unsigned char* const spend = segment + bytes - 1;
      for (s = x = *sp;  sp++ < spend;  s = x) {
        x = *sp;
        if (!(s & 0x0C)) sum++;
        if (!(s & 0x30)) sum++;
        if (!(s & 0x80) && !(x & 0x01)) sum++;
      }
      x = is_prime(seg_high+2) ? 0x00 : 0xFF;
      if (!(s & 0x0C)) sum++;
      if (!(s & 0x30)) sum++;
      if (!(s & 0x80) && !(x & 0x01)) sum++;
    }
    end_segment_primes(ctx);
  }
  return sum;
}

/* See http://numbers.computation.free.fr/Constants/Primes/twin.pdf, page 5 */
/* Upper limit is in Wu, Acta Arith 114 (2004).  4.48857*x/(log(x)*log(x) */
/* Lichtman (2021) improved the limit: https://arxiv.org/pdf/2109.02851.pdf */
UV twin_prime_count_approx(UV n)
{
  /* Best would be another estimate for n < ~ 5000 */
  if (n < 2000) return twin_prime_count(n);
  {
    /* Sebah and Gourdon 2002 */
    const long double two_C2 = 1.32032363169373914785562422L;
    const long double two_over_log_two = 2.8853900817779268147198494L;
    long double ln = (long double) n;
    long double logn = logl(ln);
    long double li2 = Ei(logn) + two_over_log_two-ln/logn;
    /* Try to minimize MSE. */
    /* We compromise to prevent discontinuities. */
    if (n < 32000000) {
      long double fm;
      if      (n <    4000) fm = 0.2952;
      else if (n <    8000) fm = 0.3102;
      else if (n <   16000) fm = 0.3090;
      else if (n <   32000) fm = 0.3096;
      else if (n <   64000) fm = 0.3097;
      else if (n <  128000) fm = 0.3094;
      else if (n <  256000) fm = 0.3099;
      else if (n <  600000) fm = .3098 + (n-256000) * (.3056-.3098) / (600000-256000);
      else if (n < 1000000) fm = .3062 + (n-600000) * (.3042-.3062) / (1000000-600000);
      else if (n < 4000000) fm = .3067 + (n-1000000) * (.3041-.3067) / (4000000-1000000);
      else if (n <16000000) fm = .3041 + (n-4000000) * (.2983-.3041) / (16000000-4000000);
      else                  fm = .2983 + (n-16000000) * (.2961-.2983) / (32000000-16000000);
      li2 *= fm * logl(12+logn);
    }
    return (UV) (two_C2 * li2 + 0.5L);
  }
}


UV nth_twin_prime(UV n)
{
  unsigned char* segment;
  double dend;
  UV nth = 0;
  UV beg, end;

  if (n < 6) {
    switch (n) {
      case 0:  nth = 0; break;
      case 1:  nth = 3; break;
      case 2:  nth = 5; break;
      case 3:  nth =11; break;
      case 4:  nth =17; break;
      case 5:
      default: nth =29; break;
    }
    return nth;
  }

  end = UV_MAX - 16;
  dend = 800.0 + 1.01L * (double)nth_twin_prime_approx(n);
  if (dend < (double)end) end = (UV) dend;

  beg = 2;
  if (n > 58980) { /* Use twin_prime_count tables to accelerate if possible */
    UV mult, exp, step = 0, base = 10000000;
    for (exp = 0; exp < twin_num_exponents && end >= base; exp++) {
      for (mult = 1; mult < 10 && n > twin_steps[step]; mult++) {
        n -= twin_steps[step++];
        beg = mult*base;
        if (exp == twin_num_exponents-1 && mult >= twin_last_mult) break;
      }
      base *= 10;
    }
  }
  if (beg == 2) { beg = 31; n -= 5; }

  {
    UV seg_base, seg_low, seg_high;
    void* ctx = start_segment_primes(beg, end, &segment);
    while (n && next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      UV p, bytes = seg_high/30 - seg_low/30 + 1;
      UV s = ((UV)segment[0]) << 8;
      for (p = 0; p < bytes; p++) {
        s >>= 8;
        if (p+1 < bytes)                    s |= (((UV)segment[p+1]) << 8);
        else if (!is_prime(seg_high+2)) s |= 0xFF00;
        if (!(s & 0x000C) && !--n) { nth=seg_base+p*30+11; break; }
        if (!(s & 0x0030) && !--n) { nth=seg_base+p*30+17; break; }
        if (!(s & 0x0180) && !--n) { nth=seg_base+p*30+29; break; }
      }
    }
    end_segment_primes(ctx);
  }
  return nth;
}

UV nth_twin_prime_approx(UV n)
{
  long double fn = (long double) n;
  long double flogn = logl(n);
  long double fnlog2n = fn * flogn * flogn;
  UV lo, hi;

  if (n < 6)
    return nth_twin_prime(n);

  /* Binary search on the TPC estimate.
   * Good results require that the TPC estimate is both fast and accurate.
   * These bounds are good for the actual nth_twin_prime values.
   */
  lo = (UV) (0.9 * fnlog2n);
  hi = (UV) ( (n >= 1e16) ? (1.04 * fnlog2n) :
              (n >= 1e13) ? (1.10 * fnlog2n) :
              (n >= 1e7 ) ? (1.31 * fnlog2n) :
              (n >= 1200) ? (1.70 * fnlog2n) :
              (2.3 * fnlog2n + 5) );
  if (hi <= lo) hi = UV_MAX;
  return inverse_interpolate(lo, hi, n, &twin_prime_count_approx, 0);
}
