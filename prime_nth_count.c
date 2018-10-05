#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FUNC_popcnt 1
#define FUNC_isqrt 1
#define FUNC_icbrt 1
#include "ptypes.h"
#include "sieve.h"
#include "cache.h"
#include "lmo.h"
#include "constants.h"
#include "prime_nth_count.h"
#include "util.h"

#include <math.h>
#if _MSC_VER || defined(__IBMC__) || defined(__IBMCPP__) || (defined(__STDC_VERSION__) && __STDC_VERSION >= 199901L)
  /* math.h should give us these as functions or macros.
   *
   *  extern long double floorl(long double);
   *  extern long double ceill(long double);
   *  extern long double sqrtl(long double);
   *  extern long double logl(long double);
   */
#else
  #define floorl(x)   (long double) floor( (double) (x) )
  #define ceill(x)    (long double) ceil( (double) (x) )
  #define sqrtl(x)    (long double) sqrt( (double) (x) )
  #define logl(x)     (long double) log( (double) (x) )
#endif

#if defined(__GNUC__)
 #define word_unaligned(m,wordsize)  ((uintptr_t)m & (wordsize-1))
#else  /* uintptr_t is part of C99 */
 #define word_unaligned(m,wordsize)  ((unsigned int)m & (wordsize-1))
#endif

/* TODO: This data is duplicated in util.c. */
static const unsigned char prime_sieve30[] =
  {0x01,0x20,0x10,0x81,0x49,0x24,0xc2,0x06,0x2a,0xb0,0xe1,0x0c,0x15,0x59,0x12,
   0x61,0x19,0xf3,0x2c,0x2c,0xc4,0x22,0xa6,0x5a,0x95,0x98,0x6d,0x42,0x87,0xe1,
   0x59,0xa9,0xa9,0x1c,0x52,0xd2,0x21,0xd5,0xb3,0xaa,0x26,0x5c,0x0f,0x60,0xfc,
   0xab,0x5e,0x07,0xd1,0x02,0xbb,0x16,0x99,0x09,0xec,0xc5,0x47,0xb3,0xd4,0xc5,
   0xba,0xee,0x40,0xab,0x73,0x3e,0x85,0x4c,0x37,0x43,0x73,0xb0,0xde,0xa7,0x8e,
   0x8e,0x64,0x3e,0xe8,0x10,0xab,0x69,0xe5,0xf7,0x1a,0x7c,0x73,0xb9,0x8d,0x04,
   0x51,0x9a,0x6d,0x70,0xa7,0x78,0x2d,0x6d,0x27,0x7e,0x9a,0xd9,0x1c,0x5f,0xee,
   0xc7,0x38,0xd9,0xc3,0x7e,0x14,0x66,0x72,0xae,0x77,0xc1,0xdb,0x0c,0xcc,0xb2,
   0xa5,0x74,0xe3,0x58,0xd5,0x4b,0xa7,0xb3,0xb1,0xd9,0x09,0xe6,0x7d,0x23,0x7c,
   0x3c,0xd3,0x0e,0xc7,0xfd,0x4a,0x32,0x32,0xfd,0x4d,0xb5,0x6b,0xf3,0xa8,0xb3,
   0x85,0xcf,0xbc,0xf4,0x0e,0x34,0xbb,0x93,0xdb,0x07,0xe6,0xfe,0x6a,0x57,0xa3,
   0x8c,0x15,0x72,0xdb,0x69,0xd4,0xaf,0x59,0xdd,0xe1,0x3b,0x2e,0xb7,0xf9,0x2b,
   0xc5,0xd0,0x8b,0x63,0xf8,0x95,0xfa,0x77,0x40,0x97,0xea,0xd1,0x9f,0xaa,0x1c,
   0x48,0xae,0x67,0xf7,0xeb,0x79,0xa5,0x55,0xba,0xb2,0xb6,0x8f,0xd8,0x2d,0x6c,
   0x2a,0x35,0x54,0xfd,0x7c,0x9e,0xfa,0xdb,0x31,0x78,0xdd,0x3d,0x56,0x52,0xe7,
   0x73,0xb2,0x87,0x2e,0x76,0xe9,0x4f,0xa8,0x38,0x9d,0x5d,0x3f,0xcb,0xdb,0xad,
   0x51,0xa5,0xbf,0xcd,0x72,0xde,0xf7,0xbc,0xcb,0x49,0x2d,0x49,0x26,0xe6,0x1e,
   0x9f,0x98,0xe5,0xc6,0x9f,0x2f,0xbb,0x85,0x6b,0x65,0xf6,0x77,0x7c,0x57,0x8b,
   0xaa,0xef,0xd8,0x5e,0xa2,0x97,0xe1,0xdc,0x37,0xcd,0x1f,0xe6,0xfc,0xbb,0x8c,
   0xb7,0x4e,0xc7,0x3c,0x19,0xd5,0xa8,0x9e,0x67,0x4a,0xe3,0xf5,0x97,0x3a,0x7e,
   0x70,0x53,0xfd,0xd6,0xe5,0xb8,0x1c,0x6b,0xee,0xb1,0x9b,0xd1,0xeb,0x34,0xc2,
   0x23,0xeb,0x3a,0xf9,0xef,0x16,0xd6,0x4e,0x7d,0x16,0xcf,0xb8,0x1c,0xcb,0xe6,
   0x3c,0xda,0xf5,0xcf};
#define NPRIME_SIEVE30 (sizeof(prime_sieve30)/sizeof(prime_sieve30[0]))

static const unsigned short primes_small[] =
  {0,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
   101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,
   193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,
   293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,
   409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499};
#define NPRIMES_SMALL (sizeof(primes_small)/sizeof(primes_small[0]))


static const unsigned char byte_zeros[256] =
  {8,7,7,6,7,6,6,5,7,6,6,5,6,5,5,4,7,6,6,5,6,5,5,4,6,5,5,4,5,4,4,3,
   7,6,6,5,6,5,5,4,6,5,5,4,5,4,4,3,6,5,5,4,5,4,4,3,5,4,4,3,4,3,3,2,
   7,6,6,5,6,5,5,4,6,5,5,4,5,4,4,3,6,5,5,4,5,4,4,3,5,4,4,3,4,3,3,2,
   6,5,5,4,5,4,4,3,5,4,4,3,4,3,3,2,5,4,4,3,4,3,3,2,4,3,3,2,3,2,2,1,
   7,6,6,5,6,5,5,4,6,5,5,4,5,4,4,3,6,5,5,4,5,4,4,3,5,4,4,3,4,3,3,2,
   6,5,5,4,5,4,4,3,5,4,4,3,4,3,3,2,5,4,4,3,4,3,3,2,4,3,3,2,3,2,2,1,
   6,5,5,4,5,4,4,3,5,4,4,3,4,3,3,2,5,4,4,3,4,3,3,2,4,3,3,2,3,2,2,1,
   5,4,4,3,4,3,3,2,4,3,3,2,3,2,2,1,4,3,3,2,3,2,2,1,3,2,2,1,2,1,1,0};
static UV count_zero_bits(const unsigned char* m, UV nbytes)
{
  UV count = 0;
#if BITS_PER_WORD == 64
  if (nbytes >= 16) {
    while ( word_unaligned(m,sizeof(UV)) && nbytes--)
      count += byte_zeros[*m++];
    if (nbytes >= 8) {
      UV* wordptr = (UV*)m;
      UV nwords = nbytes / 8;
      UV nzeros = nwords * 64;
      m += nwords * 8;
      nbytes %= 8;
      while (nwords--)
        nzeros -= popcnt(*wordptr++);
      count += nzeros;
    }
  }
#endif
  while (nbytes--)
    count += byte_zeros[*m++];
  return count;
}

/* Given a sieve of size nbytes, walk it counting zeros (primes) until:
 *
 * (1) we counted them all: return the count, which will be less than maxcount.
 *
 * (2) we hit maxcount: set position to the index of the maxcount'th prime
 *     and return count (which will be equal to maxcount).
 */
static UV count_segment_maxcount(const unsigned char* sieve, UV base, UV nbytes, UV maxcount, UV* pos)
{
  UV count = 0;
  UV byte = 0;
  const unsigned char* sieveptr = sieve;
  const unsigned char* maxsieve = sieve + nbytes;

  MPUassert(sieve != 0, "count_segment_maxcount incorrect args");
  MPUassert(pos != 0, "count_segment_maxcount incorrect args");
  *pos = 0;
  if ( (nbytes == 0) || (maxcount == 0) )
    return 0;

  /* Do fixed-length word counts to start, with possible overcounting */
  while ((count+64) < maxcount && sieveptr < maxsieve) {
    UV top = base + 3*maxcount;
    UV div = (top <       8000) ? 8 :     /* 8 cannot overcount */
             (top <    1000000) ? 4 :
             (top <   10000000) ? 3 : 2;
    UV minbytes = (maxcount-count)/div;
    if (minbytes > (UV)(maxsieve-sieveptr)) minbytes = maxsieve-sieveptr;
    count += count_zero_bits(sieveptr, minbytes);
    sieveptr += minbytes;
  }
  /* Count until we reach the end or >= maxcount */
  while ( (sieveptr < maxsieve) && (count < maxcount) )
    count += byte_zeros[*sieveptr++];
  /* If we went too far, back up. */
  while (count >= maxcount)
    count -= byte_zeros[*--sieveptr];
  /* We counted this many bytes */
  byte = sieveptr - sieve;

  MPUassert(count < maxcount, "count_segment_maxcount wrong count");

  if (byte == nbytes)
    return count;

  /* The result is somewhere in the next byte */
  START_DO_FOR_EACH_SIEVE_PRIME(sieve, 0, byte*30+1, nbytes*30-1)
    if (++count == maxcount)  { *pos = p; return count; }
  END_DO_FOR_EACH_SIEVE_PRIME;

  MPUassert(0, "count_segment_maxcount failure");
  return 0;
}

/* Given a sieve of size nbytes, counting zeros (primes) but excluding the
 * areas outside lowp and highp.
 */
static UV count_segment_ranged(const unsigned char* sieve, UV nbytes, UV lowp, UV highp)
{
  UV count, hi_d, lo_d, lo_m;

  MPUassert( sieve != 0, "count_segment_ranged incorrect args");
  if (nbytes == 0) return 0;

  count = 0;
  hi_d = highp/30;

  if (hi_d >= nbytes) {
    hi_d = nbytes-1;
    highp = hi_d*30+29;
  }

  if (highp < lowp)
    return 0;

#if 0
  /* Dead simple way */
  START_DO_FOR_EACH_SIEVE_PRIME(sieve, 0, lowp, highp)
    count++;
  END_DO_FOR_EACH_SIEVE_PRIME;
  return count;
#endif

  lo_d = lowp/30;
  lo_m = lowp - lo_d*30;
  /* Count first fragment */
  if (lo_m > 1) {
    UV upper = (highp <= (lo_d*30+29)) ? highp : (lo_d*30+29);
    START_DO_FOR_EACH_SIEVE_PRIME(sieve, 0, lowp, upper)
      count++;
    END_DO_FOR_EACH_SIEVE_PRIME;
    lowp = upper+2;
    lo_d = lowp/30;
  }
  if (highp < lowp)
    return count;

  /* Count bytes in the middle */
  {
    UV hi_m = highp - hi_d*30;
    UV count_bytes = hi_d - lo_d + (hi_m == 29);
    if (count_bytes > 0) {
      count += count_zero_bits(sieve+lo_d, count_bytes);
      lowp += 30*count_bytes;
    }
  }
  if (highp < lowp)
    return count;

  /* Count last fragment */
  START_DO_FOR_EACH_SIEVE_PRIME(sieve, 0, lowp, highp)
    count++;
  END_DO_FOR_EACH_SIEVE_PRIME;

  return count;
}


/*
 * The pi(x) prime count functions.  prime_count(x) gives an exact number,
 * but requires determining all the primes up to x, so will be much slower.
 *
 * prime_count_lower(x) and prime_count_upper(x) give lower and upper limits,
 * which will bound the exact value.  These bounds should be fairly tight.
 *
 *    pi_upper(x) - pi(x)                      pi_lower(x) - pi(x)
 *    <    10   for x <         5_371          <    10   for x <        9_437
 *    <    50   for x <       295_816          <    50   for x <      136_993
 *    <   100   for x <     1_761_655          <   100   for x <      909_911
 *    <   200   for x <     9_987_821          <   200   for x <    8_787_901
 *    <   400   for x <    34_762_891          <   400   for x <   30_332_723
 *    <  1000   for x <   372_748_528          <  1000   for x <  233_000_533
 *    <  5000   for x < 1_882_595_905          <  5000   for x <  over 4300M
 *
 * The average of the upper and lower bounds is within 9 for all x < 15809, and
 * within 50 for all x < 1_763_367.
 *
 * It is common to use the following Chebyshev inequality for x >= 17:
 *    1*x/logx <-> 1.25506*x/logx
 * but this gives terribly loose bounds.
 *
 * Rosser and Schoenfeld's bound for x >= 67 of
 *    x/(logx-1/2) <-> x/(logx-3/2)
 * is much tighter.  These bounds can be tightened even more.
 *
 * The formulas of Dusart for higher x are better yet.  I recommend the paper
 * by Burde for further information.  Dusart's thesis is also a good resource.
 *
 * I have tweaked the bounds formulas for small (under 70_000M) numbers so they
 * are tighter.  These bounds are verified via trial.  The Dusart bounds
 * (1.8 and 2.51) are used for larger numbers since those are proven.
 *
 */

#include "prime_count_tables.h"

UV segment_prime_count(UV low, UV high)
{
  const unsigned char* cache_sieve;
  unsigned char* segment;
  UV segment_size, low_d, high_d;
  UV count = 0;

  if ((low <= 2) && (high >= 2)) count++;
  if ((low <= 3) && (high >= 3)) count++;
  if ((low <= 5) && (high >= 5)) count++;
  if (low < 7)  low = 7;

  if (low > high)  return count;

#if !defined(BENCH_SEGCOUNT)
  if (low == 7 && high <= 30*NPRIME_SIEVE30) {
    count += count_segment_ranged(prime_sieve30, NPRIME_SIEVE30, low, high);
    return count;
  }

  /* If we have sparse prime count tables, use them here.  These will adjust
   * 'low' and 'count' appropriately for a value slightly less than ours.
   * This should leave just a small amount of sieving left.  They stop at
   * some point, e.g. 3000M, so we'll get the answer to that point then have
   * to sieve all the rest.  We should be using LMO or Lehmer much earlier. */
#ifdef APPLY_TABLES
  APPLY_TABLES
#endif
#endif

  low_d = low/30;
  high_d = high/30;

  /* Count full bytes only -- no fragments from primary cache */
  segment_size = get_prime_cache(0, &cache_sieve) / 30;
  if (segment_size < high_d) {
    /* Expand sieve to sqrt(n) */
    UV endp = (high_d >= (UV_MAX/30))  ?  UV_MAX-2  :  30*high_d+29;
    release_prime_cache(cache_sieve);
    segment_size = get_prime_cache( isqrt(endp) + 1 , &cache_sieve) / 30;
  }

  if ( (segment_size > 0) && (low_d <= segment_size) ) {
    /* Count all the primes in the primary cache in our range */
    count += count_segment_ranged(cache_sieve, segment_size, low, high);

    if (high_d < segment_size) {
      release_prime_cache(cache_sieve);
      return count;
    }

    low_d = segment_size;
    if (30*low_d > low)  low = 30*low_d;
  }
  release_prime_cache(cache_sieve);

  /* More primes needed.  Repeatedly segment sieve. */
  {
    void* ctx = start_segment_primes(low, high, &segment);
    UV seg_base, seg_low, seg_high;
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      segment_size = seg_high/30 - seg_low/30 + 1;
      count += count_segment_ranged(segment, segment_size, seg_low-seg_base, seg_high-seg_base);
    }
    end_segment_primes(ctx);
  }

  return count;
}

UV prime_count(UV lo, UV hi)
{
  if (lo > hi || hi < 2)
    return 0;

#if defined(BENCH_SEGCOUNT)
  return segment_prime_count(lo, hi);
#endif

  /* We use table acceleration so this is preferable for small inputs */
  if (hi < _MPU_LMO_CROSSOVER)  return segment_prime_count(lo, hi);

  { /* Rough empirical threshold for when segment faster than LMO */
    UV range_threshold = hi / (isqrt(hi)/200);
    if ( (hi-lo+1) < range_threshold )
      return segment_prime_count(lo, hi);
  }
  return LMO_prime_count(hi) - ((lo < 2) ? 0 : LMO_prime_count(lo-1));
}

UV prime_count_approx(UV n)
{
  if (n < 3000000) return segment_prime_count(2, n);
  return (UV) (RiemannR( (long double) n ) + 0.5 );
}

/* See http://numbers.computation.free.fr/Constants/Primes/twin.pdf, page 5 */
/* Upper limit is in Wu, Acta Arith 114 (2004).  4.48857*x/(log(x)*log(x) */
UV twin_prime_count_approx(UV n)
{
  /* Best would be another estimate for n < ~ 5000 */
  if (n < 2000) return twin_prime_count(3,n);
  {
    /* Sebah and Gourdon 2002 */
    const long double two_C2 = 1.32032363169373914785562422L;
    const long double two_over_log_two = 2.8853900817779268147198494L;
    long double ln = (long double) n;
    long double logn = logl(ln);
    long double li2 = Ei(logn) + two_over_log_two-ln/logn;
    /* try to minimize MSE */
    if (n < 32000000) {
      long double fm;
      if      (n <    4000) fm = 0.2952;
      else if (n <    8000) fm = 0.3152;
      else if (n <   16000) fm = 0.3090;
      else if (n <   32000) fm = 0.3096;
      else if (n <   64000) fm = 0.3100;
      else if (n <  128000) fm = 0.3089;
      else if (n <  256000) fm = 0.3099;
      else if (n <  600000) fm = .3091 + (n-256000) * (.3056-.3091) / (600000-256000);
      else if (n < 1000000) fm = .3062 + (n-600000) * (.3042-.3062) / (1000000-600000);
      else if (n < 4000000) fm = .3067 + (n-1000000) * (.3041-.3067) / (4000000-1000000);
      else if (n <16000000) fm = .3033 + (n-4000000) * (.2983-.3033) / (16000000-4000000);
      else                  fm = .2980 + (n-16000000) * (.2965-.2980) / (32000000-16000000);
      li2 *= fm * logl(12+logn);
    }
    return (UV) (two_C2 * li2 + 0.5L);
  }
}

UV prime_count_lower(UV n)
{
  long double fn, fl1, fl2, lower, a;

  if (n < 33000) return segment_prime_count(2, n);

  fn  = (long double) n;
  fl1 = logl(n);
  fl2 = fl1 * fl1;

  /* Axler 2014: https://arxiv.org/abs/1409.1780  (v7 2016), Cor 3.6
   * show variations of this. */

  if (n <= 300000) { /* Quite accurate and avoids calling Li for speed. */
    a = (n < 70200) ? 947 : (n < 176000) ? 904 : 829;
    lower = fn / (fl1 - 1 - 1/fl1 - 2.85/fl2 - 13.15/(fl1*fl2) + a/(fl2*fl2));
  } else if (n < UVCONST(4000000000)) {
    /* Loose enough that FP differences in Li(n) should be ok. */
    a = (n <     88783) ?   4.0L
      : (n <    300000) ?  -3.0L
      : (n <    303000) ?   5.0L
      : (n <   1100000) ?  -7.0L
      : (n <   4500000) ? -37.0L
      : (n <  10200000) ? -70.0L
      : (n <  36900000) ? -53.0L
      : (n <  38100000) ? -29.0L
      :                   -84.0L;
    lower = Li(fn) - (sqrtl(fn)/fl1) * (1.94L + 2.50L/fl1 + a/fl2);
  } else if (fn < 1e19) {          /* Büthe 2015 1.9      1511.02032v1.pdf */
    lower = Li(fn) - (sqrtl(fn)/fl1) * (1.94L + 3.88L/fl1 + 27.57L/fl2);
  } else {                         /* Büthe 2014 v3 7.2   1410.7015v3.pdf */
    lower = Li(fn) - fl1*sqrtl(fn)/25.132741228718345907701147L;
  }
  return (UV) ceill(lower);
}

typedef struct {
  UV thresh;
  float aval;
} thresh_t;

static const thresh_t _upper_thresh[] = {
  {     59000, 2.48 },
  {    355991, 2.54 },
  {   3550000, 2.51 },
  {   3560000, 2.49 },
  {   5000000, 2.48 },
  {   8000000, 2.47 },
  {  13000000, 2.46 },
  {  18000000, 2.45 },
  {  31000000, 2.44 },
  {  41000000, 2.43 },
  {  48000000, 2.42 },
  { 119000000, 2.41 },
  { 182000000, 2.40 },
  { 192000000, 2.395 },
  { 213000000, 2.390 },
  { 271000000, 2.385 },
  { 322000000, 2.380 },
  { 400000000, 2.375 },
  { 510000000, 2.370 },
  { 682000000, 2.367 },
  { UVCONST(2953652287), 2.362 }
};
#define NUPPER_THRESH (sizeof(_upper_thresh)/sizeof(_upper_thresh[0]))

UV prime_count_upper(UV n)
{
  int i;
  long double fn, fl1, fl2, upper, a;

  if (n < 33000) return segment_prime_count(2, n);

  fn  = (long double) n;
  fl1 = logl(n);
  fl2 = fl1 * fl1;

  /* Axler 2014: https://arxiv.org/abs/1409.1780  (v7 2016), Cor 3.5
   *
   * upper = fn/(fl1-1.0L-1.0L/fl1-3.35L/fl2-12.65L/(fl2*fl1)-89.6L/(fl2*fl2));
   * return (UV) floorl(upper);
   */

  if (BITS_PER_WORD == 32 || fn <= 821800000.0) {  /* Dusart 2010, page 2 */
    for (i = 0; i < (int)NUPPER_THRESH; i++)
      if (n < _upper_thresh[i].thresh)
        break;
    a = (i < (int)NUPPER_THRESH)  ?  _upper_thresh[i].aval  :  2.334L;
    upper = fn/fl1 * (1.0L + 1.0L/fl1 + a/fl2);
  } else if (fn < 1e19) {        /* Büthe 2015 1.10 Skewes number lower limit */
    a = (fn <   1100000000.0) ? 0.032    /* Empirical */
      : (fn <  10010000000.0) ? 0.027    /* Empirical */
      : (fn < 101260000000.0) ? 0.021    /* Empirical */
                              : 0.0;
    upper = Li(fn) - a * fl1*sqrtl(fn)/25.132741228718345907701147L;
  } else {                       /* Büthe 2014 7.4 */
    upper = Li(fn) + fl1*sqrtl(fn)/25.132741228718345907701147L;
  }
  return (UV) floorl(upper);
}

static void simple_nth_limits(UV *lo, UV *hi, long double n, long double logn, long double loglogn) {
  const long double a = (n < 228) ? .6483 : (n < 948) ? .8032 : (n < 2195) ? .8800 : (n < 39017) ? .9019 : .9484;
  *lo = n * (logn + loglogn - 1.0 + ((loglogn-2.10)/logn));
  *hi = n * (logn + loglogn - a);
  if (*hi < *lo) *hi = MPU_MAX_PRIME;
}

/* The nth prime will be less or equal to this number */
UV nth_prime_upper(UV n)
{
  long double fn, flogn, flog2n, upper;

  if (n < NPRIMES_SMALL)
    return primes_small[n];

  fn     = (long double) n;
  flogn  = logl(n);
  flog2n = logl(flogn);    /* Note distinction between log_2(n) and log^2(n) */

  if (n < 688383) {
    UV lo,hi;
    simple_nth_limits(&lo, &hi, fn, flogn, flog2n);
    while (lo < hi) {
      UV mid = lo + (hi-lo)/2;
      if (prime_count_lower(mid) < n) lo = mid+1;
      else hi = mid;
    }
    return lo;
  }

  /* Dusart 2010 page 2 */
  upper = fn * (flogn + flog2n - 1.0 + ((flog2n-2.00)/flogn));
  if        (n >= 46254381) {
     /* Axler 2017 http://arxiv.org/pdf/1706.03651.pdf Corollary 1.2 */
    upper -= fn * ((flog2n*flog2n-6*flog2n+10.667)/(2*flogn*flogn));
  } else if (n >=  8009824) {
    /* Axler 2013 page viii Korollar G */
    upper -= fn * ((flog2n*flog2n-6*flog2n+10.273)/(2*flogn*flogn));
  }

  if (upper >= (long double)UV_MAX) {
    if (n <= MPU_MAX_PRIME_IDX) return MPU_MAX_PRIME;
    croak("nth_prime_upper(%"UVuf") overflow", n);
  }

  return (UV) floorl(upper);
}

/* The nth prime will be greater than or equal to this number */
UV nth_prime_lower(UV n)
{
  double fn, flogn, flog2n, lower;

  if (n < NPRIMES_SMALL)
    return primes_small[n];

  fn     = (double) n;
  flogn  = log(n);
  flog2n = log(flogn);

  /* For small values, do a binary search on the inverse prime count */
  if (n < 2000000) {
    UV lo,hi;
    simple_nth_limits(&lo, &hi, fn, flogn, flog2n);
    while (lo < hi) {
      UV mid = lo + (hi-lo)/2;
      if (prime_count_upper(mid) < n) lo = mid+1;
      else                            hi = mid;
    }
    return lo;
  }

  { /* Axler 2017 http://arxiv.org/pdf/1706.03651.pdf Corollary 1.4 */
    double b1 = (n < 56000000)  ?  11.200  :  11.508;
    lower = fn * (flogn + flog2n-1.0 + ((flog2n-2.00)/flogn) - ((flog2n*flog2n-6*flog2n+b1)/(2*flogn*flogn)));
  }

  return (UV) ceill(lower);
}

UV nth_prime_approx(UV n)
{
  return (n < NPRIMES_SMALL)  ?  primes_small[n]  :  inverse_R(n);
}


UV nth_prime(UV n)
{
  const unsigned char* cache_sieve;
  unsigned char* segment;
  UV upper_limit, segbase, segment_size, p, count, target;

  /* If very small, return the table entry */
  if (n < NPRIMES_SMALL)
    return primes_small[n];

  /* Determine a bound on the nth prime.  We know it comes before this. */
  upper_limit = nth_prime_upper(n);
  MPUassert(upper_limit > 0, "nth_prime got an upper limit of 0");
  p = count = 0;
  target = n-3;

  /* For relatively small values, generate a sieve and count the results.
   *
   * For larger values, compute an approximate low estimate, use our fast
   * prime count, then segment sieve forwards or backwards for the rest.
   */
  if (upper_limit <= get_prime_cache(0, 0) || upper_limit <= 32*1024*30) {
    /* Generate a sieve and count. */
    segment_size = get_prime_cache(upper_limit, &cache_sieve) / 30;
    /* Count up everything in the cached sieve. */
    if (segment_size > 0)
      count += count_segment_maxcount(cache_sieve, 0, segment_size, target, &p);
    release_prime_cache(cache_sieve);
  } else {
    /* A binary search on RiemannR is nice, but ends up either often being
     * being higher (requiring going backwards) or biased and then far too
     * low.  Using the inverse Li is easier and more consistent. */
    UV lower_limit = inverse_li(n);
    /* For even better performance, add in half the usual correction, which
     * will get us even closer, so even less sieving required.  However, it
     * is now possible to get a result higher than the value, so we'll need
     * to handle that case.  It still ends up being a better deal than R,
     * given that we don't have a fast backward sieve. */
    lower_limit += inverse_li(isqrt(n))/4;
    segment_size = lower_limit / 30;
    lower_limit = 30 * segment_size - 1;
    count = prime_count(2,lower_limit);

    /* printf("We've estimated %lu too %s.\n", (count>n)?count-n:n-count, (count>n)?"FAR":"little"); */
    /* printf("Our limit %lu %s a prime\n", lower_limit, is_prime(lower_limit) ? "is" : "is not"); */

    if (count >= n) { /* Too far.  Walk backwards */
      if (is_prime(lower_limit)) count--;
      for (p = 0; p <= (count-n); p++)
        lower_limit = prev_prime(lower_limit);
      return lower_limit;
    }
    count -= 3;

    /* Make sure the segment siever won't have to keep resieving. */
    prime_precalc(isqrt(upper_limit));
  }

  if (count == target)
    return p;

  /* Start segment sieving.  Get memory to sieve into. */
  segbase = segment_size;
  segment = get_prime_segment(&segment_size);

  while (count < target) {
    /* Limit the segment size if we know the answer comes earlier */
    if ( (30*(segbase+segment_size)+29) > upper_limit )
      segment_size = (upper_limit - segbase*30 + 30) / 30;

    /* Do the actual sieving in the range */
    sieve_segment(segment, segbase, segbase + segment_size-1);

    /* Count up everything in this segment */
    count += count_segment_maxcount(segment, 30*segbase, segment_size, target-count, &p);

    if (count < target)
      segbase += segment_size;
  }
  release_prime_segment(segment);
  MPUassert(count == target, "nth_prime got incorrect count");
  return ( (segbase*30) + p );
}

/******************************************************************************/
/*                                TWIN PRIMES                                 */
/******************************************************************************/

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
static const unsigned int twin_last_mult = 4;      /* 4e19 */
#endif

UV twin_prime_count(UV beg, UV end)
{
  unsigned char* segment;
  UV sum = 0;

  /* First use the tables of #e# from 1e7 to 2e16. */
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
      unsigned char s;
      const unsigned char* sp = segment;
      const unsigned char* const spend = segment + bytes - 1;
      while (sp < spend) {
        s = *sp++;
        if (!(s & 0x0C)) sum++;
        if (!(s & 0x30)) sum++;
        if (!(s & 0x80) && !(*sp & 0x01)) sum++;
      }
      s = *sp;
      if (!(s & 0x0C)) sum++;
      if (!(s & 0x30)) sum++;
      if (!(s & 0x80) && is_prime(seg_high+2)) sum++;
    }
    end_segment_primes(ctx);
  }
  return sum;
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
  while (lo < hi) {
    UV mid = lo + (hi-lo)/2;
    if (twin_prime_count_approx(mid) < fn) lo = mid+1;
    else                                   hi = mid;
  }
  return lo;
}

/******************************************************************************/
/*                                   SUMS                                     */
/******************************************************************************/

/* The fastest way to compute the sum of primes is using a combinatorial
 * algorithm such as Deleglise 2012.  Since this code is purely native,
 * it will overflow a 64-bit result quite quickly.  Hence a relatively small
 * table plus sum over sieved primes works quite well.
 *
 * The following info is useful if we ever return 128-bit results or for a
 * GMP implementation.
 *
 * Combinatorial sum of primes < n.  Call with phisum(n, isqrt(n)).
 * Needs optimization, either caching, Lehmer, or LMO.
 * http://mathoverflow.net/questions/81443/fastest-algorithm-to-compute-the-sum-of-primes
 * http://www.ams.org/journals/mcom/2009-78-268/S0025-5718-09-02249-2/S0025-5718-09-02249-2.pdf
 * http://mathematica.stackexchange.com/questions/80291/efficient-way-to-sum-all-the-primes-below-n-million-in-mathematica
 * Deleglise 2012, page 27, simple Meissel:
 *   y = x^1/3
 *   a = Pi(y)
 *   Pi_f(x) = phisum(x,a) + Pi_f(y) - 1 - P_2(x,a)
 *   P_2(x,a) = sum prime p : y < p <= sqrt(x) of f(p) * Pi_f(x/p) -
 *              sum prime p : y < p <= sqrt(x) of f(p) * Pi_f(p-1)
 */

static const unsigned char byte_sum[256] =
  {120,119,113,112,109,108,102,101,107,106,100,99,96,95,89,88,103,102,96,95,92,
   91,85,84,90,89,83,82,79,78,72,71,101,100,94,93,90,89,83,82,88,87,81,80,77,
   76,70,69,84,83,77,76,73,72,66,65,71,70,64,63,60,59,53,52,97,96,90,89,86,85,
   79,78,84,83,77,76,73,72,66,65,80,79,73,72,69,68,62,61,67,66,60,59,56,55,49,
   48,78,77,71,70,67,66,60,59,65,64,58,57,54,53,47,46,61,60,54,53,50,49,43,42,
   48,47,41,40,37,36,30,29,91,90,84,83,80,79,73,72,78,77,71,70,67,66,60,59,74,
   73,67,66,63,62,56,55,61,60,54,53,50,49,43,42,72,71,65,64,61,60,54,53,59,58,
   52,51,48,47,41,40,55,54,48,47,44,43,37,36,42,41,35,34,31,30,24,23,68,67,61,
   60,57,56,50,49,55,54,48,47,44,43,37,36,51,50,44,43,40,39,33,32,38,37,31,30,
   27,26,20,19,49,48,42,41,38,37,31,30,36,35,29,28,25,24,18,17,32,31,25,24,21,
   20,14,13,19,18,12,11,8,7,1,0};

#if BITS_PER_WORD == 64
/* We have a much more limited range, so use a fixed interval.  We should be
 * able to get any 64-bit sum in under a half-second. */
static const UV sum_table_2e8[] =
  {1075207199997324,3071230303170813,4990865886639877,6872723092050268,8729485610396243,10566436676784677,12388862798895708,14198556341669206,15997206121881531,17783028661796383,19566685687136351,21339485298848693,23108856419719148,
   24861364231151903,26619321031799321,28368484289421890,30110050320271201,31856321671656548,33592089385327108,35316546074029522,37040262208390735,38774260466286299,40490125006181147,42207686658844380,43915802985817228,45635106002281013,
   47337822860157465,49047713696453759,50750666660265584,52449748364487290,54152689180758005,55832433395290183,57540651847418233,59224867245128289,60907462954737468,62597192477315868,64283665223856098,65961576139329367,67641982565760928,
   69339211720915217,71006044680007261,72690896543747616,74358564592509127,76016548794894677,77694517638354266,79351385193517953,81053240048141953,82698120948724835,84380724263091726,86028655116421543,87679091888973563,89348007111430334,
   90995902774878695,92678527127292212,94313220293410120,95988730932107432,97603162494502485,99310622699836698,100935243057337310,102572075478649557,104236362884241550,105885045921116836,107546170993472638,109163445284201278,
   110835950755374921,112461991135144669,114116351921245042,115740770232532531,117408250788520189,119007914428335965,120652479429703269,122317415246500401,123951466213858688,125596789655927842,127204379051939418,128867944265073217,
   130480037123800711,132121840147764197,133752985360747726,135365954823762234,137014594650995101,138614165689305879,140269121741383097,141915099618762647,143529289083557618,145146413750649432,146751434858695468,148397902396643807,
   149990139346918801,151661665434334577,153236861034424304,154885985064643097,156500983286383741,158120868946747299,159735201435796748,161399264792716319,162999489977602579,164566400448130092,166219688860475191,167836981098849796,
   169447127305804401,171078187147848898,172678849082290997,174284436375728242,175918609754056455,177525046501311788,179125593738290153,180765176633753371,182338473848291683,183966529541155489,185585792988238475,187131988176321434,
   188797837140841381,190397649440649965,191981841583560122,193609739194967419,195166830650558070,196865965063113041,198400070713177440,200057161591648721,201621899486413406,203238279253414934,204790684829891896,206407676204061001,
   208061050481364659,209641606658938873,211192088300183855,212855420483750498,214394145510853736,216036806225784861,217628995137940563,219277567478725189,220833877268454872,222430818525363309,224007307616922530,225640739533952807,
   227213096159236934,228853318075566255,230401824696558125,231961445347821085,233593317860593895,235124654760954338,236777716068869769,238431514923528303,239965003913481640,241515977959535845,243129874530821395};
#define N_SUM_TABLE  (sizeof(sum_table_2e8)/sizeof(sum_table_2e8[0]))
#endif

/* Add n to the double-word hi,lo */
#define ADD_128(hi, lo, n)  \
  do {  UV _n = n; \
        if (_n > (UV_MAX-lo)) { hi++; if (hi == 0) overflow = 1; } \
        lo += _n;   } while (0)
#define SET_128(hi, lo, n) \
  do { hi = (UV) (((n) >> 64) & UV_MAX); \
       lo = (UV) (((n)      ) & UV_MAX); } while (0)

/* Legendre method for prime sum */
int sum_primes128(UV n, UV *hi_sum, UV *lo_sum) {
#if BITS_PER_WORD == 64 && HAVE_UINT128
  uint128_t *V, *S;
  UV j, k, r = isqrt(n), r2 = r + n/(r+1);

  New(0, V, r2+1, uint128_t);
  New(0, S, r2+1, uint128_t);
  for (k = 0; k <= r2; k++) {
    uint128_t v = (k <= r)  ?  k  :  n/(r2-k+1);
    V[k] = v;
    S[k] = (v*(v+1))/2 - 1;
  }

  START_DO_FOR_EACH_PRIME(2, r) {
    uint128_t a, b, sp = S[p-1], p2 = ((uint128_t)p) * p;
    for (j = k-1; j > 1 && V[j] >= p2; j--) {
      a = V[j], b = a/p;
      if (a > r) a = r2 - n/a + 1;
      if (b > r) b = r2 - n/b + 1;
      S[a] -= p * (S[b] - sp);   /* sp = sum of primes less than p */
    }
  } END_DO_FOR_EACH_PRIME;
  SET_128(*hi_sum, *lo_sum, S[r2]);
  Safefree(V);
  Safefree(S);
  return 1;
#else
  return 0;
#endif
}

int sum_primes(UV low, UV high, UV *return_sum) {
  UV sum = 0;
  int overflow = 0;

  if ((low <= 2) && (high >= 2)) sum += 2;
  if ((low <= 3) && (high >= 3)) sum += 3;
  if ((low <= 5) && (high >= 5)) sum += 5;
  if (low < 7) low = 7;

  /* If we know the range will overflow, return now */
#if BITS_PER_WORD == 64
  if (low == 7 && high >= 29505444491)  return 0;
  if (low >= 1e10 && (high-low) >= 32e9) return 0;
  if (low >= 1e13 && (high-low) >=  5e7) return 0;
#else
  if (low == 7 && high >= 323381)  return 0;
#endif

#if 1 && BITS_PER_WORD == 64    /* Tables */
  if (low == 7 && high >= 2e8) {
    UV step;
    for (step = 1; high >= (step * 2e8) && step < N_SUM_TABLE; step++) {
      sum += sum_table_2e8[step-1];
      low = step * 2e8;
    }
  }
#endif

  if (low <= high) {
    unsigned char* segment;
    UV seg_base, seg_low, seg_high;
    void* ctx = start_segment_primes(low, high, &segment);
    while (!overflow && next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      UV bytes = seg_high/30 - seg_low/30 + 1;
      unsigned char s;
      unsigned char* sp = segment;
      unsigned char* const spend = segment + bytes - 1;
      UV i, p, pbase = 30*(seg_low/30);

      /* Clear primes before and after our range */
      p = pbase;
      for (i = 0; i < 8 && p+wheel30[i] < low; i++)
        if ( (*sp & (1<<i)) == 0 )
          *sp |= (1 << i);

      p = 30*(seg_high/30);
      for (i = 0; i < 8;  i++)
        if ( (*spend & (1<<i)) == 0 && p+wheel30[i] > high )
          *spend |= (1 << i);

      while (sp <= spend) {
        s = *sp++;
        if (sum < (UV_MAX >> 3) && pbase < (UV_MAX >> 5)) {
          /* sum block of 8 all at once */
          sum += pbase * byte_zeros[s] + byte_sum[s];
        } else {
          /* sum block of 8, checking for overflow at each step */
          for (i = 0; i < byte_zeros[s]; i++) {
            if (sum+pbase < sum) overflow = 1;
            sum += pbase;
          }
          if (sum+byte_sum[s] < sum) overflow = 1;
          sum += byte_sum[s];
          if (overflow) break;
        }
        pbase += 30;
      }
    }
    end_segment_primes(ctx);
  }
  if (!overflow && return_sum != 0)  *return_sum = sum;
  return !overflow;
}

double ramanujan_sa_gn(UV un)
{
  long double n = (long double) un;
  long double logn = logl(n);
  long double log2 = logl(2);

  return (double)( (logn + logl(logn) - log2 - 0.5) / (log2 + 0.5) );
}
