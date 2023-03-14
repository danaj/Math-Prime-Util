#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FUNC_popcnt 1
#define FUNC_isqrt 1
#include "ptypes.h"
#include "sieve.h"
#include "cache.h"
#include "lmo.h"
#include "constants.h"
#include "prime_counts.h"
#include "util.h"
#include "real.h"
#include "mathl.h"

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

UV prime_count_range(UV lo, UV hi)
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

UV prime_count(UV n)
{
  if (n < 2) return 0;

  /* We use table acceleration so this is preferable for small inputs */
  if (n < _MPU_LMO_CROSSOVER)  return segment_prime_count(0, n);

  return LMO_prime_count(n);
}

UV prime_count_approx(UV n)
{
  if (n < 3000000) return segment_prime_count(2, n);
  return (UV) (RiemannR((long double) n, 1e-6) + 0.5);
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

  if (n <= 300070) { /* Quite accurate and avoids calling Li for speed. */
    /* Based on Axler 2022, page 9, Corollary 5.1 */
    a = (n <  69720) ? 905 :
        (n <  70120) ? 961 :
        (n <  88800) ? 918.2 :
        (n < 176000) ? 887.7 :
        (n < 299270) ? 839.46 :
                       846.66;  /* Good to 300071 */
    lower = fn / (fl1 - 1 - 1/fl1 - 2.975666/fl2 - 13.024334/(fl1*fl2) + a/(fl2*fl2));
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
   *
   * Axler 2022: https://arxiv.org/pdf/2203.05917.pdf (v4 2022) improves this.
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
  long double fn, flogn, flog2n, upper, c, d;

  if (n < NPRIMES_SMALL)
    return primes_small[n];

  fn     = (long double) n;
  flogn  = logl(n);
  flog2n = logl(flogn);    /* Note distinction between log_2(n) and log^2(n) */

  /* Binary search on prime count lower.  Good but quite slow. */
  if (n < 15360) {
    UV lo,hi;
    simple_nth_limits(&lo, &hi, fn, flogn, flog2n);
    while (lo < hi) {
      UV mid = lo + (hi-lo)/2;
      if (prime_count_lower(mid) < n) lo = mid+1;
      else                            hi = mid;
    }
    return lo;
  }

  /* See: Axler 2013, Dusart 2010 */
  /*      Axler 2017: http://arxiv.org/pdf/1706.03651.pdf */

  if      (n >= 46254381) { c = 2.00; d =  10.667;  } /* Axler 2017 Cor 1.2 */
  else if (n >=  8009824) { c = 2.00; d =  10.273;  } /* Axler 2013 Kor G */
  /* This is about 3x better than Dusart (2010) for 688382-8009823:
   *
   * else if (n >=   688382) { c = 2.30; d =   0.5730; }
   *
   * but we can split the range and get another 2x improvement in MSE.
   */
  else if (n >=  5450000) { c = 2.00; d =  10.1335; }  /*5450-8009 */
  else if (n >=  3906280) { c = 1.67; d =  20.2675; }  /*3906-5450 */
  else if (n >=  2110840) { c = 2.51; d =  -5.5714; }  /*2110-3906 */
  else if (n >=   876700) { c = 2.49; d =  -4.5129; }  /* 877-2110 */
  else if (n >=   688382) { c = 3.31; d = -26.3858; }  /* 688-877 */
  /* Use the Axler framework to get good bounds for smaller inputs. */
  else if (n >=   575750) { c =-0.79; d =  83.5215; }  /* 580-688 */
  else if (n >=   467650) { c = 0.93; d =  37.1597; }  /* 467-580 */
  else if (n >=   382440) { c = 2.92; d = -15.4768; }  /* 382-467 */
  else if (n >=   301130) { c = 5.92; d = -91.3415; }  /* 301-382 */
  else if (n >=   138630) { c = 2.01; d =   7.2842; }  /* 138-301 */
  else if (n >=    85820) { c = 2.07; d =   5.2103; }  /*  86-138 */
  else if (n >=    39016) { c = 2.76; d = -11.5918; }  /*  39- 86 */
  else if (n >=    31490) { c = 1.49; d =  15.1821; }  /*  31- 39 */
  else if (n >=    25070) { c =11.89; d =-197.8951; }  /*  25- 31 */
  else if (n >=    15359) { c = 4.80; d = -51.5928; }  /*  15- 25 */
  else                    { c = 3.92; d = -33.3994; }  /*   0- 15 */

  upper = fn * ( flogn + flog2n - 1.0 + ((flog2n-c)/flogn)
                 - (flog2n*flog2n-6*flog2n+d)/(2*flogn*flogn) );

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
    double b1 = (n < 56000000)  ?  11.200  :  11.50800000002;
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
    count = prime_count(lower_limit);

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
/*                                   MISC                                     */
/******************************************************************************/


double ramanujan_axler(long double n, long double c, long double d) {
  long double res, U, c1, c2, log2 = logl(2), logn = logl(n), loglogn = logl(logn);

  c1 = 2*log2*log2 + log2 + c;
  c2 = log2*log2*log2 + 2*log2*log2 + d;

  U = (log2 * logn*loglogn*loglogn - c1*logn*loglogn + c2*logn - log2*log2*loglogn + log2*log2*log2 + log2*log2)
    / (logn*logn*logn*logn + logn*logn*logn*loglogn - logn*logn*logn*log2 - logn*logn*log2);

  res = 2*n * (1.0L + log2/logn - (log2*loglogn - log2*log2 - log2) / (logn*logn) + U);
  return res;
}
