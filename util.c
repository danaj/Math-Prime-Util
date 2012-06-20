#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* It took until C99 to get this macro */
#ifndef INFINITY
#define INFINITY (DBL_MAX + DBL_MAX)
#endif

#include "ptypes.h"
#include "util.h"
#include "sieve.h"
#include "factor.h"
#include "cache.h"

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
  while (nbytes--)
    count += byte_zeros[*m++];
  return count;
}


/* Does trial division, assuming x not divisible by 2, 3, or 5 */
static int _is_trial_prime7(UV x)
{
  UV q, i;
  i = 7;
  while (1) {   /* trial division, skipping multiples of 2/3/5 */
    q = x/i;  if (q<i) return 1;  if (x==(q*i)) return 0;   i += 4;
    q = x/i;  if (q<i) return 1;  if (x==(q*i)) return 0;   i += 2;
    q = x/i;  if (q<i) return 1;  if (x==(q*i)) return 0;   i += 4;
    q = x/i;  if (q<i) return 1;  if (x==(q*i)) return 0;   i += 2;
    q = x/i;  if (q<i) return 1;  if (x==(q*i)) return 0;   i += 4;
    q = x/i;  if (q<i) return 1;  if (x==(q*i)) return 0;   i += 6;
    q = x/i;  if (q<i) return 1;  if (x==(q*i)) return 0;   i += 2;
    q = x/i;  if (q<i) return 1;  if (x==(q*i)) return 0;   i += 6;
  }
  return 2;
}

/* Does trial division or prob tests, assuming x not divisible by 2, 3, or 5 */
static int _is_prime7(UV x)
{
  UV q, i;

  if (x > MPU_PROB_PRIME_BEST)
    return is_prob_prime(x);  /* We know this works for all 64-bit n */

  i = 7;
  while (1) {   /* trial division, skipping multiples of 2/3/5 */
    q = x/i;  if (q<i) return 1;  if (x==(q*i)) return 0;   i += 4;
    q = x/i;  if (q<i) return 1;  if (x==(q*i)) return 0;   i += 2;
    q = x/i;  if (q<i) return 1;  if (x==(q*i)) return 0;   i += 4;
    q = x/i;  if (q<i) return 1;  if (x==(q*i)) return 0;   i += 2;
    q = x/i;  if (q<i) return 1;  if (x==(q*i)) return 0;   i += 4;
    q = x/i;  if (q<i) return 1;  if (x==(q*i)) return 0;   i += 6;
    q = x/i;  if (q<i) return 1;  if (x==(q*i)) return 0;   i += 2;
    q = x/i;  if (q<i) return 1;  if (x==(q*i)) return 0;   i += 6;
  }
  return 2;
}


/* Marked bits for each n, indicating if the number is prime */
static const unsigned char prime_is_small[] =
  {0xac,0x28,0x8a,0xa0,0x20,0x8a,0x20,0x28,0x88,0x82,0x08,0x02,0xa2,0x28,0x02,
   0x80,0x08,0x0a,0xa0,0x20,0x88,0x20,0x28,0x80,0xa2,0x00,0x08,0x80,0x28,0x82,
   0x02,0x08,0x82,0xa0,0x20,0x0a,0x20,0x00,0x88,0x22,0x00,0x08,0x02,0x28,0x82,
   0x80,0x20,0x88,0x20,0x20,0x02,0x02,0x28,0x80,0x82,0x08,0x02,0xa2,0x08,0x80,
   0x80,0x08,0x88,0x20,0x00,0x0a,0x00,0x20,0x08,0x20,0x08,0x0a,0x02,0x08,0x82,
   0x82,0x20,0x0a,0x80,0x00,0x8a,0x20,0x28,0x00,0x22,0x08,0x08,0x20,0x20,0x80,
   0x80,0x20,0x88,0x80,0x20,0x02,0x22,0x00,0x08,0x20,0x00,0x0a,0xa0,0x28,0x80,
   0x00,0x20,0x8a,0x00,0x20,0x8a,0x00,0x00,0x88,0x80,0x00,0x02,0x22,0x08,0x02};
#define NPRIME_IS_SMALL (sizeof(prime_is_small)/sizeof(prime_is_small[0]))

/* Return of 2 if n is prime, 0 if not.  Do it fast. */
int is_prime(UV n)
{
  UV d, m;
  unsigned char mtab;
  const unsigned char* sieve;
  int isprime;

  if ( n < (NPRIME_IS_SMALL*8))
    return ((prime_is_small[n/8] >> (n%8)) & 1) ? 2 : 0;

  d = n/30;
  m = n - d*30;
  mtab = masktab30[m];  /* Bitmask in mod30 wheel */

  /* Return 0 if a multiple of 2, 3, or 5 */
  if (mtab == 0)
    return 0;

  isprime = (n <= get_prime_cache(0, &sieve))
            ?  2*((sieve[d] & mtab) == 0)
            :  -1;
  release_prime_cache(sieve);

  return (isprime >= 0)  ?  isprime  :  _is_prime7(n);
}

/* Shortcut, asking for a very quick response of 1 = prime, 0 = dunno.
 * No trial divisions will be done, making this useful for factoring.
 */
int is_definitely_prime(UV n)
{
  UV d, m;
  unsigned char mtab;
  const unsigned char* sieve;
  int isprime;

  if ( n < (NPRIME_IS_SMALL*8))
    return ((prime_is_small[n/8] >> (n%8)) & 1);

  d = n/30;
  m = n - d*30;
  mtab = masktab30[m];  /* Bitmask in mod30 wheel */

  /* Return 0 if a multiple of 2, 3, or 5 */
  if (mtab == 0)
    return 0;

  isprime = (n <= get_prime_cache(0, &sieve))
            ?  ((sieve[d] & mtab) == 0)
            :  -1;
  release_prime_cache(sieve);
  if (isprime >= 0)  return isprime;

  if (n > MPU_PROB_PRIME_BEST)
    return (is_prob_prime(n) == 2);

  return 0;
}


static const unsigned char prime_next_small[] =
  {2,2,3,5,5,7,7,11,11,11,11,13,13,17,17,17,17,19,19,23,23,23,23,
   29,29,29,29,29,29,31,31,37,37,37,37,37,37,41,41,41,41,43,43,47,
   47,47,47,53,53,53,53,53,53,59,59,59,59,59,59,61,61,67,67,67,67,67,67,71};
#define NPRIME_NEXT_SMALL (sizeof(prime_next_small)/sizeof(prime_next_small[0]))

UV next_trial_prime(UV n)
{
  UV d,m;

  if (n < NPRIME_NEXT_SMALL)
    return prime_next_small[n];

  d = n/30;
  m = n - d*30;
  /* Move forward one, knowing we may not be on the wheel */
  if (m == 29) { d++; m = 1; } else  { m = nextwheel30[m]; }
  while (!_is_trial_prime7(d*30+m)) {
    m = nextwheel30[m];  if (m == 1) d++;
  }
  return(d*30+m);
}


UV next_prime(UV n)
{
  UV d, m;
  const unsigned char* sieve;
  UV sieve_size;

  if (n < NPRIME_NEXT_SMALL)
    return prime_next_small[n];

  /* Overflow */
#if BITS_PER_WORD == 32
  if (n >= UVCONST(4294967291))  return 0;
#else
  if (n >= UVCONST(18446744073709551557))  return 0;
#endif

  sieve_size = get_prime_cache(0, &sieve);
  if (n < sieve_size) {
    START_DO_FOR_EACH_SIEVE_PRIME(sieve, n+1, sieve_size)
      { release_prime_cache(sieve); return p; }
    END_DO_FOR_EACH_SIEVE_PRIME;
    /* Not found, so must be larger than the cache size */
    n = sieve_size;
  }
  release_prime_cache(sieve);

  d = n/30;
  m = n - d*30;
  /* Move forward one, knowing we may not be on the wheel */
  if (m == 29) { d++; m = 1; } else  { m = nextwheel30[m]; }
  while (!_is_prime7(d*30+m)) {
    /* Move forward one, knowing we are on the wheel */
    m = nextwheel30[m];  if (m == 1) d++;
  }
  return(d*30+m);
}


UV prev_prime(UV n)
{
  UV d, m;
  const unsigned char* sieve;
  UV sieve_size;

  /* TODO: small prev prime */
  if (n <= 7)
    return (n <= 2) ? 0 : (n <= 3) ? 2 : (n <= 5) ? 3 : 5;

  d = n/30;
  m = n - d*30;

  sieve_size = get_prime_cache(0, &sieve);
  if (n < sieve_size) {
    do {
      m = prevwheel30[m];
      if (m==29) { MPUassert(d>0, "d 0 in prev_prime");  d--; }
    } while (sieve[d] & masktab30[m]);
    release_prime_cache(sieve);
  } else {
    release_prime_cache(sieve);
    do {
      m = prevwheel30[m];
      if (m==29) { MPUassert(d>0, "d 0 in prev_prime");  d--; }
    } while (!_is_prime7(d*30+m));
  }
  return(d*30+m);
}



/* Given a sieve of size nbytes, walk it counting zeros (primes) until:
 *
 * (1) we counted them all: return the count, which will be less than maxcount.
 *
 * (2) we hit maxcount: set position to the index of the maxcount'th prime
 *     and return count (which will be equal to maxcount).
 */
static UV count_segment_maxcount(const unsigned char* sieve, UV nbytes, UV maxcount, UV* pos)
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

  /* Count until we reach the end or >= maxcount */
  while ( (sieveptr < maxsieve) && (count < maxcount) )
    count += byte_zeros[*sieveptr++];
  /* If we went one too far, back up.  Count will always be < maxcount */
  if (count >= maxcount)
    count -= byte_zeros[*--sieveptr];
  /* We counted this many bytes */
  byte = sieveptr - sieve;

  MPUassert(count < maxcount, "count_segment_maxcount wrong count");

  if (byte == nbytes)
    return count;

  /* The result is somewhere in the next byte */
  START_DO_FOR_EACH_SIEVE_PRIME(sieve, byte*30+1, nbytes*30-1)
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
  UV count = 0;

  UV lo_d = lowp/30;
  UV lo_m = lowp - lo_d*30;
  UV hi_d = highp/30;
  UV hi_m = highp - hi_d*30;

  MPUassert( sieve != 0, "count_segment_ranged incorrect args");

  if (hi_d >= nbytes) {
    hi_d = nbytes-1;
    highp = hi_d*30+29;
  }

  if ( (nbytes == 0) || (highp < lowp) )
    return 0;

#if 0
  /* Dead simple way */
  START_DO_FOR_EACH_SIEVE_PRIME(sieve, lowp, highp)
    count++;
  END_DO_FOR_EACH_SIEVE_PRIME;
  return count;
#endif

  /* Count first fragment */
  if (lo_m > 1) {
    UV upper = (highp <= (lo_d*30+29)) ? highp : (lo_d*30+29);
    START_DO_FOR_EACH_SIEVE_PRIME(sieve, lowp, upper)
      count++;
    END_DO_FOR_EACH_SIEVE_PRIME;
    lowp = upper+2;
    lo_d = lowp/30;
    lo_m = lowp - lo_d*30;
  }
  if (highp < lowp)
    return count;

  /* Count bytes in the middle */
  {
    UV count_bytes = hi_d - lo_d + (hi_m == 29);
    if (count_bytes > 0) {
      count += count_zero_bits(sieve+lo_d, count_bytes);
      lowp += 30*count_bytes;
      lo_d = lowp/30;
      lo_m = lowp - lo_d*30;
    }
  }
  if (highp < lowp)
    return count;

  /* Count last fragment */
  START_DO_FOR_EACH_SIEVE_PRIME(sieve, lowp, highp)
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
 * I have tweaked the bounds formulas for small (under 4000M) numbers so they
 * are tighter.  These bounds are verified via trial.  The Dusart bounds
 * (1.8 and 2.51) are used for larger numbers since those are proven.
 *
 */

static const unsigned char prime_count_small[] =
  {0,0,1,2,2,3,3,4,4,4,4,5,5,6,6,6,6,7,7,8,8,8,8,9,9,9,9,9,9,10,10,
   11,11,11,11,11,11,12,12,12,12,13,13,14,14,14,14,15,15,15,15,15,15,
   16,16,16,16,16,16,17,17,18,18,18,18,18,18,19};
#define NPRIME_COUNT_SMALL  (sizeof(prime_count_small)/sizeof(prime_count_small[0]))

static const double F1 = 1.0;
UV prime_count_lower(UV x)
{
  double fx, flogx;
  double a = 1.80;     /* Dusart 1999, page 14 */

  if (x < NPRIME_COUNT_SMALL)
    return prime_count_small[x];

  fx = (double)x;
  flogx = log(x);

  if (x < 599)
    return (UV) (fx / (flogx-0.7));

  if      (x <     2700)  { a = 0.30; }
  else if (x <     5500)  { a = 0.90; }
  else if (x <    19400)  { a = 1.30; }
  else if (x <    32299)  { a = 1.60; }
  else if (x <   176000)  { a = 1.80; }
  else if (x <   315000)  { a = 2.10; }
  else if (x <  1100000)  { a = 2.20; }
  else if (x <  4500000)  { a = 2.31; }
  else if (x <233000000)  { a = 2.36; }
#if BITS_PER_WORD == 32
  else a = 2.32;
#else
  else if (x < UVCONST( 5433800000)) { a = 2.32; }
  else if (x < UVCONST(60000000000)) { a = 2.15; }
#endif

  return (UV) ( (fx/flogx) * (F1 + F1/flogx + a/(flogx*flogx)) );
}


UV prime_count_upper(UV x)
{
  double fx, flogx;
  double a = 2.51;    /* Dusart 1999, page 14 */

  if (x < NPRIME_COUNT_SMALL)
    return prime_count_small[x];

  fx = (double)x;
  flogx = log(x);

  /* This function is unduly complicated. */

  if (x < 1621)  return (UV) (fx / (flogx-1.048) + F1);
  if (x < 5000)  return (UV) (fx / (flogx-1.071) + F1);
  if (x < 15900) return (UV) (fx / (flogx-1.098) + F1);

  if      (x <    24000) {  a = 2.30; }
  else if (x <    59000) {  a = 2.48; }
  else if (x <   350000) {  a = 2.52; }
  else if (x <   355991) {  a = 2.54; }
  else if (x <   356000) {  a = 2.51; }
  else if (x <  3550000) {  a = 2.50; }
  else if (x <  3560000) {  a = 2.49; }
  else if (x <  5000000) {  a = 2.48; }
  else if (x <  8000000) {  a = 2.47; }
  else if (x < 13000000) {  a = 2.46; }
  else if (x < 18000000) {  a = 2.45; }
  else if (x < 31000000) {  a = 2.44; }
  else if (x < 41000000) {  a = 2.43; }
  else if (x < 48000000) {  a = 2.42; }
  else if (x <119000000) {  a = 2.41; }
  else if (x <182000000) {  a = 2.40; }
  else if (x <192000000) {  a = 2.395; }
  else if (x <213000000) {  a = 2.390; }
  else if (x <271000000) {  a = 2.385; }
  else if (x <322000000) {  a = 2.380; }
  else if (x <400000000) {  a = 2.375; }
  else if (x <510000000) {  a = 2.370; }
  else if (x <682000000) {  a = 2.367; }
#if BITS_PER_WORD == 32
  else a = 2.362;
#else
  else if (x < UVCONST(60000000000)) { a = 2.362; }
#endif

  /*
   * An alternate idea:
   *  float alog[23] = {  2.30,2.30,2.30,2.30,2.30,2.30,2.30 ,2.30,2.30,2.30,
   *                      2.47,2.49,2.53,2.50,2.49,2.49,2.456,2.44,2.40,2.370,
   *                      2.362,2.362,2.362,2.362};
   *  float clog[23] = {  0,   0,   0,   0,   0,   0,   0,    0,   0,   1,
   *                      3,   1,   2,   1,   3,   2,   5,   -6,   1,   1,
   *                      1,   1,   1,   1};
   *  if ((int)flogx < 23) {
   *    a = alog[(int)flogx];
   *    return ((UV) ( (fx/flogx) * (F1 + F1/flogx + a/(flogx*flogx)) ) + clog[(int)flogx] + 0.01);
   *  }
   *
   * Another thought is to use more terms in the Li(x) expansion along with
   * a subtraction [Li(x) > Pi(x) for x < 10^316 or so, so for our 64-bit
   * version we should be fine].
   */

  return (UV) ( (fx/flogx) * (F1 + F1/flogx + a/(flogx*flogx)) + F1 );
}


UV prime_count_approx(UV x)
{
  /*
   * A simple way:
   *     return ((prime_count_lower(x) + prime_count_upper(x)) / 2);
   * With the current bounds, this is ~131k at 10^10 and 436B at 10^19.
   *
   * The logarithmic integral works quite well, with absolute errors of
   * ~3100 at 10^10 and ~100M at 10^19.
   *
   * Riemann's R function works astoundingly well, with errors of ~1828
   * at 10^10 and 24M at 10^19.
   *
   * Getting fancier, one try using Riemann's pi formula:
   *     http://trac.sagemath.org/sage_trac/ticket/8135
   */
  double R;
  if (x < NPRIME_COUNT_SMALL)
    return prime_count_small[x];

  R = RiemannR(x);
  /* We could add the additional factor:
   *   R = R - (1.0 / log(x)) + (M_1_PI * atan(M_PI/log(x)))
   * but it's extraordinarily small, so not worth calculating here.
   */
  return (UV)(R+0.5);
}


UV prime_count(UV low, UV high)
{
  const unsigned char* cache_sieve;
  unsigned char* segment;
  UV segment_size, low_d, high_d;
  UV count = 0;

  if ( (low <= 2) && (high < NPRIME_COUNT_SMALL) )
    return prime_count_small[high];

  if ((low <= 2) && (high >= 2)) count++;
  if ((low <= 3) && (high >= 3)) count++;
  if ((low <= 5) && (high >= 5)) count++;
  if (low < 7)  low = 7;

  if (low > high)  return count;

  low_d = low/30;
  high_d = high/30;

  /* Count full bytes only -- no fragments from primary cache */
  segment_size = get_prime_cache(0, &cache_sieve) / 30;
  if (segment_size < high_d) {
    /* Expand sieve to sqrt(n) */
    UV endp = (high_d >= (UV_MAX/30))  ?  UV_MAX-2  :  30*high_d+29;
    release_prime_cache(cache_sieve);
    segment_size = get_prime_cache( sqrt(endp) + 1 , &cache_sieve) / 30;
  }

  if ( (segment_size > 0) && (low_d <= segment_size) ) {
    /* Count all the primes in the primary cache in our range */
    count += count_segment_ranged(cache_sieve, segment_size, low, high);

    if (high_d <= segment_size) {
      release_prime_cache(cache_sieve);
      return count;
    }

    low_d = segment_size;
  }
  release_prime_cache(cache_sieve);

  /* More primes needed.  Repeatedly segment sieve. */
  segment = get_prime_segment(&segment_size);
  if (segment == 0)
    croak("Could not get segment memory");

  while (low_d <= high_d)
  {
    UV seghigh_d = ((high_d - low_d) < segment_size)
                   ? high_d
                   : (low_d + segment_size-1);
    UV range_d = seghigh_d - low_d + 1;
    UV seglow  = (low_d*30 < low) ? low : low_d*30;
    UV seghigh = (seghigh_d == high_d) ? high : (seghigh_d*30+29);

    if (sieve_segment(segment, low_d, seghigh_d) == 0) {
      release_prime_segment(segment);
      croak("Could not segment sieve from %"UVuf" to %"UVuf, low_d*30+1, 30*seghigh_d+29);
    }

    count += count_segment_ranged(segment, segment_size, seglow - low_d*30, seghigh - low_d*30);

    low_d += range_d;
  }
  release_prime_segment(segment);

  return count;
}



static const unsigned short primes_small[] =
  {0,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
   101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,
   193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,
   293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,
   409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499};
#define NPRIMES_SMALL (sizeof(primes_small)/sizeof(primes_small[0]))

/* The nth prime will be greater than or equal to this number */
UV nth_prime_lower(UV n)
{
  double fn = (double) n;
  double flogn, flog2n, lower;

  if (n < NPRIMES_SMALL)
    return (n==0) ? 0 : primes_small[n];

  flogn  = log(n);
  flog2n = log(flogn);    /* Note distinction between log_2(n) and log^2(n) */

  /* Dusart 1999 page 14, for all n >= 2 */
  lower = fn * (flogn + flog2n - 1.0 + ((flog2n-2.25)/flogn));

  /* Watch out for overflow */
  if (lower >= (double)UV_MAX) {
#if BITS_PER_WORD == 32
    if (n <= UVCONST(203280221)) return UVCONST(4294967291);
#else
    if (n <= UVCONST(425656284035217743)) return UVCONST(18446744073709551557);
#endif
    croak("nth_prime_lower(%"UVuf") overflow", n);
  }

  return (UV) lower;
}


/* The nth prime will be less or equal to this number */
UV nth_prime_upper(UV n)
{
  double fn = (double) n;
  double flogn, flog2n, upper;

  if (n < NPRIMES_SMALL)
    return primes_small[n];

  flogn  = log(n);
  flog2n = log(flogn);    /* Note distinction between log_2(n) and log^2(n) */

  if (n >= 39017)
    upper = fn * ( flogn  +  flog2n - 0.9484 ); /* Dusart 1999 page 14*/
  else if (n >= 27076)
    upper = fn * (flogn + flog2n - 1.0 + ((flog2n-1.80)/flogn)); /*Dusart 1999*/
  else if (n >= 7022)
    upper = fn * ( flogn  +  0.9385 * flog2n ); /* Robin 1983 */
  else
    upper = fn * ( flogn + flog2n );

  /* Watch out for  overflow */
  if (upper >= (double)UV_MAX) {
#if BITS_PER_WORD == 32
    if (n <= UVCONST(203280221)) return UVCONST(4294967291);
#else
    if (n <= UVCONST(425656284035217743)) return UVCONST(18446744073709551557);
#endif
    croak("nth_prime_upper(%"UVuf") overflow", n);
  }

  return (UV) ceil(upper);
}


UV nth_prime_approx(UV n)
{
  double fn, flogn, flog2n, order, approx;

  if (n < NPRIMES_SMALL)
    return primes_small[n];

  /* This isn't too bad:
   *    return ((nth_prime_lower(n) + nth_prime_upper(n)) / 2);
   */

  fn = (double) n;
  flogn = log(n);
  flog2n = log(flogn);

  /* Cipolla 1902:
   *    m=0   fn * ( flogn + flog2n - 1 );
   *    m=1   + ((flog2n - 2)/flogn) );
   *    m=2   - (((flog2n*flog2n) - 6*flog2n + 11) / (2*flogn*flogn))
   *    + O((flog2n/flogn)^3)
   *
   * Shown in Dusart 1999 page 12, as well as other sources such as:
   *   http://www.emis.de/journals/JIPAM/images/153_02_JIPAM/153_02.pdf
   * where the main issue you run into is that you're doing polynomial
   * interpolation, so it oscillates like crazy with many high-order terms.
   * Hence I'm leaving it at m=2.
   */
  approx = fn * ( flogn + flog2n - 1
                  + ((flog2n - 2)/flogn)
                  - (((flog2n*flog2n) - 6*flog2n + 11) / (2*flogn*flogn))
                );

  /* Apply a correction to help keep values close */
  order = flog2n/flogn;
  order = order*order*order * fn;

  if      (n <   259) approx += 10.4 * order;
  else if (n <   775) approx += 7.52 * order;
  else if (n <  1271) approx += 5.6 * order;
  else if (n <  2000) approx += 5.2 * order;
  else if (n <  4000) approx += 4.3 * order;
  else if (n < 12000) approx += 3.0 * order;
  else if (n <150000) approx += 2.1 * order;
  else if (n <200000000) approx += 0.0 * order;
  else                approx += -0.010 * order; /* -0.025 is better */

  /* For all three analytical functions, it is possible that for a given valid
   * input, we will not be able to return an output that fits in the UV type.
   * For example, if they ask for the 203280222nd prime, we should return
   * 4294967311.  But in 32-bit, that overflows.  What we do is calculate our
   * double precision value.  If that would overflow, then we look at the input
   * and if it is <= the index of the last representable prime, then we return
   * the last representable prime.  Otherwise, we croak an overflow message.
   * This should maintain the invariant:
   *    nth_prime_lower(n)  <=  nth_prime(n)  <=  nth_prime_upper(n)
   */
  if (approx >= (double)UV_MAX) {
#if BITS_PER_WORD == 32
    if (n <= UVCONST(203280221)) return UVCONST(4294967291);
#else
    if (n <= UVCONST(425656284035217743)) return UVCONST(18446744073709551557);
#endif
    croak("nth_prime_approx(%"UVuf") overflow", n);
  }

  return (UV) (approx + 0.5);
}


UV nth_prime(UV n)
{
  const unsigned char* cache_sieve;
  unsigned char* segment;
  UV upper_limit, segbase, segment_size;
  UV p = 0;
  UV target = n-3;
  UV count = 0;

  /* If very small, return the table entry */
  if (n < NPRIMES_SMALL)
    return primes_small[n];

  /* Determine a bound on the nth prime.  We know it comes before this. */
  upper_limit = nth_prime_upper(n);
  MPUassert(upper_limit > 0, "nth_prime got an upper limit of 0");

  /* Get the primary cache, and ensure it is at least this large.  If the
   * input is small enough, get a sieve covering the range.  Otherwise, we'll
   * walk segments.  Make sure we have enough primes in the cache so the
   * segmented siever won't have to keep resieving.
   */
  if (upper_limit <= (1*1024*1024*30))
    segment_size = get_prime_cache(upper_limit, &cache_sieve) / 30;
  else
    segment_size = get_prime_cache(sqrt(upper_limit), &cache_sieve) / 30;

  /* Count up everything in the cached sieve. */
  if (segment_size > 0)
    count += count_segment_maxcount(cache_sieve, segment_size, target, &p);
  release_prime_cache(cache_sieve);
  if (count == target)
    return p;

  /* Start segment sieving.  Get memory to sieve into. */
  segbase = segment_size;
  segment = get_prime_segment(&segment_size);
  if (segment == 0)
    croak("Could not get segment memory");

  while (count < target) {
    /* Limit the segment size if we know the answer comes earlier */
    if ( (30*(segbase+segment_size)+29) > upper_limit )
      segment_size = (upper_limit - segbase*30 + 30) / 30;

    /* Do the actual sieving in the range */
    if (sieve_segment(segment, segbase, segbase + segment_size-1) == 0) {
      release_prime_segment(segment);
      croak("Could not segment sieve from %"UVuf" to %"UVuf, 30*segbase+1, 30*(segbase+segment_size)+29);
    }

    /* Count up everything in this segment */
    count += count_segment_maxcount(segment, segment_size, target-count, &p);

    if (count < target)
      segbase += segment_size;
  }
  release_prime_segment(segment);
  MPUassert(count == target, "nth_prime got incorrect count");
  return ( (segbase*30) + p );
}




/*
 * See:
 *  "Multiple-Precision Exponential Integral and Related Functions"
 *      by David M. Smith
 *  "On the Evaluation of the Complex-Valued Exponential Integral"
 *      by Vincent Pegoraro and Philipp Slusallek
 *  "Numerical Recipes" 3rd edition
 *      by William H. Press et al.
 *  "Rational Chevyshev Approximations for the Exponential Integral E_1(x)"
 *      by W. J. Cody and Henry C. Thacher, Jr.
 *
 * Any mistakes here are completely my fault.  This code has not been
 * verified for anything serious.  For better results, see:
 *    http://www.trnicely.net/pi/pix_0000.htm
 * which although the author claims are demonstration programs, will
 * undoubtedly produce more reliable results than this code does (I don't
 * know of any obvious issues with this code, but it just hasn't been used
 * by many people).
 */

static double const euler_mascheroni = 0.57721566490153286060651209008240243104215933593992;
static double const li2 = 1.045163780117492784844588889194613136522615578151;

double ExponentialIntegral(double x) {
  double const tol = 1e-16;
  double val, term, fact_n;
  double y, t;
  double sum = 0.0;
  double c = 0.0;
  int n;

  if (x == 0) croak("Invalid input to ExponentialIntegral:  x must be != 0");

  if (x < -1) {
    /* Continued fraction, good for x < -1 */
    double old;
    double lc = 0;
    double ld = 1 / (1 - x);
    val = ld * (-exp(x));
    for (n = 1; n <= 100000; n++) {
      lc = 1 / (2*n + 1 - x - n*n*lc);
      ld = 1 / (2*n + 1 - x - n*n*ld);
      old = val;
      val *= ld/lc;
      if ( fabs(val-old) <= tol*fabs(val) )
        break;
    }
  } else if (x < 0) {
    /* Rational Chebyshev approximation (Cody, Thacher), good for -1 < x < 0 */
#if 0
    static const double C2p[3] = { -4.43668255, 4.42054938, 3.16274620 };
    static const double C2q[3] = {  7.68641124, 5.65655216, 1.00000000 };
    double sumn = C2p[0] - x*(C2p[1] - x*C2p[2]);
    double sumd = C2q[0] - x*(C2q[1] - x*C2q[2]);
#else
    static const double C6p[7] = { -148151.02102575750838086,
                                    150260.59476436982420737,
                                     89904.972007457256553251,
                                     15924.175980637303639884,
                                      2150.0672908092918123209,
                                       116.69552669734461083368,
                                         5.0196785185439843791020 };
    static const double C6q[7] = {  256664.93484897117319268,
                                    184340.70063353677359298,
                                     52440.529172056355429883,
                                      8125.8035174768735759866,
                                       750.43163907103936624165,
                                        40.205465640027706061433,
                                         1.0000000000000000000000 };
    double sumn = C6p[0]-x*(C6p[1]-x*(C6p[2]-x*(C6p[3]-x*(C6p[4]-x*(C6p[5]-x*C6p[6])))));
    double sumd = C6q[0]-x*(C6q[1]-x*(C6q[2]-x*(C6q[3]-x*(C6q[4]-x*(C6q[5]-x*C6q[6])))));
#endif
    val = log(-x) - sumn/sumd;
  } else if (x < -log(tol)) {
    /* Convergent series */
    fact_n = 1;

    y = euler_mascheroni-c;  t = sum+y;  c = (t-sum)-y;  sum = t;
    y = log(x)-c;  t = sum+y;  c = (t-sum)-y;  sum = t;

    for (n = 1; n <= 200; n++) {
      fact_n *= x/n;
      term = fact_n/n;
      y = term-c;  t = sum+y;  c = (t-sum)-y;  sum = t;
      /* printf("C  after adding %.8lf, val = %.8lf\n", term, sum); */
      if (term < tol) break;
    }
    val = sum;
  } else {
    /* Asymptotic divergent series */
    double last_term;

    val = exp(x) / x;
    term = 1.0;
    y = 1.0-c;  t = sum+y;  c = (t-sum)-y;  sum = t;

    for (n = 1; n <= 200; n++) {
      last_term = term;
      term *= n/x;
      if (term < tol) break;
      if (term < last_term) {
        y = term-c;  t = sum+y;  c = (t-sum)-y;  sum = t;
        /* printf("A  after adding %.8lf, sum = %.8lf\n", term, sum); */
      } else {
        y = (-last_term/3)-c;  t = sum+y;  c = (t-sum)-y;  sum = t;
        /* printf("A  after adding %.8lf, sum = %.8lf\n", -last_term/3, sum); */
        break;
      }
    }
    val *= sum;
  }

  return val;
}

double LogarithmicIntegral(double x) {
  if (x == 0) return 0;
  if (x == 1) return -INFINITY;
  if (x == 2) return li2;
  if (x <= 0) croak("Invalid input to ExponentialIntegral:  x must be > 0");
  return ExponentialIntegral(log(x));
}

/*
 * Storing the first 10-20 Zeta values makes sense.  Past that it is purely
 * to avoid making the call to generate them ourselves.  We could cache the
 * calculated values. These all have 1 subtracted from them.  */
static const double riemann_zeta_table[] = {
  0.6449340668482264364724151666460251892,  /* zeta(2) */
  0.2020569031595942853997381615114499908,
  0.0823232337111381915160036965411679028,
  0.0369277551433699263313654864570341681,
  0.0173430619844491397145179297909205279,
  0.0083492773819228268397975498497967596,
  0.0040773561979443393786852385086524653,
  0.0020083928260822144178527692324120605,
  0.0009945751278180853371459589003190170,
  0.0004941886041194645587022825264699365,
  0.0002460865533080482986379980477396710,
  0.0001227133475784891467518365263573957,
  0.0000612481350587048292585451051353337,
  0.0000305882363070204935517285106450626,
  0.0000152822594086518717325714876367220,
  0.0000076371976378997622736002935630292,  /* zeta(17)  Past here all we're */
  0.0000038172932649998398564616446219397,  /* zeta(18)  getting is speed.   */
  0.0000019082127165539389256569577951013,
  0.0000009539620338727961131520386834493,
  0.0000004769329867878064631167196043730,
  0.0000002384505027277329900036481867530,
  0.0000001192199259653110730677887188823,
  0.0000000596081890512594796124402079358,
  0.0000000298035035146522801860637050694,
  0.0000000149015548283650412346585066307,
  0.0000000074507117898354294919810041706,
  0.0000000037253340247884570548192040184,
  0.0000000018626597235130490064039099454,
  0.0000000009313274324196681828717647350,
  0.0000000004656629065033784072989233251,
  0.0000000002328311833676505492001455976,
  0.0000000001164155017270051977592973835,
  0.0000000000582077208790270088924368599,
  0.0000000000291038504449709968692942523,
  0.0000000000145519218910419842359296322,
  0.0000000000072759598350574810145208690,
  0.0000000000036379795473786511902372363,
  0.0000000000018189896503070659475848321,
  0.0000000000009094947840263889282533118,
};
#define NPRECALC_ZETA (sizeof(riemann_zeta_table)/sizeof(riemann_zeta_table[0]))

static double evaluate_zeta(double x) {
  double const tol = 1e-16;
  double y, t;
  double sum = 0.0;
  double c = 0.0;
  double term;
  int k;

  /* Simple method.  Slow and inaccurate near x=1, but iterations and accuracy
   * improve quickly.
   */

  /* term = 1.0;          y = term-c;  t = sum+y;  c = (t-sum)-y;  sum = t; */
  term = 1.0/exp2(x);  y = term-c;  t = sum+y;  c = (t-sum)-y;  sum = t;

  for (k = 3; k <= 100000; k++) {
    if (fabs(term) < tol) break;
    if      (k ==  4) term = 1.0 / exp2(2*x);
    else if (k ==  8) term = 1.0 / exp2(4*x);
    else if (k == 16) term = 1.0 / exp2(8*x);
    else              term = 1.0 / pow(k, x);
    y = term-c;  t = sum+y;  c = (t-sum)-y;  sum = t;
  }
  return sum;
}

double RiemannR(double x) {
  double const tol = 1e-16;
  double y, t, part_term, term, flogx, zeta;
  double sum = 0.0;
  double c = 0.0;
  int k;

  if (x <= 0) croak("Invalid input to ReimannR:  x must be > 0");

  y = 1.0-c;  t = sum+y;  c = (t-sum)-y;  sum = t;

  flogx = log(x);
  part_term = 1;

  /* Do small k with zetas from table */
  for (k = 1; k <= (int)NPRECALC_ZETA; k++) {
    zeta = riemann_zeta_table[k+1-2];
    part_term *= flogx / k;
    term = part_term / (k + k * zeta);
    y = term-c;  t = sum+y;  c = (t-sum)-y;  sum = t;
    if (fabs(term) < tol) break;
    /* printf("R  after adding %.8lf, sum = %.8lf\n", term, sum); */
  }
  /* Finish with function */
  for (; k <= 10000; k++) {
    if (fabs(term) < tol) break;
    zeta = evaluate_zeta(k+1);
    part_term *= flogx / k;
    term = part_term / (k + k * zeta);
    y = term-c;  t = sum+y;  c = (t-sum)-y;  sum = t;
    /* printf("R  after adding %.8lf, sum = %.8lf\n", term, sum); */
  }

  return sum;
}
