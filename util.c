#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Use long double to get a little more precision when we're calculating the
 * math functions -- especially those calculated with a series.  Long double
 * is defined in C89 (ISO C), so it should be supported by any reasonable
 * compiler we're using (seriously is your C compiler 20+ years out of date?).
 * Noting that 'long double' on many platforms is no different than 'double'
 * so it may buy us nothing.  But it's worth trying.
 */
extern long double powl(long double, long double);
extern long double expl(long double);
extern long double logl(long double);
extern long double fabsl(long double);

/* However, standard math functions weren't defined on them until C99.  Same
 * with the macro INFINITY.  There are some reasonable platforms I've seen
 * that don't have these. */
#ifndef INFINITY
  #define INFINITY (DBL_MAX + DBL_MAX)
  #define powl(x, y)  (long double) pow( (double) (x), (double) (y) )
  #define expl(x)  (long double) exp( (double) (x) )
  #define logl(x)  (long double) log( (double) (x) )
  #define fabsl(x)  (long double) fabs( (double) (x) )
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
static int _is_trial_prime7(UV n)
{
  UV limit, i;
  limit = sqrt(n);
  i = 7;
  while (1) {   /* trial division, skipping multiples of 2/3/5 */
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 4;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 2;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 4;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 2;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 4;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 6;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 2;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 6;
  }
  return 2;
}

/* Does trial division or prob tests, assuming x not divisible by 2, 3, or 5 */
static int _is_prime7(UV n)
{
  UV limit, i;

  if (n > MPU_PROB_PRIME_BEST)
    return _XS_is_prob_prime(n);  /* We know this works for all 64-bit n */

  limit = sqrt(n);
  i = 7;
  while (1) {   /* trial division, skipping multiples of 2/3/5 */
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 4;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 2;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 4;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 2;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 4;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 6;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 2;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 6;
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
int _XS_is_prime(UV n)
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
    return (_XS_is_prob_prime(n) == 2);

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


UV _XS_next_prime(UV n)
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


UV _XS_prev_prime(UV n)
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
 * I have tweaked the bounds formulas for small (under 70_000M) numbers so they
 * are tighter.  These bounds are verified via trial.  The Dusart bounds
 * (1.8 and 2.51) are used for larger numbers since those are proven.
 *
 */

static const unsigned char prime_count_small[] =
  {0,0,1,2,2,3,3,4,4,4,4,5,5,6,6,6,6,7,7,8,8,8,8,9,9,9,9,9,9,10,10,
   11,11,11,11,11,11,12,12,12,12,13,13,14,14,14,14,15,15,15,15,15,15,
   16,16,16,16,16,16,17,17,18,18,18,18,18,18,19};
#define NPRIME_COUNT_SMALL  (sizeof(prime_count_small)/sizeof(prime_count_small[0]))

UV _XS_prime_count(UV low, UV high)
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

  if (low == 7) {
    if      (high > (1UL << 42)) { count += 156661034233-3; low = 1UL<<42; }
    else if (high > (1UL << 39)) { count +=  21151907950-3; low = 1UL<<39; }
    else if (high > (1UL << 36)) { count +=   2874398515-3; low = 1UL<<36; }
    else if (high > (1UL << 33)) { count +=    393615806-3; low = 1UL<<33; }
    else if (high > (1UL << 30)) { count +=     54400028-3; low = 1UL<<30; }
    else if (high > (1UL << 27)) { count +=      7603553-3; low = 1UL<<27; }
    else if (high > (1UL << 24)) { count +=      1077871-3; low = 1UL<<24; }
    else if (high > (1UL << 20)) { count +=        82025-3; low = 1UL<<20; }
    else if (high > (1UL << 18)) { count +=        23000-3; low = 1UL<<18; }
    else if (high > (1UL << 16)) { count +=         6542-3; low = 1UL<<16; }
    else if (high > (1UL << 14)) { count +=         1900-3; low = 1UL<<14; }
  }

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

    if (high_d < segment_size) {
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

/* Note: We're keeping this here because we use it for nth_prime */
/* The nth prime will be less or equal to this number */
static UV _XS_nth_prime_upper(UV n)
{
  double fn = (double) n;
  double flogn, flog2n, upper;

  if (n < NPRIMES_SMALL)
    return primes_small[n];

  flogn  = log(n);
  flog2n = log(flogn);    /* Note distinction between log_2(n) and log^2(n) */

  if      (n >= 688383)    /* Dusart 2010 page 2 */
    upper = fn * (flogn + flog2n - 1.0 + ((flog2n-2.00)/flogn));
  else if (n >= 178974)    /* Dusart 2010 page 7 */
    upper = fn * (flogn + flog2n - 1.0 + ((flog2n-1.95)/flogn));
  else if (n >=  39017)    /* Dusart 1999 page 14 */
    upper = fn * (flogn + flog2n - 0.9484);
  else if (n >=     6)     /* Modified from Robin 1983 for 6-39016 _only_ */
    upper = fn * ( flogn  +  0.6000 * flog2n );
  else
    upper = fn * ( flogn + flog2n );

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


UV _XS_nth_prime(UV n)
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
  upper_limit = _XS_nth_prime_upper(n);
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

static long double const euler_mascheroni = 0.57721566490153286060651209008240243104215933593992L;
static long double const li2 = 1.045163780117492784844588889194613136522615578151L;

#define KAHAN_INIT(s) \
  long double s ## _y, s ## _t; \
  long double s ## _c = 0.0; \
  long double s = 0.0;

#define KAHAN_SUM(s, term) \
  s ## _y = term - s ## _c; \
  s ## _t = s + s ## _y; \
  s ## _c = (s ## _t - s) - s ## _y; \
  s = s ## _t;

double _XS_ExponentialIntegral(double x) {
  long double const tol = 1e-16;
  long double val, term;
  unsigned int n;
  KAHAN_INIT(sum);

  if (x == 0) croak("Invalid input to ExponentialIntegral:  x must be != 0");

  if (x < -1) {
    /* Continued fraction, good for x < -1 */
    long double lc = 0;
    long double ld = 1.0L / (1.0L - (long double)x);
    long double old, t, n2;
    val = ld * (-expl(x));
    for (n = 1; n <= 100000; n++) {
      t = (long double)(2*n + 1) - (long double) x;
      n2 = n * n;
      lc = 1.0L / (t - n2 * lc);
      ld = 1.0L / (t - n2 * ld);
      old = val;
      val *= ld/lc;
      if ( fabsl(val-old) <= tol*fabsl(val) )
        break;
    }
  } else if (x < 0) {
    /* Rational Chebyshev approximation (Cody, Thacher), good for -1 < x < 0 */
    static const long double C6p[7] = { -148151.02102575750838086L,
                                    150260.59476436982420737L,
                                     89904.972007457256553251L,
                                     15924.175980637303639884L,
                                      2150.0672908092918123209L,
                                       116.69552669734461083368L,
                                         5.0196785185439843791020L };
    static const long double C6q[7] = {  256664.93484897117319268L,
                                    184340.70063353677359298L,
                                     52440.529172056355429883L,
                                      8125.8035174768735759866L,
                                       750.43163907103936624165L,
                                        40.205465640027706061433L,
                                         1.0000000000000000000000L };
    long double sumn = C6p[0]-x*(C6p[1]-x*(C6p[2]-x*(C6p[3]-x*(C6p[4]-x*(C6p[5]-x*C6p[6])))));
    long double sumd = C6q[0]-x*(C6q[1]-x*(C6q[2]-x*(C6q[3]-x*(C6q[4]-x*(C6q[5]-x*C6q[6])))));
    val = logl(-x) - sumn/sumd;
  } else if (x < -logl(tol)) {
    /* Convergent series */
    long double fact_n = x;
    for (n = 2; n <= 200; n++) {
      long double invn = 1.0L / n;
      fact_n *= (long double)x * invn;
      term = fact_n * invn;
      KAHAN_SUM(sum, term);
      /* printf("C  after adding %.8lf, val = %.8lf\n", term, sum); */
      if ( term < tol*sum) break;
    }
    KAHAN_SUM(sum, euler_mascheroni);
    KAHAN_SUM(sum, logl(x));
    KAHAN_SUM(sum, x);
    val = sum;
  } else {
    /* Asymptotic divergent series */
    long double invx = 1.0L / x;
    term = 1.0;
    for (n = 1; n <= 200; n++) {
      long double last_term = term;
      term = term * ( (long double)n * invx );
      if (term < tol*sum) break;
      if (term < last_term) {
        KAHAN_SUM(sum, term);
        /* printf("A  after adding %.8lf, sum = %.8lf\n", term, sum); */
      } else {
        KAHAN_SUM(sum, (-last_term/3) );
        /* printf("A  after adding %.8lf, sum = %.8lf\n", -last_term/3, sum); */
        break;
      }
    }
    term = expl(x) * invx;
    val = term * sum + term;
  }

  return val;
}

double _XS_LogarithmicIntegral(double x) {
  if (x == 0) return 0;
  if (x == 1) return -INFINITY;
  if (x == 2) return li2;
  if (x <= 0) croak("Invalid input to LogarithmicIntegral:  x must be > 0");
  return _XS_ExponentialIntegral(log(x));
}

/*
 * Storing the first 10-20 Zeta values makes sense.  Past that it is purely
 * to avoid making the call to generate them ourselves.  We could cache the
 * calculated values. These all have 1 subtracted from them.  */
static const long double riemann_zeta_table[] = {
  0.6449340668482264364724151666460251892L,  /* zeta(2) */
  0.2020569031595942853997381615114499908L,
  0.0823232337111381915160036965411679028L,
  0.0369277551433699263313654864570341681L,
  0.0173430619844491397145179297909205279L,
  0.0083492773819228268397975498497967596L,
  0.0040773561979443393786852385086524653L,
  0.0020083928260822144178527692324120605L,
  0.0009945751278180853371459589003190170L,
  0.0004941886041194645587022825264699365L,
  0.0002460865533080482986379980477396710L,
  0.0001227133475784891467518365263573957L,
  0.0000612481350587048292585451051353337L,
  0.0000305882363070204935517285106450626L,
  0.0000152822594086518717325714876367220L,
  0.0000076371976378997622736002935630292L,  /* zeta(17)  Past here all we're */
  0.0000038172932649998398564616446219397L,  /* zeta(18)  getting is speed.   */
  0.0000019082127165539389256569577951013L,
  0.0000009539620338727961131520386834493L,
  0.0000004769329867878064631167196043730L,
  0.0000002384505027277329900036481867530L,
  0.0000001192199259653110730677887188823L,
  0.0000000596081890512594796124402079358L,
  0.0000000298035035146522801860637050694L,
  0.0000000149015548283650412346585066307L,
  0.0000000074507117898354294919810041706L,
  0.0000000037253340247884570548192040184L,
  0.0000000018626597235130490064039099454L,
  0.0000000009313274324196681828717647350L,
  0.0000000004656629065033784072989233251L,
  0.0000000002328311833676505492001455976L,
  0.0000000001164155017270051977592973835L,
  0.0000000000582077208790270088924368599L,
  0.0000000000291038504449709968692942523L,
  0.0000000000145519218910419842359296322L,
  0.0000000000072759598350574810145208690L,
  0.0000000000036379795473786511902372363L,
  0.0000000000018189896503070659475848321L,
  0.0000000000009094947840263889282533118L,
};
#define NPRECALC_ZETA (sizeof(riemann_zeta_table)/sizeof(riemann_zeta_table[0]))

/* Riemann Zeta on the real line.  Compare to Math::Cephes::zetac */
long double ld_riemann_zeta(long double x) {
  long double const tol = 1e-17;
  long double term;
  int k;
  KAHAN_INIT(sum);

  if (x < 0.5) croak("Invalid input to RiemannZeta:  x must be >= 0.5");

  if (x == (unsigned int)x) {
    k = x - 2;
    if ((k >= 0) && (k < NPRECALC_ZETA))
      return riemann_zeta_table[k];
  }

  /* Use Cody et al. rational Chebyshev approx for small values.  This is the
   * range where the series methods take a long time and are inaccurate.  This
   * method is fast and quite accurate over the range 0.5 - 5. */
  if (x <= 5.0) {

    static const long double C8p[9] = { 1.287168121482446392809e10L,
                                   1.375396932037025111825e10L,
                                   5.106655918364406103683e09L,
                                   8.561471002433314862469e08L,
                                   7.483618124380232984824e07L,
                                   4.860106585461882511535e06L,
                                   2.739574990221406087728e05L,
                                   4.631710843183427123061e03L,
                                   5.787581004096660659109e01L };
    static const long double C8q[9] = { 2.574336242964846244667e10L,
                                   5.938165648679590160003e09L,
                                   9.006330373261233439089e08L,
                                   8.042536634283289888587e07L,
                                   5.609711759541920062814e06L,
                                   2.247431202899137523543e05L,
                                   7.574578909341537560115e03L,
                                  -2.373835781373772623086e01L,
                                   1.000000000000000000000L    };
    long double sumn = C8p[0]+x*(C8p[1]+x*(C8p[2]+x*(C8p[3]+x*(C8p[4]+x*(C8p[5]+x*(C8p[6]+x*(C8p[7]+x*C8p[8])))))));
    long double sumd = C8q[0]+x*(C8q[1]+x*(C8q[2]+x*(C8q[3]+x*(C8q[4]+x*(C8q[5]+x*(C8q[6]+x*(C8q[7]+x*C8q[8])))))));
    sum = (sumn - (x-1)*sumd) / ((x-1)*sumd);

  } else {

    /* This series needs about half the number of terms as the usual k^-x,
     * and gets slightly better numerical results.
     *   functions.wolfram.com/ZetaFunctionsandPolylogarithms/Zeta/06/04/0003
     */
    for (k = 2; k <= 1000000; k++) {
      term = powl(2*k+1, -x);
      KAHAN_SUM(sum, term);
      if (term < tol*sum) break;
    }
    KAHAN_SUM(sum, powl(3, -x) );
    term = 1.0L / (1.0L - powl(2, -x));
    sum *= term;
    sum += (term - 1.0L);
    /* printf("zeta(%lf) = %.15lf in %d iterations\n", x, sum, k); */

  }

  return sum;
}

double _XS_RiemannR(double x) {
  long double const tol = 1e-16;
  long double part_term, term, flogx;
  unsigned int k;
  KAHAN_INIT(sum);

  if (x <= 0) croak("Invalid input to ReimannR:  x must be > 0");

  KAHAN_SUM(sum, 1.0);

  flogx = logl(x);
  part_term = 1;

  for (k = 1; k <= 10000; k++) {
    part_term *= flogx / k;
    term = part_term / (k + k * ld_riemann_zeta(k+1));
    KAHAN_SUM(sum, term);
    /* printf("R  after adding %.15lg, sum = %.15lg\n", term, sum); */
    if (fabsl(term/sum) < tol) break;
  }

  return sum;
}
