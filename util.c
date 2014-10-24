#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

/* Use long double to get a little more precision when we're calculating the
 * math functions -- especially those calculated with a series.  Long double
 * is defined in C89 (ISO C), so it should be supported by any reasonable
 * compiler we're using (seriously is your C compiler 20+ years out of date?).
 * Noting that 'long double' on many platforms is no different than 'double'
 * so it may buy us nothing.  But it's worth trying.
 */

/* The C99 LD math functions are a clusterfrack.  They're defined by C99, but
 * NetBSD doesn't have them.  You need them in both the headers and libraries,
 * but there is no standard way to find out if the libraries have them.  The
 * best way (I believe) to deal with this is having the make system do test
 * compiles.  Barring that, we make limited guesses, and just give up
 * precision on any system we don't recognize.
 */
#if _MSC_VER
  /* MSVS has these as macros, and really doesn't want us defining them. */
#elif defined(__MATH_DECLARE_LDOUBLE) || \
      defined(__LONG_DOUBLE_128__) || \
      defined(__LONGDOUBLE128)
#if defined(__IBMC__) || defined(__IBMCPP__)
  /* XLC is different  */
#else
  /* GLIBC */
  extern long double powl(long double, long double);
  extern long double expl(long double);
  extern long double logl(long double);
  extern long double fabsl(long double);
  extern long double floorl(long double);
  extern long double ceill(long double);
#endif
#else
  #define powl(x, y)  (long double) pow( (double) (x), (double) (y) )
  #define expl(x)     (long double) exp( (double) (x) )
  #define logl(x)     (long double) log( (double) (x) )
  #define fabsl(x)    (long double) fabs( (double) (x) )
  #define floorl(x)   (long double) floor( (double) (x) )
  #define ceill(x)    (long double) ceil( (double) (x) )
#endif

#ifdef LDBL_INFINITY
  #undef INFINITY
  #define INFINITY LDBL_INFINITY
#elif !defined(INFINITY)
  #define INFINITY (DBL_MAX + DBL_MAX)
#endif

#ifndef LDBL_EPSILON
  #define LDBL_EPSILON 1e-16
#endif
#ifndef LDBL_MAX
  #define LDBL_MAX DBL_MAX
#endif

#define KAHAN_INIT(s) \
  long double s ## _y, s ## _t; \
  long double s ## _c = 0.0; \
  long double s = 0.0;

#define KAHAN_SUM(s, term) \
  do { \
    s ## _y = (term) - s ## _c; \
    s ## _t = s + s ## _y; \
    s ## _c = (s ## _t - s) - s ## _y; \
    s = s ## _t; \
  } while (0)


#include "ptypes.h"
#define FUNC_isqrt 1
#define FUNC_icbrt 1
#define FUNC_lcm_ui 1
#define FUNC_ctz 1
#define FUNC_log2floor 1
#define FUNC_is_perfect_square
#define FUNC_next_prime_in_sieve 1
#define FUNC_prev_prime_in_sieve 1
#include "util.h"
#include "sieve.h"
#include "primality.h"
#include "cache.h"
#include "lmo.h"
#include "factor.h"
#include "mulmod.h"
#include "constants.h"

static int _verbose = 0;
void _XS_set_verbose(int v) { _verbose = v; }
int _XS_get_verbose(void) { return _verbose; }

static int _call_gmp = 0;
void _XS_set_callgmp(int v) { _call_gmp = v; }
int  _XS_get_callgmp(void) { return _call_gmp; }

/* GCC 3.4 - 4.1 has broken 64-bit popcount.
 * GCC 4.2+ can generate awful code when it doesn't have asm (GCC bug 36041).
 * When the asm is present (e.g. compile with -march=native on a platform that
 * has them, like Nahelem+), then it is almost as fast as the direct asm. */
#if BITS_PER_WORD == 64
 #if defined(__POPCNT__) && defined(__GNUC__) && (__GNUC__> 4 || (__GNUC__== 4 && __GNUC_MINOR__> 1))
   #define popcnt(b)  __builtin_popcountll(b)
 #else
   static UV popcnt(UV b) {
     b -= (b >> 1) & 0x5555555555555555;
     b = (b & 0x3333333333333333) + ((b >> 2) & 0x3333333333333333);
     b = (b + (b >> 4)) & 0x0f0f0f0f0f0f0f0f;
     return (b * 0x0101010101010101) >> 56;
   }
 #endif
#else
 static UV popcnt(UV b) {
   b -= (b >> 1) & 0x55555555;
   b = (b & 0x33333333) + ((b >> 2) & 0x33333333);
   b = (b + (b >> 4)) & 0x0f0f0f0f;
   return (b * 0x01010101) >> 24;
 }
#endif

#if defined(__GNUC__)
 #define word_unaligned(m,wordsize)  ((uintptr_t)m & (wordsize-1))
#else  /* uintptr_t is part of C99 */
 #define word_unaligned(m,wordsize)  ((unsigned int)m & (wordsize-1))
#endif


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



/* We'll use this little static sieve to quickly answer small values of
 *   is_prime, next_prime, prev_prime, prime_count
 * for non-threaded Perl it's basically the same as getting the primary
 * cache.  It guarantees we'll have an answer with no waiting on any version.
 */
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

/* Return of 2 if n is prime, 0 if not.  Do it fast. */
int _XS_is_prime(UV n)
{
  if (n <= 10)
    return (n == 2 || n == 3 || n == 5 || n == 7) ? 2 : 0;

  if (n < UVCONST(200000000)) {
    UV d = n/30;
    UV m = n - d*30;
    unsigned char mtab = masktab30[m];  /* Bitmask in mod30 wheel */
    const unsigned char* sieve;
    int isprime;

    /* Return 0 if a multiple of 2, 3, or 5 */
    if (mtab == 0)
      return 0;

    /* Check static tiny sieve */
    if (d < NPRIME_SIEVE30)
      return (prime_sieve30[d] & mtab) ? 0 : 2;

    if (!(n%7) || !(n%11) || !(n%13)) return 0;

    /* Check primary cache */
    isprime = -1;
    if (n <= get_prime_cache(0, &sieve))
      isprime = 2*((sieve[d] & mtab) == 0);
    release_prime_cache(sieve);
    if (isprime >= 0)
      return isprime;
  }
  return is_prob_prime(n);
}


UV next_prime(UV n)
{
  UV m, sieve_size, next;
  const unsigned char* sieve;

  if (n < 30*NPRIME_SIEVE30) {
    next = next_prime_in_sieve(prime_sieve30, n,30*NPRIME_SIEVE30);
    if (next != 0) return next;
  }

  if (n >= MPU_MAX_PRIME) return 0; /* Overflow */

  sieve_size = get_prime_cache(0, &sieve);
  next = (n < sieve_size)  ?  next_prime_in_sieve(sieve, n, sieve_size)  :  0;
  release_prime_cache(sieve);
  if (next != 0) return next;

  m = n % 30;
  do { /* Move forward one. */
    n += wheeladvance30[m];
    m = nextwheel30[m];
  } while (!is_prob_prime(n));
  return n;
}


UV prev_prime(UV n)
{
  const unsigned char* sieve;
  UV m, prev;

  if (n < 30*NPRIME_SIEVE30)
    return prev_prime_in_sieve(prime_sieve30, n);

  if (n < get_prime_cache(0, &sieve)) {
    prev = prev_prime_in_sieve(sieve, n);
    release_prime_cache(sieve);
    return prev;
  }
  release_prime_cache(sieve);

  m = n % 30;
  do { /* Move back one. */
    n -= wheelretreat[m];
    m = prevwheel30[m];
  } while (!is_prob_prime(n));
  return n;
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
  START_DO_FOR_EACH_SIEVE_PRIME(sieve, lowp, highp)
    count++;
  END_DO_FOR_EACH_SIEVE_PRIME;
  return count;
#endif

  lo_d = lowp/30;
  lo_m = lowp - lo_d*30;
  /* Count first fragment */
  if (lo_m > 1) {
    UV upper = (highp <= (lo_d*30+29)) ? highp : (lo_d*30+29);
    START_DO_FOR_EACH_SIEVE_PRIME(sieve, lowp, upper)
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

#define USE_PC_TABLES 1
#if USE_PC_TABLES
/* These tables let us have fast answers up to 3000M for the cost of ~1.4k of
 * static data/code.  We can get a 4 to 100x speedup here.  We don't want to
 * push this idea too far because Lehmer's method should be faster. */
/* mpu '$step=30_000; $pc=prime_count(5); print "$pc\n", join(",", map { $spc=$pc; $pc=prime_count($_*$step); $pc-$spc; } 1..200), "\n"' */
static const unsigned short step_counts_30k[] =  /* starts at 7 */
  {3242,2812,2656,2588,2547,2494,2465,2414,2421,2355,2407,2353,2310,2323,2316,
   2299,2286,2281,2247,2279,2243,2223,2251,2214,2209,2230,2215,2207,2205,2179,
   2200,2144,2159,2193,2164,2136,2180,2152,2162,2174,2113,2131,2150,2101,2111,
   2146,2115,2123,2119,2108,2124,2097,2075,2089,2094,2119,2084,2065,2069,2101,
   2094,2083,2089,2076,2088,2027,2109,2073,2061,2033,2079,2078,2036,2025,2058,
   2083,2037,2005,2048,2048,2024,2045,2027,2025,2039,2049,2022,2034,2046,2032,
   2019,2000,2014,2069,2042,1980,2021,2014,1995,2017,1992,1985,2045,2007,1990,
   2008,2052,2033,1988,1984,2010,1943,2024,2005,2027,1937,1955,1956,1993,1976,
   2048,1940,2002,2007,1994,1954,1972,2002,1973,1993,1984,1969,1940,1960,2026,
   1966,1981,1912,1994,1971,1977,1952,1932,1977,1932,1954,1938,2018,1987,1967,
   1937,1938,1963,1973,1947,1947,1963,1959,1941,1923,1943,1957,1974,1964,1958,
   1984,1933,1935,1935,1949,1928,1943,1917,1956,1970,1932,1937,1929,1932,1947,
   1927,1944,1915,1913,1918,1925,1931,1919,1900,1952,1934,1922,1891,1926,1925,
   1903,1970,1962,1905,1905};
#define NSTEP_COUNTS_30K  (sizeof(step_counts_30k)/sizeof(step_counts_30k[0]))

/* mpu '$step=300_000; $pc=prime_count(20*$step); print "$pc\n", join(",", map { $spc=$pc; $pc=prime_count($_*$step); $pc-$spc; } 21..212), "\n"' */
static const unsigned short step_counts_300k[] =  /* starts at 6M */
  {19224,19086,19124,19036,18942,18893,18870,18853,18837,18775,18688,18674,
   18594,18525,18639,18545,18553,18424,18508,18421,18375,18366,18391,18209,
   18239,18298,18209,18294,18125,18138,18147,18115,18126,18021,18085,18068,
   18094,17963,18041,18003,17900,17881,17917,17888,17880,17852,17892,17779,
   17823,17764,17806,17762,17780,17716,17633,17758,17746,17678,17687,17613,
   17709,17628,17634,17556,17528,17598,17604,17532,17606,17548,17493,17576,
   17456,17468,17555,17452,17407,17472,17415,17500,17508,17418,17463,17240,
   17345,17351,17380,17394,17379,17330,17322,17335,17354,17113,17210,17231,
   17238,17305,17268,17219,17281,17235,17119,17292,17161,17212,17166,17277,
   17137,17260,17228,17197,17154,17097,17195,17136,17067,17058,17041,17045,
   17187,17034,17029,17037,17090,16985,17054,17017,17106,17001,17095,17125,
   17027,16948,16969,17031,16916,17031,16905,16937,16881,16952,16919,16938,
   17028,16963,16902,16922,16944,16901,16847,16969,16900,16876,16841,16874,
   16894,16861,16761,16886,16778,16820,16727,16921,16817,16845,16847,16824,
   16844,16809,16859,16783,16713,16752,16762,16857,16760,16626,16784,16784,
   16718,16745,16871,16635,16714,16630,16779,16709,16660,16730,16715,16724};
#define NSTEP_COUNTS_300K (sizeof(step_counts_300k)/sizeof(step_counts_300k[0]))

static const unsigned int step_counts_30m[] =  /* starts at 60M */
  {1654839,1624694,1602748,1585989,1571241,1559918,1549840,1540941,1533150,
   1525813,1519922,1513269,1508559,1503386,1497828,1494129,1489905,1486417,
   1482526,1478941,1475577,1472301,1469133,1466295,1464711,1461223,1458478,
   1455327,1454218,1451883,1449393,1447612,1445029,1443285,1442268,1438511,
   1437688,1435603,1433623,1432638,1431158,1429158,1427934,1426191,1424449,
   1423146,1421898,1421628,1419519,1417646,1416274,1414828,1414474,1412536,
   1412147,1410149,1409474,1408847,1406619,1405863,1404699,1403820,1402802,
   1402215,1401459,1399972,1398687,1397968,1397392,1396025,1395311,1394081,
   1393614,1393702,1391745,1390950,1389856,1389245,1388381,1387557,1387087,
   1386285,1386089,1385355,1383659,1383030,1382174,1382128,1380556,1379940,
   1379988,1379181,1378300,1378033,1376974,1376282,1375646,1374445,1373813};
#define NSTEP_COUNTS_30M  (sizeof(step_counts_30m)/sizeof(step_counts_30m[0]))
#endif

UV _XS_prime_count(UV low, UV high)
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

  if (low == 7 && high <= 30*NPRIME_SIEVE30) {
    count += count_segment_ranged(prime_sieve30, NPRIME_SIEVE30, low, high);
    return count;
  }

#if USE_PC_TABLES
  if (low == 7 && high >= 30000) {
    UV i, maxi;
    if (high < (30000*(NSTEP_COUNTS_30K+1))) {
      low = 0;
      maxi = high/30000;
      for (i = 0;  i < maxi && i < NSTEP_COUNTS_30K;  i++) {
        count += step_counts_30k[i];
        low += 30000;
      }
    } else if (high < (6000000 + 300000*(NSTEP_COUNTS_300K+1))) {
      count = 412849;
      low = 6000000;
      maxi = (high-6000000)/300000;
      for (i = 0;  i < maxi && i < NSTEP_COUNTS_300K;  i++) {
        count += step_counts_300k[i];
        low += 300000;
      }
    } else {
      count = 3562115;
      low = 60000000;
      maxi = (high-60000000)/30000000;
      for (i = 0;  i < maxi && i < NSTEP_COUNTS_30M;  i++) {
        count += step_counts_30m[i];
        low += 30000000;
      }
    }
  }
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

UV prime_count_approx(UV n)
{
  if (n < 3000000) return _XS_prime_count(2, n);
  return (UV) (_XS_RiemannR( (long double) n ) + 0.5 );
}

/* See http://numbers.computation.free.fr/Constants/Primes/twin.pdf, page 5 */
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
    long double li2 = _XS_ExponentialIntegral(logn) + two_over_log_two-ln/logn;
    /* try to minimize MSE */
    if      (n <    4000) li2 *= 1.2035 * logl(logl(logl(ln)));
    else if (n <    8000) li2 *= 0.9439 * logl(logl(logl(ln*1000)));
    else if (n <   32000) li2 *= 0.8967 * logl(logl(logl(ln*4000)));
    else if (n <  200000) li2 *= 0.8937 * logl(logl(logl(ln*4000)));
    else if (n < 1000000) li2 *= 0.8640 * logl(logl(logl(ln*16000)));
    else if (n < 4000000) li2 *= 0.8627 * logl(logl(logl(ln*16000)));
    else if (n <10000000) li2 *= 0.8536 * logl(logl(logl(ln*16000)));
    return (UV) (two_C2 * li2 + 0.5L);
  }
}

UV prime_count_lower(UV n)
{
  long double fn, flogn, lower, a;

  if (n < 33000) return _XS_prime_count(2, n);

  fn     = (long double) n;
  flogn  = logl(n);

  if      (n <   176000)  a = 1.80;
  else if (n <   315000)  a = 2.10;
  else if (n <  1100000)  a = 2.20;
  else if (n <  4500000)  a = 2.31;
  else if (n <233000000)  a = 2.36;
#if BITS_PER_WORD == 32
  else a = 2.32;
#else
  else if (n < UVCONST( 5433800000)) a = 2.32;
  else if (n < UVCONST(60000000000)) a = 2.15;
  else a = 2.00;
#endif

  lower = fn/flogn * (1.0 + 1.0/flogn + a/(flogn*flogn));
  return (UV) floorl(lower);
}

typedef struct {
  UV thresh;
  float aval;
} thresh_t;

static const thresh_t _upper_thresh[] = {
  {     59000, 2.48 },
  {    350000, 2.52 },
  {    355991, 2.54 },
  {    356000, 2.51 },
  {   3550000, 2.50 },
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
  long double fn, flogn, upper, a;

  if (n < 33000) return _XS_prime_count(2, n);

  for (i = 0; i < (int)NUPPER_THRESH; i++)
    if (n < _upper_thresh[i].thresh)
      break;

  if (i < (int)NUPPER_THRESH) a = _upper_thresh[i].aval;
  else                        a = 2.334;   /* Dusart 2010, page 2 */

  fn     = (long double) n;
  flogn  = logl(n);
  upper = fn/flogn * (1.0 + 1.0/flogn + a/(flogn*flogn));
  return (UV) ceill(upper);
}

static const unsigned short primes_small[] =
  {0,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
   101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,
   193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,
   293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,
   409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499};
#define NPRIMES_SMALL (sizeof(primes_small)/sizeof(primes_small[0]))

/* The nth prime will be less or equal to this number */
UV nth_prime_upper(UV n)
{
  long double fn, flogn, flog2n, upper;

  if (n < NPRIMES_SMALL)
    return primes_small[n];

  fn     = (long double) n;
  flogn  = logl(n);
  flog2n = logl(flogn);    /* Note distinction between log_2(n) and log^2(n) */

  if      (n >= 688383)    /* Dusart 2010 page 2 */
    upper = fn * (flogn + flog2n - 1.0 + ((flog2n-2.00)/flogn));
  else if (n >= 178974)    /* Dusart 2010 page 7 */
    upper = fn * (flogn + flog2n - 1.0 + ((flog2n-1.95)/flogn));
  else if (n >=  39017)    /* Dusart 1999 page 14 */
    upper = fn * (flogn + flog2n - 0.9484);
  else if (n >=     6)     /* Modified from Robin 1983 for 50-39016 _only_ */
    upper = fn * ( flogn  +  0.5982 * flog2n ) - 5;
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
  if (upper >= (long double)UV_MAX) {
    if (n <= MPU_MAX_PRIME_IDX) return MPU_MAX_PRIME;
    croak("nth_prime_upper(%"UVuf") overflow", n);
  }

  return (UV) ceill(upper);
}

/* The nth prime will be greater than or equal to this number */
UV nth_prime_lower(UV n)
{
  long double fn, flogn, flog2n, lower;

  if (n < NPRIMES_SMALL)
    return primes_small[n];

  fn     = (long double) n;
  flogn  = logl(n);
  flog2n = logl(flogn);    /* Note distinction between log_2(n) and log^2(n) */

  /* Dusart 2010 page 2, for all n >= 3 */
  lower = fn * (flogn + flog2n - 1.0 + ((flog2n-2.10)/flogn));
  /* Tighten small values */
  if      (n <  2679) lower = 1.003 * lower + 23;
  else if (n < 14353) lower = 1.001 * lower + 21;

  return (UV) floorl(lower);
}

UV nth_prime_approx(UV n)
{
  long double fn, flogn;
  UV lo, hi;

  if (n < NPRIMES_SMALL)
    return primes_small[n];

  /* Binary search for inverse Riemann R */
  fn    = (long double) n;
  flogn = logl(n);
  lo    = (UV) (fn * flogn);
  hi    = (UV) (fn * flogn * 2 + 2);
  if (hi <= lo) hi = UV_MAX;
  while (lo < hi) {
    UV mid = lo + (hi-lo)/2;
    if (_XS_RiemannR(mid) < fn) lo = mid+1;
    else                        hi = mid;
  }
  return lo;
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
    UV lower_limit = _XS_Inverse_Li(n);
    /* For even better performance, add in half the usual correction, which
     * will get us even closer, so even less sieving required.  However, it
     * is now possible to get a result higher than the value, so we'll need
     * to handle that case.  It still ends up being a better deal than R,
     * given that we don't have a fast backward sieve. */
    lower_limit += _XS_Inverse_Li(isqrt(n))/4;
    segment_size = lower_limit / 30;
    lower_limit = 30 * segment_size - 1;
    count = _XS_LMO_pi(lower_limit);

    /* printf("We've estimated %lu too %s.\n", (count>n)?count-n:n-count, (count>n)?"FAR":"little"); */
    /* printf("Our limit %lu %s a prime\n", lower_limit, _XS_is_prime(lower_limit) ? "is" : "is not"); */

    if (count >= n) { /* Too far.  Walk backwards */
      if (_XS_is_prime(lower_limit)) count--;
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
      if (_XS_is_prime(beg) && _XS_is_prime(beg+2) && beg <= end) sum++;
      beg += 2;
    }
    while ((end % 30) != 29) {
      if (_XS_is_prime(end) && _XS_is_prime(end+2) && beg <= end) sum++;
      end -= 2;  if (beg > end) break;
    }
  }
  if (beg <= end) {
    UV seg_base, seg_low, seg_high;
    void* ctx = start_segment_primes(beg, end, &segment);
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      UV bytes = (seg_high-seg_low+29)/30;
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
      if (!(s & 0x80) && _XS_is_prime(seg_high+2)) sum++;
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
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      UV p, bytes = (seg_high-seg_low+29)/30;
      UV s = ((UV)segment[0]) << 8;
      for (p = 0; p < bytes; p++) {
        s >>= 8;
        if (p+1 < bytes)                    s |= (((UV)segment[p+1]) << 8);
        else if (!_XS_is_prime(seg_high+2)) s |= 0xFF00;
        if (!(s & 0x000C) && !--n) { nth=seg_base+p*30+11; break; }
        if (!(s & 0x0030) && !--n) { nth=seg_base+p*30+17; break; }
        if (!(s & 0x0180) && !--n) { nth=seg_base+p*30+29; break; }
      }
      if (n == 0) break;
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
  lo = (UV) (1.2 * fnlog2n);
  hi = (UV) ( (n >= 1200) ? (1.7 * fnlog2n) : (2.3 * fnlog2n + 5) );
  if (hi <= lo) hi = UV_MAX;
  while (lo < hi) {
    UV mid = lo + (hi-lo)/2;
    if (twin_prime_count_approx(mid) < fn) lo = mid+1;
    else                                   hi = mid;
  }
  return lo;
}


/* Return a char array with lo-hi+1 elements. mu[k-lo] = µ(k) for k = lo .. hi.
 * It is the callers responsibility to call Safefree on the result. */
#define PGTLO(p,lo)  ((p) >= lo) ? (p) : ((p)*(lo/(p)) + ((lo%(p))?(p):0))
#define P2GTLO(pinit, p, lo) \
   ((pinit) >= lo) ? (pinit) : ((p)*(lo/(p)) + ((lo%(p))?(p):0))
signed char* _moebius_range(UV lo, UV hi)
{
  signed char* mu;
  UV i;
  UV sqrtn = isqrt(hi);

  /* Kuznetsov indicates that the Deléglise & Rivat (1996) method can be
   * modified to work on logs, which allows us to operate with no
   * intermediate memory at all.  Same time as the D&R method, less memory. */
  unsigned char logp;
  UV nextlog;

  Newz(0, mu, hi-lo+1, signed char);
  if (sqrtn*sqrtn != hi) sqrtn++;  /* ceil sqrtn */

  logp = 1; nextlog = 3; /* 2+1 */
  START_DO_FOR_EACH_PRIME(2, sqrtn) {
    UV p2 = p*p;
    if (p > nextlog) {
      logp += 2;   /* logp is 1 | ceil(log(p)/log(2)) */
      nextlog = ((nextlog-1)*4)+1;
    }
    for (i = PGTLO(p, lo); i <= hi; i += p)
      mu[i-lo] += logp;
    for (i = PGTLO(p2, lo); i <= hi; i += p2)
      mu[i-lo] |= 0x80;
  } END_DO_FOR_EACH_PRIME

  logp = log2floor(lo);
  nextlog = 2UL << logp;
  for (i = lo; i <= hi; i++) {
    unsigned char a = mu[i-lo];
    if (i >= nextlog) {  logp++;  nextlog *= 2;  } /* logp is log(p)/log(2) */
    if (a & 0x80)       { a = 0; }
    else if (a >= logp) { a =  1 - 2*(a&1); }
    else                { a = -1 + 2*(a&1); }
    mu[i-lo] = a;
  }
  if (lo == 0)  mu[0] = 0;

  return mu;
}

UV* _totient_range(UV lo, UV hi) {
  UV* totients;
  UV i, seg_base, seg_low, seg_high;
  unsigned char* segment;
  void* ctx;

  if (hi < lo) croak("_totient_range error hi %"UVuf" < lo %"UVuf"\n", hi, lo);
  New(0, totients, hi-lo+1, UV);

  /* Do via factoring if very small or if we have a small range */
  if (hi < 100 || (hi-lo) < 10 || hi/(hi-lo+1) > 1000) {
    for (i = lo; i <= hi; i++)
      totients[i-lo] = totient(i);
    return totients;
  }

  /* If doing a full sieve, do it monolithic.  Faster. */
  if (lo == 0) {
    UV* prime;
    double loghi = log(hi);
    UV max_index = (hi < 67)     ? 18
                 : (hi < 355991) ? 15+(hi/(loghi-1.09))
                 : (hi/loghi) * (1.0+1.0/loghi+2.51/(loghi*loghi));
    UV j, index, nprimes = 0;

    New(0, prime, max_index, UV);  /* could use prime_count_upper(hi) */
    memset(totients, 0, (hi-lo+1) * sizeof(UV));
    for (i = 2; i <= hi/2; i++) {
      index = 2*i;
      if ( !(i&1) ) {
        if (i == 2) { totients[2] = 1; prime[nprimes++] = 2; }
        totients[index] = totients[i]*2;
      } else {
        if (totients[i] == 0) {
          totients[i] = i-1;
          prime[nprimes++] = i;
        }
        for (j=0; j < nprimes && index <= hi; index = i*prime[++j]) {
          if (i % prime[j] == 0) {
            totients[index] = totients[i]*prime[j];
            break;
          } else {
            totients[index] = totients[i]*(prime[j]-1);
          }
        }
      }
    }
    Safefree(prime);
    /* All totient values have been filled in except the primes.  Mark them. */
    for (i = ((hi/2) + 1) | 1; i <= hi; i += 2)
      if (totients[i] == 0)
        totients[i] = i-1;
    totients[1] = 1;
    totients[0] = 0;
    return totients;
  }

  for (i = lo; i <= hi; i++) {
    UV v = i;
    if (i % 2 == 0)  v -= v/2;
    if (i % 3 == 0)  v -= v/3;
    if (i % 5 == 0)  v -= v/5;
    totients[i-lo] = v;
  }

  ctx = start_segment_primes(7, hi/2, &segment);
  while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
    START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_low - seg_base, seg_high - seg_base ) {
      p += seg_base;
      for (i = P2GTLO(2*p,p,lo); i <= hi; i += p)
        totients[i-lo] -= totients[i-lo]/p;
    } END_DO_FOR_EACH_SIEVE_PRIME
  }
  end_segment_primes(ctx);

  /* Fill in all primes */
  for (i = lo | 1; i <= hi; i += 2)
    if (totients[i-lo] == i)
      totients[i-lo]--;
  if (lo <= 1) totients[1-lo] = 1;

  return totients;
}

IV mertens(UV n) {
  /* See Deléglise and Rivat (1996) for O(n^2/3 log(log(n))^1/3) algorithm.
   * This implementation uses their lemma 2.1 directly, so is ~ O(n).
   * In serial it is quite a bit faster than segmented summation of mu
   * ranges, though the latter seems to be a favored method for GPUs.
   */
  UV u, i, m, nmk, maxmu;
  signed char* mu;
  IV* M;
  IV sum;

  if (n <= 1)  return n;
  u = isqrt(n);
  maxmu = (n/(u+1));              /* maxmu lets us handle u < sqrt(n) */
  if (maxmu < u) maxmu = u;
  mu = _moebius_range(0, maxmu);
  New(0, M, maxmu+1, IV);
  M[0] = 0;
  for (i = 1; i <= maxmu; i++)
    M[i] = M[i-1] + mu[i];
  sum = M[u];
  for (m = 1; m <= u; m++) {
    if (mu[m] != 0) {
      IV inner_sum = 0;
      UV lower = (u/m) + 1;
      UV last_nmk = n/(m*lower);
      UV this_k = 0;
      UV next_k = n/(m*1);
      UV nmkm = m * 2;
      for (nmk = 1; nmk <= last_nmk; nmk++, nmkm += m) {
        this_k = next_k;
        next_k = n/nmkm;
        inner_sum += M[nmk] * (this_k - next_k);
      }
      sum -= mu[m] * inner_sum;
    }
  }
  Safefree(M);
  Safefree(mu);
  return sum;
}

/* There are at least 4 ways to do this, plus hybrids.
 * 1) use a table.  Great for 32-bit, too big for 64-bit.
 * 2) Use pow() to check.  Relatively slow and FP is always dangerous.
 * 3) factor or trial factor.  Slow for 64-bit.
 * 4) Dietzfelbinger algorithm 2.3.5.  Quite slow.
 * This currently uses a hybrid of 1 and 2.
 */
int powerof(UV n) {
  int ib;
  const int iblast = (n > UVCONST(4294967295)) ? 6 : 4;
  if ((n <= 3) || (n == UV_MAX)) return 1;
  if ((n & (n-1)) == 0)          return ctz(n);  /* powers of 2    */
  if (is_perfect_square(n))      return 2 * powerof(isqrt(n));
  { UV cb = icbrt(n);  if (cb*cb*cb==n) return 3 * powerof(cb); }
  for (ib = 3; ib <= iblast; ib++) { /* prime exponents from 5 to 7-or-13 */
    UV k, pk, root, b = primes_small[ib];
    root = (UV) ( pow(n, 1.0 / b ) + 0.01 );
    pk = root * root * root * root * root;
    for (k = 5; k < b; k++)
      pk *= root;
    if (n == pk) return b * powerof(root);
  }
  if (n > 177146) {
    switch (n) { /* Check for powers of 11, 13, 17, 19 within 32 bits */
      case 177147: case 48828125: case 362797056: case 1977326743: return 11;
      case 1594323: case 1220703125: return 13;
      case 129140163: return 17;
      case 1162261467: return 19;
      default:  break;
    }
#if BITS_PER_WORD == 64
    if (n > UVCONST(4294967295)) {
    switch (n) {
      case UVCONST(762939453125):
      case UVCONST(16926659444736):
      case UVCONST(232630513987207):
      case UVCONST(100000000000000000):
      case UVCONST(505447028499293771):
      case UVCONST(2218611106740436992):
      case UVCONST(8650415919381337933):  return 17;
      case UVCONST(19073486328125):
      case UVCONST(609359740010496):
      case UVCONST(11398895185373143):
      case UVCONST(10000000000000000000): return 19;
      case UVCONST(94143178827):
      case UVCONST(11920928955078125):
      case UVCONST(789730223053602816):   return 23;
      case UVCONST(68630377364883):       return 29;
      case UVCONST(617673396283947):      return 31;
      case UVCONST(450283905890997363):   return 37;
      default:  break;
    }
    }
#endif
  }
  return 1;
}
int is_power(UV n, UV a)
{
  int ret;
  if (a > 0) {
    if (a == 1 || n <= 1) return 1;
    if ((a % 2) == 0)
      return !is_perfect_square(n) ? 0 : (a == 2) ? 1 : is_power(isqrt(n),a>>1);
    if ((a % 3) == 0)
      { UV cb = icbrt(n);
        return (cb*cb*cb != n)       ? 0 : (a == 3) ? 1 : is_power(cb, a/3); }
    if ((a % 5) == 0)
      { UV r5 = (UV)(pow(n,0.2) + 0.0001);
        return (r5*r5*r5*r5*r5 != n) ? 0 : (a == 5) ? 1 : is_power(r5, a/5); }
  }
  ret = powerof(n);
  if (a != 0) return !(ret % a);  /* Is the max power divisible by a? */
  return (ret == 1) ? 0 : ret;
}
UV rootof(UV n, UV k)
{
  switch (k) {
    case 0:  return 0;
    case 1:  return n;
    case 2:  return isqrt(n);
    case 3:  return icbrt(n);
    default: return (UV) powl((long double)n, 1.0L/(long double)k);
  }
}

UV valuation(UV n, UV k)
{
  UV v = 0;
  UV kpower = k;
  if (k < 2 || n < 2) return 0;
  if (k == 2) return ctz(n);
  while ( !(n % kpower) ) {
    kpower *= k;
    v++;
  }
  return v;
}

UV mpu_popcount(UV n) {  return popcnt(n);  }

UV mpu_popcount_string(const char* ptr, int len)
{
  int i, *s, *sptr;
  UV count = 0;

  while (len > 0 && (*ptr == '0' || *ptr == '+' || *ptr == '-'))
    {  ptr++;  len--;  }

  New(0, s, len, int);
  for (i = 0; i < len; i++)
    s[i] = ptr[i] - '0';

  while (len > 0) {
    if (s[len-1] & 1)  count++;
    /* divide by 2 */
    sptr = s;
    if (s[0] == 1) {
      if (--len == 0) break;
      *++sptr += 10;
    }
    for (i = 0; i < len; i++) {
      if ( (i+1) < len  &&  sptr[i] & 1 ) sptr[i+1] += 10;
      s[i] = sptr[i] / 2;
    }
  }
  Safefree(s);
  return count;
}


/* How many times does 2 divide n? */
#define padic2(n)  ctz(n)
#define IS_MOD8_3OR5(x)  (((x)&7)==3 || ((x)&7)==5)

static int kronecker_uu_sign(UV a, UV b, int s) {
  while (a) {
    int r = padic2(a);
    if (r) {
      if ((r&1)  &&  IS_MOD8_3OR5(b))  s = -s;
      a >>= r;
    }
    if (a & b & 2)  s = -s;
    { UV t = b % a;  b = a;  a = t; }
  }
  return (b == 1) ? s : 0;
}

int kronecker_uu(UV a, UV b) {
  int r, s;
  if (b & 1)   return kronecker_uu_sign(a, b, 1);
  if (!(a&1))  return 0;
  s = 1;
  r = padic2(b);
  if (r) {
    if ((r&1) && IS_MOD8_3OR5(a))  s = -s;
    b >>= r;
  }
  return kronecker_uu_sign(a, b, s);
}

int kronecker_su(IV a, UV b) {
  int r, s;
  if (a >= 0)  return kronecker_uu(a, b);
  if (b == 0)  return (a == 1 || a == -1) ? 1 : 0;
  s = 1;
  r = padic2(b);
  if (r) {
    if (!(a&1))  return 0;
    if ((r&1) && IS_MOD8_3OR5(a))  s = -s;
    b >>= r;
  }
  a %= (IV) b;
  if (a < 0)  a += b;
  return kronecker_uu_sign(a, b, s);
}

int kronecker_ss(IV a, IV b) {
  if (a >= 0 && b >= 0)
    return (b & 1)  ?  kronecker_uu_sign(a, b, 1)  :  kronecker_uu(a,b);
  if (b >= 0)
    return kronecker_su(a, b);
  return kronecker_su(a, -b) * ((a < 0) ? -1 : 1);
}

UV factorial(UV n) {
  UV i, r = 1;
  if ( (n > 12 && sizeof(UV) <= 4) || (n > 20 && sizeof(UV) <= 8) ) return 0;
  for (i = 2; i <= n; i++)
    r *= i;
  return r;
}

UV binomial(UV n, UV k) {    /* Thanks to MJD and RosettaCode for ideas */
  UV d, g, r = 1;
  if (k == 0) return 1;
  if (k == 1) return n;
  if (k >= n) return (k == n);
  if (k > n/2) k = n-k;
  for (d = 1; d <= k; d++) {
    if (r >= UV_MAX/n) {  /* Possible overflow */
      UV nr, dr;  /* reduced numerator / denominator */
      g = gcd_ui(n, d);  nr = n/g;  dr = d/g;
      g = gcd_ui(r, dr);  r = r/g;  dr = dr/g;
      if (r >= UV_MAX/nr) return 0;  /* Unavoidable overflow */
      r *= nr;
      r /= dr;
      n--;
    } else {
      r *= n--;
      r /= d;
    }
  }
  return r;
}

IV stirling2(UV n, UV m) {
  UV f;
  IV j, k, t, s = 0;

  if (m == n) return 1;
  if (n == 0 || m == 0 || m > n) return 0;
  if (m == 1) return 1;
  if ((f = factorial(m)) == 0) return 0;
  for (j = 1; j <= (IV)m; j++) {
    t = binomial(m, j);
    for (k = 1; k <= (IV)n; k++) {
      if (t == 0 || j >= IV_MAX/t) return 0;
      t *= j;
    }
    if ((m-j) & 1) t *= -1;
    s += t;
  }
  return s/f;
}

IV stirling1(UV n, UV m) {
  IV k, t, b1, b2, s2, s = 0;
  for (k = 1; k <= (IV)(n-m); k++) {
    b1 = binomial(k + n - 1, k + n - m);
    b2 = binomial(2 * n - m, n - k - m);
    s2 = stirling2(k - m + n, k);
    if (b1 == 0 || b2 == 0 || s2 == 0 || b1 > IV_MAX/b2) return 0;
    t = b1 * b2;
    if (s2 > IV_MAX/t) return 0;
    t *= s2;
    s += (k & 1) ? -t : t;
  }
  return s;
}

UV totient(UV n) {
  UV i, nfacs, totient, lastf, facs[MPU_MAX_FACTORS+1];
  if (n <= 1) return n;
  totient = 1;
  /* phi(2m) = 2phi(m) if m even, phi(m) if m odd */
  while ((n & 0x3) == 0) { n >>= 1; totient <<= 1; }
  if ((n & 0x1) == 0) { n >>= 1; }
  /* factor and calculate totient */
  nfacs = factor(n, facs);
  lastf = 0;
  for (i = 0; i < nfacs; i++) {
    UV f = facs[i];
    if (f == lastf) { totient *= f;               }
    else            { totient *= f-1;  lastf = f; }
  }
  return totient;
}

static const UV jordan_overflow[5] =
#if BITS_PER_WORD == 64
  {UVCONST(4294967311), 2642249, 65537, 7133, 1627};
#else
  {UVCONST(     65537),    1627,   257,   85,   41};
#endif
UV jordan_totient(UV k, UV n) {
  UV factors[MPU_MAX_FACTORS+1];
  int nfac, i;
  UV j, totient;
  if (k == 0 || n <= 1) return (n == 1);
  if (k > 6 || (k > 1 && n >= jordan_overflow[k-2])) return 0;

  totient = 1;
  /* Similar to Euler totient, shortcut even inputs */
  while ((n & 0x3) == 0) { n >>= 1; totient *= (1<<k); }
  if ((n & 0x1) == 0) { n >>= 1; totient *= ((1<<k)-1); }
  nfac = factor(n,factors);
  for (i = 0; i < nfac; i++) {
    UV p = factors[i];
    UV pk = p;
    for (j = 1; j < k; j++)  pk *= p;
    totient *= (pk-1);
    while (i+1 < nfac && p == factors[i+1]) {
      i++;
      totient *= pk;
    }
  }
  return totient;
}

UV carmichael_lambda(UV n) {
  UV fac[MPU_MAX_FACTORS+1];
  int i, nfactors;
  UV lambda = 1;

  if (n < 8) return totient(n);
  if ((n & (n-1)) == 0) return n >> 2;

  i = ctz(n);
  if (i > 0) {
    n >>= i;
    lambda <<= (i>2) ? i-2 : i-1;
  }
  nfactors = factor(n, fac);
  for (i = 0; i < nfactors; i++) {
    UV p = fac[i], pk = p-1;
    while (i+1 < nfactors && p == fac[i+1]) {
      i++;
      pk *= p;
    }
    lambda = lcm_ui(lambda, pk);
  }
  return lambda;
}

int moebius(UV n) {
  UV factors[MPU_MAX_FACTORS+1];
  UV i, nfactors;

  if (n <= 1) return (int)n;
  if ( n >= 49 && (!(n% 4) || !(n% 9) || !(n%25) || !(n%49)) )
    return 0;

  nfactors = factor(n, factors);
  for (i = 1; i < nfactors; i++)
    if (factors[i] == factors[i-1])
      return 0;
  return (nfactors % 2) ? -1 : 1;
}

UV exp_mangoldt(UV n) {
  if      (n <= 1)           return 1;
  else if ((n & (n-1)) == 0) return 2;     /* Power of 2 */
  else if ((n & 1) == 0)     return 1;     /* Even number (not 2) */
  else if (is_prob_prime(n)) return n;
  else {
    int k = powerof(n);
    if (k >= 2) {
      n = (k==2) ? isqrt(n) : (UV)(pow(n,1.0/k)+0.0000001);
      if (is_prob_prime(n)) return n;
    }
    return 1;
  }
}


UV znorder(UV a, UV n) {
  UV fac[MPU_MAX_FACTORS+1];
  UV exp[MPU_MAX_FACTORS+1];
  int i, nfactors;
  UV j, k, phi;

  if (n <= 1) return n;   /* znorder(x,0) = 0, znorder(x,1) = 1          */
  if (a <= 1) return a;   /* znorder(0,x) = 0, znorder(1,x) = 1  (x > 1) */
  if (gcd_ui(a,n) > 1)  return 0;

  /* Cohen 1.4.3 using Carmichael Lambda */
  phi = carmichael_lambda(n);
  nfactors = factor_exp(phi, fac, exp);
  k = phi;
  for (i = 0; i < nfactors; i++) {
    UV b, a1, ek, pi = fac[i], ei = exp[i];
    b = pi; for (j = 1; j < ei; j++)  b *= pi;
    k /= b;
    a1 = powmod(a, k, n);
    for (ek = 0; a1 != 1 && ek++ <= ei; a1 = powmod(a1, pi, n))
      k *= pi;
    if (ek > ei) return 0;
  }
  return k;
}

UV znprimroot(UV n) {
  UV fac[MPU_MAX_FACTORS+1];
  UV exp[MPU_MAX_FACTORS+1];
  UV a, phi;
  int i, nfactors;
  if (n <= 4) return (n == 0) ? 0 : n-1;
  if (n % 4 == 0)  return 0;
  if (is_prob_prime(n)) {
    phi = n-1;
  } else {  /* prim root exists if n is 2, 4, p^a, or 2(p^a) for odd prime p */
    if (n & 0x3) nfactors = factor_exp( (n&1) ? n : n>>1, fac, exp);
    else         nfactors = 0;
    if (nfactors != 1) return 0;
    phi = fac[0]-1;  /* n = p^a for odd prime p.  Calculate totient. */
    for (i = 1; i < (int)exp[0]; i++)
      phi *= fac[0];
  }
  nfactors = factor_exp(phi, fac, exp);
  for (i = 0; i < nfactors; i++)
    exp[i] = phi / fac[i];  /* exp[i] = phi(n) / i-th-factor-of-phi(n) */
  for (a = 2; a < n; a++) {
    if (kronecker_uu(a, n) == 0)  continue;
    for (i = 0; i < nfactors; i++)
      if (powmod(a, exp[i], n) == 1)
        break;
    if (i == nfactors) return a;
  }
  return 0;
}

IV gcdext(IV a, IV b, IV* u, IV* v, IV* cs, IV* ct) {
  IV s = 0;  IV os = 1;
  IV t = 1;  IV ot = 0;
  IV r = b;  IV or = a;
  while (r != 0) {
    IV quot = or / r;
    { IV tmp = r; r = or - quot * r;  or = tmp; }
    { IV tmp = s; s = os - quot * s;  os = tmp; }
    { IV tmp = t; t = ot - quot * t;  ot = tmp; }
  }
  if (or < 0) /* correct sign */
    { or = -or; os = -os; ot = -ot; }
  if (u  != 0) *u = os;
  if (v  != 0) *v = ot;
  if (cs != 0) *cs = s;
  if (ct != 0) *ct = t;
  return or;
}

/* Calculate 1/a mod n. */
UV modinverse(UV a, UV n) {
  IV t = 0;  UV nt = 1;
  UV r = n;  UV nr = a;
  while (nr != 0) {
    UV quot = r / nr;
    { UV tmp = nt;  nt = t - quot*nt;  t = tmp; }
    { UV tmp = nr;  nr = r - quot*nr;  r = tmp; }
  }
  if (r > 1) return 0;  /* No inverse */
  if (t < 0) t += n;
  return t;
}

UV divmod(UV a, UV b, UV n) {   /* a / b  mod n */
  UV binv = modinverse(b, n);
  if (binv == 0)  return 0;
  return mulmod(a, binv, n);
}

/* status: 1 ok, -1 no inverse, 0 overflow */
UV chinese(UV* a, UV* n, UV num, int* status) {
  UV p, gcd, i, j, lcm, sum;
  *status = 1;
  if (num == 0) return 0;

  /* Sort modulii, largest first */
  for (i = 1; i < num; i++)
    for (j = i; j > 0 && n[j-1] < n[j]; j--)
      { p=n[j-1]; n[j-1]=n[j]; n[j]=p;   p=a[j-1]; a[j-1]=a[j]; a[j]=p; }

  if (n[0] > IV_MAX) { *status = 0; return 0; }
  lcm = n[0]; sum = a[0] % n[0];
  for (i = 1; i < num; i++) {
    IV u, v, t, s;
    UV vs, ut;
    gcd = gcdext(lcm, n[i], &u, &v, &s, &t);
    if (gcd != 1 && ((sum % gcd) != (a[i] % gcd))) { *status = -1; return 0; }
    if (s < 0) s = -s;
    if (t < 0) t = -t;
    if (s > (IV)(IV_MAX/lcm)) { *status = 0; return 0; }
    lcm *= s;
    if (u < 0) u += lcm;
    if (v < 0) v += lcm;
    vs = mulmod((UV)v, (UV)s, lcm);
    ut = mulmod((UV)u, (UV)t, lcm);
    sum = addmod(  mulmod(vs, sum, lcm),  mulmod(ut, a[i], lcm),  lcm  );
  }
  return sum;
}

long double chebyshev_function(UV n, int which)
{
  long double logp, logn = logl(n);
  UV sqrtn = which ? isqrt(n) : 0;  /* for theta, p <= sqrtn always false */
  KAHAN_INIT(sum);

  if (n < primes_small[NPRIMES_SMALL-1]) {
    UV p, pi;
    for (pi = 1;  (p = primes_small[pi]) <= n; pi++) {
      logp = logl(p);
      if (p <= sqrtn) logp *= floorl(logn/logp+1e-15);
      KAHAN_SUM(sum, logp);
    }
  } else {
    UV seg_base, seg_low, seg_high;
    unsigned char* segment;
    void* ctx;
    long double logl2 = logl(2);
    long double logl3 = logl(3);
    long double logl5 = logl(5);
    if (!which) {
      KAHAN_SUM(sum,logl2); KAHAN_SUM(sum,logl3); KAHAN_SUM(sum,logl5);
    } else {
      KAHAN_SUM(sum, logl2 * floorl(logn/logl2 + 1e-15));
      KAHAN_SUM(sum, logl3 * floorl(logn/logl3 + 1e-15));
      KAHAN_SUM(sum, logl5 * floorl(logn/logl5 + 1e-15));
    }
    ctx = start_segment_primes(7, n, &segment);
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_low - seg_base, seg_high - seg_base ) {
        p += seg_base;
        logp = logl(p);
        if (p <= sqrtn) logp *= floorl(logn/logp+1e-15);
        KAHAN_SUM(sum, logp);
      } END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
  }
  return sum;
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

long double _XS_ExponentialIntegral(long double x) {
  long double val, term;
  unsigned int n;
  KAHAN_INIT(sum);

  if (x == 0) croak("Invalid input to ExponentialIntegral:  x must be != 0");
  /* Protect against messed up rounding modes */
  if (x >=  12000) return INFINITY;
  if (x <= -12000) return 0;

  if (x < -1) {
    /* Continued fraction, good for x < -1 */
    long double lc = 0;
    long double ld = 1.0L / (1.0L - (long double)x);
    val = ld * (-expl(x));
    for (n = 1; n <= 100000; n++) {
      long double old, t, n2;
      t = (long double)(2*n + 1) - (long double) x;
      n2 = n * n;
      lc = 1.0L / (t - n2 * lc);
      ld = 1.0L / (t - n2 * ld);
      old = val;
      val *= ld/lc;
      if ( fabsl(val-old) <= LDBL_EPSILON*fabsl(val) )
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
  } else if (x < -logl(LDBL_EPSILON)) {
    /* Convergent series */
    long double fact_n = x;
    for (n = 2; n <= 200; n++) {
      long double invn = 1.0L / n;
      fact_n *= (long double)x * invn;
      term = fact_n * invn;
      KAHAN_SUM(sum, term);
      /* printf("C  after adding %.20Lf, val = %.20Lf\n", term, sum); */
      if ( term < LDBL_EPSILON*sum) break;
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
      if (term < LDBL_EPSILON*sum) break;
      if (term < last_term) {
        KAHAN_SUM(sum, term);
        /* printf("A  after adding %.20llf, sum = %.20llf\n", term, sum); */
      } else {
        KAHAN_SUM(sum, (-last_term/3) );
        /* printf("A  after adding %.20llf, sum = %.20llf\n", -last_term/3, sum); */
        break;
      }
    }
    term = expl(x) * invx;
    val = term * sum + term;
  }

  return val;
}

long double _XS_LogarithmicIntegral(long double x) {
  if (x == 0) return 0;
  if (x == 1) return -INFINITY;
  if (x == 2) return li2;
  if (x < 0) croak("Invalid input to LogarithmicIntegral:  x must be >= 0");
  if (x >= LDBL_MAX) return INFINITY;

#if 0
  /* We could calculate this here directly using Ramaujan's series.
   * I did not find this any faster or more accuracte than using Ei. */
  if (x > 1) {
    long double flogx, numer, denom, factn, inner_sum, power2, term;
    unsigned int n, k;
    KAHAN_INIT(sum);

    flogx = logl(x);
    numer = flogx;
    inner_sum = 1.0L;
    factn = 1.0L;
    power2 = 2.0L;
    KAHAN_SUM(sum, flogx);
    for (n = 2, k = 1; n < 10000; n++) {
      factn *= n;
      numer *= -flogx;
      denom = factn * power2;
      power2 *= 2;
      for (; k <= (n-1) >> 1; k++)
        inner_sum += 1.0L / (2 * k + 1);
      term = (numer / denom) * inner_sum;
      KAHAN_SUM(sum, term);
      if (fabsl(term) < LDBL_EPSILON*fabsl(sum)) break;
    }
    return euler_mascheroni + logl(flogx) + sqrtl(x) * sum;
  }
#endif

  return _XS_ExponentialIntegral(logl(x));
}

/* Thanks to Kim Walisch for this idea */
UV _XS_Inverse_Li(UV x) {
  double nlogn = (double)x * log((double)x);
  UV lo = (UV) (nlogn);
  UV hi = (UV) (nlogn * 2 + 2);

  if (x == 0)  return 0;
  if (hi <= lo) hi = UV_MAX;
  while (lo < hi) {
    UV mid = lo + (hi-lo)/2;
    if (_XS_LogarithmicIntegral(mid) < x) lo = mid+1;
    else                                  hi = mid;
  }
  return lo;
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
  0.0000000000004547473783042154026799112L,
  0.0000000000002273736845824652515226821L,
  0.0000000000001136868407680227849349105L,
  0.0000000000000568434198762758560927718L,
  0.0000000000000284217097688930185545507L,
  0.0000000000000142108548280316067698343L,
  0.00000000000000710542739521085271287735L,
  0.00000000000000355271369133711367329847L,
  0.00000000000000177635684357912032747335L,
  0.000000000000000888178421093081590309609L,
  0.000000000000000444089210314381336419777L,
  0.000000000000000222044605079804198399932L,
  0.000000000000000111022302514106613372055L,
  0.0000000000000000555111512484548124372374L,
  0.0000000000000000277555756213612417258163L,
  0.0000000000000000138777878097252327628391L,
};
#define NPRECALC_ZETA (sizeof(riemann_zeta_table)/sizeof(riemann_zeta_table[0]))

/* Riemann Zeta on the real line, with 1 subtracted.
 * Compare to Math::Cephes zetac.  Also zeta with q=1 and subtracting 1.
 *
 * The Cephes zeta function uses a series (2k)!/B_2k which converges rapidly
 * and has a very wide range of values.  We use it here for some values.
 *
 * Note: Calculations here are done on long doubles and we try to generate as
 *       much accuracy as possible.  They will get returned to Perl as an NV,
 *       which is typically a 64-bit double with 15 digits.
 *
 * For values 0.5 to 5, this code uses the rational Chebyshev approximation
 * from Cody and Thacher.  This method is extraordinarily fast and very
 * accurate over its range (slightly better than Cephes for most values).  If
 * we had quad floats, we could use the 9-term polynomial.
 */
long double ld_riemann_zeta(long double x) {
  int i;

  if (x < 0)  croak("Invalid input to RiemannZeta:  x must be >= 0");
  if (x == 1) return INFINITY;

  if (x == (unsigned int)x) {
    int k = x - 2;
    if ((k >= 0) && (k < (int)NPRECALC_ZETA))
      return riemann_zeta_table[k];
  }

  /* Cody / Thacher rational Chebyshev approximation for small values */
  if (x >= 0.5 && x <= 5.0) {
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
    long double sum = (sumn - (x-1)*sumd) / ((x-1)*sumd);
    return sum;
  }

  if (x > 17000.0)
    return 0.0;

#if 0
  {
    KAHAN_INIT(sum);
    /* Simple defining series, works well. */
    for (i = 5; i <= 1000000; i++) {
      long double term = powl(i, -x);
      KAHAN_SUM(sum, term);
      if (term < LDBL_EPSILON*sum) break;
    }
    KAHAN_SUM(sum, powl(4, -x) );
    KAHAN_SUM(sum, powl(3, -x) );
    KAHAN_SUM(sum, powl(2, -x) );
    return sum;
  }
#endif

  /* The 2n!/B_2k series used by the Cephes library. */
  {
    /* gp/pari:
     *   for(i=1,13,printf("%.38g\n",(2*i)!/bernreal(2*i)))
     * MPU:
     *   use bignum;
     *   say +(factorial(2*$_)/bernreal(2*$_))->bround(38) for 1..13;
     */
    static const long double A[] = {
      12.0L,
     -720.0L,
      30240.0L,
     -1209600.0L,
      47900160.0L,
     -1892437580.3183791606367583212735166425L,
      74724249600.0L,
     -2950130727918.1642244954382084600497650L,
      116467828143500.67248729113000661089201L,
     -4597978722407472.6105457273596737891656L,
      181521054019435467.73425331153534235290L,
     -7166165256175667011.3346447367083352775L,
      282908877253042996618.18640556532523927L,
    };
    long double a, b, s, t;
    const long double w = 10.0;
    s = 0.0;
    b = 0.0;
    for (i = 2; i < 11; i++) {
      b = powl( i, -x );
      s += b;
      if (fabsl(b) < fabsl(LDBL_EPSILON * s))
        return s;
    }
    s = s + b*w/(x-1.0) - 0.5 * b;
    a = 1.0;
    for (i = 0; i < 13; i++) {
      long double k = 2*i;
      a *= x + k;
      b /= w;
      t = a*b/A[i];
      s = s + t;
      if (fabsl(t) < fabsl(LDBL_EPSILON * s))
        break;
      a *= x + k + 1.0;
      b /= w;
    }
    return s;
  }
}

long double _XS_RiemannR(long double x) {
  long double part_term, term, flogx, ki;
  unsigned int k;
  KAHAN_INIT(sum);

  if (x <= 0) croak("Invalid input to ReimannR:  x must be > 0");

  if (x > 1e19) {
    const signed char* amob = _moebius_range(0, 100);
    KAHAN_SUM(sum, _XS_ExponentialIntegral(logl(x)));
    for (k = 2; k <= 100; k++) {
      if (amob[k] == 0) continue;
      ki = 1.0L / (long double) k;
      part_term = powl(x,ki);
      if (part_term > LDBL_MAX) return INFINITY;
      term = amob[k] * ki * _XS_ExponentialIntegral(logl(part_term));
      KAHAN_SUM(sum, term);
      if (fabsl(term) < fabsl(LDBL_EPSILON*sum)) break;
    }
    Safefree(amob);
    return sum;
  }

  KAHAN_SUM(sum, 1.0);

  flogx = logl(x);
  part_term = 1;

  for (k = 1; k <= 10000; k++) {
    part_term *= flogx / k;
    if (k-1 < NPRECALC_ZETA)  term = part_term / (k+k*riemann_zeta_table[k-1]);
    else                      term = part_term / (k+k*ld_riemann_zeta(k+1));
    KAHAN_SUM(sum, term);
    /* printf("R %5d after adding %.18Lg, sum = %.19Lg\n", k, term, sum); */
    if (fabsl(term) < fabsl(LDBL_EPSILON*sum)) break;
  }

  return sum;
}

long double lambertw(long double k) {
  long double x, lastx;
  int i;

  if (k < -0.36787944118L)
    croak("Invalid input to LambertW:  k must be >= -1/e");
  /* Make first estimate */
  if (k > 1) {
    long double lk = logl(k);
    long double llk = logl(lk);
    x = lk - llk - logl(1 - llk/lk)/2;
  } else {
    x = 0.567 * k;
  }
  lastx = x;
  for (i = 0; i < 100; i++) {   /* Use Halley's method */
    long double ex = expl(x);
    long double xex = x * ex;
    long double xexk = xex - k;
    long double x1 = x + 1;
    x = x - xexk / (ex * x1 - (x+2) * xexk/(2*x1));
    /* x = x - ( (x*ex-k) / (x*ex+ex-((x+2)*(x*ex-k)/(2*x+2))) ); */
    if (fabsl(lastx-x) < fabsl(LDBL_EPSILON)) break;
    lastx = x;
  }
  return x;
}

/* 1. Perform signed integer validation on b/blen.
 * 2. Compare to a/alen using min or max based on first arg.
 * 3. Return 0 to select a, 1 to select b.
 */
int strnum_minmax(int min, char* a, STRLEN alen, char* b, STRLEN blen)
{
  int aneg, bneg;
  STRLEN i;
  /* a is checked, process b */
  if (b == 0 || blen == 0) croak("Parameter must be a positive integer");
  bneg = (b[0] == '-');
  if (b[0] == '-' || b[0] == '+') { b++; blen--; }
  while (blen > 0 && *b == '0') { b++; blen--; }
  for (i = 0; i < blen; i++)
    if (!isDIGIT(b[i]))
      break;
  if (blen == 0 || i < blen)
    croak("Parameter must be a positive integer");

  if (a == 0) return 1;

  aneg = (a[0] == '-');
  if (a[0] == '-' || a[0] == '+') { a++; alen--; }
  while (alen > 0 && *a == '0') { a++; alen--; }

  if (aneg != bneg)  return  min  ?  (bneg == 1)  :  (aneg == 1);
  if (aneg == 1)  min = !min;
  if (alen != blen)  return  min  ?  (alen > blen) :  (blen > alen);

  for (i = 0; i < blen; i++)
    if (a[i] != b[i])
      return  min  ?  (a[i] > b[i])  :  (b[i] > a[i]);
  return 0; /* equal */
}
