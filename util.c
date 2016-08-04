#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
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
  extern long double sqrtl(long double);
#endif
#else
  #define powl(x, y)  (long double) pow( (double) (x), (double) (y) )
  #define expl(x)     (long double) exp( (double) (x) )
  #define logl(x)     (long double) log( (double) (x) )
  #define fabsl(x)    (long double) fabs( (double) (x) )
  #define floorl(x)   (long double) floor( (double) (x) )
  #define ceill(x)    (long double) ceil( (double) (x) )
  #define sqrtl(x)    (long double) sqrt( (double) (x) )
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
#define FUNC_is_perfect_cube
#define FUNC_is_perfect_fifth
#define FUNC_is_perfect_seventh
#define FUNC_next_prime_in_sieve 1
#define FUNC_prev_prime_in_sieve 1
#define FUNC_is_prime_in_sieve 1
#define FUNC_ipow 1
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
    next = next_prime_in_sieve(prime_sieve30, n, 30*NPRIME_SIEVE30);
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
    n -= wheelretreat30[m];
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
  long double fn, fl1, fl2, lower, a;

  if (n < 33000) return _XS_prime_count(2, n);

  fn  = (long double) n;
  fl1 = logl(n);
  fl2 = fl1 * fl1;

  if (n < UVCONST(52600000)) {   /* Dusart 2010, tweaked */
    a = (n <    88783) ? 1.89L /* n >= 33000 */
      : (n <   176000) ? 1.99L
      : (n <   315000) ? 2.11L
      : (n <  1100000) ? 2.19L
      : (n <  4500000) ? 2.31L : 2.35;
    lower = fn/fl1 * (1.0L + 1.0L/fl1 + a/fl2);
  } else if (fn < 1e19) {          /* Büthe 2015 1.9      1511.02032v1.pdf */
    lower = _XS_LogarithmicIntegral(fn) - (sqrtl(fn)/fl1) * (1.94L + 3.88L/fl1 + 27.57L/fl2);
  } else {                         /* Büthe 2014 v3 7.2   1410.7015v3.pdf */
    lower = _XS_LogarithmicIntegral(fn) - fl1*sqrtl(fn)/25.132741228718345907701147L;
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

  if (n < 33000) return _XS_prime_count(2, n);

  fn  = (long double) n;
  fl1 = logl(n);
  fl2 = fl1 * fl1;

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
    upper = _XS_LogarithmicIntegral(fn) - a * fl1*sqrtl(fn)/25.132741228718345907701147L;
  } else {                       /* Büthe 2014 7.4 */
    upper = _XS_LogarithmicIntegral(fn) + fl1*sqrtl(fn)/25.132741228718345907701147L;
  }
  return (UV) floorl(upper);
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

  if (n < 688383) {
    UV lo = fn * (flogn + flog2n - 1.0 + ((flog2n-2.10)/flogn));
    UV hi = fn * (flogn + flog2n);
    if (hi < lo) hi = MPU_MAX_PRIME;
    while (lo < hi) {
      UV mid = lo + (hi-lo)/2;
      if (prime_count_lower(mid) < n) lo = mid+1;
      else hi = mid;
    }
    return lo;
  }

  /* Dusart 2010 page 2 */
  upper = fn * (flogn + flog2n - 1.0 + ((flog2n-2.00)/flogn));
  /* Axler 2013 page viii Korollar G */
  if (n >= 8009824) upper -= fn * ((flog2n*flog2n-6*flog2n+10.273)/(2*flogn*flogn));

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
    /* Use Dusart 2010 page 2 bounds for the endpoints */
    UV lo = fn * (flogn + flog2n - 1.0 + ((flog2n-2.10)/flogn));
    UV hi = (n < 5000)
          ?  lo + 300
          :  fn * (flogn + flog2n - 1.0 + ((flog2n-1.95)/flogn));
    if (hi < lo) hi = MPU_MAX_PRIME;
    while (lo < hi) {
      UV mid = lo + (hi-lo)/2;
      if (prime_count_upper(mid) < n) lo = mid+1;
      else                            hi = mid;
    }
    return lo;
  }

  /* Axler 2013 page vii Korollar I, for all n >= 2 */
  lower = fn * (flogn + flog2n-1.0 + ((flog2n-2.00)/flogn) - ((flog2n*flog2n-6*flog2n+11.847)/(2*flogn*flogn)));
  return (UV) ceill(lower);
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
      UV p, bytes = seg_high/30 - seg_low/30 + 1;
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

/******************************************************************************/
/*                             RAMANUJAN PRIMES                               */
/******************************************************************************/

/* For Ramanujan prime estimates:
 *  - counts are done via inverse nth, so only one thing to tune.
 *  - For nth tables, upper values ok if too high, lower values ok if too low.
 *  - both upper & lower empirically tested to 175e9 (175 thousand million),
 *    with a return value of over 10^13.
 */

/* These are playing loose with Sondow/Nicholson/Noe 2011 theorem 4.
 * The last value should be rigorously checked using actual R_n values. */
static const UV ram_upper_idx[] = {
  5215, 5223, 5261, 5271, 5553, 7431, 7451, 8582, 12589, 12620, 12762,
  18154, 18294, 18410, 25799, 28713, 40061, 45338, 63039, 65724,
  88726, 107849, 151742, 216978, 270700, 347223, 508096,
  768276, 1090285, 1568165, 2375799, 4162908, 6522443, 11406250,
  20637716, 39711166, 80161468, 174200145, 404761877, 1024431762,
  UVCONST(2868095461),
#if BITS_PER_WORD == 64  /* 1383,1382,1381, 1380,1379 */
  UVCONST(   9136430799), UVCONST(  33244053524), UVCONST( 143852101796),
  UVCONST( 760145301247), UVCONST(5136852322733)
#else
  UVCONST( 4294967295)
#endif
};
#define NRAM_UPPER_MULT 1425
#define NRAM_UPPER (sizeof(ram_upper_idx)/sizeof(ram_upper_idx[0]))

UV nth_ramanujan_prime_upper(UV n) {
  UV i, mult, res;
  if (n <= 2) return (n==0) ? 0 : (n==1) ? 2 : 11;
  /* While p_3n is a complete upper bound, Rp_n tends to p_2n, and
   * SNN(2011) theorem 4 shows how we can find (m,c) values where m < 1,
   * Rn < m*p_3n for all n > c.  Here we use various quantized m values
   * and the table gives us c values where it applies. */
  if      (n < 20) mult = 1787;
  else if (n < 98) mult = 1670;
  else if (n < 1580) mult = 1520;
  else if (n < 5214) mult = 1440;
  else {
    for (i = 0; i < NRAM_UPPER; i++)
      if (ram_upper_idx[i] > n)
        break;
    mult = NRAM_UPPER_MULT-i;
  }
  res = nth_prime_upper(3*n);
  if (res > (UV_MAX/mult)) res = (UV) (((long double) mult / 2048.0L) * res);
  else                     res = (res * mult) >> 11;
  return res;
}

static const UV ram_lower_idx[] = {
  5935, 6013, 6107, 8726, 8797, 9396, 9556, 9611, 13314, 13405, 13641,
  18655, 20457, 22665, 23744, 31015, 34432, 43298, 50564, 54168, 69194,
  75658, 97434, 110360, 149399, 164361, 224590, 265154, 337116, 386011,
  514260, 606067, 804041, 981748, 1367781, 1874883, 2448102, 3248825,
  4513600, 6240819, 8581572, 12175878, 17555268, 25999335, 39295429,
  60160729, 96043141, 153729257, 255707435, 438238581, 771748027,
  UVCONST(1413859126), UVCONST(2691176863),
#if BITS_PER_WORD == 64 /* 1055,1054,1053, 1052,1051,1050, 1049,1048 */
  UVCONST(   5351635422), UVCONST(  11116948392), UVCONST(  24348237409),
  UVCONST(  56422464192), UVCONST( 139376328492), UVCONST( 369783640662),
  UVCONST(1061320360069), UVCONST(3332007446416),
  /* These are low-biased estimates */
  UVCONST(3332007446416) *   3, UVCONST(3332007446416) *   9,
  UVCONST(3332007446416) *  27, UVCONST(3332007446416) *  81
#else
  UVCONST(   4294967295)
#endif
};
#define NRAM_LOWER_MULT 1108
#define NRAM_LOWER (sizeof(ram_lower_idx)/sizeof(ram_lower_idx[0]))

UV nth_ramanujan_prime_lower(UV n) {
  UV res, i, mult;
  if (n <= 2) return (n==0) ? 0 : (n==1) ? 2 : 11;

  res = nth_prime_lower(2*n);

  if (n < ram_lower_idx[NRAM_LOWER-1]) {
    for (i = 0; i < NRAM_LOWER; i++)
      if (ram_lower_idx[i] > n)
        break;
    mult = (NRAM_LOWER_MULT-i);
    if (res > (UV_MAX/mult)) res = (UV) (((long double) mult / 1024.0L) * res);
    else                     res = (res * mult) >> 10;
  }
  return res;
}

/* An advantage of making these binary searches on the inverse is that we
 * don't have to tune them separately, and nothing changes if the prime
 * count bounds are modified.  We do need to keep up to date with any
 * changes to nth_prime_{lower,upper} however. */

UV ramanujan_prime_count_lower(UV n) {
  UV lo, hi;
  if (n < 29) return (n < 2) ? 0 : (n < 11) ? 1 : (n < 17) ? 2 : 3;
  /* Binary search on nth_ramanujan_prime_upper */
  /* We know we're between p_2n and p_3n, probably close to the former. */
  lo = prime_count_lower(n)/3;
  hi = prime_count_upper(n) >> 1;
  while (lo < hi) {
    UV mid = lo + (hi-lo)/2;
    if (nth_ramanujan_prime_upper(mid) < n) lo = mid+1;
    else                                    hi = mid;
  }
  return lo-1;
}
UV ramanujan_prime_count_upper(UV n) {
  /* return prime_count_upper(n) >> 1; */       /* Simple bound */
  UV lo, hi;
  if (n < 29) return (n < 2) ? 0 : (n < 11) ? 1 : (n < 17) ? 2 : 3;
  /* Binary search on nth_ramanujan_prime_upper */
  /* We know we're between p_2n and p_3n, probably close to the former. */
  lo = prime_count_lower(n)/3;
  hi = prime_count_upper(n) >> 1;
  while (lo < hi) {
    UV mid = lo + (hi-lo)/2;
    if (nth_ramanujan_prime_lower(mid) < n) lo = mid+1;
    else                                    hi = mid;
  }
  return lo-1;
}

/* Return array of first n ramanujan primes.  Use Noe's algorithm. */
UV* n_ramanujan_primes(UV n) {
  UV max, k, s, *L;
  unsigned char* sieve;
  max = nth_ramanujan_prime_upper(n); /* Rn <= max, so we can sieve to there */
  if (_XS_get_verbose() >= 2) printf("sieving to %"UVuf" for first %"UVuf" Ramanujan primes\n", max, n);
  Newz(0, L, n, UV);
  L[0] = 2;
  sieve = sieve_erat30(max);
  for (s = 0, k = 7; k <= max; k += 2) {
    if (is_prime_in_sieve(sieve, k)) s++;
    if (s < n) L[s] = k+1;
    if ((k & 3) == 1 && is_prime_in_sieve(sieve, (k+1)>>1)) s--;
    if (s < n) L[s] = k+2;
  }
  Safefree(sieve);
  return L;
}

UV* n_range_ramanujan_primes(UV nlo, UV nhi) {
  UV mink, maxk, k, s, *L;
  UV seg1beg, seg1end, seg2beg, seg2end;
  unsigned char *seg1, *seg2;
  int verbose = _XS_get_verbose();

  if (nlo == 0) nlo = 1;
  if (nhi == 0) nhi = 1;

  /* If we're starting from 1, just do single monolithic sieve */
  if (nlo == 1)  return n_ramanujan_primes(nhi);

  Newz(0, L, nhi-nlo+1, UV);
  if (nlo <= 1 && nhi >= 1) L[1-nlo] =  2;
  if (nlo <= 2 && nhi >= 2) L[2-nlo] = 11;
  if (nhi < 3) return L;

  mink = nth_ramanujan_prime_lower(nlo);
  maxk = nth_ramanujan_prime_upper(nhi);
  if (mink < 15) mink = 15;
  if (mink % 2 == 0) mink--;
  s = 1 + _XS_LMO_pi(mink-2) - _XS_LMO_pi((mink-1)>>1);
  if (verbose >= 2) printf("Generate Rn[%"UVuf"] to Rn[%"UVuf"]: search %"UVuf" to %"UVuf"\n", nlo, nhi, mink, maxk);

  seg1beg = 30 * (mink/30);
  seg1end = 30 * ((maxk+29)/30);
  New(0, seg1, (seg1end - seg1beg)/30 + 1, unsigned char);
  if (verbose >= 2) printf("sieving %"UVuf" to %"UVuf"\n", seg1beg, seg1end);
  (void) sieve_segment(seg1, seg1beg/30, seg1end/30);

  seg2beg = 30 * (((mink+1)>>1)/30);
  seg2end = 30 * ((((maxk+1)>>1)+29)/30);
  New(0, seg2, (seg2end - seg2beg)/30 + 1, unsigned char);
  if (verbose >= 2) printf("sieving %"UVuf" to %"UVuf"\n", seg2beg, seg2end);
  (void) sieve_segment(seg2, seg2beg/30, seg2end/30);

  if (verbose >= 2) printf("Noe's algorithm %"UVuf" to %"UVuf"\n", mink, maxk);
  for (k = mink; k <= maxk; k += 2) {
    if (is_prime_in_sieve(seg1, k-seg1beg)) s++;
    if (s >= nlo && s <= nhi) L[s-nlo] = k+1;
    if ((k & 3) == 1 && is_prime_in_sieve(seg2, ((k+1)>>1)-seg2beg)) s--;
    if (s >= nlo && s <= nhi) L[s-nlo] = k+2;
  }
  Safefree(seg1);
  Safefree(seg2);
  if (verbose >= 2) printf("Generated %"UVuf" Ramanujan primes from %"UVuf" to %"UVuf"\n", nhi-nlo+1, L[0], L[nhi-nlo]);
  return L;
}

UV nth_ramanujan_prime(UV n) {
  UV rn, *L;
  if (n <= 2) return (n == 0) ? 0 : (n == 1) ? 2 : 11;
  L = n_range_ramanujan_primes(n, n);
  rn = L[0];
  Safefree(L);
  return rn;
}

/* Returns array of Ram primes between low and high, results from first->last */
UV* ramanujan_primes(UV* first, UV* last, UV low, UV high)
{
  UV beg, end, range, nlo, nhi, *L;

  if (high < 2 || high < low) return 0;
  if (low < 2) low = 2;

  nlo = ramanujan_prime_count_lower(low);
  nhi = ramanujan_prime_count_upper(high);
  L = n_range_ramanujan_primes(nlo, nhi);

  /* Skip out of range values at the start and end */
  range = nhi-nlo; beg = 0; end = range;

  while (beg+10000 <= range && L[beg+10000] <  low) beg += 10000;
  while (beg+100   <= range && L[beg+100  ] <  low) beg += 100;
  while (beg       <= range && L[beg      ] <  low) beg++;

  while (beg+10000 <= end   && L[end-10000] > high) end -= 10000;
  while (beg+100   <= end   && L[end-100  ] > high) end -= 100;
  while (beg       <= end   && L[end      ] > high) end--;

  *first = beg;
  *last = end;
  return L;
}

int is_ramanujan_prime(UV n) {
  UV beg, end, *L;

  if (!_XS_is_prime(n))  return 0;
  if (n < 17)            return (n == 2 || n == 11);

  /* Generate Ramanujan primes and see if we're in the list.  Slow. */
  L = ramanujan_primes(&beg, &end, n, n);
  Safefree(L);
  return (beg <= end);
}

UV ramanujan_prime_count_approx(UV n)
{
  /* TODO: overflow */
  /* We should be able to do better, but this isn't too bad. */
  /* (11*lower + 5*upper) / 16 is a bit closer */
  return (3*ramanujan_prime_count_lower(n) + 1*ramanujan_prime_count_upper(n)) >> 2;
}

#if BITS_PER_WORD == 64
#define RAMPC2 47
static const UV ramanujan_counts_pow2[RAMPC2+1] = {
   0, 1, 1, 1, 2, 4, 7, 13, 23, 42, 75, 137, 255, 463, 872, 1612,
   3030, 5706, 10749, 20387, 38635, 73584, 140336, 268216, 513705,
   985818, 1894120, 3645744, 7027290, 13561906, 26207278, 50697533,
   98182656, 190335585, 369323301, 717267167,
   UVCONST(   1394192236), UVCONST(   2712103833), UVCONST(   5279763823),
   UVCONST(  10285641777), UVCONST(  20051180846), UVCONST(  39113482639),
   UVCONST(  76344462797), UVCONST( 149100679004), UVCONST( 291354668495),
   UVCONST( 569630404447), UVCONST(1114251967767), UVCONST(2180634225768) };
#else
#define RAMPC2 31  /* input limited */
static const UV ramanujan_counts_pow2[RAMPC2+1] = {
   0, 1, 1, 1, 2, 4, 7, 13, 23, 42, 75, 137, 255, 463, 872, 1612,
   3030, 5706, 10749, 20387, 38635, 73584, 140336, 268216, 513705,
   985818, 1894120, 3645744, 7027290, 13561906, 26207278, 50697533 };
#endif
UV ramanujan_prime_count(UV lo, UV hi)
{
  UV count = 0, beg, end, inc, log2, *L;
  int addcount = 1;

  if (hi < 2 || hi < lo) return 0;

  if (lo <= 2) {
    log2 = log2floor(hi);
    /* We have some perfect powers of 2 in our table */
    if ((hi & (hi-1)) == 0 && log2 <= RAMPC2)
      return ramanujan_counts_pow2[log2];

    if (hi > 22000) {   /* Faster for small values to go direct */
      /* We can calculate Rp(n) fairly efficiently and with relatively low
       * memory if there are tight bounds available.  This lets us quickly
       * skip forward without having to calculate them all.  We'll use an
       * approximation to get close, and be prepared to go backwards. */
      UV rpcl, nthl;
      int verbose = _XS_get_verbose();
      rpcl = ramanujan_prime_count_approx(hi);
      if (verbose >= 2) printf("Calculating the %"UVuf"-th Ramanujan prime as an approximation\n", rpcl);
      nthl = nth_ramanujan_prime(rpcl);
      count = rpcl;
      if (nthl <= hi) {       /* Count forwards */
        lo = nthl+1;
      } else {                /* Count backwards */
        lo = hi+1;
        hi = nthl;
        addcount = 0;
      }
      if (verbose >= 2) printf("counting %s from %"UVuf"\n", addcount ? "forwards" : "backwards", count );
    } else {
      /* Bump up the lower limit based on our table */
      if (log2 > RAMPC2) log2 = RAMPC2;
      count = ramanujan_counts_pow2[log2];
      lo = (UVCONST(1) << log2) + 1;
    }
  }

  /* Generate all Rp from lo to hi */
  L = ramanujan_primes(&beg, &end, lo, hi);
  inc = (L && end >= beg) ? end-beg+1 : 0;
  count = addcount ? count+inc : count-inc;
  Safefree(L);
  return count;
}

/******************************************************************************/
/*                            SUMS AND PRINTING                               */
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
static int my_sprint(char* ptr, UV val) {
  int nchars;
  UV t;
  char *s = ptr;
  do {
    t = val / 10; *s++ = (char) ('0' + val - 10 * t);
  } while ((val = t));
  nchars = s - ptr + 1;  *s = '\n';
  while (--s > ptr) { char c = *s; *s = *ptr; *ptr++ = c; }
  return nchars;
}
static char* write_buf(int fd, char* buf, char* bend) {
  int res = (int) write(fd, buf, bend-buf);
  if (res == -1) croak("print_primes write error");
  return buf;
}
void print_primes(UV low, UV high, int fd) {
  char buf[8000+25];
  char* bend = buf;
  if ((low <= 2) && (high >= 2)) bend += my_sprint(bend,2);
  if ((low <= 3) && (high >= 3)) bend += my_sprint(bend,3);
  if ((low <= 5) && (high >= 5)) bend += my_sprint(bend,5);
  if (low < 7) low = 7;

  if (low <= high) {
    unsigned char* segment;
    UV seg_base, seg_low, seg_high;
    void* ctx = start_segment_primes(low, high, &segment);
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
        bend += my_sprint(bend,p);
        if (bend-buf > 8000) { bend = write_buf(fd, buf, bend); }
      END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
  }
  if (bend > buf) { bend = write_buf(fd, buf, bend); }
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
  nextlog = UVCONST(2) << logp;
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
    START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high ) {
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
  UV u, j, m, nmk, maxmu;
  signed char* mu;
  short* M;   /* 16 bits is enough range for all 32-bit M => 64-bit n */
  IV sum;

  if (n <= 1)  return n;
  u = isqrt(n);
  maxmu = (n/(u+1));              /* maxmu lets us handle u < sqrt(n) */
  if (maxmu < u) maxmu = u;
  mu = _moebius_range(0, maxmu);
  New(0, M, maxmu+1, short);
  M[0] = 0;
  for (j = 1; j <= maxmu; j++)
    M[j] = M[j-1] + mu[j];
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
      sum += (mu[m] > 0) ? -inner_sum : inner_sum;
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
  UV t;
  if ((n <= 3) || (n == UV_MAX)) return 1;
  if ((n & (n-1)) == 0)          return ctz(n);  /* powers of 2    */
  if (is_perfect_square(n))      return 2 * powerof(isqrt(n));
  if (is_perfect_cube(n))        return 3 * powerof(icbrt(n));

  /* Simple rejection filter for non-powers of 5-37.  Rejects 47.85%. */
  t = n & 511; if ((t*77855451) & (t*4598053) & 862)  return 1;

  if (is_perfect_fifth(n))       return 5 * powerof(rootof(n,5));
  if (is_perfect_seventh(n))     return 7 * powerof(rootof(n,7));

  if (n > 177146 && n <= UVCONST(1977326743)) {
    switch (n) { /* Check for powers of 11, 13, 17, 19 within 32 bits */
      case 177147: case 48828125: case 362797056: case 1977326743: return 11;
      case 1594323: case 1220703125: return 13;
      case 129140163: return 17;
      case 1162261467: return 19;
      default:  break;
    }
  }
#if BITS_PER_WORD == 64
  if (n >= UVCONST(8589934592)) {
    /* The Bloom filters reject about 90% of inputs each, about 99% for two.
     * Bach/Sorenson type sieves do about as well, but are much slower due
     * to using a powmod. */
    if ( (t = n %121, !((t*19706187) & (t*61524433) & 876897796)) &&
         (t = n % 89, !((t*28913398) & (t*69888189) & 2705511937U)) ) {
      /* (t = n % 67, !((t*117621317) & (t*48719734) & 537242019)) ) { */
      UV root = rootof(n,11);
      if (n == ipow(root,11)) return 11;
    }
    if ( (t = n %131, !((t*1545928325) & (t*1355660813) & 2771533888U)) &&
         (t = n % 79, !((t*48902028) & (t*48589927) & 404082779)) ) {
      /* (t = n % 53, !((t*79918293) & (t*236846524) & 694943819)) ) { */
      UV root = rootof(n,13);
      if (n == ipow(root,13)) return 13;
    }
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
      return !is_perfect_cube(n) ? 0 : (a == 3) ? 1 : is_power(icbrt(n),a/3);
    if ((a % 5) == 0)
      return !is_perfect_fifth(n) ? 0 : (a == 5) ? 1 :is_power(rootof(n,5),a/5);
  }
  ret = powerof(n);
  if (a != 0) return !(ret % a);  /* Is the max power divisible by a? */
  return (ret == 1) ? 0 : ret;
}

#if BITS_PER_WORD == 64
#define ROOT_MAX_3 41
static const uint32_t root_max[ROOT_MAX_3] = {0,0,0,2642245,65535,7131,1625,565,255,138,84,56,40,30,23,19,15,13,11,10,9,8,7,6,6,5,5,5,4,4,4,4,3,3,3,3,3,3,3,3,3};
#else
#define ROOT_MAX_3 21
static const uint32_t root_max[ROOT_MAX_3] = {0,0,0,1625,255,84,40,23,15,11,9,7,6,5,4,4,3,3,3,3,3};
#endif

UV rootof(UV n, UV k) {
  UV lo, hi, max;
  if (k == 0) return 0;
  if (k == 1) return n;
  if (k == 2) return isqrt(n);
  if (k == 3) return icbrt(n);

  /* Bracket between powers of 2, but never exceed max power so ipow works */
  max = 1 + ((k >= ROOT_MAX_3) ? 2 : root_max[k]);
  lo = UVCONST(1) << (log2floor(n)/k);
  hi = ((lo*2) < max) ? lo*2 : max;

  /* Binary search */
  while (lo < hi) {
    UV mid = lo + (hi-lo)/2;
    if (ipow(mid,k) <= n) lo = mid+1;
    else                  hi = mid;
  }
  return lo-1;
}

int primepower(UV n, UV* prime)
{
  int power = 0;
  if (n < 2) return 0;
  /* Check for small divisors */
  if (!(n&1)) {
    if (n & (n-1)) return 0;
    *prime = 2;
    return ctz(n);
  }
  if ((n%3) == 0) {
    do { n /= 3; power++; } while (n > 1 && (n%3) == 0);
    if (n != 1) return 0;
    *prime = 3;
    return power;
  }
  if ((n%5) == 0) {
    do { n /= 5; power++; } while (n > 1 && (n%5) == 0);
    if (n != 1) return 0;
    *prime = 5;
    return power;
  }
  if ((n%7) == 0) {
    do { n /= 7; power++; } while (n > 1 && (n%7) == 0);
    if (n != 1) return 0;
    *prime = 7;
    return power;
  }
  if (is_prob_prime(n))
    { *prime = n; return 1; }
  /* Composite.  Test for perfect power with prime root. */
  power = powerof(n);
  if (power == 1) power = 0;
  if (power) {
    UV root = rootof(n, (UV)power);
    if (is_prob_prime(root))
      *prime = root;
    else
      power = 0;
  }
  return power;
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

UV logint(UV n, UV b)
{
  /* UV e;  for (e=0; n; n /= b) e++;  return e-1; */
  UV v, e = 0;
  if (b == 2)
    return log2floor(n);
  if (n > UV_MAX/b) {
    n /= b;
    e = 1;
  }
  for (v = b; v <= n; v *= b)
    e++;
  return e;
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
  UV totient;
  if (k == 0 || n <= 1) return (n == 1);
  if (k > 6 || (k > 1 && n >= jordan_overflow[k-2])) return 0;

  totient = 1;
  /* Similar to Euler totient, shortcut even inputs */
  while ((n & 0x3) == 0) { n >>= 1; totient *= (1<<k); }
  if ((n & 0x1) == 0) { n >>= 1; totient *= ((1<<k)-1); }
  nfac = factor(n,factors);
  for (i = 0; i < nfac; i++) {
    UV p = factors[i];
    UV pk = ipow(p,k);
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

int is_carmichael(UV n) {
  UV fac[MPU_MAX_FACTORS+1];
  UV exp[MPU_MAX_FACTORS+1];
  int i, nfactors;

  /* Small or even is not a Carmichael number */
  if (n < 561 || !(n&1)) return 0;

  /* Simple pre-test for square free */
  if (!(n% 4) || !(n% 9) || !(n%25) || !(n%49) || !(n%121) || !(n%169))
    return 0;

  /* Fast check without having to factor */
  if (n > 5000000 && powmod(2,n-1,n) != 1) return 0;

  nfactors = factor_exp(n, fac, exp);
  if (nfactors < 3)
    return 0;
  for (i = 0; i < nfactors; i++) {
    if (exp[i] > 1  ||  ((n-1) % (fac[i]-1)) != 0)
      return 0;
  }
  return 1;
}

static int is_quasi_base(int nfactors, UV *fac, UV p, UV b) {
  int i;
  for (i = 0; i < nfactors; i++) {
    UV d = fac[i] - b;
    if (d == 0 || (p % d) != 0)
      return 0;
  }
  return 1;
}

/* Returns number of bases that pass */
UV is_quasi_carmichael(UV n) {
  UV nbases, fac[MPU_MAX_FACTORS+1], exp[MPU_MAX_FACTORS+1];
  UV spf, lpf, ndivisors, *divs;
  int i, nfactors;

  if (n < 35) return 0;

  /* Simple pre-test for square free */
  if (!(n% 4) || !(n% 9) || !(n%25) || !(n%49) || !(n%121) || !(n%169))
    return 0;

  nfactors = factor_exp(n, fac, exp);
  /* Must be composite */
  if (nfactors < 2)
    return 0;
  /* Must be square free */
  for (i = 0; i < nfactors; i++)
    if (exp[i] > 1)
      return 0;

  nbases = 0;
  spf = fac[0];
  lpf = fac[nfactors-1];

  /* Algorithm from Hiroaki Yamanouchi, 2015 */
  if (nfactors == 2) {
    divs = _divisor_list(n / spf - 1, &ndivisors);
    for (i = 0; i < (int)ndivisors; i++) {
      UV d = divs[i];
      UV k = spf - d;
      if (d >= spf) break;
      if (is_quasi_base(nfactors, fac, n-k, k))
        nbases++;
    }
  } else {
    divs = _divisor_list(lpf * (n / lpf - 1), &ndivisors);
    for (i = 0; i < (int)ndivisors; i++) {
      UV d = divs[i];
      UV k = lpf - d;
      if (lpf > d && k >= spf) continue;
      if (k != 0 && is_quasi_base(nfactors, fac, n-k, k))
        nbases++;
    }
  }
  Safefree(divs);
  return nbases;
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
  UV p;
  if (!primepower(n,&p)) return 1;     /* Not a prime power */
  return p;
}


UV znorder(UV a, UV n) {
  UV fac[MPU_MAX_FACTORS+1];
  UV exp[MPU_MAX_FACTORS+1];
  int i, nfactors;
  UV k, phi;

  if (n <= 1) return n;   /* znorder(x,0) = 0, znorder(x,1) = 1          */
  if (a <= 1) return a;   /* znorder(0,x) = 0, znorder(1,x) = 1  (x > 1) */
  if (gcd_ui(a,n) > 1)  return 0;

  /* Cohen 1.4.3 using Carmichael Lambda */
  phi = carmichael_lambda(n);
  nfactors = factor_exp(phi, fac, exp);
  k = phi;
  for (i = 0; i < nfactors; i++) {
    UV b, a1, ek, pi = fac[i], ei = exp[i];
    b = ipow(pi,ei);
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
  UV a, phi, on, r;
  int i, nfactors;

  if (n <= 4) return (n == 0) ? 0 : n-1;
  if (n % 4 == 0)  return 0;

  on = (n&1) ? n : (n>>1);
  a = powerof(on);
  r = rootof(on, a);
  if (!is_prob_prime(r)) return 0;        /* c^a or 2c^a */
  phi = (r-1) * (on/r);                   /* p^a or 2p^a */

  nfactors = factor_exp(phi, fac, exp);
  for (i = 0; i < nfactors; i++)
    exp[i] = phi / fac[i];  /* exp[i] = phi(n) / i-th-factor-of-phi(n) */
  for (a = 2; a < n; a++) {
    /* Skip first few perfect powers */
    if (a == 4 || a == 8 || a == 9) continue;
    /* Skip values we know can't be right: (a|n) = 0 (or 1 for odd primes) */
    if (phi == n-1) {
      if (kronecker_uu(a, n) != -1)  continue;
    } else {
      if (kronecker_uu(a, n) == 0)  continue;
    }
    for (i = 0; i < nfactors; i++)
      if (powmod(a, exp[i], n) == 1)
        break;
    if (i == nfactors) return a;
  }
  return 0;
}

int is_primitive_root(UV a, UV n, int nprime) {
  UV s, fac[MPU_MAX_FACTORS+1];
  int i, nfacs;
  if (n <= 1) return n;
  a %= n;
  if (gcd_ui(a,n) != 1) return 0;
  s = nprime ? n-1 : totient(n);
  if ((s % 2) == 0 && powmod(a, s/2, n) == 1) return 0;
  if ((s % 3) == 0 && powmod(a, s/3, n) == 1) return 0;
  if ((s % 5) == 0 && powmod(a, s/5, n) == 1) return 0;
  nfacs = factor_exp(s, fac, 0);
  for (i = 0; i < nfacs; i++) {
    if (fac[i] > 5 && powmod(a, s/fac[i], n) == 1) return 0;
  }
  return 1;
}

IV gcdext(IV a, IV b, IV* u, IV* v, IV* cs, IV* ct) {
  IV s = 0;  IV os = 1;
  IV t = 1;  IV ot = 0;
  IV r = b;  IV or = a;
  if (a == 0 && b == 0) { os = 0; t = 0; }
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

static int verify_sqrtmod(UV s, UV *rs, UV a, UV p) {
  if (p-s < s)  s = p-s;
  if (mulmod(s, s, p) != a) return 0;
  *rs = s;
  return 1;
}
int sqrtmod(UV *s, UV a, UV p) {
  if (p == 0) return 0;
  if (a >= p) a %= p;
  if (p <= 2 || a <= 1) return verify_sqrtmod(a, s,a,p);

  if ((p % 4) == 3) {
    return verify_sqrtmod(powmod(a, (p+1)>>2, p), s,a,p);
  }
  if ((p % 8) == 5) { /* Atkin's algorithm.  Faster than Legendre. */
    UV a2, alpha, beta, b;
    a2 = addmod(a,a,p);
    alpha = powmod(a2,(p-5)>>3,p);
    beta  = mulmod(a2,sqrmod(alpha,p),p);
    b     = mulmod(alpha, mulmod(a, (beta ? beta-1 : p-1), p), p);
    return verify_sqrtmod(b, s,a,p);
  }
  if ((p % 16) == 9) { /* Müller's algorithm extending Atkin */
    UV a2, alpha, beta, b, d = 1;
    a2 = addmod(a,a,p);
    alpha = powmod(a2, (p-9)>>4, p);
    beta  = mulmod(a2, sqrmod(alpha,p), p);
    if (sqrmod(beta,p) != p-1) {
      do { d += 2; } while (kronecker_uu(d,p) != -1 && d < p);
      alpha = mulmod(alpha, powmod(d,(p-9)>>3,p), p);
      beta  = mulmod(a2, mulmod(sqrmod(d,p),sqrmod(alpha,p),p), p);
    }
    b = mulmod(alpha, mulmod(a, mulmod(d,(beta ? beta-1 : p-1),p),p),p);
    return verify_sqrtmod(b, s,a,p);
  }

  /* Verify Euler condition for odd p */
  if ((p & 1) && powmod(a,(p-1)>>1,p) != 1) return 0;

  {
    UV x, q, e, t, z, r, m, b;
    q = p-1;
    e = valuation(q, 2);
    q >>= e;
    t = 3;
    while (kronecker_uu(t, p) != -1) {
      t += 2;
      if (t == 201) {           /* exit if p looks like a composite */
        if ((p % 2) == 0 || powmod(2, p-1, p) != 1 || powmod(3, p-1, p) != 1)
          return 0;
      } else if (t >= 20000) {  /* should never happen */
        return 0;
      }
    }
    z = powmod(t, q, p);
    r = e;
    b = powmod(a, q, p);
    q = (q+1) >> 1;
    x = powmod(a, q, p);
    while (b != 1) {
      t = b;
      m = 0;
      do {
        t = powmod(t, 2, p);
        m++;
      } while (m < r && t != 1);
      if (m >= r) break;
      t = powmod(z, UVCONST(1) << (r-m-1), p);
      x = mulmod(x, t, p);
      z = mulmod(t, t, p);
      b = mulmod(b, z, p);
      r = m;
    }
    return verify_sqrtmod(x, s,a,p);
  }
  return 0;
}

int sqrtmod_composite(UV *s, UV a, UV n) {
  UV fac[MPU_MAX_FACTORS+1];
  UV exp[MPU_MAX_FACTORS+1];
  UV sqr[MPU_MAX_FACTORS+1];
  UV p, j, k, gcdan;
  int i, nfactors;

  if (n == 0) return 0;
  if (a >= n) a %= n;
  if (n <= 2 || a <= 1) return verify_sqrtmod(a, s,a,n);

  /* Simple existence check.  It's still possible no solution exists.*/
  if (kronecker_uu(a, ((n%4) == 2) ? n/2 : n) == -1) return 0;

  /* if 8|n 'a' must = 1 mod 8, else if 4|n 'a' must = 1 mod 4 */
  if ((n % 4) == 0) {
    if ((n % 8) == 0) {
      if ((a % 8) != 1) return 0;
    } else {
      if ((a % 4) != 1) return 0;
    }
  }

  /* More detailed existence check before factoring.  Still possible. */
  gcdan = gcd_ui(a, n);
  if (gcdan == 1) {
    if ((n % 3) == 0 && kronecker_uu(a, 3) != 1) return 0;
    if ((n % 5) == 0 && kronecker_uu(a, 5) != 1) return 0;
    if ((n % 7) == 0 && kronecker_uu(a, 7) != 1) return 0;
  }

  /* Factor n */
  nfactors = factor_exp(n, fac, exp);

  /* If gcd(a,n)==1, this answers comclusively if a solution exists. */
  if (gcdan == 1) {
    for (i = 0; i < nfactors; i++)
      if (fac[i] > 7 && kronecker_uu(a, fac[i]) != 1) return 0;
  }

  for (i = 0; i < nfactors; i++) {

    /* Powers of 2 */
    if (fac[i] == 2) {
      if (exp[i] == 1) {
        sqr[i] = a & 1;
      } else if (exp[i] == 2) {
        sqr[i] = 1;  /* and 3 */
      } else {
        UV this_roots[256], next_roots[256];
        UV nthis = 0, nnext = 0;
        this_roots[nthis++] = 1;
        this_roots[nthis++] = 3;
        for (j = 2; j < exp[i]; j++) {
          p = UVCONST(1) << (j+1);
          nnext = 0;
          for (k = 0; k < nthis && nnext < 254; k++) {
            UV r = this_roots[k];
            if (sqrmod(r,p) == (a % p))
              next_roots[nnext++] = r;
            if (sqrmod(p-r,p) == (a % p))
              next_roots[nnext++] = p-r;
          }
          if (nnext == 0) return 0;
          /* copy next exponent's found roots to this one */
          nthis = nnext;
          for (k = 0; k < nnext; k++)
            this_roots[k] = next_roots[k];
        }
        sqr[i] = this_roots[0];
      }
      continue;
    }

    /* p is an odd prime */
    p = fac[i];
    if (!sqrtmod(&(sqr[i]), a, p))
      return 0;

    /* Lift solution of x^2 = a mod p  to  x^2 = a mod p^e */
    for (j = 1; j < exp[i]; j++) {
      UV xk2, yk, expect, sol;
      xk2 = addmod(sqr[i],sqr[i],p);
      yk = modinverse(xk2, p);
      expect = mulmod(xk2,yk,p);
      p *= fac[i];
      sol = submod(sqr[i], mulmod(submod(sqrmod(sqr[i],p), a % p, p), yk, p), p);
      if (expect != 1 || sqrmod(sol,p) != (a % p)) {
        /* printf("a %lu failure to lift to %lu^%d\n", a, fac[i], j+1); */
        return 0;
      }
      sqr[i] = sol;
    }
  }

  /* raise fac[i] */
  for (i = 0; i < nfactors; i++)
    fac[i] = ipow(fac[i], exp[i]);

  p = chinese(sqr, fac, nfactors, &i);
  return (i == 1) ? verify_sqrtmod(p, s, a, n) : 0;
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
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high ) {
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
  /* Kulsha's method.  Also not faster than our eint, though close */
  if (x > 1) {
    static long double a2[] = {1.000000000000000000L,1.333333333333333333L,1.533333333333333333L,1.676190476190476190L,1.787301587301587301L,1.878210678210678210L,1.955133755133755133L,2.021800421800421800L,2.080623951212186506L,2.133255530159554927L,2.180874577778602546L,2.224352838648167763L,2.264352838648167763L,2.301389875685204800L,2.335872634305894455L,2.368130698822023487L,2.398433729125053790L,2.427005157696482361L,2.454032184723509388L,2.479673210364535029L,2.504063454266974053L,2.527319268220462425L,2.549541490442684647L,2.570818086187365498L,2.591226249452671620L,2.610834092589926522L,2.629702017118228409L,2.647883835300046591L,2.665427694949169398L,2.682376847491542279L,2.698770290114493099L,2.714643305987508972L,2.730027921372124357L,2.744953294506452715L,2.759446048129641121L,2.773530555171894642L,2.787229185308880943L,2.800562518642214276L,2.813549531629227263L,2.826207759477328529L,2.838553438489674208L,2.850601631260758545L,2.862366337143111486L,2.873860590016674704L,2.885096545072854479L,2.896085556061865468L,2.906838244233908479L,2.917364560023382163L,2.927673838373897627L,2.937774848474907728L};
    UV n;
    long double fl1 = logl(x);
    long double e = -fl1/2;
    long double s = 0;
    long double t = 1;
    for (n = 1; n <= 100; n++) {
      t *= e/(long double)n;
      s += t*a2[(n-1)>>1];
      /* if (fabsl(t) < LDBL_EPSILON*fabsl(s)) break; */
    }
    return euler_mascheroni + logl(fl1) - 2*sqrtl(x)*s;
  }
#endif

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

int from_digit_string(UV* rn, const char* s, int base)
{
  UV max, n = 0;
  int i, len;

  /* Skip leading -/+ and zeros */
  if (s[0] == '-' || s[0] == '+') s++;
  while (s[0] == '0') s++;

  len = strlen(s);
  max = (UV_MAX-base+1)/base;

  for (i = 0, len = strlen(s); i < len; i++) {
    const char c = s[i];
    int d = !isalnum(c) ? 255 : (c <= '9') ? c-'0' : (c <= 'Z') ? c-'A'+10 : c-'a'+10;
    if (d >= base) croak("Invalid digit for base %d", base);
    if (n > max) return 0;   /* Overflow */
    n = n * base + d;
  }
  *rn = n;
  return 1;
}

int from_digit_to_UV(UV* rn, UV* r, int len, int base)
{
  UV d, n = 0;
  int i;
  if (len < 0 || len > BITS_PER_WORD)
    return 0;
  for (i = 0; i < len; i++) {
    d = r[i];
    if (n > (UV_MAX-d)/base) break;  /* overflow */
    n = n * base + d;
  }
  *rn = n;
  return (i >= len);
}


int from_digit_to_str(char** rstr, UV* r, int len, int base)
{
  char *so, *s;
  int i;

  if (len < 0 || !(base == 2 || base == 10 || base == 16)) return 0;

  if (r[0] >= (UV) base) return 0;  /* TODO: We don't apply extended carry */

  New(0, so, len + 3, char);
  s = so;
  if (base == 2 || base == 16) {
    *s++ = '0';
    *s++ = (base == 2) ? 'b' : 'x';
  }
  for (i = 0; i < len; i++) {
    UV d = r[i];
    s[i] = (d < 10) ? '0'+d : 'a'+d-10;
  }
  s[len] = '\0';
  *rstr = so;
  return 1;
}

int to_digit_array(int* bits, UV n, int base, int length)
{
  int d;

  if (base < 2 || length > 128) return -1;

  if (base == 2) {
    for (d = 0; n; n >>= 1)
      bits[d++] = n & 1;
  } else {
    for (d = 0; n; n /= base)
      bits[d++] = n % base;
  }
  if (length < 0) length = d;
  while (d < length)
    bits[d++] = 0;
  return length;
}

int to_digit_string(char* s, UV n, int base, int length)
{
  int digits[128];
  int i, len = to_digit_array(digits, n, base, length);

  if (len < 0) return -1;
  if (base > 36) croak("invalid base for string: %d", base);

  for (i = 0; i < len; i++) {
    int dig = digits[len-i-1];
    s[i] = (dig < 10) ? '0'+dig : 'a'+dig-10;
  }
  s[len] = '\0';
  return len;
}

/* Oddball primality test.
 * In this file rather than primality.c because it uses factoring (!).
 * Algorithm from Charles R Greathouse IV, 2015 */
static INLINE uint32_t _catalan_v32(uint32_t n, uint32_t p) {
  uint32_t s = 0;
  while (n /= p)  s += n % 2;
  return s;
}
static INLINE uint32_t _catalan_v(UV n, UV p) {
  uint32_t s = 0;
  while (n /= p)  s += n % 2;
  return s;
}
static UV _catalan_mult(UV m, UV p, UV n, UV a) {
  if (p > a) {
    m = mulmod(m, p, n);
  } else {
    UV pow = (n < 4294967296U) ? _catalan_v32(a<<1,p) : _catalan_v(a<<1,p);
    m = (pow == 0) ? m
      : (pow == 1) ? mulmod(m,p,n)
                   : mulmod(m,powmod(p,pow,n),n);
  }
  return m;
}
static int _catalan_vtest(UV n, UV p) {
  while (n /= p)
    if (n % 2)
      return 1;
  return 0;
}
int is_catalan_pseudoprime(UV n) {
  UV m, a;
  int i;

  if (n < 2 || ((n % 2) == 0 && n != 2)) return 0;
  if (is_prob_prime(n)) return 1;

  m = 1;
  a = n >> 1;
  {
    UV factors[MPU_MAX_FACTORS+1];
    int nfactors = factor_exp(n, factors, 0);
    /* Aebi and Cairns 2008, page 9 */
#if BITS_PER_WORD == 32
    if (nfactors == 2)
#else
    if (nfactors == 2 && (n < UVCONST(10000000000)))
#endif
      return 0;
    for (i = 0; i < nfactors; i++) {
      if (_catalan_vtest(a << 1, factors[i]))
        return 0;
    }
  }
  {
    UV seg_base, seg_low, seg_high;
    unsigned char* segment;
    void* ctx;
    m = _catalan_mult(m, 2, n, a);
    m = _catalan_mult(m, 3, n, a);
    m = _catalan_mult(m, 5, n, a);
    ctx = start_segment_primes(7, n, &segment);
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high ) {
        m = _catalan_mult(m, p, n, a);
      } END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
  }
  return (a & 1) ? (m==(n-1)) : (m==1);
}

/* If we have fast CTZ, use this GCD.  See Brent Alg V and FLINT Abhinav Baid */
UV gcdz(UV x, UV y) {
  UV f, x2, y2;

  if (x == 0) return y;
  if (y == 0) return x;

  /* Alternately:  f = ctz(x|y); x >>= ctz(x); y >>= ctz(y); */
  x2 = ctz(x);
  y2 = ctz(y);
  f = (x2 <= y2) ? x2 : y2;
  x >>= x2;
  y >>= y2;

  while (x != y) {
    if (x < y) {
      y -= x;
      y >>= ctz(y);
    } else {
      x -= y;
      x >>= ctz(x);
    }
  }
  return x << f;
}

/* The intermediate values are so large that we can only stay in 64-bit
 * up to 53 or so using the divisor_sum calculations.  So just use a table.
 * Save space by just storing the 32-bit values. */
static const int32_t tau_table[] = {
  0,1,-24,252,-1472,4830,-6048,-16744,84480,-113643,-115920,534612,-370944,-577738,401856,1217160,987136,-6905934,2727432,10661420,-7109760,-4219488,-12830688,18643272,21288960,-25499225,13865712,-73279080,24647168,128406630,-29211840,-52843168,-196706304,134722224,165742416,-80873520,167282496,-182213314,-255874080,-145589976,408038400,308120442,101267712,-17125708,-786948864,-548895690,-447438528
};
#define NTAU (sizeof(tau_table)/sizeof(tau_table[0]))
IV ramanujan_tau(UV n) {
  return (n < NTAU)  ?  tau_table[n]  :  0;
}

static UV _count_class_div(UV s, UV b2) {
  UV h = 0, i, ndivisors, *divs, lim;

  lim = isqrt(b2);
  if (lim*lim == b2) lim--;

  if ((lim-s) < 70) {  /* Iterate looking for divisors */
    for (i = s; i <= lim; i++)
      if (b2 % i == 0)
        h++;
  } else {             /* Walk through all the divisors */
    divs = _divisor_list(b2, &ndivisors);
    for (i = 0; i < ndivisors && divs[i] <= lim; i++)
      if (divs[i] >= s)
        h++;
    Safefree(divs);
  }
  return h;
}

/* Returns 12 * H(n).  See Cohen 5.3.5 or Pari/GP.
 * Pari/GP uses a different method for n > 500000, which is quite a bit
 * faster, but assumes the GRH. */
IV hclassno(UV n) {
  UV nmod4 = n % 4, b2, b, h;
  int square;

  if (n == 0) return -1;
  if (nmod4 == 1 || nmod4 == 2) return 0;
  if (n == 3) return 4;

  b = n & 1;
  b2 = (n+1) >> 2;
  square = is_perfect_square(b2);

  h = divisor_sum(b2,0) >> 1;
  if (b == 1)
    h = 1 + square + ((h - 1) << 1);
  b += 2;

  for (;  b2 = (n + b*b) >> 2, 3*b2 < n;  b += 2) {
    h += (b2 % b == 0)
      +  is_perfect_square(b2)
      +  (_count_class_div(b+1, b2) << 1);
  }
  return 12*h + ((b2*3 == n) ? 4 : square && !(n&1) ? 6 : 0);
}
