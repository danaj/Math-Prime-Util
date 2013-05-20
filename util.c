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

/* These math functions are a clusterfrack.  They're defined by C99, but
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
  /* GLIBC */
  extern long double powl(long double, long double);
  extern long double expl(long double);
  extern long double logl(long double);
  extern long double fabsl(long double);
  extern long double floorl(long double);
#else
  #define powl(x, y)  (long double) pow( (double) (x), (double) (y) )
  #define expl(x)     (long double) exp( (double) (x) )
  #define logl(x)     (long double) log( (double) (x) )
  #define fabsl(x)    (long double) fabs( (double) (x) )
  #define floorl(x)   (long double) floor( (double) (x) )
#endif

#ifndef INFINITY
  #define INFINITY (DBL_MAX + DBL_MAX)
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
#include "util.h"
#include "sieve.h"
#include "factor.h"
#include "cache.h"
#include "lehmer.h"

static int _verbose = 0;
void _XS_set_verbose(int v) { _verbose = v; }
int _XS_get_verbose(void) { return _verbose; }

static int _call_gmp = 0;
void _XS_set_callgmp(int v) { _call_gmp = v; }
int  _XS_get_callgmp(void) { return _call_gmp; }


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


/* Does trial division or prob tests, assuming x not divisible by 2, 3, or 5 */
static int _is_prime7(UV n)
{
  UV limit, i;

  if (n > MPU_PROB_PRIME_BEST)
    return _XS_is_prob_prime(n);  /* We know this works for all 64-bit n */

  limit = isqrt(n);
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
  UV d, m;
  unsigned char mtab;
  const unsigned char* sieve;
  int isprime;

  if (n <= 10) {
    switch (n) {
      case 2: case 3: case 5: case 7:   return 2;  break;
      default:                          break;
    }
    return 0;
  }
  d = n/30;
  m = n - d*30;
  mtab = masktab30[m];  /* Bitmask in mod30 wheel */

  /* Return 0 if a multiple of 2, 3, or 5 */
  if (mtab == 0)
    return 0;

  if (d < NPRIME_SIEVE30)
    return (prime_sieve30[d] & mtab) ? 0 : 2;

  isprime = (n <= get_prime_cache(0, &sieve))
            ?  2*((sieve[d] & mtab) == 0)
            :  -1;
  release_prime_cache(sieve);

  /* return (isprime >= 0)  ?  isprime  :  _is_prime7(n); */
  return (isprime >= 0)  ?  isprime  :  _XS_is_prob_prime(n);
}


UV _XS_next_prime(UV n)
{
  UV d, m;
  const unsigned char* sieve;
  UV sieve_size;

  if (n <= 10) {
    switch (n) {
      case 0: case 1:  return  2; break;
      case 2:          return  3; break;
      case 3: case 4:  return  5; break;
      case 5: case 6:  return  7; break;
      default:         return 11; break;
    }
  }
  if (n < 30*NPRIME_SIEVE30) {
    START_DO_FOR_EACH_SIEVE_PRIME(prime_sieve30, n+1, 30*NPRIME_SIEVE30)
      return p;
    END_DO_FOR_EACH_SIEVE_PRIME;
  }

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

  if (n <= 7)
    return (n <= 2) ? 0 : (n <= 3) ? 2 : (n <= 5) ? 3 : 5;

  d = n/30;
  m = n - d*30;

  if (n < 30*NPRIME_SIEVE30) {
    do {
      m = prevwheel30[m];
      if (m==29) { MPUassert(d>0, "d 0 in prev_prime");  d--; }
    } while (prime_sieve30[d] & masktab30[m]);
    return(d*30+m);
  }

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
   2008,2052,2033,1988,1984,2010,1943,2024,2005,2027,1937,1955,1956,1993,1976};
#define NSTEP_COUNTS_30K  (sizeof(step_counts_30k)/sizeof(step_counts_30k[0]))

/* mpu '$step=300_000; $pc=prime_count(10*$step); print "$pc\n", join(",", map { $spc=$pc; $pc=prime_count($_*$step); $pc-$spc; } 11..100), "\n"' */
static const unsigned short step_counts_300k[] =  /* starts at 3M */
  {20084,19826,19885,19703,19634,19491,19532,19391,19244,19243,19224,19086,
   19124,19036,18942,18893,18870,18853,18837,18775,18688,18674,18594,18525,
   18639,18545,18553,18424,18508,18421,18375,18366,18391,18209,18239,18298,
   18209,18294,18125,18138,18147,18115,18126,18021,18085,18068,18094,17963,
   18041,18003,17900,17881,17917,17888,17880,17852,17892,17779,17823,17764,
   17806,17762,17780,17716,17633,17758,17746,17678,17687,17613,17709,17628,
   17634,17556,17528,17598,17604,17532,17606,17548,17493,17576,17456,17468,
   17555,17452,17407,17472,17415,17500,17508,17418,17463,17240,17345,17351,
   17380,17394,17379,17330,17322,17335,17354,17113,17210,17231,17238,17305,
   17268,17219,17281,17235,17119,17292,17161,17212,17166,17277,17137,17260,
   17228,17197,17154,17097,17195,17136,17067,17058,17041,17045,17187,17034,
   17029,17037,17090,16985,17054,17017,17106,17001,17095,17125,17027,16948,
   16969,17031,16916,17031,16905,16937,16881,16952,16919,16938,17028,16963,
   16902,16922,16944,16901,16847,16969,16900,16876,16841,16874,16894,16861,
   16761,16886,16778,16820,16727,16921,16817,16845,16847,16824,16844,16809,
   16859,16783,16713,16752,16762,16857,16760,16626,16784,16784,16718,16745,
   16871,16635,16714,16630,16779,16709,16660,16730,16715,16724};
#define NSTEP_COUNTS_300K (sizeof(step_counts_300k)/sizeof(step_counts_300k[0]))

static const unsigned int step_counts_30m[] =  /* starts at 30M */
  {1704256,1654839,1624694,1602748,1585989,1571241,1559918,1549840,1540941,
   1533150,1525813,1519922,1513269,1508559,1503386,1497828,1494129,1489905,
   1486417,1482526,1478941,1475577,1472301,1469133,1466295,1464711,1461223,
   1458478,1455327,1454218,1451883,1449393,1447612,1445029,1443285,1442268,
   1438511,1437688,1435603,1433623,1432638,1431158,1429158,1427934,1426191,
   1424449,1423146,1421898,1421628,1419519,1417646,1416274,1414828,1414474,
   1412536,1412147,1410149,1409474,1408847,1406619,1405863,1404699,1403820,
   1402802,1402215,1401459,1399972,1398687,1397968,1397392,1396025,1395311,
   1394081,1393614,1393702,1391745,1390950,1389856,1389245,1388381,1387557,
   1387087,1386285,1386089,1385355,1383659,1383030,1382174,1382128,1380556,
   1379940,1379988,1379181,1378300,1378033,1376974,1376282,1375646,1374445};
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
    UV i = 0;
    if (high < (NSTEP_COUNTS_30K+1) * UVCONST(30000)) {
      while (i < NSTEP_COUNTS_30K && high >= (i+1) * 30000) {
        count += step_counts_30k[i++];
      }
      low = i * 30000;
    } else if (high < (NSTEP_COUNTS_300K+1) * UVCONST(300000)) {
      count = 216816;
      while (i < NSTEP_COUNTS_300K && high >= (i+11) * UVCONST(300000)) {
        count += step_counts_300k[i++];
      }
      low = (i+10) * 300000;
    } else {
      count = 1857859;
      while (i < NSTEP_COUNTS_30M && high >= (i+2) * UVCONST(30000000)) {
        count += step_counts_30m[i++];
      }
      low = (UV)(i+1) * UVCONST(30000000);
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

  /* For relatively small values, generate a sieve and count the results.
   *
   * For larger values, compute a lower bound, use Lehmer's algorithm to get
   * a fast prime count, then start segment sieving from there.
   *
   * For very large values, binary search on Riemann's R function to get a
   * good approximation, use Lehmer's algorithm to get the count, then walk
   * backwards or sieve forwards.
   */
  if (upper_limit <= 32*1024*30) {
    /* Generate a sieve and count. */
    segment_size = get_prime_cache(upper_limit, &cache_sieve) / 30;
    /* Count up everything in the cached sieve. */
    if (segment_size > 0)
      count += count_segment_maxcount(cache_sieve, segment_size, target, &p);
    release_prime_cache(cache_sieve);
  } else {
    double fn = n;
    double flogn = log(fn);
    double flog2n = log(flogn);   /* Dusart 2010, page 2, n >= 3 */
    UV lower_limit = fn * (flogn + flog2n - 1.0 + ((flog2n-2.10)/flogn));
#if BITS_PER_WORD == 32
    if (1) {
#else
    if (n <= UVCONST(20000000000)) {
#endif
      /* Calculate lower limit, get count, sieve to that */
      segment_size = lower_limit / 30;
      lower_limit = 30 * segment_size - 1;
      count = _XS_lehmer_pi(lower_limit) - 3;
      MPUassert(count <= target, "Pi(nth_prime_lower(n))) > n");
    } else {
      /* Compute approximate nth prime via binary search on R(n) */
      UV lo = lower_limit;
      UV hi = upper_limit;
      double lor = _XS_RiemannR(lo);
      double hir = _XS_RiemannR(hi);
      while (lor < hir) {
        UV mid = (UV)  ((lo + hi) / 2);
        double midr = _XS_RiemannR(mid);
        if (midr <= n) { lo = mid+1;  lor = _XS_RiemannR(lo); }
        else           { hi = mid; hir = midr; }
      }
      /* Bias toward lower, because we want to sieve up if possible */
      lower_limit = (UV) (double)(0.9999999*(lo-1));
      segment_size = lower_limit / 30;
      lower_limit = 30 * segment_size - 1;
      count = _XS_lehmer_pi(lower_limit);
      /*
        printf("We've estimated %lu too %s.\n", (count>n)?count-n:n-count, (count>n)?"FAR":"little");
        printf("Our limit %lu %s a prime\n", lower_limit, _XS_is_prime(lower_limit) ? "is" : "is not");
      */

      if (count > n) { /* Too far.  Walk backwards */
        if (_XS_is_prime(lower_limit)) count--;
        for (p = 0; p <= (count-n); p++)
          lower_limit = _XS_prev_prime(lower_limit);
        return lower_limit;
      }
      count -= 3;
    }

    /* Make sure the segment siever won't have to keep resieving. */
    prime_precalc(isqrt(upper_limit));
  }

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

/* Return an IV array with lo-hi+1 elements.  mu[k-lo] = µ(k) for k = lo .. hi.
 * It is the callers responsibility to call Safefree on the result. */
#define PGTLO(p,lo)  ((p) >= lo) ? (p) : ((p)*(lo/(p)) + ((lo%(p))?(p):0))
#define P2GTLO(pinit, p, lo) \
   ((pinit) >= lo) ? (pinit) : ((p)*(lo/(p)) + ((lo%(p))?(p):0))
char* _moebius_range(UV lo, UV hi)
{
  char* mu;
  UV i;
  UV sqrtn = isqrt(hi);
#if 1
  IV* A;

  /* This implementation follows that of Deléglise & Rivat (1996), which is
   * a segmented version of Lioen & van de Lune (1994).  Downside is that it
   * uses an array of IV's, so too much memory.  Pawleicz (2011) shows a small
   * variation but seems to just be a rearrangement -- there is no time or
   * space difference on my machines. (TODO: but maybe for hi > 2^63?)
   */

  New(0, A, hi-lo+1, IV);
  if (A == 0)
    croak("Could not get memory for %"UVuf" moebius results\n", hi-lo+1);
  for (i = lo; i <= hi; i++)
    A[i-lo] = 1;
  START_DO_FOR_EACH_PRIME(2, sqrtn) {
    UV p2 = p*p;
    for (i = PGTLO(p2, lo); i <= hi; i += p2)
      A[i-lo] = 0;
    for (i = PGTLO(p, lo); i <= hi; i += p)
      A[i-lo] *= -(IV)p;
  } END_DO_FOR_EACH_PRIME
  New(0, mu, hi-lo+1, char);
  if (mu == 0)
    croak("Could not get memory for %"UVuf" moebius results\n", hi-lo+1);
  memset(mu, 0, hi-lo+1);
  for (i = lo; i <= hi; i++) {
    IV a = A[i-lo];
    if (a != 0)
      mu[i-lo] = (a != (IV)i && -a != (IV)i)  ?  (a<0) - (a>0)
                                              :  (a>0) - (a<0);
  }
  Safefree(A);
#endif
#if 0
  UV p;
  /* Simple char method, Needs way too many primes */
  New(0, mu, hi-lo+1, char);
  if (mu == 0)
    croak("Could not get memory for %"UVuf" moebius results\n", hi-lo+1);
  memset(mu, 1, hi-lo+1);
  if (lo == 0)  mu[0] = 0;
  prime_precalc( _XS_nth_prime_upper(hi) );
  for (p = 2; p <= hi; p = _XS_next_prime(p)) {
    UV p2 = p*p;
    for (i = PGTLO(p2, lo); i <= hi; i += p2)
      mu[i-lo] = 0;
    for (i = PGTLO(p, lo); i <= hi; i += p)
      mu[i-lo] = -mu[i-lo];
  }
#endif
#if 0
  /* Kuznetsov's transform of Deléglise & Rivat (1996) into logs.
   * (1) I'm using the log function, which should be fixed (easy).
   * (2) it doesn't work.  try 64-101 vs. 64 vs. 100.
   */
  UV p;
  unsigned char* A;
  New(0, A, hi-lo+1, unsigned char);
  if (A == 0)
    croak("Could not get memory for %"UVuf" moebius results\n", hi-lo+1);
  memset(A, 0, hi-lo+1);
  if (sqrtn*sqrtn != hi) sqrtn++;  /* ceil sqrtn */
  prime_precalc(sqrtn);
  for (p = 2; p <= sqrtn; p = _XS_next_prime(p)) {
    UV p2 = p*p;
    unsigned char l = 1 | (unsigned char) ( log(p)/log(2) );
    for (i = PGTLO(p, lo); i <= hi; i += p)
      A[i-lo] += l;
    for (i = PGTLO(p2, lo); i <= hi; i += p2)
      A[i-lo] |= 0x80;
  }

  New(0, mu, hi-lo+1, char);
  if (mu == 0)
    croak("Could not get memory for %"UVuf" moebius results\n", hi-lo+1);
  for (i = lo; i <= hi; i++) {
    unsigned char a = A[i-lo];
    unsigned char log2i = (unsigned char) ( log(i)/log(2) );
    //printf("i = %lu  a = %lu  log2i = %lu\n", i, a, log2i);
    if (a & 0x80) { mu[i-lo] = 0; }
    else if (a >= log2i) { mu[i-lo] =  1 - 2*(a&1); }
    else                 { mu[i-lo] = -1 + 2*(a&1); }
  }
  if (lo == 0)  mu[0] = 0;
  Safefree(A);
#endif
  return mu;
}

UV* _totient_range(UV lo, UV hi) {
  UV* totients;
  UV i, sievehi;
  if (hi < lo) croak("_totient_range error hi %lu < lo %lu\n", hi, lo);
  New(0, totients, hi-lo+1, UV);
  if (totients == 0)
    croak("Could not get memory for %"UVuf" totients\n", hi);
  for (i = lo; i <= hi; i++)
    totients[i-lo] = i;
  sievehi = hi/2;
  for (i=P2GTLO(2*2,2,lo); i <= hi; i += 2) totients[i-lo] -= totients[i-lo]/2;
  START_DO_FOR_EACH_PRIME(3, sievehi) {
    for (i = P2GTLO(2*p,p,lo); i <= hi; i += p)
      totients[i-lo] -= totients[i-lo]/p;
  } END_DO_FOR_EACH_PRIME
  for (i = lo; i <= hi; i++)
    if (totients[i-lo] == i)
      totients[i-lo] = i-1;
  return totients;
}

IV _XS_mertens(UV n) {
#if 0
  /* Benito and Varona 2008, theorem 3.  Segment. */
  IV* mu;
  UV k;
  UV limit = 1000000;
  UV n3 = n/3;
  IV sum = 0;
  UV startk = 1;
  UV endk = limit;
  prime_precalc( (UV) (sqrt(n)+0.5) );
  while (startk <= n3) {
    if (endk > n3) endk = n3;
    mu = _moebius_range(startk, endk);
    for (k = startk; k <= endk; k++)
      if (mu[k-startk] != 0)
        sum += mu[k-startk] * ((n-k)/(2*k));
    Safefree(mu);
    startk = endk+1;
    endk += limit;
  }
  return -sum;
#else
  /* Deléglise and Rivat (1996) using u = n^1/2 and unsegmented. */
  /* Very simple, but they use u = n^1/3 and segment */
  UV u, i, m, nmk;
  char* mu;
  IV* M;
  IV sum;

  if (n <= 1)  return n;
  u = isqrt(n);
  mu = _moebius_range(0, u);
  New(0, M, u+1, IV);
  M[0] = 0;
  for (i = 1; i <= u; i++)
    M[i] = M[i-1] + mu[i];
  sum = M[u];
  for (m = 1; m <= u; m++) {
    if (mu[m] != 0) {
      IV inner_sum = 0;
      UV lower = (u/m) + 1;
      UV last_nmk = n/(m*lower);
      UV this_k = 0;
      UV next_k = n/(m*1);
      for (nmk = 1; nmk <= last_nmk; nmk++) {
        this_k = next_k;
        next_k = n/(m*(nmk+1));
        inner_sum += M[nmk] * (this_k - next_k);
      }
      sum -= mu[m] * inner_sum;
    }
  }
  Safefree(M);
  Safefree(mu);
  return sum;
#endif
}

double _XS_chebyshev_theta(UV n)
{
  KAHAN_INIT(sum);
  START_DO_FOR_EACH_PRIME(2, n) {
    KAHAN_SUM(sum, logl(p));
  } END_DO_FOR_EACH_PRIME
  return (double) sum;
}
double _XS_chebyshev_psi(UV n)
{
  UV mults_are_one = 0;
  long double logn, logp;
  KAHAN_INIT(sum);

  logn = logl(n);
  START_DO_FOR_EACH_PRIME(2, n) {
    logp = logl(p);
    if (!mults_are_one && p > (n/p))   mults_are_one = 1;
    KAHAN_SUM(sum, (mults_are_one) ? logp : logp * floorl(logn/logp + 1e-15));
  } END_DO_FOR_EACH_PRIME
  return (double) sum;
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

/* Riemann Zeta on the real line, with 1 subtracted.
 * Compare to Math::Cephes zetac.  Also zeta with q=1 and subtracting 1.
 *
 * The Cephes zeta function uses a series (2k)!/B_2k which converges rapidly
 * and has a very wide range of values.  We use it here for some values.
 *
 * Note: Calculations here are done on long doubles and we try to generate ~17
 *       digits of accuracy.  When these are returned to Perl they get put in
 *       a standard 64-bit double, so don't expect more than 15 digits.
 *
 * For values 0.5 to 5, this code uses the rational Chebyshev approximation
 * from Cody and Thacher.  This method is extraordinarily fast and very
 * accurate over its range (slightly better than Cephes for most values).  If
 * we had quad floats, we could use the 9-term polynomial.
 */
long double ld_riemann_zeta(long double x) {
  long double const tol = 1e-17;
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

  if (x > 2000.0) {
    /* 1) zeta(2000)-1 is about 8.7E-603, which is far less than a IEEE-754
     *    64-bit double can represent.  A 128-bit quad could go to ~16000.
     * 2) pow / powl start getting obnoxiously slow with values like -7500. */
    return 0.0;
  }

#if 0
  {
    KAHAN_INIT(sum);
    /* Simple defining series, works well. */
    for (i = 5; i <= 1000000; i++) {
      long double term = powl(i, -x);
      KAHAN_SUM(sum, term);
      if (term < tol*sum) break;
    }
    KAHAN_SUM(sum, powl(4, -x) );
    KAHAN_SUM(sum, powl(3, -x) );
    KAHAN_SUM(sum, powl(2, -x) );
    return sum;
  }
#endif

  /* The 2n!/B_2k series used by the Cephes library. */
  {
    /* gp/pari: factorial(2n)/bernfrac(2n) */
    static const long double A[] = {
      12.0L,
     -720.0L,
      30240.0L,
     -1209600.0L,
      47900160.0L,
     -1892437580.3183791606367583212735166426L,
      74724249600.0L,
     -2950130727918.1642244954382084600497650L,
      116467828143500.67248729113000661089202L,
     -4597978722407472.6105457273596737891657L,
      181521054019435467.73425331153534235290L,
     -7166165256175667011.3346447367083352776L,
      282908877253042996618.18640556532523927L,
    };
    long double a, b, s, t;
    const long double w = 10.0;
    s = 0.0;
    b = 0.0;
    for (i = 2; i < 11; i++) {
      b = powl( i, -x );
      s += b;
      if (fabsl(b/s) < tol)
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
      t = fabsl(t/s);
      if (t < tol)
        break;
      a *= x + k + 1.0;
      b /= w;
    }
    return s;
  }
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
