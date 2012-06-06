#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>

#include "util.h"
#include "sieve.h"
#include "factor.h"
#include "ptypes.h"

/*
 * I'm undecided as to whether we want this, or just let the functions alloc
 * and free it per call.
 */
static unsigned char* prime_segment = 0;
unsigned char* get_prime_segment(void) {
  if (prime_segment == 0)
    prime_segment = (unsigned char*) malloc( SEGMENT_CHUNK_SIZE );
  if (prime_segment == 0)
    croak("Could not allocate %"UVuf" bytes for segment sieve", SEGMENT_CHUNK_SIZE);
  return prime_segment;
}
void free_prime_segment(void) {
  if (prime_segment != 0)
    free(prime_segment);
  prime_segment = 0;
}

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

  if ( n < (NPRIME_IS_SMALL*8))
    return ((prime_is_small[n/8] >> (n%8)) & 1) ? 2 : 0;

  d = n/30;
  m = n - d*30;
  mtab = masktab30[m];  /* Bitmask in mod30 wheel */

  /* Return 0 if a multiple of 2, 3, or 5 */
  if (mtab == 0)
    return 0;

  if (n <= get_prime_cache(0, &sieve))
    return ((sieve[d] & mtab) == 0) ? 2 : 0;

  return _is_prime7(n);
}

/* Shortcut, asking for a very quick response of 1 = prime, 0 = dunno.
 * No trial divisions will be done, making this useful for factoring.
 */
int is_definitely_prime(UV n)
{
  UV d, m;
  unsigned char mtab;
  const unsigned char* sieve;

  if ( n < (NPRIME_IS_SMALL*8))
    return ((prime_is_small[n/8] >> (n%8)) & 1);

  d = n/30;
  m = n - d*30;
  mtab = masktab30[m];  /* Bitmask in mod30 wheel */

  /* Return 0 if a multiple of 2, 3, or 5 */
  if (mtab == 0)
    return 0;

  if (n <= get_prime_cache(0, &sieve))
    return ((sieve[d] & mtab) == 0);

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
  m = nextwheel30[m];  if (m == 1) d++;
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

  sieve_size = get_prime_cache(0, &sieve);
  if (n < sieve_size) {
    START_DO_FOR_EACH_SIEVE_PRIME(sieve, n+1, sieve_size)
      return p;
    END_DO_FOR_EACH_SIEVE_PRIME;
    /* Not found, so must be larger than the cache size */
    n = sieve_size;
  }

  d = n/30;
  m = n - d*30;
  m = nextwheel30[m];  if (m == 1) d++;
  while (!_is_prime7(d*30+m)) {
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
      m = prevwheel30[m];  if (m==29) { if (d == 0) return 0;  d--; }
    } while (sieve[d] & masktab30[m]);
  } else {
    do {
      m = prevwheel30[m];  if (m==29) { if (d == 0) return 0;  d--; }
    } while (!_is_prime7(d*30+m));
  }
  return(d*30+m);
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

static const float F1 = 1.0;
UV prime_count_lower(UV x)
{
  float fx, flogx;
  float a = 1.80;

  if (x < NPRIME_COUNT_SMALL)
    return prime_count_small[x];

  fx = (float)x;
  flogx = logf(x);

  if (x < 599)
    return (UV) (fx / (flogx-0.7));

  if      (x <     2700) {  a = 0.30; }
  else if (x <     5500) {  a = 0.90; }
  else if (x <    19400) {  a = 1.30; }
  else if (x <    32299) {  a = 1.60; }
  else if (x <   176000) {  a = 1.80; }
  else if (x <   315000) {  a = 2.10; }
  else if (x <  1100000) {  a = 2.20; }
  else if (x <  4500000) {  a = 2.31; }
  else if (x <233000000) {  a = 2.36; }
  else if (x <240000000) {  a = 2.32; }
  else if (x <UVCONST(0xFFFFFFFF)) {  a = 2.32; }

  return (UV) ( (fx/flogx) * (F1 + F1/flogx + a/(flogx*flogx)) );
}


UV prime_count_upper(UV x)
{
  float fx, flogx;
  float a = 2.51;

  if (x < NPRIME_COUNT_SMALL)
    return prime_count_small[x];

  fx = (float)x;
  flogx = logf(x);

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
  else if (x <UVCONST(0xFFFFFFFF)) {  a = 2.362; }

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
  /* Placeholder for fancy algorithms, like Tomás Oliveira e Silva's:
   *     http://trac.sagemath.org/sage_trac/ticket/8135
   * or an approximation to Li(x) plus a delta.
   */
  UV lower = prime_count_lower(x);
  UV upper = prime_count_upper(x);
  return ((lower + upper) / 2);
}


UV prime_count(UV n)
{
  const unsigned char* sieve;
  static UV last_bytes = 0;
  static UV last_count = 3;
  UV s, bytes;
  UV count = 3;

  if (n < NPRIME_COUNT_SMALL)
    return prime_count_small[n];

  bytes = n/30;
  s = 0;

  if (n <= get_prime_cache(0, &sieve)) {

    /* We have enough primes -- just count them. */

    /* Start from last word position if we can.  This is a big speedup when
     * calling prime_count many times with successively larger numbers. */
    if (bytes >= last_bytes) {
      s = last_bytes;
      count = last_count;
    }

    count += count_zero_bits(sieve+s, bytes-s);

    last_bytes = bytes;
    last_count = count;

    START_DO_FOR_EACH_SIEVE_PRIME(sieve, 30*bytes, n)
      count++;
    END_DO_FOR_EACH_SIEVE_PRIME;

  } else {

    /* We don't have enough primes.  Repeatedly segment sieve */
    UV const segment_size = SEGMENT_CHUNK_SIZE;
    unsigned char* segment;

    /* Get this over with once */
    prime_precalc( sqrt(n) + 2 );

    segment = get_prime_segment();
    if (segment == 0)
      return 0;

    for (s = 0; s <= bytes; s += segment_size) {
      /* We want to sieve one extra byte, to handle the last fragment */
      UV sieve_bytes = ((bytes-s) >= segment_size) ? segment_size : bytes-s+1;
      UV count_bytes = ((bytes-s) >= segment_size) ? segment_size : bytes-s;

      /* printf("sieving from %"UVuf" to %"UVuf"\n", 30*s+1, 30*(s+sieve_bytes-1)+29); */
      if (sieve_segment(segment, s, s + sieve_bytes - 1) == 0) {
        croak("Could not segment sieve from %"UVuf" to %"UVuf, 30*s+1, 30*(s+sieve_bytes)+29);
        break;
      }

      if (count_bytes > 0)
        count += count_zero_bits(segment, count_bytes);

    }
    s -= segment_size;

    /*printf("counting fragment from %"UVuf" to %"UVuf"\n", 30*bytes-30*s, n-30*s); */
    START_DO_FOR_EACH_SIEVE_PRIME(segment, 30*bytes - s*30, n - s*30)
      count++;
    END_DO_FOR_EACH_SIEVE_PRIME;

  }

  return count;
}



static const unsigned short primes_small[] =
  {0,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
   101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,
   193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,
   293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,
   409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499};
#define NPRIMES_SMALL (sizeof(primes_small)/sizeof(primes_small[0]))

/* The nth prime will be more than this number */
UV nth_prime_lower(UV n)
{
  double fn = (double) n;
  double flogn, flog2n, lower;

  if (n < NPRIMES_SMALL)
    return (n==0) ? 0 : primes_small[n]-1;

  flogn  = log(n);
  flog2n = log(flogn);

  /* Dusart 1999, for all n >= 2 */
  lower = fn * (flogn + flog2n - 1.0 + ((flog2n-2.25)/flogn));

  if (lower > (double)UV_MAX)
    return 0;

  return (UV) lower;
}


/* The nth prime will be less than this number */
UV nth_prime_upper(UV n)
{
  double fn = (double) n;
  double flogn, flog2n, upper;

  if (n < NPRIMES_SMALL)
    return primes_small[n]+1;

  flogn  = log(n);
  flog2n = log(flogn);

  if (n >= 39017)
    upper = fn * ( flogn  +  flog2n - 0.9484 ); /* Dusart 1999 */
  else if (n >= 27076)
    upper = fn * (flogn + flog2n - 1.0 + ((flog2n-1.80)/flogn)); /*Dusart 1999*/
  else if (n >= 7022)
    upper = fn * ( flogn  +  0.9385 * flog2n ); /* Robin 1983 */
  else
    upper = fn * ( flogn + flog2n );

  if (upper > (double)UV_MAX)
    return 0;

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
   */
  approx = fn * ( flogn + flog2n - 1
                  + ((flog2n - 2)/flogn)
                  - (((flog2n*flog2n) - 6*flog2n + 11) / (2*flogn*flogn))
                );

  /* Apply a correction to help keep small inputs close. */
  order = flog2n/flogn;
  order = order*order*order * fn;

  if      (n <   259) approx += 10.4 * order;
  else if (n <   775) approx += 7.52 * order;
  else if (n <  1271) approx += 5.6 * order;
  else if (n <  2000) approx += 5.2 * order;
  else if (n <  4000) approx += 4.3 * order;
  else if (n < 12000) approx += 3.0 * order;
  else if (n <150000) approx += 2.1 * order;

  return (UV) rint(approx);
}


/*
 * Given a sieve of size nbytes, walk it counting zeros (primes) until:
 *
 * (1) we counted them all: return the count, which will be less than maxcount.
 *
 * (2) we hit maxcount: set position to the index of the maxcount'th prime
 *     and return count (which will be equal to maxcount).
 */
static UV count_segment(const unsigned char* sieve, UV nbytes, UV maxcount, UV* pos)
{
  UV count = 0;
  UV bytes_left;
  UV byte = 0;

  assert(sieve != 0);
  assert(pos != 0);
  *pos = 0;
  if ( (nbytes == 0) || (maxcount == 0) )
    return 0;

  while ( (byte < nbytes) && (count < maxcount) )
    count += byte_zeros[sieve[byte++]];

  if (count >= maxcount) { /* One too far -- back up */
    count -= byte_zeros[sieve[--byte]];
  }

  assert(count < maxcount);

  if (byte == nbytes)
    return count;

  /* The result is somewhere in the next byte */
  START_DO_FOR_EACH_SIEVE_PRIME(sieve, byte*30+1, nbytes*30-1)
    if (++count == maxcount)  { *pos = p; return count; }
  END_DO_FOR_EACH_SIEVE_PRIME;

  croak("count_segment failed");
  return 0;
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
  if (upper_limit == 0) {
    croak("nth_prime(%"UVuf") would overflow", n);
    return 0;
  }

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
  count += count_segment(cache_sieve, segment_size, target, &p);
  if (count == target)
    return p;

  /* Start segment sieving.  Get memory to sieve into. */
  segbase = segment_size;
  segment_size = SEGMENT_CHUNK_SIZE;
  segment = get_prime_segment();
  if (segment == 0)
    return 0;

  while (count < target) {
    /* Limit the segment size if we know the answer comes earlier */
    if ( (30*(segbase+segment_size)+29) > upper_limit )
      segment_size = (upper_limit - segbase*30 + 30) / 30;

    /* Do the actual sieving in the range */
    if (sieve_segment(segment, segbase, segbase + segment_size-1) == 0) {
      croak("Could not segment sieve from %"UVuf" to %"UVuf, 30*segbase+1, 30*(segbase+segment_size)+29);
      break;
    }

    /* Count up everything in this segment */
    count += count_segment(segment, segment_size, target-count, &p);

    if (count < target)
      segbase += segment_size;
  }
  assert(count == target);
  return ( (segbase*30) + p );
}
