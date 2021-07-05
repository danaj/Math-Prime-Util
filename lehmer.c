#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define USE_PHI_CACHE 1
#define FUNC_isqrt 1
#define FUNC_icbrt 1
#include "lehmer.h"
#include "util.h"
#include "cache.h"
#include "sieve.h"
#include "prime_nth_count.h"
#include "prime_count_cache.h"
#include "legendre_phi.h"

/*****************************************************************************
 *
 * Counting primes with Legendre, Meissel, and Lehmer methods.
 *
 * Since we have a reasonable extended LMO, this is just for demonstration.
 *
 * The first versions of this were in 2012 and 2013, and included a novel
 * phi calculation using iterative list merging, which greatly sped up
 * the calculations compared to recursive phi calculations, even when caching
 * was added.
 *
 * Kim Walisch started his primecount project in mid-2013, which quickly
 * surpassed this in speed.  Currently (2021) his project is substantially
 * faster, as well as having support for the Deleglise-Rivat and
 * Gourdon algorithms, efficient parallelization, and big number support.
 *
 * Reference: Hans Riesel, "Prime Numbers and Computer Methods for
 * Factorization", 2nd edition, 1994.
 */

/* Below this size, just get primecount using standard methods */
#define SIEVE_LIMIT  60000000
/* Bigger prime count cache in Lehmer loop */
#define SIEVE_MULT   1

static int const verbose = 0;
#define STAGE_TIMING 0

#if STAGE_TIMING
 #include <sys/time.h>
 #define DECLARE_TIMING_VARIABLES  struct timeval t0, t1;
 #define TIMING_START   gettimeofday(&t0, 0);
 #define TIMING_END_PRINT(text) \
  { unsigned long long t; \
    gettimeofday(&t1, 0); \
    t = (t1.tv_sec-t0.tv_sec) * 1000000 + (t1.tv_usec - t0.tv_usec); \
    printf("%s: %10.5f\n", text, ((double)t) / 1000000); }
#else
 #define DECLARE_TIMING_VARIABLES
 #define TIMING_START
 #define TIMING_END_PRINT(text)
#endif


static UV P2_with_primes(UV n, UV a, UV b, const uint32_t *primes, uint32_t lastidx)
{
  UV P2, lastw, lastwpc, i;
  UV lastpc = 6 * primes[lastidx];
  void* pcache = prime_count_cache_create(lastpc);

  /* Ensure we have a large enough base sieve */
  prime_precalc(isqrt(n / primes[a+1]));

  P2 = lastw = lastwpc = 0;
  for (i = b; i > a; i--) {
    UV w = n / primes[i];
    lastwpc = (w <= lastpc) ? prime_count_cache_lookup(pcache, w)
                            : lastwpc + segment_prime_count(lastw+1, w);
    lastw = w;
    P2 += lastwpc;
  }
  P2 -= ((b+a-2) * (b-a+1) / 2) - a + 1;

  prime_count_cache_destroy(pcache);
  return P2;
}
/* b = prime_count(isqrt(n)) */
static UV P2(UV n, UV a, UV b)
{
  uint32_t lastidx, *primes;
  UV maxn, P2;

  maxn = nth_prime_upper( b );
  if (maxn > 4294967291U) maxn = 4294967291U;
  lastidx = range_prime_sieve_32(&primes, maxn, 1);
  MPUassert(lastidx >= b, "failed to generate enough primes\n");

  P2 = P2_with_primes(n, a, b, primes, lastidx);

  Safefree(primes);
  return P2;
}


/* Legendre's method.  Interesting and a good test for phi(x,a), but Lehmer's
 * method is much faster (Legendre: a = pi(n^.5), Lehmer: a = pi(n^.25)) */
UV legendre_prime_count(UV n)
{
  UV a, phina;
  if (n < SIEVE_LIMIT)
    return segment_prime_count(2, n);

  a = legendre_prime_count(isqrt(n));
  phina = legendre_phi(n, a);
  return phina + a - 1;
}


/* Meissel's method. */
UV meissel_prime_count(UV n)
{
  UV a, b, sum;
  if (n < SIEVE_LIMIT)
    return segment_prime_count(2, n);

  a = meissel_prime_count(icbrt(n));  /* a = Pi(floor(n^1/3)) [max    192725] */
  b = meissel_prime_count(isqrt(n));  /* b = Pi(floor(n^1/2)) [max 203280221] */

  sum = legendre_phi(n, a) + a - 1 - P2(n, a, b);
  return sum;
}

/* Lehmer's method.  This is basically Riesel's Lehmer function (page 22),
 * with some additional code to help optimize it.  */
UV lehmer_prime_count(UV n)
{
  UV z, a, b, c, sum, i, j, lastidx, lastpc, lastw, lastwpc;
  uint32_t* primes = 0; /* small prime cache, first b=pi(z)=pi(sqrt(n)) */
  void* pcache;  /* Prime count cache */
  DECLARE_TIMING_VARIABLES;

  if (n < SIEVE_LIMIT)
    return segment_prime_count(2, n);

  /* Protect against overflow.  2^32-1 and 2^64-1 are both divisible by 3. */
  if (n == UV_MAX) {
    if ( (n%3) == 0 || (n%5) == 0 || (n%7) == 0 || (n%31) == 0 )
      n--;
    else
      return segment_prime_count(2,n);
  }

  if (verbose > 0) printf("lehmer %lu stage 1: calculate a,b,c \n", n);
  TIMING_START;
  z = isqrt(n);
  a = lehmer_prime_count(isqrt(z));        /* a = Pi(floor(n^1/4)) [max      6542] */
  b = lehmer_prime_count(z);               /* b = Pi(floor(n^1/2)) [max 203280221] */
  c = lehmer_prime_count(icbrt(n));        /* c = Pi(floor(n^1/3)) [max    192725] */
  TIMING_END_PRINT("stage 1")

  if (verbose > 0) printf("lehmer %lu stage 2: phi(x,a) (z=%lu a=%lu b=%lu c=%lu)\n", n, z, a, b, c);
  TIMING_START;
  sum = legendre_phi(n, a) + ((b+a-2) * (b-a+1) / 2);
  TIMING_END_PRINT("phi(x,a)")

  /* The first b primes are used in stage 4.  Hence, primes to isqrt(n). */
  TIMING_START;
  lastidx = range_prime_sieve_32(&primes, isqrt(n), 1);
  MPUassert(lastidx >= b, "failed to generate enough primes\n");
  TIMING_END_PRINT("small primes")
  if (verbose > 0) printf("lehmer %lu stage 3: %lu small primes\n", n, lastidx);

  TIMING_START;
  lastpc = SIEVE_MULT * primes[lastidx];
  if (SIEVE_MULT == 1)
    pcache = prime_count_cache_create_with_primes(primes, lastidx);
  else
    pcache = prime_count_cache_create(lastpc);
  TIMING_END_PRINT("prime count cache")

  TIMING_START;
  /* Speed up all the prime counts by doing a big base sieve */
  prime_precalc( (UV) pow(n, 3.0/5.0) );
  /* Ensure we have the base sieve for big prime_count ( n/primes[i] ). */
  /* This is about 75k for n=10^13, 421k for n=10^15, 2.4M for n=10^17 */
  prime_precalc(isqrt(n / primes[a+1]));
  TIMING_END_PRINT("sieve precalc")

  if (verbose > 0) printf("lehmer %lu stage 4: loop %lu to %lu, pc to %lu\n", n, a+1, b, n/primes[a+1]);
  TIMING_START;
  /* Reverse the i loop so w increases.  Count w in segments. */
  lastw = 0;
  lastwpc = 0;
  for (i = b; i >= a+1; i--) {
    UV w = n / primes[i];
    lastwpc = (w <= lastpc) ? prime_count_cache_lookup(pcache, w)
                            : lastwpc + segment_prime_count(lastw+1, w);
    lastw = w;
    sum = sum - lastwpc;
    if (i <= c) {
      UV bi = prime_count_cache_lookup(pcache, isqrt(w));
      for (j = i; j <= bi; j++) {
        sum = sum - prime_count_cache_lookup(pcache, w / primes[j]) + j - 1;
      }
      /* We could wrap the +j-1 in:  sum += ((bi+1-i)*(bi+i))/2 - (bi-i+1); */
    }
  }
  TIMING_END_PRINT("stage 4")
  prime_count_cache_destroy(pcache);
  Safefree(primes);
  return sum;
}


/* The Lagarias-Miller-Odlyzko method.
 * Naive implementation without optimizations.
 * About the same speed as Lehmer, a bit less memory.
 * A better implementation can be 10-50x faster and much less memory.
 */
UV LMOS_prime_count(UV n)
{
  UV n13, a, b, sum, i, j, k, c, lastidx, P2, S1, S2;
  uint32_t primec;
  uint32_t* primes = 0;  /* small prime cache */
  signed char* mu = 0;   /* moebius to n^1/3 */
  uint32_t*   lpf = 0;   /* least prime factor to n^1/3 */
  void *pcache; /* Cache for recursive phi */
  DECLARE_TIMING_VARIABLES;

  if (n < SIEVE_LIMIT)
    return segment_prime_count(2, n);

  n13 = icbrt(n);                    /* n13 =  floor(n^1/3)  [max    2642245] */
  a = lehmer_prime_count(n13);       /* a = Pi(floor(n^1/3)) [max     192725] */
  b = lehmer_prime_count(isqrt(n));  /* b = Pi(floor(n^1/2)) [max  203280221] */

  TIMING_START;
  lastidx = range_prime_sieve_32(&primes, isqrt(n), 1);
  MPUassert(lastidx >= b, "failed to generate enough primes\n");
  TIMING_END_PRINT("small primes")

  TIMING_START;
  New(0, mu, n13+1, signed char);
  memset(mu, 1, sizeof(signed char) * (n13+1));
  Newz(0, lpf, n13+1, uint32_t);
  mu[0] = 0;
  for (i = 1; i <= n13; i++) {
    UV primei = primes[i];
    for (j = primei; j <= n13; j += primei) {
      mu[j] = -mu[j];
      if (lpf[j] == 0) lpf[j] = primei;
    }
    k = primei * primei;
    for (j = k; j <= n13; j += k)
      mu[j] = 0;
  }
  lpf[1] = UVCONST(4294967295);  /* Set lpf[1] to max */

  /* Remove mu[i] == 0 using lpf */
  for (i = 1; i <= n13; i++)
    if (mu[i] == 0)
      lpf[i] = 0;
  TIMING_END_PRINT("mu and lpf")

  /* Thanks to Kim Walisch for help with the S1+S2 calculations. */
  c = (a < tiny_phi_max_a()) ? a : tiny_phi_max_a();
  primec = primes[c];
  S1 = 0;
  S2 = 0;
  pcache = prepare_cached_legendre_phi(n, a);

  TIMING_START;
  for (i = 1; i <= n13; i++)
    if (lpf[i] > primec)
      S1 += mu[i] * tiny_phi(n/i, c);
  TIMING_END_PRINT("S1")

  /* TODO: This is about 1.5x slower than the old way.  Fix. */
  TIMING_START;
  for (i = c; i+1 < a; i++) {
    uint32_t p = primes[i+1];
    for (j = (n13/p)+1; j <= n13; j++)
      if (lpf[j] > p)
        S2 += -mu[j] * cached_legendre_phi(pcache, n / (j*p), i);
  }
  TIMING_END_PRINT("S2")
  destroy_cached_legendre_phi(pcache);
  Safefree(lpf);
  Safefree(mu);

  TIMING_START;
  prime_precalc( (UV) pow(n, 2.9/5.0) );
  P2 = P2_with_primes(n, a, b, primes, lastidx);
  TIMING_END_PRINT("P2")
  Safefree(primes);

  /* printf("S1 = %lu\nS2 = %lu\na  = %lu\nP2 = %lu\n", S1, S2, a, P2); */
  sum = (S1 + S2) + a - 1 - P2;
  return sum;
}
