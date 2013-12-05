#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Below this size, just sieve (with table speedup). */
#define SIEVE_LIMIT  60000000
#define MAX_PHI_MEM  (256*1024*1024)

/*****************************************************************************
 *
 * Lehmer prime counting utility.  Calculates pi(x), count of primes <= x.
 *
 * Copyright (c) 2012-2013 Dana Jacobsen (dana@acm.org).
 * This is free software; you can redistribute it and/or modify it under
 * the same terms as the Perl 5 programming language system itself.
 *
 * This file is part of the Math::Prime::Util Perl module, but also can be
 * compiled as a standalone UNIX program using the primesieve package.
 *
 *    g++ -O3 -DPRIMESIEVE_STANDALONE lehmer.c -o prime_count -lprimesieve
 *
 * For faster prime counting in stage 4 with multiprocessor machines:
 *
 *    g++ -O3 -DPRIMESIEVE_STANDALONE -DPRIMESIEVE_PARALLEL lehmer.c -o prime_count -lprimesieve -lgomp
 *
 * The phi(x,a) calculation is unique, to the best of my knowledge.  It uses
 * two lists of all x values + signed counts for the given 'a' value, and walks
 * 'a' down until it is small enough to calculate directly using Mapes' method.
 * This is relatively fast and low memory compared to many other solutions.
 * As with all Lehmer-Meissel-Legendre algorithms, memory use will be a
 * constraint with large values of x (see the table below).
 *
 * If you want something better, I highly recommend the paper "Computing
 * Pi(x): the combinatorial method" (2006) by Tomás Oliveira e Silva.  His
 * implementation is certainly much faster and lower memory than this.  I have
 * briefly run Christian Bau's LMO implementation which has the big advantage
 * of using almost no memory.  On the same machine that ran the times below,
 * I got 10^16 in 36s, 10^17 in 165s, 10^18 in 769s, 10^19 in 4162s.  All
 * using a single core.  That's 10-50x faster than this code.
 *
 * Using my sieve code with everything running in serial, calculating pi(10^12)
 * is done under 1 second on my computer.  pi(10^14) takes under 30 seconds,
 * pi(10^16) in under 20 minutes.  Compared with Thomas R. Nicely's pix4
 * program, this one is 5x faster and uses 10x less memory.  When compiled
 * with parallel primesieve it is over 10x faster.
 *   pix4(10^16) takes 124 minutes, this code + primesieve takes < 4 minutes.
 *
 * Timings with Perl + MPU with all-serial computation, no memory limit.
 * The last column is the standalone time with 12-core parallel primesieve.
 *
 *    n     phi(x,a) mem/time  |  stage 4 mem/time  | total time | pps time
 *   10^19  17884MB   1979.41  | 9144MB             |            | 7h 26m
 *   10^18   5515MB    483.46  | 2394MB             |            | 87m  0s
 *   10^17   1698MB    109.56  |  634MB   9684.1    | 163m 36  s | 17m 45s
 *   10^16    522MB     24.14  |  171MB   1066.3    |  18m 12  s |  3m 44s
 *   10^15    159MB      5.58  |   47MB    142.5    |   2m 28  s |   47.66 s
 *   10^14     48MB      1.28  |   13MB     22.26   |    23.61 s |   10.25 s
 *   10^13     14MB      0.294 |    4MB      3.82   |     4.12 s |    2.21 s
 *   10^12      4MB      0.070 |    1MB      0.707  |     0.78 s |    0.459s
 *   10^11      1MB      0.015 |             0.135  |     0.158s |    0.097s
 *   10^10               0.003 |             0.029  |     0.028s |    0.025s
 *
 * Using the standalone program with parallel primesieve speeds up stage 4
 * a lot for large values, as can be seen by the last column.  It does use
 * quite a bit more memory in stage 4 however (lowering SIEVE_MULT can reduce
 * it by as much as 10x, at the expense of some performance).
 *
 * Reference: Hans Riesel, "Prime Numbers and Computer Methods for
 * Factorization", 2nd edition, 1994.
 */

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


#ifdef PRIMESIEVE_STANDALONE

/* countPrimes can be pretty slow for small ranges, so sieve more small primes
 * and count using binary search.  Uses a lot of memory though.  For big
 * ranges, countPrimes is really fast.  If you use primesieve 4.2, the
 * crossover point is lower (better). */
#define SIEVE_MULT   10

#include <limits.h>
#include <sys/time.h>
#ifdef PRIMESIEVE_PARALLEL
 #include <primesieve/soe/ParallelPrimeSieve.h>
 ParallelPrimeSieve ps;
 #define SET_PPS_PARALLEL ps.setNumThreads(ParallelPrimeSieve::getMaxThreads())
 #define SET_PPS_SERIAL   ps.setNumThreads(1)
#else
 #include <primesieve/soe/PrimeSieve.h>
 PrimeSieve ps;
 #define SET_PPS_PARALLEL /* */
 #define SET_PPS_SERIAL   /* */
#endif

#ifdef _OPENMP
  #include <omp.h>
  int omp_threads = 8;  /* will get replaced */
#endif

/* Translations from Perl + Math::Prime::Util  to  C/C++ + primesieve */
typedef unsigned long UV;
typedef   signed long IV;
#define UV_MAX ULONG_MAX
#define UVCONST(x) ((unsigned long)x##UL)
#define New(id, mem, size, type)  mem = (type*) malloc((size)*sizeof(type))
#define Newz(id, mem, size, type) mem = (type*) calloc(size, sizeof(type))
#define Renew(mem, size, type)    mem = (type*) realloc(mem,(size)*sizeof(type))
#define Safefree(mem)             free((void*)mem)
#define _XS_prime_count(a, b)     ps.countPrimes(a, b)
#define croak(fmt,...)            { printf(fmt,##__VA_ARGS__); exit(1); }
#define prime_precalc(n)          /* */
#define BITS_PER_WORD             ((ULONG_MAX <= 4294967295UL) ? 32 : 64)

static UV isqrt(UV n) {
  UV root;
  if (sizeof(UV) == 8 && n >= 18446744065119617025UL)  return 4294967295UL;
  if (sizeof(UV) == 4 && n >= 4294836225UL)            return 65535UL;
  root = (UV) sqrt((double)n);
  while (root*root > n)  root--;
  while ((root+1)*(root+1) <= n)  root++;
  return root;
}
static UV icbrt(UV n) {
  UV b, root = 0;
  int s;
  if (sizeof(UV) == 8) {
    s = 63;  if (n >= 18446724184312856125UL)  return 2642245UL;
  } else {
    s = 30;  if (n >= 4291015625UL)            return 1625UL;
  }
  for ( ; s >= 0; s -= 3) {
    root += root;
    b = 3*root*(root+1)+1;
    if ((n >> s) >= b) {
      n -= b << s;
      root++;
    }
  }
  return root;
}

/* Callback used for creating an array of primes. */
static uint32_t* sieve_array = 0;
static UV sieve_k;
static UV sieve_n;
class stop_primesieve : public std::exception { };
void primesieve_callback(uint64_t pk) {
  if (sieve_k > sieve_n) throw stop_primesieve();
  /*
  if (pk < sieve_array[sieve_k-1])
    croak("Primes generated out of order.  Switch to serial mode.\n");
  */
  sieve_array[sieve_k++] = pk;
}

/* Generate an array of n small primes, where the kth prime is element p[k].
 * Remember to free when done. */
#define TINY_PRIME_SIZE 20000
static uint32_t* tiny_primes = 0;
static uint32_t* generate_small_primes(UV n)
{
  uint32_t* primes;
  double fn = (double)n;
  double flogn  = log(fn);
  double flog2n  = log(flogn);
  UV nth_prime =  /* Dusart 2010 for > 179k, custom for 18-179k */
     (n >= 688383) ? (UV) ceil(fn*(flogn+flog2n-1.0+((flog2n-2.00)/flogn))) :
     (n >= 178974) ? (UV) ceil(fn*(flogn+flog2n-1.0+((flog2n-1.95)/flogn))) :
     (n >= 18)     ? (UV) ceil(fn*(flogn+flog2n-1.0+((flog2n+0.30)/flogn)))
                   : 59;
  New(0, primes, n+1, uint32_t);
  if (primes == 0)
    croak("Can not allocate small primes\n");
  if (n < TINY_PRIME_SIZE) {
    if (tiny_primes == 0)
      tiny_primes = generate_small_primes(TINY_PRIME_SIZE+1);
    memcpy(primes, tiny_primes, (n+1) * sizeof(uint32_t));
    return primes;
  }
  primes[0] = 0;
  sieve_array = primes;
  sieve_n = n;
  sieve_k = 1;
  {
    PrimeSieve sps;
    try {
      sps.generatePrimes(2, nth_prime, primesieve_callback);
    } catch (stop_primesieve&) { }
  }
  sieve_array = 0;
  return primes;
}

#else

/* We will use pre-sieving to speed up counting for small ranges */
#define SIEVE_MULT   1

#define FUNC_icbrt 1
#include "lehmer.h"
#include "util.h"
#include "cache.h"
#include "sieve.h"

/* Generate an array of n small primes, where the kth prime is element p[k].
 * Remember to free when done. */
static uint32_t* generate_small_primes(UV n)
{
  uint32_t* primes;
  UV  i = 0;
  double fn = (double)n;
  double flogn  = log(fn);
  double flog2n  = log(flogn);
  UV nth_prime =  /* Dusart 2010 for > 179k, custom for 18-179k */
     (n >= 688383) ? (UV) ceil(fn*(flogn+flog2n-1.0+((flog2n-2.00)/flogn))) :
     (n >= 178974) ? (UV) ceil(fn*(flogn+flog2n-1.0+((flog2n-1.95)/flogn))) :
     (n >= 18)     ? (UV) ceil(fn*(flogn+flog2n-1.0+((flog2n+0.30)/flogn)))
                   : 59;

  if (n > 203280221)
    croak("generate small primes with argument too large: %lu\n", (unsigned long)n);
  New(0, primes, n+1, uint32_t);
  if (primes == 0)
    croak("Can not allocate small primes\n");
  primes[0] = 0;
  START_DO_FOR_EACH_PRIME(2, nth_prime) {
    if (i >= n) break;
    primes[++i] = p;
  } END_DO_FOR_EACH_PRIME
  if (i < n)
    croak("Did not generate enough small primes.\n");
  if (verbose > 1) printf("generated %lu small primes, from 2 to %lu\n", i, (unsigned long)primes[i]);
  return primes;
}
#endif


/* Given an array of primes[1..lastprime], return Pi(n) where n <= lastprime.
 * This is actually quite fast, and definitely faster than sieving.  By using
 * this we can avoid caching prime counts and also skip most calls to the
 * segment siever.
 */
static UV bs_prime_count(uint32_t n, uint32_t const* const primes, uint32_t lastidx)
{
  UV i, j;
  if (n <= 2)  return (n == 2);
  /* if (n > primes[lastidx])  return _XS_prime_count(2, n); */
  if (n >= primes[lastidx]) {
    if (n == primes[lastidx]) return lastidx;
    croak("called bspc(%u) with counts up to %u\n", n, primes[lastidx]);
  }
  j = lastidx;
  if (n < 8480) {
    i = 1 + (n>>4);
    if (j > 1060) j = 1060;
  } else if (n < 25875000) {
    i = 793 + (n>>5);
    if (j > (n>>3)) j = n>>3;
  } else {
    i = 1617183;
    if (j > (n>>4)) j = n>>4;
  }
  while (i < j) {
    UV mid = (i+j)/2;
    if (primes[mid] <= n)  i = mid+1;
    else                   j = mid;
  }
  return i-1;
}

#define FAST_DIV(x,y) \
  ( ((x) <= 4294967295U) ? (uint32_t)(x)/(uint32_t)(y) : (x)/(y) )

/* Use Mapes' method to calculate phi(x,a) for small a.  This is really
 * convenient and a little Perl script will spit this code out for whatever
 * limit we select.  It gets unwieldy with large a values.
 */
static UV mapes(UV x, UV a)
{
  IV val;
  if (a == 0)  return x;
  if (a == 1)  return x-x/2;
  val = x-x/2-x/3+x/6;
  if (a >= 3) val += 0-x/5+x/10+x/15-x/30;
  if (a >= 4) val += 0-x/7+x/14+x/21-x/42+x/35-x/70-x/105+x/210;
  if (a >= 5) val += 0-x/11+x/22+x/33-x/66+x/55-x/110-x/165+x/330+x/77-x/154-x/231+x/462-x/385+x/770+x/1155-x/2310;
  if (a >= 6) val += 0-x/13+x/26+x/39-x/78+x/65-x/130-x/195+x/390+x/91-x/182-x/273+x/546-x/455+x/910+x/1365-x/2730+x/143-x/286-x/429+x/858-x/715+x/1430+x/2145-x/4290-x/1001+x/2002+x/3003-x/6006+x/5005-x/10010-x/15015+x/30030;
  if (a >= 7) val += 0-x/17+x/34+x/51-x/102+x/85-x/170-x/255+x/510+x/119-x/238-x/357+x/714-x/595+x/1190+x/1785-x/3570+x/187-x/374-x/561+x/1122-x/935+x/1870+x/2805-x/5610-x/1309+x/2618+x/3927-x/7854+x/6545-x/13090-x/19635+x/39270+x/221-x/442-x/663+x/1326-x/1105+x/2210+x/3315-x/6630-x/1547+x/3094+x/4641-x/9282+x/7735-x/15470-x/23205+x/46410-x/2431+x/4862+x/7293-x/14586+x/12155-x/24310-x/36465+x/72930+x/17017-x/34034-x/51051+x/102102-x/85085+x/170170+x/255255-x/510510;
  return (UV) val;
}

#define mapes7(x) (((x) <= 4294967295U) ? mapes7_32(x) : mapes7_64(x))
static UV mapes7_64(UV x) {    /* A tiny bit faster setup for a=7 */
  UV val = x-x/2-x/3-x/5+x/6-x/7+x/10-x/11-x/13+x/14+x/15-x/17+x/21+x/22+x/26
          -x/30+x/33+x/34+x/35+x/39-x/42+x/51+x/55+x/65-x/66-x/70+x/77-x/78
          +x/85+x/91-x/102-x/105-x/110+x/119-x/130+x/143-x/154-x/165-x/170
          -x/182+x/187-x/195+x/210+x/221-x/231-x/238-x/255-x/273-x/286+x/330
          -x/357-x/374-x/385+x/390-x/429-x/442-x/455+x/462+x/510+x/546-x/561
          -x/595-x/663+x/714;
  if (x >= 715) {
    val += 0-x/715+x/770+x/858+x/910-x/935-x/1001-x/1105+x/1122+x/1155+x/1190
            -x/1309+x/1326+x/1365+x/1430-x/1547+x/1785+x/1870+x/2002+x/2145
            +x/2210-x/2310-x/2431+x/2618-x/2730+x/2805+x/3003+x/3094+x/3315
            -x/3570+x/3927-x/4290+x/4641+x/4862+x/5005-x/5610-x/6006+x/6545
            -x/6630+x/7293+x/7735-x/7854;
    if (x >= 9282)
      val += 0-x/9282-x/10010+x/12155-x/13090-x/14586-x/15015-x/15470+x/17017
              -x/19635-x/23205-x/24310+x/30030-x/34034-x/36465+x/39270+x/46410
              -x/51051+x/72930-x/85085+x/102102+x/170170+x/255255-x/510510;
  }
  return val;
}

static uint32_t mapes7_32(uint32_t x) {
  uint32_t val = x-x/2-x/3-x/5+x/6-x/7+x/10-x/11-x/13+x/14+x/15-x/17+x/21
          +x/22+x/26-x/30+x/33+x/34+x/35+x/39-x/42+x/51+x/55+x/65-x/66-x/70
          +x/77-x/78+x/85+x/91-x/102-x/105-x/110+x/119-x/130+x/143-x/154
          -x/165-x/170-x/182+x/187-x/195+x/210+x/221-x/231-x/238-x/255-x/273
          -x/286+x/330-x/357-x/374-x/385+x/390-x/429-x/442-x/455+x/462+x/510
          +x/546-x/561-x/595-x/663+x/714;
  if (x >= 715) {
    val += 0-x/715+x/770+x/858+x/910-x/935-x/1001-x/1105+x/1122+x/1155+x/1190
            -x/1309+x/1326+x/1365+x/1430-x/1547+x/1785+x/1870+x/2002+x/2145
            +x/2210-x/2310-x/2431+x/2618-x/2730+x/2805+x/3003+x/3094+x/3315
            -x/3570+x/3927-x/4290+x/4641+x/4862+x/5005-x/5610-x/6006+x/6545
            -x/6630+x/7293+x/7735-x/7854;
    if (x >= 9282)
      val += 0-x/9282-x/10010+x/12155-x/13090-x/14586-x/15015-x/15470+x/17017
              -x/19635-x/23205-x/24310+x/30030-x/34034-x/36465+x/39270+x/46410
              -x/51051+x/72930-x/85085+x/102102+x/170170+x/255255-x/510510;
  }
  return val;
}

/* Max memory = 2*A*X bytes, e.g. 2*400*24000 = 18.3 MB */
#define PHICACHEA 257
#define PHICACHEX 32769
typedef struct
{
  uint32_t max[PHICACHEA];
  int16_t* val[PHICACHEA];
} cache_t;
static void phicache_init(cache_t* cache) {
  int a;
  for (a = 0; a < PHICACHEA; a++) {
    cache->val[a] = 0;
    cache->max[a] = 0;
  }
}
static void phicache_free(cache_t* cache) {
  int a;
  for (a = 0; a < PHICACHEA; a++) {
    if (cache->val[a] != 0)
      Safefree(cache->val[a]);
    cache->val[a] = 0;
    cache->max[a] = 0;
  }
}

#define PHI_CACHE_POPULATED(x, a) \
  ((a) < PHICACHEA && (UV) cache->max[a] > (x) && cache->val[a][x] != 0)

static void phi_cache_insert(uint32_t x, uint32_t a, IV sum, cache_t* cache) {
  uint32_t cap = ( (x+32) >> 5) << 5;
  if (cache->val[a] == 0) {
    Newz(0, cache->val[a], cap, int16_t);
    cache->max[a] = cap;
  } else if (cache->max[a] < cap) {
    uint32_t i;
    Renew(cache->val[a], cap, int16_t);
    for (i = cache->max[a]; i < cap; i++)
      cache->val[a][i] = 0;
    cache->max[a] = cap;
  }
  if (sum < SHRT_MIN || sum > SHRT_MAX)
    croak("phi(%u,%u) 16-bit overflow: sum = %ld\n", x, a, sum);
  if (cache->val[a] == 0)
    croak("phi cache allocation failure");
  cache->val[a][x] = sum;
}

static IV _phi3(UV x, UV a, int sign, const uint32_t* const primes, const uint32_t lastidx, cache_t* cache)
{
  IV sum;

  if (a < 3)
    return sign * ((a==0) ? x : (a==1) ? x-x/2 : x-x/2-x/3+x/6);
  else if (PHI_CACHE_POPULATED(x, a))
    return sign * cache->val[a][x];
  else if (x < primes[a+1])
    sum = sign;
  else if (1 && x <= primes[lastidx] && x < primes[a]*primes[a])
    sum = sign * (bs_prime_count(x, primes, lastidx) - a + 1);
  else {
    UV a2;
    if (a*a > x) {
      UV ix = isqrt(x);
      UV iters = bs_prime_count(ix, primes, a+1);
      sum = sign * (x - a + iters);
      for (a2 = 1; a2 <= iters; a2++)
        sum += _phi3( FAST_DIV(x, primes[a2]), a2-1, -sign, primes, lastidx, cache);
    } else if (a >= 7) {
      if (PHI_CACHE_POPULATED(x, 7))
        sum = sign * cache->val[7][x];
      else
        sum = sign * mapes7(x);
      for (a2 = 8; a2 <= a; a2++)
        sum += _phi3( FAST_DIV(x,primes[a2]), a2-1, -sign, primes, lastidx, cache);
    } else {
      if (PHI_CACHE_POPULATED(x, a))
        sum = sign * cache->val[a][x];
      else
        sum = sign * mapes(x, a);
    }
  }
  if (a < PHICACHEA && x < PHICACHEX)
    phi_cache_insert(x, a, sign * sum, cache);
  return sum;
}
#define phi_small(x, a, primes, lastidx, cache)  _phi3(x, a, 1, primes, lastidx, cache)

/******************************************************************************/
/*   In-order lists for manipulating our UV value / IV count pairs            */
/******************************************************************************/

typedef struct {
  UV v;
  IV c;
} vc_t;

typedef struct {
  vc_t* a;
  UV size;
  UV n;
} vcarray_t;

static vcarray_t vcarray_create(void)
{
  vcarray_t l;
  l.a = 0;
  l.size = 0;
  l.n = 0;
  return l;
}
static void vcarray_destroy(vcarray_t* l)
{
  if (l->a != 0) {
    if (verbose > 2) printf("FREE list %p\n", l->a);
    Safefree(l->a);
  }
  l->size = 0;
  l->n = 0;
}
/* Insert a value/count pair.  We do this indirection because about 80% of
 * the calls result in a merge with the previous entry. */
static void vcarray_insert(vcarray_t* l, UV val, IV count)
{
  UV n = l->n;
  if (n > 0 && l->a[n-1].v < val)
    croak("Previous value was %lu, inserting %lu out of order\n", l->a[n-1].v, val);
  if (n >= l->size) {
    UV new_size;
    if (l->size == 0) {
      new_size = 20000;
      if (verbose>2) printf("ALLOCing list, size %lu (%luk)\n", new_size, new_size*sizeof(vc_t)/1024);
      New(0, l->a, new_size, vc_t);
    } else {
      new_size = (UV) (1.5 * l->size);
      if (verbose>2) printf("REALLOCing list %p, new size %lu (%luk)\n",l->a,new_size, new_size*sizeof(vc_t)/1024);
      Renew( l->a, new_size, vc_t );
    }
    if (l->a == 0) croak("could not allocate list\n");
    l->size = new_size;
  }
  /* printf(" inserting %lu %ld\n", val, count); */
  l->a[n].v = val;
  l->a[n].c = count;
  l->n++;
}

/* Merge the two sorted lists A and B into A.  Each list has no duplicates,
 * but they may have duplications between the two.  We're quite interested
 * in saving memory, so first remove all the duplicates, then do an in-place
 * merge. */
static void vcarray_merge(vcarray_t* a, vcarray_t* b)
{
  long ai, bi, bj, k, kn;
  long an = a->n;
  long bn = b->n;
  vc_t* aa = a->a;
  vc_t* ba = b->a;

  /* Merge anything in B that appears in A. */
  for (ai = 0, bi = 0, bj = 0; bi < bn; bi++) {
    /* Skip forward in A until empty or aa[ai].v <= ba[bi].v */
    UV bval = ba[bi].v;
    while (ai < an && aa[ai].v > bval)
      ai++;
    /* if A empty then copy the remaining elements */
    if (ai >= an) {
      if (bi == bj)
        bj = bn;
      else
        while (bi < bn)
          ba[bj++] = ba[bi++];
      break;
    }
    if (aa[ai].v == bval)
      aa[ai].c += ba[bi].c;
    else
      ba[bj++] = ba[bi];
  }
  if (verbose>3) printf("  removed %lu duplicates from b\n", bn - bj);
  bn = bj;

  if (bn == 0) {  /* In case they were all duplicates */
    b->n = 0;
    return;
  }

  /* kn = the final merged size.  All duplicates are gone, so this is exact. */
  kn = an+bn;
  if ((long)a->size < kn) {  /* Make A big enough to hold kn elements */
    UV new_size = (UV) (1.2 * kn);
    if (verbose>2) printf("REALLOCing list %p, new size %lu (%luk)\n", a->a, new_size, new_size*sizeof(vc_t)/1024);
    Renew( a->a, new_size, vc_t );
    aa = a->a;  /* this could have been changed by the realloc */
    a->size = new_size;
  }

  /* merge A and B.  Very simple using reverse merge. */
  ai = an-1;
  bi = bn-1;
  for (k = kn-1; k >= 0; k--) {
    if (ai < 0) { /* A is exhausted, just filling in B */
      if (bi < 0) croak("ran out of data during merge");
      aa[k] = ba[bi--];
    } else if (bi < 0) { /* We've caught up with A */
      break;
    } else if (aa[ai].v < ba[bi].v) {
      aa[k] = aa[ai--];
    } else {
      if (aa[ai].v == ba[bi].v) croak("deduplication error");
      aa[k] = ba[bi--];
    }
  }
  a->n = kn;    /* A now has this many items */
  b->n = 0;     /* B is marked empty */
}

static void vcarray_remove_zeros(vcarray_t* a)
{
  long ai = 0;
  long aj = 0;
  long an = a->n;
  vc_t* aa = a->a;

  while (aj < an) {
    if (aa[aj].c != 0) {
      if (ai != aj)
        aa[ai] = aa[aj];
      ai++;
    }
    aj++;
  }
  a->n = ai;
}


/*
 * The main phi(x,a) algorithm.  In this implementation, it takes under 10%
 * of the total time for the Lehmer algorithm, but is a big memory consumer.
 */
#define NTHRESH (MAX_PHI_MEM/16)

static UV phi(UV x, UV a)
{
  UV i, val, sval, lastidx, lastprime;
  UV sum = 0;
  IV count;
  const uint32_t* primes;
  vcarray_t a1, a2;
  vc_t* arr;
  cache_t pcache; /* Cache for recursive phi */

  if (a == 1)  return ((x+1)/2);
  if (a <= 7)  return mapes(x, a);

  lastidx = a+1;
  primes = generate_small_primes(lastidx);
  if (primes == 0)
    croak("Could not generate primes for phi(%lu,%lu)\n", x, a);
  lastprime = primes[lastidx];
  if (x < lastprime)  { Safefree(primes); return (x > 0) ? 1 : 0; }
  phicache_init(&pcache);

  a1 = vcarray_create();
  a2 = vcarray_create();
  vcarray_insert(&a1, x, 1);

  while (a > 7) {
    UV primea = primes[a];
    UV sval_last = 0;
    IV sval_count = 0;
    arr = a1.a;
    for (i = 0; i < a1.n; i++) {
      count = arr[i].c;
      val  = arr[i].v;
      sval = FAST_DIV(val, primea);
      if (sval < primea) break;      /* stop inserting into a2 if small */
      if (sval != sval_last) {       /* non-merged value.  Insert into a2 */
        if (sval_last != 0) {
          if (sval_last <= lastprime && sval_last < primes[a-1]*primes[a-1])
            sum += sval_count*(bs_prime_count(sval_last,primes,lastidx)-a+2);
          else
            vcarray_insert(&a2, sval_last, sval_count);
        }
        sval_last = sval;
        sval_count = 0;
      }
      sval_count -= count;           /* Accumulate count for this sval */
    }
    if (sval_last != 0) {            /* Insert the last sval */
      if (sval_last <= lastprime && sval_last < primes[a-1]*primes[a-1])
        sum += sval_count*(bs_prime_count(sval_last,primes,lastidx)-a+2);
      else
        vcarray_insert(&a2, sval_last, sval_count);
    }
    /* For each small sval, add up the counts */
    for ( ; i < a1.n; i++)
      sum -= arr[i].c;
    /* Merge a1 and a2 into a1.  a2 will be emptied. */
    vcarray_merge(&a1, &a2);
    /* If we've grown too large, use recursive phi to clip. */
    if ( a1.n > NTHRESH ) {
      arr = a1.a;
      if (verbose > 0) printf("clipping small values at a=%lu a1.n=%lu \n", a, a1.n);
#ifdef _OPENMP
      /* #pragma omp parallel for reduction(+: sum) firstprivate(pcache) schedule(dynamic, 16) */
#endif
      for (i = 0; i < a1.n-NTHRESH+NTHRESH/50; i++) {
        UV j = a1.n - 1 - i;
        IV count = arr[j].c;
        if (count != 0) {
          sum += count * phi_small( arr[j].v, a-1, primes, lastidx, &pcache );
          arr[j].c = 0;
        }
      }
    }
    vcarray_remove_zeros(&a1);
    a--;
  }
  phicache_free(&pcache);
  vcarray_destroy(&a2);
  arr = a1.a;
#ifdef _OPENMP
  #pragma omp parallel for reduction(+: sum) schedule(dynamic, 16)
#endif
  for (i = 0; i < a1.n; i++)
    sum += arr[i].c * mapes7( arr[i].v );
  vcarray_destroy(&a1);
  Safefree(primes);
  return (UV) sum;
}


extern UV _XS_meissel_pi(UV n);
/* b = prime_count(isqrt(n)) */
static UV Pk_2_p(UV n, UV a, UV b, const uint32_t* primes, uint32_t lastidx)
{
  UV lastw, lastwpc, i, P2;
  UV lastpc = primes[lastidx];

  /* Ensure we have a large enough base sieve */
  prime_precalc(isqrt(n / primes[a+1]));

  P2 = lastw = lastwpc = 0;
  for (i = b; i > a; i--) {
    UV w = n / primes[i];
    lastwpc = (w <= lastpc) ? bs_prime_count(w, primes, lastidx)
                            : lastwpc + _XS_prime_count(lastw+1, w);
    lastw = w;
    P2 += lastwpc;
  }
  P2 -= ((b+a-2) * (b-a+1) / 2) - a + 1;
  return P2;
}
static UV Pk_2(UV n, UV a, UV b)
{
  UV lastprime = ((b*SIEVE_MULT+1) > 203280221) ? 203280221 : b*SIEVE_MULT+1;
  const uint32_t* primes = generate_small_primes(lastprime);
  UV P2 = Pk_2_p(n, a, b, primes, lastprime);
  Safefree(primes);
  return P2;
}


/* Legendre's method.  Interesting and a good test for phi(x,a), but Lehmer's
 * method is much faster (Legendre: a = pi(n^.5), Lehmer: a = pi(n^.25)) */
UV _XS_legendre_pi(UV n)
{
  UV a, phina;
  if (n < SIEVE_LIMIT)
    return _XS_prime_count(2, n);

  a = _XS_legendre_pi(isqrt(n));
  phina = phi(n, a);
  return phina + a - 1;
}


/* Meissel's method. */
UV _XS_meissel_pi(UV n)
{
  UV a, b, sum;
  if (n < SIEVE_LIMIT)
    return _XS_prime_count(2, n);

  a = _XS_meissel_pi(icbrt(n));       /* a = Pi(floor(n^1/3)) [max    192725] */
  b = _XS_meissel_pi(isqrt(n));       /* b = Pi(floor(n^1/2)) [max 203280221] */

  sum = phi(n, a) + a - 1 - Pk_2(n, a, b);
  return sum;
}

/* Lehmer's method.  This is basically Riesel's Lehmer function (page 22),
 * with some additional code to help optimize it.  */
UV _XS_lehmer_pi(UV n)
{
  UV z, a, b, c, sum, i, j, lastprime, lastpc, lastw, lastwpc;
  const uint32_t* primes = 0; /* small prime cache, first b=pi(z)=pi(sqrt(n)) */
  DECLARE_TIMING_VARIABLES;

  if (n < SIEVE_LIMIT)
    return _XS_prime_count(2, n);

  /* Protect against overflow.  2^32-1 and 2^64-1 are both divisible by 3. */
  if (n == UV_MAX) {
    if ( (n%3) == 0 || (n%5) == 0 || (n%7) == 0 || (n%31) == 0 )
      n--;
    else
      return _XS_prime_count(2,n);
  }

  if (verbose > 0) printf("lehmer %lu stage 1: calculate a,b,c \n", n);
  TIMING_START;
  z = isqrt(n);
  a = _XS_lehmer_pi(isqrt(z));        /* a = Pi(floor(n^1/4)) [max      6542] */
  b = _XS_lehmer_pi(z);               /* b = Pi(floor(n^1/2)) [max 203280221] */
  c = _XS_lehmer_pi(icbrt(n));        /* c = Pi(floor(n^1/3)) [max    192725] */
  TIMING_END_PRINT("stage 1")

  if (verbose > 0) printf("lehmer %lu stage 2: phi(x,a) (z=%lu a=%lu b=%lu c=%lu)\n", n, z, a, b, c);
  TIMING_START;
  sum = phi(n, a) + ((b+a-2) * (b-a+1) / 2);
  TIMING_END_PRINT("phi(x,a)")

  /* We get an array of the first b primes.  This is used in stage 4.  If we
   * get more than necessary, we can use them to speed up some.
   */
  lastprime = b*SIEVE_MULT+1;
  if (lastprime > 203280221) lastprime = 203280221;
  if (verbose > 0) printf("lehmer %lu stage 3: %lu small primes\n", n, lastprime);
  TIMING_START;
  primes = generate_small_primes(lastprime);
  if (primes == 0) croak("Error generating primes.\n");
  lastpc = primes[lastprime];
  TIMING_END_PRINT("small primes")

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
    lastwpc = (w <= lastpc) ? bs_prime_count(w, primes, lastprime)
                            : lastwpc + _XS_prime_count(lastw+1, w);
    lastw = w;
    sum = sum - lastwpc;
    if (i <= c) {
      UV bi = bs_prime_count( isqrt(w), primes, lastprime );
      for (j = i; j <= bi; j++) {
        sum = sum - bs_prime_count(w / primes[j], primes, lastprime) + j - 1;
      }
      /* We could wrap the +j-1 in:  sum += ((bi+1-i)*(bi+i))/2 - (bi-i+1); */
    }
  }
  TIMING_END_PRINT("stage 4")
  Safefree(primes);
  return sum;
}


/* The Lagarias-Miller-Odlyzko method.
 * Naive implementation without optimizations.
 * About the same speed as Lehmer, a bit less memory.
 * A better implementation can be 10-50x faster and much less memory.
 */
UV _XS_LMOS_pi(UV n)
{
  UV n13, a, b, sum, i, j, k, lastprime, P2, S1, S2;
  const uint32_t* primes = 0;  /* small prime cache */
  signed char* mu = 0;   /* moebius to n^1/3 */
  uint32_t*   lpf = 0;   /* least prime factor to n^1/3 */
  cache_t pcache; /* Cache for recursive phi */
  DECLARE_TIMING_VARIABLES;

  if (n < SIEVE_LIMIT)
    return _XS_prime_count(2, n);

  n13 = icbrt(n);                    /* n13 =  floor(n^1/3)  [max    2642245] */
  a = _XS_lehmer_pi(n13);            /* a = Pi(floor(n^1/3)) [max     192725] */
  b = _XS_lehmer_pi(isqrt(n));       /* b = Pi(floor(n^1/2)) [max  203280221] */

  lastprime = b*SIEVE_MULT+1;
  if (lastprime > 203280221) lastprime = 203280221;
  if (lastprime < n13) lastprime = n13;
  primes = generate_small_primes(lastprime);
  if (primes == 0) croak("Error generating primes.\n");

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

  /* Thanks to Kim Walisch for help with the S1+S2 calculations. */
  k = (a < 7) ? a : 7;
  S1 = 0;
  S2 = 0;
  phicache_init(&pcache);
  TIMING_START;
  for (i = 1; i <= n13; i++)
    if (lpf[i] > primes[k])
      /* S1 += mu[i] * phi_small(n/i, k, primes, lastprime, &pcache); */
      S1 += mu[i] * phi(n/i, k);
  TIMING_END_PRINT("S1")

  TIMING_START;
  for (i = k; i+1 < a; i++) {
    uint32_t p = primes[i+1];
    /* TODO: #pragma omp parallel for reduction(+: S2) firstprivate(pcache) schedule(dynamic, 16) */
    for (j = (n13/p)+1; j <= n13; j++)
      if (lpf[j] > p)
        S2 += -mu[j] * phi_small(n / (j*p), i, primes, lastprime, &pcache);
  }
  TIMING_END_PRINT("S2")
  phicache_free(&pcache);
  Safefree(lpf);
  Safefree(mu);

  TIMING_START;
  prime_precalc( (UV) pow(n, 2.9/5.0) );
  P2 = Pk_2_p(n, a, b, primes, lastprime);
  TIMING_END_PRINT("P2")
  Safefree(primes);

  /* printf("S1 = %lu\nS2 = %lu\na  = %lu\nP2 = %lu\n", S1, S2, a, P2); */
  sum = (S1 + S2) + a - 1 - P2;
  return sum;
}

static const unsigned char primes_small[] =
  {0,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
   101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191};
#define NPRIMES_SMALL (sizeof(primes_small)/sizeof(primes_small[0]))

UV _XS_legendre_phi(UV x, UV a) {
  /* For small values, calculate directly */
  if (a < 3)  return (a == 0) ? x : (a == 1) ? x-x/2 : x-x/2-x/3+x/6;
  if (a <= 7) return mapes(x, a);
  /* For large values, do our non-recursive phi */
  if (a > NPRIMES_SMALL) return phi(x,a);
  /* Otherwise, recurse */
  {
    UV i;
    UV sum = mapes7(x);
    for (i = 8; i <= a; i++) {
      uint32_t p = primes_small[i];
      UV xp = x/p;
      if (xp < p) {
        while (x < primes_small[a])
          a--;
        return (sum - a + i - 1);
      }
      sum -= _XS_legendre_phi(xp, i-1);
    }
    return sum;
  }
}


#ifdef PRIMESIEVE_STANDALONE
int main(int argc, char *argv[])
{
  UV n, pi;
  double t;
  const char* method;
  struct timeval t0, t1;

  if (argc <= 1) { printf("usage: %s  <n>  [<method>]\n", argv[0]); return(1); }
  n = strtoul(argv[1], 0, 10);
  if (n < 2) { printf("Pi(%lu) = 0\n", n); return(0); }

  if (argc > 2)
    method = argv[2];
  else
    method = "lehmer";

  gettimeofday(&t0, 0);

  SET_PPS_PARALLEL;
  if      (!strcasecmp(method, "lehmer"))   { pi = _XS_lehmer_pi(n);      }
  else if (!strcasecmp(method, "meissel"))  { pi = _XS_meissel_pi(n);     }
  else if (!strcasecmp(method, "legendre")) { pi = _XS_legendre_pi(n);    }
  else if (!strcasecmp(method, "lmo"))      { pi = _XS_LMOS_pi(n);  }
  else if (!strcasecmp(method, "sieve"))    { pi = _XS_prime_count(2, n); }
  else {
    printf("method must be one of: lehmer, meissel, legendre, lmo, or sieve\n");
    return(2);
  }
  gettimeofday(&t1, 0);
  t = (t1.tv_sec-t0.tv_sec);  t *= 1000000.0;  t += (t1.tv_usec - t0.tv_usec);
  printf("%8s Pi(%lu) = %lu   in %10.5fs\n", method, n, pi, t / 1000000.0);
  return(0);
}
#endif
