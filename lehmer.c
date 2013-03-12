#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Below this size, just sieve. */
#define SIEVE_LIMIT  60000000

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
 * Timings with Perl + MPU with all-serial computation.
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

static UV isqrt(UV n)
{
  UV root;
  if (sizeof(UV) == 8 && n >= 18446744065119617025UL)  return 4294967295UL;
  if (sizeof(UV) == 4 && n >= 4294836225UL)            return 65535UL;
  root = (UV) sqrt((double)n);
  while (root*root > n)  root--;
  while ((root+1)*(root+1) <= n)  root++;
  return root;
}

/* Callback used for creating an array of primes. */
static UV* sieve_array = 0;
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
static UV* generate_small_primes(UV n)
{
  UV* primes;
  double fn = (double)n;
  double flogn  = log(fn);
  double flog2n  = log(flogn);
  UV nth_prime =  /* Dusart 2010 for > 179k, custom for 18-179k */
     (n >= 688383) ? (UV) ceil(fn*(flogn+flog2n-1.0+((flog2n-2.00)/flogn))) :
     (n >= 178974) ? (UV) ceil(fn*(flogn+flog2n-1.0+((flog2n-1.95)/flogn))) :
     (n >= 18)     ? (UV) ceil(fn*(flogn+flog2n-1.0+((flog2n+0.30)/flogn)))
                   : 59;
  New(0, primes, n+1, UV);
  if (primes == 0)
    croak("Can not allocate small primes\n");
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

#include "lehmer.h"
#include "util.h"
#include "cache.h"
#include "sieve.h"

/* Generate an array of n small primes, where the kth prime is element p[k].
 * Remember to free when done. */
static UV* generate_small_primes(UV n)
{
  UV* primes;
  UV  i = 0;
  double fn = (double)n;
  double flogn  = log(fn);
  double flog2n  = log(flogn);
  UV nth_prime =  /* Dusart 2010 for > 179k, custom for 18-179k */
     (n >= 688383) ? (UV) ceil(fn*(flogn+flog2n-1.0+((flog2n-2.00)/flogn))) :
     (n >= 178974) ? (UV) ceil(fn*(flogn+flog2n-1.0+((flog2n-1.95)/flogn))) :
     (n >= 18)     ? (UV) ceil(fn*(flogn+flog2n-1.0+((flog2n+0.30)/flogn)))
                   : 59;

  New(0, primes, n+1, UV);
  if (primes == 0)
    croak("Can not allocate small primes\n");
  primes[0] = 0;
  START_DO_FOR_EACH_PRIME(2, nth_prime) {
    if (i >= n) break;
    primes[++i] = p;
  } END_DO_FOR_EACH_PRIME
  if (i < n)
    croak("Did not generate enough small primes.\n");
  if (verbose > 1) printf("generated %lu small primes, from 2 to %lu\n", i, primes[i]);
  return primes;
}

#endif

static UV icbrt(UV n)
{
  UV root = 0;
  /* int s = BITS_PER_WORD - (BITS_PER_WORD % 3); */
#if BITS_PER_WORD == 32
  int s = 30;
  if (n >= UVCONST(4291015625)) return UVCONST(1625);
#else
  int s = 63;
  if (n >= UVCONST(18446724184312856125)) return UVCONST(2642245);
#endif
#if 0
  /* The integer cube root code is about 30% faster for me */
  root = (UV) pow(n, 1.0/3.0);
  if (root*root*root > n) {
    root--;
    while (root*root*root > n)  root--;
  } else {
    while ((root+1)*(root+1)*(root+1) <= n)  root++;
  }
#else
  for ( ; s >= 0; s -= 3) {
    UV b;
    root += root;
    b = 3*root*(root+1)+1;
    if ((n >> s) >= b) {
      n -= b << s;
      root++;
    }
  }
#endif
  return root;
}


/* Given an array of primes[1..lastprime], return Pi(n) where n <= lastprime.
 * This is actually quite fast, and definitely faster than sieving.  By using
 * this we can avoid caching prime counts and also skip most calls to the
 * segment siever.
 */
static UV bs_prime_count(UV n, UV const* const primes, UV lastprime)
{
  UV i, j;
  if (n < 2)  return 0;
  /* if (n > primes[lastprime])  return _XS_prime_count(2, n); */
  if (n >= primes[lastprime]) {
    if (n == primes[lastprime]) return lastprime;
    croak("called bspc(%lu) with counts up to %lu\n", n, primes[lastprime]);
  }
  i = 1;
  j = lastprime;
  while (i < j) {
    UV mid = (i+j)/2;
    if (primes[mid] <= n)  i = mid+1;
    else                   j = mid;
  }
  return i-1;
}


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

static UV mapes7(UV x) {    /* A tiny bit faster setup for a=7 */
  IV val = x-x/2-x/3-x/5+x/6-x/7+x/10-x/11-x/13+x/14+x/15-x/17+x/21+x/22+x/26
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
  return (UV) val;
}

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
      if (verbose>2) printf("ALLOCing list, size %lu\n", new_size);
      New(0, l->a, new_size, vc_t);
    } else {
      new_size = (UV) (1.5 * l->size);
      if (verbose>2) printf("REALLOCing list %p, new size %lu\n",l->a,new_size);
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
  if (verbose>2) printf("  removed %lu duplicates from b\n", bn - bj);
  bn = bj;

  if (bn == 0) {  /* In case they were all duplicates */
    b->n = 0;
    return;
  }

  /* kn = the final merged size.  All duplicates are gone, so this is exact. */
  kn = an+bn;
  if ((long)a->size < kn) {  /* Make A big enough to hold kn elements */
    UV new_size = (UV) (1.2 * kn);
    if (verbose>2) printf("REALLOCing list %p, new size %lu\n", a->a, new_size);
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


/*
 * The main phi(x,a) algorithm.  In this implementation, it takes under 10%
 * of the total time for the Lehmer algorithm, but is a big memory consumer.
 */

static UV phi(UV x, UV a)
{
  UV i, val, sval;
  UV sum = 0;
  IV count;
  const UV* primes;
  vcarray_t a1, a2;
  vc_t* arr;

  if (a == 1)  return ((x+1)/2);
  if (a <= 7)  return mapes(x, a);

  primes = generate_small_primes(a+1);
  if (primes == 0)
    croak("Could not generate primes for phi(%lu,%lu)\n", x, a);
  if (x < primes[a+1])  { Safefree(primes); return (x > 0) ? 1 : 0; }

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
      if (count == 0) continue;      /* Skip if count = 0 */
      val  = arr[i].v;
      sval = val / primea;
      if (sval < primea) break;      /* stop inserting into a2 if small */
      if (sval != sval_last) {       /* non-merged value.  Insert into a2 */
        if (sval_last != 0)
          vcarray_insert(&a2, sval_last, sval_count);
        sval_last = sval;
        sval_count = 0;
      }
      sval_count -= count;           /* Accumulate count for this sval */
    }
    if (sval_last != 0)              /* Insert the last sval */
      vcarray_insert(&a2, sval_last, sval_count);
    /* For each small sval, add up the counts */
    for ( ; i < a1.n; i++) 
      sum -= arr[i].c;
    /* Merge a1 and a2 into a1.  a2 will be emptied. */
    vcarray_merge(&a1, &a2);
    a--;
  }
  vcarray_destroy(&a2);
  if (a != 7) croak("final loop is set for a=7, a = %lu\n", a);
  arr = a1.a;
  for (i = 0; i < a1.n; i++) {
    count = arr[i].c;
    if (count != 0)
      sum += count * mapes7( arr[i].v );
  }
  vcarray_destroy(&a1);
  Safefree(primes);
  return (UV) sum;
}





/* Legendre's method.  Interesting and a good test for phi(x,a), but Lehmer's
 * method is much faster (Legendre: a = pi(n^.5), Lehmer: a = pi(n^.25)) */
UV _XS_legendre_pi(UV n)
{
  UV a;
  if (n < SIEVE_LIMIT)
    return _XS_prime_count(2, n);

  a = _XS_legendre_pi(isqrt(n));

  return phi(n, a) + a - 1;
}


/* Meissel's method. */
UV _XS_meissel_pi(UV n)
{
  UV a, b, c, sum, i, lastprime, lastpc, lastw, lastwpc;
  const UV* primes = 0;  /* small prime cache */
  DECLARE_TIMING_VARIABLES;
  if (n < SIEVE_LIMIT)
    return _XS_prime_count(2, n);

  if (verbose > 0) printf("meissel %lu stage 1: calculate a,b,c \n", n);
  TIMING_START;
  a = _XS_meissel_pi(icbrt(n));        /* a = floor(n^1/3) */
  b = _XS_meissel_pi(isqrt(n));        /* b = floor(n^1/2) */
  c = a;                               /* c = a            */
  TIMING_END_PRINT("stage 1")

  if (verbose > 0) printf("meissel %lu stage 2: phi(x,a) (a=%lu b=%lu c=%lu)\n", n, a, b, c);
  TIMING_START;
  sum = phi(n, a) + ((b+a-2) * (b-a+1) / 2);
  if (verbose > 0) printf("phi(%lu,%lu) = %lu.  sum = %lu\n", n, a, sum - ((b+a-2) * (b-a+1) / 2), sum);
  TIMING_END_PRINT("phi(x,a)")

  lastprime = b*SIEVE_MULT+1;
  if (verbose > 0) printf("meissel %lu stage 3: %lu small primes\n", n, lastprime);
  TIMING_START;
  primes = generate_small_primes(lastprime);
  if (primes == 0) croak("Error generating primes.\n");
  lastpc = primes[lastprime];
  TIMING_END_PRINT("small primes")

  prime_precalc(isqrt(n / primes[a+1]));
  prime_precalc( (UV) pow(n, 2.0/5.0) );  /* Sieve more for speed */

  if (verbose > 0) printf("meissel %lu stage 4: loop %lu to %lu, pc to %lu\n", n, a+1, b, n/primes[a+1]);
  TIMING_START;
  /* Reverse the i loop so w increases.  Count w in segments. */
  lastw = 0;
  lastwpc = 0;
  for (i = b; i > a; i--) {
    UV w = n / primes[i];
    lastwpc = (w <= lastpc) ? bs_prime_count(w, primes, lastprime)
                            : lastwpc + _XS_prime_count(lastw+1, w);
    lastw = w;
    sum = sum - lastwpc;
  }
  TIMING_END_PRINT("stage 4")
  Safefree(primes);
  return sum;
}

/* Lehmer's method.  This is basically Riesel's Lehmer function (page 22),
 * with some additional code to help optimize it.  */
UV _XS_lehmer_pi(UV n)
{
  UV z, a, b, c, sum, i, j, lastprime, lastpc, lastw, lastwpc;
  const UV* primes = 0; /* small prime cache, first b=pi(z)=pi(sqrt(n)) */
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
  a = _XS_lehmer_pi(isqrt(z));         /* a = floor(n^1/4) */
  b = _XS_lehmer_pi(z);                /* b = floor(n^1/2) */
  c = _XS_lehmer_pi(icbrt(n));         /* c = floor(n^1/3) */
  TIMING_END_PRINT("stage 1")

  if (verbose > 0) printf("lehmer %lu stage 2: phi(x,a) (z=%lu a=%lu b=%lu c=%lu)\n", n, z, a, b, c);
  TIMING_START;
  sum = phi(n, a) + ((b+a-2) * (b-a+1) / 2);
  TIMING_END_PRINT("phi(x,a)")

  /* We get an array of the first b primes.  This is used in stage 4.  If we
   * get more than necessary, we can use them to speed up some.
   */
  lastprime = b*SIEVE_MULT+1;
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


UV _XS_LMO_pi(UV n)
{
  UV a, b, sum, i, lastprime, lastpc, lastw, lastwpc;
  UV n13, n12, n23;
  IV S1;
  UV S2, P2;
  const UV* primes = 0;  /* small prime cache */
  char* mu = 0;          /* moebius to n^1/3 */
  UV*   lpf = 0;         /* least prime factor to n^1/3 */
  DECLARE_TIMING_VARIABLES;
  if (n < SIEVE_LIMIT)
    return _XS_prime_count(2, n);

  if (verbose > 0) printf("LMO %lu stage 1: calculate pi(n^1/3) \n", n);
  TIMING_START;
  n13 = icbrt(n);
  n12 = isqrt(n);
  n23 = (UV) (pow(n, 2.0/3.0)+0.01);
  a = _XS_lehmer_pi(n13);
  b = _XS_lehmer_pi(n12);
  TIMING_END_PRINT("stage 1")

  lastprime = b*SIEVE_MULT+1;
  if (verbose > 0) printf("LMO %lu stage 2: %lu small primes\n", n, lastprime);
  TIMING_START;
  primes = generate_small_primes(lastprime);
  if (primes == 0) croak("Error generating primes.\n");
  lastpc = primes[lastprime];
  TIMING_END_PRINT("small primes")

  if (verbose > 0) printf("LMO %lu stage 3: calculate mu/lpf to %lu\n", n, a);
  TIMING_START;
  /* We could call MPU's:
   *    mu = _moebius_range(0, n13+1)
   * but (1) it's a bit slower (something to be addressed), and (2) we will
   * do the least prime factor calculation at the same time.
   */
  New(0, mu, n13+1, char);
  memset(mu, 1, sizeof(char) * (n13+1));
  New(0, lpf, n13+1, UV);
  memset(lpf, 0, sizeof(UV) * (n13+1));
  mu[0] = 0;
  for (i = 1; i <= a; i++) {
    UV primei = primes[i];
    UV j, isquared;
    for (j = primei; j <= n13; j += primei) {
      mu[j] = -mu[j];
      if (lpf[j] == 0) lpf[j] = primei;
    }
    isquared = primei * primei;
    for (j = isquared; j <= n13; j += isquared)
      mu[j] = 0;
  }
  /* for (i = 0; i <= n13; i++) { printf("mu %lu %ld\n", i, (IV)mu[i]); } */
  TIMING_END_PRINT("mu")

  if (verbose > 0) printf("LMO %lu stage 4: calculate S1 (%lu)\n", n, n13);
  TIMING_START;
  S1 = 0;
  for (i = 1; i <= n13; i++)
    if (mu[i] != 0)
      S1 += mu[i] * (IV) (n/i);
  TIMING_END_PRINT("S1")
  if (verbose > 0) printf("LMO %lu stage 4: S1 = %ld\n", n, S1);

  S2 = 0;
  /* TODO... */

  Safefree(mu);
  Safefree(lpf);

  prime_precalc(isqrt(n / primes[a+1]));
  if (verbose > 0) printf("LMO %lu stage 5: P2 loop %lu to %lu, pc to %lu\n", n, a+1, b, n/primes[a+1]);
  TIMING_START;
  P2 = 0;
  /* Reverse the i loop so w increases.  Count w in segments. */
  lastw = 0;
  lastwpc = 0;
  for (i = b; i > a; i--) {
    UV w = n / primes[i];
    lastwpc = (w <= lastpc) ? bs_prime_count(w, primes, lastprime)
                            : lastwpc + _XS_prime_count(lastw+1, w);
    lastw = w;
    P2 += lastwpc;
  }
  P2 -= ((b+a-2) * (b-a+1) / 2) - a + 1;
  TIMING_END_PRINT("P2")
  if (verbose > 0) printf("LMO %lu stage 5: P2 = %lu\n", n, P2);
  Safefree(primes);
  sum = P2 + S1 + S2;
  return sum;
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
  else if (!strcasecmp(method, "lmo"))      { pi = _XS_LMO_pi(n);    }
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
