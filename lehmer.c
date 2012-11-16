#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lehmer.h"
#include "util.h"
#include "cache.h"
#include "sieve.h"

/*****************************************************************************
 *
 * Lehmer prime counting utility.  Calculates pi(x), count of primes <= x.
 *
 * Copyright (c) 2012 Dana Jacobsen (dana@acm.org).
 * This is free software; you can redistribute it and/or modify it under
 * the same terms as the Perl 5 programming language system itself.
 *
 * This file is part of the Math::Prime::Util Perl module.  It relies on two
 * features: (1) a sieve to generate very small primes which could be as
 * simple as a three-line SoE, and (2) _XS_prime_count(low, high) is expected
 * to return the prime count within the segment low to high inclusive.  These
 * are used for the relatively small segments in the final stage, but making
 * it fast is important for the final performance.  The primesieve package is
 * one source of excellent routines for either task.
 *
 * The phi(x,a) calculation is unique, to the best of my knowledge.  It keeps
 * a heap of all x values + signed counts for the given 'a' value, and walks
 * 'a' down until it is small enough to calculate directly (either with Mapes
 * or using a calculated table using the primorial/totient method).  This
 * is relatively fast and low memory compared to many other solutions.  As with
 * all Lehmer-Meissel-Legendre algorithms, memory use will be a constraint
 * with large values of x.
 *
 * Calculating pi(10^11) is done in under 1 second on my computer.  pi(10^14)
 * takes about 1 minute, pi(10^16) in under an hour.  Compared to Thomas
 * R. Nicely's pix4 program, this one is about 3-4x faster and uses 3x less
 * memory.  However those comparisons were done without the big pre-computed
 * prime count tables pix4 can use to speed up the final stage.
 *
 *    n     phi(x,a) mem/time  |  stage 4 mem/time  | total time
 *   10^17   5988MB  1244.11   |   2688MB
 *   10^16   1877MB   254.83   |   945MB  1478.9    |  28m 51.4s
 *   10^15    754MB    50.50   |   392MB   224.2    |   4m 33.8s
 *   10^14    260MB     9.879  |   217MB    38.06   |    47.83s
 *   10^13              1.852  |   162MB     6.987  |     8.864s
 *   10^12              0.363  |             1.358  |     1.768s
 *   10^11              0.064  |             0.276  |     0.381s
 *   10^10  xxxxxMB     0.011  | xxxxxMB     0.056  |     0.105s
 *
 * Reference: "Prime numbers and computer methods for factorization", Riesel.
 */

static int const verbose = 0;
#define DO_TIMING 0
#ifdef DO_TIMING
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

/* Use Mapes' method to calculate phi(x,a) for small a.  This is really
 * convenient and a little Perl script will spit this code out for whatever
 * limit we select.  It gets unwieldy however, and using the primorial/totient
 * method looks faster for the final output.
 */
static IV mapes(IV x, UV a)
{
  IV val;
  if (a == 0)  return x;
  if (a == 1)  return x-x/2;
  val = x-x/2-x/3+x/6;
  if (a >= 3) val += -x/5+x/10+x/15-x/30;
  if (a >= 4) val += -x/7+x/14+x/21-x/42+x/35-x/70-x/105+x/210;
  if (a >= 5) val += -x/11+x/22+x/33-x/66+x/55-x/110-x/165+x/330+x/77-x/154-x/231+x/462-x/385+x/770+x/1155-x/2310;
  if (a >= 6) val += -x/13+x/26+x/39-x/78+x/65-x/130-x/195+x/390+x/91-x/182-x/273+x/546-x/455+x/910+x/1365-x/2730+x/143-x/286-x/429+x/858-x/715+x/1430+x/2145-x/4290-x/1001+x/2002+x/3003-x/6006+x/5005-x/10010-x/15015+x/30030;
  if (a >= 7) val += -x/17+x/34+x/51-x/102+x/85-x/170-x/255+x/510+x/119-x/238-x/357+x/714-x/595+x/1190+x/1785-x/3570+x/187-x/374-x/561+x/1122-x/935+x/1870+x/2805-x/5610-x/1309+x/2618+x/3927-x/7854+x/6545-x/13090-x/19635+x/39270+x/221-x/442-x/663+x/1326-x/1105+x/2210+x/3315-x/6630-x/1547+x/3094+x/4641-x/9282+x/7735-x/15470-x/23205+x/46410-x/2431+x/4862+x/7293-x/14586+x/12155-x/24310-x/36465+x/72930+x/17017-x/34034-x/51051+x/102102-x/85085+x/170170+x/255255-x/510510;
  return val;
}

/******************************************************************************/
/*   Modified heap for manipulating our UV value / IV count pairs             */
/******************************************************************************/

/* TODO: This should be cleaned up */

/* heap of values and counts, stored by value */
typedef struct {
  UV v;
  IV c;
} vc_t;

typedef struct {
  vc_t* array;
  UV    array_size;
  UV    N;
  UV    small_limit;
  IV*   small_array;
  UV    Nsmall;
  int   ptr_small;
} heap_t;

static void heap_insert(heap_t* h, vc_t v)
{
  UV k;
  UV val = v.v;
  if (v.c == 0)
    return;
  if (val < h->small_limit) {
    if (h->small_array[val] == 0) ++(h->Nsmall);
    h->small_array[val] += v.c;
    if (h->small_array[val] == 0) --(h->Nsmall);
    else if (h->ptr_small < (int)val) h->ptr_small = (int) val;
// printf("small array insert %lu count %ld, coalesce to %ld, Nsmall = %lu, ptr %d\n", v.v, v.c, h->small_array[val], h->Nsmall, h->ptr_small);
    return;
  }
  UV n = ++(h->N);
  vc_t* a = h->array;
  if (n >= h->array_size) {
    UV new_size = (h->array_size == 0) ? 20000 : 1.5 * h->array_size;
    if (verbose>2) printf("REALLOCing %p, new size %lu\n", h->array, new_size);
    h->array = realloc( h->array, new_size * sizeof(vc_t) );
    if (h->array == 0)
      croak("could not allocate heap\n");
    a = h->array;
    h->array_size = new_size-1;
    a[0].v = UV_MAX;  a[0].c = 0;
  }
  a[n] = v;
  k = n;
  while (a[k/2].v <= v.v) {  /* upheap */
    a[k] = a[k/2];
    k /= 2;
  }
  a[k] = v;
}

static vc_t heap_remove(heap_t* h)
{
  UV k;
  UV n = h->N;
  vc_t* a = h->array;
  vc_t v, top;
  if (n == 0)
    croak("removing from empty heap\n");
  top = a[1];
  v = a[1] = a[n];
  n = --(h->N);                     /* downheap */
  k = 1;
  while (k <= n/2) {
    UV j = k+k;
    if (j < n && a[j].v < a[j+1].v)  j++;
    if (v.v >= a[j].v) break;
    a[k] = a[j];
    k = j;
  }
  a[k] = v;
  return top;
}
static vc_t heap_remove_coalesce(heap_t* h)
{
  vc_t v;
  if (h->N == 0) {
    int i;
    if (h->Nsmall == 0) croak("remove from empty heap\n");
    for (i = h->ptr_small; i >= 0; i--) {
      if (h->small_array[i] != 0) {
        v.v = (UV) i;
        v.c = h->small_array[i];
        h->small_array[i] = 0;
        h->ptr_small = i-1;
        --(h->Nsmall);
 //printf("small array returning %lu count %ld\n", v.v, v.c);
        return v;
      }
    }
    croak("walked off small array\n");
  }
  v = heap_remove(h);
  /* get rest of entries with same value */
  while (h->N > 0 && h->array[1].v == v.v) {
    v.c += h->array[1].c;
    (void) heap_remove(h);
  }
  return v;
}
static heap_t heap_create(UV small_size)
{
  heap_t h;
  h.array = 0;
  h.array_size = 0;
  h.N = 0;
  h.small_limit = small_size;
  h.small_array = 0;
  if (small_size > 0) {
    if (verbose>1)printf("creating small array of size %lu\n", small_size);
    h.small_array = (IV*) calloc( small_size, sizeof(IV) );
  }
  h.Nsmall = 0;
  h.ptr_small = -1;
  return h;
}
static void heap_destroy(heap_t* h)
{
  if (h->array != 0) {
    if (verbose > 2) printf("FREE %p\n", h->array);
    free(h->array);
  }
  h->array = 0;
  h->array_size = 0;
  h->N = 0;
  if (h->small_array != 0) {
    free(h->small_array);
  }
  h->small_array = 0;
  h->Nsmall = 0;
  h->ptr_small = -1;
}
static void heap_empty(heap_t* h)
{
  h->N = 0;
  h->Nsmall = 0;
  h->ptr_small = -1;
  if (h->small_limit > 0)
    memset(h->small_array, 0, h->small_limit * sizeof(IV));
}

#define heap_not_empty(h) ((h).N > 0 || (h).Nsmall > 0)

/******************************************************************************/



/*
 * The main phi(x,a) algorithm.  In this implementation, it takes about 25%
 * of the total time for the Lehmer algorithm, but it is by far the most memory
 * consuming section.
 */

static UV phi(UV x, UV a, UV* primes)
{
  UV i;
  IV  sum = 0;
  UV  primea = primes ? primes[a+1] : _XS_nth_prime(a+1);
  if (x < primea)  return (x > 0) ? 1 : 0;
  if (a == 1)      return ((x+1)/2);
  if (a <= 7)      return mapes(x, a);

  primea = primes ? primes[a] : _XS_prev_prime(primea);

  heap_t h1 = heap_create(20 * x/primea/primea/primea);
  heap_t h2 = heap_create(20 * x/primea/primea/primea);
  vc_t v;

  v.v = x;  v.c = 1;
  heap_insert(&h1, v);
  while (a > 6) {
    heap_empty(&h2);
    while ( heap_not_empty(h1) ) {
      UV sval;
      v = heap_remove_coalesce(&h1);
      if (v.c == 0)
        continue;
      heap_insert(&h2, v);
      sval = v.v / primea;
      if (sval >= primea) {
        v.v = sval;
        v.c = -v.c;
        heap_insert(&h2, v);
      } else {
        sum -= v.c;
      }
    }
    { heap_t t = h1; h1 = h2; h2 = t; }
    a--;
    primea = primes ? primes[a] : _XS_prev_prime(primea);
  }
  heap_destroy(&h2);
  while ( heap_not_empty(h1) ) {
    v = heap_remove_coalesce(&h1);
    if (v.c != 0)
      sum += v.c * mapes(v.v, a);  /* This could be faster */
  }
  heap_destroy(&h1);
  return sum;
}


/* Legendre's method.  Interesting and a good test for phi(x,a), but Lehmer's
 * method is much faster (Legendre: a = pi(n^.5), Lehmer: a = pi(n^.25)) */
UV _XS_legendre_pi(UV n)
{
  if (n < 30000)
    return _XS_prime_count(2, n);

  UV a = _XS_legendre_pi(sqrt(n));
  prime_precalc(a);
  return phi(n, a, 0) + a - 1;
}

/* Lehmer's method.  See Riesel for basic program. */

UV _XS_lehmer_pi(UV n)
{
  DECLARE_TIMING_VARIABLES;
  UV z, a, b, c, sum, i, j, lastw, lastwpc;
  UV* primea = 0;   /* small prime cache, first b=pi(z)=pi(sqrt(n)) primes */
  UV* pca = 0;      /* small prime count cache, first z=sqrt(n) numbers */
  if (n < 1000000)
    return _XS_prime_count(2, n);

  if (verbose > 0) printf("lehmer %lu stage 1: calculate a,b,c \n", n);
  TIMING_START;
  z = sqrt((double)n+0.5);
  a = _XS_lehmer_pi(sqrt((double)z)+0.5);
  b = _XS_lehmer_pi(z);
  c = _XS_lehmer_pi(pow((double)n, 1.0/3.0)+0.5);
  TIMING_END_PRINT("stage 1")


  if (verbose > 0) printf("lehmer %lu stage 2: %lu small primes, %lu counts\n", n, b, z);
  TIMING_START;
  { /* Calculate the first b primes */
    const unsigned char* sieve;
    if (get_prime_cache(z, &sieve) < z) {
      release_prime_cache(sieve);
      croak("Could not generate sieve for %"UVuf, z);
    }
    primea = (UV*) malloc( (b+4) * sizeof(UV) );
    if (primea == 0)
      croak("Can not allocate small primes\n");
    primea[0] = 0;
    primea[1] = 2;
    primea[2] = 3;
    primea[3] = 5;
    i = 3;
    START_DO_FOR_EACH_SIEVE_PRIME( sieve, 7, z ) {
      if (i > b) break;
      primea[++i] = p;
    } END_DO_FOR_EACH_SIEVE_PRIME
    release_prime_cache(sieve);
    if (i < b)
      croak("Did not generate enough small primes.  Bug in Lehmer code.\n");
    if (verbose > 1) printf("generated first %lu small primes, from 2 to %lu\n", i, z);
  }
  TIMING_END_PRINT("small primes")


  if (verbose > 0) printf("lehmer %lu stage 3: phi(x,a) (z=%lu a=%lu b=%lu c=%lu)\n", n, z, a, b, c);
  TIMING_START;
  sum = phi(n, a, primea) + ( (b+a-2) * (b-a+1) / 2);
  TIMING_END_PRINT("phi(x,a)")


  TIMING_START;
  { /* Generate primecounts up to z used by stage 4. */
    UV count = 0;
    pca = (UV*) malloc( (z+4) * sizeof(UV) );
    if (pca == 0)
      croak("Can not allocate small prime counts\n");
    for (i = 0; i <= z && count+1 <= b; i++) {
      if (i >= primea[count+1])
        count++;
      pca[i] = count;
    }
    if (verbose > 1) printf("generated first %lu prime counts, from 2 to %lu\n", count, z);
  }
 TIMING_END_PRINT("prime counts")
printf("precalculated counts to %lu\n", z);

  /* Ensure we have the base sieve for prime_count ( n/primea[i] ). */
  /* This is about 75k for n=10^13, 421k for n=10^15, 2.4M for n=10^17 */
  prime_precalc( sqrt( n / primea[a+1] ) );

 /* TODO:
  *  (1) generate 2b primes instead of b
  *  (2) generate counts for all 2b primes
  *  (3) store counts on odd numbers only.
  *  (4) in loop below, check if we have the result in pca and use it.
  *  This will cut in half the prime_count calls below, at a little memory
  *  cost (almost entirely mitigated by change #3)
  */

  if (verbose > 0) printf("lehmer %lu stage 4: loop %lu to %lu, pc to %lu\n", n, a+1, b, n/primea[a+1]);
  TIMING_START;
  /* Reverse the i loop so w increases.  Count w in segments. */
  lastw = 0;
  lastwpc = 0;
  for (i = b; i >= a+1; i--) {
    UV w = n / primea[i];
    lastwpc += _XS_prime_count(lastw+1, w);
    lastw = w;
    sum = sum - lastwpc;
    if (i <= c) {
      UV bi = pca[(UV) (sqrt(w) + 0.5)];
      for (j = i; j <= bi; j++) {
        sum = sum - pca[w / primea[j]] + j - 1;
      }
    }
  }
  TIMING_END_PRINT("stage 4")
  if (pca != 0)
    free(pca);
  if (primea != 0)
    free(primea);
  return sum;
}
