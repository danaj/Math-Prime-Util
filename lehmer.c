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
 * takes under 1 minute, pi(10^16) in a half hour.  Compared with Thomas
 * R. Nicely's pix4 program, this one is 3-5x faster and uses 2-3x less memory.
 *
 *    n     phi(x,a) mem/time  |  stage 4 mem/time  | total time
 *   10^17   4953MB    871.14   |   2988MB  9911.9    | 179m 37.5s
 *   10^16   1436MB    168.02   |   901MB   1195.7    |  22m 45.6s
 *   10^15    432MB     31.34   |   394MB    165.6    |   3m 17.5s
 *   10^14    203MB      5.509  |   223MB     25.96   |    31.69s
 *   10^13               0.949  |   165MB      4.284  |     5.336s
 *   10^12               0.174  |              0.755  |     0.990s
 *   10^11               0.034  |              0.138  |     0.213s
 *   10^10               0.007  |              0.025  |     0.064s
 *
 * Reference: Hans Riesel, "Prime Numbers and Computer Methods for
 * Factorization", 2nd edition, 1994.
 */

static int const verbose = 0;
#define DO_TIMING 0
#if DO_TIMING
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
 * limit we select.  It gets unwieldy with large a values.
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

static IV mapes7(IV x) {    /* A tiny bit faster setup for a=7 */
  IV val = x-x/2-x/3-x/5+x/6-x/7+x/10-x/11-x/13+x/14+x/15-x/17+x/21+x/22+x/26
          -x/30+x/33+x/34+x/35+x/39-x/42+x/51+x/55+x/65-x/66-x/70+x/77-x/78
          +x/85+x/91-x/102-x/105-x/110+x/119-x/130+x/143-x/154-x/165-x/170
          -x/182+x/187-x/195+x/210+x/221-x/231-x/238-x/255-x/273-x/286+x/330
          -x/357-x/374-x/385+x/390-x/429-x/442-x/455+x/462+x/510+x/546-x/561
          -x/595-x/663+x/714;
  if (x >= 715) {
    val += -x/715+x/770+x/858+x/910-x/935-x/1001-x/1105+x/1122+x/1155+x/1190
           -x/1309+x/1326+x/1365+x/1430-x/1547+x/1785+x/1870+x/2002+x/2145
           +x/2210-x/2310-x/2431+x/2618-x/2730+x/2805+x/3003+x/3094+x/3315
           -x/3570+x/3927-x/4290+x/4641+x/4862+x/5005-x/5610-x/6006+x/6545
           -x/6630+x/7293+x/7735-x/7854;
    if (x >= 9282)
      val += -x/9282-x/10010+x/12155-x/13090-x/14586-x/15015-x/15470+x/17017
             -x/19635-x/23205-x/24310+x/30030-x/34034-x/36465+x/39270+x/46410
             -x/51051+x/72930-x/85085+x/102102+x/170170+x/255255-x/510510;
  }
  return val;
}

/******************************************************************************/
/*   Modified heap for manipulating our UV value / IV count pairs             */
/******************************************************************************/

/*
 * This is a heap augmented with a small array.  We store values and signed
 * counts, where all counts for the same value are summed.  An easy way to
 * do this in Perl/Python is a hash.  In plain C, I don't believe this is the
 * best solution.  A heap can be implemented both easily and very efficiently
 * (using a linear array), and as we pull items off the heap we can combine
 * all similar values.
 *
 * Below some threshold value ('small_limit') the items become dense.  That is,
 * not only are the values small but we have many items in that range.  Hence
 * the small array augmentation.  All values below the threshold are just put
 * directly into an array.  This not only handles them a little faster but
 * helps reduce the heap size a bit, as we don't put any repeated values in
 * the heap for the small items.  Since they're dense in this range, we can
 * do a linear scan to find the next non-zero count.
 *
 * An ideal data structure for our purpose would coalesce values on insertion,
 * and would allow operating in place (so we could retrieve all our items and
 * add new items as we go, without them appearing on this scan).  The former
 * is possible using an ordered list or a balanced tree.  I don't know how we
 * would achieve the latter.  The point being that we're pulling items off of
 * h1 and adding it (plus possibly a new item) to h2, so ideally we'd manage
 * to use that space freed up by h1.
 */

typedef struct {
  UV v;
  IV c;
} vc_t;

typedef struct {
  UV    small_limit;    /* small count array: size */
  UV    small_N;        /* small count array: number of non-zero elements */
  int   small_ptr;      /* small count array: index of largest non-zero value */
  UV    array_size;     /* heap: allocated size in elements */
  UV    N;              /* heap: number of elements */
  IV*   small_array;    /* small count array data */
  vc_t* array;          /* heap data */
} heap_t;

static void heap_insert(heap_t* h, UV val, IV count)
{
  UV n;
  vc_t* a;
  if (val < h->small_limit) {
    IV* saptr = h->small_array + val;
    if (*saptr == 0)
      ++(h->small_N);
    *saptr += count;
    if (*saptr == 0)
      --(h->small_N);
    else if (h->small_ptr < (int)val)
      h->small_ptr = (int) val;
    return;
  }
  n = ++(h->N);
  a = h->array;
  if (n >= h->array_size) {
    UV new_size;
    if (h->array_size == 0) {
      new_size = (h->small_limit <= (20000/2)) ? 20000 : 2*h->small_limit;
      if (verbose>2) printf("ALLOCing heap, size %lu\n", new_size);
      New(0, h->array, new_size, vc_t);
    } else {
      new_size = 1.5 * h->array_size;
      if (verbose>2) printf("REALLOCing heap %p, new size %lu\n", h->array, new_size);
      Renew( h->array, new_size, vc_t );
    }
    if (h->array == 0) croak("could not allocate heap\n");
    a = h->array;
    h->array_size = new_size-1;
    a[0].v = UV_MAX;  a[0].c = 0;
  }
  while (a[n/2].v <= val) {  /* upheap */
    a[n] = a[n/2];
    n /= 2;
  }
  a[n].v = val;
  a[n].c = count;
}
static void heap_remove(heap_t* h, UV* val, IV* count)
{
  if (h->N == 0) {
    /* Search small_array for a non-zero count from small_ptr down */
    IV* saptr = h->small_array + h->small_ptr;
    if (h->small_N == 0) croak("remove from empty heap\n");
    while (!*saptr)
      saptr--;
    if (saptr < h->small_array) croak("walked off small array\n");
    *val = (saptr - h->small_array);
    *count = *saptr;
    *saptr = 0;
    h->small_ptr = *val-1;
    --(h->small_N);
    return;
  } else {
    vc_t* a = h->array;
    UV ival, k, n = h->N;
    *val = a[1].v;
    *count = 0;
    do {
      *count += a[1].c;
      /* remove top element */
      ival = a[n--].v;
      k = 1;
      while (k <= n/2) {
        UV j = k+k;
        if (j < n && a[j].v < a[j+1].v)  j++;
        if (ival >= a[j].v) break;
        a[k] = a[j];
        k = j;
      }
      a[k] = a[n+1];
    } while (n > 0 && a[1].v == *val);
    h->N = n;
  }
}
static heap_t heap_create(UV small_size)
{
  heap_t h;
  h.array = 0;
  h.array_size = 0;
  h.small_limit = small_size;
  h.small_array = 0;
  if (small_size > 0) {
    if (verbose>1)printf("creating small array of size %lu\n", small_size);
    Newz(0, h.small_array, small_size, IV);
  }
  h.N = 0;
  h.small_N = 0;
  h.small_ptr = -1;
  return h;
}
static void heap_destroy(heap_t* h)
{
  if (h->array != 0) {
    if (verbose > 2) printf("FREE heap %p\n", h->array);
    Safefree(h->array);
  }
  h->array = 0;
  h->array_size = 0;
  h->N = 0;
  if (h->small_array != 0)
    Safefree(h->small_array);
  h->small_array = 0;
  h->small_N = 0;
  h->small_ptr = -1;
}
#define heap_not_empty(h)  ((h).N > 0 || (h).small_N > 0)

/******************************************************************************/



/*
 * The main phi(x,a) algorithm.  In this implementation, it takes about 15%
 * of the total time for the Lehmer algorithm, but it is by far the most
 * memory consuming part.
 */

static UV phi(UV x, UV a, UV const* const primes)
{
  heap_t h1, h2;
  UV val;
  IV count, sum = 0;
  UV primea = primes ? primes[a+1] : _XS_nth_prime(a+1);
  if (x < primea)  return (x > 0) ? 1 : 0;
  if (a == 1)      return ((x+1)/2);
  if (a <= 7)      return mapes(x, a);

  primea = primes ? primes[a] : _XS_prev_prime(primea);

  h1 = heap_create(a * 1000);
  h2 = heap_create(a * 1000);

  heap_insert(&h1, x, 1);

  while (a > 7) {
    if (heap_not_empty(h2)) croak("h2 heap isn't empty.");
    while ( heap_not_empty(h1) ) {
      UV sval;
      heap_remove(&h1, &val, &count);
      if (count == 0)
        continue;
      heap_insert(&h2, val, count);
      sval = val / primea;
      if (sval >= primea) {
        heap_insert(&h2, sval, -count);
      } else {
        sum -= count;
      }
    }
    { heap_t t = h1; h1 = h2; h2 = t; }
    a--;
    primea = primes ? primes[a] : _XS_prev_prime(primea);
  }
  heap_destroy(&h2);
  if (a != 7) croak("final loop is set for a=7, a = %lu\n", a);
  while ( heap_not_empty(h1) ) {
    heap_remove(&h1, &val, &count);
    if (count != 0)
      sum += count * mapes7(val);
  }
  heap_destroy(&h1);
  return (UV) sum;
}


static UV* generate_small_primes(UV *n)
{
  const unsigned char* sieve;
  UV* primea;
  UV  nth_prime;
  UV  i;

  /* Dusart 1999 bound */
  nth_prime = (*n <= 10) ? 29 : *n * ( log(*n) + log(log(*n)) ) + 1;

  if (get_prime_cache(nth_prime, &sieve) < nth_prime) {
    release_prime_cache(sieve);
    croak("Could not generate sieve for %"UVuf, nth_prime);
  }
  New(0, primea, *n+4, UV);
  if (primea == 0)
    croak("Can not allocate small primes\n");
  primea[0] = 0; primea[1] = 2; primea[2] = 3; primea[3] = 5;
  i = 3;
  START_DO_FOR_EACH_SIEVE_PRIME( sieve, 7, nth_prime ) {
    if (i >= *n) break;
    primea[++i] = p;
  } END_DO_FOR_EACH_SIEVE_PRIME
  release_prime_cache(sieve);
  if (i < *n)
    croak("Did not generate enough small primes.\n");
  if (verbose > 1) printf("generated %lu small primes, from 2 to %lu\n", i, primea[i]);
  return primea;
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
  if (n > primes[lastprime])
    croak("called bspc(%lu) with counts up to %lu\n", n, primes[lastprime]);
  i = 1;
  j = lastprime;
  while (i < j) {
    UV mid = (i+j)/2;
    if (primes[mid] <= n)  i = mid+1;
    else                   j = mid;
  }
  return i-1;
}


/* Legendre's method.  Interesting and a good test for phi(x,a), but Lehmer's
 * method is much faster (Legendre: a = pi(n^.5), Lehmer: a = pi(n^.25)) */
UV _XS_legendre_pi(UV n)
{
  UV a;
  if (n < 30000)
    return _XS_prime_count(2, n);

  a = _XS_legendre_pi(sqrt(n));
  prime_precalc(a);
  return phi(n, a, 0) + a - 1;
}

/* Lehmer's method.  This is basically Riesel's Lehmer function (page 22),
 * with some additional code to help optimize it.
 */

UV _XS_lehmer_pi(UV n)
{
  UV z, a, b, c, sum, i, j, lastprimea, lastpc, lastw, lastwpc;
  UV* primea = 0;   /* small prime cache, first b=pi(z)=pi(sqrt(n)) primes */
  DECLARE_TIMING_VARIABLES;

  if (n < 1000000)
    return _XS_prime_count(2, n);

  if (verbose > 0) printf("lehmer %lu stage 1: calculate a,b,c \n", n);
  TIMING_START;
  z = sqrt((double)n+0.5);
  a = _XS_lehmer_pi(sqrt((double)z)+0.5);
  b = _XS_lehmer_pi(z);
  c = _XS_lehmer_pi(pow((double)n, 1.0/3.0)+0.5);
  TIMING_END_PRINT("stage 1")


  /* We're going to sieve for small primes twice, in the interest of saving
   * memory.  For phi(x,a), get just as many as are needed (a+1) because it
   * is our memory bottleneck.  Then sieve for more afterwards. */
  if (verbose > 0) printf("lehmer %lu stage 2: phi(x,a) (z=%lu a=%lu b=%lu c=%lu)\n", n, z, a, b, c);
  TIMING_START;
  lastprimea = a+1;
  primea = generate_small_primes(&lastprimea);
  if (primea == 0 || lastprimea < a+1) croak("Error generating primes.\n");

  sum = phi(n, a, primea) + ( (b+a-2) * (b-a+1) / 2);

  Safefree(primea);
  TIMING_END_PRINT("phi(x,a)")


  /* Sieve to get small primes.  Get more than the minimum needed (b) to allow
   * fast prime counts.  Using a higher value here will mean more memory but
   * faster operation.  A lower value saves memory at the expense of more
   * segment sieving.*/
  lastprimea = b*16;
  if (verbose > 0) printf("lehmer %lu stage 3: %lu small primes\n", n, lastprimea);
  TIMING_START;
  primea = generate_small_primes(&lastprimea);
  if (primea == 0 || lastprimea < b) croak("Error generating primes.\n");
  lastpc = primea[lastprimea];
  TIMING_END_PRINT("small primes")


  /* Ensure we have the base sieve for big prime_count ( n/primea[i] ). */
  /* This is about 75k for n=10^13, 421k for n=10^15, 2.4M for n=10^17 */
  prime_precalc( sqrt( n / primea[a+1] ) );

  if (verbose > 0) printf("lehmer %lu stage 4: loop %lu to %lu, pc to %lu\n", n, a+1, b, n/primea[a+1]);
  TIMING_START;
  /* Reverse the i loop so w increases.  Count w in segments. */
  lastw = 0;
  lastwpc = 0;
  for (i = b; i >= a+1; i--) {
    UV w = n / primea[i];
    lastwpc = (w <= lastpc) ? bs_prime_count(w, primea, lastprimea)
                            : lastwpc + _XS_prime_count(lastw+1, w);
    lastw = w;
    sum = sum - lastwpc;
    if (i <= c) {
      UV bi = bs_prime_count( (UV) (sqrt(w) + 0.5), primea, lastprimea );
      for (j = i; j <= bi; j++) {
        sum = sum - bs_prime_count(w / primea[j], primea, lastprimea) + j - 1;
      }
    }
  }
  TIMING_END_PRINT("stage 4")
  Safefree(primea);
  return sum;
}
