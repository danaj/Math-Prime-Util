#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lehmer.h"
#include "util.h"
#include "cache.h"
#include "sieve.h"

static int const verbose = 0;

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

/* heap of values and counts, stored by value */
typedef struct {
  UV v;
  IV c;
} vc_t;

typedef struct {
  vc_t* array;
  UV    array_size;
  UV    N;
} heap_t;

static void heap_insert(heap_t* h, vc_t v)
{
  UV k;
  UV n = ++(h->N);
  vc_t* a = h->array;
  if (n >= h->array_size) {
    UV new_size = (h->array_size == 0) ? 20000 : 2 * h->array_size;
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
  vc_t v = heap_remove(h);
  /* get rest of entries with same value */
  while (h->N > 0 && h->array[1].v == v.v) {
    v.c += h->array[1].c;
    (void) heap_remove(h);
  }
  return v;
}
static heap_t heap_create(void)
{
  heap_t h;
  h.array = 0;
  h.array_size = 0;
  h.N = 0;
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
}
static void heap_empty(heap_t* h)
{
  h->N = 0;
}


static UV phi(UV x, UV a)
{
  UV i;
  IV  primea;
  IV  sum = 0;
  if (x < _XS_nth_prime(a+1))  return (x > 0) ? 1 : 0;
  if (a == 1) return ((x+1)/2);
  if (a <= 7) return mapes(x, a);

  primea = _XS_nth_prime(a);

  heap_t h1 = heap_create();
  heap_t h2 = heap_create();
  vc_t v;

  v.v = x;  v.c = 1;
  heap_insert(&h1, v);
  while (a > 5) {
    heap_empty(&h2);
    /* Walk heap */
    while (h1.N > 0) {
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
    primea = _XS_prev_prime(primea);
  }
  heap_destroy(&h2);
  while (h1.N > 0) {
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
  return phi(n, a) + a - 1;
}

UV _XS_lehmer_pi(UV n)
{
  UV z, a, b, c, sum, i, j, lastw, lastwpc;
  UV* primea = 0;   /* small prime cache, first b=pi(z)=pi(sqrt(n)) primes */
  UV* pca = 0;      /* small prime count cache, first z=sqrt(n) numbers */
  if (n < 1000000)
    return _XS_prime_count(2, n);

  if (verbose > 0) printf("lehmer stage 1 - %lu\n", n);
  z = sqrt((double)n+0.5);
  a = _XS_lehmer_pi(sqrt((double)z)+0.5);
  b = _XS_lehmer_pi(z);
  c = _XS_lehmer_pi(pow((double)n, 1.0/3.0)+0.5);
  if (verbose > 0) printf("lehmer stage 2 - %lu (z=%lu a=%lu b=%lu c=%lu)\n", n, z, a, b, c);
  sum = phi(n, a) + ( (b+a-2) * (b-a+1) / 2);

  if (verbose > 0) printf("lehmer stage 3 - n=%lu a=%lu\n", n, a);
  /* Calculate the first b primes and first z primecounts */
  {
    UV count = 0;
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

  /* Ensure we have the base sieve for prime_count ( n/primea[i] ) */
  prime_precalc( sqrt( n / primea[a+1] ) );


  if (verbose > 0) printf("lehmer stage 4 - n=%lu a=%lu\n", n, a);
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
  if (pca != 0)
    free(pca);
  if (primea != 0)
    free(primea);
  return sum;
}
