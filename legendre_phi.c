#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FUNC_isqrt 1
#define FUNC_icbrt 1
#include "util.h"
#include "sieve.h"
#include "lmo.h"
#include "prime_nth_count.h"
#include "prime_count_cache.h"
#include "legendre_phi.h"
#ifdef _OPENMP
  #include <omp.h>
#endif


/*
 * tablephi / tiny_phi
 *    a must be very small (e.g. 6, 7)
 *    direct answer
 *
 * phi_recurse_small
 *    memoryless recursive
 *    calls tablephi, nth_prime, prev_prime, next_prime
 *    good for very small a (e.g. less than 25)
 *
 * phi_recurse
 *    recursive with a small cache
 *    calls tablephi, prime_count_cache, phi_recurse internal
 *    generates primes to nth_prime(a)
 *
 * phi_walk
 *    iterative using list merges
 *    calls tablephi, prime_count_cache, phi_recurse internal
 *    generates primes to nth_prime(a)
 *    complicated, can be much faster than the others, but uses a lot of memory
 *
 * legendre_phi
 *    decides what to do, including handling some special cases
 */

/*============================================================================*/

#define FAST_DIV(x,y) \
  ( ((x) <= 4294967295U) ? (uint32_t)(x)/(uint32_t)(y) : (x)/(y) )

#define PHIC  6U  /* phi(x,a) with a <= PHIC can go to tablephi */

#define PHIS 15U  /* phi(x,a) with a <= PHIS can go to phi_small */
#define PHIS_XMIN (_snth[PHIS+1]-1)   /* nth_prime(PHIS+1)-1 */

#define PHIR 20U  /* phi(x,a) with a <= PHIR is faster with phi_recurse_small */

/*============================================================================*/


/* static const uint8_t _s0[ 1] = {0};
   static const uint8_t _s1[ 2] = {0,1};
   static const uint8_t _s2[ 6] = {0,1,1,1,1,2}; */
static const uint8_t _s3[30] = {0,1,1,1,1,1,1,2,2,2,2,3,3,4,4,4,4,5,5,6,6,6,6,7,7,7,7,7,7,8};
static const uint8_t _s4[210]= {0,1,1,1,1,1,1,1,1,1,1,2,2,3,3,3,3,4,4,5,5,5,5,6,6,6,6,6,6,7,7,8,8,8,8,8,8,9,9,9,9,10,10,11,11,11,11,12,12,12,12,12,12,13,13,13,13,13,13,14,14,15,15,15,15,15,15,16,16,16,16,17,17,18,18,18,18,18,18,19,19,19,19,20,20,20,20,20,20,21,21,21,21,21,21,21,21,22,22,22,22,23,23,24,24,24,24,25,25,26,26,26,26,27,27,27,27,27,27,27,27,28,28,28,28,28,28,29,29,29,29,30,30,30,30,30,30,31,31,32,32,32,32,33,33,33,33,33,33,34,34,35,35,35,35,35,35,36,36,36,36,36,36,37,37,37,37,38,38,39,39,39,39,40,40,40,40,40,40,41,41,42,42,42,42,42,42,43,43,43,43,44,44,45,45,45,45,46,46,47,47,47,47,47,47,47,47,47,47,48};
static UV tablephi(UV x, uint32_t a)
{
  switch (a) {
    case 0: return x;
    case 1: return x-x/2;
    case 2: return x-x/2-x/3+x/6;
    case 3: return (x/    30U) *     8U + _s3[x %     30U];
    case 4: return (x/   210U) *    48U + _s4[x %    210U];
    case 5: {
              UV xp  = x / 11U;
              return ((x /210) * 48 + _s4[x  % 210]) -
                     ((xp/210) * 48 + _s4[xp % 210]);
             }
    case 6:
#if PHIC == 6
    default:
#endif
            {
              UV xp  = x / 11U;
              UV x2  = x / 13U;
              UV x2p = x2 / 11U;
              return ((x  /210) * 48 + _s4[x  % 210]) -
                     ((xp /210) * 48 + _s4[xp % 210]) -
                     ((x2 /210) * 48 + _s4[x2 % 210]) +
                     ((x2p/210) * 48 + _s4[x2p% 210]);
            }
#if PHIC == 7
    case 7:
    default:return tablephi(x,a-1) - tablephi(x/17,a-1);   /* Hacky */
#endif
  }
}
/*============================================================================*/

/* Iterate with simple arrays */

static const unsigned char _snth[25+1] = {0,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};

static UV phi_small(UV x, uint32_t a) {
  UV sum = 0, xpos[1025], xneg[1025];  /* For 32-bit x, 848 is enough */
  uint32_t i, npos, nneg;

  if (a < 4) {
    return (a==0) ? x :
           (a==1) ? x-x/2 :
           (a==2) ? x-x/2-x/3+x/6
                  : (x/30U) * 8U + _s3[x % 30U];
  }
  MPUassert(a <= PHIS, "phi_small: a too large");
  if (x < _snth[a+1]) return (x>0);

  for (npos = nneg = 0, xpos[npos++] = x;  a > 4U;  a--) {
    uint32_t oneg = nneg,  opos = npos;
    for (i = 0; i < opos; i++)
      if (xpos[i] >= _snth[a])
        xneg[nneg++] = xpos[i]/_snth[a];
    for (i = 0; i < oneg; i++)
      if (xneg[i] >= _snth[a])
        xpos[npos++] = xneg[i]/_snth[a];
  }
  for (i = 0; i < npos; i++)
    sum += (xpos[i]/210U)*48U + _s4[xpos[i] % 210U];
  for (i = 0; i < nneg; i++)
    sum -= (xneg[i]/210U)*48U + _s4[xneg[i] % 210U];
  return sum;
}

/*============================================================================*/

/* Recurse until a <= PHIS */

static UV phi_recurse_small(UV x, UV a, UV npa) {
  UV sum, i, xp, p, lp;

  if (x < 1 || a >= x) return (x > 0);
  if (a <= PHIS || x <= PHIS_XMIN)  return phi_small(x, a);

  sum = phi_small(x, PHIS);
  p = _snth[PHIS];
  for (i = PHIS+1; i <= a; i++) {
    lp = p;
    p = next_prime(p);
    xp = FAST_DIV(x,p);
    if (xp < p) {
      while (x < npa) {
        a--;
        npa = prev_prime(npa);
      }
      return (sum - a + i - 1);
    }
    sum -= phi_recurse_small(xp, i-1, lp);
  }
  return sum;
}

/*============================================================================*/
/*============================================================================*/

/* Cache for phi(x,a) */

#define PHICACHEA 512
typedef struct
{
  uint32_t  siz[PHICACHEA];   /* how many entries we have allocated */
  uint16_t *val[PHICACHEA];
  uint32_t  xlim;
} phi_cache_t;

static phi_cache_t* phi_cache_create(uint32_t xlim) {
  phi_cache_t *cache;
  int a;
  New(0, cache, 1, phi_cache_t);
  for (a = 0; a < PHICACHEA; a++) {
    cache->val[a] = 0;
    cache->siz[a] = 0;
  }
  cache->xlim = (xlim < 0xFFFFFFFFU) ? xlim : xlim-1;  /* Reserve 0xFFFFFFFF */
  return cache;
}

static void phi_cache_destroy(phi_cache_t* cache) {
  int a;
  for (a = 0; a < PHICACHEA; a++) {
    if (cache->val[a] != 0)
      Safefree(cache->val[a]);
  }
  Safefree(cache);
}

static void phi_cache_insert(uint32_t x, uint32_t a, IV sum, phi_cache_t* cache) {
  uint32_t i, newsize;
  if (sum < 0) sum = -sum;
  if (sum > 65535) return;  /* If sum is too large for the cache, ignore it. */
  if (x >= cache->siz[a]) {
    newsize = (x >= 0xFFFFFFFFUL-32)  ?  0xFFFFFFFFUL-1  :  x+32;
    if (cache->val[a] == 0) {
      Newz(0, cache->val[a], newsize, uint16_t);
    } else {
      Renew(cache->val[a], newsize, uint16_t);
      for (i = cache->siz[a]; i < newsize; i++)  /* Zero the new entries */
        cache->val[a][i] = 0;
    }
    cache->siz[a] = newsize;
  }
  cache->val[a][x] = (uint16_t) sum;
}

#define PHI_CACHE_HAS_VAL(x, a, pcache) \
  ((a) < PHICACHEA && (x) < (pcache)->siz[a] && (pcache)->val[a][x] != 0)

#define PHI_CACHE_INSERT(x, a, val, pcache) \
  if ((a) < PHICACHEA && (x) <= (pcache)->xlim) \
    { phi_cache_insert((x), (a), (val), pcache); }

#define PHI_CACHE_GET_VAL(x, a, pcache) \
  (pcache)->val[a][x]

/* End of Phi cache definitions */

/* Struct of everything needed for recursive phi call */

typedef struct {
  const uint32_t* primes;
  uint32_t lastidx;
  void* cachepc;
  phi_cache_t* cachephi;
} phidata_t;

static phidata_t* phidata_create(const uint32_t* primes, uint32_t lastidx, UV x, UV a)
{
  phidata_t *d;
  uint32_t xlim = isqrt(x);  /* Not ideal */

  New(0, d, 1, phidata_t);
  d->primes = primes;
  d->lastidx = lastidx;
  d->cachepc = prime_count_cache_create_with_primes(primes, lastidx);
  d->cachephi = phi_cache_create(xlim);
  return d;
}
static void phidata_destroy(phidata_t *d)
{
  phi_cache_destroy(d->cachephi);
  prime_count_cache_destroy(d->cachepc);
  /* They own the primes */
  Safefree(d);
}

#define PHI_IS_X_SMALL(x, a) \
  ( ((x) <= primes[d->lastidx]) && ((x) < (UV)primes[a+1] * primes[a+1]) )
#define PHI_PRIMECOUNT(x) \
  prime_count_cache_lookup(d->cachepc, (x))

/* The recursive cached phi routine, given the struct with primes and cache */

static IV _phi3(UV x, UV a, int sign, phidata_t *d)
{
  const uint32_t* const primes = d->primes;
  phi_cache_t* pcache = d->cachephi;

  if (x < primes[a+1])
    return sign;
  else if (a <= PHIC)
    return sign * tablephi(x,a);
  else if (PHI_IS_X_SMALL(x,a))
    return sign * (PHI_PRIMECOUNT(x) - a + 1);
  else if (PHI_CACHE_HAS_VAL(x, a, pcache))
    return sign * PHI_CACHE_GET_VAL(x, a, pcache);
  else {
    UV ai, iters = ((UV)a*a > x)  ?  PHI_PRIMECOUNT(isqrt(x))  :  a;
    UV c = (iters > PHIC) ? PHIC : iters;
    IV sum = sign * (iters - a + tablephi(x,c));

    for (ai = c+1; ai <= iters; ai++)
      sum += _phi3(FAST_DIV(x,primes[ai]), ai-1, -sign, d);

    PHI_CACHE_INSERT(x, a, sum, pcache);
    return sum;
  }
}

static UV phi_recurse(UV x, UV a, UV primes_to_n)
{
  uint32_t* primes;
  uint32_t lastidx;
  UV sum = 1;

  if (x < 1 || a >= x) return (x > 0);
  if (a <= PHIS || x <= PHIS_XMIN)  return phi_small(x, a);
  if (a > 203280221) croak("64-bit phi out of range");

  lastidx = range_prime_sieve_32(&primes, primes_to_n, 1);

  if (primes[a] < x) {
    phidata_t *d = phidata_create(primes, lastidx, x, a);
    /* Ensure testing with legendre_phi(1e13, 203280221) +/- 2 */
    /* sum = (UV) _phi3(x, a, 1, d); */
    sum = (UV) _phi3(x, a-1, 1, d)  -  (UV) _phi3(x/primes[a], a-1, 1, d);
    phidata_destroy(d);
  }

  Safefree(primes);
  return sum;
}

/*============================================================================*/
/*============================================================================*/

static int const verbose = 0;
#define MAX_PHI_MEM  (896*1024*1024)
#define NTHRESH (MAX_PHI_MEM/16)

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
    l->size = new_size;
  }
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
    UV bval = ba[bi].v;
    /* Skip forward in A until empty or aa[ai].v <= ba[bi].v */
    while (ai+8 < an && aa[ai+8].v > bval)  ai += 8;
    while (ai   < an && aa[ai  ].v > bval)  ai++;
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
  for (k = kn-1; k >= 0 && bi >= 0; k--) {
    UV bval = ba[bi].v;
    long startai = ai;
    while (ai >= 15 && aa[ai-15].v < bval) ai -= 16;
    while (ai >=  3 && aa[ai- 3].v < bval) ai -= 4;
    while (ai >=  0 && aa[ai   ].v < bval) ai--;
    if (startai > ai) {
      k = k - (startai - ai) + 1;
      memmove(aa+k, aa+ai+1, (startai-ai) * sizeof(vc_t));
    } else {
      if (ai >= 0 && aa[ai].v == bval) croak("deduplication error");
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

/* phi(x,a) non-recursive, using list merging.   Memory intensive. */

static UV phi_walk(UV x, UV a, UV primes_to_n)
{
  UV i, val, sval, lastidx, lastprime;
  UV sum = 0;
  IV count;
  uint32_t* primes;
  vcarray_t a1, a2;
  vc_t* arr;
  phidata_t *d;

  if (x < 1 || a >= x) return (x > 0);
  if (x <= PHIC || a <= PHIC)  return tablephi(x, (a > PHIC) ? PHIC : a);
  if (a > 203280221) croak("64-bit phi out of range");

  lastidx = range_prime_sieve_32(&primes, primes_to_n, 1);
  lastprime = primes[lastidx];
  if (x < lastprime)  { Safefree(primes); return 1; }

  d = phidata_create(primes, lastidx, x, a);

  a1 = vcarray_create();
  a2 = vcarray_create();
  vcarray_insert(&a1, x, 1);

  while (a > PHIC) {
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
          if (PHI_IS_X_SMALL(sval_last, a-1))
            sum += sval_count*(PHI_PRIMECOUNT(sval_last)-a+2);
          else
            vcarray_insert(&a2, sval_last, sval_count);
        }
        sval_last = sval;
        sval_count = 0;
      }
      sval_count -= count;           /* Accumulate count for this sval */
    }
    if (sval_last != 0) {            /* Insert the last sval */
      if (PHI_IS_X_SMALL(sval_last, a-1))
        sum += sval_count*(PHI_PRIMECOUNT(sval_last)-a+2);
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
          sum += count * _phi3( arr[j].v, a-1, 1, d );
          arr[j].c = 0;
        }
      }
    }
    vcarray_remove_zeros(&a1);
    a--;
  }
  phidata_destroy(d);
  Safefree(primes);
  vcarray_destroy(&a2);
  arr = a1.a;
#ifdef _OPENMP
  #pragma omp parallel for reduction(+: sum) schedule(dynamic, 16)
#endif
  for (i = 0; i < a1.n; i++)
    sum += arr[i].c * tablephi( arr[i].v, PHIC );
  vcarray_destroy(&a1);
  return (UV) sum;
}

/*============================================================================*/
/*============================================================================*/

uint32_t tiny_phi_max_a(void) { return PHIC; }

UV tiny_phi(UV n, uint32_t a) {
  return (a <= PHIC) ? tablephi(n, a)
       : (a <= PHIS) ? phi_small(n, a)
                     : phi_recurse_small(n, a, nth_prime(a));
}

uint32_t small_phi_max_a(void) { return PHIS; }

UV small_phi(UV n, uint32_t a) {
  UV npa;
  if (a <= PHIS) return phi_small(n, a);
  npa = nth_prime(a);
  if (n <= npa) return 1;
  if (npa < isqrt(n)) npa = isqrt(n);
  return phi_recurse(n, a, npa);
}

/*============================================================================*/
/*============================================================================*/

void* prepare_cached_legendre_phi(UV x, UV a)
{
  uint32_t npa, lastidx, *primes;
  UV sum;

  if (a > 203280221)  a = 203280221;
  npa = nth_prime_upper(a);
  if (npa < isqrt(x)) npa = isqrt(x);
  lastidx = range_prime_sieve_32(&primes, npa, 1);
  return (void*) phidata_create(primes, lastidx, x, a);
}
UV cached_legendre_phi(void* cache, UV x, UV a)
{
  phidata_t *d = (phidata_t*) cache;

  if (x < 1 || a >= x) return (x > 0);
  if (x <= PHIC || a <= PHIC)  return tablephi(x, (a > PHIC) ? PHIC : a);
  if (a > (x >> 1))  return 1;

  /* Make the function work even if x,a outside of cached conditions */
  if (a > 203280221) {  /* prime_count(2**32) */
    UV pc = LMO_prime_count(x);
    return (a >= pc)  ?  1  :  pc - a + 1;
  }
  if (a > d->lastidx)
    return legendre_phi(x, a);

  return (UV) _phi3(x, a-1, 1, d)  -  (UV) _phi3(x/d->primes[a], a-1, 1, d);
}

void destroy_cached_legendre_phi(void* cache)
{
  phidata_t *d = (phidata_t*) cache;
  Safefree(d->primes);
  phidata_destroy(d);
}

/* static UV phi_stupid(UV x, UV a) {
  if (a <= PHIC) return tablephi(x,a);
  return phi_stupid(x, a-1) - phi_stupid(x/nth_prime(a), a-1);
} */

/*============================================================================*/
/*============================================================================*/

UV legendre_phi(UV x, UV a)
{
  UV npa;
  /* If 'x' is very small, give a quick answer with any 'a' */
  if (x < 1 || a >= x) return (x > 0);
  if (x <= PHIC || a <= PHIC)  return tablephi(x, (a > PHIC) ? PHIC : a);

  /* Two shortcuts for large values, from R. Andrew Ohana */
  if (a > (x >> 1))  return 1;

  /* If a > prime_count(2^32), then we need not be concerned with composite
   * x values with all factors > 2^32, as x is limited to 64-bit. */
  if (a > 203280221) {  /* prime_count(2**32) */
    UV pc = LMO_prime_count(x);
    return (a >= pc)  ?  1  :  pc - a + 1;
  }
  /* After the two checks above,  7 <= a <= MIN(203280221, 2*x) */

  if (a <= PHIS)  return phi_small(x, a);

  if (a > PHIR  &&  a >= prime_count_upper(x)) return 1;
  npa = nth_prime(a);
  if (x <= npa) return 1;  /* Technically x < next_prime(npa) */
  if (a <= PHIR) return phi_recurse_small(x, a, npa); /* npa must = P[a] */

  /* For large x, it is faster to compute more primes. */
  if (npa < isqrt(x)) npa = isqrt(x);

  /* The best crossover between recurse and walk is complicated */
  if (x <= 1e10)
    return phi_recurse(x, a, npa);

  /* TODO: More tuning of these, or just improve them. */

  if (x <= 1e11) return (a <   5000) ? phi_walk(x,a,npa) : phi_recurse(x,a,npa);
  if (x <= 1e12) return (a <  15000) ? phi_walk(x,a,npa) : phi_recurse(x,a,npa);
  if (x <= 1e13) return (a <  31000) ? phi_walk(x,a,npa) : phi_recurse(x,a,npa);
  if (x <= 1e14) return (a <  70000) ? phi_walk(x,a,npa) : phi_recurse(x,a,npa);
  if (x <= 1e15) return (a < 120000) ? phi_walk(x,a,npa) : phi_recurse(x,a,npa);

  return (a < 200000) ? phi_walk(x,a,npa) : phi_recurse(x,a,npa);
}
