#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FUNC_isqrt 1
#include "util.h"
#include "sieve.h"
#include "lmo.h"
#include "prime_counts.h"
#include "prime_count_cache.h"
#include "legendre_phi.h"


/*
 * tablephi / tiny_phi
 *    a must be very small (e.g. 6, 7)
 *    direct answer
 *
 * phi_small
 *    a must be very small (e.g. 15)
 *    calls tablephi
 *    simple iteration using fixed size lists
 *
 * phi_recurse_small
 *    memoryless recursive
 *    calls phi_small, nth_prime (if a > 25), prev_prime, next_prime
 *    good for very small a (e.g. less than 25)
 *
 * phi_recurse
 *    recursive with a cache
 *    calls tablephi, prime_count_cache, phi_recurse internal
 *    generates primes to max(nth_prime(a),isqrt(x))
 *
 * phi_walk
 *    iterative using list merges
 *    calls tablephi, prime_count_cache, phi_recurse internal
 *    generates primes to max(nth_prime(a),isqrt(x))
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

/* For x >= 1 and a >= 4,  phi(x,a) = phi(x-_pred7[x%210],a)
 * This allows us to collapse multiple x values, useful for caching. */
static const unsigned char _pred7[210] = {1,0,1,2,3,4,5,6,7,8,9,0,1,0,1,2,3,0,1,0,1,2,3,0,1,2,3,4,5,0,1,0,1,2,3,4,5,0,1,2,3,0,1,0,1,2,3,0,1,2,3,4,5,0,1,2,3,4,5,0,1,0,1,2,3,4,5,0,1,2,3,0,1,0,1,2,3,4,5,0,1,2,3,0,1,2,3,4,5,0,1,2,3,4,5,6,7,0,1,2,3,0,1,0,1,2,3,0,1,0,1,2,3,0,1,2,3,4,5,6,7,0,1,2,3,4,5,0,1,2,3,0,1,2,3,4,5,0,1,0,1,2,3,0,1,2,3,4,5,0,1,0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,0,1,0,1,2,3,0,1,2,3,4,5,0,1,0,1,2,3,4,5,0,1,2,3,0,1,0,1,2,3,0,1,0,1,2,3,4,5,6,7,8,9,0};

/* Maps x to value <= x not divisible by first 4 primes */
/* mpu 'say join(",",map { legendre_phi($_,4)-1 } 0..209);' */
static const int8_t _coprime_idx210[210]={-1,0,0,0,0,0,0,0,0,0,0,1,1,2,2,2,2,3,3,4,4,4,4,5,5,5,5,5,5,6,6,7,7,7,7,7,7,8,8,8,8,9,9,10,10,10,10,11,11,11,11,11,11,12,12,12,12,12,12,13,13,14,14,14,14,14,14,15,15,15,15,16,16,17,17,17,17,17,17,18,18,18,18,19,19,19,19,19,19,20,20,20,20,20,20,20,20,21,21,21,21,22,22,23,23,23,23,24,24,25,25,25,25,26,26,26,26,26,26,26,26,27,27,27,27,27,27,28,28,28,28,29,29,29,29,29,29,30,30,31,31,31,31,32,32,32,32,32,32,33,33,34,34,34,34,34,34,35,35,35,35,35,35,36,36,36,36,37,37,38,38,38,38,39,39,39,39,39,39,40,40,41,41,41,41,41,41,42,42,42,42,43,43,44,44,44,44,45,45,46,46,46,46,46,46,46,46,46,46,47};
static UV _toindex210(UV x) {
  UV q = x / 210, r = x % 210;
  return 48 * q + _coprime_idx210[r];
}

/* Small table of nth primes */
static const unsigned char _snth[25+1] = {0,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};

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

/* Iterate with simple arrays, no merging or cleverness. */

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

static UV phi_recurse_small(UV x, UV a) {
  UV sum, i, xp, p, npa;

  if (x < 1 || a >= x) return (x > 0);
  if (a <= PHIS || x <= PHIS_XMIN)  return phi_small(x, a);

  npa = (a <= 25) ? _snth[a] : nth_prime(a);
  sum = phi_small(x, PHIS);
  p = _snth[PHIS];
  for (i = PHIS+1; i <= a; i++) {
    p = next_prime(p);
    xp = FAST_DIV(x,p);
    if (xp < p) {
      while (x < npa) {
        a--;
        npa = prev_prime(npa);
      }
      return (sum - a + i - 1);
    }
    sum -= phi_recurse_small(xp, i-1);
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
  uint32_t xlim = (UV) pow(x, 1.0/2.70);
  if (xlim < 256) xlim = 256;

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
  UV mapx;

  if (x < primes[a+1])
    return sign;
  else if (a <= PHIC)
    return sign * tablephi(x,a);
  else if (PHI_IS_X_SMALL(x,a))
    return sign * (PHI_PRIMECOUNT(x) - a + 1);

  /* Choose a mapping:   x,  (x+1)>>1,  _toindex30(x),  _toindex210(x) */
  mapx = (a < PHICACHEA)  ?  _toindex210(x)  :  0;

  if (a < PHICACHEA && mapx < pcache->siz[a]) {
    IV v = pcache->val[a][mapx];
    if (v != 0)
      return sign * v;
  }
  {
    UV xp, i, iters = ((UV)a*a > x)  ?  PHI_PRIMECOUNT(isqrt(x))  :  a;
    UV c = (iters > PHIC) ? PHIC : iters;
    IV sum = sign * (iters - a + tablephi(x,c));

    /* for (i=c; i<iters; i++)  sum += _phi3(x/primes[i+1], i, -sign, d); */

    if (c < iters)
      sum += -sign * tablephi(FAST_DIV(x,primes[c+1]), c);
    for (i = c+1; i < iters; i++) {
      xp = FAST_DIV(x,primes[i+1]);
      if (PHI_IS_X_SMALL(xp,i))
        break;
      sum += _phi3(xp, i, -sign, d);
    }
    for (; i < iters; i++) {
      xp = FAST_DIV(x,primes[i+1]);
      if (xp < primes[i+1])
        break;
      sum += -sign * (PHI_PRIMECOUNT(xp) - i + 1);
    }
    if (i < iters)
      sum += -sign * (iters - i);

    if (a < PHICACHEA && mapx <= pcache->xlim)
      phi_cache_insert(mapx, a, sum, pcache);
    return sum;
  }
}

static UV phi_recurse(UV x, UV a)
{
  uint32_t* primes;
  uint32_t lastidx;
  UV primes_to_n, sum = 1;

  if (x < 1 || a >= x) return (x > 0);
  if (a <= PHIS || x <= PHIS_XMIN)  return phi_small(x, a);
  if (a > 203280221) croak("64-bit phi out of range");

  primes_to_n = nth_prime_upper(a);
  if (isqrt(x) > primes_to_n) primes_to_n = isqrt(x);
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

/* Insert a value/count pair.  Must be done in decreasing size order. */
static void vcarray_insert(vcarray_t* l, UV val, IV count)
{
  UV n = l->n;
  vc_t* arr = l->a;

  if (n > 0 && arr[n-1].v <= val) {
    if (arr[n-1].v == val) {
      arr[n-1].c += count;
      return;
    }
    croak("Previous value was %lu, inserting %lu out of order\n", arr[n-1].v, val);
  }
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
    arr = l->a;
  }
  arr[n].v = val;
  arr[n].c = count;
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

static UV phi_walk(UV x, UV a)
{
  UV i, sval, lastidx, lastprime, primes_to_n;
  UV sum = 0;
  uint32_t* primes;
  vcarray_t a1, a2;
  vc_t* arr;
  phidata_t *d;

  if (x < 1 || a >= x) return (x > 0);
  if (x <= PHIC || a <= PHIC)  return tablephi(x, (a > PHIC) ? PHIC : a);
  if (a > 203280221) croak("64-bit phi out of range");

  primes_to_n = nth_prime_upper(a);
  if (isqrt(x) > primes_to_n) primes_to_n = isqrt(x);

  lastidx = range_prime_sieve_32(&primes, primes_to_n, 1);
  lastprime = primes[lastidx];
  if (x < lastprime)  { Safefree(primes); return 1; }

  d = phidata_create(primes, lastidx, x, a);

  a1 = vcarray_create();
  a2 = vcarray_create();
  vcarray_insert(&a1, x, 1);

  while (a > PHIC) {
    UV primea = primes[a];
    arr = a1.a;

    for (i = 0; i < a1.n; i++) {
      sval = FAST_DIV(arr[i].v, primea);
      sval -= _pred7[sval % 210];     /* Reduce to lower value if possible */
      if (sval < primea || PHI_IS_X_SMALL(sval, a-1))
        break;
      vcarray_insert(&a2, sval, -arr[i].c);
    }
    for ( ; i < a1.n; i++) {
      sval = FAST_DIV(arr[i].v, primea);
      if (sval < primea)
        break;
      sum -= arr[i].c * (PHI_PRIMECOUNT(sval)-a+2);
    }
    for ( ; i < a1.n; i++)
      sum -= arr[i].c;

    /* Merge a1 and a2 into a1.  a2 will be emptied. */
    vcarray_merge(&a1, &a2);

    /* If we've grown too large, use recursive phi to clip. */
    if ( a1.n > NTHRESH ) {
      arr = a1.a;
      if (verbose > 0) printf("clipping small values at a=%lu a1.n=%lu \n", a, a1.n);
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
                     : phi_recurse_small(n, a);
}

uint32_t small_phi_max_a(void) { return PHIS; }

UV small_phi(UV n, uint32_t a) {
  return (a <= PHIS) ? phi_small(n, a) : phi_recurse(n, a);
}

/*============================================================================*/
/*============================================================================*/

void* prepare_cached_legendre_phi(UV x, UV a)
{
  uint32_t npa, lastidx, *primes;

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
  UV sqrtx = isqrt(x);

  /* If 'x' is very small, give a quick answer with any 'a' */
  if (x < 1 || a >= x) return (x > 0);
  if (x <= PHIC || a <= PHIC)  return tablephi(x, (a > PHIC) ? PHIC : a);

  /* Very fast shortcuts for large values */
  if (a > (x >> 1))
    return 1;
  if (a >= sqrtx || a > 203280221) {   /* 203280221 = prime_count(2^32) */
    UV pc = LMO_prime_count(x);
    return (a >= pc)  ?  1  :  pc - a + 1;
  }
  /* After this:  7 <= a <= MIN(203280221, sqrtx) */

  /* For very small a, calculate now. */
  if (a <= PHIS)  return phi_small(x, a);
  if (a <= PHIR)  return phi_recurse_small(x, a);

  /* Better shortcuts, slightly more time */
  if (prime_count_upper(x) <= a)
    return 1;
  /* Use 'a' instead of 'a+1' to ensure Legendre Pi doesn't call here */
  if (prime_count_upper(sqrtx) < a) {
    UV pc = LMO_prime_count(x);
    return (a >= pc)  ?  1  :  pc - a + 1;
  }
  /* Because we used the fast bounds, there are still a few easy cases. */

  /* The best crossover between recurse and walk is complicated */
  /* TODO: More tuning of the crossovers, or just improve the algorithms. */

  if (x < 1e10)
    return phi_recurse(x, a);

  if ( (x >= 1e10 && x < 1e11 && a <   2000) ||
       (x >= 1e11 && x < 1e12 && a <   4000) ||
       (x >= 1e12 && x < 1e13 && a <  10000) ||
       (x >= 1e13 && x < 1e14 && a <  24000) ||
       (x >= 1e14 && x < 1e15 && a <  80000) ||
       (x >  1e15             && a < 150000) )
    return phi_walk(x, a);

  return phi_recurse(x, a);
}




/*============================================================================*/

#if 0
  // TODO: setup with initial function.  optimize.  export.
IV phi_sum(UV x, UV a, int sign) {
  IV sum = 0;
  //if (x < 1) return 0;
  for (; a > 0; a--) {
    UV p = nth_prime(a);
    if (x <= p) {
      return sum + (long)sign;
    }
    sum += p * phi_sum(x / p, a-1, -sign);
  }
  if (sign > 0) sum += (x*(x+1))/2; else sum -= (x*(x+1))/2;
  return sum;
}
#endif
