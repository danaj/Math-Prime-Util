/******************************************************************************/
/*                              OMEGA PRIMES                                  */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ptypes.h"
#include "constants.h"
#include "cache.h"
#include "sieve.h"
#include "util.h"
#include "sort.h"
#include "prime_counts.h"
#include "prime_powers.h"
#include "factor.h"
#include "inverse_interpolate.h"
#include "omega_primes.h"

bool is_omega_prime(uint32_t k, UV n) {
  if (n <= 1) return (n == 1 && k == 0);

  if (k > 0 && !(n& 1)) { k--; do { n >>= 1; } while (!(n& 1)); }
  if (k > 0 && !(n% 3)) { k--; do { n /=  3; } while (!(n% 3)); }
  if (k > 0 && !(n% 5)) { k--; do { n /=  5; } while (!(n% 5)); }
  if (k > 0 && !(n% 7)) { k--; do { n /=  7; } while (!(n% 7)); }
  if (k > 0 && !(n%11)) { k--; do { n /= 11; } while (!(n%11)); }

  if (n == 1) return (k == 0);
  if (k == 0) return (n == 1);
  if (k == 1) return is_prime_power(n);
  if (n < ipowsafe(13,k)) return 0;

  return ((UV)prime_omega(n) == k);
}

/* See https://arxiv.org/pdf/2006.16491.pdf page 12 for a brief note */

/* For the interpolation */
static UV opce(UV mid, UV k) { return omega_prime_count(k, mid); }


/*********************************  Construction  *****************************/

static void _omega_prime_push(UV** kop, UV* skop, UV* nkop, UV v)
{
  UV *l = *kop, lsize = *skop, n = *nkop;

  if (n >= lsize) {
    UV nsize = lsize + lsize/5 + 1;
    if (nsize <= lsize || nsize > (UV)(MAX_SIZET / sizeof(UV)))
      croak("omega prime list is too large");
    Renew(l, nsize, UV);
    *kop = l;
    *skop = nsize;
  }
  l[n++] = v;
  *nkop = n;
}

static void _omega_prime_gen_one(UV** kop, UV* skop, UV* nkop,
                                 uint32_t e, UV lo, UV hi, UV m, UV p)
{
  UV ppow, v;

  if ((m % p) == 0) return;
  ppow = (e == 1) ? p : ipowsafe(p, e);
  if (ppow > hi/m) return;
  v = m * ppow;
  if (v >= lo)
    _omega_prime_push(kop, skop, nkop, v);
}

static void _omega_prime_gen_pow(UV** kop, UV* skop, UV* nkop, uint32_t e,
                                 UV lo, UV hi, UV m, UV pmin, UV pmax)
{
  const UV cache_limit = UVCONST(100000000);

  if (pmax < pmin) return;

  if (pmax <= cache_limit || pmax <= get_prime_cache(0, 0)) {
    START_DO_FOR_EACH_PRIME(pmin, pmax) {
      _omega_prime_gen_one(kop, skop, nkop, e, lo, hi, m, p);
    } END_DO_FOR_EACH_PRIME
  } else {
    if (pmin <= 2 && pmax >= 2)
      _omega_prime_gen_one(kop, skop, nkop, e, lo, hi, m, 2);
    if (pmin <= 3 && pmax >= 3)
      _omega_prime_gen_one(kop, skop, nkop, e, lo, hi, m, 3);
    if (pmin <= 5 && pmax >= 5)
      _omega_prime_gen_one(kop, skop, nkop, e, lo, hi, m, 5);
    if (pmin < 7) pmin = 7;

    if (pmin <= pmax) {
      unsigned char *segment;
      UV seg_base, seg_low, seg_high;
      void *ctx = start_segment_primes(pmin, pmax, &segment);

      while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
        START_DO_FOR_EACH_SIEVE_PRIME(segment, seg_base, seg_low, seg_high) {
          _omega_prime_gen_one(kop, skop, nkop, e, lo, hi, m, p);
        } END_DO_FOR_EACH_SIEVE_PRIME;
      }
      end_segment_primes(ctx);
    }
  }
}

static void _omega_prime_gen_last(UV** kop, UV* skop, UV* nkop,
                                  UV lo, UV hi, UV m, UV pstart)
{
  UV minpow = lo/m + ((lo % m) != 0);
  UV maxpow = hi/m;
  uint32_t e;

  for (e = 1; ; e++) {
    UV pmin, pmax = rootint(maxpow, e);
    if (pmax <= pstart) break;
    pmin = crootint(minpow, e);
    if (pmin < pstart) pmin = pstart;
    _omega_prime_gen_pow(kop, skop, nkop, e, lo, hi, m, pmin, pmax);
  }
}

static void _omega_prime_gen_rec(UV** kop, UV* skop, UV* nkop, uint32_t k, UV lo, UV hi, UV m, UV pstart) {
  UV v;

  if (k > 1) {
    SIMPLE_FOR_EACH_PRIME(pstart, rootint(hi/m, k)) {
      if ((m % p) == 0) continue;
      if (UV_MAX/m < p) break;
      for (v = m*p; UV_MAX/v >= p && v*p <= hi; v *= p)
        _omega_prime_gen_rec(kop, skop, nkop, k-1, lo, hi, v, p);
    } END_SIMPLE_FOR_EACH_PRIME
    return;
  }

  _omega_prime_gen_last(kop, skop, nkop, lo, hi, m, pstart);
}

static UV rec_omega_primes(UV** ret, uint32_t k, UV lo, UV hi) {
  UV nkop, skop;

  if (hi < lo) croak("range_omega_prime_sieve error hi %"UVuf" < lo %"UVuf"\n",hi,lo);

  nkop = 0;
  skop = 256;
  New(0, *ret, skop, UV);
  _omega_prime_gen_rec(ret, &skop, &nkop, k, lo, hi, 1, 2);
  sort_uv_array(*ret, nkop);
  return nkop;
}

static UV sieve_omega_primes(UV** ret, uint32_t k, UV lo, UV hi)
{
  const UV chunk_size = UVCONST(5000000);
  UV *l, lsize = 1000, n = 0;

  New(0, l, lsize, UV);
  while (lo <= hi) {
    UV i, chi, range;
    unsigned char *nf;

    chi = (hi-lo >= chunk_size) ? lo+chunk_size-1 : hi;
    range = chi-lo+1;
    nf = range_nfactor_sieve(lo, chi, 0);
    for (i = 0; i < range; i++)
      if (nf[i] == k)
        _omega_prime_push(&l, &lsize, &n, lo+i);
    Safefree(nf);
    if (chi == hi) break;
    lo = chi+1;
  }
  *ret = l;
  return n;
}


UV range_omega_prime_sieve(UV** ret, uint32_t k, UV lo, UV hi) {
  UV min, *l;

  if (hi < lo) croak("range_omega_prime_sieve error hi %"UVuf" < lo %"UVuf"\n",hi,lo);

  if (ret == 0)
    return omega_prime_count(k,hi) - (lo <= 1 ? 0 : omega_prime_count(k,lo-1));

  *ret = 0;

  if (k == 0) {
    if (lo <= 1 && hi >= 1) {
      New(0, l, 1, UV);
      l[0] = 1;
      *ret = l;
      return 1;
    }
    return 0;
  }

  min = pn_primorial(k);
  if (min == 0 || min > hi) return 0;
  if (lo < min) lo = min;

  if (k == 1) return prime_power_sieve(ret, lo, hi);

  /* Recursion efficiently constructs a prefix.  It does substantial work below
   * lo, so use bounded factor-sieve chunks for arbitrary ranges. */
  if (lo == min &&
      (((hi-lo) > 100000000UL) || (k >= 10 && (hi-lo) > 5000000UL)))
    return rec_omega_primes(ret, k, lo, hi);

  return sieve_omega_primes(ret, k, lo, hi);
}


/*********************************  Counting  *********************************/

UV max_omega_prime_count(uint32_t k) {
#if BITS_PER_WORD == 32
  static const UV max[10] = {1,203287168,838888926,1389246717,1178725572,540561553,129357524,14327954,567659,4221};
  if (k >= 10) return 0;
#else
  static const UV max[16] = {1, UVCONST(425656284140516112), /* prime powers */
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* TODO: find these */
    UVCONST(5512134903353),UVCONST(87133769732),UVCONST(446745559),299178
  };
  if (k >= 16) return 0;
#endif
  if (k >= BITS_PER_WORD) return 0;
  if (max[k] == 0) return UV_MAX;
  return max[k];
}

UV max_nth_omega_prime(uint32_t k) {
#if BITS_PER_WORD == 32
  static const UV offset[10] = {0,4,1,8,5,0,34,3,1305,46665};
  if (k >= 10) return 0;
#else
  static const UV offset[16] = {0,58,7,2,3,5,25,0,48,255,1155,46017,15,
  UVCONST(125585475),UVCONST(522131625),UVCONST(338362334325)};
  if (k >= 16) return 0;
#endif
  if (k >= BITS_PER_WORD) return 0;
  return UV_MAX - offset[k];
}


#define RECURSIVE_OMEGA_COUNT(k,n,pr,npr) \
  _omega_prime_count_rec2(k, n, 1, 2, rootint(n,k), 1, pr, npr)

/*  Initial call:   m = 1, p = 2, s = sqrtn(n), j = 1  */
static UV _omega_prime_count_rec2(uint32_t k, UV n, UV m, UV p, UV s, UV j, uint32_t* pr, UV numprimes) {
  UV t, r, count = 0;

  if (k == 2) {
    UV r2, w;
    uint32_t e;
    for (;  p <= s;  j++, p = r) {
      r = (j < numprimes)  ?  pr[j]  :  next_prime(p);
      for (t = m*p, w = n/t;  t <= n && w >= r;  t *= p, w = n/t) {
        count += prime_count(w) - j;
        for (e = 2; ; e++) {
          r2 = rootint(w, e);
          if (r2 < r) break;
          count += prime_count(r2) - j;
        }
        if (t > n/p) break;
      }
    }
    return count;
  }

  for (;  p <= s;  j++, p = r) {
    r = (j < numprimes)  ?  pr[j]  :  next_prime(p);
    for (t = m*p; t <= n; t *= p) {
      UV S = rootint(n/t, k-1);
      if (r > S) break;
      count += _omega_prime_count_rec2(k-1, n, t, r, S, j+1, pr, numprimes);
      if (t > n/p) break;
    }
  }
  return count;
}

UV omega_prime_count(uint32_t k, UV n)
{
  uint32_t* pr;
  UV maxpr, npr, sum, lo;

  if (k == 0) return (n >= 1);
  if (k == 1) return prime_power_count(n);

  /* The first k-omega-prime is primorial(p_k) (ignoring zero for k=1) */
  lo = pn_primorial(k);
  if (lo == 0 || n < lo) return 0;

  maxpr = rootint(n, (k > 10)  ?  4  :  (k > 6)  ?  3  :  2);
  npr = range_prime_sieve_32(&pr, maxpr, 0);  /* p[0]=2, p[1]=3,... */
  sum = RECURSIVE_OMEGA_COUNT(k, n, pr, npr);
  Safefree(pr);
  return sum;
}

/* An upper bound for the omega prime count, when n >= 10^12 is shown in
 * Bayless,Kinlaw,Klyve 2019, page 4
 * https://www.researchgate.net/profile/Paul-Kinlaw/publication/329788487_Sums_over_primitive_sets_with_a_fixed_number_of_prime_factors/links/5c44103d92851c22a3825286/Sums-over-primitive-sets-with-a-fixed-number-of-prime-factors.pdf
 * double logn = log(n), loglogn = log(logn);
 * double lim = (1.0989 * n * pow(loglogn + 1.1174, k-1)) / (factorial(k-1)*logn);
 */


/************************************  nth  ***********************************/

UV nth_omega_prime(uint32_t k, UV n) {
  UV lo, hi;

  if (n == 0) return 0;
  if (k == 0) return (n == 1) ? 1 : 0;

  if (k > 15 || n > max_omega_prime_count(k)) return 0;

  lo = pn_primorial(k);
  if (lo == 0) return 0;
  if (n == 1) return lo;

  if (k == 1) {
    hi = nth_prime(n);
    if (hi == 0) hi = max_nth_omega_prime(1);
    lo = hi >> 1;  /* We could do better */
  } else {
    hi = 0;        /* Let the interpolation routine find it */
  }
  hi = inverse_interpolate_k(lo, hi, n, k, &opce, 600);
  if (hi == 0) return 0;

  while (!is_omega_prime(k,hi))
    hi--;
  /* if (omega_prime_count(k,hi) != n) croak("bad nth"); */
  return hi;
}
