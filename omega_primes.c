/******************************************************************************/
/*                              ALMOST PRIMES                                 */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ptypes.h"
#include "constants.h"
#define FUNC_isqrt 1
#define FUNC_ctz 1
#include "util.h"
#include "cache.h"
#include "sieve.h"
#include "prime_counts.h"
#include "prime_powers.h"
#include "factor.h"
#include "inverse_interpolate.h"
#include "omega_primes.h"

/******************************************************************************/
/*                              OMEGA PRIMES                                  */
/******************************************************************************/

int is_omega_prime(uint32_t k, UV n) {
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

static void _omega_prime_gen_rec(UV** kop, UV* skop, UV* nkop, uint32_t k, UV lo, UV hi, UV m, UV pstart) {
  UV v, *l = *kop, lsize = *skop, n = *nkop;

  if (k > 1) {
    SIMPLE_FOR_EACH_PRIME(pstart, rootint(hi/m, k)) {
      if ((m % p) == 0) continue;
      for (v = m*p; v <= hi; v *= p)
        if ((v*p) <= hi)
          _omega_prime_gen_rec(kop, skop, nkop, k-1, lo, hi, v, p);
    } END_SIMPLE_FOR_EACH_PRIME
    return;
  }

  START_DO_FOR_EACH_PRIME(pstart, rootint(hi/m, k)) {
    if ((m % p) == 0) continue;
    for (v = m*p; v <= hi; v *= p) {
      if (v >= lo) { // Add v to kop list
        if (n >= lsize) {
          lsize = 1 + lsize * 1.2;
          Renew(l, lsize, UV);
        }
        l[n++] = v;
      }
    }
  } END_DO_FOR_EACH_PRIME
  *kop = l;  *skop = lsize;  *nkop = n;
}

UV rec_omega_primes(UV** ret, uint32_t k, UV lo, UV hi) {
  UV min, nkop, skop, *kop;

  min = pn_primorial(k);
  if (lo < min) lo = min;

  if (hi < lo) croak("range_omega_prime_sieve error hi %"UVuf" < lo %"UVuf"\n",hi,lo);

  nkop = 0;
  skop = 256;
  New(0, *ret, skop, UV);
  _omega_prime_gen_rec(ret, &skop, &nkop, k, lo, hi, 1, 2);
  qsort(*ret, nkop, sizeof(UV), _numcmp);
  return nkop;
}


UV range_omega_prime_sieve(UV** ret, uint32_t k, UV lo, UV hi) {
  UV i, lmax = 0, n = 0;
  UV* l = 0;
  unsigned char *nf;

  if (hi < lo) croak("range_omega_prime_sieve error hi %"UVuf" < lo %"UVuf"\n",hi,lo);

  /* TODO: The recursive routine should compute primes like the count does */
  if ( ((hi-lo) > 1000000000UL) || (k >= 10 && (hi-lo) > 10000000UL) )
    return rec_omega_primes(ret, k, lo, hi);

  nf = range_nfactor_sieve(lo, hi, 0);
  if (ret != 0) {
    lmax = 1000;
    New(0, l, lmax, UV);
  }
  for (i = 0; i < hi-lo+1; i++) {
    if (nf[i] != k) continue;
    if (l != 0) {
      if (n >= lmax)  { lmax = 1 + lmax * 1.2;  Renew(l, lmax, UV); }
      l[n] = lo+i;
    }
    n++;
  }
  Safefree(nf);
  if (ret != 0)  *ret = l;
  return n;
}

/* TODO: Should make a single construct routine that calls sieve or recurse */


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
    UV r2, w, u, k, rlim;
    for (;  p <= s;  j++, p = r) {
      r = (j < numprimes)  ?  pr[j]  :  next_prime(p);
      for (t = m*p, w = n/t;  t <= n && w >= r;  t *= p, w = n/t) {
#if 1
        count += prime_count(w) - j;
        for (k = j, r2 = r, rlim = isqrt(w);
             r2 <= rlim;
             r2 = (++k < numprimes) ? pr[k] : rlim+1) {
          u = t * r2;
          do {  u *= r2;  count++;  } while (n/r2 >= u);
        }
#else
        /* This is the basic method from the definition, before optimizing */
        UV q;
        count += prime_power_count(w);
        rlim = prev_prime(r);
        for (k = 1, q = 2;
             q <= rlim;
             q = (++k < numprimes) ? pr[k-1] : nth_prime(k)) {
          count -= logint(w, q);
        }
#endif
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

  while (!is_omega_prime(k,hi))
    hi--;
  /* if (omega_prime_count(k,hi) != n) croak("bad nth"); */
  return hi;
}
