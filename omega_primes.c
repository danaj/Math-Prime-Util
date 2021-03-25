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
#include "lmo.h"
#include "prime_nth_count.h"
#include "factor.h"
#include "inverse_interpolate.h"
#include "omega_primes.h"

/******************************************************************************/
/*                              OMEGA PRIMES                                  */
/******************************************************************************/

/* See https://arxiv.org/pdf/2006.16491.pdf page 12 for a brief note */

/* For the interpolation */
static UV opce(UV mid, UV k) { return omega_prime_count(k, mid); }

UV range_omega_prime_sieve(UV** ret, uint32_t k, UV lo, UV hi) {
  UV i, lmax = 0, n = 0;
  UV* l = 0;
  unsigned char *nf;

  if (hi < lo) croak("range_omega_prime_sieve error hi %"UVuf" < lo %"UVuf"\n",hi,lo);
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


/* Recursive nested summation over prime powers */
static UV _opcn(uint32_t k, UV n, UV qi, UV PQ, uint32_t* pr, UV numprimes) {
  UV i, maxi, sum, nPQ = n / PQ;

  if (k == 1) {
    sum = prime_power_count(nPQ);
    for (i = 0; i < qi; i++)
      sum -= logint(nPQ, pr[i]);
  } else {
    sum = 0;
    maxi = prime_count(0, rootint(nPQ,k));
    if (maxi > numprimes) croak("not enough primes in opcn");
    for (i = qi;  i < maxi;  i++) {
      UV R, r = pr[i];
      for (R = r;  R < n / (PQ*r);  R *= r)    /* TODO: overflow */
        sum += _opcn(k-1, n, i+1, PQ*R, pr, numprimes);
    }
  }
  return sum;
}

UV omega_prime_count(uint32_t k, UV n)
{
#if 0
  UV const incr = 500000;
  UV i, lo, hi, count;

  if (k == 0) return (n >= 1);
  if (k == 1) return prime_power_count(n);

  /* The first k-omega-prime is primorial(p_k) (ignoring zero for k=1) */
  lo = pn_primorial(k);
  if (lo == 0) return 0;

  prime_precalc(isqrt(n));

  for (count = 0;  lo <= n;  lo = hi+1) {
    hi = lo + incr - 1;
    if (hi > n) hi = n;
    {
      unsigned char* nf = range_nfactor_sieve(lo, hi, 0);
      for (i = 0; i < hi-lo+1; i++)
        if (nf[i] == k)
          count++;
      Safefree(nf);
    }
  }
  /* Do not count prime_omega(0) = 1 */
  return count;
#else
  uint32_t* pr;
  UV npr, sum, lo;

  if (k == 0) return (n >= 1);
  if (k == 1) return prime_power_count(n);

  /* The first k-omega-prime is primorial(p_k) (ignoring zero for k=1) */
  lo = pn_primorial(k);
  if (lo == 0 || n < lo) return 0;

  npr = range_prime_sieve_32(&pr, isqrt(n), 0);  /* p[0]=2, p[1]=3,... */
  sum = _opcn(k, n, 0, 1, pr, npr);
  Safefree(pr);
  return sum;
#endif
}

/* An upper bound for the omega prime count, when n >= 10^12 is shown in
 * Bayless,Kinlaw,Klyve 2019, page 4
 * https://www.researchgate.net/profile/Paul-Kinlaw/publication/329788487_Sums_over_primitive_sets_with_a_fixed_number_of_prime_factors/links/5c44103d92851c22a3825286/Sums-over-primitive-sets-with-a-fixed-number-of-prime-factors.pdf
 * double logn = log(n), loglogn = log(logn);
 * double lim = (1.0989 * n * pow(loglogn + 1.1174, k-1)) / (factorial(k-1)*logn);
 */

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
    if (hi == 0) hi = UV_MAX;
    lo = hi >> 1;  /* We could do better */
  } else {
    hi = lo;  /* TODO better bounds */
    while (omega_prime_count(k,hi) < n) {
      if (hi >= UV_MAX/2) { hi = UV_MAX;  break; }
      hi *= 2;
    }
  }
  if (hi > max_nth_omega_prime(k))  hi = max_nth_omega_prime(k);
  hi = inverse_interpolate(lo, hi, n, k, &opce, 600);

  while (!is_omega_prime(k,hi))
    hi--;
  /* if (omega_prime_count(k,hi) != n) croak("bad nth"); */
  return hi;
}
