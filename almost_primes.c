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
#include "prime_count_cache.h"
#include "semi_primes.h"
#include "factor.h"   /* for omega primes */
#include "inverse_interpolate.h"
#include "almost_primes.h"

/******************************************************************************/
/*                                KAP UTILITY                                 */
/******************************************************************************/

#if BITS_PER_WORD == 32
static uint32_t const _pow3[21] = {1,3,9,27,81,243,729,2187,6561,19683,59049,177147,531441,1594323,4782969,14348907,43046721,129140163,387420489,1162261467,3486784401U};
#else
static UV const _pow3[41] = {1,3,9,27,81,243,729,2187,6561,19683,59049,177147,531441,1594323,4782969,14348907,43046721,129140163,387420489,1162261467,3486784401U,UVCONST(10460353203),UVCONST(31381059609),UVCONST(94143178827),UVCONST(282429536481),UVCONST(847288609443),UVCONST(2541865828329),UVCONST(7625597484987),UVCONST(22876792454961),UVCONST(68630377364883),UVCONST(205891132094649),UVCONST(617673396283947),UVCONST(1853020188851841),UVCONST(5559060566555523),UVCONST(16677181699666569),UVCONST(50031545098999707),UVCONST(150094635296999121),UVCONST(450283905890997363),UVCONST(1350851717672992089),UVCONST(4052555153018976267),UVCONST(12157665459056928801)};
#endif
#define A078843_MAX_K 49
static const uint32_t _first_3[A078843_MAX_K+1] = {1, 2, 3, 5, 8, 14, 23, 39, 64, 103, 169, 269, 427, 676, 1065, 1669, 2628, 4104, 6414, 10023, 15608, 24281, 37733, 58503, 90616, 140187, 216625, 334527, 516126, 795632, 1225641, 1886570, 2901796, 4460359, 6851532, 10518476, 16138642, 24748319, 37932129, 58110457, 88981343, 136192537, 208364721, 318653143, 487128905, 744398307, 1137129971, 1736461477, 2650785552U, 4045250962U};

/* For all n <= hi, we can get the same results using 2*result with lower k */
static uint32_t reduce_k_for_n(uint32_t k, UV n) {
  uint32_t r = 0;
  if (k <= 1 || k >= BITS_PER_WORD) return 0;
  if (k > MPU_MAX_POW3)    /* Larger n would not fit in a UV type */
    r = k-MPU_MAX_POW3;
  while ((k-r) > 1 && (n>>r) < _pow3[k-r])
    r++;
  return r;
}

/* Least r s.t.  almost_prime_count(k, n)  =  almost_prime_count(k-r, n >> r) */
static void reduce_prime_count_factor(uint32_t *pk, UV *n) {
  uint32_t k = *pk, r = 0;
  if (k > MPU_MAX_POW3)    /* Larger n would not fit in a UV type */
    r = k-MPU_MAX_POW3;
  while (k >= r && ((*n)>>r) < _pow3[k-r])
    r++;
  /* Reduce */
  if (r > 0) {
    *pk -= r;
    *n >>= r;
  }
}

/* Least r s.t.  nth_almost_prime(k,n)  =  nth_almost_prime(k-r,n) << r */
static uint32_t reduce_nth_factor(uint32_t k, UV n) {
  uint32_t r = 0;
  if (k <= 1 || k >= BITS_PER_WORD) return 0;
  if (k > A078843_MAX_K) {
    if (n >= _first_3[A078843_MAX_K])
      return 0;
    r = k-A078843_MAX_K+1;
  }
  while (n < _first_3[k-r])
    r++;
  return r;
}


/* This could be easily extended to 16 or 32 */
static UV _fast_small_nth_almost_prime(uint32_t k, UV n) {
  static const uint8_t semi[8] = {0, 4,  6,  9, 10, 14, 15, 21};
  static const uint8_t mult[8] = {0, 8, 12, 18, 20, 27, 28, 30};
  MPUassert(n < 8 && k >= 2, "Fast small nth almost prime out of range");
  if (k == 2) return semi[n];
  return mult[n] * (UVCONST(1) << (k-3));
}

static void _almost_prime_count_bounds(UV *lower, UV *upper, uint32_t k, UV n);

#if 0   /* Not currently used */
/* Somewhere around k=20 it is faster to do:
 *    return nth_almost_prime(h, 1+almost_prime_count(k,n));
 */
static UV _next_almost_prime(uint32_t k, UV n) {
  while (!is_almost_prime(k, ++n))
    ;
  return n;
}
static UV _prev_almost_semiprime(uint32_t k, UV n) {
  while (!is_almost_prime(k, --n))
    ;
  return n;
}
#endif


/******************************************************************************/
/*                                KAP COUNT                                   */
/******************************************************************************/


/* almost_prime_count(k,n) is the main interface, it will call the recursive
 * function _cs(), with the terminal function _final_sum(). */

#define CACHED_PC(cache,n) prime_count_cache_lookup(cache,n)

/* for Pi from Pi to isqrt(N/Pi) [pc[n/Pi]-idx(Pi)+1] */
/* semiprime count = _final_sum(n, 1, 2, cache); */
/* 3-almost prime count = sum(Pj < icbrt(n) of _final_sum(n, Pj, Pj, cache); */
static UV _final_sum(UV n, UV pdiv, UV lo, void *cache) {
  UV s = 0, hi = isqrt(n/pdiv);
  UV j = CACHED_PC(cache, lo) - 1;  /* IDX(Pi) */

  if (hi-lo < 500) {
    SIMPLE_FOR_EACH_PRIME(lo, hi) {
      s += CACHED_PC(cache,n/(pdiv*p)) - j++;
    } END_SIMPLE_FOR_EACH_PRIME
    return s;
  }

  if ((lo <= 2) && (hi >= 2)) s += CACHED_PC(cache,n/(pdiv*2)) - j++;
  if ((lo <= 3) && (hi >= 3)) s += CACHED_PC(cache,n/(pdiv*3)) - j++;
  if ((lo <= 5) && (hi >= 5)) s += CACHED_PC(cache,n/(pdiv*5)) - j++;
  if (lo < 7) lo = 7;
  if (lo <= hi) {
    unsigned char* segment;
    UV seg_base, seg_low, seg_high;
    void* ctx = start_segment_primes(lo, hi, &segment);
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
        s += CACHED_PC(cache,n/(pdiv*p)) - j++;
      END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
  }
  return s;
}

static UV _cs(UV n, UV pdiv, UV lo, uint32_t k, void *cache) {
  UV count = 0;

  if (k == 2)
    return _final_sum(n, pdiv, lo, cache);

  SIMPLE_FOR_EACH_PRIME(lo, rootint(n/pdiv,k)) {
    if (k == 3)  count += _final_sum(n, pdiv*p, p, cache);
    else         count += _cs(n, pdiv*p, p, k-1, cache);
  } END_SIMPLE_FOR_EACH_PRIME

  return count;
}

UV almost_prime_count(uint32_t k, UV n)
{
  void* cache;
  UV count, csize;

  if (k == 0) return (n >= 1);
  if (k >= BITS_PER_WORD || (n >> k) == 0) return 0;
  reduce_prime_count_factor(&k, &n); /* Reduce to lower k,n if possible */

  if (k == 1) return LMO_prime_count(n);
  if (k == 2) return semiprime_count(0,n);
  if (n <  3*(UVCONST(1) << (k-1))) return 1;
  if (n <  9*(UVCONST(1) << (k-2))) return 2;
  if (n < 10*(UVCONST(1) << (k-2))) return 3;

  /* Decide how much we will cache prime counts.
   *
   * n/(1UL << (k+M)) has 0,1,2,7,15,37,84,187,... lookups for M=-2,-1,0,...
   * The number of non-cached counts performed follows OEIS A052130. */

  csize = n / (1UL << (k-2));
  if (csize < 32) csize = 32;
  if (csize >       16UL*1024)  csize = n / (1UL << (k+2));  /*  15 */
  if (csize >      128UL*1024)  csize = n / (1UL << (k+4));  /*  84 */
  if (csize >   1UL*1024*1024)  csize = n / (1UL << (k+6));  /* 421 */
  if (((csize >> 16) >> 16) >= 3)  csize >>= 1;

  cache = prime_count_cache_create( csize );
  count = _cs(n, 1, 2, k, cache);
  prime_count_cache_destroy(cache);
  return count;
}

#if 0
UV almost_prime_count(uint32_t k, UV lo, UV hi) {
  if (k == 0) return (lo <= 1 || hi >= 1);
  if (k == 1) return prime_count(lo, hi);
  if (k == 2) return semiprime_count(lo, hi);
  /* See semiprime_count.  Possibly clever solutions for small ranges. */
  return _almost_prime_count(k, hi) - _almost_prime_count(k, lo-1);
}
#endif

UV almost_prime_count_approx(uint32_t k, UV n) {
  UV lo, hi;

  if (k == 0) return (n >= 1);
  if (k >= BITS_PER_WORD || (n >> k) == 0) return 0;
  reduce_prime_count_factor(&k, &n); /* Reduce to lower k,n if possible */

  if (k == 1) return prime_count_approx(n);
  if (k == 2) return semiprime_count_approx(n);
  if (n <  3*(UVCONST(1) << (k-1))) return 1;
  if (n <  9*(UVCONST(1) << (k-2))) return 2;
  if (n < 10*(UVCONST(1) << (k-2))) return 3;

  if (k == 3 && n < 102) {
    unsigned char const sm3[19] = {27,28,30,42,44,45,50,52,63,66,68,70,75,76,78,92,98,99};
    for (lo=0; lo < 19; lo++)
      if (n < sm3[lo])
        break;
    return 4+lo;
  }

  _almost_prime_count_bounds(&lo, &hi, k, n);

  if (k == 3) {  /* Much better fit for k=3. */
    double x = n, logx = log(x), loglogx = log(logx);
    double a = 1.0, s = 2.0;
    UV est;
    if      (x <=      638) { s = 1.554688; a = 0.865814; }
    else if (x <=     1544) { s = 1.050000; a = 0.822256; }
    else if (x <=     1927) { s = 0.625000; a = 0.791747; }
    else if (x <=   486586) { s = 2.865611; a = 1.004090; }
    else if (x <=  1913680) { s = 2.790963; a = 0.999618; }
    else if (x <= 22347532) { s = 2.719238; a = 0.995635; }
    else if (x <= 2.95e8)   { s = 2.584473; a = 0.988802; }
    else if (x <= 4.20e9)   { s = 2.457108; a = 0.983098; }
    else if (x <= 7.07e10)  { s = 2.352818; a = 0.978931; }
    else if (x <= 1.36e12)  { s = 2.269745; a = 0.975953; }
    else if (x <= 4.1e13)   { s = 2.203002; a = 0.973796; }
    else if (x <= 9.2e14)   { s = 2.148463; a = 0.972213; }
    else                    { s = 2.119279; a = 0.971438; }
    est = 0.5 * a * x * ((loglogx+0.26153)*(loglogx+0.26153)) / (logx+s)+0.5;
    if      (est < lo) est = lo;
    else if (est > hi) est = hi;
    return est;
  }

  /* We should look at (1) better bounds, (2) better weighting here */
  /* return lo + (hi-lo)/2; */
  return lo + (hi-lo) * 0.8;
}


/******************************************************************************/
/*                                  NTH KAP                                   */
/******************************************************************************/

#if 0
/* Asymptotic estimate for the nth k-almost prime */
static double _asymptotic_nth(uint32_t k, UV n) {
  uint32_t i;  double x, logn, loglogn;
  if (k == 0 || n == 0) return 0;
  if (n == 1) return UVCONST(1) << k;
  logn = log(n);
  loglogn = log(logn);
  x = n * logn;
  for (i = 1; i < k; i++)
    x *= (double)i / loglogn;
  return x;
}
#endif

static UV apcu(UV mid, UV k) { return almost_prime_count_upper(k, mid); }
static UV apcl(UV mid, UV k) { return almost_prime_count_lower(k, mid); }
static UV apca(UV mid, UV k) { return almost_prime_count_approx(k, mid); }
static UV apce(UV mid, UV k) { return almost_prime_count(k, mid); }

UV nth_almost_prime_upper(uint32_t k, UV n) {
  UV r, max, lo, hi;

  if (n == 0) return 0;
  if (k == 0) return (n == 1) ? 1 : 0;
  if (k == 1) return nth_prime_upper(n);
  if (n < 8) return _fast_small_nth_almost_prime(k, n);

  max = max_almost_prime_count(k);
  if (n > max) return 0;
  r = reduce_nth_factor(k,n);
  if (r > 0)  return nth_almost_prime_upper(k-r, n) << r;

  /* We start out with the literal min and max because we have NO idea. */
  lo = 5 * (UVCONST(1) << k);  /* For k >= 1 and n >= 8 */
  hi = max_nth_almost_prime(k);

  return inverse_interpolate(lo, hi, n, k, &apcl, 0);
}

UV nth_almost_prime_lower(uint32_t k, UV n) {
  UV r, max, lo, hi;

  if (n == 0) return 0;
  if (k == 0) return (n == 1) ? 1 : 0;
  if (k == 1) return nth_prime_lower(n);
  if (n < 8) return _fast_small_nth_almost_prime(k, n);

  max = max_almost_prime_count(k);
  if (n > max) return 0;
  r = reduce_nth_factor(k,n);
  if (r > 0)  return nth_almost_prime_lower(k-r, n) << r;

  /* We start out with the literal min and max because we have NO idea. */
  lo = 5 * (UVCONST(1) << k);  /* For k >= 1 and n >= 8 */
  hi = max_nth_almost_prime(k);

  return inverse_interpolate(lo, hi, n, k, &apcu, 0);
}

UV nth_almost_prime_approx(uint32_t k, UV n) {
  UV max, lo, hi;

  if (n == 0) return 0;
  if (k == 0) return (n == 1) ? 1 : 0;
  if (k == 1) return nth_prime_approx(n);
  if (k == 2) return nth_semiprime_approx(n);

  max = max_almost_prime_count(k);
  if (n > max) return 0;
  /* We could reduce but really no reason to do it */

  if (n < 8) return _fast_small_nth_almost_prime(k,n);
  lo = nth_almost_prime_lower(k,n);
  hi = nth_almost_prime_upper(k,n);

  return inverse_interpolate(lo, hi, n, k, &apca, 0);
}

UV nth_almost_prime(uint32_t k, UV n) {
  UV r, max, lo, hi;

  if (n == 0) return 0;
  if (k == 0) return (n == 1) ? 1 : 0;
  if (k == 1) return nth_prime(n);
  if (k == 2) return nth_semiprime(n);

  max = max_almost_prime_count(k);
  if (n > max) return 0;

  /* For k >= 3 and small n we can answer this quickly. */
  if (n < 8) return _fast_small_nth_almost_prime(k,n);
  r = reduce_nth_factor(k,n);
  if (r > 0)  return nth_almost_prime(k-r,n) << r;
  /* NOTE: given n a 64-bit integer, k always <= 40 after reduction */

  lo = nth_almost_prime_lower(k,n);
  hi = nth_almost_prime_upper(k,n);
  hi = inverse_interpolate(lo, hi, n, k, &apce, 60000);
  while (!is_almost_prime(k,hi))
    hi--;
  return hi;
}


/******************************************************************************/
/*                                  Bounds                                    */
/******************************************************************************/

/* Bayless et al. (2018) and Kinlaw (2019) are main references.
 * Note that they use Pi_k(x) to mean square-free numbers, and
 * Tau_k(x) to mean the general count like we use.
 * They also have results for k = 2,3,4 only.
 * Also see https://archimede.mat.ulaval.ca/MAINE-QUEBEC/19/Kinlaw19.pdf.
 *
 * We split into three ranges:
 *    1    - 2^20  complete computations
 *    2^20 - 2^32  complete computations
 *    2^32 - 2^64  correct upper for k=2,3,4.  correct lower for k=2.
 *                 estimates elsewhere.
 *
 */

static const double _upper_20[13] = {0,0, 1.006,0.7385,0.6830,0.5940,0.3596,0.2227,0.1439, 0.09786,0.07017,0.05304,0.04205};
static const double _upper_32[21] = {0,0, 1.013,0.8095,0.7486,
  /*  5-12 */ 0.6467,0.3984,0.2465,0.1573,0.105,0.07374,0.05462,0.04275,
  /* 13-20 */ 0.03552,0.03093,0.02814,0.02665,0.02619,0.02667,0.02808,0.03054};
static const double _upper_64[41] = {0,0, 1.028, 1.028, 1.3043,/* <--corrrect */
  /* 5-12 */
  0.70,0.45,0.30,0.18,0.12,0.09,0.07,0.05,
  /* 13-20 */
  0.05, 0.03007, 0.02711, 0.02555, 0.02505, 0.02554, 0.02699, 0.02955,
  /* 21-28 */
  0.03331, 0.03879, 0.04643, 0.05696, 0.07167, 0.09293, 0.1234, 0.1671,
  /* 29-31 */
  0.2318,0.3288,0.477,
  /* 32-40 */
  0.7071,1.071,1.654,2.607,4.188,6.857,11.45,19.43,33.58,
};

static const double _lower_20[13] = {0,0, 0.7716, 0.3994, 0.4999,
  /*  5-12 */ 0.4778,0.2699,0.1631,0.1050,0.0718,0.05185,0.03943,0.03149};
static const double _lower_32[21] = {0,0, 1.004,  0.7383, 0.6828,
  /*  5-12 */ 0.5939,0.3594,0.2222,0.1438,0.09754,0.06981,0.05245,0.04151,
  /* 13-20 */ 0.03426,0.0290,0.02617,0.02344,0.02183,0.01972,0.02073,0.02252 };
static const double _lower_64[41] = {0,0, 1.011,  0.8093, 0.7484,
  /* 5-12 */
  0.6465,0.3982,0.2463,0.1571,0.1048,0.07363,0.0545,0.042,
  /* 13-20 */
  0.033,0.025,0.024,0.021,0.020,0.019,0.019,0.021,
  /* 21-28 */
  0.022,0.027,0.029,0.035,0.040,0.057,0.0115,0.10,
  /* 29-31 */
  0.13,0.20,0.29,
  /* 32-40 */
  0.44,0.67,1.06,1.71,2.8,4.7,8.0,13.89,23.98,
};

static void _almost_prime_count_bounds(UV *lower, UV *upper, uint32_t k, UV n) {
  double x, logx, loglogx, logplus, multl, multu, boundl, boundu;
  UV max;
  uint32_t i;

  if (k >= BITS_PER_WORD || (n >> k) == 0) { *lower = *upper = 0; return; }
  reduce_prime_count_factor(&k, &n); /* Reduce to lower k,n if possible */
  if (k == 0) { *lower = *upper = (n >= 1); return; }
  if (k == 1) { *lower = prime_count_lower(n); *upper = prime_count_upper(n); return; }
  if (n <  3*(UVCONST(1) << (k-1))) { *lower = *upper = 1; return; }
  if (n <  9*(UVCONST(1) << (k-2))) { *lower = *upper = 2; return; }
  if (n < 10*(UVCONST(1) << (k-2))) { *lower = *upper = 3; return; }

  max = max_almost_prime_count(k);
  if (n >= max_nth_almost_prime(k) && max != UV_MAX)
    { *lower = *upper = max; return; }

  x = (double) n;
  logx = log(x);
  loglogx = log(logx);
  logplus = loglogx + 0.26153;

  /* Select the appropriate table for n's range.
   * 20/32/64-bit n will always reduce k to these limits. */
  if (n <= 1048575U) {
    MPUassert(k <= 12, "almost prime count: 20-bit n doesn't exceed k 12");
    multu = _upper_20[k];  multl = _lower_20[k];
  } else if (n <= 4294967295U) {
    MPUassert(k <= 20, "almost prime count: 32-bit n doesn't exceed k 20");
    multu = _upper_32[k];  multl = _lower_32[k];
  } else {
    MPUassert(k <= 40, "almost prime count: after reduction, k <= 40");
    multu = _upper_64[k];  multl = _lower_64[k];
  }

  if (k == 2) {
    boundl = boundu = x * (loglogx + 0.261536) / logx;
    if (x >= 1e12) {
      boundl = x*(loglogx+0.1769)/logx * (1+0.4232/logx);
      multl = 1;
    }
  } else if (k == 3) {
    boundu = x * (logplus*logplus + 1.055852) / (2*logx);
    /* Kinlaw (2019) Theorem 1 (with 1.000) */
    boundl = x * loglogx * loglogx / (2*logx);
    multl = (x <=      500194)  ?  0.8418
          : (x <= 3184393786U)  ?  1.0000
          :                        1.04;
  } else if (k == 4) {
    /* Bayless Theorem 5.4 part 1 (with multu = 1.3043) */
    boundl = boundu = x * logplus*logplus*logplus / (6*logx);
    /* Bayless Theorem 5.4 part 2 */
    if (x > 1e12) {
      boundu += 0.511977 * x * (log(log(x/4)) + 0.261536) / logx;
      multu = 1.028;
    }
  } else {
    /* Completely empirical and by no means optimal.
     * It is easy and seems fairly reasonable through k=20 or so.
     */
    boundl = x / logx;
    logplus = loglogx + (log(k)*log(log(k))-0.504377); /* k=5 => 0.26153 */
    for (i = 1; i < k; i++)
      boundl *= logplus / (double)i;
    boundu = boundl;
  }

#if 0
  printf("  lower: %lf * %lf = %lf\n", boundl, multl, boundl*multl);
  printf("  upper: %lf * %lf = %lf\n", boundu, multu, boundu*multu);
  printf("  max: %lu\n", max);
#endif
  boundl *= multl;
  boundu *= multu;

  *lower = (boundl >= UV_MAX || (max > 0 && boundl > max)) ? max : (UV)boundl;
  *upper = (boundu >= UV_MAX || (max > 0 && boundu > max)) ? max : (UV)(boundu+1.0);
}

UV almost_prime_count_upper(uint32_t k, UV n) {
  UV l, u;
  _almost_prime_count_bounds(&l, &u, k, n);
  return u;
}

UV almost_prime_count_lower(uint32_t k, UV n) {
  UV l, u;
  _almost_prime_count_bounds(&l, &u, k, n);
  return l;
}

UV max_nth_almost_prime(uint32_t k) {
#if BITS_PER_WORD == 32
  static const UV offset[32] = {0,4,1,9,5,0,7,47,31,3,15,511,1263,5119,1023,255,23295,2559,4095,126975,16383,262143,2359295,2097151,5767167,1048575,33554431,16777215,100663295,67108863,268435455,1073741823};
#else
  static const UV offset[64] = {0,58,14,2,4,
    /*  5-12 */  3,17,0,1,195,51,127,63,
    /* 13-22 */  767,1535,511,255,15,8191,1023,83967,16383,111615,
    /* 23-32 */  557055,2097151,524287,65535,1048575,6553599,33554431,4194303,671088639,16777215,
    /* 33-63 */  UVCONST(536870911),UVCONST(2684354559),UVCONST(2147483647),
    UVCONST(25769803775),UVCONST(4294967295),UVCONST(268435455),
    UVCONST(206158430207),UVCONST(137438953471),UVCONST(17179869183),
    UVCONST(68719476735),UVCONST(2199023255551),UVCONST(5428838662143),
    UVCONST(21990232555519),UVCONST(4398046511103),UVCONST(1099511627775),
    UVCONST(100055558127615),UVCONST(10995116277759),UVCONST(17592186044415),
    UVCONST(545357767376895),UVCONST(70368744177663),UVCONST(1125899906842623),
    UVCONST(10133099161583615),UVCONST(9007199254740991),
    UVCONST(24769797950537727),UVCONST(4503599627370495),
    UVCONST(144115188075855871),UVCONST(72057594037927935),
    UVCONST(432345564227567615),UVCONST(288230376151711743),
    UVCONST(1152921504606846975),UVCONST(4611686018427387903)
  };
#endif
  if (k >= BITS_PER_WORD) return 0;
  return UV_MAX - offset[k];
}

UV max_almost_prime_count(uint32_t k) {
#if BITS_PER_WORD == 32
  static const UV max[32] = {1,203280221,658662065,967785236,916899721,665533848,410630859,229679168,121092503,61600699,30653019,15043269,7315315,3535071,1700690,814699,389357,185245,87964,41599,19611,9184,4283,2001,914,421,187,84,37,15,7,2,0};
#else
  static const UV max[64] = {1,
    UVCONST(425656284035217743),  /* max prime count */
    UVCONST(1701500000000000000), /* max semiprime count (estimate) */
    UVCONST(3738241841329689241), /* max 3-almost-prime count (upper limit) */
    UVCONST(5626791887645067628), /* max 4-almost-prime count (upper limit) */
    /* 5-12 */  0,0,0,0,0,UVCONST(263157990621533964),UVCONST(135128109904869290),UVCONST(68283616225825256),
    /* 13-22 */ UVCONST(34151861008771016),UVCONST(16967424859951587),UVCONST(8393048221327186),UVCONST(4139595949113890),UVCONST(2037655246635364),UVCONST(1001591348315641),UVCONST(491808604962296),UVCONST(241293656953012),UVCONST(118304122014405),UVCONST(57968649799947),
    /* 23-32 */ UVCONST(28388662714236),UVCONST(13895161400556),UVCONST(6797526392535),UVCONST(3323560145881),UVCONST(1624109166018),UVCONST(793189260998),UVCONST(387148515886),UVCONST(188844769357),UVCONST(92054377509),UVCONST(44841620426),
    /* 33-63 */ UVCONST(21827124353),UVCONST(10616326552),UVCONST(5159281045),UVCONST(2505087309),1215204383,588891145,285076316,137840686,66567488,32103728,15460810,7433670,3567978,1709640,817053,389954,185387,87993,41604,19611,9184,4283,2001,914,421,187,84,37,15,7,2
  };
#endif
  if (k >= BITS_PER_WORD) return 0;
  if (max[k] == 0) return UV_MAX;  /* TODO: find these */
  return max[k];
}


/******************************************************************************/
/*                              Construction                                  */
/******************************************************************************/

/* There are a few options for constructing KAPs without sieving/factoring.
 *
 *   1) we can make an iterator that recursively constructs them using
 *      a prime list and a k-1 iterator.  This is a generalization of
 *      Dijkstra's Hamming number algorithm.
 *
 *   2) Given a range [lo,hi], We can ask for all k-1 kaps less than hi/2,
 *      then multiply through by primes to see which fall in our range.
 *
 * Each of these (and sieving) is limited in some ways.  For example, #1
 * can output 500-almost-primes quite rapidly which some other methods have
 * trouble with, even with all the calculation in Perl.  But it rapidly
 * slows down with increasing n.
 *
 * I suspect there are far more efficient methods.
 */

static void _tidy_list(UV **list, UV *Lsize, UV *count, int minimal) {
  UV *L = *list;

  if (*count > 1) {
    UV i, j;
    qsort(L, *count, sizeof(UV), _numcmp);
    for (j = 0, i = 1; i < *count; i++) {
      if (L[i] != L[j])
        L[++j] = L[i];
    }
    *count = j+1;
  }
  if (minimal) {
    *Lsize = *count;
    Renew(*list, *Lsize, UV);
  } else if (*count * 1.5 > *Lsize) {
    *Lsize = *count * 2 + 100;
    Renew(*list, *Lsize, UV);
  }
}

UV range_construct_almost_prime(UV** list, uint32_t k, UV lo, UV hi) {
  UV *L, minkap1, lastprime, count = 0;

  if (k == 0 || k >= BITS_PER_WORD) { *list = 0; return 0; }
  if ((lo >> k) == 0) lo = UVCONST(1) << k;
  if (hi > max_nth_almost_prime(k)) hi = max_nth_almost_prime(k);
  if (lo > hi) { *list = 0; return 0; }

  if (k == 1) return range_prime_sieve(list, lo, hi);
  if (k == 2) return range_semiprime_sieve(list, lo, hi);
  /* if (k <= 5) return range_almost_prime_sieve(list, k, lo, hi); */

  minkap1 = 1 << (k-1);
  lastprime = hi / minkap1;  /* lastprime = prev_prime(lastprime+1); */

  {
    UV i, Lsize;
    UV *lkap1, nkap1=range_construct_almost_prime(&lkap1, k-1, minkap1, hi>>1);

    /* Now multiply through exhaustively. */
    Lsize = nkap1*4 + 100;
    New(0, L, Lsize, UV);

    START_DO_FOR_EACH_PRIME(2, lastprime) {
      for (i = 0; i < nkap1; i++) {
        UV prod = p * lkap1[i];
        if (prod < lo) continue;
        if (prod > hi) break;
        if (count >= Lsize)
          _tidy_list(&L, &Lsize, &count, 0);
        L[count++] = prod;
      }
    } END_DO_FOR_EACH_PRIME
    _tidy_list(&L, &Lsize, &count, 1);
    Safefree(lkap1);
  }
  *list = L;
  return count;
}

UV range_almost_prime_sieve(UV** list, uint32_t k, UV slo, UV shi)
{
  UV *S, Ssize, i, j, count;
  const UV thresh_pred = 40;

  if (k == 0 || k >= BITS_PER_WORD) { *list = 0; return 0; }
  if ((slo >> k) == 0) slo = UVCONST(1) << k;
  if (shi > max_nth_almost_prime(k)) shi = max_nth_almost_prime(k);
  if (slo > shi) { *list = 0; return 0; }

#if 1
  if (shi-slo+1 < thresh_pred) {
    Ssize = 3 + (thresh_pred >> 1);
    New(0, S, Ssize, UV);
    for (i = 0, j = 0; i < shi-slo+1; i++)
      if (is_almost_prime(k, slo+i))
        S[j++] = slo+i;
    *list = S;
    return j;
  }
#endif

  if (k == 1) return range_prime_sieve(list, slo, shi);
  if (k == 2) return range_semiprime_sieve(list, slo, shi);

  /* See if we can reduce k.
   * If for all possible kap from 1 to shi, ap(k,n) = 2*ap(k-1,n), then
   * sieve for k-1 from lo/2 to hi/2+1.
   * For large k this can continue even further so we might reduce a lot.
   */
  {
    uint32_t r = reduce_k_for_n(k, shi);
    if (r > 0) {
      UV lo = (slo >> r) + (((slo >> r) << r) < slo);
      UV hi = shi >> r;
      count = range_almost_prime_sieve(&S, k-r, lo, hi);
      for (i = 0; i < count; i++)
        S[i] <<= r;
      *list = S;
      return count;
    }
  }

  Ssize = (almost_prime_count_approx(k,shi) - almost_prime_count_approx(k,slo) + 1) * 1.2 + 100;
  New(0, S, Ssize, UV);

  /* Do a range nfactor sieve in small windows, with one optimization.
   *
   * We know that we are looking for numbers with k factors, hence after
   * looking for small factors we could get a remainder R as large as:
   * 2 x 2 x ... x R where R could be prime or semiprime.  Hence we can
   * reduce the sieve limit somewhat.  Effectively we are sieving to the
   * maximum possible *second largest* factor for a k-almost-prime,
   * allowing us to correctly decide whether R is prime or semiprime
   * (if it has >= k factors).
   *
   * This isn't a big deal for small k, but it's a big impact for high k.
   *
   * I still think there should be a better way to do this for high k.
   * Is there any way to do this just sieving to rootint(hi,k+1)?
   * Given hi=10^6:
   *   k=3 =>         97^3,   2x701^2,   2x2x249989
   *   k=4 => 31^4, 2x79^3, 2x2x499^2, 2x2x2x724991
   */
  {
    unsigned char* nf;
    UV const segsize = 65536*4;
    UV *N, lo, hi, range;
    UV kdiv = (k < 3) ? UVCONST(1) : (UVCONST(1) << (k-2));

    New(0, nf, segsize+1, unsigned char);
    New(0, N, segsize+1, UV);
    prime_precalc(isqrt(shi/kdiv));
    count = 0;

    for (lo = slo; lo <= shi && lo >= slo; lo = hi+1) {
      hi = lo+segsize-1;
      if (hi > shi  ||  hi < lo)  hi = shi;
      range = hi - lo + 1;
      memset(nf, 0, range);
      for (i = lo; i <= hi; i++) {
        if (!(i&1) && i >= 2) {
          unsigned char nz = ctz(i);
          nf[i-lo] = nz;
          N[i-lo] = UVCONST(1) << nz;
        } else
          N[i-lo] = 1;
      }
      START_DO_FOR_EACH_PRIME(3, isqrt(hi/kdiv)) {
        UV pk, maxpk = UV_MAX/p; \
        for (i = P_GT_LO_0(p,p,lo); i < range; i += p)
          { N[i] *= p;  nf[i]++; }
        for (pk = p*p; pk <= hi; pk *= p) {
          for (i = P_GT_LO_0(pk,pk,lo); i < range; i += pk)
            { N[i] *= p;  nf[i]++; }
          if (pk >= maxpk) break;  /* Overflow protection */
        }
      } END_DO_FOR_EACH_PRIME
      for (i = 0; i < range; i++) {
        if (N[i] < (lo+i))
          nf[i]++;
        if (nf[i] == k) {
          if (count >= Ssize)
            Renew(S, Ssize += 10000, UV);
          S[count++] = i+lo;
        }
      }
    }
    Safefree(N);
    Safefree(nf);
  }
  *list = S;
  return count;
}

/******************************************************************************/
/*                              OMEGA PRIMES                                  */
/******************************************************************************/

/* See https://arxiv.org/pdf/2006.16491.pdf page 12 for a brief note */

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


/* Recursive nested summation over prime powers */
static UV _opcn(uint32_t k, UV n, UV qi, UV PQ, uint32_t* pr, UV numprimes) {
  UV nPQ = n / PQ;
  UV i, maxi, sum;

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
      for (R = r;  R < n / (PQ*r);  R *= r) {    /* TODO: overflow */
        sum += _opcn(k-1, n, i+1, PQ*R, pr, numprimes);
      }
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
  lo = primorial(nth_prime(k));
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
  lo = primorial(nth_prime(k));
  if (lo == 0 || n < lo) return 0;

  npr = range_prime_sieve_32(&pr, isqrt(n), 0);  /* p[0]=2, p[1]=3,... */
  sum = _opcn(k, n, 0, 1, pr, npr);
  Safefree(pr);
  return sum;
#endif
}
