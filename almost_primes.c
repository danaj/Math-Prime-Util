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
#include "prime_count_cache.h"
#include "semi_primes.h"
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

#define CACHED_PC(cache,n) prime_count_cache_lookup(cache,n)

/* Debatably useful.  Slightly faster for small n, the same for larger. */
static UV almost3prime_count(UV n) {
  UV sum = 0, cbrtn = prev_prime(rootint(n,3)+1);
  void *cache = prime_count_cache_create( (UV)pow(n,0.72) );

  SIMPLE_FOR_EACH_PRIME(2, cbrtn) {
    UV pdiv = p, lo = p, hi = isqrt(n/pdiv);
    UV j = CACHED_PC(cache, lo) - 1;  /* IDX(Pi) */
    if ((lo <= 2) && (hi >= 2)) sum += CACHED_PC(cache,n/(pdiv*2)) - j++;
    if ((lo <= 3) && (hi >= 3)) sum += CACHED_PC(cache,n/(pdiv*3)) - j++;
    if ((lo <= 5) && (hi >= 5)) sum += CACHED_PC(cache,n/(pdiv*5)) - j++;
    if (lo < 7) lo = 7;
    if (lo <= hi) {
      unsigned char* segment;
      UV seg_base, seg_low, seg_high;
      void* ctx = start_segment_primes(lo, hi, &segment);
      while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
        START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
          sum += CACHED_PC(cache,n/(pdiv*p)) - j++;
        END_DO_FOR_EACH_SIEVE_PRIME
      }
      end_segment_primes(ctx);
    }
  } END_SIMPLE_FOR_EACH_PRIME
  prime_count_cache_destroy(cache);
  return sum;
}

/* almost_prime_count(k,n) is the main interface, it will call the recursive
 * function _cs(), with the terminal function _final_sum(). */

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

  if (n >= max_nth_almost_prime(k))
    return max_almost_prime_count(k);

  if (k == 0) return n;
  if (k == 1) return prime_count(n);
  if (k == 2) return semiprime_count(n);
  if (k == 3) return almost3prime_count(n);
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

UV almost_prime_count_range(uint32_t k, UV lo, UV hi) {
  if (k == 0) return (lo <= 1 || hi >= 1);
  if (k == 1) return prime_count_range(lo, hi);
  if (k == 2) return semiprime_count_range(lo, hi);
  /* See semiprime_count.  Possibly clever solutions for small ranges. */
  if (k >= BITS_PER_WORD || (hi >> k) == 0 || hi < lo) return 0;
  return almost_prime_count(k, hi)
         - (((lo >> k) == 0) ? 0 : almost_prime_count(k,lo-1));
}

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

#if 0     /* Equation 6 from https://arxiv.org/pdf/2103.09866v3.pdf */
  {
    const double nu[21] = {
      1.0, 2.61497e-1, -5.62153e-1, 3.05978e-1, 2.62973e-2, -6.44501e-2,
      3.64064e-2, -4.70865e-3, -4.33984e-4, 1.50850e-3, -1.83548e-4,
      1.49365e-4, 4.99174e-5, 1.82657e-5, 1.30241e-5, 5.52779e-6,
      2.90194e-6, 1.45075e-6, 7.19861e-7, 3.61606e-7, 1.80517e-7 };
    double sum = 0, x = n, logx = log(x), loglogx = log(logx);
    uint32_t i, j;
    for (j = 0; j < k; j++) {
      uint32_t idx = k-1-j;
      double v = (idx <= 20) ? nu[idx] : 0.1893475 * powl(2.0, -(double)k);
      for (i = 1; i <= j; i++)
        v = v * loglogx / i;
      sum += v;
    }
    sum = (x / logx) * sum;
    return (UV) (sum+0.5);
  }
#endif

  /* We should look at (1) better bounds, (2) better weighting here */
  /* return lo + (hi-lo)/2; */
  /* Consider two variables to control our weight: k and n */
  if (k > 11) return lo + (hi-lo) * 0.20;
  return lo + (hi-lo) * 0.76;
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
  UV r, max, lo;

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

  return inverse_interpolate_k(lo, 0, n, k, &apcl, 0);
}

UV nth_almost_prime_lower(uint32_t k, UV n) {
  UV r, max, lo;

  if (n == 0) return 0;
  if (k == 0) return (n == 1) ? 1 : 0;
  if (k == 1) return nth_prime_lower(n);
  if (n < 8) return _fast_small_nth_almost_prime(k, n);

  max = max_almost_prime_count(k);
  if (n > max) return 0;
  r = reduce_nth_factor(k,n);
  if (r > 0)  return nth_almost_prime_lower(k-r, n) << r;

  /* We start out with the literal min and max because we have NO idea. */
  /*  \_/ note 3 instead of 5!  TODO:  apcu is not tight enough, so reduce */
  lo = 3 * (UVCONST(1) << k);  /* For k >= 1 and n >= 8 */

  return inverse_interpolate_k(lo, 0, n, k, &apcu, 0);
}

UV nth_almost_prime_approx(uint32_t k, UV n) {
  UV max, lo;

  if (n == 0) return 0;
  if (k == 0) return (n == 1) ? 1 : 0;
  if (k == 1) return nth_prime_approx(n);
  if (k == 2) return nth_semiprime_approx(n);

  max = max_almost_prime_count(k);
  if (n > max) return 0;
  /* We could reduce but really no reason to do it */

  if (n < 8) return _fast_small_nth_almost_prime(k,n);
  lo = nth_almost_prime_lower(k,n);

  return inverse_interpolate_k(lo, 0, n, k, &apca, 0);
}

static UV   _cb_nth3(UV n) { return nth_almost_prime_approx(3,n); }
static UV   _cb_cnt3(UV n) { return almost_prime_count(3,n); }
static bool _cb_is3(UV n)  { return is_almost_prime(3,n); }

static UV   _cb_nth4(UV n) { return nth_almost_prime_approx(4,n); }
static UV   _cb_cnt4(UV n) { return almost_prime_count(4,n); }
static bool _cb_is4(UV n)  { return is_almost_prime(4,n); }

UV nth_almost_prime(uint32_t k, UV n) {
  UV r, lo, hi;

  if (n == 0) return 0;
  if (k == 0) return (n == 1) ? 1 : 0;
  if (k == 1) return nth_prime(n);
  if (k == 2) return nth_semiprime(n);

  if (n > max_almost_prime_count(k)) return 0;

  /* For k >= 3 and small n we can answer this quickly. */
  if (n < 8) return _fast_small_nth_almost_prime(k,n);
  r = reduce_nth_factor(k,n);
  if (r > 0)  return nth_almost_prime(k-r,n) << r;
  /* NOTE: given n a 64-bit integer, k always <= 40 after reduction */

  /* Using the approximation to narrow in is *much* more efficient.  But
   * there is no good way to make it generic without closures (GCC extension)
   * or statics (not thread-safe). */
  if (k == 3)
    return interpolate_with_approx(n, 0, 20000, &_cb_nth3, &_cb_cnt3, &_cb_is3);
  if (k == 4)
    return interpolate_with_approx(n, 0, 20000, &_cb_nth4, &_cb_cnt4, &_cb_is4);

  lo = nth_almost_prime_lower(k,n);
  hi = nth_almost_prime_upper(k,n);
  hi = inverse_interpolate_k(lo, hi, n, k, &apce, 60000);
  while (!is_almost_prime(k,hi))
    hi--;
  return hi;
}


/******************************************************************************/
/*                                  Bounds                                    */
/******************************************************************************/

/* Bayless et al. (2018) and Kinlaw (2019) are main references.
 *
 * https://www.researchgate.net/publication/329788487_Sums_over_primitive_sets_with_a_fixed_number_of_prime_factors
 * http://math.colgate.edu/~integers/t22/t22.pdf
 * https://arxiv.org/pdf/2103.09866v3.pdf
 *
 * Note that they use Pi_k(x) to mean square-free numbers, and
 * Tau_k(x) to mean the general count like we use.
 * They also have results for k = 2,3,4 only.
 * Also see https://archimede.mat.ulaval.ca/MAINE-QUEBEC/19/Kinlaw19.pdf.
 *
 * We split into three ranges:
 *    1    - 2^20  complete computations
 *    2^20 - 2^32  complete computations
 *    2^32 - 2^64  correct upper for k=2,3,4.  correct lower for k=2.
 *                 empirical for other k.
 *
 */

static const double _upper_20[13] = {0,0, 1.006,0.7385,0.6830,0.5940,0.3596,0.2227,0.1439, 0.09785,0.07016,0.05303,0.04202};
static const double _upper_32[21] = {0,0, 1.013,0.8094,0.7485,
  /*  5-12 */ 0.6467,0.3984,0.2464,0.1572,0.1049,0.07364,0.05452,0.04266,
  /* 13-20 */ 0.03542,0.03082,0.02798,0.02642,0.02585,0.02615,0.02808,0.03054};
static const double _upper_64[41] = {0,0, 1.028, 1.028, 1.3043,/* <--corrrect */
  /* 5-12 */
  0.72208, 0.46609, 0.29340,0.18571,0.12063,0.0815,0.0575,0.0427,
  /* 13-20 */
  0.03490, 0.03007, 0.02710, 0.02554, 0.02504, 0.02554, 0.02699, 0.02954,
  /* 21-28 */
  0.03294, 0.03779, 0.04453, 0.05393, 0.06703, 0.08543, 0.1117, 0.1494,
  /* 29-31 */
  0.205,0.287,0.410,
  /* 32-40 */
  0.60,0.90,1.36,2.12,3.35,5.38,8.83,14.75,25.07,
};

static const double _lower_20[13] = {0,0, 0.8197, 0.8418, 0.5242,
  /*  5-12 */ 0.5154,0.3053,0.1901,0.1253,0.0892,0.06551,0.05082,0.04101};
static const double _lower_32[21] = {0,0, 1.004,  0.7383, 0.6828,
  /*  5-12 */ 0.5939,0.3594,0.2222,0.1438,0.09754,0.06981,0.05245,0.04151,
  /* 13-20 */ 0.03461,0.03006,0.02709,0.02553,0.02502,0.02552,0.02697,0.02945 };
static const double _lower_64[41] = {0,0, 1.011,  0.8093, 0.7484,
  /* 5-12 */
  0.6465,0.3982,0.2463,0.1571,0.1048,0.07363,0.0545,0.0422,
  /* 13-20 */
  0.0331,0.0270,0.0232,0.0208,0.0194,0.0190,0.0193,0.0203,
  /* 21-28 */
  0.0222,0.0252,0.0295,0.0356,0.0444,0.0570,0.0753,0.102,
  /* 29-31 */
  0.14,0.20,0.297,
  /* 32-40 */
  0.44,0.68,1.07,1.71,2.8,4.7,8.0,13.89,23.98,
};

/*
   k,count n <= 2^64-1
   1,425656284035217743
   2,1701748900850019777                 10 hours
   3,3167597434038354478                320 hours
   4,3787884015050788482                322 hours
   5,3378907169603895030                294 hours
   6,2466706950238087748                209 hours
   7,1571012171387856192                123 hours
   8,913164427599983727                  82 hours
   9,499840874923678341                  42 hours
  10,263157990621533964                  20 hours
  11,135128109904869290                  12 hours
  12,68283616225825256                    7 hours
  13,34151861008771016                    4 hours
  14,16967424859951587                    2 hours
  15,8393048221327186
  16,4139595949113890
  17,2037655246635364
  18,1001591348315641
  19,491808604962296
  20,241293656953012
  21,118304122014405
  22,57968649799947
  23,28388662714236
  24,13895161400556
  25,6797526392535
  26,3323560145881
  27,1624109166018
  28,793189260998
  29,387148515886
  30,188844769357
  31,92054377509
  32,44841620426
  33,21827124353
  34,10616326552
  35,5159281045
  36,2505087309
  37,1215204383
  38,588891145
  39,285076316
  40,137840686
*/

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
  if (n >= max_nth_almost_prime(k))
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
    /* Bayless (2018) Theorem 5.3 proves that multu=1.028 is a correct bound
     * for all x >= 10^12.  However it is not a tight bound for the range
     * 2^32 to 2^64.  We tighten it a lot for the reduced range.
     */
    if (n > 4294967295U)  multu = 0.8711;
  } else if (k == 4) {
    /* Bayless Theorem 5.4 part 1 (with multu = 1.3043) */
    boundl = boundu = x * logplus*logplus*logplus / (6*logx);
    /* Bayless Theorem 5.4 part 2 */
    if (x > 1e12) {
      boundu += 0.511977 * x * (log(log(x/4)) + 0.261536) / logx;
      multu = 1.028;
    }
    /* As with k=3, adjust to tighten in the finite range. */
    if (n > 4294967295U)  multu = 0.780;
    if (x > 1e12)         multu = 0.6921;
  } else {
    /* Completely empirical and by no means optimal.
     * It is easy and seems fairly reasonable through k=20 or so.
     *
     * For high k, this follows the lower bound well but upper grows too fast.
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
  if (k == 2 && n < 256) return semiprime_count(n);
  _almost_prime_count_bounds(&l, &u, k, n);
  return u;
}

UV almost_prime_count_lower(uint32_t k, UV n) {
  UV l, u;
  if (k == 2 && n < 256) return semiprime_count(n);
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
  static const UV max[32] = {1,203280221,658662065,967785236,916899721,665533848,410630859,229679168,121092503,61600699,30653019,15043269,7315315,3535071,1700690,814699,389357,185245,87964,41599,19611,9184,4283,2001,914,421,187,84,37,15,7,2};
#else
  static const UV max[64] = {1,
    UVCONST( 425656284035217743), /* max prime count */
    UVCONST(1701748900850019777), /* max semiprime count */
    UVCONST(3167597434038354478), /* max 3-almost-prime count */
    UVCONST(3787884015050788482), /* max 4-almost-prime count */
    /* 5-12 */  UVCONST(3378907169603895030),UVCONST(2466706950238087748),UVCONST(1571012171387856192),UVCONST(913164427599983727),UVCONST(499840874923678341),UVCONST(263157990621533964),UVCONST(135128109904869290),UVCONST(68283616225825256),
    /* 13-22 */ UVCONST(34151861008771016),UVCONST(16967424859951587),UVCONST(8393048221327186),UVCONST(4139595949113890),UVCONST(2037655246635364),UVCONST(1001591348315641),UVCONST(491808604962296),UVCONST(241293656953012),UVCONST(118304122014405),UVCONST(57968649799947),
    /* 23-32 */ UVCONST(28388662714236),UVCONST(13895161400556),UVCONST(6797526392535),UVCONST(3323560145881),UVCONST(1624109166018),UVCONST(793189260998),UVCONST(387148515886),UVCONST(188844769357),UVCONST(92054377509),UVCONST(44841620426),
    /* 33-63 */ UVCONST(21827124353),UVCONST(10616326552),UVCONST(5159281045),UVCONST(2505087309),1215204383,588891145,285076316,137840686,66567488,32103728,15460810,7433670,3567978,1709640,817053,389954,185387,87993,41604,19611,9184,4283,2001,914,421,187,84,37,15,7,2
  };
#endif
  if (k >= BITS_PER_WORD) return 0;
  /* if (max[k] == 0) return UV_MAX;  All filled in so no need */
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
    sort_uv_array(L, *count);
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
  if (Ssize > 10000000UL) Ssize = 10000000UL;
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

/* Algorithm from Trizen, May 2022 */
static void _genkap(UV lo, UV hi, uint32_t k, UV m, UV begp, UV **List, UV *Lpos, UV *Lsize) {
  if (k == 1) {

    UV pos = *Lpos,  size = *Lsize, *L = *List;
    UV start = lo/m + (lo % m != 0),  endp = hi/m;

    if (start > begp) begp = start;

    if (endp < 10000000U) {
      START_DO_FOR_EACH_PRIME(begp, endp) {
        if (L != 0) {
          if (pos >= size)  Renew(L, size += 100000, UV);
          L[pos] = m*p;
        }
        pos++;
      } END_DO_FOR_EACH_PRIME
    } else {
      UV i, count, *list;
      count = range_prime_sieve(&list, begp, endp);
      if (L == 0) {
        pos += count;
      } else {
        if ((pos + count - 1) >= size)  Renew(L, size += (count + 100000), UV);
        for (i = 0; i < count; i++)
          L[pos++] = m * list[i];
      }
      Safefree(list);
    }
    *Lpos = pos;
    *Lsize = size;
    *List = L;

  } else {

    UV p, s;
    for (s = rootint(hi/m, k), p = begp;  p <= s;  p = next_prime(p)) {
      UV t = m * p;
      if ((lo/t + (lo % t != 0)) <= (hi/t))
        _genkap(lo, hi, k-1, t, p, List, Lpos, Lsize);
    }

  }
}

UV generate_almost_primes(UV** list, uint32_t k, UV lo, UV hi) {
  UV *L, Lpos = 0, Lsize, countest;

  if (k == 0 || k >= BITS_PER_WORD) { *list = 0; return 0; }
  if ((lo >> k) == 0) lo = UVCONST(1) << k;
  if (hi > max_nth_almost_prime(k)) hi = max_nth_almost_prime(k);
  if (lo > hi) { *list = 0; return 0; }

  /* For these small k values, these are typically faster */
  if (k == 1) return range_prime_sieve(list, lo, hi);
  if (k == 2) return range_semiprime_sieve(list, lo, hi);

  /* Large base with small range: better to sieve */
  if ( (k >= 3 && hi >= 1e12 && (hi-lo) <= 5e6) ||
       (k >= 3 && hi >= 1e13 && (hi-lo) <= 2e8) ||
       (k >= 3 && hi >= 1e14 && (hi-lo) <= 4e8)    )
    return range_almost_prime_sieve(list, k, lo, hi);

  /* Optional:  we could try reduce_k_for_n() here. */

  prime_precalc(10000000U);
  countest = almost_prime_count_approx(k,hi) - almost_prime_count_approx(k,lo-1);
  Lsize = (countest > 10000000U) ? 10000000U : countest+1000;

  New(0, L, Lsize, UV);
  _genkap(lo, hi, k, 1, 2, &L, &Lpos, &Lsize);
  sort_uv_array(L, Lpos);
  *list = L;
  return Lpos;
}


/******************************************************************************/
/*                                CHEN PRIMES                                 */
/******************************************************************************/

/* consider Chen(h,k) where p prime and bigomega(p+h) <= k */

#if BITS_PER_WORD == 64
  #define MAX_CHEN_PRIME UVCONST(18446744073709551437)
#else
  #define MAX_CHEN_PRIME UVCONST(4294967291)
#endif

bool is_chen_prime(UV n) {
  if (n < 2 || n > MAX_CHEN_PRIME) return 0;
  return (is_prime(n) && (is_prime(n+2) || is_semiprime(n+2)));
}

UV next_chen_prime(UV n) {
  for ( n = next_prime(n);  n != 0 && n < MAX_CHEN_PRIME;  n = next_prime(n+2) )
    if (is_prime(n+2) || is_semiprime(n+2))
      return n;
  return 0;
}
