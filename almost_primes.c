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
#include "semi_primes.h"
#include "almost_primes.h"

#define HCPC 101
static const uint32_t _pcnt[1 + ((HCPC-1)>>1)] = {
1,2,3,4,4,5,6,6,7,8,8,9,9,9,10,11,11,11,12,12,13,14,14,15,15,15,16,16,16,17,18,18,18,19,19,20,21,21,21,22,22,23,23,23,24,24,24,24,25,25,26
};

/******************************************************************************/
/*                              ALMOST PRIMES                                 */
/******************************************************************************/

static UV _next_almost_prime(UV n, UV k) {
  while (!is_almost_prime(++n, k))
    ;
  return n;
}
static UV _prev_almost_semiprime(UV n, UV k) {
  while (!is_almost_prime(--n, k))
    ;
  return n;
}

typedef struct {
  uint32_t *count;
  uint32_t lastidx;
} prime_array_t;

static UV _cached_count(UV n, prime_array_t precalc) {
  if (n < 2) return 0;
  return (n <= precalc.lastidx)  ?  precalc.count[(n-1)>>1]
                                 :  LMO_prime_count(n);
}
static void _create_precalc(UV limit, prime_array_t *precalc) {
  uint32_t *counts = 0, idx = 1, cnt = 1;

  if (limit <= HCPC) {
    limit = HCPC;
    counts = (uint32_t*) _pcnt;
    if (_pcnt[(limit-1)>>1] != LMO_prime_count(limit)) croak("fail limit");
  } else {
    limit |= 1;
    New(0, counts, limit+1, uint32_t);
    counts[0] = 1;
    START_DO_FOR_EACH_PRIME(3, limit) {
      while (idx < ((p-1)>>1)) counts[idx++] = cnt;
      counts[idx++] = ++cnt;
    } END_DO_FOR_EACH_PRIME
    while (idx <= ((limit-1)>>1)) counts[idx++] = cnt;
  }
  precalc->lastidx = limit;
  precalc->count = counts;
}
static void _destroy_precalc(prime_array_t *precalc) {
  if (precalc->count != 0 && precalc->count != _pcnt)
    Safefree(precalc->count);
}

/* for Pi from Pi to isqrt(N/Pi) [pc[n/Pi]-idx(Pi)+1] */
/* semiprime count = _final_sum(n, 1, 2, cache); */
/* 3-almost prime count = sum(Pj < icbrt(n) of _final_sum(n, Pj, Pj, cache); */
static UV _final_sum(UV n, UV pdiv, UV lo, prime_array_t cache) {
  UV s = 0, hi = isqrt(n/pdiv);
  UV j = _cached_count(lo, cache) - 1;  /* IDX(Pi) */

  if (hi-lo < 500) {
    SIMPLE_FOR_EACH_PRIME(lo, hi) {
      s += _cached_count(n/(pdiv*p),cache) - j++;
    } END_SIMPLE_FOR_EACH_PRIME
    return s;
  }

  if ((lo <= 2) && (hi >= 2)) s += _cached_count(n/(pdiv*2),cache) - j++;
  if ((lo <= 3) && (hi >= 3)) s += _cached_count(n/(pdiv*3),cache) - j++;
  if ((lo <= 5) && (hi >= 5)) s += _cached_count(n/(pdiv*5),cache) - j++;
  if (lo < 7) lo = 7;
  if (lo <= hi) {
    unsigned char* segment;
    UV seg_base, seg_low, seg_high;
    void* ctx = start_segment_primes(lo, hi, &segment);
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
        s += _cached_count(n/(pdiv*p),cache) - j++;
      END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
  }
  return s;
}

static UV _cs(UV n, UV pdiv, UV lo, UV k, prime_array_t cache) {
  UV count = 0;

  if (k == 2)
    return _final_sum(n, pdiv, lo, cache);

  SIMPLE_FOR_EACH_PRIME(lo, rootof(n/pdiv,k)) {
    if (k == 3)  count += _final_sum(n, pdiv*p, p, cache);
    else         count += _cs(n, pdiv*p, p, k-1, cache);
  } END_SIMPLE_FOR_EACH_PRIME

  return count;
}

UV almost_prime_count(UV n, UV k)
{
  prime_array_t cache;
  UV count, csize;

  if (k == 0) return (n >= 1);
  if (k == 1) return LMO_prime_count(n);
  if (k == 2) return semiprime_count(0,n);

  if ((n >> k) == 0) return 0;

  csize = isqrt(n);
  if (csize > 128*1024*1024) csize = 128*1024*1024;
  _create_precalc(csize, &cache);
  count = _cs(n, 1, 2, k, cache);
  _destroy_precalc(&cache);
  return count;
}

#if 0
UV almost_prime_count(UV lo, UV hi, UV k) {
  if (k == 0) return (lo <= 1 || hi >= 1);
  if (k == 1) return prime_count(lo, hi);
  if (k == 2) return semiprime_count(lo, hi);
  return _almost_prime_count(hi, k) - _almost_prime_count(lo-1, k);
}
#endif

/* Bayless et al. (2018) give very nice results.
 * Note that they use Pi_k(x) to mean square-free numbers, and
 * Tau_k(x) to mean the general count like we use. */
UV almost_prime_count_upper(UV n, UV k) {
  double x, bound, logx, loglogx;

  if (k == 0) return (n >= 1);
  if (k == 1) return prime_count_lower(n);
  if ((n >> k) == 0) return 0;

  x = (double) n;
  logx = log(x);
  loglogx = log(logx);

  if (k == 2) {
    bound = 1.028 * x * (loglogx + 0.261536) / logx;
  } else if (k == 3) {
    /* See http://math.colgate.edu/~integers/t22/t22.pdf */
    /* Bound for x >= 10^12 */
    bound = 1.028 * x * ((loglogx+0.26153)*(loglogx+0.26153) + 1.055852) / (2*logx);
    /* Verified ok for x < 1e9, which leaves 1e9 to 1e12  */
    /* TODO: For x < 10^12, ....? */
  } else if (k == 4) {
    bound = 1.3043 * x * (pow(loglogx+0.26153,3)) / (6*logx);
    /* TODO: For x >= 10^12 can be improved */
  } else {
    /* TODO ???? */
  }
  if (bound >= UV_MAX) croak("bound overflow\n");
  return (UV) (bound + 1e-8);
}

UV almost_prime_count_approx(UV n, UV k) {
  double loglogn, est;
  UV i;

  if (k == 0) return (n >= 1);
  if (k == 1) return prime_count_approx(n);
  if (k == 2) return semiprime_count_approx(n);
  if ((n >> k) == 0) return 0;

  /* Landau:  Pi_k(x) ~ x/log(x) * (loglogx)^(k-1) / (k-1)!
   * This is really bad for "small" x values however.
   * Selberg shows a complicated improvement. */

  loglogn = log(log(n));

  est = n * log(n);
  for (i = 1; i < k; i++)
    est *= loglogn/(double)i;

  return (UV) (est+0.5);
}

UV nth_almost_prime_approx(UV n, UV k) {
  UV lo, hi;

  if (n == 0) return 0;

  /* TODO: Overflow is a concern */
  /* I have no idea what an actual upper bound is. Even the asymptotic
   * has an enormous error term.  So we'll guess rather low, then test
   * it.  This isn't so bad as if it is low it isn't expensive, and gets
   * us a better lower bound. */

  lo = 1 << k;
  hi = lo + nth_prime_upper(n);
  if (k > 5) hi *= 2;
  if (k > 10) hi *= 2;

  while (almost_prime_count_approx(hi, k) <= n) {
    lo = hi+1;
    hi = hi*2;
  }

  while (lo < hi) {
    UV mid = lo + (hi-lo)/2;
    if (almost_prime_count_approx(mid, k) <= n) lo = mid+1;
    else                                 hi = mid;
  }
  return lo-1;
}

UV nth_almost_prime(UV n, UV k) {
  UV lo, hi, apc;

  if (n == 0) return 0;
  if (k == 0) return (n == 1) ? 1 : 0;
  if (k == 1) return nth_prime(n);
  if (k == 2) return nth_semiprime(n);

  /* TODO: Overflow is a big concern */
  /* I have no idea what an actual upper bound is. Even the asymptotic
   * has an enormous error term.  So we'll guess rather low, then test
   * it.  This isn't so bad as if it is low it isn't expensive, and gets
   * us a better lower bound. */

  /* TODO: Idea -- use the approx count to get an initial value. */
  /* TODO: See nth_ramanujan_prime_approx for interpolation. */
  /* TODO: See nth_semiprime for another way, including sieving. */

  lo = 1 << k;
  hi = nth_almost_prime_approx(n,k);  /* This is a *terrible* estimate */

  apc = almost_prime_count(hi,k);
  while (apc <= n) {
    lo = hi+1;
    /* hi->apc.  We got one kap for each (hi/apc).  We need n-apc more. */
    double f = (double) hi / apc;
    hi = hi + 1.02 * f * (n-apc) + 100;
    apc = almost_prime_count(hi,k);
  }

  while (lo < hi) {
    UV mid = lo + (hi-lo)/2;
    if (almost_prime_count(mid, k) <= n) lo = mid+1;
    else                                 hi = mid;
  }
  /* So ... this is really expensive for big k. */
  /* We need the exact point where we go from count = n-1 to count = n */
  return _prev_almost_semiprime(lo-1+1, k);
}

UV max_nth_almost_prime(UV k) {
#if BITS_PER_WORD == 32
  static const UV offset[32] = {0,4,1,9,5,0,7,47,31,3,15,511,1263,5119,1023,255,23295,2559,4095,126975,16383,262143,2359295,2097151,5767167,1048575,33554431,16777215,100663295,67108863,268435455,1073741823};
#else
  static const UV offset[64] = {0,58,14,2,4,3,17,0,1,195,51,127,63,767,1535,511,255,15,8191,1023,83967,16383,111615,557055,2097151,524287,65535,1048575,6553599,33554431,
    4194303, 671088639, 16777215, 536870911,
    UVCONST(2684354559),UVCONST(2147483647),
    UVCONST(25769803775),UVCONST(4294967295),
    UVCONST(268435455),UVCONST(206158430207),
    UVCONST(137438953471),UVCONST(17179869183),
    UVCONST(68719476735),UVCONST(2199023255551),
    UVCONST(5428838662143),UVCONST(21990232555519),
    UVCONST(4398046511103),UVCONST(1099511627775),
    UVCONST(100055558127615),UVCONST(10995116277759),
    UVCONST(17592186044415),UVCONST(545357767376895),
    UVCONST(70368744177663),UVCONST(1125899906842623),
    UVCONST(10133099161583615), UVCONST(9007199254740991),
    UVCONST(24769797950537727),UVCONST(4503599627370495),
    UVCONST(144115188075855871),UVCONST(72057594037927935),
    UVCONST(432345564227567615),UVCONST(288230376151711743),
    UVCONST(1152921504606846975),UVCONST(4611686018427387903)
  };
#endif
  if (k >= BITS_PER_WORD) return 0;
  return UV_MAX - offset[k];
}

UV max_almost_prime_count(UV k) {
#if BITS_PER_WORD == 32
  static const UV max[32] = {1,203280221,658662065,967785236,916899721,665533848,410630859,229679168,121092503,61600699,30653019,15043269,7315315,3535071,1700690,814699,389357,185245,87964,41599,19611,9184,4283,2001,914,421,187,84,37,15,7,2,0};
#else
  static const UV max[64] = {1,
    UVCONST(425656284035217743),  /* max prime count */
    UVCONST(1701500000000000000), /* max semiprime count (estimate) */
    UVCONST(3738241841329689241), /* max 3-almost-prime count (upper limit) */
    UVCONST(5626791887645067628), /* max 4-almost-prime count (upper limit) */
    /* 5-12 */  0,0,0,0,0,0,0,0,
    /* 13-22 */ 0,0,UVCONST(8393048221327186),UVCONST(4139595949113890),UVCONST(2037655246635364),UVCONST(1001591348315641),UVCONST(491808604962296),UVCONST(241293656953012),UVCONST(118304122014405),UVCONST(57968649799947),
    /* 23-32 */ UVCONST(28388662714236),UVCONST(13895161400556),UVCONST(6797526392535),UVCONST(332356014588),UVCONST(1624109166018),UVCONST(793189260998),UVCONST(387148515886),UVCONST(188844769357),UVCONST(92054377509),UVCONST(44841620426),
    /* 33-63 */ UVCONST(21827124353),UVCONST(10616326552),UVCONST(5159281045),UVCONST(2505087309),1215204383,588891145,285076316,137840686,66567488,32103728,15460810,7433670,3567978,1709640,817053,389954,185387,87993,41604,19611,9184,4283,2001,914,421,187,84,37,15,7,2
  };
#endif
  if (k >= BITS_PER_WORD) return 0;
  if (max[k] == 0) return UV_MAX;  /* TODO: find these */
  return max[k];
}
