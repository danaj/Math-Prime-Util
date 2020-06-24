#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ptypes.h"
#include "constants.h"
#define FUNC_isqrt 1
#define FUNC_icbrt 1
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
  UV j = _cached_count((lo < 2) ? 2 : lo, cache) - 1;  /* IDX(Pi) */

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
  if (csize > 32*1024*1024) csize = 32*1024*1024;
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


UV almost_prime_count_approx(UV n, UV k) {
  double loglogn, pc, est;
  UV i;

  if (k == 0) return (n >= 1);
  if (k == 1) return prime_count_approx(n);
  if (k == 2) return semiprime_count_approx(n);
  if ((n >> k) == 0) return 0;

  pc = est = prime_count_approx(n);
  loglogn = log(log(n));

  for (i = 1; i < k; i++)
    est *= loglogn/i;

  if      (k == 3) est -= .050 * pc * loglogn * loglogn;
  else if (k == 4) est -= .025  * pc * loglogn * loglogn * loglogn;

  return (UV) (est+0.5);
}
