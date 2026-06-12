/***********    inverse_sigma0_list   and   inverse_sigma0_count    ***********/

#include <string.h>
#define FUNC_isqrt 1
#define FUNC_is_perfect_square 1
#define FUNC_log2floor 1
#include "ptypes.h"
#include "cache.h"
#include "constants.h"
#include "factor.h"
#include "inverse_sigma0.h"
#include "prime_counts.h"
#include "primality.h"
#include "sieve.h"
#include "sort.h"
#include "util.h"

/*
 * Given a range [lo,hi] and a k value, returns or counts all the values
 * in the range that have a divisor count equal to k.
 *
 * We could loop over the range computing scalar(divisors(n)), equivalently
 * divisor_sum(n,0).  This is fine for a few values but not very fast for
 * larger ranges.  The MPU level 'forfactored { }' iterator helps.
 *
 * 1) SIEVE
 *    A lightweight segmented sieve.  It tracks only the divisor count and
 *    residual value for each integer in the segment, not the full factor list.
 *    When sieving to sqrt(hi) is cheap, residual values are prime.  Otherwise
 *    we sieve through cuberoot(hi), after which each residual has at most two
 *    prime factors, so classifying it as 1, prime, prime square, or semiprime
 *    is enough.
 *
 * 2) GENERATE
 *    Here rather than sieving the range, we generate values in the range
 *    that have the desired divisor count.  Can be very efficient for large
 *    ranges, but can also spend too long for some narrow range / k
 *    combinations.
 */

#ifndef INVERSE_SIGMA0_TAU_MAX_WIDTH
#define INVERSE_SIGMA0_TAU_MAX_WIDTH UVCONST(10000000)
#endif

#ifndef INVERSE_SIGMA0_TAU_MAX_CHUNKED_WIDTH
#define INVERSE_SIGMA0_TAU_MAX_CHUNKED_WIDTH UVCONST(100000000)
#endif

#if BITS_PER_WORD == 64
#ifndef INVERSE_SIGMA0_TAU_CHUNK_MIN_LO
#define INVERSE_SIGMA0_TAU_CHUNK_MIN_LO UVCONST(100000000000)
#endif
#endif

static bool full_factor_sieve(UV sqrtn, UV range)  /* range = hi-lo */
{
  if (sqrtn <    1000000U) return (range > 200000);  /* Below 10^12 */
  if (sqrtn <   10000000U) return (range > 300000);  /* Below 10^14 */
  if (sqrtn <   35000000U) return (range > 350000);  /* Below 10^15 */
  if (sqrtn <  100000000U) return (range > 500000);  /* Below 10^16 */
  if (sqrtn <  350000000U) return (range > 750000);  /* Below 10^17 */
  if (sqrtn < 1000000000U) return (range > 1000000); /* Below 10^18 */
  if (sqrtn < 3500000000U) return (range > 1500000); /* Below 10^19 */
  return (range > 2000000);
}

static int first_multiple_in_range(UV *A, UV lo, UV hi, UV p)
{
  UV a;

  a = (lo / p) * p;
  if (a < lo) {
    if (a > UV_MAX - p) return 0;
    a += p;
  }
  if (a > hi) return 0;

  *A = a;
  return 1;
}

static void tau_update(uint16_t *tau, UV k, UV mult)
{
  UV t;

  if (*tau == 0) return;
  t = *tau;
  if (mult == 0 || t > k / mult) {
    *tau = 0;
  } else {
    t *= mult;
    if (k % t != 0)
      *tau = 0;
    else
      *tau = (uint16_t)t;
  }
}

static void tau_update_residual(uint16_t *tau, UV k, UV rem, bool full_sieve)
{
  UV need;

  if (*tau == 0 || rem <= 1) return;
  if (k % *tau != 0) { *tau = 0; return; }
  need = k / *tau;

  if (full_sieve) {
    *tau = (need == 2) ? (uint16_t)k : 0;
  } else if (need == 2) {
    *tau = is_def_prime(rem) ? (uint16_t)k : 0;
  } else if (need == 3) {
    *tau = is_perfect_square(rem) ? (uint16_t)k : 0;
  } else if (need == 4) {
    *tau = (!is_perfect_square(rem) && !is_def_prime(rem)) ? (uint16_t)k : 0;
  } else {
    *tau = 0;
  }
}

static UV* tau_sieve(UV *count, UV lo, UV hi, UV k, bool count_only)
{
  UV range, sqrthi, sievelim, j, len, alloc, *rem, *list;
  uint16_t *tau;
  bool full_sieve;

  range = hi - lo + 1;
  sqrthi = isqrt(hi);
  full_sieve = full_factor_sieve(sqrthi, hi-lo);
  sievelim = full_sieve ? sqrthi : icbrt(hi);

  len = alloc = 0;
  list = 0;
  New(0, rem, range, UV);
  New(0, tau, range, uint16_t);
  for (j = 0; j < range; j++) {
    rem[j] = lo + j;
    tau[j] = 1;
  }

  START_DO_FOR_EACH_PRIME(2, sievelim) {
    UV A, idx;
    if (first_multiple_in_range(&A, lo, hi, p)) {
      for (idx = A - lo; idx < range; idx += p) {
        UV e = 0;
        if (tau[idx] == 0) continue;
        while ((rem[idx] % p) == 0) {
          e++;
          rem[idx] /= p;
        }
        if (e > 0)
          tau_update(&tau[idx], k, e+1);
      }
    }
  } END_DO_FOR_EACH_PRIME

  for (j = 0; j < range; j++) {
    if (tau[j] != 0 && rem[j] > 1)
      tau_update_residual(&tau[j], k, rem[j], full_sieve);
    if (tau[j] == k) {
      if (count_only) {
        len++;
      } else {
        if (len >= alloc) {
          if (alloc == 0) {
            alloc = 32;
            New(0, list, alloc, UV);
          } else {
            alloc *= 2;
            Renew(list, alloc, UV);
          }
        }
        list[len++] = lo + j;
      }
    }
  }

  Safefree(tau);
  Safefree(rem);
  *count = len;

  return list;
}

static UV* tau_sieve_chunked(UV *count, UV lo, UV hi, UV k, bool count_only)
{
  UV len, alloc, seglo, *list;

  len = alloc = 0;
  list = 0;
  seglo = lo;

  while (seglo <= hi) {
    UV segcount, seghi, *seglist;

    seghi = (hi - seglo + 1 > INVERSE_SIGMA0_TAU_MAX_WIDTH)
          ? seglo + INVERSE_SIGMA0_TAU_MAX_WIDTH - 1
          : hi;
    seglist = tau_sieve(&segcount, seglo, seghi, k, count_only);

    if (count_only) {
      len += segcount;
    } else if (segcount > 0) {
      if (len + segcount > alloc) {
        alloc = (alloc == 0) ? 32 : alloc;
        while (len + segcount > alloc)
          alloc *= 2;
        if (list == 0) New(0, list, alloc, UV);
        else           Renew(list, alloc, UV);
      }
      Copy(seglist, list+len, segcount, UV);
      len += segcount;
    }

    Safefree(seglist);
    if (seghi == hi) break;
    seglo = seghi + 1;
  }

  *count = len;
  return list;
}

static bool use_tau_sieve(UV lo, UV hi, UV k)
{
  UV width = hi - lo + 1;

  if (k > 65535) return 0;
  if (width <= INVERSE_SIGMA0_TAU_MAX_WIDTH) return 1;
#if BITS_PER_WORD == 64
  if (lo >= INVERSE_SIGMA0_TAU_CHUNK_MIN_LO &&
      width <= INVERSE_SIGMA0_TAU_MAX_CHUNKED_WIDTH)
    return 1;
#endif
  return 0;
}

typedef struct {
  UV lo;
  UV hi;
  UV count;
  UV alloc;
  UV *list;
  UV *divs;
  UV ndivs;
  int nexps;
  uint8_t exps[MPU_MAX_FACTORS];
  bool count_only;
} invsig0_ctx_t;

static int uv_in_list(UV v, const UV *list, int len)
{
  int i;
  for (i = 0; i < len; i++)
    if (list[i] == v)
      return 1;
  return 0;
}

static UV ceil_div_uv(UV n, UV d)
{
  return n / d + ((n % d) != 0);
}

static UV ceil_root_uv(UV n, uint32_t k)
{
  UV r;

  if (n <= 1 || k <= 1) return n;
  r = rootint(n, k);
  return r + (ipowsafe(r, k) < n);
}

static UV mul_pow_limit(UV n, UV p, uint8_t e, UV limit)
{
  while (e-- > 0) {
    if (p == 0 || n > limit / p) return 0;
    n *= p;
  }
  return n;
}

static void record(invsig0_ctx_t *ctx, UV n)
{
  if (n < ctx->lo || n > ctx->hi) return;

  if (ctx->count_only) {
    ctx->count++;
  } else {
    if (ctx->count >= ctx->alloc) {
      if (ctx->alloc == 0) {
        ctx->alloc = 32;
        New(0, ctx->list, ctx->alloc, UV);
      } else {
        ctx->alloc *= 2;
        Renew(ctx->list, ctx->alloc, UV);
      }
    }
    ctx->list[ctx->count++] = n;
  }
}

static int min_signature_le_hi(const uint8_t *exps, int nexps, UV hi)
{
  UV n = 1, p = 0;
  int i;

  for (i = 0; i < nexps; i++) {
    p = (i+1 < NPRIMES_SMALL) ? primes_small[i+1] : next_prime(p);
    n = mul_pow_limit(n, p, exps[i], hi);
    if (n == 0) return 0;
  }
  return 1;
}

static UV prime_upper_bound(invsig0_ctx_t *ctx, int idx, UV maxpow)
{
  UV denom = 1, base;
  uint32_t root, lower_sum = 0;
  int j, same;
  uint8_t e = ctx->exps[idx];

  /* Include the current run of equal exponents in the bound.  Without this,
   * signatures like [1,1,1] try far too many leading primes. */
  same = 1;
  while (idx+same < ctx->nexps && ctx->exps[idx+same] == e)
    same++;

  for (j = idx+same; j < ctx->nexps; j++)
    lower_sum += ctx->exps[j];

  if (lower_sum >= BITS_PER_WORD)
    return 0;
  denom = UVCONST(1) << lower_sum;
  base = maxpow / denom;
  if (base == 0)
    return 0;

  root = (uint32_t)e * (uint32_t)same;
  return rootint(base, root);
}

static void recurse(invsig0_ctx_t *ctx, int idx, UV n, UV *used)
{
  UV minpow, maxpow, pl, ph, *P, pcnt, i;
  uint8_t e = ctx->exps[idx];

  if (idx == ctx->nexps-1) {
    minpow = (n >= ctx->lo) ? 1 : ceil_div_uv(ctx->lo, n);
    maxpow = ctx->hi / n;
    pl = (e == 1) ? minpow : ceil_root_uv(minpow, e);
    ph = (e == 1) ? maxpow : rootint(maxpow, e);
    if (pl < 2) pl = 2;
    if (idx > 0 && e == ctx->exps[idx-1] && pl <= used[idx-1])
      pl = used[idx-1] + 1;
    if (pl > ph) return;

    if (ctx->count_only) {
      pcnt = prime_count_range(pl, ph);
      for (i = 0; i < (UV)idx; i++)
        if (used[i] >= pl && used[i] <= ph)
          pcnt--;
      ctx->count += pcnt;
      return;
    }

    pcnt = range_prime_sieve(&P, pl, ph);
    for (i = 0; i < pcnt; i++) {
      UV p = P[i];
      if (!uv_in_list(p, used, idx)) {
        UV v = (e == 1) ? n * p : mul_pow_limit(n, p, e, ctx->hi);
        if (v != 0) record(ctx, v);
      }
    }
    Safefree(P);
    return;
  }

  maxpow = ctx->hi / n;
  pl = 2;
  ph = rootint(maxpow, e);
  {
    UV sigph = prime_upper_bound(ctx, idx, maxpow);
    if (sigph < ph) ph = sigph;
  }
  if (idx > 0 && e == ctx->exps[idx-1] && pl <= used[idx-1])
    pl = used[idx-1] + 1;
  if (pl > ph) return;

  pcnt = range_prime_sieve(&P, pl, ph);
  for (i = 0; i < pcnt; i++) {
    UV p = P[i];
    if (!uv_in_list(p, used, idx)) {
      UV v = mul_pow_limit(n, p, e, ctx->hi);
      if (v != 0) {
        used[idx] = p;
        recurse(ctx, idx+1, v, used);
      }
    }
  }
  Safefree(P);
}

static void process_signature(invsig0_ctx_t *ctx, UV *parts, int len)
{
  UV used[MPU_MAX_FACTORS];
  UV maxexp = log2floor(ctx->hi);
  int i;

  if (len > MPU_MAX_FACTORS) return;
  for (i = 0; i < len; i++) {
    UV e = parts[len-1-i] - 1;
    if (e == 0 || e > maxexp) return;
    ctx->exps[i] = (uint8_t)e;
  }
  if (!min_signature_le_hi(ctx->exps, len, ctx->hi))
    return;

  ctx->nexps = len;
  recurse(ctx, 0, 1, used);
}

static void partitions(invsig0_ctx_t *ctx, UV rem, UV min_d, UV *parts, int len)
{
  UV i;

  if (len >= MPU_MAX_FACTORS) return;
  for (i = 0; i < ctx->ndivs; i++) {
    UV d = ctx->divs[i];
    if (d < 2 || d < min_d) continue;
    if (d > rem / d) break;
    if (rem % d == 0) {
      parts[len] = d;
      partitions(ctx, rem/d, d, parts, len+1);
    }
  }
  if (rem >= 2 && rem >= min_d) {
    parts[len] = rem;
    process_signature(ctx, parts, len+1);
  }
}

static void by_signature(invsig0_ctx_t *ctx, UV k)
{
  UV parts[MPU_MAX_FACTORS];

  ctx->divs = divisor_list(k, &ctx->ndivs, k);
  partitions(ctx, k, 2, parts, 0);
  Safefree(ctx->divs);
}


UV inverse_sigma0_count(UV lo, UV hi, UV k)
{
  invsig0_ctx_t ctx;

  if (k == 0 || hi < 1 || k > hi) return 0;
  if (k == 1) return (lo <= 1 && hi >= 1) ? 1 : 0;
  if (k == 2) return prime_count_range(lo, hi);
  if (lo < 1) lo = 1;
  if (lo > hi) return 0;
  if (use_tau_sieve(lo, hi, k)) {
    UV count;
    tau_sieve_chunked(&count, lo, hi, k, 1);
    return count;
  }

  memset(&ctx, 0, sizeof(ctx));
  ctx.lo = lo;
  ctx.hi = hi;
  ctx.count_only = 1;
  by_signature(&ctx, k);
  return ctx.count;
}

UV* inverse_sigma0_list(UV *count, UV lo, UV hi, UV k)
{
  invsig0_ctx_t ctx;

  *count = 0;
  if (k == 0 || hi < 1 || k > hi) return 0;
  if (k == 1) {
    UV *list;
    if (lo <= 1 && hi >= 1) {
      New(0, list, 1, UV);
      list[0] = 1;
      *count = 1;
      return list;
    }
    return 0;
  }
  if (k == 2) {
    UV *list;
    if (hi < 2) return 0;
    *count = range_prime_sieve(&list, (lo < 2) ? 2 : lo, hi);
    return list;
  }
  if (lo < 1) lo = 1;
  if (lo > hi) return 0;
  if (use_tau_sieve(lo, hi, k))
    return tau_sieve_chunked(count, lo, hi, k, 0);

  memset(&ctx, 0, sizeof(ctx));
  ctx.lo = lo;
  ctx.hi = hi;
  by_signature(&ctx, k);
  if (ctx.count > 1) {
    size_t len = ctx.count;
    sort_dedup_uv_array(ctx.list, 0, &len);
    ctx.count = len;
  }
  *count = ctx.count;
  return ctx.list;
}
