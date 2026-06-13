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
 * 1) SIEVE
 *    A lightweight segmented sieve.  It tracks only the divisor count and
 *    residual value for each integer in the segment, not the full factor list.
 *    When sieving to sqrt(hi) is cheap, residual values are prime.  Otherwise
 *    we sieve through cuberoot(hi), after which each residual has at most two
 *    prime factors, so classifying it as 1, prime, prime square, or semiprime
 *    is enough.
 *
 * 2) GENERATE
 *    We enumerate all prime signatures: sorted-descending exponent lists
 *    (e0 >= e1 >= ... >= er) with (e0+1)*(e1+1)*...*(er+1) == k, factors
 *    drawn from divisors of k.  For each signature we count or list all
 *    integers in the target range with that multiplicative structure.
 *
 *    The counting recursion is a direct translation of the Perl function
 *    count_prime_signature_numbers() including all its optimizations:
 *
 *    - rem_sum bound: hi = rootint(n/m, rem_sum) where rem_sum is the sum
 *      of ALL remaining exponents.  Every subsequent prime is also >= p,
 *      so p^rem_sum <= n/m is necessary, giving a tighter hi than using
 *      only the current exponent.
 *
 *    - @seen deduplication: repeated exponent values in the remaining
 *      signature are skipped (only the first occurrence is processed),
 *      eliminating symmetric subtrees.
 *
 *    - local_j optimization: j = prime_count(lo - 1) is threaded through
 *      the recursion, avoiding prime_count(lo-1) recomputation.
 *        k==1 base: count += prime_count(hi) - j
 *        k==2 base: count += prime_count(u) - (++local_j) for each outer p
 *        general:   recurse with j = ++local_j (= prime_count(p))
 *
 *    The list path retains the original positional recursion which was
 *    already correct.
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

/******************************************************************************/
/* Sieve path (unchanged from original)                                        */
/******************************************************************************/

static bool full_factor_sieve(UV sqrtn, UV range)
{
  if (sqrtn <    1000000U) return (range > 200000);
  if (sqrtn <   10000000U) return (range > 300000);
  if (sqrtn <   35000000U) return (range > 350000);
  if (sqrtn <  100000000U) return (range > 500000);
  if (sqrtn <  350000000U) return (range > 750000);
  if (sqrtn < 1000000000U) return (range > 1000000);
  if (sqrtn < 3500000000U) return (range > 1500000);
  return (range > 2000000);
}

static int first_multiple_in_range(UV *A, UV lo, UV hi, UV p)
{
  UV a = (lo / p) * p;
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
  if (mult == 0 || t > k / mult) { *tau = 0; return; }
  t *= mult;
  *tau = (k % t != 0) ? 0 : (uint16_t)t;
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
    *tau = is_def_prime(rem)                               ? (uint16_t)k : 0;
  } else if (need == 3) {
    *tau = is_perfect_square(rem)                          ? (uint16_t)k : 0;
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
  full_sieve = full_factor_sieve(sqrthi, hi - lo);
  sievelim = full_sieve ? sqrthi : icbrt(hi);

  len = alloc = 0;
  list = 0;
  New(0, rem, range, UV);
  New(0, tau, range, uint16_t);
  for (j = 0; j < range; j++) { rem[j] = lo + j; tau[j] = 1; }

  START_DO_FOR_EACH_PRIME(2, sievelim) {
    UV A, idx;
    if (first_multiple_in_range(&A, lo, hi, p)) {
      for (idx = A - lo; idx < range; idx += p) {
        UV e = 0;
        if (tau[idx] == 0) continue;
        while ((rem[idx] % p) == 0) { e++; rem[idx] /= p; }
        if (e > 0) tau_update(&tau[idx], k, e + 1);
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
          if (alloc == 0) { alloc = 32; New(0, list, alloc, UV); }
          else            { alloc *= 2; Renew(list, alloc, UV);  }
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
  UV len = 0, alloc = 0, seglo = lo, *list = 0;

  while (seglo <= hi) {
    UV segcount, seghi, *seglist;
    seghi = (hi - seglo + 1 > INVERSE_SIGMA0_TAU_MAX_WIDTH)
          ? seglo + INVERSE_SIGMA0_TAU_MAX_WIDTH - 1 : hi;
    seglist = tau_sieve(&segcount, seglo, seghi, k, count_only);
    if (count_only) {
      len += segcount;
    } else if (segcount > 0) {
      if (len + segcount > alloc) {
        alloc = (alloc == 0) ? 32 : alloc;
        while (len + segcount > alloc) alloc *= 2;
        if (list == 0) New(0, list, alloc, UV);
        else           Renew(list, alloc, UV);
      }
      Copy(seglist, list + len, segcount, UV);
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
      width <= INVERSE_SIGMA0_TAU_MAX_CHUNKED_WIDTH) return 1;
#endif
  return 0;
}

/******************************************************************************/
/* Shared helpers                                                              */
/******************************************************************************/

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

static UV mul_pow_limit(UV m, UV p, UV e, UV limit)
{
  while (e-- > 0) {
    if (p == 0 || m > limit / p) return 0;
    m *= p;
  }
  return m;
}

/******************************************************************************/
/* Partition enumeration (tau_partitions)                                      */
/*                                                                             */
/* Direct translation of tau_partitions() from the Perl scripts.              */
/* Factors k as a product of its own divisors >= 2, in non-decreasing order,  */
/* then reverses to get descending exponents.                                  */
/*                                                                             */
/* Key: the loop is  for d in divs while d <= target  (not d <= sqrt(target)) */
/* so d == target is naturally handled as the "last factor" case.             */
/******************************************************************************/

#define MAX_SIG_LEN MPU_MAX_FACTORS

typedef void (*sig_cb_t)(void *data, const UV *sig, int nsig);

typedef struct {
  UV       max_sum_e;
  UV      *divs;        /* divisors of k >= 2, ascending */
  int      ndivs;
  UV       path[MAX_SIG_LEN];
  int      pathlen;
  sig_cb_t callback;
  void    *cb_data;
} part_ctx_t;

static void partition_recurse(part_ctx_t *ctx, UV target, int min_idx,
                              UV curr_sum_e)
{
  int i;

  if (target == 1) {
    /* Reverse path[] to get descending exponents */
    UV sig[MAX_SIG_LEN];
    int len = ctx->pathlen;
    for (i = 0; i < len; i++) sig[i] = ctx->path[len - 1 - i];
    ctx->callback(ctx->cb_data, sig, len);
    return;
  }

  if (ctx->pathlen >= MAX_SIG_LEN) return;

  for (i = min_idx; i < ctx->ndivs; i++) {
    UV d = ctx->divs[i];
    UV e = d - 1;
    if (d > target)                      break;  /* sorted, none larger useful */
    if (curr_sum_e + e > ctx->max_sum_e) break;  /* e only grows, budget gone */
    if (target % d == 0) {
      ctx->path[ctx->pathlen++] = e;
      partition_recurse(ctx, target / d, i, curr_sum_e + e);
      ctx->pathlen--;
    }
  }
}

static void enumerate_signatures(UV k, UV max_sum_e,
                                 sig_cb_t callback, void *cb_data)
{
  part_ctx_t ctx;
  UV ndivs_all, skip, *all_divs;

  if (k == 0) return;
  if (k == 1) { callback(cb_data, NULL, 0); return; }

  all_divs = divisor_list(k, &ndivs_all, k);
  skip = 0;
  while (skip < ndivs_all && all_divs[skip] < 2) skip++;

  ctx.max_sum_e = max_sum_e;
  ctx.divs      = all_divs + skip;
  ctx.ndivs     = (int)(ndivs_all - skip);
  ctx.pathlen   = 0;
  ctx.callback  = callback;
  ctx.cb_data   = cb_data;

  partition_recurse(&ctx, k, 0, 0);
  Safefree(all_divs);
}

/******************************************************************************/
/* Count path                                                                  */
/*                                                                             */
/* count_sig(n, m, lo, sig, nsig, j):                                         */
/*   Count x <= n of the form  x = m * p0^s0 * p1^s1 * ...                   */
/*   with p0 < p1 < ... all primes >= lo.                                     */
/*   j = prime_count(lo - 1)  [invariant maintained by all callers]           */
/******************************************************************************/

static UV count_sig(UV n, UV m, UV lo, const UV *sig, int nsig, UV j);

static UV count_sig(UV n, UV m, UV lo, const UV *sig, int nsig, UV j)
{
  UV rem_sum, hi, count, i;

  if (nsig <= 0 || m == 0) return 0;

  rem_sum = 0;
  for (i = 0; i < (UV)nsig; i++) rem_sum += sig[i];
  if (rem_sum == 0) return 0;
  if (n < m) return 0;
  hi = rootint(n / m, rem_sum);

  if (lo > hi) return 0;

  count = 0;

  for (i = 0; i < (UV)nsig; i++) {
    UV e = sig[i];
    UV p, new_sig[MAX_SIG_LEN];
    int nk, ni;

    if (i > 0 && sig[i] == sig[i - 1]) continue;

    nk = 0;
    for (ni = 0; ni < nsig; ni++)
      if (ni != (int)i) new_sig[nk++] = sig[ni];

    if (nsig == 1) {
      UV pc_hi = prime_count(hi);
      if (pc_hi > j) count += pc_hi - j;

    } else if (nsig == 2) {
      UV e2 = new_sig[0];
      UV local_j = j;
      for (p = lo; p <= hi; p = next_prime(p)) {
        UV t = mul_pow_limit(m, p, e, n);
        UV u, pc_u;
        if (t == 0) break;
        u    = (e2 == 1) ? n / t : rootint(n / t, e2);
        pc_u = prime_count(u);
        ++local_j;
        if (pc_u > local_j) count += pc_u - local_j;
      }

    } else {
      UV local_j = j;
      for (p = lo; p <= hi;) {
        UV t = mul_pow_limit(m, p, e, n);
        UV r;
        if (t == 0) break;
        r = next_prime(p);
        ++local_j;
        count += count_sig(n, t, r, new_sig, nk, local_j);
        p = r;
      }
    }
  }

  return count;
}

/******************************************************************************/
/* Count callback and per-signature count-in-range                            */
/******************************************************************************/

typedef struct {
  UV n;
  UV count;
} count_cb_t;

static void count_callback(void *data, const UV *sig, int nsig)
{
  count_cb_t *d = (count_cb_t *)data;
  UV sum_e, i;

  if (nsig == 0) { if (d->n >= 1) d->count++; return; }

  sum_e = 0;
  for (i = 0; i < (UV)nsig; i++) sum_e += sig[i];
  if (sum_e == 0 || sum_e > log2floor(d->n)) return;

  /* Initial call: m=1, lo=2, j=prime_count(1)=0 */
  d->count += count_sig(d->n, 1, 2, sig, nsig, 0);
}

/******************************************************************************/
/* List path                                                                   */
/******************************************************************************/

typedef struct {
  UV      lo;
  UV      hi;
  UV      count;
  UV      alloc;
  UV     *list;
  UV     *divs;
  UV      ndivs;
  int     nexps;
  uint8_t exps[MPU_MAX_FACTORS];
} invsig0_ctx_t;

static int uv_in_list(UV v, const UV *list, int len)
{
  int i;
  for (i = 0; i < len; i++)
    if (list[i] == v) return 1;
  return 0;
}

static int min_signature_le_hi(const uint8_t *exps, int nexps, UV hi)
{
  UV n = 1, p = 0;
  int i;
  for (i = 0; i < nexps; i++) {
    p = (i + 1 < NPRIMES_SMALL) ? primes_small[i + 1] : next_prime(p);
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

  same = 1;
  while (idx + same < ctx->nexps && ctx->exps[idx + same] == e) same++;
  for (j = idx + same; j < ctx->nexps; j++)
    lower_sum += ctx->exps[j];
  if (lower_sum >= BITS_PER_WORD) return 0;
  denom = UVCONST(1) << lower_sum;
  base  = maxpow / denom;
  if (base == 0) return 0;
  root = (uint32_t)e * (uint32_t)same;
  return rootint(base, root);
}

static void record(invsig0_ctx_t *ctx, UV n)
{
  if (n < ctx->lo || n > ctx->hi) return;
  if (ctx->count >= ctx->alloc) {
    if (ctx->alloc == 0) { ctx->alloc = 32; New(0, ctx->list, ctx->alloc, UV); }
    else                 { ctx->alloc *= 2; Renew(ctx->list, ctx->alloc, UV);  }
  }
  ctx->list[ctx->count++] = n;
}

static void recurse(invsig0_ctx_t *ctx, int idx, UV n, UV *used)
{
  UV minpow, maxpow, pl, ph, *P, pcnt, i;
  uint8_t e = ctx->exps[idx];

  if (idx == ctx->nexps - 1) {
    minpow = (n >= ctx->lo) ? 1 : ceil_div_uv(ctx->lo, n);
    maxpow = ctx->hi / n;
    pl = (e == 1) ? minpow : ceil_root_uv(minpow, e);
    ph = (e == 1) ? maxpow : rootint(maxpow, e);
    if (pl < 2) pl = 2;
    if (idx > 0 && e == ctx->exps[idx - 1] && pl <= used[idx - 1])
      pl = used[idx - 1] + 1;
    if (pl > ph) return;

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
  if (idx > 0 && e == ctx->exps[idx - 1] && pl <= used[idx - 1])
    pl = used[idx - 1] + 1;
  if (pl > ph) return;

  pcnt = range_prime_sieve(&P, pl, ph);
  for (i = 0; i < pcnt; i++) {
    UV p = P[i];
    if (!uv_in_list(p, used, idx)) {
      UV v = mul_pow_limit(n, p, e, ctx->hi);
      if (v != 0) {
        used[idx] = p;
        recurse(ctx, idx + 1, v, used);
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
    UV e = parts[len - 1 - i] - 1;
    if (e == 0 || e > maxexp) return;
    ctx->exps[i] = (uint8_t)e;
  }
  if (!min_signature_le_hi(ctx->exps, len, ctx->hi)) return;
  ctx->nexps = len;
  recurse(ctx, 0, 1, used);
}

static void partitions(invsig0_ctx_t *ctx, UV rem, UV min_d,
                       UV *parts, int len)
{
  UV i;
  if (len >= MPU_MAX_FACTORS) return;
  for (i = 0; i < ctx->ndivs; i++) {
    UV d = ctx->divs[i];
    if (d < 2 || d < min_d) continue;
    if (d > rem / d) break;
    if (rem % d == 0) {
      parts[len] = d;
      partitions(ctx, rem / d, d, parts, len + 1);
    }
  }
  if (rem >= 2 && rem >= min_d) {
    parts[len] = rem;
    process_signature(ctx, parts, len + 1);
  }
}

static void by_signature(invsig0_ctx_t *ctx, UV k)
{
  UV parts[MPU_MAX_FACTORS];
  ctx->divs = divisor_list(k, &ctx->ndivs, k);
  partitions(ctx, k, 2, parts, 0);
  Safefree(ctx->divs);
}

/******************************************************************************/
/* Public API                                                                  */
/******************************************************************************/

UV inverse_sigma0_count(UV lo, UV hi, UV k)
{
  count_cb_t cbd_hi, cbd_lo;
  UV count_hi, count_lo;

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

  /* count(lo,hi,k) = count_upto(hi,k) - count_upto(lo-1,k) */
  cbd_hi.n = hi;     cbd_hi.count = 0;
  enumerate_signatures(k, log2floor(hi), count_callback, &cbd_hi);
  count_hi = cbd_hi.count;

  if (lo <= 1) return count_hi;

  cbd_lo.n = lo - 1; cbd_lo.count = 0;
  enumerate_signatures(k, log2floor(lo - 1), count_callback, &cbd_lo);
  count_lo = cbd_lo.count;

  return (count_hi > count_lo) ? count_hi - count_lo : 0;
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
      list[0] = 1; *count = 1;
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
