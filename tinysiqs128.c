/*============================================================================

  Fixed-width self-initializing quadratic sieve for inputs through 128 bits.

  This is the no-GMP specialization of Math::Prime::Util::GMP's SIQS.  It
  uses the same normalized polynomial

      q(x) = ((D*x + B)^2 - kN) / D,       D = d*A, d in {1,2},

  the same byte sieve and candidate resieve, one-large-prime pairing, and
  exact dense linear algebra.  Every over-width setup expression is either
  rearranged exactly or checked before it is formed.

  Copyright (c) 2026 Dana Jacobsen
  Written by Dana Jacobsen, September 2026, with assistance from OpenAI Codex.

============================================================================*/

#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "tinysiqs128.h"

#if HAVE_FACTOR128

#ifndef UINT8_MAX
# define UINT8_MAX 255U
#endif
#ifndef UINT16_MAX
# define UINT16_MAX 65535U
#endif
#ifndef UINT32_MAX
# define UINT32_MAX 4294967295U
#endif
#ifndef UINT64_C
# define UINT64_C(value) ((uint64_t)(value))
#endif
#ifndef M_LN2
# define M_LN2 0.69314718055994530942
#endif

#define TS_U128_MAX             ((uint128_t)-1)
#define TS_I128_MAX_MAG         (TS_U128_MAX >> 1)
#define TS_NO_INDEX             UINT32_MAX
#define TS_NO_ROOT              UINT32_MAX
#define TS_SIEVE_ALIGN          256U
#define TS_MAX_EXTRA_RELS       512U
#define TS_RETRY_BATCH_MAX      128U
#define TS_A_TOLERANCE_DEFAULT  8U
#define TS_MULTIPLIER_MAX       255U
#define TS_MULTIPLIER_CAPACITY  128U
#define TS_POSTFILTER_MAX_SMALL 256U
#define TS_RESIEVE_COEFFICIENT  8U
#define TS_RESIEVE_FLOOR        1024U
#define TS_MAX_Q                6U
#define TS_MAX_FB               768U
#define TS_MAX_RELATIONS        2048U
/* A signed-128 q value has at most 25 distinct prime factors.  Adding sign,
 * d, and six A primes keeps a raw relation below 40 entries; the union of
 * two raw relations remains below 72. */
#define TS_MAX_RAW_FACTORS      40U
#define TS_MAX_FULL_FACTORS     72U
#define TS_MAX_FAMILIES_Q2      4096U
#define TS_MAX_FAMILIES_OTHER   131072U
#define TS_HASH_EMPTY           UINT64_C(0)

#if defined(DEBUGGING) || defined(TINYSIQS128_DEBUG)
# define TS_DEBUG 1
#else
# define TS_DEBUG 0
#endif
#ifdef TINYSIQS128_DIAGNOSTIC
# define TS_DIAGNOSTIC 1
#else
# define TS_DIAGNOSTIC 0
#endif

typedef struct {
  uint32_t p;
  uint32_t sqrt_kn;
  uint8_t logp;
  uint8_t in_a;
} ts_fb_t;

typedef struct {
  uint16_t row;
  uint16_t exponent;
} ts_factor_t;

typedef struct {
  uint128_t y;
  ts_factor_t factors[TS_MAX_RAW_FACTORS];
  uint16_t nfactors;
  uint32_t lp;
} ts_raw_relation_t;

typedef struct {
  uint128_t y;
  uint32_t factor_offset;
  uint16_t nfactors;
} ts_full_relation_t;

typedef struct {
  uint32_t label;
  ts_raw_relation_t relation;
} ts_lp_anchor_t;

typedef struct {
  ts_lp_anchor_t *anchors;
  uint32_t count;
  uint32_t alloc;
} ts_lp_state_t;

typedef struct {
  uint64_t *slots;
  uint32_t count;
  uint32_t alloc;
} ts_hashset_t;

typedef struct {
  uint64_t state;
} ts_rng_t;

typedef struct {
  const char *name;
  uint8_t first_bits;
  uint8_t last_bits;
  uint8_t q_count;
  uint8_t lp_multiplier;
  uint8_t bias_base;
  uint8_t bias_step_bits;
  uint8_t bias_origin;
  uint8_t sieve_free_units;
  uint16_t fb_floor;
  uint16_t relation_extra;
  uint8_t a_tolerance;
  uint32_t fixed_M;
  double fb_base;
  double fb_step;
  uint8_t fb_origin;
  double interval_base;
  double interval_step;
  uint8_t interval_origin;
  double sieve_start_exponent;
} ts_policy_band_t;

typedef struct {
  const char *name;
  uint32_t bits;
  uint32_t fb_size;
  uint32_t M;
  uint32_t q_count;
  uint32_t d;
  uint32_t lp_multiplier;
  uint32_t relation_extra;
  uint32_t target_relations;
  uint32_t a_tolerance;
  uint32_t stage1_bias;
  uint32_t sieve_free_units;
  uint32_t sieve_start;
  double sieve_start_exponent;
  uint8_t sieve_initial;
  double log_scale;
  double sieve_hit_bound;
  uint32_t max_families;
} ts_parameters_t;

typedef struct {
  uint64_t A;
  uint64_t D;
  int128_t B;
  int128_t C;
  uint64_t target_A;
  uint64_t H[TS_MAX_Q];
  uint16_t a_index[TS_MAX_Q];
  uint32_t *corrections;
  uint32_t q_count;
  uint32_t b_index;
  uint32_t b_limit;
} ts_poly_t;

typedef struct {
  int32_t x;
  uint32_t first_hit;
  uint8_t sieve_score;
} ts_candidate_t;

typedef struct {
  uint16_t fb_index;
  uint32_t next;
} ts_hit_t;

typedef struct {
  uint128_t n;
  uint128_t ninv;
  uint128_t r2;
  uint128_t one;
} ts_mont_t;

typedef struct {
  unsigned long *data;
  unsigned long weight;
  unsigned long orig;
} ts_la_col_t;

typedef struct {
  uint128_t n;
  uint128_t kn;
  uint32_t multiplier;
  ts_parameters_t params;
  ts_fb_t *fb;
  uint32_t *prime;
  uint32_t *root1;
  uint32_t *root2;
  uint32_t *reciprocal;
  uint32_t largest_fb_prime;
  uint32_t resieve_one_subtract_index;
  uint8_t active_sieve_initial;
  uint8_t *sieve;
  uint32_t sieve_length;
  uint16_t *candidate_at;
  uint32_t *candidate_at_wide;
  int candidate_wide;
  ts_candidate_t *candidates;
  uint32_t candidate_count;
  uint32_t candidate_alloc;
  ts_hit_t *hits;
  uint32_t hit_count;
  uint32_t hit_alloc;
  ts_full_relation_t *full;
  uint32_t full_count;
  uint32_t full_alloc;
  ts_factor_t *full_factors;
  uint32_t full_factor_count;
  uint32_t full_factor_alloc;
  ts_lp_state_t one_lp;
  ts_hashset_t a_hashes;
  ts_rng_t poly_rng;
  ts_mont_t mont;
  int mont_ready;
  uint32_t *factor_counts;
  uint16_t *factor_touched;
  uint32_t factor_touched_count;
  uint32_t matrix_next_target;
  uint32_t matrix_last_count;
  uint32_t matrix_retry_batch;
  uint32_t matrix_target_limit;
  uint128_t factor;
  uint64_t total_candidates;
  uint64_t accepted_smooth;
  uint64_t accepted_partial;
} ts_ctx_t;

/*----------------------------------------------------------------------------
 * Allocation and scalar helpers
 *----------------------------------------------------------------------------*/

static size_t ts_array_bytes(size_t count, size_t size) {
  if (size != 0 && count > (size_t)-1 / size)
    croak("tinysiqs128: allocation size overflow");
  return count * size;
}

static void *ts_malloc_array(size_t count, size_t size) {
  size_t bytes = ts_array_bytes(count, size);
  void *p = malloc(bytes ? bytes : 1);
  if (p == NULL)
    croak("tinysiqs128: unable to allocate memory");
  return p;
}

static void *ts_calloc_array(size_t count, size_t size) {
  void *p;
  (void)ts_array_bytes(count, size);
  p = calloc(count ? count : 1, size ? size : 1);
  if (p == NULL)
    croak("tinysiqs128: unable to allocate memory");
  return p;
}

static void *ts_realloc_array(void *old, size_t count, size_t size) {
  size_t bytes = ts_array_bytes(count, size);
  void *p = realloc(old, bytes ? bytes : 1);
  if (p == NULL)
    croak("tinysiqs128: unable to grow allocation");
  return p;
}

static uint32_t ts_bits128(uint128_t n) {
  uint64_t hi = (uint64_t)(n >> 64);
  uint64_t lo = (uint64_t)n;
  uint32_t bits = 0;
  uint64_t v = hi ? hi : lo;
  if (v == 0)
    return 0;
#if defined(__GNUC__) || defined(__clang__)
  bits = 64U - (uint32_t)__builtin_clzll(v);
#else
  while (v != 0) { bits++; v >>= 1; }
#endif
  return hi ? bits + 64U : bits;
}

static double ts_log128(uint128_t n) {
  uint64_t hi = (uint64_t)(n >> 64);
  uint64_t lo = (uint64_t)n;
  if (hi == 0)
    return log((double)lo);
  return log((double)hi + ldexp((double)lo, -64)) + 64.0 * M_LN2;
}

static uint128_t ts_abs_i128(int128_t value) {
  uint128_t u = (uint128_t)value;
  return value < 0 ? (uint128_t)0 - u : u;
}

static uint32_t ts_signed_mod_u32(int128_t value, uint32_t p) {
  uint32_t r = (uint32_t)(ts_abs_i128(value) % p);
  return value < 0 && r != 0 ? p - r : r;
}

static uint128_t ts_signed_mod_n(int128_t value, uint128_t n) {
  uint128_t r = ts_abs_i128(value) % n;
  return value < 0 && r != 0 ? n - r : r;
}

static uint64_t ts_mix64(uint64_t x) {
  x ^= x >> 30;
  x *= UINT64_C(0xbf58476d1ce4e5b9);
  x ^= x >> 27;
  x *= UINT64_C(0x94d049bb133111eb);
  x ^= x >> 31;
  return x ? x : UINT64_C(0x9e3779b97f4a7c15);
}

static uint64_t ts_rand64(ts_rng_t *rng) {
  uint64_t x = rng->state;
  x ^= x >> 12;
  x ^= x << 25;
  x ^= x >> 27;
  rng->state = x;
  return x * UINT64_C(2685821657736338717);
}

static uint32_t ts_rand_range(ts_rng_t *rng, uint32_t limit) {
  return limit == 0 ? 0 : (uint32_t)(ts_rand64(rng) % limit);
}

static uint32_t ts_ctz32(uint32_t n) {
  if (n == 0)
    croak("tinysiqs128: ctz of zero");
#if defined(__GNUC__) || defined(__clang__)
  return (uint32_t)__builtin_ctz(n);
#else
  uint32_t count = 0;
  while (!(n & 1U)) { n >>= 1; count++; }
  return count;
#endif
}

static unsigned int ts_ctz64(uint64_t n) {
  if (n == 0)
    croak("tinysiqs128: ctz of zero");
#if defined(__GNUC__) || defined(__clang__)
  return (unsigned int)__builtin_ctzll(n);
#else
  unsigned int count = 0;
  while (!(n & 1U)) { n >>= 1; count++; }
  return count;
#endif
}

static uint64_t ts_isqrt128(uint128_t n) {
  uint128_t r, next;
  uint32_t bits;
  if (n == 0)
    return 0;
  bits = ts_bits128(n);
  r = (uint128_t)1 << ((bits + 1U) / 2U);
  for (;;) {
    next = (r + n / r) >> 1;
    if (next >= r)
      break;
    r = next;
  }
  while (r > n / r)
    r--;
  return (uint64_t)r;
}

static uint128_t ts_gcd128(uint128_t a, uint128_t b) {
  while (b != 0) {
    uint128_t t = a % b;
    a = b;
    b = t;
  }
  return a;
}

static uint128_t ts_addmod128(uint128_t a, uint128_t b, uint128_t n) {
  uint128_t nb = n - b;
  return a >= nb ? a - nb : a + b;
}

static uint128_t ts_submod128(uint128_t a, uint128_t b, uint128_t n) {
  return a >= b ? a - b : n - (b - a);
}

static uint128_t ts_pow2mod128(unsigned int exponent, uint128_t n) {
  uint128_t r = 1;
  while (exponent-- != 0)
    r = ts_addmod128(r, r, n);
  return r;
}

/* Low 128 bits of a full product. */
static uint128_t ts_lo128(uint128_t a, uint128_t b) {
  uint64_t a0 = (uint64_t)a, a1 = (uint64_t)(a >> 64);
  uint64_t b0 = (uint64_t)b, b1 = (uint64_t)(b >> 64);
  uint128_t lo = (uint128_t)a0 * b0;
  lo += ((uint128_t)a0 * b1) << 64;
  lo += ((uint128_t)a1 * b0) << 64;
  return lo;
}

static void ts_mul256(uint128_t a, uint128_t b,
                      uint128_t *lo_out, uint128_t *hi_out) {
  uint64_t a0 = (uint64_t)a, a1 = (uint64_t)(a >> 64);
  uint64_t b0 = (uint64_t)b, b1 = (uint64_t)(b >> 64);
  uint128_t p00 = (uint128_t)a0 * b0;
  uint128_t p01 = (uint128_t)a0 * b1;
  uint128_t p10 = (uint128_t)a1 * b0;
  uint128_t p11 = (uint128_t)a1 * b1;
  uint128_t mid = p01 + p10;
  int middle_carry = mid < p01;
  uint128_t lo = p00 + (mid << 64);
  int low_carry = lo < p00;
  *lo_out = lo;
  *hi_out = p11 + (mid >> 64)
          + ((uint128_t)middle_carry << 64) + (uint128_t)low_carry;
}

static uint128_t ts_mulhi128(uint128_t a, uint128_t b) {
  uint128_t lo, hi;
  ts_mul256(a, b, &lo, &hi);
  (void)lo;
  return hi;
}

static uint128_t ts_mont_mul(uint128_t a, uint128_t b,
                             const ts_mont_t *ctx) {
  uint128_t lo, hi, m, mn_hi, t;
  int carry, overflow;
  ts_mul256(a, b, &lo, &hi);
  m = ts_lo128(lo, ctx->ninv);
  carry = lo != 0;
  mn_hi = ts_mulhi128(m, ctx->n);
  t = hi + mn_hi;
  overflow = t < hi;
  {
    uint128_t old = t;
    t += (uint128_t)carry;
    overflow += t < old;
  }
  return overflow || t >= ctx->n ? t - ctx->n : t;
}

static uint128_t ts_mont_enter(uint128_t a, const ts_mont_t *ctx) {
  return ts_mont_mul(a % ctx->n, ctx->r2, ctx);
}

static uint128_t ts_mont_exit(uint128_t a, const ts_mont_t *ctx) {
  return ts_mont_mul(a, 1, ctx);
}

static void ts_mont_setup(ts_mont_t *ctx, uint128_t n) {
  uint128_t x;
  if (n <= 1 || !(n & 1U))
    croak("tinysiqs128: invalid Montgomery modulus");
  ctx->n = n;
  x = (3 * n) ^ 2;
  x *= 2 - n * x;
  x *= 2 - n * x;
  x *= 2 - n * x;
  x *= 2 - n * x;
  x *= 2 - n * x;
  ctx->ninv = (uint128_t)0 - x;
  ctx->r2 = ts_pow2mod128(256, n);
  ctx->one = ts_mont_enter(1, ctx);
}

static uint128_t ts_mont_pow_u32(uint128_t a, uint32_t exponent,
                                 const ts_mont_t *ctx) {
  uint128_t r = ctx->one;
  while (exponent != 0) {
    if (exponent & 1U)
      r = ts_mont_mul(r, a, ctx);
    exponent >>= 1;
    if (exponent != 0)
      a = ts_mont_mul(a, a, ctx);
  }
  return r;
}

/* Coefficients are kept modulo n in Montgomery form. */
static uint128_t ts_modinv128(uint128_t a, const ts_mont_t *ctx) {
  uint128_t r0 = ctx->n, r1 = a % ctx->n;
  uint128_t t0 = 0, t1 = ctx->one;
  if (r1 == 0)
    return 0;
  while (r1 != 0) {
    uint128_t q = r0 / r1;
    uint128_t r = r0 - q * r1;
    uint128_t qm = ts_mont_enter(q, ctx);
    uint128_t t = ts_submod128(t0, ts_mont_mul(qm, t1, ctx), ctx->n);
    r0 = r1;
    r1 = r;
    t0 = t1;
    t1 = t;
  }
  return r0 == 1 ? ts_mont_exit(t0, ctx) : 0;
}

/* Compare base^exponent with limit without overflowing. */
static int ts_power_compare(uint128_t base, uint32_t exponent,
                            uint128_t limit) {
  uint128_t product = 1;
  uint32_t i;
  for (i = 0; i < exponent; i++) {
    if (base != 0 && product > limit / base)
      return 1;
    product *= base;
  }
  return product < limit ? -1 : product > limit ? 1 : 0;
}

static uint128_t ts_exact_root(uint128_t n, uint32_t exponent) {
  uint32_t bits = ts_bits128(n);
  uint32_t root_bits = (bits + exponent - 1U) / exponent;
  uint128_t lo = 2, hi;
  if (root_bits >= 128)
    hi = n;
  else
    hi = (uint128_t)1 << root_bits;
  while (lo <= hi) {
    uint128_t mid = lo + ((hi - lo) >> 1);
    int cmp = ts_power_compare(mid, exponent, n);
    if (cmp == 0)
      return mid;
    if (cmp < 0)
      lo = mid + 1;
    else {
      if (mid == 0)
        break;
      hi = mid - 1;
    }
  }
  return 0;
}

static uint128_t ts_trial_factor32(uint32_t n) {
  static const uint8_t steps[8] = {4, 2, 4, 2, 4, 6, 2, 6};
  uint32_t p = 7, step_index = 0;
  while ((uint64_t)p * p <= n) {
    if (n % p == 0)
      return p;
    p += steps[step_index];
    step_index = (step_index + 1U) & 7U;
  }
  return 0;
}

/* Exhausted HOLF32 is faster on known-composite, post-trial inputs near
 * 32 bits, but cannot replace either this guaranteed trial path or SQUFOF.
 * Keep the low end to two splitters. */

/* Racing SQUFOF is a much better compact splitter than setting up SIQS at
 * the bottom of its range.  Keeping the multiplied input below 2^62 makes
 * every continued-fraction state fit in uint32_t. */
typedef struct {
  int valid;
  uint32_t p;
  uint32_t bn;
  uint32_t qn;
  uint32_t q0;
  uint32_t b0;
  uint32_t iteration;
  uint32_t batch;
} ts_squfof_state_t;

static uint64_t ts_gcd64(uint64_t a, uint64_t b) {
  while (b != 0) {
    uint64_t t = a % b;
    a = b;
    b = t;
  }
  return a;
}

static uint32_t ts_isqrt62(uint64_t n) {
  uint64_t root = (uint64_t)(sqrt((double)n) + 1e-6);
  if (root * root > n)
    root--;
  return (uint32_t)root;
}

static int ts_square_u32(uint32_t n, uint32_t *root_out) {
  uint32_t root;
  /* Quadratic residues modulo 64 reject 81.25% of nonsquares here. */
  if ((UINT64_C(1) << (n & 63U)) &
      UINT64_C(0xfdfdfdedfdfcfdec))
    return 0;
  root = ts_isqrt62(n);
  if ((uint64_t)root * root != n)
    return 0;
  *root_out = root;
  return 1;
}

static uint64_t ts_squfof_unit(uint64_t n, ts_squfof_state_t *state) {
  uint32_t p = state->p, bn = state->bn;
  uint32_t qn = state->qn, q0 = state->q0, b0 = state->b0;
  uint32_t iteration = state->iteration;
  uint32_t limit = iteration + state->batch;

  for (;;) {
    uint32_t root, ro, so, s, bbn;
    uint32_t symmetry_iterations = 0;
    if (iteration & 1U) {
      uint32_t old_p = p, old_qn = qn;
      p = bn * qn - p;
      qn = q0 + bn * (old_p - p);
      q0 = old_qn;
      bn = (b0 + p) / qn;
      iteration++;
    }
    for (;;) {
      uint32_t old_p, old_qn;
      if (iteration >= limit) {
        state->p = p;
        state->bn = bn;
        state->qn = qn;
        state->q0 = q0;
        state->iteration = iteration;
        return 0;
      }
      old_p = p;
      old_qn = qn;
      p = bn * qn - p;
      qn = q0 + bn * (old_p - p);
      q0 = old_qn;
      bn = (b0 + p) / qn;
      iteration++;
      if (ts_square_u32(qn, &root))
        break;

      old_p = p;
      old_qn = qn;
      p = bn * qn - p;
      qn = q0 + bn * (old_p - p);
      q0 = old_qn;
      bn = (b0 + p) / qn;
      iteration++;
    }

    s = root;
    state->iteration = iteration;
    ro = p + s * ((b0 - p) / s);
    so = (uint32_t)((n - (uint64_t)ro * ro) / s);
    bbn = (b0 + ro) / so;

    for (;;) {
      uint32_t old_ro = ro, old_so = so;
      ro = bbn * so - ro;
      so = s + bbn * (old_ro - ro);
      s = old_so;
      bbn = (b0 + ro) / so;
      if (ro == old_ro)
        break;
      if (++symmetry_iterations > 2000000U) {
        state->valid = 0;
        return 0;
      }
    }

    {
      uint64_t factor = ts_gcd64(ro, n);
      if (factor > 1)
        return factor;
    }
  }
}

static uint64_t ts_squfof64(uint64_t n, uint32_t rounds) {
  /* Gower and Wagstaff 2008, section 5.3, with MPU's 13, 17, and 19
   * extensions in analytic-score order.  Large tuning corpora found this
   * complete race both faster and safer than reduced portfolios. */
  static const uint32_t multipliers[] = {
    3U*5U*7U*11U, 3U*5U*7U, 3U*5U*7U*11U*13U, 3U*5U*7U*13U,
    3U*5U*7U*11U*17U, 3U*5U*11U, 3U*5U*7U*17U, 3U*5U,
    3U*5U*7U*11U*19U, 3U*5U*11U*13U, 3U*5U*7U*19U,
    3U*5U*7U*13U*17U, 3U*5U*13U, 3U*7U*11U, 3U*7U, 5U*7U*11U,
    3U*7U*13U, 5U*7U, 3U*5U*17U, 5U*7U*13U, 3U*5U*19U, 3U*11U,
    3U*7U*17U, 3U, 3U*11U*13U, 5U*11U, 3U*7U*19U, 3U*13U, 5U,
    5U*11U*13U, 5U*7U*19U, 5U*13U, 7U*11U, 7U, 3U*17U,
    7U*13U, 11U, 1U
  };
  ts_squfof_state_t states[sizeof(multipliers) / sizeof(multipliers[0])];
  uint32_t rounds_done = 0, racing = (uint32_t)(sizeof(multipliers) /
                                                sizeof(multipliers[0]));
  uint32_t i;
  for (i = 0; i < racing; i++) {
    states[i].valid = -1;
    states[i].iteration = 0;
  }

  while (racing != 0 && rounds_done < rounds) {
    uint32_t count = (uint32_t)(sizeof(multipliers) /
                                sizeof(multipliers[0]));
    for (i = 0; i < count && rounds_done < rounds; i++) {
      ts_squfof_state_t *state = &states[i];
      uint32_t multiplier = multipliers[i];
      uint64_t nn, factor;
      if (state->valid == 0)
        continue;
      if ((UINT64_MAX >> 2) / multiplier < n) {
        state->valid = 0;
        racing--;
        continue;
      }
      nn = n * multiplier;
      if (state->valid < 0) {
        uint32_t sqrt_nn = ts_isqrt62(nn);
        state->valid = 1;
        state->q0 = 1;
        state->b0 = sqrt_nn;
        state->p = sqrt_nn;
        state->qn = (uint32_t)(nn - (uint64_t)sqrt_nn * sqrt_nn);
        if (state->qn == 0) {
          factor = ts_gcd64(sqrt_nn, n);
          if (factor > 1 && factor < n)
            return factor;
          state->valid = 0;
          racing--;
          continue;
        }
        state->bn = (2U * sqrt_nn) / state->qn;
        state->batch = (uint32_t)(sqrt((double)sqrt_nn) / 16.0);
        if (state->batch < 20U)
          state->batch = 20U;
        if (state->batch > rounds)
          state->batch = rounds;
      }
      if (racing == 1U)
        state->batch = rounds - rounds_done;
      factor = ts_squfof_unit(nn, state);
      if (factor > 1) {
        uint64_t reduced = ts_gcd64(factor, n);
        if (reduced > 1 && reduced < n && n % reduced == 0)
          return reduced;
        state->valid = 0;
      }
      if (state->valid == 0)
        racing--;
      rounds_done += state->batch;
    }
  }
  return 0;
}

/*----------------------------------------------------------------------------
 * Small-prime arithmetic
 *----------------------------------------------------------------------------*/

static uint32_t ts_inverse_u32(uint32_t a, uint32_t p) {
  uint32_t dividend = p, divisor = a;
  uint32_t ps1 = 1, ps2 = 0, parity = 0;
  if (a == 0)
    return 0;
  while (divisor > 1) {
    uint32_t q = dividend / divisor;
    uint32_t rem = dividend - q * divisor;
    uint32_t next = q * ps1 + ps2;
    dividend = divisor;
    divisor = rem;
    ps2 = ps1;
    ps1 = next;
    parity ^= 1U;
  }
  return parity ? p - ps1 : ps1;
}

static uint32_t ts_powmod_u32(uint32_t a, uint32_t exponent, uint32_t p) {
  uint32_t r = 1;
  while (exponent != 0) {
    if (exponent & 1U)
      r = (uint32_t)((uint64_t)r * a % p);
    exponent >>= 1;
    if (exponent != 0)
      a = (uint32_t)((uint64_t)a * a % p);
  }
  return r;
}

static int ts_legendre_u32(uint32_t a, uint32_t p) {
  uint32_t value;
  if (a == 0)
    return 0;
  value = ts_powmod_u32(a, (p - 1U) >> 1, p);
  return value == 1 ? 1 : value == p - 1U ? -1 : 0;
}

static int ts_jacobi_odd_u32(uint32_t a, uint32_t n) {
  int sign = 1;
  a %= n;
  while (a != 0) {
    while (!(a & 1U)) {
      a >>= 1;
      if ((n & 7U) == 3U || (n & 7U) == 5U)
        sign = -sign;
    }
    {
      uint32_t t = a;
      a = n;
      n = t;
    }
    if ((a & 3U) == 3U && (n & 3U) == 3U)
      sign = -sign;
    a %= n;
  }
  return n == 1 ? sign : 0;
}

static uint32_t ts_sqrtmod_u32(uint32_t n, uint32_t p) {
  uint32_t q, s, z, c, x, t, m;
  if (n == 0)
    return 0;
  if ((p & 3U) == 3U)
    return ts_powmod_u32(n, (p + 1U) / 4U, p);
  s = ts_ctz32(p - 1U);
  q = (p - 1U) >> s;
  for (z = 2; ts_legendre_u32(z, p) != -1; z++)
    ;
  c = ts_powmod_u32(z, q, p);
  x = ts_powmod_u32(n, (q + 1U) >> 1, p);
  t = ts_powmod_u32(n, q, p);
  m = s;
  while (t != 1) {
    uint32_t i = 1, b = (uint32_t)((uint64_t)t * t % p);
    while (i < m && b != 1) {
      b = (uint32_t)((uint64_t)b * b % p);
      i++;
    }
    if (i >= m)
      croak("tinysiqs128: Tonelli-Shanks failed");
    b = c;
    {
      uint32_t j;
      for (j = m - i - 1U; j != 0; j--)
        b = (uint32_t)((uint64_t)b * b % p);
    }
    x = (uint32_t)((uint64_t)x * b % p);
    b = (uint32_t)((uint64_t)b * b % p);
    t = (uint32_t)((uint64_t)t * b % p);
    c = b;
    m = i;
  }
  return x;
}

static int ts_is_prime_u32(uint32_t n) {
  uint32_t p;
  if (n < 2)
    return 0;
  if (!(n & 1U))
    return n == 2;
  if (n % 3U == 0)
    return n == 3;
  for (p = 5; (uint64_t)p * p <= n; p += 6)
    if (n % p == 0 || n % (p + 2U) == 0)
      return 0;
  return 1;
}

static uint32_t ts_next_prime_u32(uint32_t n) {
  if (n < 2)
    return 2;
  n++;
  if (!(n & 1U))
    n++;
  while (!ts_is_prime_u32(n))
    n += 2;
  return n;
}

/*----------------------------------------------------------------------------
 * Policies, multiplier, factor base, and hash tables
 *----------------------------------------------------------------------------*/

static const ts_policy_band_t ts_primary_policies[] = {
  /* Inputs through 31 bits use the wheel splitter, while 32 bits first use
   * compact SQUFOF.  The smaller FB32 q=1 primary keeps 33--36 bits compact;
   * the floor-40 interval recoveries below absorb its rare underfills. */
  { "smooth_k1_q1_33_36_4k_fb32", 33, 36, 1, 1, 0, 0, 0,
    8, 32, 2, 8, 4096,
    0.315, 0.0, 65, 0.0, 0.0, 1, 0.0 },
  { "smooth_k1_q2_low", 37, 41, 2, 1, 0, 0, 0,
    8, 40, 2, 8, 0,
    0.315, 0.0, 65, 0.5, 0.0, 65, 0.0 },
  { "smooth_k1_q2_fb36_bias3", 42, 49, 2, 1, 3, 0, 0,
    8, 36, 4, 32, 4096,
    0.315, 0.0, 42, 0.0, 0.0, 42, 0.0 },
  { "smooth_k1_q3_floor40", 50, 64, 3, 1, 0, 0, 0,
    8, 40, 4, 8, 0,
    0.315, 0.0, 65, 0.0, 0.041666666666666667, 50, 0.0 },
  { "smooth_k1_q3_fb_low", 65, 80, 3, 1, 6, 0, 0,
    8, 48, 4, 8, 0,
    0.306956091, 0.0000655453, 65,
    0.0, 0.041666666666666667, 50, 0.0 },
  { "smooth_k1_q4_fb_low", 81, 95, 4, 1, 6, 0, 0,
    8, 48, 4, 8, 0,
    0.311589139, 0.0000229912, 81,
    0.5, 0.05, 65, 0.0 },
  { "one_lp_k2_q5_geometry_ramp", 96, 106, 5, 2, 10, 3, 96,
    8, 160, 2, 8, 0,
    0.315, -0.0005929708, 96,
    1.7, 0.08, 96, 0.45 },
  { "one_lp_k2_q5_geometry", 107, 116, 5, 2, 14, 0, 0,
    8, 160, 2, 8, 0,
    0.309104504, 0.0000317423, 107,
    2.5375, 0.0375, 107, 0.45 },
  { "one_lp_k4_q6_geometry_ramp", 117, 128, 6, 4, 16, 0, 0,
    8, 160, 0, 8, 0,
    0.309419667, 0.00046502775, 117,
    3.125, 0.0, 117, 0.45 }
};

static const ts_policy_band_t ts_recovery_q1_8k = {
  "smooth_k1_q1_low_recovery_8k", 33, 36, 1, 1, 0, 0, 0,
  8, 40, 2, 8, 8192,
  0.315, 0.0, 65, 0.0, 0.0, 1, 0.0
};
static const ts_policy_band_t ts_recovery_q1_16k = {
  "smooth_k1_q1_low_recovery_16k", 33, 36, 1, 1, 0, 0, 0,
  8, 40, 2, 8, 16384,
  0.315, 0.0, 65, 0.0, 0.0, 1, 0.0
};
static const ts_policy_band_t ts_recovery_q1_legacy = {
  "smooth_k1_q1_low_recovery_legacy", 33, 36, 1, 1, 0, 0, 0,
  8, 40, 2, 8, 0,
  0.315, 0.0, 65, 4.0, 0.0, 1, 0.0
};
static const ts_policy_band_t ts_recovery_42_96k = {
  "one_lp_k60_q1_low_recovery_96k", 42, 49, 1, 60, 0, 0, 0,
  8, 40, 4, 8, 96U * 1024U,
  0.315, 0.0, 42, 0.0, 0.0, 42, 0.0
};
static const ts_policy_band_t ts_recovery_42_192k = {
  "one_lp_k60_q1_low_recovery_192k", 42, 49, 1, 60, 0, 0, 0,
  8, 40, 4, 8, 192U * 1024U,
  0.315, 0.0, 42, 0.0, 0.0, 42, 0.0
};
static const ts_policy_band_t ts_recovery_42_384k = {
  "one_lp_k60_q1_low_recovery_384k", 42, 49, 1, 60, 0, 0, 0,
  8, 40, 4, 8, 384U * 1024U,
  0.315, 0.0, 42, 0.0, 0.0, 42, 0.0
};
static const ts_policy_band_t ts_recovery_42_1m = {
  "one_lp_k60_q1_low_recovery_1m", 42, 49, 1, 60, 0, 0, 0,
  8, 40, 4, 8, 1U << 20,
  0.315, 0.0, 42, 0.0, 0.0, 42, 0.0
};

static const ts_policy_band_t *ts_primary_policy(uint32_t bits) {
  uint32_t i;
  for (i = 0; i < (uint32_t)(sizeof(ts_primary_policies) /
                              sizeof(ts_primary_policies[0])); i++)
    if (bits >= ts_primary_policies[i].first_bits &&
        bits <= ts_primary_policies[i].last_bits)
      return &ts_primary_policies[i];
  return NULL;
}

static int ts_resolve_parameters(ts_parameters_t *p, uint128_t n,
                                 const ts_policy_band_t *band) {
  double ln_n, ln_term, fb, interval;
  int32_t distance;
  memset(p, 0, sizeof(*p));
  p->bits = ts_bits128(n);
  if (band == NULL)
    band = ts_primary_policy(p->bits);
  if (band == NULL || p->bits < band->first_bits ||
      p->bits > band->last_bits)
    return 0;
  p->name = band->name;
  p->q_count = band->q_count;
  p->lp_multiplier = band->lp_multiplier;
  p->relation_extra = band->relation_extra;
  p->a_tolerance = band->a_tolerance;
  p->stage1_bias = band->bias_base;
  if (band->bias_step_bits != 0)
    p->stage1_bias += (p->bits - band->bias_origin) /
                      band->bias_step_bits;
  p->sieve_free_units = band->sieve_free_units;
  p->sieve_start_exponent = band->sieve_start_exponent;

  ln_n = (double)p->bits * M_LN2;
  ln_term = sqrt(ln_n * log(ln_n));
  distance = (int32_t)p->bits - (int32_t)band->fb_origin;
  /* The exponent applies to sqrt(log N log log N), not directly to e. */
  fb = exp((band->fb_base + band->fb_step * (double)distance) * ln_term);
  if (fb < band->fb_floor)
    fb = band->fb_floor;
  if (fb > TS_MAX_FB)
    return 0;
  p->fb_size = (uint32_t)fb;

  if (band->fixed_M != 0) {
    p->M = band->fixed_M;
  } else {
    double scale;
    distance = (int32_t)p->bits - (int32_t)band->interval_origin;
    scale = band->interval_base +
            band->interval_step * (double)distance;
    interval = (6144.0 + exp(0.37 * ln_term)) * scale;
    p->M = interval > 0.0 ? (uint32_t)interval : 0;
    p->M = (p->M + TS_SIEVE_ALIGN - 1U) & ~(TS_SIEVE_ALIGN - 1U);
    if (p->M < 4096)
      p->M = 4096;
  }
  if (p->M < 4096 || (p->M & (TS_SIEVE_ALIGN - 1U)) != 0)
    return 0;
  p->target_relations = p->fb_size + 1U + p->relation_extra;
  p->max_families = p->q_count == 1 ? 1U
                  : p->q_count == 2 ? TS_MAX_FAMILIES_Q2
                                    : TS_MAX_FAMILIES_OTHER;
  return 1;
}

static int ts_squarefree_small(uint32_t n) {
  uint32_t p;
  for (p = 3; p * p <= n; p += 2)
    if (n % (p * p) == 0)
      return 0;
  return 1;
}

static uint32_t ts_choose_multiplier(uint128_t n, uint32_t fb_size) {
  uint32_t kval[TS_MULTIPLIER_CAPACITY];
  uint16_t accepted[TS_MULTIPLIER_CAPACITY];
  double scores[TS_MULTIPLIER_CAPACITY];
  uint128_t quotient = TS_U128_MAX / n;
  uint32_t max_k = quotient > TS_MULTIPLIER_MAX
                 ? TS_MULTIPLIER_MAX : (uint32_t)quotient;
  uint32_t kcount = 0, wanted = fb_size / 20U;
  uint32_t k, i, unfinished, p = 2;
  uint32_t best = 1;
  double best_score;
  if (wanted < 1)
    wanted = 1;
  if (wanted > fb_size)
    wanted = fb_size;
  memset(accepted, 0, sizeof(accepted));
  for (k = 1; k <= max_k; k += 2) {
    uint32_t mod8;
    if (!ts_squarefree_small(k))
      continue;
    kval[kcount] = k;
    mod8 = (uint32_t)((n & 7U) * k) & 7U;
    scores[kcount] = -0.5 * M_LN2 * log((double)k);
    scores[kcount] += mod8 == 1 ? 2.0 * M_LN2
                       : mod8 == 5 ? M_LN2 : 0.5 * M_LN2;
    if (mod8 == 1)
      scores[kcount] += 3.0 * M_LN2 / 16.0;
    kcount++;
  }
  if (kcount == 0)
    return 1;
  unfinished = kcount;
  while (unfinished != 0) {
    uint32_t nmod;
    int nsymbol;
    p = ts_next_prime_u32(p);
    nmod = (uint32_t)(n % p);
    nsymbol = nmod == 0 ? 0 : ts_jacobi_odd_u32(nmod, p);
    for (i = 0; i < kcount; i++) {
      uint32_t km;
      int match = 0;
      if (accepted[i] == wanted)
        continue;
      km = kval[i] % p;
      if (km == 0) {
        scores[i] += log((double)p) / p;
        match = 1;
      } else if (nsymbol != 0 &&
                 ts_jacobi_odd_u32(km, p) == nsymbol) {
        scores[i] += 2.0 * log((double)p) / (p - 1U);
        match = 1;
      }
      if (match && ++accepted[i] == wanted)
        unfinished--;
    }
  }
  best_score = scores[0];
  for (i = 1; i < kcount; i++)
    if (scores[i] > best_score) {
      best_score = scores[i];
      best = kval[i];
    }
  return best;
}

static int ts_build_factor_base(ts_ctx_t *ctx) {
  uint32_t count = 0, p = 2;
  ctx->fb = (ts_fb_t *)ts_calloc_array(ctx->params.fb_size,
                                       sizeof(*ctx->fb));
  ctx->fb[0].p = 2;
  ctx->fb[0].sqrt_kn = (uint32_t)(ctx->kn & 1U);
  count = 1;
  while (count < ctx->params.fb_size) {
    uint32_t r;
    int symbol;
    p = ts_next_prime_u32(p);
    r = (uint32_t)(ctx->kn % p);
    symbol = ts_legendre_u32(r, p);
    if (symbol < 0)
      continue;
    if (r == 0 && ctx->n % p == 0) {
      ctx->factor = p;
      return 0;
    }
    ctx->fb[count].p = p;
    ctx->fb[count].sqrt_kn = r == 0 ? 0 : ts_sqrtmod_u32(r, p);
    count++;
  }
  ctx->largest_fb_prime = ctx->fb[ctx->params.fb_size - 1U].p;
  return 1;
}

static void ts_hashset_init(ts_hashset_t *set, uint32_t initial) {
  uint32_t alloc = 16;
  while (alloc < initial) {
    if (alloc > UINT32_MAX / 2U)
      croak("tinysiqs128: hash allocation overflow");
    alloc *= 2U;
  }
  set->slots = (uint64_t *)ts_calloc_array(alloc, sizeof(*set->slots));
  set->alloc = alloc;
  set->count = 0;
}

static void ts_hashset_clear(ts_hashset_t *set) {
  free(set->slots);
  memset(set, 0, sizeof(*set));
}

static void ts_hashset_grow(ts_hashset_t *set) {
  uint64_t *old = set->slots;
  uint32_t old_alloc = set->alloc, i;
  if (old_alloc > UINT32_MAX / 2U)
    croak("tinysiqs128: hash is too large");
  set->alloc *= 2U;
  set->slots = (uint64_t *)ts_calloc_array(set->alloc,
                                           sizeof(*set->slots));
  set->count = 0;
  for (i = 0; i < old_alloc; i++) {
    uint64_t value = old[i];
    if (value != TS_HASH_EMPTY) {
      uint32_t pos = (uint32_t)ts_mix64(value) & (set->alloc - 1U);
      while (set->slots[pos] != TS_HASH_EMPTY)
        pos = (pos + 1U) & (set->alloc - 1U);
      set->slots[pos] = value;
      set->count++;
    }
  }
  free(old);
}

static int ts_hashset_insert(ts_hashset_t *set, uint64_t value) {
  uint32_t pos;
  if (value == TS_HASH_EMPTY)
    value = UINT64_C(0x9e3779b97f4a7c15);
  if (set->count * 10U >= set->alloc * 7U)
    ts_hashset_grow(set);
  pos = (uint32_t)ts_mix64(value) & (set->alloc - 1U);
  while (set->slots[pos] != TS_HASH_EMPTY) {
    if (set->slots[pos] == value)
      return 0;
    pos = (pos + 1U) & (set->alloc - 1U);
  }
  set->slots[pos] = value;
  set->count++;
  return 1;
}

/*----------------------------------------------------------------------------
 * Polynomial construction
 *----------------------------------------------------------------------------*/

static int ts_u16_cmp(const void *va, const void *vb) {
  uint16_t a = *(const uint16_t *)va;
  uint16_t b = *(const uint16_t *)vb;
  return a < b ? -1 : a > b ? 1 : 0;
}

static uint64_t ts_target_A(uint128_t kn, uint32_t M, uint32_t d) {
  uint64_t L = (uint64_t)M * d;
  uint64_t S;
  uint128_t a, W;
  uint64_t r;
  if (L < 2 || L > UINT64_MAX / L)
    croak("tinysiqs128: invalid target-A denominator");
  S = L * L;
  a = kn / S;
  r = (uint64_t)(kn % S);
  if (a > TS_U128_MAX / 2U)
    croak("tinysiqs128: target-A quotient overflow");
  W = 2U * a;
  if (r >= S / 2U + (S & 1U))
    W++;
  return ts_isqrt128(W);
}

static uint32_t ts_nearest_fb_index(const ts_ctx_t *ctx, uint32_t wanted) {
  uint32_t lo = 1, hi = ctx->params.fb_size;
  while (lo < hi) {
    uint32_t mid = lo + (hi - lo) / 2U;
    if (ctx->fb[mid].p < wanted)
      lo = mid + 1U;
    else
      hi = mid;
  }
  if (lo >= ctx->params.fb_size)
    return ctx->params.fb_size - 1U;
  if (lo > 1U && wanted - ctx->fb[lo - 1U].p <=
                 ctx->fb[lo].p - wanted)
    return lo - 1U;
  return lo;
}

static uint32_t ts_nearest_available_fb(ts_ctx_t *ctx, uint32_t wanted) {
  uint32_t center = ts_nearest_fb_index(ctx, wanted);
  uint32_t step;
  for (step = 0; step < ctx->params.fb_size; step++) {
    uint32_t up = center + step;
    if (up < ctx->params.fb_size && !ctx->fb[up].in_a &&
        ctx->fb[up].sqrt_kn != 0)
      return up;
    if (step != 0 && center >= step) {
      uint32_t down = center - step;
      if (down > 0 && !ctx->fb[down].in_a &&
          ctx->fb[down].sqrt_kn != 0)
        return down;
    }
  }
  return TS_NO_INDEX;
}

static int ts_choose_A(ts_ctx_t *ctx, ts_poly_t *poly) {
  uint32_t attempt, i;
  uint32_t attempt_limit = poly->q_count <= 3 ? 30000U : 10000U;
  double ideal_d = pow((double)poly->target_A, 1.0 / poly->q_count);
  uint32_t ideal = ideal_d < 3.0 ? 3U : (uint32_t)ideal_d;
  uint32_t center = ts_nearest_fb_index(ctx, ideal);
  uint32_t variance = (uint32_t)(0.75 * sqrt((double)ctx->params.fb_size));
  if (variance < 8U)
    variance = 8U;

  for (i = 0; i < poly->q_count; i++) {
    if (poly->a_index[i] < ctx->params.fb_size)
      ctx->fb[poly->a_index[i]].in_a = 0;
    poly->a_index[i] = UINT16_MAX;
  }

  for (attempt = 0; attempt < attempt_limit; attempt++) {
    uint32_t tolerance = attempt < 10000U ? 2U
                       : attempt < 20000U ? 4U
                       : ctx->params.a_tolerance;
    uint64_t product = 1;
    int acceptable = poly->q_count == 1;

    for (i = 0; i + 1U < poly->q_count; i++) {
      int64_t wanted_index;
      uint32_t search_variance, index;
      if (attempt >= 20000U) {
        wanted_index = 1 + (int64_t)ts_rand_range(
            &ctx->poly_rng, ctx->params.fb_size - 1U);
      } else {
        int64_t offset;
        search_variance = attempt < 10000U ? variance : 2U * variance;
        offset = (int64_t)ts_rand_range(
            &ctx->poly_rng, 2U * search_variance + 1U) - search_variance;
        wanted_index = (int64_t)center + offset;
        if (wanted_index < 1)
          wanted_index = 1;
        if (wanted_index >= (int64_t)ctx->params.fb_size)
          wanted_index = ctx->params.fb_size - 1U;
      }
      index = ts_nearest_available_fb(
          ctx, ctx->fb[(uint32_t)wanted_index].p);
      if (index == TS_NO_INDEX || product > UINT64_MAX / ctx->fb[index].p)
        goto reject_A;
      poly->a_index[i] = (uint16_t)index;
      ctx->fb[index].in_a = 1;
      product *= ctx->fb[index].p;
    }
    {
      uint64_t remaining = product == 0 ? 0 : poly->target_A / product;
      uint32_t wanted = remaining > UINT32_MAX
                      ? UINT32_MAX : (uint32_t)remaining;
      uint32_t index = ts_nearest_available_fb(ctx, wanted);
      if (index == TS_NO_INDEX || product > UINT64_MAX / ctx->fb[index].p)
        goto reject_A;
      poly->a_index[poly->q_count - 1U] = (uint16_t)index;
      ctx->fb[index].in_a = 1;
      product *= ctx->fb[index].p;
    }
    qsort(poly->a_index, poly->q_count, sizeof(poly->a_index[0]),
          ts_u16_cmp);
    if (!acceptable) {
      uint128_t scaled_product = (uint128_t)product * tolerance;
      uint128_t scaled_target = (uint128_t)poly->target_A * tolerance;
      acceptable = scaled_product >= poly->target_A &&
                   product <= scaled_target;
    }
    if (acceptable) {
      /* A is square-free, so its exact product uniquely identifies the
       * selected factor-base set without a fingerprint collision. */
      if (ts_hashset_insert(&ctx->a_hashes, product)) {
        poly->A = product;
        if (product > UINT64_MAX / ctx->params.d)
          goto reject_A;
        poly->D = product * ctx->params.d;
        return 1;
      }
    }

reject_A:
    for (i = 0; i < poly->q_count; i++) {
      if (poly->a_index[i] < ctx->params.fb_size)
        ctx->fb[poly->a_index[i]].in_a = 0;
      poly->a_index[i] = UINT16_MAX;
    }
    if (poly->q_count == 1)
      break;
  }
  return 0;
}

static int ts_compute_C(const ts_ctx_t *ctx, ts_poly_t *poly) {
  uint128_t bmag = ts_abs_i128(poly->B);
  uint128_t b2, difference, magnitude;
  int negative;
  if (bmag != 0 && bmag > TS_U128_MAX / bmag)
    return 0;
  b2 = bmag * bmag;
  if (b2 >= ctx->kn) {
    difference = b2 - ctx->kn;
    negative = 0;
  } else {
    difference = ctx->kn - b2;
    negative = 1;
  }
  if (poly->D == 0 || difference % poly->D != 0)
    return 0;
  magnitude = difference / poly->D;
  if (magnitude > TS_I128_MAX_MAG)
    return 0;
  poly->C = negative ? -(int128_t)magnitude : (int128_t)magnitude;
  return 1;
}

static int ts_check_family_bounds(const ts_ctx_t *ctx,
                                  const ts_poly_t *poly) {
  uint128_t bmax = ctx->params.d == 2 ? poly->A : 0;
  uint128_t dm, twice_b, y_bound, linear_bound, cmax, qbound;
  uint32_t i;
  for (i = 0; i < poly->q_count; i++) {
    if (poly->H[i] > UINT64_MAX / 2U)
      return 0;
    bmax += poly->H[i];
  }
  if (bmax > TS_I128_MAX_MAG ||
      poly->D > TS_U128_MAX / ctx->params.M)
    return 0;
  dm = (uint128_t)poly->D * ctx->params.M;
  if (bmax > TS_U128_MAX / 2U)
    return 0;
  twice_b = 2U * bmax;
  if (dm > TS_U128_MAX - bmax || dm > TS_U128_MAX - twice_b)
    return 0;
  y_bound = dm + bmax;
  linear_bound = dm + twice_b;
  if (bmax != 0 && bmax > TS_U128_MAX / bmax)
    return 0;
  cmax = bmax * bmax;
  if (ctx->kn > cmax)
    cmax = ctx->kn;
  cmax /= poly->D;
  if (linear_bound > TS_U128_MAX / ctx->params.M)
    return 0;
  qbound = linear_bound * ctx->params.M;
  if (qbound > TS_U128_MAX - cmax)
    return 0;
  qbound += cmax;
  return y_bound <= TS_I128_MAX_MAG && qbound <= TS_I128_MAX_MAG;
}

static void ts_set_special_roots(ts_ctx_t *ctx, const ts_poly_t *poly) {
  uint32_t i, M = ctx->params.M;
  if (ctx->params.d == 2) {
    ctx->root1[0] = 0;
    ctx->root2[0] = 1;
  } else {
    ctx->root1[0] = (ts_signed_mod_u32(poly->C, 2) + (M & 1U)) & 1U;
    ctx->root2[0] = TS_NO_ROOT;
  }
  for (i = 0; i < poly->q_count; i++) {
    uint32_t index = poly->a_index[i];
    uint32_t p = ctx->fb[index].p;
    uint32_t b = ts_signed_mod_u32(poly->B, p);
    uint32_t c = ts_signed_mod_u32(poly->C, p);
    uint32_t denom = (uint32_t)((2ULL * b) % p);
    uint32_t inv = ts_inverse_u32(denom, p);
    uint32_t x;
    if (inv == 0)
      croak("tinysiqs128: singular polynomial root");
    x = (uint32_t)((uint64_t)(c ? p - c : 0) * inv % p);
    ctx->root1[index] = (x + M % p) % p;
    ctx->root2[index] = TS_NO_ROOT;
  }
}

#if TS_DEBUG
static uint32_t ts_debug_polynomial_at_root(const ts_ctx_t *ctx,
                                            const ts_poly_t *poly,
                                            uint32_t index,
                                            uint32_t root) {
  uint32_t p = ctx->fb[index].p;
  uint32_t mmod = ctx->params.M % p;
  uint32_t x = root >= mmod ? root - mmod : root + p - mmod;
  uint32_t dmod = (uint32_t)(poly->D % p);
  uint32_t bmod = ts_signed_mod_u32(poly->B, p);
  uint32_t cmod = ts_signed_mod_u32(poly->C, p);
  uint32_t linear = (uint32_t)((2ULL * bmod) % p);
  uint32_t value = (uint32_t)(((uint64_t)dmod * x + linear) % p);
  return (uint32_t)(((uint64_t)value * x + cmod) % p);
}

static void ts_verify_polynomial(const ts_ctx_t *ctx,
                                 const ts_poly_t *poly) {
  uint64_t product = 1;
  uint32_t i;
  for (i = 0; i < poly->q_count; i++) {
    uint32_t index = poly->a_index[i];
    if (index >= ctx->params.fb_size ||
        product > UINT64_MAX / ctx->fb[index].p)
      croak("tinysiqs128: invalid polynomial A factors");
    product *= ctx->fb[index].p;
  }
  if (product != poly->A ||
      poly->D != (uint64_t)ctx->params.d * poly->A)
    croak("tinysiqs128: polynomial A invariant failed");
  if (ctx->params.d == 2 &&
      ((ctx->kn & 7U) != 1U || !(ts_abs_i128(poly->B) & 1U) ||
       ts_signed_mod_u32(poly->C, 4U) != 0))
    croak("tinysiqs128: d=2 polynomial invariant failed");
  for (i = 0; i < ctx->params.fb_size; i++) {
    uint32_t p = ctx->fb[i].p;
    int expect_two = i == 0 ? ctx->params.d == 2
                   : !ctx->fb[i].in_a && ctx->fb[i].sqrt_kn != 0;
    if (ctx->root1[i] >= p ||
        ts_debug_polynomial_at_root(ctx, poly, i, ctx->root1[i]) != 0)
      croak("tinysiqs128: first polynomial root invariant failed");
    if (expect_two) {
      if (ctx->root2[i] >= p || ctx->root2[i] == ctx->root1[i] ||
          ts_debug_polynomial_at_root(ctx, poly, i,
                                      ctx->root2[i]) != 0)
        croak("tinysiqs128: second polynomial root invariant failed");
    } else if (ctx->root2[i] != TS_NO_ROOT) {
      croak("tinysiqs128: unexpected second polynomial root");
    }
  }
}
#endif

static int ts_first_B_and_roots(ts_ctx_t *ctx, ts_poly_t *poly) {
  uint32_t i, j;
  poly->B = 0;
  for (i = 0; i < poly->q_count; i++) {
    uint32_t index = poly->a_index[i];
    uint32_t p = ctx->fb[index].p;
    uint64_t adiv = poly->A / p;
    uint32_t amod = (uint32_t)(adiv % p);
    uint32_t inv = ts_inverse_u32(amod, p);
    uint32_t gamma = (uint32_t)((uint64_t)ctx->fb[index].sqrt_kn * inv % p);
    if (gamma > p / 2U)
      gamma = p - gamma;
    if (gamma != 0 && adiv > UINT64_MAX / gamma)
      return 0;
    poly->H[i] = adiv * gamma;
    poly->B += (int128_t)poly->H[i];
  }
  if (ctx->params.d == 2 && !(ts_abs_i128(poly->B) & 1U))
    poly->B += (int128_t)poly->A;
  if (!ts_check_family_bounds(ctx, poly) || !ts_compute_C(ctx, poly))
    return 0;

  for (j = 1; j < ctx->params.fb_size; j++) {
    uint32_t p = ctx->fb[j].p;
    if (ctx->fb[j].in_a)
      continue;
    {
      uint32_t amod = (uint32_t)(poly->D % p);
      uint32_t inva = ts_inverse_u32(amod, p);
      uint32_t bmod = ts_signed_mod_u32(poly->B, p);
      uint32_t negb = bmod ? p - bmod : 0;
      uint32_t s = ctx->fb[j].sqrt_kn;
      uint32_t x1 = (uint32_t)((uint64_t)((negb + s) % p) * inva % p);
      uint32_t x2 = (uint32_t)((uint64_t)((negb + p - s) % p) * inva % p);
      ctx->root1[j] = (x1 + ctx->params.M % p) % p;
      ctx->root2[j] = x1 == x2 ? TS_NO_ROOT
                    : (x2 + ctx->params.M % p) % p;
      for (i = 0; i + 1U < poly->q_count; i++) {
        uint32_t hmod = (uint32_t)(poly->H[i] % p);
        poly->corrections[(size_t)i * ctx->params.fb_size + j] =
          (uint32_t)((uint64_t)((2ULL * hmod) % p) * inva % p);
      }
    }
  }
  ts_set_special_roots(ctx, poly);
#if TS_DEBUG
  ts_verify_polynomial(ctx, poly);
#endif
  poly->b_index = 0;
  return 1;
}

static void ts_set_family_sieve_initial(ts_ctx_t *ctx) {
  double expected = ctx->params.sieve_free_units;
  uint32_t i, initial;
  if (ctx->params.d == 2)
    expected += ctx->fb[0].logp;
  for (i = 0; i < ctx->params.sieve_start; i++) {
    uint32_t roots = ctx->root2[i] == TS_NO_ROOT ? 1U : 2U;
    expected += roots * (double)ctx->fb[i].logp / ctx->fb[i].p;
  }
  initial = (uint32_t)(expected + 0.5);
  if (initial > 127U)
    initial = 127U;
  ctx->active_sieve_initial = (uint8_t)initial;
}

static int ts_new_family(ts_ctx_t *ctx, ts_poly_t *poly) {
  uint32_t i;
  for (i = 0; i < poly->q_count; i++)
    if (poly->a_index[i] < ctx->params.fb_size)
      ctx->fb[poly->a_index[i]].in_a = 0;
  if (!ts_choose_A(ctx, poly)) {
#if TS_DIAGNOSTIC
    fprintf(stderr, "# tinysiqs128 choose-A exhausted target_A=%llu\n",
            (unsigned long long)poly->target_A);
#endif
    return 0;
  }
  if (!ts_first_B_and_roots(ctx, poly)) {
#if TS_DIAGNOSTIC
    fprintf(stderr, "# tinysiqs128 rejected family A=%llu D=%llu\n",
            (unsigned long long)poly->A, (unsigned long long)poly->D);
#endif
    return 0;
  }
  ts_set_family_sieve_initial(ctx);
  return 1;
}

static void ts_update_roots(ts_ctx_t *ctx, const ts_poly_t *poly,
                            uint32_t bit, int subtract,
                            uint32_t first, uint32_t end) {
  const uint32_t *correction = poly->corrections +
      (size_t)bit * ctx->params.fb_size;
  uint32_t j;
  for (j = first; j < end; j++) {
    uint32_t p = ctx->prime[j], corr = correction[j];
    if (subtract) {
      ctx->root1[j] += corr;
      if (ctx->root1[j] >= p) ctx->root1[j] -= p;
      if (ctx->root2[j] != TS_NO_ROOT) {
        ctx->root2[j] += corr;
        if (ctx->root2[j] >= p) ctx->root2[j] -= p;
      }
    } else {
      ctx->root1[j] = ctx->root1[j] >= corr
                    ? ctx->root1[j] - corr : ctx->root1[j] + p - corr;
      if (ctx->root2[j] != TS_NO_ROOT)
        ctx->root2[j] = ctx->root2[j] >= corr
                      ? ctx->root2[j] - corr : ctx->root2[j] + p - corr;
    }
  }
}

static int ts_next_B(ts_ctx_t *ctx, ts_poly_t *poly) {
  uint32_t bit, gray, i, first;
  int subtract;
  if (poly->b_index + 1U >= poly->b_limit)
    return 0;
  poly->b_index++;
  bit = ts_ctz32(poly->b_index);
  gray = poly->b_index ^ (poly->b_index >> 1);
  subtract = (gray & (1U << bit)) != 0;
  if (subtract)
    poly->B -= (int128_t)(2U * poly->H[bit]);
  else
    poly->B += (int128_t)(2U * poly->H[bit]);
  first = 1;
  for (i = 0; i < poly->q_count; i++) {
    uint32_t a_index = poly->a_index[i];
    ts_update_roots(ctx, poly, bit, subtract, first, a_index);
    first = a_index + 1U;
  }
  ts_update_roots(ctx, poly, bit, subtract, first, ctx->params.fb_size);
  if (!ts_compute_C(ctx, poly))
    croak("tinysiqs128: Gray polynomial is not integral");
  ts_set_special_roots(ctx, poly);
#if TS_DEBUG
  ts_verify_polynomial(ctx, poly);
#endif
  return 1;
}

static int ts_poly_init(ts_ctx_t *ctx, ts_poly_t *poly) {
  size_t correction_count;
  uint32_t i;
  memset(poly, 0, sizeof(*poly));
  poly->q_count = ctx->params.q_count;
  poly->b_limit = 1U << (poly->q_count - 1U);
  for (i = 0; i < TS_MAX_Q; i++)
    poly->a_index[i] = UINT16_MAX;
  correction_count = (size_t)(poly->q_count - 1U) * ctx->params.fb_size;
  poly->corrections = (uint32_t *)ts_calloc_array(
      correction_count, sizeof(*poly->corrections));
  poly->target_A = ts_target_A(ctx->kn, ctx->params.M, ctx->params.d);
  if (poly->target_A == 0) {
    free(poly->corrections);
    memset(poly, 0, sizeof(*poly));
    return 0;
  }
  return 1;
}

static void ts_poly_clear(ts_ctx_t *ctx, ts_poly_t *poly) {
  uint32_t i;
  for (i = 0; i < poly->q_count; i++)
    if (poly->a_index[i] < ctx->params.fb_size)
      ctx->fb[poly->a_index[i]].in_a = 0;
  free(poly->corrections);
  memset(poly, 0, sizeof(*poly));
}

static void ts_set_log_weights(ts_ctx_t *ctx) {
  double log_q = 0.5 * ts_log128(ctx->kn)
               + log((double)ctx->params.M) + 0.5 * M_LN2;
  double nominal = exp(0.15 * (double)ctx->params.bits * M_LN2);
  double bound = 2.0 * (double)ctx->params.lp_multiplier *
                       ctx->largest_fb_prime;
  double denom, scale, expected;
  uint32_t i, initial;
  if (bound < nominal)
    bound = nominal;
  ctx->params.sieve_hit_bound = bound;
  denom = log_q - log(bound);
  if (denom < 8.0)
    denom = 8.0;
  scale = (128.0 - ctx->params.sieve_free_units) / denom;
  ctx->params.log_scale = scale;
  for (i = 0; i < ctx->params.fb_size; i++) {
    uint32_t weight = (uint32_t)(log((double)ctx->fb[i].p) * scale + 0.5);
    if (weight == 0) weight = 1;
    if (weight > 127U) weight = 127U;
    ctx->fb[i].logp = (uint8_t)weight;
  }
  ctx->params.sieve_start = (uint32_t)cbrt((double)ctx->params.fb_size);
  if (ctx->params.sieve_start < 1U)
    ctx->params.sieve_start = 1U;
  if (ctx->params.sieve_start_exponent > 0.0) {
    uint32_t index = (uint32_t)pow((double)ctx->params.fb_size,
                                   ctx->params.sieve_start_exponent);
    if (index > ctx->params.sieve_start)
      ctx->params.sieve_start = index;
  }
  if (ctx->params.sieve_start >= ctx->params.fb_size)
    croak("tinysiqs128: invalid sieve start");
  expected = ctx->params.sieve_free_units;
  for (i = 0; i < ctx->params.sieve_start; i++) {
    double contribution = (double)ctx->fb[i].logp / ctx->fb[i].p;
    expected += contribution *
      (i == 0 ? (double)ctx->params.d
              : ctx->fb[i].sqrt_kn == 0 ? 1.0 : 2.0);
  }
  if (ctx->params.d == 2)
    expected += ctx->fb[0].logp;
  initial = (uint32_t)(expected + 0.5);
  if (initial > 127U) initial = 127U;
  ctx->params.sieve_initial = (uint8_t)initial;
  ctx->active_sieve_initial = (uint8_t)initial;
}

/*----------------------------------------------------------------------------
 * Relations and one-large-prime anchors
 *----------------------------------------------------------------------------*/

static void ts_raw_set(ts_raw_relation_t *relation, uint128_t y,
                       const ts_factor_t *factors,
                       uint32_t nfactors, uint32_t lp) {
  if (nfactors > TS_MAX_RAW_FACTORS)
    croak("tinysiqs128: too many relation factors");
  memcpy(relation->factors, factors,
         (size_t)nfactors * sizeof(*relation->factors));
  relation->y = y;
  relation->nfactors = (uint16_t)nfactors;
  relation->lp = lp;
}

static void ts_store_full(ts_ctx_t *ctx, uint128_t y,
                          const ts_factor_t *factors,
                          uint32_t nfactors) {
  ts_full_relation_t *relation;
  uint32_t factor_alloc;
  if (ctx->full_count >= TS_MAX_RELATIONS)
    croak("tinysiqs128: full relation capacity exceeded");
  if (nfactors > TS_MAX_FULL_FACTORS ||
      nfactors > UINT32_MAX - ctx->full_factor_count)
    croak("tinysiqs128: full-relation factor capacity exceeded");
  if (ctx->full_count == ctx->full_alloc) {
    uint32_t alloc = ctx->full_alloc ? 2U * ctx->full_alloc : 128U;
    if (alloc > TS_MAX_RELATIONS)
      alloc = TS_MAX_RELATIONS;
    ctx->full = (ts_full_relation_t *)ts_realloc_array(
        ctx->full, alloc, sizeof(*ctx->full));
    ctx->full_alloc = alloc;
  }
  if (ctx->full_factor_count + nfactors > ctx->full_factor_alloc) {
    factor_alloc = ctx->full_factor_alloc
                 ? 2U * ctx->full_factor_alloc : 4096U;
    if (factor_alloc < ctx->full_factor_count + nfactors)
      factor_alloc = ctx->full_factor_count + nfactors;
    if (factor_alloc < ctx->full_factor_alloc)
      croak("tinysiqs128: full-factor pool overflow");
    ctx->full_factors = (ts_factor_t *)ts_realloc_array(
        ctx->full_factors, factor_alloc, sizeof(*ctx->full_factors));
    ctx->full_factor_alloc = factor_alloc;
  }
  relation = &ctx->full[ctx->full_count++];
  relation->y = y;
  relation->factor_offset = ctx->full_factor_count;
  relation->nfactors = (uint16_t)nfactors;
  memcpy(ctx->full_factors + ctx->full_factor_count, factors,
         (size_t)nfactors * sizeof(*factors));
  ctx->full_factor_count += nfactors;
}

static void ts_materialize_smooth(ts_ctx_t *ctx,
                                  const ts_raw_relation_t *raw) {
  ts_store_full(ctx, raw->y, raw->factors, raw->nfactors);
  ctx->accepted_smooth++;
}

static uint16_t ts_merge_relation_factors(
    const ts_raw_relation_t *left, const ts_raw_relation_t *right,
    ts_factor_t merged[TS_MAX_FULL_FACTORS]) {
  uint32_t li = 0, ri = 0, count = 0, unique = 0;
  while (li < left->nfactors || ri < right->nfactors) {
    if (ri == right->nfactors ||
        (li < left->nfactors &&
         left->factors[li].row < right->factors[ri].row))
      li++;
    else if (li == left->nfactors ||
             right->factors[ri].row < left->factors[li].row)
      ri++;
    else {
      li++;
      ri++;
    }
    unique++;
  }
  if (unique > TS_MAX_FULL_FACTORS)
    croak("tinysiqs128: merged relation is too large");
  li = ri = 0;
  while (li < left->nfactors || ri < right->nfactors) {
    if (ri == right->nfactors ||
        (li < left->nfactors &&
         left->factors[li].row < right->factors[ri].row)) {
      merged[count++] = left->factors[li++];
    } else if (li == left->nfactors ||
               right->factors[ri].row < left->factors[li].row) {
      merged[count++] = right->factors[ri++];
    } else {
      uint32_t exponent = (uint32_t)left->factors[li].exponent +
                          right->factors[ri].exponent;
      if (exponent > UINT16_MAX)
        croak("tinysiqs128: relation exponent overflow");
      merged[count].row = left->factors[li].row;
      merged[count].exponent = (uint16_t)exponent;
      count++;
      li++;
      ri++;
    }
  }
  return (uint16_t)count;
}

static void ts_lp_init(ts_lp_state_t *state, uint32_t initial) {
  uint32_t alloc = 64;
  while (alloc < initial)
    alloc *= 2U;
  state->anchors = (ts_lp_anchor_t *)ts_calloc_array(
      alloc, sizeof(*state->anchors));
  state->alloc = alloc;
  state->count = 0;
}

static void ts_lp_clear(ts_lp_state_t *state) {
  free(state->anchors);
  memset(state, 0, sizeof(*state));
}

static uint32_t ts_lp_hash(uint32_t label, uint32_t mask) {
  return (uint32_t)ts_mix64(label) & mask;
}

static void ts_lp_grow(ts_lp_state_t *state) {
  ts_lp_anchor_t *old = state->anchors;
  uint32_t old_alloc = state->alloc, i;
  if (old_alloc > UINT32_MAX / 2U)
    croak("tinysiqs128: LP table is too large");
  state->alloc *= 2U;
  state->anchors = (ts_lp_anchor_t *)ts_calloc_array(
      state->alloc, sizeof(*state->anchors));
  state->count = 0;
  for (i = 0; i < old_alloc; i++) {
    if (old[i].label != 0) {
      uint32_t pos = ts_lp_hash(old[i].label, state->alloc - 1U);
      while (state->anchors[pos].label != 0)
        pos = (pos + 1U) & (state->alloc - 1U);
      state->anchors[pos] = old[i];
      state->count++;
    }
  }
  free(old);
}

static ts_raw_relation_t *ts_lp_anchor(ts_lp_state_t *state,
                                       uint32_t label,
                                       const ts_raw_relation_t *relation) {
  uint32_t pos;
  if (label <= 1)
    croak("tinysiqs128: invalid LP label");
  if (state->count * 10U >= state->alloc * 7U)
    ts_lp_grow(state);
  pos = ts_lp_hash(label, state->alloc - 1U);
  while (state->anchors[pos].label != 0) {
    if (state->anchors[pos].label == label)
      return &state->anchors[pos].relation;
    pos = (pos + 1U) & (state->alloc - 1U);
  }
  state->anchors[pos].label = label;
  state->anchors[pos].relation.y = relation->y;
  state->anchors[pos].relation.nfactors = relation->nfactors;
  state->anchors[pos].relation.lp = relation->lp;
  memcpy(state->anchors[pos].relation.factors, relation->factors,
         (size_t)relation->nfactors * sizeof(*relation->factors));
  state->count++;
  return NULL;
}

static void ts_ensure_mont(ts_ctx_t *ctx) {
  if (!ctx->mont_ready) {
    ts_mont_setup(&ctx->mont, ctx->n);
    ctx->mont_ready = 1;
  }
}

#if TS_DEBUG
static void ts_verify_relation(ts_ctx_t *ctx, uint128_t y,
                               const ts_factor_t *factors,
                               uint32_t nfactors, uint32_t lp) {
  uint128_t lhs, rhs;
  uint32_t i;
  ts_ensure_mont(ctx);
  lhs = ts_mont_enter(y, &ctx->mont);
  lhs = ts_mont_mul(lhs, lhs, &ctx->mont);
  rhs = ctx->mont.one;
  for (i = 0; i < nfactors; i++) {
    uint32_t row = factors[i].row;
    uint32_t exponent = factors[i].exponent;
    if (row > ctx->params.fb_size || exponent == 0)
      croak("tinysiqs128: invalid relation factor");
    if (row == 0) {
      if (exponent & 1U)
        rhs = rhs == 0 ? 0 : ctx->n - rhs;
    } else {
      uint128_t base = ts_mont_enter(ctx->fb[row - 1U].p,
                                    &ctx->mont);
      rhs = ts_mont_mul(rhs,
                        ts_mont_pow_u32(base, exponent, &ctx->mont),
                        &ctx->mont);
    }
  }
  if (lp != 0) {
    uint64_t pmax = ctx->largest_fb_prime;
    uint64_t bound = (uint64_t)ctx->params.lp_multiplier * pmax;
    if (lp <= pmax || lp > bound ||
        (uint64_t)lp >= pmax * pmax || !ts_is_prime_u32(lp))
      croak("tinysiqs128: invalid large-prime residual");
    rhs = ts_mont_mul(rhs, ts_mont_enter(lp, &ctx->mont), &ctx->mont);
  }
  if (lhs != rhs)
    croak("tinysiqs128: relation identity failed");
}
#endif

static int ts_materialize_lp_pair(ts_ctx_t *ctx,
                                  const ts_raw_relation_t *anchor,
                                  const ts_raw_relation_t *relation,
                                  uint32_t label) {
  ts_factor_t factors[TS_MAX_FULL_FACTORS];
  uint16_t nfactors;
  uint128_t inverse, y;
  ts_ensure_mont(ctx);
  inverse = ts_modinv128(label, &ctx->mont);
  if (inverse == 0) {
    uint128_t factor = ts_gcd128(label, ctx->n);
    if (factor > 1 && factor < ctx->n)
      ctx->factor = factor;
    return 0;
  }
  y = ts_mont_mul(ts_mont_enter(anchor->y, &ctx->mont),
                  ts_mont_enter(relation->y, &ctx->mont), &ctx->mont);
  y = ts_mont_mul(y, ts_mont_enter(inverse, &ctx->mont), &ctx->mont);
  y = ts_mont_exit(y, &ctx->mont);
  nfactors = ts_merge_relation_factors(anchor, relation, factors);
#if TS_DEBUG
  ts_verify_relation(ctx, y, factors, nfactors, 0);
#endif
  ts_store_full(ctx, y, factors, nfactors);
  return 1;
}

static void ts_accept_relation(ts_ctx_t *ctx,
                               const ts_raw_relation_t *relation) {
#if TS_DEBUG
  ts_verify_relation(ctx, relation->y, relation->factors,
                     relation->nfactors, relation->lp);
#endif
  if (relation->lp == 0) {
    ts_materialize_smooth(ctx, relation);
  } else {
    ts_raw_relation_t *anchor = ts_lp_anchor(
        &ctx->one_lp, relation->lp, relation);
    ctx->accepted_partial++;
    if (anchor != NULL) {
      (void)ts_materialize_lp_pair(ctx, anchor, relation, relation->lp);
    }
  }
}

/*----------------------------------------------------------------------------
 * Byte sieve, candidate filter, resieve, and relation evaluation
 *----------------------------------------------------------------------------*/

static void ts_sieve_add(uint8_t *cell, uint8_t logp) {
#if TS_DEBUG
  uint32_t value = (uint32_t)*cell + logp;
  if (value > UINT8_MAX)
    croak("tinysiqs128: byte sieve overflow");
  *cell = (uint8_t)value;
#else
  *cell += logp;
#endif
}

static void ts_sieve_one_root(uint8_t *sieve, uint32_t length,
                              uint32_t root, uint32_t p, uint8_t logp) {
  uint32_t pos = root;
  for (; pos + 4U * p < length; pos += 4U * p) {
    ts_sieve_add(sieve + pos, logp);
    ts_sieve_add(sieve + pos + p, logp);
    ts_sieve_add(sieve + pos + 2U * p, logp);
    ts_sieve_add(sieve + pos + 3U * p, logp);
  }
  for (; pos < length; pos += p)
    ts_sieve_add(sieve + pos, logp);
}

static void ts_sieve_two_roots(uint8_t *sieve, uint32_t length,
                               uint32_t root1, uint32_t root2,
                               uint32_t p, uint8_t logp) {
  uint32_t pos, gap1, gap2, p4 = 4U * p;
  if (root1 > root2) {
    uint32_t t = root1;
    root1 = root2;
    root2 = t;
  }
  pos = root1;
  gap1 = root2 - root1;
  gap2 = p - gap1;
  while (pos + p4 < length) {
    ts_sieve_add(sieve + pos, logp); pos += gap1;
    ts_sieve_add(sieve + pos, logp); pos += gap2;
    ts_sieve_add(sieve + pos, logp); pos += gap1;
    ts_sieve_add(sieve + pos, logp); pos += gap2;
    ts_sieve_add(sieve + pos, logp); pos += gap1;
    ts_sieve_add(sieve + pos, logp); pos += gap2;
    ts_sieve_add(sieve + pos, logp); pos += gap1;
    ts_sieve_add(sieve + pos, logp); pos += gap2;
  }
  while (pos < length) {
    ts_sieve_add(sieve + pos, logp);
    pos += gap1;
    if (pos >= length)
      break;
    ts_sieve_add(sieve + pos, logp);
    pos += gap2;
  }
}

static void ts_sieve_large(uint8_t *sieve, uint32_t length,
                           uint32_t root1, uint32_t root2,
                           uint32_t p, uint8_t logp, uint32_t count) {
  uint32_t i, pos = root1;
  for (i = 0; i < count; i++, pos += p) {
#if TS_DEBUG
    if (pos < length)
#endif
      ts_sieve_add(sieve + pos, logp);
  }
  if (root2 != TS_NO_ROOT) {
    pos = root2;
    for (i = 0; i < count; i++, pos += p) {
#if TS_DEBUG
      if (pos < length)
#endif
        ts_sieve_add(sieve + pos, logp);
    }
  }
#if !TS_DEBUG
  (void)length;
#endif
}

static void ts_run_sieve(ts_ctx_t *ctx) {
  uint32_t i, length = ctx->sieve_length;
  uint32_t initial = ctx->active_sieve_initial + ctx->params.stage1_bias;
  if (initial > UINT8_MAX)
    croak("tinysiqs128: sieve initializer overflow");
  memset(ctx->sieve, (uint8_t)initial, length);
  for (i = ctx->params.sieve_start;
       i < ctx->params.fb_size && ctx->prime[i] <= length / 6U; i++) {
    if (ctx->root2[i] == TS_NO_ROOT)
      ts_sieve_one_root(ctx->sieve, length, ctx->root1[i],
                        ctx->prime[i], ctx->fb[i].logp);
    else
      ts_sieve_two_roots(ctx->sieve, length, ctx->root1[i], ctx->root2[i],
                         ctx->prime[i], ctx->fb[i].logp);
  }
#define TS_LARGE_RANGE(bound, count) do {                                  \
    for (; i < ctx->params.fb_size && ctx->prime[i] <= (bound); i++)       \
      ts_sieve_large(ctx->sieve, length, ctx->root1[i], ctx->root2[i],     \
                     ctx->prime[i], ctx->fb[i].logp, (count));             \
  } while (0)
  TS_LARGE_RANGE(length / 5U, 6U);
  TS_LARGE_RANGE(length / 4U, 5U);
  TS_LARGE_RANGE(length / 3U, 4U);
  TS_LARGE_RANGE(length / 2U, 3U);
  TS_LARGE_RANGE(length, 2U);
  TS_LARGE_RANGE(UINT32_MAX, 1U);
#undef TS_LARGE_RANGE
}

static void ts_clear_candidates(ts_ctx_t *ctx) {
  uint32_t i;
  for (i = 0; i < ctx->candidate_count; i++) {
    uint32_t pos = (uint32_t)(ctx->candidates[i].x +
                              (int32_t)ctx->params.M);
    ctx->candidate_at[pos] = 0;
    if (ctx->candidate_wide)
      ctx->candidate_at_wide[pos] = 0;
  }
  ctx->candidate_wide = 0;
  ctx->candidate_count = 0;
  ctx->hit_count = 0;
}

static void ts_add_candidate(ts_ctx_t *ctx, uint32_t pos) {
  uint32_t physical = ctx->sieve[pos];
  uint32_t logical;
  if (physical < 128U || physical < ctx->params.stage1_bias)
    croak("tinysiqs128: invalid candidate score");
  logical = physical - ctx->params.stage1_bias;
  if (ctx->candidate_count == UINT16_MAX) {
    uint32_t i;
    if (ctx->candidate_at_wide == NULL)
      ctx->candidate_at_wide = (uint32_t *)ts_calloc_array(
          ctx->sieve_length, sizeof(*ctx->candidate_at_wide));
    for (i = 0; i < ctx->candidate_count; i++) {
      uint32_t old_pos = (uint32_t)(ctx->candidates[i].x +
                                    (int32_t)ctx->params.M);
      ctx->candidate_at_wide[old_pos] = i + 1U;
    }
    ctx->candidate_wide = 1;
  }
  if (ctx->candidate_count == ctx->candidate_alloc) {
    uint32_t alloc = ctx->candidate_alloc ? 2U * ctx->candidate_alloc : 1024U;
    if (alloc < ctx->candidate_alloc || alloc > ctx->sieve_length)
      alloc = ctx->sieve_length;
    if (alloc == ctx->candidate_alloc)
      croak("tinysiqs128: candidate capacity exceeded");
    ctx->candidates = (ts_candidate_t *)ts_realloc_array(
        ctx->candidates, alloc, sizeof(*ctx->candidates));
    ctx->candidate_alloc = alloc;
  }
  ctx->candidates[ctx->candidate_count].x =
      (int32_t)pos - (int32_t)ctx->params.M;
  ctx->candidates[ctx->candidate_count].first_hit = TS_NO_INDEX;
  ctx->candidates[ctx->candidate_count].sieve_score = (uint8_t)logical;
  if (ctx->candidate_wide)
    ctx->candidate_at_wide[pos] = ctx->candidate_count + 1U;
  else
    ctx->candidate_at[pos] = (uint16_t)(ctx->candidate_count + 1U);
  ctx->candidate_count++;
}

static void ts_find_candidates(ts_ctx_t *ctx) {
  static const uint64_t high_bits = UINT64_C(0x8080808080808080);
  uint32_t pos, b, bulk = ctx->sieve_length & ~31U;
  ts_clear_candidates(ctx);
  for (pos = 0; pos < bulk; pos += 32U) {
    uint64_t w0, w1, w2, w3;
    memcpy(&w0, ctx->sieve + pos, sizeof(w0));
    memcpy(&w1, ctx->sieve + pos + 8U, sizeof(w1));
    memcpy(&w2, ctx->sieve + pos + 16U, sizeof(w2));
    memcpy(&w3, ctx->sieve + pos + 24U, sizeof(w3));
    if (!((w0 | w1 | w2 | w3) & high_bits))
      continue;
    for (b = 0; b < 32U; b++)
      if (ctx->sieve[pos + b] & 0x80U)
        ts_add_candidate(ctx, pos + b);
  }
  for (; pos < ctx->sieve_length; pos++)
    if (ctx->sieve[pos] & 0x80U)
      ts_add_candidate(ctx, pos);
  ctx->total_candidates += ctx->candidate_count;
}

static void ts_add_hit(ts_ctx_t *ctx, uint32_t candidate,
                       uint32_t fb_index) {
  ts_candidate_t *c = &ctx->candidates[candidate];
  if (fb_index > UINT16_MAX)
    croak("tinysiqs128: factor-base index overflow");
  if (ctx->hit_count == ctx->hit_alloc) {
    uint32_t alloc = ctx->hit_alloc ? 2U * ctx->hit_alloc : 8192U;
    if (alloc < ctx->hit_alloc)
      croak("tinysiqs128: hit capacity overflow");
    ctx->hits = (ts_hit_t *)ts_realloc_array(
        ctx->hits, alloc, sizeof(*ctx->hits));
    ctx->hit_alloc = alloc;
  }
  ctx->hits[ctx->hit_count].fb_index = (uint16_t)fb_index;
  ctx->hits[ctx->hit_count].next = c->first_hit;
  c->first_hit = ctx->hit_count++;
}

static uint32_t ts_reduce_u32(uint32_t n, uint32_t p,
                              uint32_t reciprocal) {
  uint32_t q = (uint32_t)(((uint64_t)n * reciprocal) >> 32);
  uint32_t r = n - q * p;
  return r >= p ? r - p : r;
}

static int128_t ts_eval_q(const ts_poly_t *poly, int32_t x,
                          int128_t *y_out) {
  int128_t sx = (int128_t)x;
  int128_t y = (int128_t)poly->D * sx + poly->B;
  int128_t q = (y + poly->B) * sx + poly->C;
  if (y_out != NULL)
    *y_out = y;
  return q;
}

static uint32_t ts_remove_u128(uint128_t *value, uint32_t p) {
  uint32_t exponent = 0;
  if (*value == 0)
    croak("tinysiqs128: valuation of zero");
  while (*value % p == 0) {
    *value /= p;
    exponent++;
  }
  return exponent;
}

static void ts_filter_candidates(ts_ctx_t *ctx, const ts_poly_t *poly) {
  uint32_t original = ctx->candidate_count, read, out = 0;
  uint16_t small_hits[TS_POSTFILTER_MAX_SMALL];
  double log_smooth = log((double)ctx->params.lp_multiplier *
                          ctx->largest_fb_prime);
  if (ctx->params.sieve_start > TS_POSTFILTER_MAX_SMALL)
    croak("tinysiqs128: too many postfilter primes");
  for (read = 0; read < original; read++) {
    ts_candidate_t candidate = ctx->candidates[read];
    uint32_t pos = (uint32_t)(candidate.x + (int32_t)ctx->params.M);
    uint32_t small_count = 0, i;
    int128_t signed_q = ts_eval_q(poly, candidate.x, NULL);
    int pass;
    ctx->candidate_at[pos] = 0;
    if (ctx->candidate_wide)
      ctx->candidate_at_wide[pos] = 0;
    if (signed_q == 0) {
      pass = 1;
    } else {
      uint128_t q = ts_abs_i128(signed_q);
      for (i = 0; i < ctx->params.sieve_start; i++) {
        uint32_t rem = ts_reduce_u32(pos, ctx->prime[i], ctx->reciprocal[i]);
        if (rem == ctx->root1[i] || rem == ctx->root2[i]) {
          if (ts_remove_u128(&q, ctx->prime[i]) == 0)
            croak("tinysiqs128: postfilter root is not a factor");
          small_hits[small_count++] = (uint16_t)i;
        }
      }
      pass = (double)candidate.sieve_score - ctx->active_sieve_initial -
             ctx->params.log_scale * (ts_log128(q) - log_smooth) >= 0.0;
    }
    if (pass) {
      ctx->candidates[out] = candidate;
      ctx->candidates[out].first_hit = TS_NO_INDEX;
      if (ctx->candidate_wide)
        ctx->candidate_at_wide[pos] = out + 1U;
      else
        ctx->candidate_at[pos] = (uint16_t)(out + 1U);
      for (i = 0; i < small_count; i++)
        ts_add_hit(ctx, out, small_hits[i]);
      out++;
    }
  }
  ctx->candidate_count = out;
}

static void ts_resieve_root16(ts_ctx_t *ctx, uint32_t fb_index,
                              uint32_t root, uint32_t p) {
  uint32_t pos;
  for (pos = root; pos < ctx->sieve_length; pos += p) {
    uint16_t candidate = ctx->candidate_at[pos];
    if (candidate != 0)
      ts_add_hit(ctx, candidate - 1U, fb_index);
  }
}

static void ts_resieve_root32(ts_ctx_t *ctx, uint32_t fb_index,
                              uint32_t root, uint32_t p) {
  uint32_t pos;
  for (pos = root; pos < ctx->sieve_length; pos += p) {
    uint32_t candidate = ctx->candidate_at_wide[pos];
    if (candidate != 0)
      ts_add_hit(ctx, candidate - 1U, fb_index);
  }
}

static void ts_resieve_candidates(ts_ctx_t *ctx) {
  uint32_t begin = ctx->params.sieve_start, cutoff = begin;
  uint32_t cutoff_prime, barrett_end, i, c;
  if (ctx->candidate_count == 0)
    return;
  {
    uint64_t proposed = (uint64_t)TS_RESIEVE_COEFFICIENT *
                        ctx->sieve_length / ctx->candidate_count;
    cutoff_prime = proposed > UINT32_MAX ? UINT32_MAX : (uint32_t)proposed;
  }
  if (cutoff_prime < TS_RESIEVE_FLOOR)
    cutoff_prime = TS_RESIEVE_FLOOR;
  while (cutoff < ctx->params.fb_size &&
         ctx->prime[cutoff] <= cutoff_prime)
    cutoff++;
  barrett_end = cutoff < ctx->resieve_one_subtract_index
              ? cutoff : ctx->resieve_one_subtract_index;
  for (c = 0; c < ctx->candidate_count; c++) {
    uint32_t pos = (uint32_t)(ctx->candidates[c].x +
                              (int32_t)ctx->params.M);
    for (i = begin; i < barrett_end; i++) {
      uint32_t rem = ts_reduce_u32(pos, ctx->prime[i], ctx->reciprocal[i]);
      if (rem == ctx->root1[i] || rem == ctx->root2[i])
        ts_add_hit(ctx, c, i);
    }
    for (; i < cutoff; i++) {
      uint32_t p = ctx->prime[i];
      uint32_t rem = pos >= p ? pos - p : pos;
      if (rem == ctx->root1[i] || rem == ctx->root2[i])
        ts_add_hit(ctx, c, i);
    }
  }
  for (i = cutoff; i < ctx->params.fb_size; i++) {
    if (ctx->candidate_wide) {
      ts_resieve_root32(ctx, i, ctx->root1[i], ctx->prime[i]);
      if (ctx->root2[i] != TS_NO_ROOT)
        ts_resieve_root32(ctx, i, ctx->root2[i], ctx->prime[i]);
    } else {
      ts_resieve_root16(ctx, i, ctx->root1[i], ctx->prime[i]);
      if (ctx->root2[i] != TS_NO_ROOT)
        ts_resieve_root16(ctx, i, ctx->root2[i], ctx->prime[i]);
    }
  }
}

static uint32_t ts_merge_candidate_factors(
    ts_factor_t *factors, uint32_t hit_offset, uint32_t hit_count,
    const ts_factor_t *extra, uint32_t extra_count) {
  uint32_t hit = 0, add = 0, out = 0;
  while (hit < hit_count || add < extra_count) {
    ts_factor_t next;
    if (hit == hit_count ||
        (add < extra_count &&
         extra[add].row < factors[hit_offset + hit].row)) {
      next = extra[add++];
    } else if (add == extra_count ||
               factors[hit_offset + hit].row < extra[add].row) {
      next = factors[hit_offset + hit++];
    } else {
      uint32_t exponent = (uint32_t)extra[add].exponent +
                          factors[hit_offset + hit].exponent;
      if (exponent > UINT16_MAX)
        croak("tinysiqs128: candidate exponent overflow");
      next = extra[add++];
      next.exponent = (uint16_t)exponent;
      hit++;
    }
    factors[out++] = next;
  }
  return out;
}

static void ts_evaluate_candidate(ts_ctx_t *ctx, const ts_poly_t *poly,
                                  uint32_t candidate_index) {
  const ts_candidate_t *candidate = &ctx->candidates[candidate_index];
  ts_factor_t extra[TS_MAX_Q + 2U];
  ts_factor_t factors[TS_MAX_RAW_FACTORS];
  ts_raw_relation_t relation;
  uint32_t hit_count = 0, hit, required, hit_offset;
  uint32_t extra_count = 0, count, lp = 0;
  int128_t signed_y, signed_q = ts_eval_q(poly, candidate->x, &signed_y);
  uint128_t rest;
  if (signed_q == 0) {
    uint128_t factor = ts_gcd128(ts_signed_mod_n(signed_y, ctx->n), ctx->n);
    if (factor > 1 && factor < ctx->n)
      ctx->factor = factor;
    return;
  }
  for (hit = candidate->first_hit; hit != TS_NO_INDEX;
       hit = ctx->hits[hit].next)
    hit_count++;
  hit_offset = poly->q_count + 2U;
  required = hit_offset + hit_count;
  if (required > TS_MAX_RAW_FACTORS)
    croak("tinysiqs128: candidate factor capacity exceeded");
  rest = ts_abs_i128(signed_q);
  if (signed_q < 0) {
    extra[extra_count].row = 0;
    extra[extra_count].exponent = 1;
    extra_count++;
  }
  {
    uint32_t write = hit_offset + hit_count;
    for (hit = candidate->first_hit; hit != TS_NO_INDEX;
         hit = ctx->hits[hit].next) {
      uint32_t fb_index = ctx->hits[hit].fb_index;
      uint32_t exponent = ts_remove_u128(&rest, ctx->prime[fb_index]);
      if (exponent == 0 || exponent > UINT16_MAX)
        croak("tinysiqs128: resieve hit is not a factor");
      write--;
      factors[write].row = (uint16_t)(fb_index + 1U);
      factors[write].exponent = (uint16_t)exponent;
    }
    if (write != hit_offset)
      croak("tinysiqs128: hit count changed");
  }
  if (ctx->params.d == 2) {
    extra[extra_count].row = 1;
    extra[extra_count].exponent = 1;
    extra_count++;
  }
  {
    uint32_t i;
    for (i = 0; i < poly->q_count; i++) {
      extra[extra_count].row = (uint16_t)(poly->a_index[i] + 1U);
      extra[extra_count].exponent = 1;
      extra_count++;
    }
  }
  count = ts_merge_candidate_factors(
      factors, hit_offset, hit_count, extra, extra_count);
  if (rest != 1) {
    uint64_t pmax = ctx->largest_fb_prime;
    uint64_t bound = (uint64_t)ctx->params.lp_multiplier * pmax;
    if (rest > UINT32_MAX || rest <= pmax || rest > bound ||
        rest >= (uint64_t)pmax * pmax)
      return;
    lp = (uint32_t)rest;
  }
  ts_raw_set(&relation, ts_signed_mod_n(signed_y, ctx->n),
             factors, count, lp);
  ts_accept_relation(ctx, &relation);
}

static int ts_solve(ts_ctx_t *ctx);

static void ts_sieve_polynomial(ts_ctx_t *ctx, const ts_poly_t *poly) {
  uint32_t i;
  ts_run_sieve(ctx);
  ts_find_candidates(ctx);
  ts_filter_candidates(ctx, poly);
  ts_resieve_candidates(ctx);
  for (i = 0; i < ctx->candidate_count && ctx->factor == 0 &&
              ctx->full_count < TS_MAX_RELATIONS; i++) {
    ts_evaluate_candidate(ctx, poly, i);
    if (ctx->factor == 0 &&
        ctx->full_count >= ctx->matrix_next_target &&
        ctx->full_count > ctx->matrix_last_count) {
      ctx->matrix_last_count = ctx->full_count;
      if (ts_solve(ctx))
        break;
      if (ctx->full_count >
          ctx->matrix_target_limit - ctx->matrix_retry_batch) {
        ctx->matrix_next_target = ctx->matrix_target_limit + 1U;
        break;
      }
      ctx->matrix_next_target =
          ctx->full_count + ctx->matrix_retry_batch;
    }
  }
  ts_clear_candidates(ctx);
}

/*----------------------------------------------------------------------------
 * Exact dense linear algebra and square-root extraction
 *----------------------------------------------------------------------------*/

#define TS_BIT64(bit) (UINT64_C(1) << (bit))

static int ts_column_cmp(const void *va, const void *vb) {
  const ts_la_col_t *a = (const ts_la_col_t *)va;
  const ts_la_col_t *b = (const ts_la_col_t *)vb;
  return a->weight < b->weight ? -1 : a->weight > b->weight ? 1 : 0;
}

static ts_la_col_t *ts_build_matrix(const ts_ctx_t *ctx,
                                    unsigned long *nrows,
                                    unsigned long *ncols) {
  ts_la_col_t *columns;
  uint32_t i, j;
  *nrows = (unsigned long)ctx->params.fb_size + 1UL;
  *ncols = (unsigned long)ctx->full_count;
  columns = (ts_la_col_t *)ts_calloc_array(ctx->full_count,
                                           sizeof(*columns));
  for (i = 0; i < ctx->full_count; i++) {
    const ts_full_relation_t *relation = &ctx->full[i];
    const ts_factor_t *factors =
        ctx->full_factors + relation->factor_offset;
    unsigned long weight = 0;
    for (j = 0; j < relation->nfactors; j++) {
      if (factors[j].row >= *nrows)
        croak("tinysiqs128: relation matrix row is out of range");
      if (factors[j].exponent & 1U)
        weight++;
    }
    columns[i].data = (unsigned long *)ts_malloc_array(
        weight, sizeof(*columns[i].data));
    columns[i].weight = 0;
    columns[i].orig = i;
    for (j = 0; j < relation->nfactors; j++)
      if (factors[j].exponent & 1U)
        columns[i].data[columns[i].weight++] =
            factors[j].row;
  }
  qsort(columns, ctx->full_count, sizeof(*columns), ts_column_cmp);
  return columns;
}

/* Remove singleton-bearing columns to convergence.  Keep every surviving
 * column: the dense solver enumerates the complete basis in 64-wide batches,
 * including zero and duplicate columns that can be immediately useful. */
static void ts_reduce_matrix(unsigned long nrows, unsigned long *ncols,
                             ts_la_col_t *columns) {
  unsigned long previous_rows, previous_cols, reduced_rows, reduced_cols;
  unsigned long *counts;
  unsigned long c, i, j, k;
  counts = (unsigned long *)ts_calloc_array(nrows, sizeof(*counts));
  for (i = 0; i < *ncols; i++)
    for (j = 0; j < columns[i].weight; j++)
      counts[columns[i].data[j]]++;
  reduced_rows = nrows;
  reduced_cols = *ncols;
  do {
    previous_rows = reduced_rows;
    previous_cols = reduced_cols;
    do {
      c = reduced_cols;
      for (i = j = 0; i < reduced_cols; i++) {
        ts_la_col_t *column = columns + i;
        for (k = 0; k < column->weight; k++)
          if (counts[column->data[k]] < 2UL)
            break;
        if (k < column->weight) {
          for (k = 0; k < column->weight; k++)
            counts[column->data[k]]--;
          free(column->data);
          column->data = NULL;
        } else {
          if (j != i) {
            columns[j] = columns[i];
            columns[i].data = NULL;
          }
          j++;
        }
      }
      reduced_cols = j;
    } while (c != reduced_cols);
    reduced_rows = 0;
    for (i = 0; i < nrows; i++)
      if (counts[i] != 0)
        reduced_rows++;
  } while (previous_rows != reduced_rows ||
           previous_cols != reduced_cols);
  free(counts);
  *ncols = reduced_cols;
}

/* Return one successive batch of at most 64 exact nullspace basis vectors.
 * skip_dependencies makes later calls enumerate genuinely new free columns. */
static uint64_t *ts_dense_nullspace_batch(
    unsigned long nrows, unsigned long ncols, const ts_la_col_t *columns,
    unsigned long skip_dependencies, uint64_t *mask,
    unsigned long *total_dependencies) {
  const unsigned long no_row = ~0UL;
  uint8_t *row_active;
  unsigned long *row_map, *pivot_row;
  uint64_t *matrix, *used_rows, *result, *check;
  unsigned long active_rows = 0, row_words, dependencies = 0, packed = 0;
  unsigned long c, i, j, word;
  size_t matrix_words;

  *mask = 0;
  *total_dependencies = 0;
  if (ncols == 0)
    return NULL;
  row_active = (uint8_t *)ts_calloc_array(nrows, sizeof(*row_active));
  row_map = (unsigned long *)ts_malloc_array(nrows, sizeof(*row_map));
  pivot_row = (unsigned long *)ts_malloc_array(ncols, sizeof(*pivot_row));
  result = (uint64_t *)ts_calloc_array(ncols, sizeof(*result));

  for (c = 0; c < ncols; c++)
    for (i = 0; i < columns[c].weight; i++)
      row_active[columns[c].data[i]] = 1;
  for (i = 0; i < nrows; i++)
    if (row_active[i])
      row_map[i] = active_rows++;
  row_words = (active_rows + 63UL) / 64UL;
  matrix_words = ts_array_bytes((size_t)ncols, (size_t)row_words);
  matrix = (uint64_t *)ts_calloc_array(
      matrix_words, sizeof(*matrix));
  used_rows = (uint64_t *)ts_calloc_array(row_words,
                                           sizeof(*used_rows));

  for (c = 0; c < ncols; c++) {
    uint64_t *dense = matrix + (size_t)c * row_words;
    pivot_row[c] = no_row;
    for (i = 0; i < columns[c].weight; i++) {
      unsigned long row = row_map[columns[c].data[i]];
      dense[row >> 6] ^= TS_BIT64(row & 63UL);
    }
  }
  for (c = 0; c < ncols; c++) {
    uint64_t *pivot_column = matrix + (size_t)c * row_words;
    unsigned long pivot = no_row, pivot_word;
    uint64_t pivot_bit;
    for (word = 0; word < row_words; word++) {
      uint64_t available = pivot_column[word] & ~used_rows[word];
      if (available != 0) {
        pivot = word * 64UL + ts_ctz64(available);
        break;
      }
    }
    if (pivot == no_row)
      continue;
    pivot_row[c] = pivot;
    pivot_word = pivot >> 6;
    pivot_bit = TS_BIT64(pivot & 63UL);
    used_rows[pivot_word] |= pivot_bit;
    pivot_column[pivot_word] &= ~pivot_bit;
    for (i = c + 1UL; i < ncols; i++) {
      uint64_t *other = matrix + (size_t)i * row_words;
      if (other[pivot_word] & pivot_bit)
        for (j = 0; j < row_words; j++)
          other[j] ^= pivot_column[j];
    }
    pivot_column[pivot_word] |= pivot_bit;
  }

  for (c = 0; c < ncols; c++) {
    const uint64_t *column;
    uint64_t dependency_bit;
    if (pivot_row[c] != no_row)
      continue;
    if (dependencies++ < skip_dependencies)
      continue;
    if (packed >= 64UL)
      continue;
    dependency_bit = TS_BIT64(packed++);
    result[c] |= dependency_bit;
    column = matrix + (size_t)c * row_words;
    for (i = 0; i < c; i++) {
      unsigned long pivot = pivot_row[i];
      if (pivot != no_row &&
          (column[pivot >> 6] & TS_BIT64(pivot & 63UL)))
        result[i] |= dependency_bit;
    }
  }
  *total_dependencies = dependencies;
  if (packed == 0) {
    free(result);
    result = NULL;
    goto done;
  }

  /* Verify against the untouched sparse matrix, not the elimination copy. */
  check = (uint64_t *)ts_calloc_array(nrows, sizeof(*check));
  for (c = 0; c < ncols; c++) {
    uint64_t value = result[c];
    if (value == 0)
      continue;
    for (i = 0; i < columns[c].weight; i++)
      check[columns[c].data[i]] ^= value;
  }
  for (i = 0; i < nrows; i++)
    if (check[i] != 0)
      croak("tinysiqs128: dense nullspace verification failed");
  free(check);
  *mask = packed == 64UL ? (uint64_t)-1 : TS_BIT64(packed) - 1U;

done:
  free(used_rows);
  free(matrix);
  free(pivot_row);
  free(row_map);
  free(row_active);
  return result;
}

static void ts_reset_touched_factors(ts_ctx_t *ctx) {
  uint32_t i;
  for (i = 0; i < ctx->factor_touched_count; i++)
    ctx->factor_counts[ctx->factor_touched[i]] = 0;
  ctx->factor_touched_count = 0;
}

static void ts_touch_factor(ts_ctx_t *ctx, uint32_t row,
                            uint32_t exponent) {
  uint32_t old;
  if (row > ctx->params.fb_size)
    croak("tinysiqs128: dependency factor row is out of range");
  old = ctx->factor_counts[row];
  if (old == 0) {
    if (ctx->factor_touched_count > ctx->params.fb_size)
      croak("tinysiqs128: touched-factor capacity exceeded");
    ctx->factor_touched[ctx->factor_touched_count++] = (uint16_t)row;
  }
  if (exponent > UINT32_MAX - old)
    croak("tinysiqs128: dependency exponent overflow");
  ctx->factor_counts[row] = old + exponent;
}

static uint128_t ts_mont_neg(uint128_t value, const ts_mont_t *mont) {
  return value == 0 ? 0 : mont->n - value;
}

static int ts_test_dependencies(ts_ctx_t *ctx,
                                const ts_la_col_t *columns,
                                unsigned long ncols,
                                const uint64_t *nullrows,
                                uint64_t mask) {
  uint32_t dependency;
  ts_ensure_mont(ctx);
  for (dependency = 0; dependency < 64U && ctx->factor == 0;
       dependency++) {
    uint128_t lhs = ctx->mont.one, rhs = ctx->mont.one;
    uint128_t normal_lhs, normal_rhs, divisor;
    unsigned long i;
    uint32_t j;
    if (!(mask & TS_BIT64(dependency)))
      continue;
    ts_reset_touched_factors(ctx);
    for (i = 0; i < ncols; i++) {
      const ts_full_relation_t *relation;
      if (!(nullrows[i] & TS_BIT64(dependency)))
        continue;
      relation = &ctx->full[columns[i].orig];
      lhs = ts_mont_mul(lhs, ts_mont_enter(relation->y, &ctx->mont),
                        &ctx->mont);
      {
        const ts_factor_t *factors =
            ctx->full_factors + relation->factor_offset;
        for (j = 0; j < relation->nfactors; j++)
          ts_touch_factor(ctx, factors[j].row, factors[j].exponent);
      }
    }
    for (j = 0; j < ctx->factor_touched_count; j++) {
      uint32_t row = ctx->factor_touched[j];
      uint32_t exponent = ctx->factor_counts[row];
      if (exponent & 1U)
        croak("tinysiqs128: dependency is not a square");
      if (row == 0) {
        if ((exponent >> 1) & 1U)
          rhs = ts_mont_neg(rhs, &ctx->mont);
      } else if (exponent != 0) {
        uint128_t base = ts_mont_enter(ctx->fb[row - 1U].p,
                                      &ctx->mont);
        rhs = ts_mont_mul(rhs,
                          ts_mont_pow_u32(base, exponent >> 1, &ctx->mont),
                          &ctx->mont);
      }
    }
    normal_lhs = ts_mont_exit(lhs, &ctx->mont);
    normal_rhs = ts_mont_exit(rhs, &ctx->mont);
    divisor = ts_gcd128(normal_lhs >= normal_rhs
                      ? normal_lhs - normal_rhs
                      : normal_rhs - normal_lhs, ctx->n);
    if (divisor <= 1 || divisor >= ctx->n) {
      uint128_t sum = ts_addmod128(normal_lhs, normal_rhs, ctx->n);
      divisor = ts_gcd128(sum, ctx->n);
    }
    if (divisor > 1 && divisor < ctx->n && ctx->n % divisor == 0)
      ctx->factor = divisor;
    ts_reset_touched_factors(ctx);
  }
  return ctx->factor != 0;
}

static int ts_solve(ts_ctx_t *ctx) {
  unsigned long original_cols, nrows, ncols, i;
  unsigned long skip = 0, total = 0;
  ts_la_col_t *columns;
  if (ctx->full_count == 0)
    return 0;
  columns = ts_build_matrix(ctx, &nrows, &ncols);
  original_cols = ncols;
  ts_reduce_matrix(nrows, &ncols, columns);
#if TS_DIAGNOSTIC
  fprintf(stderr, "# tinysiqs128 matrix %lu/%lu columns, %lu rows\n",
          ncols, original_cols, nrows);
#endif
  while (ncols != 0 && ctx->factor == 0) {
    uint64_t mask = 0;
    uint64_t *nullrows = ts_dense_nullspace_batch(
        nrows, ncols, columns, skip, &mask, &total);
#if TS_DIAGNOSTIC
    fprintf(stderr, "# tinysiqs128 dependency batch skip=%lu total=%lu\n",
            skip, total);
#endif
    if (nullrows == NULL)
      break;
    (void)ts_test_dependencies(ctx, columns, ncols, nullrows, mask);
    free(nullrows);
    if (total <= skip + 64UL)
      break;
    skip += 64UL;
  }
  for (i = 0; i < original_cols; i++)
    free(columns[i].data);
  free(columns);
  return ctx->factor != 0;
}

#undef TS_BIT64

/*----------------------------------------------------------------------------
 * Collection driver and public entry point
 *----------------------------------------------------------------------------*/

static int ts_ctx_init(ts_ctx_t *ctx, uint128_t n,
                       const ts_policy_band_t *policy) {
  uint64_t seed;
  memset(ctx, 0, sizeof(*ctx));
  ctx->n = n;
  if (!ts_resolve_parameters(&ctx->params, n, policy))
    return 0;
  ctx->multiplier = ts_choose_multiplier(n, ctx->params.fb_size);
  if (ctx->multiplier == 0 || n > TS_U128_MAX / ctx->multiplier)
    return 0;
  ctx->kn = n * ctx->multiplier;
  ctx->params.d = (ctx->kn & 7U) == 1U ? 2U : 1U;
  seed = (uint64_t)(n % UINT32_C(4294967291));
  seed ^= (uint64_t)(n % UINT32_C(4294967279)) << 32;
  seed ^= (uint64_t)ctx->params.bits << 17;
  ctx->poly_rng.state = ts_mix64(seed);
  return 1;
}

static int ts_ctx_allocate(ts_ctx_t *ctx) {
  uint32_t rows, i;
  size_t sieve_alloc, pmax_alloc;
  uint64_t lp_bound, pmax2;
  if (!ts_build_factor_base(ctx))
    return 0;
  ctx->prime = (uint32_t *)ts_malloc_array(ctx->params.fb_size,
                                           sizeof(*ctx->prime));
  ctx->root1 = (uint32_t *)ts_malloc_array(ctx->params.fb_size,
                                           sizeof(*ctx->root1));
  ctx->root2 = (uint32_t *)ts_malloc_array(ctx->params.fb_size,
                                           sizeof(*ctx->root2));
  ctx->reciprocal = (uint32_t *)ts_malloc_array(
      ctx->params.fb_size, sizeof(*ctx->reciprocal));
  for (i = 0; i < ctx->params.fb_size; i++) {
    ctx->prime[i] = ctx->fb[i].p;
    ctx->root1[i] = 0;
    ctx->root2[i] = TS_NO_ROOT;
    ctx->reciprocal[i] =
        (uint32_t)(UINT64_C(0x100000000) / ctx->prime[i]);
  }

  lp_bound = (uint64_t)ctx->params.lp_multiplier *
             ctx->largest_fb_prime;
  pmax2 = (uint64_t)ctx->largest_fb_prime * ctx->largest_fb_prime;
  if (lp_bound > UINT32_MAX || lp_bound >= pmax2)
    croak("tinysiqs128: invalid one-large-prime bound");
  ts_set_log_weights(ctx);

  if (ctx->params.M > UINT32_MAX / 2U)
    croak("tinysiqs128: sieve interval is too large");
  ctx->sieve_length = 2U * ctx->params.M;
  sieve_alloc = 2U * (size_t)ctx->sieve_length;
  pmax_alloc = (size_t)ctx->largest_fb_prime + 1U;
  if (sieve_alloc < pmax_alloc)
    sieve_alloc = pmax_alloc;
  ctx->sieve = (uint8_t *)ts_calloc_array(sieve_alloc + 8U, 1U);
  ctx->candidate_at = (uint16_t *)ts_calloc_array(
      ctx->sieve_length, sizeof(*ctx->candidate_at));

  ctx->resieve_one_subtract_index = 0;
  while (ctx->resieve_one_subtract_index < ctx->params.fb_size &&
         ctx->prime[ctx->resieve_one_subtract_index] <= ctx->params.M)
    ctx->resieve_one_subtract_index++;

  rows = ctx->params.fb_size + 1U;
  ctx->factor_counts = (uint32_t *)ts_calloc_array(
      rows, sizeof(*ctx->factor_counts));
  ctx->factor_touched = (uint16_t *)ts_malloc_array(
      rows, sizeof(*ctx->factor_touched));
  if (ctx->params.lp_multiplier > 1U)
    ts_lp_init(&ctx->one_lp, ctx->params.fb_size * 2U);
  ts_hashset_init(&ctx->a_hashes, 1024U);
  return 1;
}

static void ts_ctx_clear(ts_ctx_t *ctx) {
  free(ctx->full);
  free(ctx->full_factors);
  free(ctx->fb);
  free(ctx->prime);
  free(ctx->root1);
  free(ctx->root2);
  free(ctx->reciprocal);
  free(ctx->sieve);
  free(ctx->candidate_at);
  free(ctx->candidate_at_wide);
  free(ctx->candidates);
  free(ctx->hits);
  free(ctx->factor_counts);
  free(ctx->factor_touched);
  ts_lp_clear(&ctx->one_lp);
  ts_hashset_clear(&ctx->a_hashes);
  memset(ctx, 0, sizeof(*ctx));
}

static int ts_collect_relations(ts_ctx_t *ctx, ts_poly_t *poly,
                                uint32_t *families,
                                uint32_t *polynomials) {
  while (ctx->factor == 0 &&
         ctx->matrix_next_target <= ctx->matrix_target_limit &&
         ctx->full_count < ctx->matrix_next_target &&
         ctx->full_count < TS_MAX_RELATIONS &&
         *families < ctx->params.max_families) {
    if (!ts_new_family(ctx, poly))
      return 0;
    (*families)++;
    for (;;) {
      ts_sieve_polynomial(ctx, poly);
      (*polynomials)++;
      if (ctx->factor != 0 ||
          ctx->matrix_next_target > ctx->matrix_target_limit ||
          ctx->full_count >= TS_MAX_RELATIONS)
        break;
      if (ctx->full_count >= ctx->matrix_next_target ||
          ctx->full_count >= TS_MAX_RELATIONS)
        break;
      if (!ts_next_B(ctx, poly))
        break;
    }
  }
  return ctx->factor != 0;
}

static int ts_run(ts_ctx_t *ctx) {
  ts_poly_t poly;
  uint32_t retry_batch = ctx->params.fb_size / 16U;
  uint32_t families = 0, polynomials = 0;
  uint32_t target_limit = ctx->params.fb_size + 1U + TS_MAX_EXTRA_RELS;
  if (retry_batch < 8U)
    retry_batch = 8U;
  if (retry_batch > TS_RETRY_BATCH_MAX)
    retry_batch = TS_RETRY_BATCH_MAX;
  ctx->matrix_next_target = ctx->params.target_relations;
  ctx->matrix_last_count = 0;
  ctx->matrix_retry_batch = retry_batch;
  ctx->matrix_target_limit = target_limit;
  if (!ts_ctx_allocate(ctx))
    return ctx->factor != 0;
#if TS_DIAGNOSTIC
  fprintf(stderr, "# tinysiqs128 policy=%s bits=%u k=%u FB=%u pmax=%u "
          "M=%u q=%u d=%u K=%u target=%u\n",
          ctx->params.name, ctx->params.bits, ctx->multiplier,
          ctx->params.fb_size, ctx->largest_fb_prime, ctx->params.M,
          ctx->params.q_count, ctx->params.d, ctx->params.lp_multiplier,
          ctx->params.target_relations);
#endif
  if (!ts_poly_init(ctx, &poly))
    return 0;

  (void)ts_collect_relations(ctx, &poly, &families, &polynomials);
  if (ctx->factor == 0 && ctx->full_count > ctx->matrix_last_count) {
    ctx->matrix_last_count = ctx->full_count;
    (void)ts_solve(ctx);
  }
#if TS_DIAGNOSTIC
  fprintf(stderr, "# tinysiqs128 collect next=%u full=%u families=%u "
          "polynomials=%u\n", ctx->matrix_next_target, ctx->full_count,
          families, polynomials);
#endif
  ts_poly_clear(ctx, &poly);
  (void)families;
  (void)polynomials;
  return ctx->factor != 0;
}

static uint128_t ts_try_policy(uint128_t n,
                               const ts_policy_band_t *policy) {
  ts_ctx_t ctx;
  uint128_t factor = 0;
  uint64_t root;
  if (!ts_ctx_init(&ctx, n, policy))
    return 0;
  root = ts_isqrt128(ctx.kn);
  if ((uint128_t)root * root == ctx.kn) {
    factor = ts_gcd128(root, n);
    if (factor <= 1 || factor >= n)
      factor = 0;
  }
  if (factor == 0) {
    (void)ts_run(&ctx);
    factor = ctx.factor;
  }
  ts_ctx_clear(&ctx);
  return factor > 1 && factor < n && n % factor == 0 ? factor : 0;
}

uint128_t tinysiqs128(uint128_t n) {
  static const ts_policy_band_t *const q1_recoveries[] = {
    &ts_recovery_q1_8k,
    &ts_recovery_q1_16k,
    &ts_recovery_q1_legacy
  };
  static const ts_policy_band_t *const q2_recoveries[] = {
    &ts_recovery_42_96k,
    &ts_recovery_42_192k,
    &ts_recovery_42_384k,
    &ts_recovery_42_1m
  };
  const ts_policy_band_t *primary;
  uint128_t factor, root;
  uint32_t bits, i, exponent;
  if (n < 4)
    return 0;
  if (!(n & 1U))
    return 2;
  if (n % 3U == 0)
    return 3;
  if (n % 5U == 0)
    return n == 5 ? 0 : 5;
  root = ts_isqrt128(n);
  if (root > 1 && root < n && root * root == n)
    return root;
  bits = ts_bits128(n);
  if (bits <= 31U)
    return ts_trial_factor32((uint32_t)n);

  root = ts_exact_root(n, 3U);
  if (root > 1 && root < n)
    return root;
  if (is_prime128(n))
    return 0;

  /* Compact SQUFOF avoids the SIQS setup cliff through 51 bits.  Failure is
   * harmless; the normal SIQS policy and its recoveries remain authoritative. */
  if (bits <= 51U) {
    factor = ts_squfof64((uint64_t)n, 100000U);
    if (factor != 0)
      return factor;
    if (bits == 32U)
      return ts_trial_factor32((uint32_t)n);
  }

  primary = ts_primary_policy(bits);
  factor = ts_try_policy(n, primary);
  if (factor != 0)
    return factor;
  if (bits <= 36U) {
    for (i = 0; i < (uint32_t)(sizeof(q1_recoveries) /
                                sizeof(q1_recoveries[0])); i++) {
      factor = ts_try_policy(n, q1_recoveries[i]);
      if (factor != 0)
        return factor;
    }
  } else if (bits >= 42U && bits <= 49U) {
    for (i = 0; i < (uint32_t)(sizeof(q2_recoveries) /
                                sizeof(q2_recoveries[0])); i++) {
      factor = ts_try_policy(n, q2_recoveries[i]);
      if (factor != 0)
        return factor;
    }
  }

  /* Prime powers can make every QS dependency trivial.  Cubes were screened
   * before collection; test the remaining odd prime exponents on exhaustion. */
  for (exponent = 5U; exponent <= bits; exponent += 2U) {
    if (!ts_is_prime_u32(exponent))
      continue;
    root = ts_exact_root(n, exponent);
    if (root > 1 && root < n)
      return root;
  }
  return 0;
}

#endif /* HAVE_FACTOR128 */
