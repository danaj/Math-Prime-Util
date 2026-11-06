/* Fixed-width modular arithmetic for an odd modulus n > 1.
 *
 * Adapted from the tecm64 helpers in Math::Prime::Util's factor.c.
 * Uses the shared types/inlining in ptypes.h.
 * No GMP dependency.
 * With unsigned 128-bit arithmetic, values are in Montgomery form (R = 2^64),
 * otherwise the same interface uses ordinary residues and portable
 * double-and-add multiplication.  Set MONT64_HAVE_UINT128=0 to exercise that
 * fallback on any compiler.
 *
 * Add/subtract and Montgomery multiply operands must be reduced modulo n.
 * Enter/exit convert representations; ctx->one is the represented value 1.
 * Inline definitions keep the hot arithmetic in the caller's loops.
 *
 * Copyright (c) 2026 Dana Jacobsen.
 */
#ifndef MPU_MONT64_H
#define MPU_MONT64_H

#include "ptypes.h"

#if HAVE_UINT64

#ifndef MONT64_HAVE_UINT128
# define MONT64_HAVE_UINT128 HAVE_UINT128
#endif

typedef struct {
  uint64_t n, one;
#if MONT64_HAVE_UINT128
  uint64_t ninv, r2;
#endif
} mont64_t;

static INLINE MAYBE_UNUSED uint64_t mont64_add(uint64_t a, uint64_t b, uint64_t n) {
  uint64_t gap = n - b;
  return a >= gap ? a - gap : a + b;
}
static INLINE MAYBE_UNUSED uint64_t mont64_sub(uint64_t a, uint64_t b, uint64_t n) {
  return a >= b ? a - b : n - (b - a);
}
/* Ordinary modular multiplication; inputs need not already be reduced. */
static INLINE MAYBE_UNUSED uint64_t mont64_mulmod(uint64_t a, uint64_t b, uint64_t n) {
#if MONT64_HAVE_UINT128
  return (uint64_t)(((uint128_t)a * b) % n);
#else
  uint64_t r = 0;
  if (a >= n) a %= n;
  if (b >= n) b %= n;
  if (a < b) { uint64_t t = a; a = b; b = t; }
  while (b != 0) {
    if (b & 1U) r = mont64_add(r, a, n);
    b >>= 1;
    if (b != 0) a = mont64_add(a, a, n);
  }
  return r;
#endif
}
static INLINE MAYBE_UNUSED void mont64_init(mont64_t *ctx, uint64_t n) {
  ctx->n = n;
#if MONT64_HAVE_UINT128
  {
    uint64_t x = (3 * n) ^ 2U;
    x *= (uint64_t)2 - n * x;
    x *= (uint64_t)2 - n * x;
    x *= (uint64_t)2 - n * x;
    x *= (uint64_t)2 - n * x;
    ctx->ninv = (uint64_t)0 - x;
    /* R mod n uses just one 64-bit remainder. For odd n > 1 it is nonzero. */
    ctx->one = ((uint64_t)-1) % n + 1;
    ctx->r2 = mont64_mulmod(ctx->one, ctx->one, n);
  }
#else
  ctx->one = 1;
#endif
}
static INLINE MAYBE_UNUSED uint64_t mont64_mul(uint64_t a, uint64_t b, const mont64_t *ctx) {
#if MONT64_HAVE_UINT128
  uint128_t ab = (uint128_t)a * b;
  uint64_t lo = (uint64_t)ab, hi = (uint64_t)(ab >> 64);
  uint64_t m = lo * ctx->ninv;
  uint64_t mn_hi = (uint64_t)(((uint128_t)m * ctx->n) >> 64);
  uint64_t u = hi + mn_hi;
  int carry = u < hi;
  uint64_t v = u + (lo != 0);
  /* The low-word sum is zero; its carry is exactly (lo != 0).
   * Preserve the carry above bit 127 when the modulus exceeds 2^63. */
  if (carry || v < u || v >= ctx->n) v -= ctx->n;
  return v;
#else
  return mont64_mulmod(a, b, ctx->n);
#endif
}
static INLINE MAYBE_UNUSED uint64_t mont64_enter(uint64_t a, const mont64_t *ctx) {
  a %= ctx->n;
#if MONT64_HAVE_UINT128
  return mont64_mul(a, ctx->r2, ctx);
#else
  return a;
#endif
}
static INLINE MAYBE_UNUSED uint64_t mont64_exit(uint64_t a, const mont64_t *ctx) {
#if MONT64_HAVE_UINT128
  return mont64_mul(a, 1, ctx);
#else
  (void)ctx;
  return a;
#endif
}

#endif /* HAVE_UINT64 */

#endif /* MPU_MONT64_H */
