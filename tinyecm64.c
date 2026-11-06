/* tinyecm64 1.0.1
 *
 * ECM for uint64_t inputs.  Portable C, no GMP, no heap allocation.
 *
 * In 2026 I wrote a 64-bit ECM implementation (tecm64) from scratch for
 * Math::Prime::Util.  This substantial update retains the same portable
 * uint64_t concept but includes several enhancements, making it much faster.
 *
 * Representative CPU microseconds/input for balanced semiprimes:
 *
 *   Input bits       32    40    48    56    64
 *   -------------- ----- ----- ----- ----- -----
 *   Apple M1 Pro     1.7   3.1   6.5  13.8  29.8
 *   AMD EPYC 9B45    1.8   4.6   9.6  20.0  43.1
 *
 * One thread, try_small=0.  Clang 21 (ARM), GCC 13 (AMD),
 * -O3 -march=native -falign-loops=32.  Each machine had a total of 11
 * misses in 1.5 million attempts.  Misses included; no fallback timed.
 *
 * Scalar and four-curve code share one caller-owned RNG interface.
 * Montgomery's PRAC chains are generated offline and packed into 32-bit words.
 * No assembly, explicit SIMD intrinsics, GMP, or heap allocation.
 *
 * Copyright (c) 2026 Dana Jacobsen.  See LICENSE for redistribution terms.
 */

#include "tinyecm64.h"

#if HAVE_UINT64

#include "mont64.h"
#include <string.h>

#ifndef UINT64_C
#define UINT64_C(value) ((uint64_t)(value))
#endif

/* Stage-2 products share a GCD across eight windows. Failed products are
 * recovered using every unchecked window's own giant point. */
enum { ECM64_GCD_WINDOWS = 8, ECM64_LANES = 4 };

#if defined(__GNUC__) || defined(__clang__)
#define ECM64_INLINE static __inline__ __attribute__((always_inline))
#else
#define ECM64_INLINE static INLINE
#endif

/* Subtractive REDC avoids needing a carry beyond the double-width product.
 * Domain entry/exit retain the shared mont64.h representation. */
static INLINE uint64_t multiply(uint64_t a, uint64_t b, const mont64_t *ctx) {
#if MONT64_HAVE_UINT128
  uint128_t product = (uint128_t)a * b;
  uint64_t high = (uint64_t)(product >> 64);
  uint64_t m = (uint64_t)product * ((uint64_t)0 - ctx->ninv);
  uint64_t q = (uint64_t)(((uint128_t)m * ctx->n) >> 64);
  uint64_t r = high - q;
  return high < q ? r + ctx->n : r;
#else
  return mont64_mul(a, b, ctx);
#endif
}
static INLINE uint64_t leading_bit(uint64_t k) {
#if defined(__GNUC__) || defined(__clang__)
  return UINT64_C(1) << (63 - (unsigned)__builtin_clzll(k));
#else
  uint64_t bit = UINT64_C(1) << 63;
  while (!(k & bit))
    bit >>= 1;
  return bit;
#endif
}

/* Point/inverse arithmetic from MPU's factor. */
static INLINE uint64_t addmod(uint64_t a, uint64_t b, uint64_t n) {
  uint64_t gap = n - b;
  return a >= gap ? a - gap : a + b;
}

static INLINE uint64_t submod(uint64_t a, uint64_t b, uint64_t n) {
  return (a >= b) ? a - b : n - (b - a);
}

static INLINE uint64_t halfmod(uint64_t a, uint64_t n) {
  return (a & 1) ? (a >> 1) + (n >> 1) + 1 : a >> 1;
}

static unsigned trailing_zeros(uint64_t n) {
#if defined(__GNUC__) || defined(__clang__)
  return (unsigned)__builtin_ctzll(n);
#else
  unsigned k = 0;
  while (!(n & 1)) {
    n >>= 1;
    ++k;
  }
  return k;
#endif
}
static uint64_t gcd64(uint64_t a, uint64_t b) {

  if (!a) return b;
  if (!b) return a;
  unsigned shift = trailing_zeros(a | b);
  a >>= trailing_zeros(a);
  b >>= trailing_zeros(b);
  while (a != b) {
    uint64_t ab = a - b, ba = b - a;
    uint64_t small = a >= b ? b : a;
    uint64_t delta = a >= b ? ab : ba;
    /* ab and ba have the same trailing zeros even when one wraps. */
    unsigned zeros = trailing_zeros(ab);
    a = small;
    b = delta >> zeros;
  }
  return a << shift;
}

/* Binary extended modular inverse for odd n.  Returns 0 if gcd(a,n) > 1. */
static uint64_t inverse64(uint64_t a, uint64_t n) {

  uint64_t u, v, x1, x2;
  u = a % n;
  if (u == 0) return 0;
  v = n;
  x1 = 1;
  x2 = 0;
  while (u != 1 && v != 1) {
    while ((u & 1) == 0) {
      u >>= 1;
      x1 = (x1 & 1) ? halfmod(x1, n) : x1 >> 1;
    }
    while ((v & 1) == 0) {
      v >>= 1;
      x2 = (x2 & 1) ? halfmod(x2, n) : x2 >> 1;
    }
    if (u >= v) {
      u -= v;
      x1 = submod(x1, x2, n);
    } else {
      v -= u;
      x2 = submod(x2, x1, n);
    }
    if (u == 0 || v == 0) return 0;
  }
  return (u == 1) ? x1 : x2;
}

/* Use the project's shared full-width Montgomery implementation. */
typedef struct {
  uint64_t X, Z;
} ecm_point;

static INLINE void point_double(ecm_point *R, const ecm_point *P, uint64_t A24,
                                const mont64_t *ctx) {

  uint64_t n = ctx->n;
  uint64_t u = multiply(submod(P->X, P->Z, n), submod(P->X, P->Z, n), ctx);
  uint64_t v = multiply(addmod(P->X, P->Z, n), addmod(P->X, P->Z, n), ctx);
  uint64_t w = submod(v, u, n);
  R->X = multiply(u, v, ctx);
  R->Z = multiply(w, addmod(u, multiply(A24, w, ctx), n), ctx);
}

static INLINE void point_add(ecm_point *R, const ecm_point *P, const ecm_point *Q,
                             const ecm_point *Pm, const mont64_t *ctx) {

  uint64_t n = ctx->n;
  uint64_t u = multiply(submod(P->X, P->Z, n), addmod(Q->X, Q->Z, n), ctx);
  uint64_t v = multiply(addmod(P->X, P->Z, n), submod(Q->X, Q->Z, n), ctx);
  uint64_t s = addmod(u, v, n);
  uint64_t d = submod(u, v, n);
  /* PRAC also aliases output with the difference point; preserve both inputs. */
  uint64_t x = multiply(multiply(s, s, ctx), Pm->Z, ctx);
  uint64_t z = multiply(multiply(d, d, ctx), Pm->X, ctx);
  R->X = x;
  R->Z = z;
}

static void point_multiply(ecm_point *R, const ecm_point *P, uint64_t k, uint64_t A24,
                           const mont64_t *ctx) {
  ecm_point R0, R1;
  uint64_t bit;
  if (k == 1) {
    *R = *P;
    return;
  }
  R0 = *P;
  point_double(&R1, &R0, A24, ctx);
  if (k == 2) {
    *R = R1;
    return;
  }
  bit = leading_bit(k);
  for (bit >>= 1; bit; bit >>= 1) {
    if (k & bit) {
      point_add(&R0, &R0, &R1, P, ctx);
      point_double(&R1, &R1, A24, ctx);
    } else {
      point_add(&R1, &R0, &R1, P, ctx);
      point_double(&R0, &R0, A24, ctx);
    }
  }
  *R = R0;
}

/* R is coprime to odd N, so a Montgomery residue has the same GCD. */
static INLINE uint64_t gcd_residue(uint64_t value, const mont64_t *ctx) {
  return gcd64(value, ctx->n);
}

/* Scalar binary stage 1; PRAC is retained for the four-curve tail. */

static void scalar_stage1(ecm_point *out, const ecm_point *in, uint32_t k, uint64_t a24,
                          const mont64_t *ctx) {
  ecm_point a = *in, b;
  /* Binary scalar ladder, omitting the unused final register update. */
  uint64_t bit;
  if (!(k & (k - 1))) {
    for (; k > 1; k >>= 1)
      point_double(&a, &a, a24, ctx);
    *out = a;
    return;
  }
  point_double(&b, &a, a24, ctx);
  bit = leading_bit(k);
  for (bit >>= 1; bit > 1; bit >>= 1) {
    if (k & bit) {
      point_add(&a, &a, &b, in, ctx);
      point_double(&b, &b, a24, ctx);
    } else {
      point_add(&b, &a, &b, in, ctx);
      point_double(&a, &a, a24, ctx);
    }
  }
  /* Only the returned register is needed for the final bit. */
  if (k & 1)
    point_add(out, &a, &b, in, ctx);
  else
    point_double(out, &a, a24, ctx);
}

/* Odd multipliers use three-bit PRAC opcodes, indexed by k >> 1.
 * A zero high part terminates the word; powers of two use doubling. */
static const uint32_t chain_words[] = {
    0,        0,      0,     18,       0,       19,       434,      0, 186,     3474,
    0,        155,    27794, 1490,     3475,    28050,    0,        0, 11411,   0,
    27795,    11954,  0,     28051,    222642,  0,        222355,   0, 0,       11923,
    27803,    0,      0,     224403,   0,       222674,   97714,    0, 0,       1781170,
    1797522,  11155,  0,     0,        1797554, 0,        0,        0, 1666195, 0,
    97683,    765106, 0,     781714,   1781139, 0,        765362,   0, 0,       0,
    14378386, 0,      81555, 73107,    0,       14380434, 0,        0, 652690,  1795227,
    0,        0,      0,     0,        5991571, 1779666,  0,        0, 781498,  0,
    0,        781715, 0,     13315475, 6253746, 0,        6251954,  0, 0,       106524050,
    13331859, 0,      0,     0,        0,       6254002,  14253490, 0, 6253716, 14380435,
};

static void scalar_stage1_prac(ecm_point *out, const ecm_point *in, uint32_t k, uint64_t a24,
                               const mont64_t *ctx) {
  if (!(k & 1) || k == 1) {
    scalar_stage1(out, in, k, a24, ctx);
    return;
  }
  ecm_point a, b = *in, c = *in, t, u;
  /* k comes from the immutable plans; every such odd multiplier has a
   * verified nonempty word. Its highest unused bits terminate decoding. */
  uint32_t word = chain_words[k >> 1];
  point_double(&a, &b, a24, ctx);
  for (; word; word >>= 3) {
    switch (word & 7) {
      /* Three-bit opcode, with any swap folded into the operation. */
    case 1:
      point_add(&t, &a, &b, &c, ctx);
      point_add(&u, &t, &a, &b, ctx);
      point_add(&b, &b, &t, &a, ctx);
      a = u;
      break;
    case 2:
      point_add(&c, &b, &a, &c, ctx);
      t = b;
      b = c;
      c = t;
      break;
    case 3:
      point_add(&b, &b, &a, &c, ctx);
      point_double(&a, &a, a24, ctx);
      break;
    case 4:
      point_add(&c, &c, &a, &b, ctx);
      point_double(&a, &a, a24, ctx);
      break;
    case 5:
      t = a;
      a = b;
      b = t;
      point_add(&t, &a, &b, &c, ctx);
      point_add(&u, &t, &a, &b, ctx);
      point_add(&b, &b, &t, &a, ctx);
      a = u;
      break;
    case 6:
      t = a;
      a = b;
      b = t;
      point_add(&c, &b, &a, &c, ctx);
      t = b;
      b = c;
      c = t;
      break;
    case 7:
      t = a;
      a = b;
      b = t;
      point_add(&b, &b, &a, &c, ctx);
      point_double(&a, &a, a24, ctx);
      break;
    }
  }
  point_add(out, &a, &b, &c, ctx);
}

typedef struct {
  uint32_t b1, b2, npowers, d, windows, babies;
  const uint32_t *powers;
  const uint8_t *slots, *indices;
  const uint16_t *offsets;
} ecm64_plan;

/* Offline wheel schedules cover every stage-2 prime in (B1,B2]. */
static const uint32_t powers_27_675[] = {
    16, 27, 25, 7, 11, 13, 17, 19, 23,
};
static const uint8_t slots_27_675[] = {
    255, 0, 255, 255, 255, 255, 255, 1, 255, 255, 255, 2, 255, 3, 255, 255,
};
static const uint8_t indices_27_675[] = {
    0, 1, 2, 3, 3, 1, 0, 2, 2, 1, 0, 3, 3, 2, 1, 3, 2, 0, 1, 3, 1, 0, 2, 3, 2, 0,
    3, 2, 1, 0, 3, 1, 0, 2, 1, 2, 3, 3, 0, 1, 3, 2, 1, 0, 2, 1, 0, 2, 0, 3, 2, 1,
    0, 3, 3, 0, 1, 2, 2, 1, 0, 3, 0, 1, 3, 1, 0, 3, 1, 0, 3, 2, 0, 3, 1, 0,
};
static const uint16_t offsets_27_675[] = {
    0, 4, 8, 12, 15, 19, 23, 26, 30, 34, 37, 40, 44, 47, 50, 54, 58, 62, 64, 67, 70, 73, 76,
};
static const uint32_t powers_70_1750[] = {
    64, 27, 25, 49, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
};
static const uint8_t slots_70_1750[] = {
    255, 0,   255, 255, 255, 255, 255, 1,   255, 255, 255, 2,   255, 3,   255, 255,
    255, 4,   255, 5,   255, 255, 255, 6,   255, 255, 255, 255, 255, 7,   255, 8,
    255, 255, 255, 255, 255, 9,   255, 255, 255, 10,  255, 11,  255, 255,
};
static const uint8_t indices_70_1750[] = {
    5,  4,  2,  1, 0,  3,  6,  9, 10, 11, 10, 8,  7, 6,  4,  3,  1,  0,  2,  5, 11, 10, 9,
    8,  7,  5,  3, 1,  0,  2,  6, 11, 7,  6,  3,  2, 1,  0,  5,  9,  10, 10, 8, 7,  5,  4,
    2,  1,  0,  3, 9,  10, 9,  8, 5,  4,  0,  1,  6, 7,  11, 9,  8,  7,  6,  4, 3,  2,  0,
    11, 9,  7,  5, 2,  0,  1,  3, 6,  8,  10, 10, 9, 6,  3,  0,  2,  4,  5,  7, 11, 11, 10,
    9,  6,  5,  4, 3,  1,  2,  7, 11, 9,  6,  5,  3, 1,  0,  7,  8,  10, 10, 8, 7,  5,  4,
    2,  1,  3,  6, 9,  11, 10, 5, 4,  1,  0,  2,  6, 8,  11, 11, 9,  8,  7,  6, 2,  0,  4,
    5,  10, 11, 8, 7,  6,  2,  4, 10, 8,  4,  3,  2, 1,  0,  5,  11, 11, 10, 9, 8,  5,  1,
    0,  3,  6,  7, 10, 9,  6,  5, 3,  2,  1,  0,  4, 11, 11, 10, 4,  3,  2,  0, 6,  8,  9,
};
static const uint16_t offsets_70_1750[] = {
    0, 9, 20, 31, 41, 51, 60, 69, 80, 90, 100, 110, 121, 130, 140, 146, 155, 165, 175, 184,
};
static const uint32_t powers_125_3125[] = {
    64, 81, 125, 49, 121, 13, 17, 19, 23, 29, 31,  37,  41,  43,  47,
    53, 59, 61,  67, 71,  73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
};
static const uint8_t slots_125_3125[] = {
    255, 0,   255, 255, 255, 255, 255, 1,   255, 255, 255, 2,   255, 3,   255, 255, 255, 4,   255,
    5,   255, 255, 255, 6,   255, 255, 255, 255, 255, 7,   255, 8,   255, 255, 255, 255, 255, 9,
    255, 255, 255, 10,  255, 11,  255, 255, 255, 12,  255, 13,  255, 255, 255, 14,  255, 255, 255,
    255, 255, 15,  255, 16,  255, 255, 255, 255, 255, 17,  255, 255, 255, 18,  255, 19,  255, 255,
};
static const uint8_t indices_125_3125[] = {
    6,  5,  3,  2,  0,  1,  4,  7,  8,  10, 11, 12, 13, 16, 19, 19, 18, 17, 16, 15, 13, 11, 9,
    8,  7,  6,  5,  4,  1,  2,  3,  12, 14, 18, 17, 16, 14, 13, 10, 8,  7,  5,  4,  2,  1,  0,
    3,  9,  15, 19, 15, 14, 11, 9,  8,  7,  6,  3,  1,  0,  4,  5,  10, 12, 16, 19, 19, 17, 15,
    13, 10, 8,  6,  4,  2,  1,  0,  5,  9,  12, 16, 18, 19, 18, 16, 12, 11, 10, 9,  6,  5,  4,
    3,  1,  2,  7,  14, 17, 19, 17, 15, 14, 10, 9,  8,  7,  5,  4,  2,  0,  3,  11, 12, 18, 13,
    12, 9,  7,  5,  3,  1,  0,  4,  6,  8,  15, 19, 18, 17, 16, 15, 14, 13, 12, 11, 8,  7,  6,
    2,  4,  19, 18, 17, 16, 14, 13, 12, 10, 7,  5,  4,  3,  2,  1,  0,  6,  8,  11, 15, 18, 17,
    14, 13, 11, 10, 9,  8,  7,  6,  3,  1,  4,  5,  12, 15, 19, 17, 15, 14, 12, 10, 6,  4,  3,
    2,  0,  8,  16, 18, 19, 19, 18, 16, 13, 11, 9,  5,  4,  0,  6,  7,  12, 14, 17, 19, 18, 16,
    12, 9,  8,  5,  4,  3,  2,  0,  7,  10, 11, 14, 18, 12, 11, 9,  7,  3,  2,  1,  0,  4,  5,
    6,  8,  15, 16, 17, 16, 15, 14, 13, 11, 7,  6,  5,  4,  2,  1,  0,  9,  10, 12, 19, 19, 12,
    7,  5,  2,  1,  0,  10, 11, 15, 17, 18, 17, 14, 11, 10, 9,  7,  6,  4,  3,  2,  1,  0,  5,
    8,  13, 19, 16, 15, 14, 13, 12, 8,  4,  3,  1,  0,  2,  7,  9,  17, 19, 16, 12, 11, 9,  8,
    7,  0,  2,  5,  6,  10, 13, 17, 18, 17, 16, 10, 8,  7,
};
static const uint16_t offsets_125_3125[] = {
    0,   15,  33,  50,  66,  82,  98,  113, 126, 140, 159,
    176, 190, 204, 219, 234, 251, 263, 278, 293, 307, 313,
};
static const uint32_t powers_165_4125[] = {
    128, 81, 125, 49, 121, 13, 17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,  61,  67,
    71,  73, 79,  83, 89,  97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163,
};
static const uint8_t slots_165_4125[] = {
    255, 0,   255, 255, 255, 255, 255, 1,   255, 255, 255, 2,   255, 3,   255, 255, 255, 4,   255,
    5,   255, 255, 255, 6,   255, 255, 255, 255, 255, 7,   255, 8,   255, 255, 255, 255, 255, 9,
    255, 255, 255, 10,  255, 11,  255, 255, 255, 12,  255, 13,  255, 255, 255, 14,  255, 255, 255,
    255, 255, 15,  255, 16,  255, 255, 255, 255, 255, 17,  255, 255, 255, 18,  255, 19,  255, 255,
};
static const uint8_t indices_165_4125[] = {
    4,  6,  7,  8,  10, 11, 12, 13, 16, 19, 19, 18, 17, 16, 15, 13, 11, 9,  8,  7,  6,  5,  4,  1,
    2,  3,  12, 14, 18, 17, 16, 14, 13, 10, 8,  7,  5,  4,  2,  1,  0,  3,  9,  15, 19, 15, 14, 11,
    9,  8,  7,  6,  3,  1,  0,  4,  5,  10, 12, 16, 19, 19, 17, 15, 13, 10, 8,  6,  4,  2,  1,  0,
    5,  9,  12, 16, 18, 19, 18, 16, 12, 11, 10, 9,  6,  5,  4,  3,  1,  2,  7,  14, 17, 19, 17, 15,
    14, 10, 9,  8,  7,  5,  4,  2,  0,  3,  11, 12, 18, 13, 12, 9,  7,  5,  3,  1,  0,  4,  6,  8,
    15, 19, 18, 17, 16, 15, 14, 13, 12, 11, 8,  7,  6,  2,  4,  19, 18, 17, 16, 14, 13, 12, 10, 7,
    5,  4,  3,  2,  1,  0,  6,  8,  11, 15, 18, 17, 14, 13, 11, 10, 9,  8,  7,  6,  3,  1,  4,  5,
    12, 15, 19, 17, 15, 14, 12, 10, 6,  4,  3,  2,  0,  8,  16, 18, 19, 19, 18, 16, 13, 11, 9,  5,
    4,  0,  6,  7,  12, 14, 17, 19, 18, 16, 12, 9,  8,  5,  4,  3,  2,  0,  7,  10, 11, 14, 18, 12,
    11, 9,  7,  3,  2,  1,  0,  4,  5,  6,  8,  15, 16, 17, 16, 15, 14, 13, 11, 7,  6,  5,  4,  2,
    1,  0,  9,  10, 12, 19, 19, 12, 7,  5,  2,  1,  0,  10, 11, 15, 17, 18, 17, 14, 11, 10, 9,  7,
    6,  4,  3,  2,  1,  0,  5,  8,  13, 19, 16, 15, 14, 13, 12, 8,  4,  3,  1,  0,  2,  7,  9,  17,
    19, 16, 12, 11, 9,  8,  7,  0,  2,  5,  6,  10, 13, 17, 18, 17, 16, 10, 8,  7,  3,  4,  5,  9,
    14, 15, 18, 13, 12, 11, 10, 7,  0,  1,  3,  5,  6,  8,  15, 16, 19, 16, 15, 11, 9,  4,  0,  1,
    2,  3,  5,  10, 13, 17, 19, 18, 17, 16, 15, 14, 11, 10, 7,  5,  4,  1,  3,  6,  8,  9,  19, 15,
    14, 13, 10, 8,  6,  4,  2,  5,  7,  11, 12, 18, 17, 14, 13, 12, 9,  6,  5,  2,  1,  4,  7,  8,
    11, 16, 13, 12, 11, 9,  8,  7,  6,  0,  1,  10,
};
static const uint16_t offsets_165_4125[] = {
    0,   10,  28,  45,  61,  77,  93,  108, 121, 135, 154, 171, 185, 199,
    214, 229, 246, 258, 273, 288, 302, 314, 329, 342, 358, 372, 385, 396,
};
static const uint32_t powers_205_5125[] = {
    128, 81,  125, 49,  121, 169, 17,  19,  23,  29,  31,  37,  41,  43,  47,  53,
    59,  61,  67,  71,  73,  79,  83,  89,  97,  101, 103, 107, 109, 113, 127, 131,
    137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199,
};
static const uint8_t slots_205_5125[] = {
    255, 0,   255, 255, 255, 255, 255, 1,   255, 255, 255, 2,   255, 3,   255, 255, 255, 4,
    255, 5,   255, 255, 255, 6,   255, 255, 255, 255, 255, 7,   255, 8,   255, 255, 255, 255,
    255, 9,   255, 255, 255, 10,  255, 11,  255, 255, 255, 12,  255, 13,  255, 255, 255, 14,
    255, 255, 255, 255, 255, 15,  255, 16,  255, 255, 255, 255, 255, 17,  255, 255, 255, 18,
    255, 19,  255, 255, 255, 20,  255, 21,  255, 255, 255, 22,  255, 255, 255, 255, 255, 23,
    255, 24,  255, 255, 255, 255, 255, 25,  255, 255, 255, 26,  255, 27,  255, 255,
};
static const uint8_t indices_205_5125[] = {
    0,  3,  4,  5,  6,  7,  8,  10, 12, 14, 15, 16, 17, 18, 19, 22, 25, 26, 27, 27, 23, 22, 19, 18,
    17, 16, 14, 12, 10, 9,  8,  6,  5,  2,  0,  3,  7,  11, 15, 21, 26, 23, 22, 19, 17, 16, 15, 14,
    11, 9,  8,  7,  6,  4,  3,  2,  0,  12, 18, 21, 25, 27, 26, 25, 23, 22, 21, 18, 17, 14, 11, 8,
    7,  5,  4,  3,  2,  0,  6,  9,  10, 12, 27, 25, 22, 21, 19, 17, 15, 14, 10, 9,  8,  7,  5,  4,
    2,  0,  3,  11, 12, 26, 25, 23, 21, 19, 17, 15, 12, 11, 9,  8,  7,  6,  2,  0,  4,  5,  10, 16,
    26, 27, 25, 23, 18, 16, 12, 11, 10, 9,  8,  6,  5,  4,  2,  0,  3,  7,  14, 19, 21, 22, 26, 26,
    25, 22, 21, 19, 18, 17, 16, 15, 14, 11, 6,  4,  3,  2,  5,  7,  10, 27, 27, 26, 23, 21, 17, 15,
    11, 7,  6,  5,  4,  3,  2,  0,  10, 16, 22, 25, 27, 26, 25, 23, 22, 19, 18, 16, 12, 9,  8,  5,
    4,  3,  2,  0,  7,  10, 11, 14, 21, 27, 25, 23, 19, 18, 17, 15, 11, 10, 9,  7,  6,  4,  3,  0,
    8,  12, 16, 21, 22, 26, 27, 25, 22, 21, 19, 16, 14, 12, 11, 4,  0,  2,  5,  6,  7,  8,  9,  15,
    18, 23, 26, 25, 22, 19, 18, 17, 15, 14, 12, 11, 10, 9,  8,  6,  5,  4,  2,  0,  16, 23, 27, 27,
    25, 23, 22, 21, 16, 14, 11, 9,  8,  6,  3,  0,  4,  7,  15, 18, 26, 26, 23, 22, 18, 17, 16, 10,
    8,  7,  3,  4,  5,  9,  14, 15, 21, 27, 27, 26, 23, 16, 15, 14, 12, 10, 9,  8,  7,  4,  3,  0,
    2,  19, 25, 27, 26, 21, 18, 15, 14, 11, 10, 9,  8,  7,  6,  3,  2,  0,  12, 16, 17, 19, 23, 27,
    23, 22, 21, 18, 16, 14, 12, 10, 5,  3,  2,  0,  4,  6,  11, 17, 19, 25, 26, 26, 22, 21, 19, 18,
    17, 16, 15, 12, 11, 6,  0,  2,  3,  4,  7,  8,  9,  23, 27, 26, 23, 19, 18, 17, 16, 12, 11, 10,
    6,  0,  2,  4,  5,  7,  8,  14, 15, 22, 25, 22, 19, 18, 16, 14, 12, 9,  5,  3,  0,  2,  8,  10,
    25, 27, 27, 26, 25, 19, 18, 15, 14, 9,  7,  6,  4,  0,  5,  8,  11, 22, 26, 25, 21, 18, 12, 11,
    10, 9,  8,  7,  4,  3,  0,  15, 19, 23, 27, 27, 25, 23, 22, 19, 18, 17, 14, 12, 10, 9,  8,  7,
    5,  4,  0,  2,  15, 16, 21,
};
static const uint16_t offsets_205_5125[] = {
    0,   19,  41,  62,  82,  102, 121, 143, 162, 180, 201, 222, 243,
    263, 281, 298, 315, 335, 355, 375, 395, 410, 426, 443, 463,
};
/* The optional tiny-factor probe is scalar and uses the binary ladder.
 * Its separate plan leaves main-plan lookup unchanged by this extra bound. */
static const uint32_t prefix15_powers[] = {
    8, 9, 5, 7, 11, 13,
};
static const uint8_t prefix15_slots[] = {
    255, 0, 255, 255, 255, 255, 255, 1, 255, 255, 255, 2, 255, 3, 255, 255,
};
static const uint8_t prefix15_indices[] = {
    3, 2, 1, 0, 3, 1, 0, 2, 2, 1, 0, 3, 3, 2, 1, 3, 2, 0, 1, 3, 1, 0, 2, 3,
    2, 0, 3, 2, 1, 0, 3, 1, 0, 2, 1, 2, 3, 3, 0, 1, 3, 2, 1, 0,
};
static const uint16_t prefix15_offsets[] = {
    0, 4, 8, 12, 15, 19, 23, 26, 30, 34, 37, 40, 44,
};
static const ecm64_plan prefix15 = {
    15, 375, 6, 15, 12, 4,
    prefix15_powers, prefix15_slots, prefix15_indices, prefix15_offsets,
};

/* Small main plan: extend the powers of 2, 3 and 5 while omitting
 * the larger stage-1 primes. B2 remains 25*B1. */
static const uint32_t low16_powers[] = {
    32, 27, 25, 7, 11, 13,
};
static const uint8_t low16_slots[] = {
    255, 0, 255, 255, 255, 255, 255, 1, 255, 255, 255, 2, 255, 3, 255, 255,
};
static const uint8_t low16_indices[] = {
    3, 2, 1, 0, 3, 1, 0, 2, 2, 1, 0, 3, 3, 2, 1, 3,
    2, 0, 1, 3, 1, 0, 2, 3, 2, 0, 3, 2, 1, 0, 3, 1,
    0, 2, 1, 2, 3, 3, 0, 1, 3, 2, 1, 0, 2, 1, 0,
};
static const uint16_t low16_offsets[] = {
    0, 4, 8, 12, 15, 19, 23, 26, 30, 34, 37, 40, 44, 47,
};

/* Small main plan: extend the powers of 2, 3 and 5 while omitting
 * the larger stage-1 primes. B2 remains 25*B1. */
static const uint32_t low19_powers[] = {
    32, 27, 25, 7, 11, 13, 17, 19,
};
static const uint8_t low19_slots[] = {
    255, 0, 255, 255, 255, 255, 255, 1, 255, 255, 255, 2, 255, 3, 255, 255,
};
static const uint8_t low19_indices[] = {
    1, 0, 2, 3, 3, 1, 0, 2, 2, 1, 0, 3, 3, 2, 1, 3,
    2, 0, 1, 3, 1, 0, 2, 3, 2, 0, 3, 2, 1, 0, 3, 1,
    0, 2, 1, 2, 3, 3, 0, 1, 3, 2, 1, 0, 2, 1, 0, 2,
    0, 3, 2, 1, 0, 3, 3,
};
static const uint16_t low19_offsets[] = {
    0, 4, 8, 12, 15, 19, 23, 26, 30, 34, 37, 40, 44, 47, 50, 54,
    55,
};

static const ecm64_plan plans[] = {
    /* Reuse the existing small-prefix schedule for the lowest main rung. */
    {15,375,6,15,12,4,prefix15_powers,prefix15_slots,prefix15_indices,prefix15_offsets},
    {16,400,6,15,13,4,low16_powers,low16_slots,low16_indices,low16_offsets},
    {19,475,8,15,16,4,low19_powers,low19_slots,low19_indices,low19_offsets},

    {27, 675, 9, 15, 22, 4, powers_27_675, slots_27_675, indices_27_675, offsets_27_675},
    {70, 1750, 19, 45, 19, 12, powers_70_1750, slots_70_1750, indices_70_1750, offsets_70_1750},
    {125, 3125, 30, 75, 21, 20, powers_125_3125, slots_125_3125, indices_125_3125,
     offsets_125_3125},
    {165, 4125, 38, 75, 27, 20, powers_165_4125, slots_165_4125, indices_165_4125,
     offsets_165_4125},
    {205, 5125, 46, 105, 24, 28, powers_205_5125, slots_205_5125, indices_205_5125,
     offsets_205_5125},
};
#define ECM64_WORK_POINTS 52

typedef struct {
  ecm_point nq[ECM64_WORK_POINTS];
} scalar_workspace;

/* Stage 2 uses projective cross-products, avoiding affine inversion. Check all
 * denominators first: a nonunit Z may already expose a factor. When products
 * hide different factors, scan individual Z values before abandoning a curve. */
static int check_point_denominators(uint64_t *factor, const ecm_point *points,
                                     unsigned count, const mont64_t *ctx) {
  unsigned i;
  uint64_t product = ctx->one, g;
  *factor = 0;
  for (i = 0; i < count; i++) product = multiply(product, points[i].Z, ctx);
  g = gcd64(product, ctx->n);
  if (g == 1) return 1;
  if (g < ctx->n) *factor = g;
  else for (i = 0; i < count; i++) {
    g = gcd64(points[i].Z, ctx->n);
    if (g > 1 && g < ctx->n) { *factor = g; break; }
  }
  return 0;
}

static const ecm64_plan *find_plan(uint32_t b1, uint32_t b2) {
  size_t i;
  for (i = 0; i < sizeof(plans) / sizeof(plans[0]); i++)
    if (plans[i].b1 == b1 && plans[i].b2 == b2) return &plans[i];
  return NULL;
}

static void scalar_geometry(const ecm_point *q, uint64_t a24, const mont64_t *ctx,
                            const ecm64_plan *p, scalar_workspace *ws) {
  unsigned i, w;
  ecm_point q2, q3, q5, q6, q7, step, current, previous, next;
  ecm_point *babies = ws->nq, *giants = babies + p->babies;
  point_double(&q2, q, a24, ctx);
  point_add(&q3, &q2, q, q, ctx);
  point_double(&q6, &q3, a24, ctx);
  point_add(&q5, &q3, &q2, q, ctx);
  point_add(&q7, &q6, q, &q5, ctx);
  babies[p->slots[1]] = *q;
  previous = *q;
  current = q7;
  for (i = 7; i <= p->d; i += 6) {
    if (i % 5) babies[p->slots[i]] = current;
    if (i + 6 <= p->d) {
      point_add(&next, &current, &q6, &previous, ctx);
      previous = current;
      current = next;
    }
  }
  previous = *q;
  current = q5;
  for (i = 5; i <= p->d; i += 6) {
    if (i % 5) babies[p->slots[i]] = current;
    if (i + 6 <= p->d) {
      point_add(&next, &current, &q6, &previous, ctx);
      previous = current;
      current = next;
    }
  }
  point_multiply(&step, &q6, 2 * p->d / 6, a24, ctx);
  giants[0] = step;
  if (p->windows > 1) point_double(&giants[1], &step, a24, ctx);
  for (w = 2; w < p->windows; w++)
    point_add(&giants[w], &step, &giants[w - 1], &giants[w - 2], ctx);
}

/* Xg/Zg - Xb/Zb has the same factor GCD as Xg*Zb - Xb*Zg
 * after the denominator check. Products stay in the Montgomery domain. */
static INLINE uint64_t stage2_diff(const ecm_point *giant, const ecm_point *baby,
                                   const mont64_t *ctx) {
  return submod(multiply(giant->X, baby->Z, ctx),
                multiply(baby->X, giant->Z, ctx), ctx->n);
}

static uint64_t scalar_recover(const ecm64_plan *p, const ecm_point *points,
                               unsigned first, unsigned last, const mont64_t *ctx) {
  unsigned w, i;
  for (w = first; w <= last; w++)
    for (i = p->offsets[w]; i < p->offsets[w + 1]; i++) {
      uint64_t f = gcd64(stage2_diff(points + p->babies + w, points + p->indices[i], ctx), ctx->n);
      if (f > 1 && f < ctx->n) return f;
    }
  return 0;
}

static uint64_t scalar_stage2(const ecm_point *q, uint64_t a24, const mont64_t *ctx,
                              const ecm64_plan *p, scalar_workspace *ws) {
  unsigned w, i;
  const uint64_t n = ctx->n;
  ecm_point *babies = ws->nq;
  uint64_t f;

  scalar_geometry(q, a24, ctx, p, ws);

  if (!check_point_denominators(&f, babies, p->babies + p->windows, ctx)) return f;

  uint64_t product = ctx->one;
  unsigned checked_window = 0;
  for (w = 0; w < p->windows; w++) {
    for (i = p->offsets[w]; i < p->offsets[w + 1]; i++) {
      uint64_t diff = stage2_diff(babies + p->babies + w, babies + p->indices[i], ctx);
      product = multiply(product, diff, ctx);
    }
    if ((((w) + 1) % ECM64_GCD_WINDOWS == 0 || (w) + 1 == (p)->windows) &&
        p->offsets[w + 1] > p->offsets[checked_window]) {
      f = gcd_residue(product, ctx);
      if (f > 1) {
        if (f < n) return f;
        return scalar_recover(p, babies, checked_window, w, ctx);
      }
      checked_window = w + 1;
    }
  }
  /* The final nonempty group is checked even when shorter than the cadence. */
  return 0;
}

static void setup_terms(const mont64_t *ctx, uint16_t sigma, ecm_point *point, uint64_t *num,
                        uint64_t *den) {
  unsigned i;
  uint64_t n = ctx->n, ui = (uint64_t)sigma * sigma - 5, vi = 4 * (uint64_t)sigma;
  uint64_t u = mont64_enter(((ui) % (n)), ctx), v = mont64_enter(((vi) % (n)), ctx);
  uint64_t delta = mont64_enter(((ui - vi) % (n)), ctx);
  uint64_t u2 = multiply(u, u, ctx), v2 = multiply(v, v, ctx);
  point->X = multiply(u2, u, ctx);
  point->Z = multiply(v2, v, ctx);
  uint64_t d2 = multiply(delta, delta, ctx), d3 = multiply(d2, delta, ctx);
  *num = multiply(d3, mont64_enter(((3 * ui + vi) % (n)), ctx), ctx);
  *num = *num ? n - *num : 0;
  *den = multiply(point->X, v, ctx);
  for (i = 0; i < 4; i++)
    *den = mont64_add(*den, *den, n);
}

static uint64_t setup_curve(const mont64_t *ctx, uint16_t sigma, ecm_point *point, uint64_t *a24,
                            int *ready) {
  uint64_t numerator, denominator;
  setup_terms(ctx, sigma, point, &numerator, &denominator);
  denominator = mont64_exit(denominator, ctx);
  uint64_t inverse = inverse64(denominator, ctx->n);
  *ready = 0;
  if (!inverse) {
    uint64_t f = gcd64(denominator, ctx->n);
    return f > 1 && f < ctx->n ? f : 0;
  }
  *a24 = multiply(numerator, mont64_enter(inverse, ctx), ctx);
  *ready = 1;
  return 0;
}

static uint64_t scalar_curve(const mont64_t *ctx, const ecm64_plan *p, uint16_t sigma,
                             scalar_workspace *ws, int use_prac) {
  uint32_t j;

  const uint64_t n = ctx->n;
  uint64_t a24;
  ecm_point point;
  int ready;
  uint64_t initial = setup_curve(ctx, sigma, &point, &a24, &ready);
  if (initial || !ready) return initial;

  for (j = 0; j < p->npowers; j++) {
    if (use_prac)
      scalar_stage1_prac(&point, &point, p->powers[j], a24, ctx);
    else
      scalar_stage1(&point, &point, p->powers[j], a24, ctx);
    if ((j % 64) == 0) {
      uint64_t f = gcd_residue(point.Z, ctx);
      if (f > 1) return f < n ? f : 0;
    }
  }
  uint64_t f = gcd_residue(point.Z, ctx);
  if (f > 1 && f < n) return f;
  return scalar_stage2(&point, a24, ctx, p, ws);
}

static uint64_t splitmix(uint64_t *state) {
  uint64_t z = (*state += UINT64_C(0x9e3779b97f4a7c15));
  z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
  return z ^ (z >> 31);
}

static uint16_t next_sigma(uint16_t *sigmas, unsigned i, uint64_t *state) {
  unsigned j;
  uint16_t s;
  int duplicate;
  do {
    s = (uint16_t)(6 + splitmix(state) % 65530);
    duplicate = 0;
    for (j = 0; j < i; j++)
      if (sigmas[j] == s) {
        duplicate = 1;
        break;
      }
  } while (duplicate);
  sigmas[i] = s;
  return s;
}

#if MONT64_HAVE_UINT128
/* Four independent curves in structure-of-arrays layout. */

typedef struct {
  uint64_t n, inv;
} four_modulus;
typedef struct {
  uint64_t X[ECM64_LANES], Z[ECM64_LANES];
} four_point;
ECM64_INLINE uint64_t high_product(uint64_t a, uint64_t b) {
#if defined(__x86_64__) || defined(_M_X64)
  uint64_t a0 = (uint32_t)a, a1 = a >> 32, b0 = (uint32_t)b, b1 = b >> 32;
  uint64_t w0 = a0 * b0, t = a1 * b0 + (w0 >> 32);
  uint64_t w1 = (uint32_t)t, w2 = t >> 32;
  w1 += a0 * b1;
  return a1 * b1 + w2 + (w1 >> 32);
#else
  return (uint64_t)(((uint128_t)a * b) >> 64);
#endif
}
ECM64_INLINE uint64_t four_multiply(uint64_t a, uint64_t b, const four_modulus *ctx) {

  uint64_t high = high_product(a, b), m = (a * b) * ctx->inv;
  uint64_t q = high_product(m, ctx->n), r = high - q;
  return high < q ? r + ctx->n : r;
}
ECM64_INLINE uint64_t four_addmod(uint64_t a, uint64_t b, uint64_t n) {
  uint64_t gap = n - b;
  return a >= gap ? a - gap : a + b;
}
ECM64_INLINE uint64_t four_submod(uint64_t a, uint64_t b, uint64_t n) {
  return a >= b ? a - b : n - (b - a);
}
ECM64_INLINE void four_double(four_point *R, const four_point *P, const uint64_t *a24,
                              const four_modulus *ctx) {
  unsigned i;

  four_point result;
  for (i = 0; i < ECM64_LANES; i++) {
    uint64_t x = P->X[i], z = P->Z[i], a = four_submod(x, z, ctx->n), b = four_addmod(x, z, ctx->n);
    uint64_t u = four_multiply(a, a, ctx), v = four_multiply(b, b, ctx),
             w = four_submod(v, u, ctx->n);
    result.X[i] = four_multiply(u, v, ctx);
    result.Z[i] = four_multiply(w, four_addmod(u, four_multiply(a24[i], w, ctx), ctx->n), ctx);
  }
  *R = result; /* Local output preserves all permitted in-place aliases. */
}
ECM64_INLINE void four_add(four_point *R, const four_point *P, const four_point *Q,
                           const four_point *D, const four_modulus *ctx) {
  unsigned i;

  four_point result;
  for (i = 0; i < ECM64_LANES; i++) {
    uint64_t u = four_multiply(four_submod(P->X[i], P->Z[i], ctx->n),
                               four_addmod(Q->X[i], Q->Z[i], ctx->n), ctx);
    uint64_t v = four_multiply(four_addmod(P->X[i], P->Z[i], ctx->n),
                               four_submod(Q->X[i], Q->Z[i], ctx->n), ctx);
    uint64_t sum = four_addmod(u, v, ctx->n), diff = four_submod(u, v, ctx->n);
    result.X[i] = four_multiply(four_multiply(sum, sum, ctx), D->Z[i], ctx);
    result.Z[i] = four_multiply(four_multiply(diff, diff, ctx), D->X[i], ctx);
  }
  *R = result;
}
static void four_ladder(four_point *R, const four_point *P, uint64_t k, const uint64_t *a24,
                        const four_modulus *ctx) {
  four_point r0 = *P, r1;
  uint64_t bit;
  if (k == 1) {
    *R = *P;
    return;
  }
  four_double(&r1, &r0, a24, ctx);
  if (k == 2) {
    *R = r1;
    return;
  }
  bit = leading_bit(k);
  for (bit >>= 1; bit; bit >>= 1) {
    if (k & bit) {
      four_add(&r0, &r0, &r1, P, ctx);
      four_double(&r1, &r1, a24, ctx);
    } else {
      four_add(&r1, &r0, &r1, P, ctx);
      four_double(&r0, &r0, a24, ctx);
    }
  }
  *R = r0;
}

static void four_stage1_clean(four_point *out, const four_point *in, uint32_t k,
                              const uint64_t *a24, const four_modulus *ctx) {
  four_point a = *in, b;
  /* Binary scalar ladder, omitting the unused final register update. */
  uint64_t bit;
  if (!(k & (k - 1))) {
    for (; k > 1; k >>= 1)
      four_double(&a, &a, a24, ctx);
    *out = a;
    return;
  }
  four_double(&b, &a, a24, ctx);
  bit = leading_bit(k);
  for (bit >>= 1; bit > 1; bit >>= 1) {
    if (k & bit) {
      four_add(&a, &a, &b, in, ctx);
      four_double(&b, &b, a24, ctx);
    } else {
      four_add(&b, &a, &b, in, ctx);
      four_double(&a, &a, a24, ctx);
    }
  }
  /* Only the returned register is needed for the final bit. */
  if (k & 1)
    four_add(out, &a, &b, in, ctx);
  else
    four_double(out, &a, a24, ctx);
}

static void four_stage1(four_point *out, const four_point *in, uint32_t k, const uint64_t *a24,
                        const four_modulus *ctx) {
  if (!(k & 1) || k == 1) {
    four_stage1_clean(out, in, k, a24, ctx);
    return;
  }
  four_point a, b = *in, c = *in, t, u;
  /* k comes from the immutable plans; every such odd multiplier has a
   * verified nonempty word. Its highest unused bits terminate decoding. */
  uint32_t word = chain_words[k >> 1];
  four_double(&a, &b, a24, ctx);
  for (; word; word >>= 3) {
    switch (word & 7) {
      /* Three-bit opcode, with any swap folded into the operation. */
    case 1:
      four_add(&t, &a, &b, &c, ctx);
      four_add(&u, &t, &a, &b, ctx);
      four_add(&b, &b, &t, &a, ctx);
      a = u;
      break;
    case 2:
      four_add(&c, &b, &a, &c, ctx);
      t = b;
      b = c;
      c = t;
      break;
    case 3:
      four_add(&b, &b, &a, &c, ctx);
      four_double(&a, &a, a24, ctx);
      break;
    case 4:
      four_add(&c, &c, &a, &b, ctx);
      four_double(&a, &a, a24, ctx);
      break;
    case 5:
      t = a;
      a = b;
      b = t;
      four_add(&t, &a, &b, &c, ctx);
      four_add(&u, &t, &a, &b, ctx);
      four_add(&b, &b, &t, &a, ctx);
      a = u;
      break;
    case 6:
      t = a;
      a = b;
      b = t;
      four_add(&c, &b, &a, &c, ctx);
      t = b;
      b = c;
      c = t;
      break;
    case 7:
      t = a;
      a = b;
      b = t;
      four_add(&b, &b, &a, &c, ctx);
      four_double(&a, &a, a24, ctx);
      break;
    }
  }
  four_add(out, &a, &b, &c, ctx);
}

typedef struct {
  four_point points[ECM64_WORK_POINTS];
} four_workspace;

static void four_geometry(const four_point *q, const uint64_t *a24, const four_modulus *ctx,
                          const ecm64_plan *p, four_workspace *ws) {
  unsigned i, w;
  four_point q2, q3, q5, q6, q7, step, current, previous, next;
  four_point *babies = ws->points, *giants = babies + p->babies;
  four_double(&q2, q, a24, ctx);
  four_add(&q3, &q2, q, q, ctx);
  four_double(&q6, &q3, a24, ctx);
  four_add(&q5, &q3, &q2, q, ctx);
  four_add(&q7, &q6, q, &q5, ctx);
  babies[p->slots[1]] = *q;
  previous = *q;
  current = q7;
  for (i = 7; i <= p->d; i += 6) {
    if (i % 5) babies[p->slots[i]] = current;
    if (i + 6 <= p->d) {
      four_add(&next, &current, &q6, &previous, ctx);
      previous = current;
      current = next;
    }
  }
  previous = *q;
  current = q5;
  for (i = 5; i <= p->d; i += 6) {
    if (i % 5) babies[p->slots[i]] = current;
    if (i + 6 <= p->d) {
      four_add(&next, &current, &q6, &previous, ctx);
      previous = current;
      current = next;
    }
  }
  four_ladder(&step, &q6, 2 * p->d / 6, a24, ctx);
  giants[0] = step;
  if (p->windows > 1) four_double(&giants[1], &step, a24, ctx);
  for (w = 2; w < p->windows; w++)
    four_add(&giants[w], &step, &giants[w - 1], &giants[w - 2], ctx);
}

/* Input and output arrays are distinct; all four values share one modulus. */
static int invert_four(uint64_t out[4], const uint64_t in[4], const mont64_t *ctx) {
  unsigned i;
  uint64_t prefix[4], inv;
  prefix[0] = in[0];
  for (i = 1; i < 4; i++)
    prefix[i] = multiply(prefix[i - 1], in[i], ctx);
  inv = inverse64(mont64_exit(prefix[3], ctx), ctx->n);
  if (!inv) return 0;
  inv = mont64_enter(inv, ctx);
  for (i = 3; i > 0; i--) {
    out[i] = multiply(inv, prefix[i - 1], ctx);
    inv = multiply(inv, in[i], ctx);
  }
  out[0] = inv;
  return 1;
}

static uint64_t setup_four(const mont64_t *ctx, const uint16_t *sigmas, four_point *point,
                           uint64_t a24[4], unsigned *active) {
  unsigned lane;
  uint64_t num[4], den[4], inv[4];

  for (lane = 0; lane < 4; lane++) {
    ecm_point p;
    setup_terms(ctx, sigmas[lane], &p, &num[lane], &den[lane]);
    point->X[lane] = p.X;
    point->Z[lane] = p.Z;
  }
  if (invert_four(inv, den, ctx)) {
    for (lane = 0; lane < 4; lane++)
      a24[lane] = multiply(num[lane], inv[lane], ctx);
    *active = 15;
    return 0;
  }
  /* A product can hide factors from different lanes, or one denominator
   * can be zero. Reuse the original setup/factor/rejected-curve handling. */
  *active = 0;
  for (lane = 0; lane < 4; lane++) {
    ecm_point p;
    int ready;
    uint64_t f = setup_curve(ctx, sigmas[lane], &p, &a24[lane], &ready);
    if (f) return f;
    if (ready)
      *active |= 1U << lane;
    else
      p.X = p.Z = a24[lane] = ctx->one;
    point->X[lane] = p.X;
    point->Z[lane] = p.Z;
  }
  return 0;
}

static uint64_t check_four_denominators(four_workspace *ws, unsigned count,
                                        const mont64_t *ctx, unsigned *active) {
  unsigned lane, i;
  four_modulus mod = {ctx->n, (uint64_t)0 - ctx->ninv};
  uint64_t acc[4] = {ctx->one, ctx->one, ctx->one, ctx->one};
  uint64_t product = ctx->one;
  for (lane = 0; lane < 4; lane++)
    if (!(*active & (1U << lane)))
      for (i = 0; i < count; i++)
        ws->points[i].X[lane] = ws->points[i].Z[lane] = ctx->one;
  for (i = 0; i < count; i++)
    for (lane = 0; lane < 4; lane++)
      acc[lane] = four_multiply(acc[lane], ws->points[i].Z[lane], &mod);
  for (lane = 0; lane < 4; lane++) product = multiply(product, acc[lane], ctx);
  if (gcd64(product, ctx->n) == 1) return 0;
  /* Match scalar recovery order when different lanes hide each other's factors. */
  for (lane = 0; lane < 4; lane++) if (*active & (1U << lane)) {
    uint64_t g = gcd64(acc[lane], ctx->n);
    if (g > 1 && g < ctx->n) return g;
    if (g == ctx->n) {
      for (i = 0; i < count; i++) {
        g = gcd64(ws->points[i].Z[lane], ctx->n);
        if (g > 1 && g < ctx->n) return g;
      }
      *active &= ~(1U << lane);
    }
  }
  return 0;
}

/* Combine active lanes; scan them separately if their factors hide each other. */
static uint64_t gcd_four(const uint64_t values[4], unsigned active, const mont64_t *ctx,
                         unsigned *lost) {
  unsigned lane;
  uint64_t prod = ctx->one, g;
  *lost = 0;
  for (lane = 0; lane < 4; lane++)
    if (active & (1U << lane)) prod = multiply(prod, values[lane], ctx);
  g = gcd_residue(prod, ctx);
  if (g > 1 && g < ctx->n) return g;
  if (g == ctx->n)
    for (lane = 0; lane < 4; lane++)
      if (active & (1U << lane)) {
        g = gcd_residue(values[lane], ctx);
        if (g > 1 && g < ctx->n) return g;
        if (g == ctx->n) *lost |= 1U << lane;
      }
  return 0;
}

static INLINE uint64_t stage2_diff_four(const four_workspace *ws, unsigned g, unsigned b,
                                         unsigned lane, const four_modulus *mod) {
  return four_submod(four_multiply(ws->points[g].X[lane], ws->points[b].Z[lane], mod), four_multiply(ws->points[b].X[lane], ws->points[g].Z[lane], mod), mod->n);
}

static uint64_t recover_lane(const ecm64_plan *p, unsigned first, unsigned last,
                             const four_workspace *ws, unsigned lane, const mont64_t *ctx) {
  unsigned w, i;
  four_modulus mod = {ctx->n, (uint64_t)0 - ctx->ninv};
  for (w = first; w <= last; w++)
    for (i = p->offsets[w]; i < p->offsets[w + 1]; i++) {
      uint64_t diff = stage2_diff_four(ws, p->babies + w, p->indices[i], lane, &mod);
      uint64_t f = gcd64(diff, ctx->n);
      if (f > 1 && f < ctx->n) return f;
    }
  return 0;
}

static uint64_t four_curve(const mont64_t *ctx, const ecm64_plan *p, const uint16_t *sigmas) {
  uint32_t j;
  unsigned w, i, lane;
  four_point point;
  uint64_t a24[4], f, n = ctx->n;
  unsigned active = 0;
  /* Factors in the first stage-1 Z persist through later multiplications.
   * Save it instead of taking an early GCD; recover it if the final product
   * hides factors from different primes. Scalar prefixes keep their early check. */
  uint64_t first_z[4] = {ctx->one, ctx->one, ctx->one, ctx->one};
  four_modulus mod = {n, (uint64_t)0 - ctx->ninv};
  f = setup_four(ctx, sigmas, &point, a24, &active);
  if (f) return f;
  if (!active) return 0;

  for (j = 0; j < p->npowers; j++) {
    four_stage1(&point, &point, p->powers[j], a24, &mod);
    if (j == 0) memcpy(first_z, point.Z, sizeof(first_z));
    else if ((j % 64) == 0) {
      unsigned lost;
      f = gcd_four(point.Z, active, ctx, &lost);
      if (f) return f;
      if (lost) {
        unsigned ignored;
        f = gcd_four(first_z, active & lost, ctx, &ignored);
        if (f) return f;
      }
      active &= ~lost;
    }
    if (!active) return 0;
  }
  unsigned stage1_lost;
  f = gcd_four(point.Z, active, ctx, &stage1_lost);
  if (f) return f;
  if (stage1_lost) {
    unsigned ignored;
    f = gcd_four(first_z, active & stage1_lost, ctx, &ignored);
    if (f) return f;
  }
  four_workspace ws;

  four_geometry(&point, a24, &mod, p, &ws);

  unsigned points = p->babies + p->windows;
  f = check_four_denominators(&ws, points, ctx, &active);
  if (f) return f;
  if (!active) return 0;

  uint64_t product[4] = {ctx->one, ctx->one, ctx->one, ctx->one};
  unsigned checked_window = 0;
  for (w = 0; w < p->windows; w++) {
    for (i = p->offsets[w]; i < p->offsets[w + 1]; i++)
      for (lane = 0; lane < 4; lane++) {
        uint64_t diff = stage2_diff_four(&ws, p->babies + w, p->indices[i], lane, &mod);
        product[lane] = four_multiply(product[lane], diff, &mod);
      }
    if ((((w) + 1) % ECM64_GCD_WINDOWS == 0 || (w) + 1 == (p)->windows) &&
        p->offsets[w + 1] > p->offsets[checked_window]) {
      unsigned lost;
      f = gcd_four(product, active, ctx, &lost);
      if (f) return f;
      for (lane = 0; lane < 4; lane++)
        if (lost & (1U << lane)) {
          f = recover_lane(p, checked_window, w, &ws, lane, ctx);
          if (f) return f;
          active &= ~(1U << lane);
        }
      checked_window = w + 1;
    }
    if (!active) return 0;
  }
  return 0;
}

#endif

/* Select one main schedule by input size. try_small prepends its extra
 * curves independently of this budget. A main miss gets fresh four-curve batches
 * at the next original B1 setting (at least 70), or at the same bounds
 * for the final setting. New low rungs retain this conservative rescue.
 * Try the second batch only if the first one fails.
 * All main, prefix and rescue plans use B2 = 25 * B1. */
static const struct {
  unsigned max_bits, b1, curves;
} policy[] = {{28, 15, 44}, {37, 16, 44}, {40, 19, 44}, {44, 27, 32}, {53, 70, 32}, {58, 125, 32},
              {60, 165, 32}, {64, 205, 46}};
#if defined(__x86_64__) || defined(_M_X64)
enum { FOUR_CURVE_MIN_BITS = 40 };
#else
enum { FOUR_CURVE_MIN_BITS = 34 };
#endif
enum { MAX_MAIN_CURVES = 46, MAX_PREFIX_CURVES = 2, RESCUE_CURVES = 8 };

/* Keep the uncommon continuation out of the main routine. Copy its short
 * sigma history here so ordinary calls need no extra rescue workspace.
 * Include prefix sigmas and continue the existing RNG state: every rescue
 * curve is new. In particular, append after all 46 main curves without
 * regrouping their existing final two scalar curves. */
static NOINLINE uint64_t rescue_curves(const mont64_t *ctx, unsigned b1,
                                      const uint16_t *seen, unsigned nseen,
                                      uint64_t *state, int four) {
  unsigned j;
  unsigned next_b1 = b1;
  for (j = 0; j < sizeof(policy) / sizeof(policy[0]); j++) {
    if (policy[j].b1 >= 70 && policy[j].b1 > b1) {
      next_b1 = policy[j].b1;
      break;
    }
  }
  const ecm64_plan *plan = find_plan(next_b1, 25 * next_b1);
  uint16_t sigmas[MAX_MAIN_CURVES + MAX_PREFIX_CURVES + RESCUE_CURVES];
  memcpy(sigmas, seen, nseen * sizeof(*sigmas));
#if MONT64_HAVE_UINT128
  if (four) {
    unsigned i;
    for (i = 0; i < RESCUE_CURVES; i += ECM64_LANES) {
      for (j = 0; j < ECM64_LANES; j++)
        next_sigma(sigmas, nseen + i + j, state);
      uint64_t f = four_curve(ctx, plan, sigmas + nseen + i);
      if (f) return f;
    }
    return 0;
  }
#endif
  scalar_workspace ws;
  for (j = 0; j < RESCUE_CURVES; j++) {
    uint16_t sigma = next_sigma(sigmas, nseen + j, state);
    uint64_t f = scalar_curve(ctx, plan, sigma, &ws, four);
    if (f) return f;
  }
  return 0;
}

static uint64_t factor_with_policy(uint64_t n, unsigned b1, unsigned count,
                                   int try_small, int four, uint64_t *rng_state) {
  unsigned i;
  uint64_t state = splitmix(rng_state);
  uint16_t sigmas[MAX_MAIN_CURVES + MAX_PREFIX_CURVES];
  mont64_t ctx;
  mont64_init(&ctx, n);
  scalar_workspace ws;
  /* Probe cheaply for small factors before paying for a main batch. Above
   * B1=27, add one scalar B1=70 curve. Both are additional
   * to the main and rescue budgets; try_small=0 skips them entirely. */
  unsigned prefix = try_small ? 1 + (b1 > 27) : 0;
  for (i = 0; i < prefix; i++) {
    next_sigma(sigmas, i, &state);
    const ecm64_plan *small_plan = i == 0 ? &prefix15 : find_plan(70, 1750);
    uint64_t f = scalar_curve(&ctx, small_plan, sigmas[i], &ws, 0);
    if (f) return f;
  }
  const ecm64_plan *plan = find_plan(b1, 25 * b1);
  i = 0;
#if MONT64_HAVE_UINT128
  if (four)
    for (; i + ECM64_LANES <= count; i += ECM64_LANES) {
      unsigned lane;
      for (lane = 0; lane < ECM64_LANES; lane++)
        next_sigma(sigmas, prefix + i + lane, &state);
      uint64_t f = four_curve(&ctx, plan, sigmas + prefix + i);
      if (f) return f;
    }
#else
  four = 0;
#endif
  for (; i < count; i++) {
    next_sigma(sigmas, prefix + i, &state);
    /* A short four-curve tail retains its scalar PRAC stage-1 schedule. */
    uint64_t f = scalar_curve(&ctx, plan, sigmas[prefix + i], &ws, four);
    if (f) return f;
  }
  return rescue_curves(&ctx, b1, sigmas, prefix + count, &state, four);
}
uint64_t tinyecm64_rng(uint64_t n, int try_small, uint64_t *rng_state) {
  uint64_t v;
  if (n < 4) return 0;
  if (!(n & 1)) return 2;
  unsigned bits = 0;
  for (v = n; v; v >>= 1)
    bits++;
  unsigned rung = 0;
  while (bits > policy[rung].max_bits)
    rung++;
#if MONT64_HAVE_UINT128 && (defined(__x86_64__) || defined(_M_X64))
  /* The scalar x86 path favors the original B1=27 plan from 34 bits until
   * four curves start at 40. ARM retains the common low plans above. */
  unsigned b1 = policy[rung].b1, curves = policy[rung].curves;
  if (bits >= 34 && bits < 40) { b1 = 27; curves = 32; }
  return factor_with_policy(n, b1, curves, try_small,
                            bits >= FOUR_CURVE_MIN_BITS, rng_state);
#else
  return factor_with_policy(n, policy[rung].b1, policy[rung].curves, try_small,
                            bits >= FOUR_CURVE_MIN_BITS, rng_state);
#endif
}

uint64_t tinyecm64(uint64_t n, int try_small) {
  uint64_t state = n;
  return tinyecm64_rng(n, try_small, &state);
}

#endif /* HAVE_UINT64 */
