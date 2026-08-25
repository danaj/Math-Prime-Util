#ifndef MPU_UTIL_MATH_H
#define MPU_UTIL_MATH_H

#include "ptypes.h"
#include <math.h>

/* Inputs are assumed to be unsigned. */
#define is_divis_2_3(n)     ( (n)%2 == 0 || (n) % 3 == 0 )

#if defined(__arm64__)
#define is_divis_2_3_5(n)   ( (n)%2 == 0 || (0x1669>>((n)%15))&1 )
#else
#define is_divis_2_3_5(n)   ( (n)%2 == 0 || (n) % 3 == 0 || (n) % 5 == 0 )
#endif
/* 2,3,5 could use the single test:  (0x1f75d77d >> (n % 30)) & 1 */

#define is_divis_2_3_5_7(n) ( is_divis_2_3_5(n) || (n) % 7 == 0 )

static MAYBE_UNUSED UV ipow(UV n, UV k) {
  UV p = 1;
  while (k) {
    if (k & 1) p *= n;
    k >>= 1;
    if (k)     n *= n;
  }
  return p;
}

/* If we have a very fast ctz, then use the fast FLINT version of gcd. */
#if defined(__GNUC__) && 100*__GNUC__ + __GNUC_MINOR >= 304
#include "util.h"
#define gcd_ui(x,y) gcdz(x,y)
#else
static MAYBE_UNUSED UV gcd_ui(UV x, UV y) {
  UV t;
  if (y > x) { t = x; x = y; y = t; }
  while (y > 0) {
    t = y;  y = x % y;  x = t;  /* y1 <- x0 % y0 ; x1 <- y0 */
  }
  return x;
}
#endif

static MAYBE_UNUSED UV lcm_ui(UV x, UV y) {
  /* Can overflow if lcm(x,y) > UV_MAX. */
  if (x == 0 || y == 0) return 0;
  return x * (y / gcd_ui(x,y));
}

/* Correct for all 64-bit inputs and all FP rounding modes. */
static MAYBE_UNUSED uint32_t isqrt(UV n) {
  /* The small addition means we only need to check for fixing downwards. */
  IV r = (IV)(sqrt((double)n) + 1e-6f);
  IV diff = n - (UV)r*r;
  return (uint32_t)(r - (diff < 0));
}

static MAYBE_UNUSED bool is_perfect_square_ret(UV n, uint32_t *root)
{
  uint32_t r;
  /* Fast filters reject 95.0% of non-squares. */
#if BITS_PER_WORD == 64
  if ((UVCONST(1) << (n&63)) & UVCONST(0xfdfdfdedfdfcfdec)) return 0;
  /* if ((UVCONST(1) << (n%45)) & UVCONST(0xfffffeeb7df6f9ec)) return 0; */
#else
  /* uint32_t m; */
  if ((1U << (n&31)) & 0xfdfcfdec) return 0;
  /* m = n % 105; if ((m*0xd24554cd) & (m*0x0929579a) & 0x38020141) return 0; */
#endif
  r = isqrt(n);
  if (root != 0) *root = r;
  return ((UV)r*r == n);
}
#define is_perfect_square(n)  is_perfect_square_ret(n,0)

#endif
