#ifndef MPU_UTIL_H
#define MPU_UTIL_H

#include "ptypes.h"

extern int  _XS_get_verbose(void);
extern void _XS_set_verbose(int v);
extern int  _XS_get_callgmp(void);
extern void _XS_set_callgmp(int v);

extern int _XS_is_prime(UV x);
extern UV  _XS_next_prime(UV x);
extern UV  _XS_prev_prime(UV x);

extern UV  _XS_prime_count(UV low, UV high);
extern UV  _XS_nth_prime(UV x);

extern signed char* _moebius_range(UV low, UV high);
extern UV*    _totient_range(UV low, UV high);
extern IV     _XS_mertens(UV n);
extern double _XS_chebyshev_theta(UV n);
extern double _XS_chebyshev_psi(UV n);

extern double _XS_ExponentialIntegral(double x);
extern double _XS_LogarithmicIntegral(double x);
extern long double ld_riemann_zeta(long double x);
extern double _XS_RiemannR(double x);

/* Above this value, is_prime will do deterministic Miller-Rabin */
/* With 64-bit math, we can do much faster mulmods from 2^16-2^32 */
#if (BITS_PER_WORD == 64) || HAVE_STD_U64
  #define MPU_PROB_PRIME_BEST  UVCONST(100000)
#else
  #define MPU_PROB_PRIME_BEST  UVCONST(100000000)
#endif

static UV isqrt(UV n) {
  UV root;
#if BITS_PER_WORD == 32
  if (n >= UVCONST(4294836225)) return UVCONST(65535);
#else
  if (n >= UVCONST(18446744065119617025)) return UVCONST(4294967295);
#endif
  root = (UV) sqrt((double)n);
  while (root*root > n)  root--;
  while ((root+1)*(root+1) <= n)  root++;
  return root;
}

#ifdef FUNC_gcd_ui
static UV gcd_ui(UV x, UV y) {
  UV t;
  if (y < x) { t = x; x = y; y = t; }
  while (y > 0) {
    t = y;  y = x % y;  x = t;  /* y1 <- x0 % y0 ; x1 <- y0 */
  }
  return x;
}
#endif

#endif
