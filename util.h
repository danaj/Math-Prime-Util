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
extern IV     mertens(UV n);
extern double _XS_chebyshev_theta(UV n);
extern double _XS_chebyshev_psi(UV n);

extern double _XS_ExponentialIntegral(double x);
extern double _XS_LogarithmicIntegral(double x);
extern long double ld_riemann_zeta(long double x);
extern double _XS_RiemannR(double x);
extern UV _XS_Inverse_Li(UV x);

extern int kronecker_uu(UV a, UV b);
extern int kronecker_su(IV a, UV b);
extern int kronecker_ss(IV a, IV b);

extern UV totient(UV n);
extern int moebius(UV n);
extern UV exp_mangoldt(UV n);
extern UV carmichael_lambda(UV n);
extern UV znprimroot(UV n);
extern UV znorder(UV a, UV n);

/* Above this value, is_prime will do deterministic Miller-Rabin */
/* With 64-bit math, we can do much faster mulmods from 2^16-2^32 */
#if (BITS_PER_WORD == 64) || HAVE_STD_U64
  #define MPU_PROB_PRIME_BEST  UVCONST(100000)
#else
  #define MPU_PROB_PRIME_BEST  UVCONST(100000000)
#endif

#ifdef FUNC_isqrt
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
#endif

#ifdef FUNC_icbrt
static UV icbrt(UV n) {
  UV b, root = 0;
#if BITS_PER_WORD == 32
  int s = 30;
  if (n >= UVCONST(4291015625)) return UVCONST(1625);
#else
  int s = 63;
  if (n >= UVCONST(18446724184312856125)) return UVCONST(2642245);
#endif
  for ( ; s >= 0; s -= 3) {
    root += root;
    b = 3*root*(root+1)+1;
    if ((n >> s) >= b) {
      n -= b << s;
      root++;
    }
  }
  return root;
}
#endif

#if defined(FUNC_gcd_ui) || defined(FUNC_lcm_ui)
static UV gcd_ui(UV x, UV y) {
  UV t;
  if (y < x) { t = x; x = y; y = t; }
  while (y > 0) {
    t = y;  y = x % y;  x = t;  /* y1 <- x0 % y0 ; x1 <- y0 */
  }
  return x;
}
#endif

#ifdef FUNC_lcm_ui
static UV lcm_ui(UV x, UV y) {
  /* Can overflow if lcm(x,y) > 2^64 (e.g. two primes each > 2^32) */
  return x * (y / gcd_ui(x,y));
}
#endif

#endif
