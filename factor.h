#ifndef MPU_FACTOR_H
#define MPU_FACTOR_H

#include "ptypes.h"

#define MPU_MAX_FACTORS 64

extern int factor(UV n, UV *factors);

extern int trial_factor(UV n, UV *factors, UV maxtrial);

extern int fermat_factor(UV n, UV *factors, UV rounds);
extern int holf_factor(UV n, UV *factors, UV rounds);
extern int pbrent_factor(UV n, UV *factors, UV maxrounds, UV a);
extern int prho_factor(UV n, UV *factors, UV maxrounds);
extern int pminus1_factor(UV n, UV *factors, UV B1, UV B2);
extern int squfof_factor(UV n, UV *factors, UV rounds);
extern int racing_squfof_factor(UV n, UV *factors, UV rounds);


extern int _XS_is_pseudoprime(UV n, UV a);
extern int _XS_miller_rabin(UV n, const UV *bases, int nbases);
extern int _SPRP2(UV n);
extern int _XS_is_prob_prime(UV n);
extern int _XS_is_strong_lucas_pseudoprime(UV n);
extern int _XS_is_frobenius_underwood_pseudoprime(UV n);

extern UV _XS_divisor_sum(UV n);

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
