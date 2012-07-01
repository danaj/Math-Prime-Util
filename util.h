#ifndef MPU_UTIL_H
#define MPU_UTIL_H

#include "ptypes.h"

extern int is_prime(UV x);
extern int is_definitely_prime(UV x);
extern UV  next_trial_prime(UV x);
extern UV  next_prime(UV x);
extern UV  prev_prime(UV x);

extern UV  prime_count(UV low, UV high);
extern UV  nth_prime(UV x);

/* These have been moved into the main Util.pm */
extern UV  _XS_prime_count_lower(UV x);
extern UV  _XS_prime_count_upper(UV x);
extern UV  _XS_prime_count_approx(UV x);
extern UV  _XS_nth_prime_lower(UV n);
extern UV  _XS_nth_prime_upper(UV x);
extern UV  _XS_nth_prime_approx(UV x);

extern double _XS_ExponentialIntegral(double x);
extern double _XS_LogarithmicIntegral(double x);
extern double _XS_RiemannR(double x);

/* Above this value, is_prime will do deterministic Miller-Rabin */
/* With 64-bit math, we can do much faster mulmods from 2^16-2^32 */
#if (BITS_PER_WORD == 64) || HAVE_STD_U64
  #define MPU_PROB_PRIME_BEST  UVCONST(100000)
#else
  #define MPU_PROB_PRIME_BEST  UVCONST(100000000)
#endif

#endif
