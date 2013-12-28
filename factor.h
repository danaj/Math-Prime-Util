#ifndef MPU_FACTOR_H
#define MPU_FACTOR_H

#include "ptypes.h"

#define MPU_MAX_FACTORS 64

extern int factor(UV n, UV *factors);
extern int factor_exp(UV n, UV *factors, UV* exponents);

extern int trial_factor(UV n, UV *factors, UV maxtrial);

extern int fermat_factor(UV n, UV *factors, UV rounds);
extern int holf_factor(UV n, UV *factors, UV rounds);
extern int pbrent_factor(UV n, UV *factors, UV maxrounds, UV a);
extern int prho_factor(UV n, UV *factors, UV maxrounds);
extern int pminus1_factor(UV n, UV *factors, UV B1, UV B2);
extern int pplus1_factor(UV n, UV *factors, UV B);
extern int squfof_factor(UV n, UV *factors, UV rounds);
extern int racing_squfof_factor(UV n, UV *factors, UV rounds);

extern UV  _XS_divisor_sum(UV n, UV k);
extern UV* _divisor_list(UV n, UV *num_divisors);

#endif
