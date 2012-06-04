#ifndef MPU_FACTOR_H
#define MPU_FACTOR_H

#include "EXTERN.h"
#include "perl.h"

#define MPU_MAX_FACTORS 64

extern int trial_factor(UV n, UV *factors, UV maxtrial);

extern int fermat_factor(UV n, UV *factors);

extern int squfof_factor(UV n, UV *factors, UV rounds);

extern int pbrent_factor(UV n, UV *factors, UV maxrounds);

extern int prho_factor(UV n, UV *factors, UV maxrounds);

extern int pminus1_factor(UV n, UV *factors, UV maxrounds);

#endif
