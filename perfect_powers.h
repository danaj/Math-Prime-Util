#ifndef MPU_PERFECT_POWERS_H
#define MPU_PERFECT_POWERS_H

#include "ptypes.h"

extern int is_perfect_power(UV n);

extern UV  next_perfect_power(UV x);
extern UV  prev_perfect_power(UV x);

extern UV  perfect_power_count_range(UV lo, UV hi);

extern UV  perfect_power_count(UV n);
extern UV  perfect_power_count_lower(UV n);
extern UV  perfect_power_count_upper(UV n);
extern UV  perfect_power_count_approx(UV n);

extern UV  nth_perfect_power(UV n);
extern UV  nth_perfect_power_lower(UV n);
extern UV  nth_perfect_power_upper(UV n);
extern UV  nth_perfect_power_approx(UV n);

#endif
