#ifndef MPU_PRIME_COUNTS_H
#define MPU_PRIME_COUNTS_H

#include "ptypes.h"


extern UV  segment_prime_count(UV low, UV high);

extern UV  prime_count_range(UV low, UV high);

extern UV  prime_count(UV x);
extern UV  prime_count_upper(UV x);
extern UV  prime_count_lower(UV x);
extern UV  prime_count_approx(UV x);

extern UV  nth_prime(UV x);
extern UV  nth_prime_upper(UV x);
extern UV  nth_prime_lower(UV x);
extern UV  nth_prime_approx(UV x);

/* Used for Ramanujan prime bounds */
extern double ramanujan_axler(long double n, long double c, long double d);

#endif
