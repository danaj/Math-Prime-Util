#ifndef MPU_PRIME_NTH_COUNT_H
#define MPU_PRIME_NTH_COUNT_H

#include "ptypes.h"

extern UV  _XS_prime_count(UV low, UV high);
extern UV  nth_prime(UV x);
extern UV  nth_prime_upper(UV x);
extern UV  nth_prime_lower(UV x);
extern UV  nth_prime_approx(UV x);
extern UV  prime_count_upper(UV x);
extern UV  prime_count_lower(UV x);
extern UV  prime_count_approx(UV x);
extern UV  twin_prime_count(UV low, UV high);
extern UV  twin_prime_count_approx(UV n);
extern UV  nth_twin_prime(UV n);
extern UV  nth_twin_prime_approx(UV n);
extern int sum_primes(UV low, UV high, UV *sum);
extern int sum_primes128(UV n, UV *hisum, UV *losum);

#endif
