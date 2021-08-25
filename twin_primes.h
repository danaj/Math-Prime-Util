#ifndef MPU_TWIN_PRIMES_H
#define MPU_TWIN_PRIMES_H

#include "ptypes.h"

extern UV  twin_prime_count(UV low, UV high);
extern UV  twin_prime_count_approx(UV n);
extern UV  nth_twin_prime(UV n);
extern UV  nth_twin_prime_approx(UV n);

#endif
