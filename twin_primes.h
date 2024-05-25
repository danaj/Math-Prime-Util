#ifndef MPU_TWIN_PRIMES_H
#define MPU_TWIN_PRIMES_H

#include "ptypes.h"

extern UV  twin_prime_count_range(UV beg, UV end);

extern UV  twin_prime_count(UV n);
extern UV  twin_prime_count_approx(UV n);

extern UV  nth_twin_prime(UV n);
extern UV  nth_twin_prime_approx(UV n);

extern UV  range_twin_prime_sieve(UV** list, UV lo, UV hi);

#endif
