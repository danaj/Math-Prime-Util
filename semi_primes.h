#ifndef MPU_SEMI_PRIMES_H
#define MPU_SEMI_PRIMES_H

#include "ptypes.h"

extern UV  range_semiprime_sieve(UV** semis, UV lo, UV hi);

extern UV  semiprime_count_range(UV lo, UV hi);

extern UV  semiprime_count(UV n);
extern UV  semiprime_count_approx(UV n);

extern UV  nth_semiprime(UV n);
extern UV  nth_semiprime_approx(UV n);

#endif
