#ifndef MPU_ALMOST_PRIMES_H
#define MPU_ALMOST_PRIMES_H

#include "ptypes.h"

#if 0
extern UV  nth_almost_prime(UV n, UV k);
extern UV  nth_almost_prime_approx(UV n, UV k);
extern UV  almost_prime_count(UV low, UV high, UV k);
extern UV  almost_prime_count_approx(UV n, UV k);
extern UV range_almost_prime_sieve(UV** semis, UV lo, UV hi, UV k);
#endif

extern UV  almost_prime_count(UV n, UV k);
extern UV  almost_prime_count_approx(UV n, UV k);

#endif
