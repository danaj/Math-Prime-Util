#ifndef MPU_PRIME_SUMS_H
#define MPU_PRIME_SUMS_H

#include "ptypes.h"

extern bool sum_primes(UV low, UV high, UV *sum);
extern UV  sum_primes64(UV n);
#define HAVE_SUM_PRIMES128 (BITS_PER_WORD == 64 && HAVE_UINT128)
extern bool sum_primes128(UV n, UV *hisum, UV *losum);

#endif
