#ifndef MPU_PRIME_SUMS_H
#define MPU_PRIME_SUMS_H

#include "ptypes.h"

extern bool sum_primes(UV low, UV high, UV *sum);

#define HAVE_SUM_PRIMES128 (HAVE_UINT64 && HAVE_UINT128)
#if HAVE_SUM_PRIMES128
extern bool sum_primes128(uint64_t lo, uint64_t hi, uint128_t *sum);
#endif

#endif
