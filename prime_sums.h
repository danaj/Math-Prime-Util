#ifndef MPU_PRIME_SUMS_H
#define MPU_PRIME_SUMS_H

#include "ptypes.h"

extern int sum_primes(UV low, UV high, UV *sum);
extern UV  sum_primes64(UV n);
#define HAVE_SUM_PRIMES128 (BITS_PER_WORD == 64 && HAVE_UINT128)
extern int sum_primes128(UV n, UV *hisum, UV *losum);

extern int sum_primes128_lmos(UV n, UV *hi_sum, UV *lo_sum);
#endif
