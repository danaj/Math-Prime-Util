#ifndef MPU_RAMANUJAN_PRIMES_H
#define MPU_RAMANUJAN_PRIMES_H

#include "ptypes.h"

extern UV* n_ramanujan_primes(UV n);
extern UV* n_range_ramanujan_primes(UV nlo, UV nhi);
extern UV* ramanujan_primes(UV* first, UV* last, UV low, UV high);
extern int is_ramanujan_prime(UV n);
extern UV  ramanujan_prime_count(UV lo, UV hi);
extern UV  ramanujan_prime_count_upper(UV n);
extern UV  ramanujan_prime_count_lower(UV n);
extern UV  ramanujan_prime_count_approx(UV n);
extern UV  nth_ramanujan_prime(UV n);
extern UV  nth_ramanujan_prime_upper(UV n);
extern UV  nth_ramanujan_prime_lower(UV n);
extern UV  nth_ramanujan_prime_approx(UV n);

#endif
