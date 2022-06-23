#ifndef MPU_OMEGA_PRIMES_H
#define MPU_OMEGA_PRIMES_H

#include "ptypes.h"

/* Almost primes use prime_bigomega, Omega primes use prime_omega */

extern int is_omega_prime(uint32_t k, UV n);

extern UV max_nth_omega_prime(uint32_t k);
extern UV max_omega_prime_count(uint32_t k);

extern UV range_omega_prime_sieve(UV** ret, uint32_t k, UV lo, UV hi);

extern UV omega_prime_count(uint32_t k, UV n);
extern UV nth_omega_prime(uint32_t k, UV n);

#endif
