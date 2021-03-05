#ifndef MPU_ALMOST_PRIMES_H
#define MPU_ALMOST_PRIMES_H

#include "ptypes.h"

extern UV  almost_prime_count(uint32_t k, UV n);
extern UV  almost_prime_count_approx(uint32_t k, UV n);

extern UV  nth_almost_prime(uint32_t k, UV n);
extern UV  nth_almost_prime_approx(uint32_t k, UV n);
extern UV  nth_almost_prime_lower(uint32_t k, UV n);
extern UV  nth_almost_prime_upper(uint32_t k, UV n);

/* The largest k-almost-prime that fits in a UV */
extern UV max_nth_almost_prime(uint32_t k);
/* The k-almost-prime count for 2^64-1 */
extern UV max_almost_prime_count(uint32_t k);

/* Reasonably tight bounds on the counts */
extern UV almost_prime_count_upper(uint32_t k, UV n);
extern UV almost_prime_count_lower(uint32_t k, UV n);

extern UV range_construct_almost_prime(UV** list, uint32_t k, UV lo, UV hi);
extern UV range_almost_prime_sieve(UV** list, uint32_t k, UV lo, UV hi);

/* Using prime_omega rather than prime_bigomega */
extern UV range_omega_prime_sieve(UV** ret, uint32_t k, UV lo, UV hi);
extern UV omega_prime_count(uint32_t k, UV n);

#endif
