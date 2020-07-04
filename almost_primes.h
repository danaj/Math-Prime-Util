#ifndef MPU_ALMOST_PRIMES_H
#define MPU_ALMOST_PRIMES_H

#include "ptypes.h"

#if 0
extern UV  nth_almost_prime_approx(UV n, UV k);
extern UV  almost_prime_count(UV low, UV high, UV k);
extern UV  almost_prime_count_approx(UV n, UV k);
extern UV range_almost_prime_sieve(UV** semis, UV lo, UV hi, UV k);
#endif

extern UV  almost_prime_count(UV k, UV n);
extern UV  almost_prime_count_approx(UV k, UV n);

extern UV  nth_almost_prime(UV k, UV n);
extern UV  nth_almost_prime_approx(UV k, UV n);
extern UV  nth_almost_prime_lower(UV k, UV n);
extern UV  nth_almost_prime_upper(UV k, UV n);

/* The largest k-almost-prime that fits in a UV */
extern UV max_nth_almost_prime(UV k);
/* The k-almost-prime count for 2^64-1 */
extern UV max_almost_prime_count(UV k);

/* Reasonably tight bounds on the counts */
extern UV almost_prime_count_upper(UV k, UV n);
extern UV almost_prime_count_lower(UV k, UV n);

extern UV range_construct_almost_prime(UV** list, UV k, UV lo, UV hi);

#endif
