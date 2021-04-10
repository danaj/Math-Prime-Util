#ifndef MPU_PRIME_COUNT_CACHE_H
#define MPU_PRIME_COUNT_CACHE_H

#include "ptypes.h"

/*
 * Create will attempt to cache up to n, but possibly less.
 *
 * Destroy is necessary to call to free the memory.
 *
 * Lookup will return prime_count(n) (it looks in cache first, then calls LMO).
 */

extern void* prime_count_cache_create(UV n);
extern void  prime_count_cache_destroy(void* cache);
extern UV    prime_count_cache_lookup(void* cache, UV n);

/* Alternate creation method using an array of primes:
 *     P[0]=0, P[1]=2, P[2]=3, P[3]=5, ..., P[lastidx] = <nth_prime(lastidx)>
 * If the caller already has this, this can be faster.
 */
extern void* prime_count_cache_create_with_primes(const uint32_t *primes, uint32_t lastidx);

#endif
