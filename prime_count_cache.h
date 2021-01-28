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

#endif
