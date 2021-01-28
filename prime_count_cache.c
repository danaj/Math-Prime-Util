#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "prime_nth_count.h"
#include "sieve.h"
#include "cache.h"
#include "lmo.h"

/*
 * Cache small counts directly using a uint32_t array.
 * Very fast, but space intensive for larger sizes.
 * We limit this to 16MB, e.g. counts to 8M.
 *
 * If caching past this is requested, also make a uint32_t list of primes.
 * We use this to do a binary search.
 *
 */

typedef struct {
  uint32_t *count;
  uint32_t *primes;
  uint32_t last_count_n;
  uint32_t last_prime_idx;
  /* Statistics counting */
  unsigned long nl_total;
  unsigned long nl_small;
  unsigned long nl_large;
  unsigned long nl_lmo;
} pc_cache_t;

UV prime_count_cache_lookup(void* cobj, UV n) {
  pc_cache_t *cache = (pc_cache_t*)cobj;
  uint32_t *primes;
  UV i, j;

  cache->nl_total++;

  if (n <= 2)
    return (n==2);

  /* Look in the small direct cache. */
  if (n <= cache->last_count_n) {
    cache->nl_small++;
    return cache->count[(n-1)>>1];
  }

  /* Binary search if within our limits, otherwise call LMO to calculate. */
  primes = cache->primes;
  if (!primes || n >= primes[cache->last_prime_idx]) {
    cache->nl_lmo++;
    return LMO_prime_count(n);
  }

  cache->nl_large++;

  /* Select starting and ending limits for the binary search */
  j = cache->last_prime_idx;
#if 0
  if      (n <     8480) { i = 1 + (n>>4);    if (j > 1060) j = 1060; }
  else if (n < 25875000) { i = 793 + (n>>5);  if (j > (n>>3)) j = n>>3; }
  else                   { i = 1617183;       if (j > (n>>4)) j = n>>4; }
#else
  if (n < 8480) {
    i =      1 + ((n*31)>>8);
    j =      3 + ((n*68)>>8) + 1;
    if (j > 1060) j = 1060;
  } else if (n < 8400000) {
    i =    495 + ((n*17)>>8);
    j =   4620 + ((n*19)>>8) + 1;
    if (j > (n>>3)) j = n>>3;
  } else if (n < 25875000) {
    i =  72670 + ((n*15)>>8);
    j = 101170 + ((n*15)>>8) + 1;
    if (j > (n>>3)) j = n>>3;
  } else {
#if BITS_PER_WORD == 64
    i = 404290 + ((n*12)>>8);
    j = 890170 + ((n*14)>>8) + 1;
#else
    i = 1617183;
#endif
    if (j > (n>>4)) j = n>>4;
  }
  if (j > cache->last_prime_idx) j = cache->last_prime_idx;
#endif

  while (i < j) {
    UV mid = i + (j-i)/2;  
    if (primes[mid] <= n)  i = mid+1;
    else                   j = mid;
  }
  return i-1;
}

void prime_count_cache_destroy(void* cobj) {
  pc_cache_t *cache = (pc_cache_t*)cobj;

  MPUverbose(2, "  Prime Count Cache:  Small: %u (%uk)  Large: %u (%uk)\n  Lookups  Total %lu  Small %lu  Large %lu  LMO %lu\n",
    cache->last_count_n, (((cache->last_count_n-1)>>1)+1)*4/1024,
    cache->primes[cache->last_prime_idx], (cache->last_prime_idx+1)*4/1024,
    cache->nl_total, cache->nl_small, cache->nl_large, cache->nl_lmo);

  if (cache->primes != 0)
    Safefree(cache->primes);
  if (cache->count != 0)
    Safefree(cache->count);
  Safefree(cache);
}

#if 1 /* About 16MB */
#define LIM_SMALL  8388607
#else /* About 8MB */
#define LIM_SMALL  4194303
#endif

void* prime_count_cache_create(UV n) {
  pc_cache_t *cache;
  uint32_t *counts = 0,  *primes = 0;
  uint32_t count_last_n = 0,  prime_last_idx = 0;

  if (n < 4)  n = 4;
  if (n > 4294967295U)   n = 4294967295U;   /* Limit to uint32_t */
  prime_precalc(LIM_SMALL);

  { /* First make the counts */
    uint32_t idx = 1, cnt = 1, count_last_idx;

    count_last_n   = ((n <= LIM_SMALL)  ?  n  :  LIM_SMALL) | 1;
    count_last_idx = (count_last_n-1) >> 1;
    New(0, counts, count_last_idx+1, uint32_t);
    counts[0] = 1;
    START_DO_FOR_EACH_PRIME(3, count_last_n) {
      while (idx < ((p-1)>>1)) counts[idx++] = cnt;
      counts[idx++] = ++cnt;
    } END_DO_FOR_EACH_PRIME
    while (idx <= count_last_idx) counts[idx++] = cnt;
  }
  /* If more needed, generate small primes. */
  if (n > count_last_n)
    prime_last_idx = range_prime_sieve_32(&primes, n, /* Offset */ 1);

  New(0, cache, 1, pc_cache_t);
  cache->last_count_n = count_last_n;
  cache->count = counts;
  cache->last_prime_idx = prime_last_idx;
  cache->primes = primes;
  cache->nl_total = cache->nl_small = cache->nl_large = cache->nl_lmo = 0;
  return cache;
}
