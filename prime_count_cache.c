#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FUNC_popcnt 1
#include "ptypes.h"
#include "prime_nth_count.h"
#include "sieve.h"
#include "cache.h"
#include "lmo.h"
#include "util.h"
#include "prime_count_cache.h"

/*
 * Cache small counts directly using a uint32_t array.
 * Very fast, but space intensive.
 *
 * Cache larger counts using a base count + a single-word bit count.
 *
 * We used to use a binary search on a prime list, which is reasonable,
 * but the bit mask uses less memory and is faster on average.  It also
 * easily allows larger sizes.  Note: in 32-bit this isn't very efficient.
 *
 * If memory is a concern, we could switch to a base count every two words.
 */

typedef struct {
  uint16_t *count;
  uint32_t *bm_count;
  UV       *bm_mask;

  UV last_n;
  UV last_count_n;
  UV last_bmask_n;

  /* Statistics counting */
  unsigned long nl_small;
  unsigned long nl_bmask;
  unsigned long nl_lmo;
} pc_cache_t;


UV prime_count_cache_lookup(void* cobj, UV n) {
  pc_cache_t *cache = (pc_cache_t*)cobj;

  if (n <= 2)
    return (n==2);
  /* Look in the small direct cache. */
  if (n <= cache->last_count_n) {
    cache->nl_small++;
    return cache->count[(n-1)>>1];
  }
  /* Look in bitmask */
  if (n <= cache->last_bmask_n) {
    UV m = (n-1) >> 1;
    uint32_t idx = m / BITS_PER_WORD;
    uint32_t rem = m % BITS_PER_WORD;
    cache->nl_bmask++;
    return (UV)cache->bm_count[idx]  +  popcnt(cache->bm_mask[idx] >> (BITS_PER_WORD - 1 - rem));
  }
  /* OK, call LMO/segment */
  cache->nl_lmo++;
  return LMO_prime_count(n);
}


#if 0
static void _checkn(pc_cache_t *cache, UV n, UV count) {
  UV pc = prime_count_cache_lookup(cache, n);
  if (pc != count)
    croak("  pc cache [%lu] returned %lu instead of %lu\n", n, pc, count);
}
static void verify_cache(pc_cache_t *cache) {
  UV n = 3, c = 1, lastn = cache->last_n;

  _checkn(cache, 0, 0);
  _checkn(cache, 1, 0);
  _checkn(cache, 2, 1);
  START_DO_FOR_EACH_PRIME(3, next_prime(lastn)) {
    while (n < p)  _checkn(cache, n++, c);
    _checkn(cache, n++, ++c);
  } END_DO_FOR_EACH_PRIME
  printf("  prime count cache verified to %lu complete\n", lastn);
}
static UV _bm_lookup(pc_cache_t *cache, UV n) {
  uint32_t m = (n-1) >> 1;
  uint32_t idx = m / BITS_PER_WORD;
  uint32_t rem = m % BITS_PER_WORD;
  return cache->bm_count[idx]  +  popcnt(cache->bm_mask[idx] >> (BITS_PER_WORD - 1 - rem));
}
#else
#define verify_cache(cache)  /* nothing */
#define _bm_lookup(cache,n) prime_count_cache_lookup(cache,n)
#endif


void prime_count_cache_destroy(void* cobj) {
  pc_cache_t *cache = (pc_cache_t*)cobj;

  MPUverbose(2, "  Prime Count Cache (max %lu):\n", (UV)cache->last_n);
  MPUverbose(2, "    Small: %lu (%luk)   Mask: %lu (%luk)\n",
    (unsigned long)cache->last_count_n,
    cache->last_count_n ? (unsigned long)(((cache->last_count_n-1)>>1)+1)*4/1024 : 0,
    (unsigned long)cache->last_bmask_n,
    (unsigned long) (sizeof(UV)+sizeof(uint32_t)) * (cache->last_bmask_n/(2*BITS_PER_WORD) + 1) / 1024);
  MPUverbose(2, "    Lookups  Small %lu  Mask %lu  LMO %lu\n",
    cache->nl_small, cache->nl_bmask, cache->nl_lmo);

  if (cache->count != 0)
    Safefree(cache->count);
  if (cache->bm_count != 0)
    Safefree(cache->bm_count);
  if (cache->bm_mask != 0)
    Safefree(cache->bm_mask);
  Safefree(cache);
}


/* prime_count(LIM_SMALL) <= 65535 */
#define LIM_SMALL 821640

void* prime_count_cache_create(UV n) {
  pc_cache_t *cache;
  uint32_t i, idx, cnt;

  if (n < 5)  n = 5;
#if BITS_PER_WORD == 64
  /* The prime count has to fit in a uint32_t, so must be < 104484802057 */
  /* Further limit to ~ 3GB. */
  if (n > UVCONST( 34359738367))   n = UVCONST( 34359738367);
#endif
  prime_precalc(LIM_SMALL);

  Newz(0, cache, 1, pc_cache_t);  /* Allocate cache object, everything zero */
  cache->last_n = n;

  /* Fill in small counts */
  {
    uint16_t *counts;
    uint32_t count_last_n   = (n <= LIM_SMALL) ? n : LIM_SMALL;
    uint32_t count_last_idx = (count_last_n-1) >> 1;
    New(0, counts, count_last_idx+1, uint16_t);
    counts[0] = 1;
    idx = cnt = 1;
    START_DO_FOR_EACH_PRIME(3, count_last_n) {
      while (idx < ((p-1)>>1))  counts[idx++] = cnt;
      counts[idx++] = ++cnt;
    } END_DO_FOR_EACH_PRIME
    MPUassert(cnt <= 65535, "small count overflow");
    while (idx <= count_last_idx)  counts[idx++] = cnt;

    cache->count = counts;
    cache->last_count_n = count_last_n;
  }

  /* Fill in bitmask and base counts */
  if (n > cache->last_count_n) {
    UV       *mask;
    uint32_t *count;
    uint32_t words = (n / (2*BITS_PER_WORD)) + 1;  /* 0-127=1, 128-255=2 */
    Newz(0, count, words, uint32_t);
    Newz(0,  mask, words, UV);

    mask[0] = UVCONST(15) << (BITS_PER_WORD-4);
    {
      unsigned char* segment;
      UV seg_base, seg_low, seg_high;
      void* ctx = start_segment_primes(7, n, &segment);
      while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
        START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
          UV m = (p-1)>>1;
          uint32_t midx = m / BITS_PER_WORD;
          uint32_t mrem = m % BITS_PER_WORD;
          mask[midx] |= UVCONST(1) << (BITS_PER_WORD-1-mrem);
        END_DO_FOR_EACH_SIEVE_PRIME
      }
      end_segment_primes(ctx);
    }
    for (i = 1; i < words; i++)
      count[i] = count[i-1] + popcnt(mask[i-1]);

    cache->bm_mask  = mask;
    cache->bm_count = count;
    cache->last_bmask_n = n;
  }
  verify_cache(cache);
  return cache;
}

void* prime_count_cache_create_with_primes(const uint32_t *primes, uint32_t lastidx) {
#if 0  /* Slower  */
  return prime_count_cache_create(primes[lastidx]);
#else  /* Faster, but so much code duplication.... */
  pc_cache_t *cache;
  uint32_t i, idx, cnt, n;

  MPUassert(primes != 0, "prime_count_cache_create called with null pointer");
  if (lastidx <= 1) return prime_count_cache_create(5);
  if (lastidx > 203280221) lastidx = 203280221;

  Newz(0, cache, 1, pc_cache_t);  /* Allocate cache object, everything zero */
  cache->last_n = n = primes[lastidx];

  /* Fill in small counts */
  {
    uint16_t *counts;
    uint32_t count_last_n   = (n <= LIM_SMALL) ? n : LIM_SMALL;
    uint32_t count_last_idx = (count_last_n-1) >> 1;
    New(0, counts, count_last_idx+1, uint16_t);
    counts[0] = 1;
    idx = cnt = 1;
    for (i = 2;  i <= lastidx;  i++) {
      uint32_t p = primes[i];
      if (p > count_last_n) break;
      while (idx < ((p-1)>>1))  counts[idx++] = cnt;
      counts[idx++] = ++cnt;
    }
    MPUassert(cnt <= 65535, "small count overflow");
    while (idx <= count_last_idx)  counts[idx++] = cnt;

    cache->count = counts;
    cache->last_count_n = count_last_n;
  }

  /* Fill in bitmask and base counts */
  if (n > cache->last_count_n) {
    UV       *mask;
    uint32_t *count;
    uint32_t words = (n / (2*BITS_PER_WORD)) + 1;  /* 0-127=1, 128-255=2 */
    Newz(0, count, words, uint32_t);
    Newz(0,  mask, words, UV);

    mask[0] = UVCONST(1) << (BITS_PER_WORD-1);
    for (i = 2;  i <= lastidx;  i++) {
      uint32_t p = primes[i];
      uint32_t m = (p-1)>>1;
      uint32_t midx = m / BITS_PER_WORD;
      uint32_t mrem = m % BITS_PER_WORD;
      mask[midx] |= UVCONST(1) << (BITS_PER_WORD-1-mrem);
    }

    for (i = 1; i < words; i++)
      count[i] = count[i-1] + popcnt(mask[i-1]);

    cache->bm_mask  = mask;
    cache->bm_count = count;
    cache->last_bmask_n = n;
  }
  verify_cache(cache);
  return cache;
#endif
}
