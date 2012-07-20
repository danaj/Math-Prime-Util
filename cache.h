#ifndef MPU_CACHE_H
#define MPU_CACHE_H

#include "EXTERN.h"
#include "perl.h"

  /* Sieve from 0 to x and store in primary cache */
extern void  prime_precalc(UV x);
  /* Release all extra memory -- go back to initial amounts */
extern void  prime_memfree(void);
  /* Seriously shut everything down, including destroying mutexes.
   * This should ONLY be called when we're leaving for good.
   */
extern void  _prime_memfreeall(void);

  /* Get the primary cache (mod-30 wheel sieve).
   * Try to make sure it contains n.
   * Returns the maximum value in the cache.
   * Sets sieve* to the cache, unless given 0.
   * If you get a pointer back, you MUST call release when you're done.
   *
   * Ex: just give me the current size:
   *   UV cache_size = get_prime_cache(0, 0);
   *
   * Ex: give me the current cache and size:
   *   UV cache_size = get_prime_cache(0, &sieve);
   *
   * Ex: give me the cache at least size n:
   *   UV cache_size = get_prime_cache(n, &sieve);
   */
extern UV   get_prime_cache(UV n, const unsigned char** sieve);
  /* Inform the system we're done using the primary cache if we got a ptr. */
extern void release_prime_cache(const unsigned char* sieve);

  /* Get the segment cache.  Set size to its size. */
extern unsigned char* get_prime_segment(UV* size);
  /* Inform the system we're done using the segment cache. */
extern void release_prime_segment(unsigned char* segment);

#endif
