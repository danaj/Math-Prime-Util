#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "cache.h"
#include "sieve.h"
#include "constants.h"   /* _MPU_FILL_EXTRA_N and _MPU_INITIAL_CACHE_SIZE */

#include "threadlock.h"

/*
 * These functions are used internally by the .c and .xs files.
 * They handle a cached primary set of primes, as well as a segment
 * area for use by all the functions that want to do segmented operation.
 *
 * We must be thread-safe, and we want to allow a good deal of concurrency.
 * It is imperative these be used correctly.  After calling the get method,
 * use the sieve or segment, then release.  You MUST call release before you
 * return or croak.  You ought to release as soon as you're done using the
 * sieve or segment.
 */

static int mutex_init = 0;
MUTEX_DECL(segment);
READ_WRITE_LOCK_DECL(primary_cache);

static unsigned char* prime_cache_sieve = 0;
static UV             prime_cache_size = 0;

/* Erase the primary cache and fill up to n. */
/* Note: You must have a write lock before calling this! */
static void _erase_and_fill_prime_cache(UV n) {
  UV padded_n;

  if (n >= (UV_MAX-_MPU_FILL_EXTRA_N))
    padded_n = UV_MAX;
  else
    padded_n = ((n + _MPU_FILL_EXTRA_N)/30)*30;

  /* If new size isn't larger or smaller, then we're done. */
  if (prime_cache_size == padded_n)
    return;

  if (prime_cache_sieve != 0)
    Safefree(prime_cache_sieve);
  prime_cache_sieve = 0;
  prime_cache_size = 0;

  if (n > 0) {
    prime_cache_sieve = sieve_erat30(padded_n);
    MPUassert(prime_cache_sieve != 0, "sieve returned null");
    prime_cache_size = padded_n;
  }
}

/*
 * Get the size and a pointer to the cached prime sieve.
 * Returns the maximum sieved value available.
 * Allocates and sieves if needed.
 *
 * The sieve holds 30 numbers per byte, using a mod-30 wheel.
 */
UV get_prime_cache(UV n, const unsigned char** sieve)
{
#ifdef USE_ITHREADS
  if (sieve == 0) {
    if (prime_cache_size < n) {
      WRITE_LOCK_START(primary_cache);
        _erase_and_fill_prime_cache(n);
      WRITE_LOCK_END(primary_cache);
    }
    return prime_cache_size;
  }

  /* This could be done more efficiently if we converted a write lock to a
   * reader after doing the expansion.  But I think this solution is less
   * error prone (though could lead to starvation in pathological cases).
   */
  READ_LOCK_START(primary_cache);
  while (prime_cache_size < n) {
    /* The cache isn't big enough.  Expand it. */
    READ_LOCK_END(primary_cache);
    /* thread reminder: the world can change right here */
    WRITE_LOCK_START(primary_cache);
      if (prime_cache_size < n)
        _erase_and_fill_prime_cache(n);
    WRITE_LOCK_END(primary_cache);
    /* thread reminder: the world can change right here */
    READ_LOCK_START(primary_cache);
  }
  MPUassert(prime_cache_size >= n, "prime cache is too small!");

  *sieve = prime_cache_sieve;
  return prime_cache_size;
#else
  if (prime_cache_size < n)
    _erase_and_fill_prime_cache(n);
  MPUassert(prime_cache_size >= n, "prime cache is too small!");
  if (sieve != 0)
    *sieve = prime_cache_sieve;
  return prime_cache_size;
#endif
}

#ifdef USE_ITHREADS
void release_prime_cache(const unsigned char* mem) {
  (void)mem; /* We don't currently care about the pointer */
  READ_LOCK_END(primary_cache);
}
#endif



/* The segment everyone is trying to share */
#define PRIMARY_SEGMENT_CHUNK_SIZE    UVCONST(256*1024-16)
static unsigned char* prime_segment = 0;
static int prime_segment_is_available = 1;
/* If that's in use, malloc a new one of this size */
#define SECONDARY_SEGMENT_CHUNK_SIZE  UVCONST(128*1024-16)

unsigned char* get_prime_segment(UV *size) {
  unsigned char* mem;
  int use_prime_segment = 0;

  MPUassert(size != 0, "get_prime_segment given null size pointer");
  MPUassert(mutex_init == 1, "segment mutex has not been initialized");

  MUTEX_LOCK(&segment_mutex);
    if (prime_segment_is_available) {
      prime_segment_is_available = 0;
      use_prime_segment = 1;
    }
  MUTEX_UNLOCK(&segment_mutex);

  if (use_prime_segment) {
    if (prime_segment == 0)
      New(0, prime_segment, PRIMARY_SEGMENT_CHUNK_SIZE, unsigned char);
    *size = PRIMARY_SEGMENT_CHUNK_SIZE;
    mem = prime_segment;
  } else {
    New(0, mem, SECONDARY_SEGMENT_CHUNK_SIZE, unsigned char);
    *size = SECONDARY_SEGMENT_CHUNK_SIZE;
  }
  MPUassert(mem != 0, "get_prime_segment allocation failure");

  return mem;
}

void release_prime_segment(unsigned char* mem) {
  MUTEX_LOCK(&segment_mutex);
    if (mem == prime_segment) {
      prime_segment_is_available = 1;
      mem = 0;
    }
  MUTEX_UNLOCK(&segment_mutex);
  if (mem)
    Safefree(mem);
}



void prime_precalc(UV n)
{
  if (!mutex_init) {
    MUTEX_INIT(&segment_mutex);
    MUTEX_INIT(&primary_cache_mutex);
    COND_INIT(&primary_cache_turn);
    mutex_init = 1;
  }

  /* On initialization, make a few primes (30k per 1k memory) */
  if (n == 0)
    n = _MPU_INITIAL_CACHE_SIZE;
  get_prime_cache(n, 0);   /* Sieve to n */

  /* TODO: should we prealloc the segment here? */
}


void prime_memfree(void)
{
  unsigned char* old_segment = 0;

  /* This can happen in global destructor, and PL_dirty has porting issues */
  /* MPUassert(mutex_init == 1, "cache mutexes have not been initialized"); */
  if (mutex_init == 0) return;

  MUTEX_LOCK(&segment_mutex);
  /* Don't free if another thread is using it */
  if ( (prime_segment != 0) && (prime_segment_is_available) ) {\
    unsigned char* new_segment = old_segment;
    old_segment = prime_segment;
    prime_segment = new_segment; /* Exchanged old_segment / prime_segment */
  }
  MUTEX_UNLOCK(&segment_mutex);
  if (old_segment) Safefree(old_segment);

  WRITE_LOCK_START(primary_cache);
    /* Put primary cache back to initial state */
    _erase_and_fill_prime_cache(_MPU_INITIAL_CACHE_SIZE);
  WRITE_LOCK_END(primary_cache);
}


void _prime_memfreeall(void)
{
  /* No locks.  We're shutting everything down. */
  if (mutex_init) {
    mutex_init = 0;
    MUTEX_DESTROY(&segment_mutex);
    MUTEX_DESTROY(&primary_cache_mutex);
    COND_DESTROY(&primary_cache_turn);
  }
  if (prime_cache_sieve != 0)
    Safefree(prime_cache_sieve);
  prime_cache_sieve = 0;
  prime_cache_size = 0;

  if (prime_segment != 0)
    Safefree(prime_segment);
  prime_segment = 0;
}
