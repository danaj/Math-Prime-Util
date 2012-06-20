#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "cache.h"
#include "sieve.h"

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

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
#ifdef USE_ITHREADS
static perl_mutex segment_mutex;

/* See: http://en.wikipedia.org/wiki/Readers-writers_problem */
static perl_mutex primary_mutex_no_waiting;
static perl_mutex primary_mutex_no_accessing;
static perl_mutex primary_mutex_counter;
static int        primary_number_of_readers = 0;
#endif

static unsigned char* prime_cache_sieve = 0;
static UV             prime_cache_size = 0;

/* Fill the primary cache up to n */
static void _fill_prime_cache(UV n, int nowait) {
  if (!nowait) { MUTEX_LOCK(&primary_mutex_no_waiting); }
  MUTEX_LOCK(&primary_mutex_no_accessing);
  if (!nowait) { MUTEX_UNLOCK(&primary_mutex_no_waiting); }
  if (prime_cache_size < n) {
    if (prime_cache_sieve != 0)
      Safefree(prime_cache_sieve);
    prime_cache_sieve = 0;
    prime_cache_size = 0;

    /* Sieve a bit more than asked, to mitigate thrashing */
    if (n >= (UV_MAX-3840))
      n = UV_MAX;
    else
      n = ((n + 3840)/30)*30;
    /* TODO: testing near 2^32-1 */

    prime_cache_sieve = sieve_erat30(n);

    if (prime_cache_sieve != 0)
      prime_cache_size = n;
    /* printf("size to %llu\n", prime_cache_size); fflush(stdout); */
  }
  MUTEX_UNLOCK(&primary_mutex_no_accessing);
}

/*
 * Get the size and a pointer to the cached prime sieve.
 * Returns the maximum sieved value in the sieve.
 * Allocates and sieves if needed.
 *
 * The sieve holds 30 numbers per byte, using a mod-30 wheel.
 */
UV get_prime_cache(UV n, const unsigned char** sieve)
{
  int prev_readers;

  if (sieve == 0) {
    if (prime_cache_size < n) {
      _fill_prime_cache(n, 0);
    }
    return prime_cache_size;
  }

  if (prime_cache_size < n)
    _fill_prime_cache(n, 0);

  /* TODO: We've got a problem.  If another thread does a memfree right here,
   * then we'll return a size less than n.  Everything technically works, but
   * there will be sieve croaks because they can't get enough primes.
   */

  MUTEX_LOCK(&primary_mutex_no_waiting);
    MUTEX_LOCK(&primary_mutex_counter);
      prev_readers = primary_number_of_readers;
      primary_number_of_readers++;
    MUTEX_UNLOCK(&primary_mutex_counter);
    if (prev_readers == 0) { MUTEX_LOCK(&primary_mutex_no_accessing); }
  MUTEX_UNLOCK(&primary_mutex_no_waiting);

  *sieve = prime_cache_sieve;
  return prime_cache_size;
}
void release_prime_cache(const unsigned char* mem) {
  int current_readers;
  MUTEX_LOCK(&primary_mutex_counter);
    primary_number_of_readers--;
    current_readers = primary_number_of_readers;
  MUTEX_UNLOCK(&primary_mutex_counter);
  if (current_readers == 0) { MUTEX_UNLOCK(&primary_mutex_no_accessing); }
}



#define SEGMENT_CHUNK_SIZE  UVCONST(256*1024*1024-8)
static unsigned char* prime_segment = 0;
static int prime_segment_is_available = 1;

unsigned char* get_prime_segment(UV *size) {
  unsigned char* mem;
  int use_prime_segment;

  MPUassert(size != 0, "get_prime_segment given null size pointer");
  MPUassert(mutex_init == 1, "segment mutex has not been initialized");

  MUTEX_LOCK(&segment_mutex);
    if (prime_segment_is_available) {
      prime_segment_is_available = 0;
      use_prime_segment = 1;
    } else {
      use_prime_segment = 0;
    }
  MUTEX_UNLOCK(&segment_mutex);

  if (use_prime_segment) {
    if (prime_segment == 0)
      New(0, prime_segment, SEGMENT_CHUNK_SIZE, unsigned char);
    *size = SEGMENT_CHUNK_SIZE;
    mem = prime_segment;
  } else {
    UV chunk_size = 64*1024*1024-8;
    New(0, mem, chunk_size, unsigned char);
    *size = chunk_size;
  }

  if (mem == 0)
    croak("Could not allocate %"UVuf" bytes for segment sieve", *size);

  return mem;
}
void release_prime_segment(unsigned char* mem) {
  MUTEX_LOCK(&segment_mutex);
    if (mem == prime_segment) {
      prime_segment_is_available = 1;
    } else {
      Safefree(mem);
    }
  MUTEX_UNLOCK(&segment_mutex);
}



void prime_precalc(UV n)
{
  if (!mutex_init) {
    MUTEX_INIT(&segment_mutex);
    MUTEX_INIT(&primary_mutex_no_waiting);
    MUTEX_INIT(&primary_mutex_no_accessing);
    MUTEX_INIT(&primary_mutex_counter);
    mutex_init = 1;
  }

  /* On initialization, make a few primes (2-30k using 1k memory) */
  if (n == 0)
    n = (1024-16)*30;
  get_prime_cache(n, 0);   /* Sieve to n */

  /* TODO: should we prealloc the segment here? */
}


void prime_memfree(void)
{
  MPUassert(mutex_init == 1, "segment mutex has not been initialized");

  MUTEX_LOCK(&primary_mutex_no_waiting);
  MUTEX_LOCK(&primary_mutex_no_accessing);
  MUTEX_UNLOCK(&primary_mutex_no_waiting);
  if ( (prime_cache_sieve != 0) ) {
    /* printf("size to 0  nreaders: %d\n", primary_number_of_readers); fflush(stdout); */
    Safefree(prime_cache_sieve);
    prime_cache_sieve = 0;
    prime_cache_size = 0;
  }
  MUTEX_UNLOCK(&primary_mutex_no_accessing);

  MUTEX_LOCK(&segment_mutex);
  /* Don't free if another thread is using it */
  if ( (prime_segment != 0) && (prime_segment_is_available) ) {
    Safefree(prime_segment);
    prime_segment = 0;
  }
  MUTEX_UNLOCK(&segment_mutex);

  prime_precalc(0);
}


void _prime_memfreeall(void)
{
  /* No locks.  We're shutting everything down. */
  if (mutex_init) {
    MUTEX_DESTROY(&segment_mutex);
    MUTEX_DESTROY(&primary_mutex_no_waiting);
    MUTEX_DESTROY(&primary_mutex_no_accessing);
    MUTEX_DESTROY(&primary_mutex_counter);
    mutex_init = 0;
  }
  if (prime_cache_sieve != 0)
    Safefree(prime_cache_sieve);
  prime_cache_sieve = 0;
  prime_cache_size = 0;

  if (prime_segment != 0)
    Safefree(prime_segment);
  prime_segment = 0;
}
