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

/* To avoid thrashing, sieve a little farther than we absolutely need to. */
#define FILL_EXTRA_N (128*30)

/* Erase the primary cache and fill up to n. */
static void _erase_and_fill_prime_cache(UV n) {
  UV padded_n;
  /* Note: You need to handle mutexes around this.
   *   MUTEX_LOCK(&primary_mutex_no_waiting);
   *   MUTEX_LOCK(&primary_mutex_no_accessing);
   *   MUTEX_UNLOCK(&primary_mutex_no_waiting);
   *   _fill_prime_cache(n);
   *   MUTEX_UNLOCK(&primary_mutex_no_accessing);
   */

  if (n >= (UV_MAX-FILL_EXTRA_N))
    padded_n = UV_MAX;
  else
    padded_n = ((n + FILL_EXTRA_N)/30)*30;

  /* If new size isn't larger or smaller, then we're done. */
  if (prime_cache_size == padded_n)
    return;

  if (prime_cache_sieve != 0)
    Safefree(prime_cache_sieve);
  prime_cache_sieve = 0;
  prime_cache_size = 0;

  if (n > 0) {
    prime_cache_sieve = sieve_erat30(padded_n);
    if (prime_cache_sieve != 0)
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
  int prev_readers;
  int i_hold_access_lock = 0;

  if (sieve == 0) {
    if (prime_cache_size < n) {
      MUTEX_LOCK(&primary_mutex_no_waiting);
      MUTEX_LOCK(&primary_mutex_no_accessing);
      MUTEX_UNLOCK(&primary_mutex_no_waiting);
      _erase_and_fill_prime_cache(n);
      MUTEX_UNLOCK(&primary_mutex_no_accessing);
    }
    return prime_cache_size;
  }

  MUTEX_LOCK(&primary_mutex_no_waiting);
    /* If we need more primes, get the access lock right now */
    if (prime_cache_size < n) {
      MUTEX_LOCK(&primary_mutex_no_accessing);
      i_hold_access_lock = 1;
    }

    MUTEX_LOCK(&primary_mutex_counter);
      prev_readers = primary_number_of_readers;
      primary_number_of_readers++;
    MUTEX_UNLOCK(&primary_mutex_counter);

    if ( (prev_readers == 0) && (!i_hold_access_lock) ) {
      MUTEX_LOCK(&primary_mutex_no_accessing);
      i_hold_access_lock = 1;
    }
    if (prime_cache_size < n) {
      _erase_and_fill_prime_cache(n);
    }
  MUTEX_UNLOCK(&primary_mutex_no_waiting);

  MPUassert(prime_cache_size >= n, "prime cache is too small!");

  *sieve = prime_cache_sieve;
  return prime_cache_size;
#else
  if (prime_cache_size < n)
    _erase_and_fill_prime_cache(n);
  if (sieve != 0)
    *sieve = prime_cache_sieve;
  return prime_cache_size;
#endif
}

void release_prime_cache(const unsigned char* mem) {
#ifdef USE_ITHREADS
  int current_readers;
  MUTEX_LOCK(&primary_mutex_counter);
    primary_number_of_readers--;
    current_readers = primary_number_of_readers;
  MUTEX_UNLOCK(&primary_mutex_counter);
  if (current_readers == 0) { MUTEX_UNLOCK(&primary_mutex_no_accessing); }
#endif
}



/* The segment everyone is trying to share */
#define PRIMARY_SEGMENT_CHUNK_SIZE    UVCONST(256*1024*1024-16)
static unsigned char* prime_segment = 0;
static int prime_segment_is_available = 1;
/* If that's in use, malloc a new one of this size */
#define SECONDARY_SEGMENT_CHUNK_SIZE  UVCONST( 64*1024*1024-16)

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
      New(0, prime_segment, PRIMARY_SEGMENT_CHUNK_SIZE, unsigned char);
    *size = PRIMARY_SEGMENT_CHUNK_SIZE;
    mem = prime_segment;
  } else {
    New(0, mem, SECONDARY_SEGMENT_CHUNK_SIZE, unsigned char);
    *size = SECONDARY_SEGMENT_CHUNK_SIZE;
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



#define INITIAL_CACHE_SIZE ((1024-16)*30 - FILL_EXTRA_N)
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
    n = INITIAL_CACHE_SIZE;
  get_prime_cache(n, 0);   /* Sieve to n */

  /* TODO: should we prealloc the segment here? */
}


void prime_memfree(void)
{
  MPUassert(mutex_init == 1, "cache mutexes have not been initialized");

  MUTEX_LOCK(&segment_mutex);
  /* Don't free if another thread is using it */
  if ( (prime_segment != 0) && (prime_segment_is_available) ) {
    Safefree(prime_segment);
    prime_segment = 0;
  }
  MUTEX_UNLOCK(&segment_mutex);

  MUTEX_LOCK(&primary_mutex_no_waiting);
  MUTEX_LOCK(&primary_mutex_no_accessing);
  MUTEX_UNLOCK(&primary_mutex_no_waiting);
    /* Put primary cache back to initial state */
    _erase_and_fill_prime_cache(INITIAL_CACHE_SIZE);
  MUTEX_UNLOCK(&primary_mutex_no_accessing);
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
