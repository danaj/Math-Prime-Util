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
 * Since we're trying to be thread-safe (and ideally allow a good deal
 * of concurrency), it is imperative these be used correctly.  You need
 * to call the get method, do your stuff, then call free.  Do *not* return
 * out of your function or croak without calling free.
 */

static int mutex_init = 0;
#ifdef USE_ITHREADS
static perl_mutex segment_mutex;
static perl_mutex primary_mutex;
#endif

static unsigned char* prime_cache_sieve = 0;
static UV             prime_cache_size = 0;

/* Get the maximum sieved value of the cached prime sieve. */
UV get_prime_cache_size(void) { return prime_cache_size; }

/*
 * Get the size and a pointer to the cached prime sieve.
 * Returns the maximum sieved value in the sieve.
 * Allocates and sieves if needed.
 *
 * The sieve holds 30 numbers per byte, using a mod-30 wheel.
 */
UV get_prime_cache(UV n, const unsigned char** sieve)
{
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
  }

  if (sieve != 0)
    *sieve = prime_cache_sieve;
  return prime_cache_size;
}



#define SEGMENT_CHUNK_SIZE  UVCONST(262144)
static unsigned char* prime_segment = 0;

unsigned char* get_prime_segment(UV *size) {
  MPUassert(size != 0, "get_prime_segment given null size pointer");
  MPUassert(mutex_init == 1, "segment mutex has not been initialized");
  MUTEX_LOCK(&segment_mutex);
  if (prime_segment == 0)
    New(0, prime_segment, SEGMENT_CHUNK_SIZE, unsigned char);
  if (prime_segment == 0) {
    MUTEX_UNLOCK(&segment_mutex);
    croak("Could not allocate %"UVuf" bytes for segment sieve", SEGMENT_CHUNK_SIZE);
  }
  *size = SEGMENT_CHUNK_SIZE;
  return prime_segment;
}
void free_prime_segment(unsigned char* mem) {
  /* Thanks for letting us know you're done. */
  MUTEX_UNLOCK(&segment_mutex);
}



void prime_precalc(UV n)
{
  if (!mutex_init) {
    MUTEX_INIT(&segment_mutex);
    MUTEX_INIT(&primary_mutex);
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

  if (prime_cache_sieve != 0) {
    MUTEX_LOCK(&primary_mutex);
    Safefree(prime_cache_sieve);
    prime_cache_sieve = 0;
    prime_cache_size = 0;
    MUTEX_UNLOCK(&primary_mutex);
  }

  if (prime_segment != 0) {
    MUTEX_LOCK(&segment_mutex);
    Safefree(prime_segment);
    prime_segment = 0;
    MUTEX_UNLOCK(&segment_mutex);
  }

  prime_precalc(0);
}


void _prime_memfreeall(void)
{
  /* No locks.  We're shutting everything down. */
  if (mutex_init) {
    MUTEX_DESTROY(&segment_mutex);
    MUTEX_DESTROY(&primary_mutex);
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
