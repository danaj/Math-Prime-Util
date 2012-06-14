#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "cache.h"
#include "sieve.h"


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



static unsigned char* prime_segment = 0;
unsigned char* get_prime_segment(void) {
  if (prime_segment == 0)
    New(0, prime_segment, SEGMENT_CHUNK_SIZE, unsigned char);
  if (prime_segment == 0)
    croak("Could not allocate %"UVuf" bytes for segment sieve", SEGMENT_CHUNK_SIZE);
  return prime_segment;
}
void free_prime_segment(void) {
  if (prime_segment != 0)
    Safefree(prime_segment);
  prime_segment = 0;
}



void _prime_memfreeall(void)
{
  if (prime_cache_sieve != 0)
    Safefree(prime_cache_sieve);
  prime_cache_sieve = 0;
  prime_cache_size = 0;

  free_prime_segment();
}

void prime_memfree(void)
{
  _prime_memfreeall();

  prime_precalc(0);
}

