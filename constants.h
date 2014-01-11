#ifndef MPU_CONSTANTS_H
#define MPU_CONSTANTS_H

#include "EXTERN.h"
#include "perl.h"

#if BITS_PER_WORD == 32
 #define MPU_MAX_PRIME      UVCONST(4294967291)
 #define MPU_MAX_PRIME_IDX  UVCONST(203280221)
#else
 #define MPU_MAX_PRIME      UVCONST(18446744073709551557)
 #define MPU_MAX_PRIME_IDX  UVCONST(425656284035217743)
#endif


/****************************************************************************/
/* Configuration */

/* To avoid thrashing, sieve a little farther than needed */
#define _MPU_FILL_EXTRA_N (128*30)
/* The initial cache size.  30k primes per 1k of cache. */
#define _MPU_INITIAL_CACHE_SIZE ((4096-16)*30 - _MPU_FILL_EXTRA_N)

#endif
