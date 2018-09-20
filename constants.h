#ifndef MPU_CONSTANTS_H
#define MPU_CONSTANTS_H

#include "ptypes.h"

#if BITS_PER_WORD == 32
 #define MPU_MAX_PRIME          UVCONST(4294967291)
 #define MPU_MAX_PRIME_IDX      UVCONST(203280221)
 #define MPU_MAX_TWIN_PRIME     UVCONST(4294965839)
 #define MPU_MAX_TWIN_PRIME_IDX UVCONST(12739574)
 #define MPU_MAX_RMJN_PRIME     UVCONST(4294967279)
 #define MPU_MAX_RMJN_PRIME_IDX UVCONST(98182656)
 #define MPU_MAX_SEMI_PRIME     UVCONST(4294967294)
 #define MPU_MAX_SEMI_PRIME_IDX UVCONST(658662065)
#else
 #define MPU_MAX_PRIME          UVCONST(18446744073709551557)
 #define MPU_MAX_PRIME_IDX      UVCONST(425656284035217743)
 #define MPU_MAX_TWIN_PRIME     UVCONST(18446744073709550771)
 #define MPU_MAX_TWIN_PRIME_IDX UVCONST(12975810317986308)    /* Approx */
 #define MPU_MAX_RMJN_PRIME     UVCONST(18446744073709550771) /* Not correct */
 #define MPU_MAX_RMJN_PRIME_IDX UVCONST(12975810317986308)    /* Not correct */
 #define MPU_MAX_SEMI_PRIME     UVCONST(18446744073709551601)
 #define MPU_MAX_SEMI_PRIME_IDX UVCONST(1600000000000000000)  /* Max est */
#endif


/****************************************************************************/
/* Configuration */

/* To avoid thrashing, sieve a little farther than needed */
#define _MPU_FILL_EXTRA_N (128*30)
/* The initial cache size.  30k primes per 1k of cache. */
#define _MPU_INITIAL_CACHE_SIZE ((4096-16)*30 - _MPU_FILL_EXTRA_N)
/* Where to start using LMO instead of segment sieve */
#define _MPU_LMO_CROSSOVER 66000000

#endif
