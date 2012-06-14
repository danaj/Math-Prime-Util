#ifndef MPU_CACHE_H
#define MPU_CACHE_H

#include "EXTERN.h"
#include "perl.h"

extern UV  get_prime_cache_size(void);
extern UV  get_prime_cache(UV n, const unsigned char** sieve);

extern void  prime_memfree(void);

#define SEGMENT_CHUNK_SIZE  UVCONST(262144)
extern unsigned char* get_prime_segment(void);

#endif
