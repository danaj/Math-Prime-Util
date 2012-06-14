#ifndef MPU_CACHE_H
#define MPU_CACHE_H

#include "EXTERN.h"
#include "perl.h"

extern UV  get_prime_cache_size(void);
extern UV  get_prime_cache(UV n, const unsigned char** sieve);
extern void free_prime_cache(unsigned char* sieve);

extern void  prime_precalc(UV x);
extern void  prime_memfree(void);

extern unsigned char* get_prime_segment(UV* size);
extern void free_prime_segment(unsigned char* segment);

#endif
