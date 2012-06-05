#ifndef MPU_UTIL_H
#define MPU_UTIL_H

#include "EXTERN.h"
#include "perl.h"

extern int is_prime(UV x);
extern UV  next_trial_prime(UV x);
extern UV  next_prime(UV x);
extern UV  prev_prime(UV x);

extern UV  prime_count_lower(UV x);
extern UV  prime_count_upper(UV x);
extern UV  prime_count_approx(UV x);
extern UV  prime_count(UV x);

extern UV  nth_prime_lower(UV n);
extern UV  nth_prime_upper(UV x);
extern UV  nth_prime_approx(UV x);
extern UV  nth_prime(UV x);

#define SEGMENT_CHUNK_SIZE 262144
extern unsigned char* get_prime_segment(void);
extern void           free_prime_segment(void);

#endif
