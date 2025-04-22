#ifndef MPU_RANDOM_PRIME_H
#define MPU_RANDOM_PRIME_H

#include "ptypes.h"

extern UV random_safe_prime(void* ctx, UV bits);
extern UV random_nbit_prime(void* ctx, UV bits);
extern UV random_ndigit_prime(void* ctx, UV digits);
extern UV random_prime(void* ctx, UV lo, UV hi);

extern bool is_mr_random(void* ctx, UV n, UV k);

extern UV random_semiprime(void* ctx, UV bits);
extern UV random_unrestricted_semiprime(void* ctx, UV bits);

#endif
