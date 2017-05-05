#ifndef MPU_RANDOM_PRIME_H
#define MPU_RANDOM_PRIME_H

#include "ptypes.h"

extern UV random_nbit_prime(UV b);
extern UV random_ndigit_prime(UV d);
extern UV random_prime(UV lo, UV hi);

extern int is_mr_random(UV n, UV k);

extern UV random_semiprime(UV b, int type);

#endif
