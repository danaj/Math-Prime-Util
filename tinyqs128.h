#ifndef MPU_TINYQS128_H
#define MPU_TINYQS128_H

#include "factor128.h"

#if HAVE_FACTOR128
/* Return a factor of n up to 128 bits, or zero on failure. */
extern uint128_t tinyqs128(uint128_t n);
#endif

#endif
