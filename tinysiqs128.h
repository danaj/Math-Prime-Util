#ifndef MPU_TINYSIQS128_H
#define MPU_TINYSIQS128_H

#include "factor128.h"

#if HAVE_FACTOR128
/* Return a proper nontrivial factor of n, or zero when none was found. */
extern uint128_t tinysiqs128(uint128_t n);
#endif

#endif
