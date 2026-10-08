#ifndef MPU_TINYSIQS128_H
#define MPU_TINYSIQS128_H

#include "factor128.h"

#if HAVE_FACTOR128
/* Factor using SIQS (or shortcuts for small inputs).
 *
 * 0 and 1 return themselves.
 * For other inputs, the returned value is:
 *    N      input was prime (BPSW when needed)
 *    0      failed to factor, we don't expect this
 *    F      a proper divisor of n, 1 < F < N
 * F divides N but is not guaranteed to be prime.
 */
extern uint128_t tinysiqs128(uint128_t n);

/* Factor as above, but never runs the primality test.
 *
 * This is meant for callers who have already done a primality test
 * and don't want the factoring routine to run it again.  If primes
 * do come in, they may or may not run the full SIQS, and the return
 * value doesn't distinguish between prime and failure to factor.
 *
 * 0 and 1 return themselves.
 * For other inputs, the returned value is:
 *    N      input was prime
 *    0      either failed to factor or input was prime
 *    F      a proper divisor of n, 1 < F < N
 * F divides N but is not guaranteed to be prime.
 */
extern uint128_t tinysiqs128_composite(uint128_t n);
#endif

#endif
