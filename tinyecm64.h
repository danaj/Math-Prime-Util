/* tinyecm64 1.0.1
 *
 * ECM for uint64_t inputs.  Portable C, no GMP, no heap allocation.
 */
#ifndef TINYECM64_H
#define TINYECM64_H
#include "ptypes.h"

#if HAVE_UINT64

/* Return a proper divisor or zero after a finite main budget and rescue.
 * For odd n>3, advance *rng_state once and use a private stream for the
 * curves.  Zero is a valid seed; the pointer must be nonnull.
 * Nonzero try_small adds a short small-factor prefix, not trial division.
 * All mutable state is local.  Callers retain their own fallback policy. */
uint64_t tinyecm64_rng(uint64_t n, int try_small, uint64_t *rng_state);

/* As above but uses n as the rng state.  Simpler and deterministic. */
uint64_t tinyecm64(uint64_t n, int try_small);

#endif

#endif
