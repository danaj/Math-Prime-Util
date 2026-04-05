#ifndef MPU_FACTOR128_H
#define MPU_FACTOR128_H

#include "ptypes.h"

#if BITS_PER_WORD == 64 && HAVE_UINT128

#define MPU_MAX_DFACTORS128  27

typedef struct {
  uint128_t n;
  UV        f[MPU_MAX_DFACTORS128];  /* small prime factors (<= UV_MAX) */
  uint8_t   e[MPU_MAX_DFACTORS128];  /* exponents of f[] */
  uint16_t  nfactors;                /* number of entries in f[]/e[] */
  uint128_t flarge;                  /* a prime > UV_MAX, exponent always 1;
                                        0 if no such factor */
} factored128_t;

/* Factor n into nf.  n must be > 1.  Returns false if cannot factor. */
extern bool factorintp128(factored128_t *nf, uint128_t n);

static INLINE factored128_t factorint128(uint128_t n)
  { factored128_t nf; factorintp128(&nf, n); return nf; }

#endif /* BITS_PER_WORD == 64 && HAVE_UINT128 */

#endif /* MPU_FACTOR128_H */
