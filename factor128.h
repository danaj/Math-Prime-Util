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

/* Factor n into nf.  Croaks if an internal split fails. */
extern void factorintp128(factored128_t *nf, uint128_t n);

MAYBE_UNUSED static INLINE factored128_t factorint128(uint128_t n)
  { factored128_t nf; factorintp128(&nf, n); return nf; }

extern signed char moebius128(uint128_t n);
extern bool        is_semiprime128(uint128_t n);

/* Simple helper functions */
extern uint32_t    factored128p_total_factors(const factored128_t *nf);
extern uint32_t    factored128p_distinct_factors(const factored128_t *nf);
extern bool        factored128p_is_square_free(const factored128_t *nf);
extern signed char factored128p_moebius(const factored128_t *nf);

#endif /* BITS_PER_WORD == 64 && HAVE_UINT128 */

#endif /* MPU_FACTOR128_H */
