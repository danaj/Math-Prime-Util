#ifndef MPU_FACTOR128_H
#define MPU_FACTOR128_H

#include "ptypes.h"

/* factor128 works if we have uint64_t and uint128_t.
 * Caution: callers that use 64-bit or 128-bit factor results must not
 * assume they fit in a UV. */
#define HAVE_FACTOR128 (HAVE_UINT64 && HAVE_UINT128)

#if HAVE_FACTOR128

/* Conversions of uint128_t <=> string */
extern int  u128_to_str(char str[40], uint128_t n);
extern bool str_to_u128(uint128_t *out, const char *s, size_t len);

#define MPU_MAX_DFACTORS128  27

typedef struct {
  uint128_t n;
  uint64_t  f[MPU_MAX_DFACTORS128];  /* 64-bit prime factors */
  uint8_t   e[MPU_MAX_DFACTORS128];  /* exponents of f[] */
  uint16_t  nfactors;                /* number of entries in f[]/e[] */
  uint128_t flarge;                  /* a prime > UINT64_MAX with exp 1;
                                        0 if no such factor */
} factored128_t;

/* Factor n into nf.  Croaks if an internal split fails. */
extern void factorintp128(factored128_t *nf, uint128_t n);

MAYBE_UNUSED static INLINE factored128_t factorint128(uint128_t n)
  { factored128_t nf; factorintp128(&nf, n); return nf; }

extern bool        is_prime128(uint128_t n);
extern bool        is_bpsw128(uint128_t n);
extern signed char moebius128(uint128_t n);
extern bool        is_semiprime128(uint128_t n);
extern uint128_t   muladdmod128_s(uint128_t a, uint128_t b, uint128_t c,
                                  uint128_t n, int sub);

/* Simple helper functions */
extern uint32_t    factored128p_total_factors(const factored128_t *nf);
extern uint32_t    factored128p_distinct_factors(const factored128_t *nf);
extern bool        factored128p_is_square_free(const factored128_t *nf);
extern signed char factored128p_moebius(const factored128_t *nf);

#endif /* HAVE_FACTOR128 */

#endif /* MPU_FACTOR128_H */
