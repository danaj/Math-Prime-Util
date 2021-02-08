#ifndef MPU_LEGENDRE_PHI_H
#define MPU_LEGENDRE_PHI_H

#include "ptypes.h"

extern UV legendre_phi(UV n, UV a);

/* These have almost no overhead, so make good stopping points */

/* Extremely fast, a <= ~6 */
extern uint32_t tiny_phi_max_a(void);
extern UV tiny_phi(UV n, uint32_t a);
/* Fast, a <= ~15 */
extern uint32_t small_phi_max_a(void);
extern UV small_phi(UV n, uint32_t a);

/* For doing many calls with small values */
extern void* prepare_cached_legendre_phi(UV n, UV a);
extern UV cached_legendre_phi(void* cache, UV n, UV a);
extern void destroy_cached_legendre_phi(void* cache);

#endif
