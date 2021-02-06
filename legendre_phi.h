#ifndef MPU_LEGENDRE_PHI_H
#define MPU_LEGENDRE_PHI_H

#include "ptypes.h"

extern UV legendre_phi(UV n, UV a);

extern uint32_t tiny_phi_max_a(void);
extern UV tiny_phi(UV n, uint32_t a);

#endif
