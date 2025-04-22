#ifndef MPU_LUCAS_SEQ_H
#define MPU_LUCAS_SEQ_H

#include "ptypes.h"

extern void lucasuvmod(UV* U, UV* V, UV P, UV Q, UV k, UV n);

/* If caller has IV P and Q, send through ivmod(P,n) first */
extern UV lucasvmod(UV P, UV Q, UV k, UV n);
extern UV lucasumod(UV P, UV Q, UV k, UV n);

/* No modulo, returns false if it overflows. */
extern bool lucasuv(IV* U, IV* V, IV P, IV Q, UV k);

/* Legacy version, also returns Qk */
extern void lucas_seq(UV* U, UV* V, UV* Qk,  UV n, IV P, IV Q, UV k);

#endif
