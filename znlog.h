#ifndef MPU_ZNLOG_H
#define MPU_ZNLOG_H

#include "ptypes.h"

/* Generic znlog returns k that solves a = g^k mod p */
extern UV znlog(UV a, UV g, UV p);

/* znlog given prime gorder = znorder(g,p) */
extern UV znlog_solve(UV a, UV g, UV p, UV gorder);

#endif
