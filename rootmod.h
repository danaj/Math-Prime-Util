#ifndef MPU_ROOTMOD_H
#define MPU_ROOTMOD_H

#include "ptypes.h"

extern bool sqrtmodp(UV *r, UV a, UV p);        /* sqrt(a) mod p */
extern bool sqrtmod(UV *r, UV a, UV n);         /* sqrt(a) mod n */

extern bool rootmodp(UV *r, UV a, UV k, UV p);  /* a^(1/k) mod p */
extern bool rootmod(UV *r, UV a, UV k, UV n);   /* a^(1/k) mod n */

extern UV* allsqrtmod(UV* nroots, UV a, UV n);        /* all results */
extern UV* allrootmod(UV* nroots, UV a, UV k, UV n);  /* all results */

#endif
