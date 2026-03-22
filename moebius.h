#ifndef MPU_MOEBIUS_H
#define MPU_MOEBIUS_H

#include "ptypes.h"

extern int          moebius(UV n);
extern signed char* range_moebius(UV lo, UV hi);
extern IV           mertens(UV n);
extern int          liouville(UV n);
extern IV           sumliouville(UV n);
extern signed char* range_liouville(UV lo, UV hi);

#endif
