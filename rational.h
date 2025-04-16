#ifndef MPU_RATIONAL_H
#define MPU_RATIONAL_H

#include "ptypes.h"

extern int contfrac(UV** cfrac, UV *rem, UV num, UV den);

extern int next_calkin_wilf(UV* num, UV* den);
extern int next_stern_brocot(UV* num, UV* den);

extern UV calkin_wilf_n(UV num, UV den);
extern UV stern_brocot_n(UV num, UV den);

extern int nth_calkin_wilf(UV* num, UV* den, UV n);
extern int nth_stern_brocot(UV* num, UV* den, UV n);

extern UV nth_stern_diatomic(UV n);

#endif
