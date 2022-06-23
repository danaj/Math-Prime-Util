#ifndef MPU_FINVERSE_H
#define MPU_FINVERSE_H

#include "ptypes.h"

extern UV inverse_interpolate(UV lo, UV hi, UV n, UV (*func)(UV mid), UV threshold);

extern UV inverse_interpolate_k(UV lo, UV hi, UV n, UV k, UV (*func)(UV mid, UV k), UV threshold);


extern UV interpolate_with_approx(UV n, UV *gcount, UV tol, UV (*fnth)(UV n), UV (*fcnt)(UV n), int (*fis)(UV n));

#endif
