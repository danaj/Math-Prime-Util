#ifndef MPU_FINVERSE_H
#define MPU_FINVERSE_H

#include "ptypes.h"

extern UV inverse_interpolate(UV lo, UV hi, UV n, UV (*func)(UV mid), UV threshold);

extern UV inverse_interpolate_k(UV lo, UV hi, UV n, UV k, UV (*func)(UV mid, UV k), UV threshold);


/*
 * We need to be given:
 *   n          we're trying to find nth_xxx(n)
 *   gcount     we'll fill this in with the exact count at the return value
 *   tol        how close we care about getting
 *   fnth(n)    a callback for nth_xxx_approx(n) (approx nth)
 *   fcnt(g)    a callback for xxx_count(g)      (exact count at g)
 *   fis(g)     a callback for is_xxx(g)
 * returns g, the nth_xxx(n) value
 *
 * if fis(g) is set, then the result is exact.  After we narrow in using the
 * count and approximations, we step to the correct spot using the is_xxx(g)
 * function.  If it is not set, then we will return once within the tolerance,
 * with the exact count returned as well.  The caller may then want to use
 * fast next / prev functions to quickly step to the right spot.
 *
 *
 * Each caller can have slightly different optimization considerations.
 * For example, perfect_power has a super fast count function, so we want
 * to narrow down quickly.  It also has fast next/prev functions.
 */
extern UV interpolate_with_approx(UV n, UV *gcount, UV tol, UV (*fnth)(UV n), UV (*fcnt)(UV n), bool (*fis)(UV n));

#endif
