#ifndef MPU_REAL_H
#define MPU_REAL_H

#include "ptypes.h"

extern NV Ei(NV x);
extern NV Li(NV x);
extern long double ld_inverse_li(long double x);
extern long double ld_riemann_zeta(long double x);
extern long double RiemannR(long double x, long double eps);
extern NV lambertw(NV k);
extern UV inverse_li(UV x);
extern UV inverse_R(UV x);

extern NV  chebyshev_psi(UV n);
extern NV  chebyshev_theta(UV n);

#endif
