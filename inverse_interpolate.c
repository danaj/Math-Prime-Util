#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "inverse_interpolate.h"

#define LINEAR_INTERP(n, lo, hi, rlo, rhi) \
  (lo + (UV) (((double)(n-rlo) * (double)(hi-lo) / (double)(rhi-rlo))+0.5))

UV inverse_interpolate(UV lo, UV hi, UV n, UV k,
                       UV (*func)(UV mid, UV k),
                       UV threshold) {
  UV i, mid, rlo, rhi;

  rlo = (*func)(lo,k);
  rhi = (*func)(hi,k);

  mid = LINEAR_INTERP(n,lo,hi,rlo,rhi);

  for (i = 1; i <= 12 && lo < hi && rhi != n; i++) {
    if (mid > lo && mid < hi) {
      UV rmid = (*func)(mid,k);
      if (rmid < n) { lo = mid; rlo = rmid; }
      else          { hi = mid; rhi = rmid; }
    }
    mid = LINEAR_INTERP(n,lo,hi,rlo,rhi);
  }
  /* Try to pull in */
  if (rlo < n-1) {
    mid = LINEAR_INTERP(n-1,lo,hi,rlo,rhi);
    if (mid > lo && mid < hi) {
      if ((*func)(mid,k) < n) lo = mid+1;
      else                    hi = mid;
    }
  }

  while (lo < hi) {
    if (rhi == n && (hi-lo < threshold)) break;
    mid = lo + ((hi-lo)>>1);
    if ((*func)(mid,k) < n) lo = mid+1;
    else                    hi = mid;
  }
  return hi;
}
