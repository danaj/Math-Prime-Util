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
  UV mid, rlo, rhi, rmid;

  rlo = (*func)(lo,k);
  if (rlo == n)  return lo;  /* Possible bad limit */
  rhi = (*func)(hi,k);

  MPUassert(rlo <= n && rhi >= n, "interpolation: bad initial limits");

  if ((hi-lo) == 1)
    return (rlo == n) ? lo : hi;

  mid = LINEAR_INTERP(n,lo,hi,rlo,rhi);
  if (mid == lo) mid++;  else if (mid == hi) mid--;
  rmid = (*func)(mid,k);

  while (lo < mid && mid < hi && rmid != n && rlo != n) {
    if (rmid > n) {
      hi = mid;  rhi = rmid;
      mid -= (UV)((double)(mid-lo) * (double)(rmid-n)/(double)(rmid-rlo) + 0.5);
    } else {
      lo = mid;  rlo = rmid;
      mid += (UV)((double)(hi-mid) * (double)(n-rmid)/(double)(rhi-rmid) + 0.5);
    }
    if (mid == lo) mid++;  else if (mid == hi) mid--;
    rmid = (*func)(mid,k);
  }

  MPUassert(lo <= mid && mid <= hi, "interpolation: range error");
  MPUassert(rlo <= n && rhi >= n, "interpolation: bad initial interpolation");

  if (rlo == n) return lo;

  if (hi-lo <= 2) {
    if ((hi-lo) == 2 && (hi == mid || lo == mid)) croak("bad mid assumption");
    return (rlo == n) ? lo : (mid <= hi && rmid == n) ? mid : hi;
  }

  MPUassert(rmid == n, "interpolation: bad rmid");

  /* mid / rmid are correct, but it coulud be a one to many mapping
   * and they want the least value.
   *
   * We will return x such that x <= n and x+1+threshold > n
   */

  if ((mid-lo) <= threshold)
    return mid;

  /* hi and rhi will always remain a correct value, we're pulling lo in */
  hi = mid; rhi = rmid;

  /* It's possible 'lo' is very far away.  Try to pull 'lo' in a lot. */
  if (rlo < n-1) {
    mid = mid - (UV)((hi-lo) * 0.10 + 1);
    if (mid == lo) mid++;  else if (mid == hi) mid--;
    rmid = (*func)(mid,k);
    if (rmid < n) lo = mid;
    else          hi = mid;
  }

  /* Binary search until within threshold */
  while (lo < hi && (hi-lo) >= threshold) {
    mid = lo + ((hi-lo)>>1);
    if ((*func)(mid,k) < n) lo = mid+1;
    else                    hi = mid;
  }
  return hi;
}
