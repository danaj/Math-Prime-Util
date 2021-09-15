#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "inverse_interpolate.h"
#include "util.h"

/* Return x with v(x)=(*func)(x,k) s.t. either of:
 *    1.  v(x) == n  and v(x-1-threshold) < n
 *    2.  v(x) < n   and v(x+1) > n
 */

#define LINEAR_INTERP(n, lo, hi, rlo, rhi) \
  (lo + (UV) (((double)(n-rlo) * (double)(hi-lo) / (double)(rhi-rlo))+0.5))

#define CALLBACK(n)  ((funck) ? (*funck)(n,k) : (*func)(n))

#if 0   /* Debugging return, checking the conditions above. */
#define RETURNI(x) \
  { \
    UV v = x; \
    UV rv = CALLBACK(v); \
    /* printf("v  %lu    rv %lu   n %lu\n",v,rv,n); */\
    MPUassert( rv <= n, "BAD INTERP  v > n" ); \
    if (rv == n) { \
      if (v > threshold) { \
        /* printf("threshold %lu v %lu    func(%lu) = %lu\n", threshold, v, v-1-threshold, CALLBACK(v-1-threshold)); */\
        MPUassert( CALLBACK(v-1-threshold) < n, "BAD INTERP  v-1-thresh >= n" ); \
      } \
    } else { \
      MPUassert( CALLBACK(v+1) > n, "BAD INTERP  v+1 <= n" ); \
    } \
    return v; \
  }
#else
  #define RETURNI(x) { return x; }
#endif

static UV _inverse_interpolate(UV lo, UV hi, UV n,
                               UV k, UV (*funck)(UV mid, UV k),
                               UV (*func)(UV mid),
                               UV threshold) {
  UV mid, rlo, rhi, rmid;

  rlo = CALLBACK(lo);
  if (rlo == n)  RETURNI(lo);  /* Possible bad limit */
  rhi = CALLBACK(hi);

  /* printf("n %lu\n   lo %lu  rlo  %lu\n   hi %lu  rhi  %lu\n",n,lo,rlo,hi,rhi);*/
  MPUassert(rlo <= n && rhi >= n, "interpolation: bad initial limits");

  if ((hi-lo) == 1)
    RETURNI( (rlo == n) ? lo : hi );

  mid = LINEAR_INTERP(n,lo,hi,rlo,rhi);
  if (mid == lo) mid++;  else if (mid == hi) mid--;
  rmid = CALLBACK(mid);

  while (lo < mid && mid < hi && rmid != n && rlo != n) {
    if (rmid > n) {
      hi = mid;  rhi = rmid;
      mid -= (UV)((double)(mid-lo) * (double)(rmid-n)/(double)(rmid-rlo) + 0.5);
    } else {
      lo = mid;  rlo = rmid;
      mid += (UV)((double)(hi-mid) * (double)(n-rmid)/(double)(rhi-rmid) + 0.5);
    }
    if (mid == lo) mid++;  else if (mid == hi) mid--;
    rmid = CALLBACK(mid);
  }

  MPUassert(lo <= mid && mid <= hi, "interpolation: range error");
  MPUassert(rlo <= n && rhi >= n, "interpolation: bad initial interpolation");

  if (rlo == n)
    RETURNI(lo);

  if ((hi-lo) == 1)
    RETURNI( (rlo == n || (rlo < n && rhi > n)) ? lo : hi );
  if ((hi-lo) == 2) {
    if (hi == mid || lo == mid) croak("bad mid assumption");
    if      (rlo  == n || (rlo  < n && rmid > n)) { RETURNI(lo);  }
    else if (rmid == n || (rmid < n && rhi  > n)) { RETURNI(mid); }
    else                                          { RETURNI(hi);  }
  }

  /* printf("   lo %lu  rlo  %lu\n  mid %lu  rmid %lu\n   hi %lu  rhi  %lu\n",lo,rlo,mid,rmid,hi,rhi); */
  MPUassert(rmid == n, "interpolation: bad rmid");

  /* mid / rmid are correct, but we need to find the least value */

  if ((mid-lo) <= threshold)
    RETURNI(mid);

  /* hi and rhi will always remain a correct value, we're pulling lo in */
  hi = mid; rhi = rmid;

  /* It's possible 'lo' is very far away.  Try to pull 'lo' in a lot. */
  if (rlo < n-1) {
    mid = mid - (UV)((hi-lo) * 0.10 + 1);
    if (mid == lo) mid++;  else if (mid == hi) mid--;
    rmid = CALLBACK(mid);
    if (rmid < n) lo = mid;
    else          hi = mid;
  }

  /* Binary search until within threshold */
  while (lo < hi && (hi-lo) >= threshold) {
    mid = lo + ((hi-lo)>>1);
    if (CALLBACK(mid) < n) lo = mid+1;
    else                    hi = mid;
  }
  RETURNI(hi);
}


UV inverse_interpolate(UV lo, UV hi, UV n, UV (*func)(UV mid), UV threshold) {
  return _inverse_interpolate(lo,hi,n,0,0,func,threshold);
}

UV inverse_interpolate_k(UV lo, UV hi, UV n, UV k, UV (*funck)(UV mid, UV k), UV threshold) {
  return _inverse_interpolate(lo,hi,n,k,funck,0,threshold);
}


/******************************************************************************/


UV interpolate_with_approx(UV n,
                           UV *gcount,
                           UV tol,
                           UV (*fnth)(UV n),
                           UV (*fcnt)(UV n),
                           int (*fis)(UV n)   /* optional */
                          ) {
  UV guess, gn, count, ming = 0, maxg = UV_MAX;

  guess = (*fnth)(n);
  for (gn = 2; gn < 20; gn++) {
    IV adjust;
    MPUverbose(2, "  interp  %"UVuf"-th is around %"UVuf" ... ", n, guess);
    count = (*fcnt)(guess);
    MPUverbose(2, "(%"IVdf")\n", (IV)(n-count));
    /* Stop guessing if within our tolerance */
    if (n==count || (n>count && n-count < tol) || (n<count && count-n < tol)) break;
    /* Determine how far off we think we are */
    adjust = (IV) ((fnth)(n) - (fnth)(count));
    /* When computing new guess, ensure we don't overshoot.  Rarely used. */
    if (count <= n && guess > ming) ming = guess;   /* Previous guesses */
    if (count >= n && guess < maxg) maxg = guess;
    guess += adjust;
    if (guess <= ming || guess >= maxg) MPUverbose(2, "  fix min/max for %"UVuf"\n",n);
    if (guess <= ming) guess = ming + tol - 1;
    if (guess >= maxg) guess = maxg - tol + 1;
    /* TODO: if min/max dist is small, split the difference. */
  }
  if (gn == 20) count = (*fcnt)(guess);

  if (fis) {
    if (count < n) {

      /* Increase count one at a time if needed */
      for ( ; count < n; count++)
        while (!(*fis)(++guess))
          ;

    } else if (count >= n) {

      /* Make sure this is the least value at this count */
      while (!(*fis)(guess))  guess--;
      /* Reduce count one at a time if needed */
      for ( ; count > n; count--)
        while (!(*fis)(--guess))
          ;

    }
  }

  if (gcount) *gcount = count;
  return guess;
}
