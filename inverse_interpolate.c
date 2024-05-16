#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "inverse_interpolate.h"
#include "util.h"
#include "mathl.h"

static const int _dbgprint = 0;

/* TODO:  Consider Brent's method. */


/* Return x with v(x)=func(x,k) s.t. either of:
 *    1.  v(x) == n  and v(x-1-threshold) < n
 *    2.  v(x) < n   and v(x+1) > n
 */

#define LINEAR_INTERP(n, lo, hi, rlo, rhi) \
  (lo + (UV) (((double)(n-rlo) * (double)(hi-lo) / (double)(rhi-rlo))+0.5))

#define CALLBACK(n)  ((funck) ? funck(n,k) : func(n))

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
  UV mid, rlo, rhi, rmid, iloopc;

  if (hi != 0) {
    /* Given both lo and hi, halve the range on start. */
    mid = lo + ((hi-lo)>>1);
    rmid = CALLBACK(mid);
    if(_dbgprint)printf("  01 lo %lu  mid %lu  hi %lu\n", lo, mid, hi);
    if (rmid >= n) {
      hi = mid;  rhi = rmid;
      rlo = CALLBACK(lo);
      if (rlo == n)  RETURNI(lo);  /* Possible bad limit */
    } else {
      lo = mid;  rlo = rmid;
      rhi = CALLBACK(hi);
    }
  } else {
    /* They don't know what hi might be, so estimate something. */
    rlo = CALLBACK(lo);
    if (rlo == n)  RETURNI(lo);  /* Possible bad limit */
    while (hi == 0) {
      double estf = (double)n/(double)rlo - 0.004;
      if      (estf <= 1.004) estf = 1.004;
      else if (estf > 8.0)    estf = 8.0;
      mid =  ((double)UV_MAX/(double)lo <= estf)  ?  UV_MAX
          :  (UV) (estf * (double)lo + 1);
      if(_dbgprint)printf("  0s lo %lu  mid %lu  hi %lu\n", lo, mid, hi);
      rmid = CALLBACK(mid);
      if (rmid >= n) {  hi = mid;  rhi = rmid;  }
      else           {  lo = mid;  rlo = rmid;  }
      if (lo == UV_MAX)  break;  /* Overflow */
    }
  }

  MPUassert(rlo <= n && rhi >= n, "interpolation: bad initial limits");
  if ((hi-lo) <= 1)   RETURNI( (rlo == n || (rlo < n && rhi > n)) ? lo : hi );

  /* Step 1.  Linear interpolation until rhi is correct. */
  if(_dbgprint)printf("  1  lo %lu hi %lu\n", lo, hi);

  mid = (n == rhi)  ?  hi-1  :  LINEAR_INTERP(n,lo,hi,rlo,rhi);
  if (mid == lo) mid++;  else if (mid == hi) mid--;

  for (iloopc = 1;  (hi-lo) > 1 && rhi > n;  iloopc++) {
    MPUassert(lo < mid && mid < hi, "interpolation: assume 3 unique points");
    rmid = CALLBACK(mid);
    if (rmid >= n) { hi = mid;  rhi = rmid; }
    else           { lo = mid;  rlo = rmid; }
    if (rhi == n) break;
    mid += (IV)(((double)n-(double)rmid)*(double)(hi-lo) / (double)(rhi-rlo));
    /* Sometimes we get stuck getting closer and closer but not bracketing.
     * We could do Ridder's method of alternating bisection, or using a
     * multiplier on mid on alternate iterations to reflect about n.
     * What we're going to do instead is, every few loops, check if we're
     * very close to one of the edges and try to pull in the other edge.
     */
    if ((iloopc % 6) == 0) {
      UV close = .003*(hi-lo) + 1.0;
      if      (lo+close > mid) mid = lo+close;
      else if (hi-close < mid) mid = hi-close;
    }
    /* Alternately:
        if      (mid == lo) { mid = lo + .01*(hi-lo); }
        else if (mid == hi) { mid = hi - .01*(hi-lo); }
     */
    if (mid <= lo) mid=lo+1;  else if (mid >= hi) mid=hi-1;
    MPUassert(lo <= mid && mid <= hi, "interpolation: range error");
    if(_dbgprint)printf("  1s lo %lu  mid %lu  hi %lu  (%lu)\n", lo, mid, hi, rhi-n);
  }

  if (rlo == n)       RETURNI(lo);
  if ((hi-lo) <= 1)   RETURNI((rlo == n || (rlo < n && rhi > n)) ? lo : hi);

  MPUassert(rlo < n && rhi == n, "interpolation: bad step 1 interpolation");

  /* Step 2.  Ridder's method until we're very close. */

  MPUassert(rlo < n && rhi >= n, "interpolation: Ridder initial assumption");
  if(_dbgprint)printf("  2  lo %lu  mid %lu  hi %lu\n", lo, mid, hi);

  while ((hi-lo) > 8 && ((hi-lo) > threshold || rhi > n)) {
    UV x0 = lo,  x1 = lo + ((hi-lo)>>1);   /* x2 = hi */
    UV rx1 = CALLBACK(x1);
    IV fx0 = rlo-n,  fx1 = rx1-n,  fx2=rhi-n+1;

    double pos = ((double)(x1-x0) * (double)fx1)
               / sqrtl((double)fx1 * (double)fx1 - (double)fx0 * (double)fx2);
    UV x3 = x1 - (IV)(pos+0.5);

    if (x3 >= hi || x3 <= lo || x3 == x1) {
      /* We got nothing from the new point.  Just use the bisection. */
      if (rx1 >= n) { hi = x1; rhi = rx1; }
      else          { lo = x1; rlo = rx1; }
    } else {
      UV rx3 = CALLBACK(x3);
      /* Swap if needed to have:   [lo  x1  x3  hi]  */
      if (rx1 > rx3) { UV t=x1; x1=x3; x3=t;  t=rx1; fx1=rx3; rx3=t; }
      if      (rx1 >= n) {                      hi = x1; rhi = rx1; }
      else if (rx3 >= n) { lo = x1; rlo = rx1;  hi = x3; rhi = rx3; }
      else               { lo = x3; rlo = rx3; }
    }
    MPUassert(rlo < n && rhi >= n, "interpolation: Ridder step error");
    if(_dbgprint)printf("  2s lo %lu  mid %lu  hi %lu  (%lu)\n", lo, mid, hi, rhi-n);
  }

  /* Step 3.  Binary search. */

  /* Binary search until within threshold */
  while ((hi-lo) > 1 && ((hi-lo) > threshold || rhi > n)) {
    mid = lo + ((hi-lo)>>1);
    if (CALLBACK(mid) < n) lo = mid;   /* Keeps invariant f(lo) < n */
    else                   hi = mid;
  }
  if(_dbgprint)printf("final %lu - %lu threshold %lu\n", lo, hi, threshold);
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
  UV approx_nth_n, guess, gn, count, ming = 0, maxg = UV_MAX;

  approx_nth_n = guess = fnth(n);
  for (gn = 2; gn < 20; gn++) {
    IV adjust;
    MPUverbose(2, "  interp  %"UVuf"-th is around %"UVuf" ... ", n, guess);
    count = fcnt(guess);
    MPUverbose(2, "(%"IVdf")\n", (IV)(n-count));
    /* Stop guessing if within our tolerance */
    if (n==count || (n>count && n-count < tol) || (n<count && count-n < tol)) break;
    /* Determine how far off we think we are */
    adjust = (IV) (approx_nth_n - fnth(count));
    /* When computing new guess, ensure we don't overshoot.  Rarely used. */
    if (count <= n && guess > ming) ming = guess;   /* Previous guesses */
    if (count >= n && guess < maxg) maxg = guess;
    guess += adjust;
    if (guess <= ming || guess >= maxg) MPUverbose(2, "  fix min/max for %"UVuf"\n",n);
    if (guess <= ming) guess = ming + tol - 1;
    if (guess >= maxg) guess = maxg - tol + 1;
    /* TODO: if min/max dist is small, split the difference. */
  }
  if (gn == 20) count = fcnt(guess);

  if (fis) {
    if (count < n) {

      /* Increase count one at a time if needed */
      for ( ; count < n; count++)
        while (!fis(++guess))
          ;

    } else if (count >= n) {

      /* Make sure this is the least value at this count */
      while (!fis(guess))  guess--;
      /* Reduce count one at a time if needed */
      for ( ; count > n; count--)
        while (!fis(--guess))
          ;

    }
  }

  if (gcount) *gcount = count;
  return guess;
}
