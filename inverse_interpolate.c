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
 *    1.  v(x) == n  and v(x-1-threshold) < n  [within supplied bounds]
 *    2.  v(x) < n   and v(x+1) > n
 */

#define MPU_CALLBACK(n)  ((funck) ? funck(n,k) : func(n))

#if 0   /* Debugging return, checking the conditions above. */
#define RETURNI(x) \
  { \
    UV v = x; \
    UV rv = MPU_CALLBACK(v); \
    /* printf("v  %"UVuf"    rv %"UVuf"   n %"UVuf"\n",v,rv,n); */\
    MPUassert( rv <= n, "BAD INTERP  v > n" ); \
    if (rv == n) { \
      /* Approximation callbacks may not be monotonic, so we must only    */ \
      /* check inside the caller-supplied bounds.                         */ \
      /* We lose the ability to detect lower bounds that were too high.   */ \
      if (v > bound_lo && v-bound_lo > threshold) { \
        /* printf("threshold %"UVuf" v %"UVuf"    func(%"UVuf") = %"UVuf"\n", threshold, v, v-1-threshold, MPU_CALLBACK(v-1-threshold)); */\
        MPUassert( MPU_CALLBACK(v-1-threshold) < n, "BAD INTERP  v-1-thresh >= n" ); \
      } \
    } else { \
      MPUassert( MPU_CALLBACK(v+1) > n, "BAD INTERP  v+1 <= n" ); \
    } \
    return v; \
  }
#else
  #define RETURNI(x) { return x; }
#endif

/* Return a point strictly inside [lo,hi], interpolating the endpoint values
 * [rlo,rhi].  Subtract before converting to double: large UVs that differ by
 * only a few units can become the same double. */
static UV _linear_interpolate(UV lo, UV hi, UV rlo, UV rhi, UV n)
{
  UV rdiff, rspan, xspan, move;
  double dmove;

  MPUassert(hi-lo > 1, "interpolation: no interior point");
  MPUassert(rlo <= n && n <= rhi, "interpolation: target outside range");

  if (n == rlo) return lo+1;
  if (n == rhi) return hi-1;

  rdiff = n-rlo;
  rspan = rhi-rlo;
  xspan = hi-lo;
  dmove = (double)rdiff * (double)xspan / (double)rspan + 0.5;
  move = (dmove >= (double)xspan) ? xspan : (UV)dmove;

  if (move == 0)     return lo+1;
  if (move >= xspan) return hi-1;
  return lo+move;
}

static UV _inverse_interpolate(UV bound_lo, UV bound_hi, UV n,
                               UV k, UV (*funck)(UV mid, UV k),
                               UV (*func)(UV mid),
                               UV threshold) {
  UV mid, rlo, rhi, rmid, iloopc, lo = bound_lo, hi = bound_hi;

  if (hi != 0) {
    /* Given both lo and hi, halve the range on start. */
    mid = lo + ((hi-lo)>>1);
    rmid = MPU_CALLBACK(mid);
    if(_dbgprint)printf("  01 lo %"UVuf"  mid %"UVuf"  hi %"UVuf"\n", lo, mid, hi);
    if (rmid >= n) {
      hi = mid;  rhi = rmid;
      rlo = MPU_CALLBACK(lo);
      if (rlo == n)  RETURNI(lo);  /* Possible bad limit */
    } else {
      lo = mid;  rlo = rmid;
      rhi = MPU_CALLBACK(hi);
    }
  } else {
    /* They don't know what hi might be, so estimate something. */
    rlo = MPU_CALLBACK(lo);
    if (rlo == n)  RETURNI(lo);  /* Possible bad limit */
    for (;;) {
      double estf = (double)n/(double)rlo - 0.004;
      if      (estf <= 1.004) estf = 1.004;
      else if (estf > 8.0)    estf = 8.0;
      mid =  ((double)UV_MAX/(double)lo <= estf)  ?  UV_MAX
          :  (UV) (estf * (double)lo + 1);
      if(_dbgprint)printf("  0s lo %"UVuf"  mid %"UVuf"  hi %"UVuf"\n", lo, mid, hi);
      rmid = MPU_CALLBACK(mid);
      if (rmid >= n) {
        hi = mid;
        rhi = rmid;
        break;
      }
      lo = mid;
      rlo = rmid;
      if (lo == UV_MAX) return 0;  /* No representable upper bracket */
    }
  }

  MPUassert(rlo <= n && rhi >= n, "interpolation: bad initial limits");
  if ((hi-lo) <= 1)   RETURNI( (rlo == n || (rlo < n && rhi > n)) ? lo : hi );

  /* Step 1.  Linear interpolation until rhi is correct. */
  if(_dbgprint)printf("  1  lo %"UVuf" hi %"UVuf"\n", lo, hi);

  mid = _linear_interpolate(lo, hi, rlo, rhi, n);

  for (iloopc = 1;  (hi-lo) > 1 && rhi > n;  iloopc++) {
    MPUassert(lo < mid && mid < hi, "interpolation: assume 3 unique points");
    rmid = MPU_CALLBACK(mid);
    if (rmid >= n) { hi = mid;  rhi = rmid; }
    else           { lo = mid;  rlo = rmid; }
    if (rhi == n || (hi-lo) <= 1) break;
    mid = _linear_interpolate(lo, hi, rlo, rhi, n);
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
    if(_dbgprint)printf("  1s lo %"UVuf"  mid %"UVuf"  hi %"UVuf"  (%"UVuf")\n", lo, mid, hi, rhi-n);
  }

  if (rlo == n)       RETURNI(lo);
  if ((hi-lo) <= 1)   RETURNI((rlo == n || (rlo < n && rhi > n)) ? lo : hi);

  MPUassert(rlo < n && rhi == n, "interpolation: bad step 1 interpolation");

  /* Step 2.  Ridder's method until we're very close. */

  MPUassert(rlo < n && rhi >= n, "interpolation: Ridder initial assumption");
  if(_dbgprint)printf("  2  lo %"UVuf"  mid %"UVuf"  hi %"UVuf"\n", lo, mid, hi);

  while ((hi-lo) > 8 && ((hi-lo) > threshold || rhi > n)) {
    UV x0 = lo,  x1 = lo + ((hi-lo)>>1);   /* x2 = hi */
    UV rx1 = MPU_CALLBACK(x1);
    double d0 = (double) (n-rlo);     /* > 0 */
    double d1 = (double) (rx1 >= n ? rx1-n : n-rx1);
    double d2 = (double) (rhi-n+1);   /* > 0 */
    double den = d1*d1 + d0*d2;       /* f1^2 - f0*f2 */
    double sgn = (rx1 >= n) ? 1.0L : -1.0L;
    double pos = sgn * ((double)(x1 - x0) * d1) / sqrtl(den);
    UV x3 = x1 - (IV)(pos + (pos >= 0 ? 0.5L : -0.5L));

    if(_dbgprint)printf("  2s lo %"UVuf"  mid %"UVuf"  hi %"UVuf"  (%"UVuf")\n", lo, x1, hi, (rx1>n) ? rx1-n : n-rx1);

    if (x3 >= hi || x3 <= lo || x3 == x1) {
      /* We got nothing from the new point.  Just use the bisection. */
      if (rx1 >= n) { hi = x1; rhi = rx1; }
      else          { lo = x1; rlo = rx1; }
    } else {
      UV rx3 = MPU_CALLBACK(x3);
      if(_dbgprint)printf("  2S lo %"UVuf"  mid %"UVuf"  hi %"UVuf"  (%"UVuf")\n", lo, x3, hi, (rx3>n) ? rx3-n : n-rx3);
      /* Swap if needed to have:   [lo  x1  x3  hi]  */
      if (rx1 > rx3) { UV t=x1; x1=x3; x3=t;  t=rx1; rx1=rx3; rx3=t; }
      if      (rx1 >= n) {                      hi = x1; rhi = rx1; }
      else if (rx3 >= n) { lo = x1; rlo = rx1;  hi = x3; rhi = rx3; }
      else               { lo = x3; rlo = rx3; }
    }
    MPUassert(rlo < n && rhi >= n, "interpolation: Ridder step error");
  }

  /* Step 3.  Binary search. */

  /* Binary search until within threshold */
  while ((hi-lo) > 1 && ((hi-lo) > threshold || rhi > n)) {
    mid = lo + ((hi-lo)>>1);
    if (MPU_CALLBACK(mid) < n) lo = mid;   /* Keeps invariant f(lo) < n */
    else                   hi = mid;
  }
  if(_dbgprint)printf("final %"UVuf" - %"UVuf" threshold %"UVuf"\n", lo, hi, threshold);
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
                           bool (*fis)(UV n)   /* optional */
                          ) {
  UV approx_nth_n, guess, gn, count, ming = 0, maxg = UV_MAX;

  MPUassert(tol > 0, "interpolate_with_approx: tolerance must be positive");

  approx_nth_n = guess = fnth(n);
  for (gn = 2; gn < 20; gn++) {
    UV adjust, approx_nth_count;
    MPUverbose(2, "  interp  %"UVuf"-th is around %"UVuf" ... ", n, guess);
    count = fcnt(guess);
    if (n >= count) { MPUverbose(2, "(%"UVuf")\n", n-count); }
    else            { MPUverbose(2, "(-%"UVuf")\n", count-n); }
    /* Stop guessing if within our tolerance */
    if (n==count || (n>count && n-count < tol) || (n<count && count-n < tol)) break;
    /* Determine how far off we think we are */
    approx_nth_count = fnth(count);
    /* When computing new guess, ensure we don't overshoot.  Rarely used. */
    if (count <= n && guess > ming) ming = guess;   /* Previous guesses */
    if (count >= n && guess < maxg) maxg = guess;
    if (approx_nth_n >= approx_nth_count) {
      adjust = approx_nth_n - approx_nth_count;
      guess = (adjust > UV_MAX-guess) ? UV_MAX : guess+adjust;
    } else {
      adjust = approx_nth_count - approx_nth_n;
      guess = (adjust > guess) ? 0 : guess-adjust;
    }
    /* Pull an outlying estimate back inside the known bounds when possible. */
    if (guess <= ming || guess >= maxg) {
      UV span = maxg-ming;
      UV step = (tol > 1) ? tol-1 : 1;
      bool below = (guess <= ming);
      MPUverbose(2, "  fix min/max for %"UVuf"\n",n);
      if (span <= 1)      guess = below ? maxg : ming;
      else if (step >= span) guess = ming + (span >> 1);
      else                guess = below ? ming+step : maxg-step;
    }
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
