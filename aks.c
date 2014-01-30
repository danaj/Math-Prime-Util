#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/*
 * The AKS v6 algorithm, for native integers.  Based on the AKS v6 paper.
 * As with most AKS implementations, it's really slow.
 *
 * When n < 2^(wordbits/2)-1, we can do a straightforward intermediate:
 *      r = (r + a * b) % n
 * If n is larger, then these are replaced with:
 *      r = addmod( r, mulmod(a, b, n), n)
 * which is a lot more work, but keeps us correct.
 *
 * Software that does polynomial convolutions followed by a modulo can be
 * very fast, but will fail when n >= (2^wordbits)/r.
 *
 * This is all much easier in GMP.
 *
 * Copyright 2012, Dana Jacobsen.
 */

#define SQRTN_SHORTCUT 1

#include "ptypes.h"
#include "aks.h"
#define FUNC_isqrt 1
#include "util.h"
#include "cache.h"
#include "mulmod.h"

/* Naive znorder.  Works well here because limit will be very small. */
static UV order(UV r, UV n, UV limit) {
  UV j;
  UV t = 1;
  for (j = 1; j <= limit; j++) {
    t = mulmod(t, n, r);
    if (t == 1)
      break;
  }
  return j;
}

#if 0
static void poly_print(UV* poly, UV r)
{
  int i;
  for (i = r-1; i >= 1; i--) {
    if (poly[i] != 0)
      printf("%lux^%d + ", poly[i], i);
  }
  if (poly[0] != 0) printf("%lu", poly[0]);
  printf("\n");
}
#endif

static void poly_mod_mul(UV* px, UV* py, UV* res, UV r, UV mod)
{
  UV degpx, degpy;
  UV i, j, pxi, pyj, rindex;

  memset(res, 0, r * sizeof(UV));

  /* Determine max degree of px and py */
  for (degpx = r-1; degpx > 0 && !px[degpx]; degpx--) ; /* */
  for (degpy = r-1; degpy > 0 && !py[degpy]; degpy--) ; /* */
  /* We can sum at least j values at once */
  j = (mod >= HALF_WORD) ? 0 : (UV_MAX / ((mod-1)*(mod-1)));

  if (j >= degpx || j >= degpy) {
    for (rindex = 0; rindex < r; rindex++) {
      UV sum = 0;
      j = rindex;
      for (i = 0; i <= degpx; i++) {
        if (j <= degpy)
          sum += px[i] * py[j];
        j = (j == 0) ? r-1 : j-1;
      }
      res[rindex] = sum % mod;
    }
  } else {
    for (i = 0; i <= degpx; i++) {
      pxi = px[i];
      if (pxi == 0)  continue;
      if (mod < HALF_WORD) {
        for (j = 0; j <= degpy; j++) {
          pyj = py[j];
          rindex = i+j;   if (rindex >= r)  rindex -= r;
          res[rindex] = (res[rindex] + (pxi*pyj) ) % mod;
        }
      } else {
        for (j = 0; j <= degpy; j++) {
          pyj = py[j];
          rindex = i+j;   if (rindex >= r)  rindex -= r;
          res[rindex] = muladdmod(pxi, pyj, res[rindex], mod);
        }
      }
    }
  }
  memcpy(px, res, r * sizeof(UV)); /* put result in px */
}
static void poly_mod_sqr(UV* px, UV* res, UV r, UV mod)
{
  UV c, d, s, sum, rindex, maxpx;
  UV degree = r-1;

  memset(res, 0, r * sizeof(UV)); /* zero out sums */
  /* Discover index of last non-zero value in px */
  for (s = degree; s > 0; s--)
    if (px[s] != 0)
      break;
  maxpx = s;
  /* 1D convolution */
  for (d = 0; d <= 2*degree; d++) {
    UV *pp1, *pp2, *ppend;
    UV s_beg = (d <= degree) ? 0 : d-degree;
    UV s_end = ((d/2) <= maxpx) ? d/2 : maxpx;
    if (s_end < s_beg) continue;
    sum = 0;
    pp1 = px + s_beg;
    pp2 = px + d - s_beg;
    ppend = px + s_end;
    while (pp1 < ppend)
      sum += 2 * *pp1++  *  *pp2--;
    /* Special treatment for last point */
    c = px[s_end];
    sum += (s_end*2 == d)  ?  c*c  :  2*c*px[d-s_end];
    rindex = (d < r) ? d : d-r;  /* d % r */
    res[rindex] = (res[rindex] + sum) % mod;
  }
  memcpy(px, res, r * sizeof(UV)); /* put result in px */
}

static UV* poly_mod_pow(UV* pn, UV power, UV r, UV mod)
{
  UV* res;
  UV* temp;
  int use_sqr = (mod > isqrt(UV_MAX/r)) ? 0 : 1;

  Newz(0, res, r, UV);
  New(0, temp, r, UV);
  if ( (res == 0) || (temp == 0) )
    croak("Couldn't allocate space for polynomial of degree %lu\n", (unsigned long) r);

  res[0] = 1;

  while (power) {
    if (power & 1)  poly_mod_mul(res, pn, temp, r, mod);
    power >>= 1;
    if (power) {
      if (use_sqr)  poly_mod_sqr(pn, temp, r, mod);
      else          poly_mod_mul(pn, pn, temp, r, mod);
    }
  }
  Safefree(temp);
  return res;
}

static int test_anr(UV a, UV n, UV r)
{
  UV* pn;
  UV* res;
  UV i;
  int retval = 1;

  Newz(0, pn, r, UV);
  if (pn == 0)
    croak("Couldn't allocate space for polynomial of degree %lu\n", (unsigned long) r);
  a %= r;
  pn[0] = a;
  pn[1] = 1;
  res = poly_mod_pow(pn, n, r, n);
  res[n % r] = addmod(res[n % r], n - 1, n);
  res[0]     = addmod(res[0],     n - a, n);

  for (i = 0; i < r; i++)
    if (res[i] != 0)
      retval = 0;
  Safefree(res);
  Safefree(pn);
  return retval;
}

int _XS_is_aks_prime(UV n)
{
  UV sqrtn, limit, r, rlimit, a;
  double log2n;
  int verbose;

  if (n < 2)
    return 0;
  if (n == 2)
    return 1;

  if (is_power(n, 0))
    return 0;

  sqrtn = isqrt(n);
  log2n = log(n) / log(2);   /* C99 has a log2() function */
  limit = (UV) floor(log2n * log2n);

  verbose = _XS_get_verbose();
  if (verbose) { printf("# aks limit is %lu\n", (unsigned long) limit); }

  for (r = 2; r < n; r++) {
    if ((n % r) == 0)
      return 0;
#if SQRTN_SHORTCUT
    if (r > sqrtn)
      return 1;
#endif
    if (order(r, n, limit) > limit)
      break;
  }

  if (r >= n)
    return 1;

  rlimit = (UV) floor(sqrt(r-1) * log2n);

  if (verbose) { printf("# aks r = %lu  rlimit = %lu\n", (unsigned long) r, (unsigned long) rlimit); }

  for (a = 1; a <= rlimit; a++) {
    if (! test_anr(a, n, r) )
      return 0;
    if (verbose>1) { printf("."); fflush(stdout); }
  }
  if (verbose>1) { printf("\n"); }
  return 1;
}
