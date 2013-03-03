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
#include "util.h"
#include "sieve.h"
#include "factor.h"
#include "cache.h"
#include "mulmod.h"

static UV log2floor(UV n) {
  UV log2n = 0;
  while (n >>= 1)
    log2n++;
  return log2n;
}

/* Bach and Sorenson (1993) would be better */
static int is_perfect_power(UV n) {
  UV b, last;
  if ((n <= 3) || (n == UV_MAX)) return 0;
  if ((n & (n-1)) == 0)  return 1;          /* powers of 2    */
  last = log2floor(n-1) + 1;
#if (BITS_PER_WORD == 32) || (DBL_DIG > 19)
  if (1) {
#elif DBL_DIG == 10
  if (n < UVCONST(10000000000)) {
#elif DBL_DIG == 15
  if (n < UVCONST(1000000000000000)) {
#else
  if ( n < (UV) pow(10, DBL_DIG) ) {
#endif
    /* Simple floating point method.  Fast, but need enough mantissa. */
    b = isqrt(n); if (b*b == n)  return 1; /* perfect square */
    for (b = 3; b < last; b = _XS_next_prime(b)) {
      UV root = pow(n, 1.0 / (double)b) + 0.5;
      if ( ((UV)(pow(root, b)+0.5)) == n)  return 1;
    }
  } else {
    /* Dietzfelbinger, algorithm 2.3.5 (without optimized exponential) */
    for (b = 2; b <= last; b++) {
      UV a = 1;
      UV c = n;
      while (c >= HALF_WORD) c = (1+c)>>1;
      while ((c-a) >= 2) {
        UV m, maxm, p, i;
        m = (a+c) >> 1;
        maxm = UV_MAX / m;
        p = m;
        for (i = 2; i <= b; i++) {
          if (p > maxm) { p = n+1; break; }
          p *= m;
        }
        if (p == n)  return 1;
        if (p < n)
          a = m;
        else
          c = m;
      }
    }
  }
  return 0;
}

static UV order(UV r, UV n, UV limit) {
  UV j;
  UV t = 1;
  for (j = 1; j <= limit; j++) {
    t = (t * n) % r;
    if (t == 1)
      break;
  }
  return j;
}

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

static void poly_mod_mul(UV* px, UV* py, UV* res, UV r, UV mod)
{
  UV i, j, pxi, pyj, rindex;

  memset(res, 0, r * sizeof(UV));
  for (i = 0; i < r; i++) {
    pxi = px[i];
    if (pxi == 0)  continue;
    for (j = 0; j < r; j++) {
      pyj = py[j];
      if (pyj == 0)  continue;
      rindex = (i+j) < r ? i+j : i+j-r; /* (i+j) % r */
      if (mod < HALF_WORD) {
        res[rindex] = (res[rindex] + (pxi*pyj) ) % mod;
      } else {
        res[rindex] = muladdmod(pxi, pyj, res[rindex], mod);
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
    UV s_beg = (d <= degree) ? 0 : d-degree;
    UV s_end = ((d/2) <= maxpx) ? d/2 : maxpx;
    if (s_end < s_beg) continue;
    sum = 0;
    for (s = s_beg; s < s_end; s++) {
      c = px[s];
      sum += 2*c * px[d-s];
    }
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
  int verbose = _XS_get_verbose();

  if (n < 2)
    return 0;
  if (n == 2)
    return 1;

  if (is_perfect_power(n))
    return 0;

  sqrtn = isqrt(n);
  log2n = log(n) / log(2);   /* C99 has a log2() function */
  limit = (UV) floor(log2n * log2n);

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
