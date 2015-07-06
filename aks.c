#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* The AKS primality algorithm for native integers.
 *
 * There are two versions here.  The v6 algorithm from the latest AKS paper,
 * as well as one with improvements from Bernstein and Voloch and better r/s
 * selection derived from Folkmar Bornemann's 2002 Pari implementation.
 *
 * Note that AKS is very, very slow compared to other methods.  It is, however,
 * polynomial in log(N), and log-log performance graphs show nice straight
 * lines for both implementations.  However APR-CL and ECPP both start out
 * much faster and the slope will be less for any sizes of N that we're
 * interested in.
 *
 * For native 64-bit integers this is purely a coding exercise, as BPSW is
 * a million times faster and gives proven results.
 *
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
 * Copyright 2012-2014, Dana Jacobsen.
 */

#define SQRTN_SHORTCUT 1

/* Use improvements from Bornemann's 2002 implementation */
#define IMPL_BORNEMANN 1

#include "ptypes.h"
#include "aks.h"
#define FUNC_isqrt 1
#define FUNC_gcd_ui 1
#include "util.h"
#include "cache.h"
#include "mulmod.h"
#include "factor.h"

#if IMPL_BORNEMANN
static int is_primitive_root(UV n, UV r)
{
  UV fac[MPU_MAX_FACTORS+1];
  int i, nfacs;
  /* if ( (r&1) & powmod(n, (r-1)>>1, r) == 1 )  return 0; */
  nfacs = factor_exp(r-1, fac, 0);
  for (i = 0; i < nfacs; i++) {
    UV m = powmod(n, (r-1)/fac[i], r);
    if (m == 1) return 0;
  }
  return (gcd_ui(n,r) == 1);
}
/* We could use lgamma, but it isn't in MSVC and not in pre-C99.  The only
 * sure way to find if it is available is test compilation (ala autoconf).
 * Instead, we'll just use our own implementation.
 * See http://mrob.com/pub/ries/lanczos-gamma.html for alternates. */
static double lanczos_coef[8+1] =
{ 0.99999999999980993, 676.5203681218851, -1259.1392167224028,
  771.32342877765313, -176.61502916214059, 12.507343278686905,
  -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 };
static double log_sqrt_two_pi =  0.91893853320467274178;
static double log_gamma(double x)
{
  double base = x + 7 + 0.5;
  double sum = 0;
  int i;
  for (i = 8; i >= 1; i--)
    sum += lanczos_coef[i] / (x + (double)i);
  sum += lanczos_coef[0];
  sum = log_sqrt_two_pi + log(sum/x) + ( (x+0.5)*log(base) - base );
  return sum;
}
#undef lgamma
#define lgamma(x) log_gamma(x)
#endif

#if 0
/* Naive znorder.  Works well if limit is small.  Note arguments.  */
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

  /* Determine max degree of px and py */
  for (degpx = r-1; degpx > 0 && !px[degpx]; degpx--) ; /* */
  for (degpy = r-1; degpy > 0 && !py[degpy]; degpy--) ; /* */
  /* We can sum at least j values at once */
  j = (mod >= HALF_WORD) ? 0 : (UV_MAX / ((mod-1)*(mod-1)));

  if (j >= degpx || j >= degpy) {
    /* res will be written completely, so no need to set */
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
    memset(res, 0, r * sizeof(UV));  /* Zero result accumulator */
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
  int native_sqr = (mod > isqrt(UV_MAX/(2*r))) ? 0 : 1;

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
    if (native_sqr) {
      while (pp1 < ppend)
        sum += 2 * *pp1++  *  *pp2--;
      /* Special treatment for last point */
      c = px[s_end];
      sum += (s_end*2 == d)  ?  c*c  :  2*c*px[d-s_end];
      rindex = (d < r) ? d : d-r;  /* d % r */
      res[rindex] = (res[rindex] + sum) % mod;
#if HAVE_UINT128
    } else {
      uint128_t max = ((uint128_t)1 << 127) - 1;
      uint128_t c128, sum128 = 0;

      while (pp1 < ppend) {
        c128 = ((uint128_t)*pp1++)  *  ((uint128_t)*pp2--);
        if (c128 > max) c128 %= mod;
        c128 <<= 1;
        if (c128 > max) c128 %= mod;
        sum128 += c128;
        if (sum128 > max) sum128 %= mod;
      }
      c128 = px[s_end];
      if (s_end*2 == d) {
        c128 *= c128;
      } else {
        c128 *= px[d-s_end];
        if (c128 > max) c128 %= mod;
        c128 <<= 1;
      }
      if (c128 > max) c128 %= mod;
      sum128 += c128;
      if (sum128 > max) sum128 %= mod;
      rindex = (d < r) ? d : d-r;  /* d % r */
      res[rindex] = ((uint128_t)res[rindex] + sum128) % mod;
#else
    } else {
      while (pp1 < ppend) {
        UV p1 = *pp1++;
        UV p2 = *pp2--;
        sum = addmod(sum, mulmod(2, mulmod(p1, p2, mod), mod), mod);
      }
      c = px[s_end];
      if (s_end*2 == d)
        sum = addmod(sum, sqrmod(c, mod), mod);
      else
        sum = addmod(sum, mulmod(2, mulmod(c, px[d-s_end], mod), mod), mod);
      rindex = (d < r) ? d : d-r;  /* d % r */
      res[rindex] = addmod(res[rindex], sum, mod);
#endif
    }
  }
  memcpy(px, res, r * sizeof(UV)); /* put result in px */
}

static UV* poly_mod_pow(UV* pn, UV power, UV r, UV mod)
{
  UV *res, *temp;

  Newz(0, res, r, UV);
  New(0, temp, r, UV);
  res[0] = 1;

  while (power) {
    if (power & 1)  poly_mod_mul(res, pn, temp, r, mod);
    power >>= 1;
    if (power)      poly_mod_sqr(pn, temp, r, mod);
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

/*
 * Avanzi and Mihǎilescu, 2007
 * http://www.uni-math.gwdg.de/preda/mihailescu-papers/ouraks3.pdf
 * "As a consequence, one cannot expect the present variants of AKS to
 *  compete with the earlier primality proving methods like ECPP and
 *  cyclotomy." - conclusion regarding memory consumption
 */
int is_aks_prime(UV n)
{
  UV r, s, a;
  int verbose;

  if (n < 2)
    return 0;
  if (n == 2)
    return 1;

  if (is_power(n, 0))
    return 0;

  if (n > 11 && ( !(n%2) || !(n%3) || !(n%5) || !(n%7) || !(n%11) )) return 0;
  /* if (!is_prob_prime(n)) return 0; */

  verbose = _XS_get_verbose();
#if IMPL_BORNEMANN == 0
  {
    UV sqrtn = isqrt(n);
    double log2n = log(n) / log(2);   /* C99 has a log2() function */
    UV limit = (UV) floor(log2n * log2n);

    if (verbose) { printf("# aks limit is %lu\n", (unsigned long) limit); }

    for (r = 2; r < n; r++) {
      if ((n % r) == 0)
        return 0;
#if SQRTN_SHORTCUT
      if (r > sqrtn)
        return 1;
#endif
      if (znorder(n, r) > limit)
        break;
    }

    if (r >= n)
      return 1;

    s = (UV) floor(sqrt(r-1) * log2n);
  }
#else
  {
    UV fac[MPU_MAX_FACTORS+1];
    UV slim;
    double c1, c2, x;
    double const t = 48;
    double const t1 = (1.0/((t+1)*log(t+1)-t*log(t)));
    double const dlogn = log(n);
    r = next_prime( (UV) (t1*t1 * dlogn*dlogn) );
    while (!is_primitive_root(n,r))
      r = next_prime(r);

    slim = (UV) (2*t*(r-1));
    c1 = lgamma(r-1);
    c2 = dlogn * floor(sqrt(r));
    { /* Binary search for first s in [1,slim] where x >= 0 */
      UV i = 1;
      UV j = slim;
      while (i < j) {
        s = i + (j-i)/2;
        x = (lgamma(r-1+s) - c1 - lgamma(s+1)) / c2 - 1.0;
        if (x < 0)  i = s+1;
        else        j = s;
      }
      s = i-1;
    }
    s = (s+3) >> 1;
    /* Bornemann checks factors up to (s-1)^2, we check to max(r,s) */
    /* slim = (s-1)*(s-1); */
    slim = (r > s) ? r : s;
    if (verbose > 1) printf("# aks trial to %lu\n", slim);
    if (trial_factor(n, fac, slim) > 1)
      return 0;
    if (slim >= HALF_WORD || (slim*slim) >= n)
      return 1;
  }
#endif

  if (verbose) { printf("# aks r = %lu  s = %lu\n", (unsigned long) r, (unsigned long) s); }

  /* Almost every composite will get recognized by the first test.
   * However, we need to run 's' tests to have the result proven for all n
   * based on the theorems we have available at this time. */
  for (a = 1; a <= s; a++) {
    if (! test_anr(a, n, r) )
      return 0;
    if (verbose>1) { printf("."); fflush(stdout); }
  }
  if (verbose>1) { printf("\n"); }
  return 1;
}
