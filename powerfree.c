#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FUNC_isqrt 1
#define FUNC_ipow 1
#define FUNC_ctz 1
#include "ptypes.h"
#include "constants.h"
#include "powerfree.h"
#include "util.h"
#include "factor.h"

static INLINE UV T(UV n) {
  return (n*(n+1)) >> 1;
}
static UV fprod(UV n, UV r) {
  UV P, fac[MPU_MAX_FACTORS+1];
  int i, nfactors;

  P = 1;
  nfactors = factor_exp(n, fac, 0);
  for (i = 0; i < nfactors; i++)
    P *= 1 - ipow(fac[i], r);
  return P;
}

int is_powerfree(UV n, uint32_t k)
{
  UV fac[MPU_MAX_FACTORS+1];
  UV exp[MPU_MAX_FACTORS+1];
  int i, nfactors;

  if (k < 2 || n <= 1)   return (n==1);

  if (k >= BITS_PER_WORD || n < (UVCONST(1) << (k-1)))
    return 1;

  /* Try to quickly find common powers so we don't have to factor */
  if (k == 2) {
    if (!(n %   4) || !(n %   9) || !(n %  25) || !(n %  49) ||
        !(n % 121) || !(n % 169) || !(n % 289) || !(n % 361))
      return 0;
    if (n < 529)
      return 1;
  } else if (k == 3) {
    if ( !(n % 8) || !(n % 27) || !(n % 125) )
      return 0;
  }

  /* A factor iterator would be good to use here */
  nfactors = factor_exp(n, fac, exp);
  for (i = 0; i < nfactors; i++) {
    if (exp[i] >= k)
      return 0;
  }

  return 1;
}

/* Basic method from https://arxiv.org/pdf/1107.4890.pdf */
static UV squarefree_count(UV n)
{
  signed char* mu;
  IV *M, *Mx, Mxisum;
  UV sqrtn, I, D, i, j, S1 = 0, S2 = 0;

  if (n < 4) return n;

  sqrtn = isqrt(n);
  I = rootint(n, 5);   /* times loglogn ^ (4/5) */
  D = isqrt(n / I);
  mu = range_moebius(0, D);

  S1 += n;
  for (i = 2; i <= D; i++)
    if (mu[i] != 0)
      S1 += mu[i] * (n/(i*i));

  New(0, M, D+1, IV);
  M[0] = 0;
  for (i = 1; i <= D; i++)
    M[i] = M[i-1] + mu[i];
  Safefree(mu);

  Newz(0, Mx, I+1, IV);
  Mxisum = 0;
  for (i = I-1; i > 0; i--) {
    IV Mxi = 1;
    UV xi = isqrt(n/i);
    UV L = isqrt(xi);
    for (j = 1; j <= xi/(L+1); j++)
      Mxi -= M[j] * (xi/j - xi/(j+1));
    for (j = 2; j <= L; j++)
      Mxi -=  (xi/j <= D)  ?  M[xi/j]  :  Mx[j*j*i];
    Mx[i] = Mxi;
    Mxisum += Mxi;
  }
  S2 = Mxisum - (I - 1) * M[D];
  Safefree(Mx);
  Safefree(M);

  return S1 + S2;
}

UV powerfree_count(UV n, uint32_t k)
{
  UV i, nk, count;

  if (k < 2) return (n >= 1);
  if (n < 4) return n;
  if (k == 2) return squarefree_count(n);

  count = n;
  nk = rootint(n, k);

  if (nk <= 100) {
    for (i = 2; i <= nk; i++) {
      int m = moebius(i);
      if (m != 0)
        count += m * (n / ipow(i, k));
    }
  } else {
    signed char* mu = range_moebius(0, nk);
    for (i = 2; i <= nk; i++)
      if (mu[i] != 0)
        count += mu[i] * (n/ipow(i,k));
    Safefree(mu);
  }
  return count;
}

UV powerfree_sum(UV n, uint32_t k)
{
  UV i, nk, sum;

  if (k < 2) return (n >= 1);

  if (n >= (UVCONST(1) << (BITS_PER_WORD/2)))  return 0;  /* Overflow */

  sum = T(n);
  nk = rootint(n, k);

  for (i = 2; i <= nk; i++) {
    int m = moebius(i);
    if (m != 0) {
      UV ik = (k==2) ? i*i : ipow(i,k);
      UV nik = n / ik;
      sum += m * ik * T(nik);
    }
  }
  return sum;
}


UV powerfree_part(UV n, uint32_t k)
{
  UV fac[MPU_MAX_FACTORS+1];
  UV exp[MPU_MAX_FACTORS+1];
  UV t, P;
  int i, nfactors;

  if (k < 2 || n <= 1)
    return (n==1);

  if (k >= BITS_PER_WORD || n < (UVCONST(1) << (k-1)))
    return n;

  /* Pull all powers of two out */
  t = ctz(n);
  P = n >> t;
  if ((t % k))  P <<= (t % k);

  nfactors = factor_exp(P, fac, exp);
  for (i = 0; i < nfactors; i++)
    if (exp[i] >= k)
      P /= ipow(fac[i], exp[i] - (exp[i] % k));

  return P;
}


UV powerfree_part_sum(UV n, uint32_t k)
{
  UV j, nk, sum = 0;

  if (k < 2 || n <= 1) return (n >= 1);

  if (n >= (UVCONST(1) << (BITS_PER_WORD/2)))  return 0;  /* Overflow */

  sum = T(n);
  nk = rootint(n,k);

  /* Using the factor iterator is overkill because of the limited range. */

  if (nk <= 100) {
    for (j = 2; j <= nk; j++)
      sum += fprod(j,k) * T(n/ipow(j,k));
  } else {
    UV P, *factors;
    factor_range_context_t fctx;
    int i, nfactors;

    fctx = factor_range_init(2, nk, 0);
    for (j = 2; j <= nk; j++) {
      nfactors = factor_range_next(&fctx);
      factors = fctx.factors;
      for (P = 1, i = 0; i < nfactors; i++)
        if (i == 0 || factors[i] != factors[i-1])
          P *= 1 - ipow(factors[i], k);
      sum += P * T(n/ipow(j,k));
    }
    factor_range_destroy(&fctx);
  }

  return sum;
}
