#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

UV powerfree_count(UV n, uint32_t k)
{
  UV i, nk, count;

  if (k < 2) return (n >= 1);

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
        count += mu[i] * ( (k==2) ? n/(i*i) : n/ipow(i,k) );
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
