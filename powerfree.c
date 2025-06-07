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
#include "real.h"

static INLINE UV T(UV n) {
  return (n+1)/2 * (n|1);
}
static UV fprod(UV n, UV r) {
  factored_t nf;
  UV P;
  uint32_t i;

  nf = factorint(n);
  for (P = 1, i = 0; i < nf.nfactors; i++)
    P *= 1 - ipow(nf.f[i], r);
  return P;
}

bool is_powerfree(UV n, uint32_t k)
{
  factored_t nf;
  uint32_t i;

  if (k < 2 || n <= 1)           return (n==1);

  if (k >= BITS_PER_WORD)        return 1;
  if (n < (UVCONST(1) << (k-1))) return 1;
  if (n == ((n >> k) << k))      return 0;
  if (k == 2)                    return is_square_free(n);

  /* Try to quickly find common powers so we don't have to factor */
  if (k == 3) {
    if ( !(n % 27) || !(n % 125) || !(n % 343) || !(n%1331) || !(n%2197) )
      return 0;
    if (n < 4913) return 1;
  }

  /* A factor iterator would be good to use here */
  nf = factorint(n);
  for (i = 0; i < nf.nfactors; i++) {
    if (nf.e[i] >= k)
      return 0;
  }

  return 1;
}

/* Basic method from https://arxiv.org/pdf/1107.4890.pdf */
static UV squarefree_count(UV n)
{
  signed char* mu;
  IV *M, *Mx, Mxisum, mert;
  UV I, D, i, j, S1 = 0, S2 = 0;

  if (n < 4) return n;

  I = rootint(n, 5);   /* times loglogn ^ (4/5) */
  D = isqrt(n / I);
  mu = range_moebius(0, D);

  S1 += n;
  New(0, M, D+1, IV);
  M[0] = 0;
  M[1] = 1;
  mert = 1;
  for (i = 2; i <= D; i++) {
    if (mu[i] != 0) {
      S1 += mu[i] * (n/(i*i));
      mert += mu[i];
    }
    M[i] = mert;
  }
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
  factored_t nf;
  UV t, P;
  uint32_t i;

  if (k < 2 || n <= 1)
    return (n==1);

  if (k >= BITS_PER_WORD || n < (UVCONST(1) << (k-1)))
    return n;

  /* Pull all powers of two out */
  t = ctz(n);
  P = n >> t;
  if ((t % k))  P <<= (t % k);

  nf = factorint(P);
  for (i = 0; i < nf.nfactors; i++)
    if (nf.e[i] >= k)
      P /= ipow(nf.f[i], nf.e[i] - (nf.e[i] % k));

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

#if BITS_PER_WORD == 64
  #define MAX_PFC2 UVCONST(11214275663373200251)
  #define MAX_PFC3 UVCONST(15345982395028449439)
  #define MAX_PFC4 UVCONST(17043655258566511333)
#else
  #define MAX_PFC2 UVCONST(2611027094)
  #define MAX_PFC3 UVCONST(3573014938)
  #define MAX_PFC4 UVCONST(3968285222)
#endif

UV nth_powerfree(UV n, uint32_t k)
{
  long double zm;
  UV qk, count, diff, thresh, i;

  if (k < 2) return 0;
  if (n < 4) return n;

  /* Check for overflow. */
  if (k == 2 && n > MAX_PFC2) return 0;
  if (k == 3 && n > MAX_PFC3) return 0;
  if (k >= 4 && n > MAX_PFC4) {
    if (k == 4) return 0;
    if (n > powerfree_count(UV_MAX,k)) return 0;
  }

  /* Step 1:  Density ZM and expected value QK. */
  zm = 1.0 + ld_riemann_zeta(k);
  qk = (UV)(zm * (long double) n + 0.5);
  thresh = (k <= 2) ? 200 : (k == 3) ? 60 : (k == 4) ? 2 : 1;

  for (i = 0; i < 10; i++) {
    /* Step 2: Initial count at QK and difference from goal. */
    count = powerfree_count(qk, k);
    diff = (count >= n) ? count-n : n-count;
    /* Step 3: Update estimate using expected density. */
    if (diff <= thresh) break;
    if (count > n)  qk -= (UV)((long double)diff * zm);
    else            qk += (UV)((long double)diff * zm);
  }

  /* Step 4: Get ourselves onto a powerfree number */
  while (!is_powerfree(qk,k))  qk--;

  /* Step 5: Walk forwards or backwards until we get to the goal. */
  while (count != n) {
    do { qk += (count < n) ? 1 : -1; } while (!is_powerfree(qk,k));
    count += (count < n) ? 1 : -1;
  }
  return qk;
}

/******************************************************************************/

UV squarefree_kernel(UV n)
{
  factored_t nf;
  UV P;
  uint32_t i;

  nf = factorint(n);
  for (P = 1, i = 0; i < nf.nfactors; i++)
    P *= nf.f[i];
  return P;
}
