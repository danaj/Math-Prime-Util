#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#define FUNC_isqrt 1
#define FUNC_ctz 1
#define FUNC_gcd_ui 1
#define FUNC_ipow 1
#include "util.h"
#include "sieve.h"
#include "cache.h"
#include "constants.h"
#include "inverse_interpolate.h"

int is_powerful(UV n, UV k) {
  UV pk;
  int res;

  if (n <= 1 || k <= 1) return 1;

  if (!(n&1)) { /* Check and remove all multiples of 2 */
    if (n & ((UVCONST(1) << k)-1)) return 0;
    n >>= ctz(n);
    if (n == 1) return 1;
  }

  if (k > MPU_MAX_POW3) return 0;
  /* if (k > logint(n,3)) return 0; */

  /* Quick checks */
  if (k == 2) {
    if (   (!(n %  3) && (n %     9))
        || (!(n %  5) && (n %    25))
        || (!(n %  7) && (n %    49))
        || (!(n % 11) && (n %   121))
        || (!(n % 13) && (n %   169)) ) return 0;
  } else if (k == 3) {
    if (   (!(n %  3) && (n %    27))
        || (!(n %  5) && (n %   125))
        || (!(n %  7) && (n %   343))
        || (!(n % 11) && (n %  1331)) ) return 0;
  } else {
    if (   (!(n %  3) && (n %    81))
        || (!(n %  5) && (n %   625))
        || (!(n %  7) && (n %  2401))
        || (!(n % 11) && (n % 14641)) ) return 0;
  }

#if 0   /* Full factoring.  Very simple and reasonably efficient. */
  {
    UV fac[MPU_MAX_FACTORS+1];
    UV exp[MPU_MAX_FACTORS+1];
    int i, nfactors;
    nfactors = factor_exp(n, fac, exp);
    for (i = 0; i < nfactors; i++)
      if (exp[i] < k)
        return 0;
    return 1;
  }
#endif

  /* Rather than full factoring, we'll use trial division.  For k=2, we
   * only need to check up to the fourth root of n, and k=3 to the sixth.
   * Even for k=2 this is faster than full factoring on average. */

  /* At every checkpoint (new prime p) for k=2 either:
   *   1) N < p^4  and N=1, p^2, q^2, p^3, or q^3. (q>p).  Return 1.
   *   2) N < p^4  otherwise, N cannot be powerful.  Return 0;
   *   3) N = p^4  is_square caught this and returned 1.
   * So the next possibility is p^2 * q^2 where q = next_prime(p).
   * Check n < p^4 before each new prime, and condition 1 after modifying n.
   */

  /* We could save 5-10% by special casing k=2 and pulling 3 from the loop */

  if (n == 1 || powerof(n) >= k) return 1;
  res = 0;
  START_DO_FOR_EACH_PRIME(3, rootint(n, 2*k)) {
    pk = ipow(p,k);
    if (n < pk*pk) break;
    if (!(n%p)) {
      if (n%pk) break;
      for (n /= pk; (n%p) == 0; n /= p) ;
      if (n == 1 || powerof(n) >= k) { res = 1; break; }
    }
  } END_DO_FOR_EACH_PRIME
  return res;
}


static unsigned char* _squarefree_range(UV lo, UV hi) {
  unsigned char* isf;
  UV i, i2, j, range = hi-lo+1, sqrthi = isqrt(hi);

  if (hi < lo) return 0;

  New(0, isf, range, unsigned char);
  memset(isf, 1, range);
  if (lo == 0) isf[0] = 0;
  for (i = 2; i <= sqrthi; i++)
    for (i2 = i*i, j = P_GT_LO(i2, i2, lo); j <= hi; j += i2)
      isf[j] = 0;
  return isf;
}

static UV _pcr(UV n, UV k, unsigned char* isf, UV m, UV r) {
  UV i, sum = 0, lim = rootint(n/m, r);

  if (r <= k) return lim;

  if (r-1 == k) {
    for (i = 1; i <= lim; i++)
      if (isf[i] && gcd_ui(m,i) == 1)
        sum += rootint(n/(m*ipow(i,r)),k);
  } else {
    for (i = 1; i <= lim; i++)
      if (isf[i] && gcd_ui(m,i) == 1)
        sum += _pcr(n, k, isf,  m * ipow(i,r), r-1);
  }
  return sum;
}

UV powerful_count(UV n, UV k) {
  UV i, r, lim, sum = 0;
  unsigned char *isf;

  if (k == 0) return 0;
  if (k == 1 || n <= 1) return n;
  if (k >= BITS_PER_WORD) return 1;

  lim = rootint(n, k+1);
  isf = _squarefree_range(0, lim);

  if (k == 2) {
    for (i = 1; i <= lim; i++)
      if (isf[i])
        sum += isqrt(n/(i*i*i));
  } else {
    /* sum = _pcr(n, k, isf,  1, 2*k-1); */
    r = 2*k-1;
    lim = rootint(n, r);
    for (i = 1; i <= lim; i++)
      if (isf[i])
        sum += _pcr(n, k, isf,  ipow(i,r), r-1);
  }

  Safefree(isf);
  return sum;
}

/* We want:
 *    k=0 turned into k=2 in XS (0 here ok)
 *    n=0 undef in XS (0 here ok)
 *    k=1 => n
 *    n=1 => 1
 *    n=2 => 1<<k
 *    overflow here should return 0
 */
UV nth_powerful(UV n, UV k) {
  static unsigned char const mink[20+1] = {0,0,1,2,4,6,7,9,11,12,14,16,18,19,21,23,24,26,28,30,31};
#if BITS_PER_WORD == 64
  static UV const maxpow[11] = {0,UV_MAX,9330124695,11938035,526402,85014,25017,10251,5137,2903,1796};
#else
  static UV const maxpow[11] = {0,UV_MAX,140008,6215,1373,536,281,172,115,79,57};
#endif
  UV lo, hi, max;
  double nc, npow, nest, dlo, dhi;

  if (k == 0 || k >= BITS_PER_WORD) return 0;
  if (k == 1 || n <= 1) return n;

  max = (k <= 10) ? maxpow[k] : powerful_count(UV_MAX,k);
  if (n > max) return 0;

  if (n <= 20 && k >= mink[n]) return UVCONST(1) << (k+(n-2));
  /* Now k >= 2, n >= 4 */

  nc = pow(n, 2) / pow(2.1732543125195541, 2);

  if (k == 2) { /* From Mincu and Panaitopol 2009 */
    npow = pow(n, 5.0/3.0);
    dlo = nc + 0.3 * npow;
    dhi = nc + 0.5 * npow;
    lo = (UV) dlo;
    hi = (n < 170) ? 8575 : (dhi >= UV_MAX) ? UV_MAX : 1 + (UV) dhi;
  } else if (k == 3) {
    /* Splitting the range is hacky but overall this isn't bad */
    if (n < 84000) {
      nest = .06003 * pow(n, 2.865);
      dlo = 0.96 * (nc + nest);
      dhi = 1.08 * (nc + nest);
    } else {
      nest = .02209 * pow(n, 2.955);
      dlo = 0.987 * (nc + nest);
      dhi = 1.020 * (nc + nest);
    }
    lo = (UV) dlo;
    if (n < 900) dhi *= 1.3;
    if (n < 160) dhi = 1.3 * dhi + 600;
    hi = (dhi >= UV_MAX) ? UV_MAX : 1 + (UV) dhi;
  } else if (k <= 10) {
    /* Slopppy but better than linear.  4 <= k <= 10. */
    if (n < 200) {
      npow = pow(n, 3.031 + 0.460*(k-4));
      nest = (.5462 / pow(1.15, k-4)) * npow;
      dlo = 0.51 * (nc + nest);
      dhi = 1.86 * (nc + nest);
    } else {
      npow = pow(n, 3.690 + 0.665*(k-4));
      nest = (.01275 / pow(4.11, k-4)) * npow;
      dlo = 0.76 * (nc + nest);
      dhi = 4.3 * (nc + nest);
    }
    lo = (UV) dlo;
    hi = (dhi >= UV_MAX) ? UV_MAX : 1 + (UV) dhi;
  } else {
    lo = (UVCONST(1) << (k+1))+1;
    hi = UV_MAX;
    /* Linear from min to max rather than a nice power fit as above */
    if (n < max)  hi = lo + (((double)n / (double)max) * (UV_MAX-lo) + 1);
  }

  return inverse_interpolate_k(lo, hi, n, k, &powerful_count, 0);
}


static UV _sumpowerful(UV m, UV r, UV n, UV k, unsigned char* isf) {
  UV i, rootdiv, sum;

  if (r < k)  return m;

  rootdiv = rootint(n/m, r);

  if (r == k)  return m * powersum(rootdiv, r);

  for (sum = 0, i = 1; i <= rootdiv; i++)
    if (isf[i] && gcd_ui(i,m) == 1)
      sum += _sumpowerful(m * ipow(i,r), r-1, n, k, isf);

  return sum;
}
UV sumpowerful(UV n, UV k)
{
  UV lim, sum;
  unsigned char *isf;

#if BITS_PER_WORD == 64
  static UV const maxpow[41] = {0,6074000999,8676161894447,263030040471727,
  1856371767674975,6768543273131775,17199267839999999,35098120384607510,
  62985599999999999,104857599999999999,157641800537109374,246512345193381887,
  312499999999999999,406381963906121727,499999999999999999,592297667290202111,
  701982420492091391,935976560656121855,1184595334580404223,1350851717672992088,
  1579460446107205631,2105947261476274175,2369190669160808447,
  4052555153018976266,4738381338321616895,7450580596923828124,
  7450580596923828124,7450580596923828124,
  UVCONST(9223372036854775807),UVCONST(9223372036854775807),
  UVCONST(9223372036854775807),UVCONST(9223372036854775807),
  UVCONST(9223372036854775807),UVCONST(9223372036854775807),
  UVCONST(9223372036854775807),UVCONST(9223372036854775807),
  UVCONST(9223372036854775807),UVCONST(9223372036854775807),
  UVCONST(9223372036854775807),UVCONST(9223372036854775807),
  UVCONST(12157665459056928800)};

  if (k == 0 || n == 0) return 0;
  if (k >= 64) return 1;
  if (k < 41 && n > maxpow[k]) return 0;
#else
  static UV const maxpow[21] = {0,92681,3367224,18224999,48599999,102036671,161243135,244140624,362797055,536870911,725594111,1088391167,1220703124,1220703124,2147483647,2147483647,2147483647,2147483647,2147483647,2147483647,3486784400};

  if (k == 0 || n == 0) return 0;
  if (k >= 32) return 1;
  if (k < 21 && n > maxpow[k]) return 0;
#endif

  if (k == 1)  return (n & 1)  ?  n*((n+1)>>1)  :  (n>>1)*(n+1);

  lim = rootint(n, k+1);
  isf = _squarefree_range(0, lim);
  sum = _sumpowerful(1, 2*k-1, n, k, isf);
  Safefree(isf);
  return sum;
}
