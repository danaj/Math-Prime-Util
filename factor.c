#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define FUNC_isqrt  1
#define FUNC_icbrt  1
#define FUNC_gcd_ui 1
#define FUNC_is_perfect_square 1
#define FUNC_clz 1
#include "ptypes.h"
#include "factor.h"
#include "sieve.h"
#include "util.h"
#include "mulmod.h"
#include "cache.h"
#include "primality.h"
#include "montmath.h"
static int holf32(uint32_t n, UV *factors, uint32_t rounds);

/*
 * You need to remember to use UV for unsigned and IV for signed types that
 * are large enough to hold our data.
 *   If you use int, that's 32-bit on LP64 and LLP64 machines.  You lose.
 *   If you use long, that's 32-bit on LLP64 machines.  You lose.
 *   If you use long long, you may be too large which isn't so bad, but some
 *                         compilers may not understand the type at all.
 * perl.h already figured all this out, and provided us with these types which
 * match the native integer type used inside our Perl, so just use those.
 */

static const unsigned short primes_small[] =
  {0,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
   101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,
   193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,
   293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,
   409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,
   521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,
   641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,
   757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,
   881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009,
   1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,1087,1091,1093,
   1097,1103,1109,1117,1123,1129,1151,1153,1163,1171,1181,1187,1193,1201,
   1213,1217,1223,1229,1231,1237,1249,1259,1277,1279,1283,1289,1291,1297,
   1301,1303,1307,1319,1321,1327,1361,1367,1373,1381,1399,1409,1423,1427,
   1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,1489,1493,1499,
   1511,1523,1531,1543,1549,1553,1559,1567,1571,1579,1583,1597,1601,1607,
   1609,1613,1619,1621,1627,1637,1657,1663,1667,1669,1693,1697,1699,1709,
   1721,1723,1733,1741,1747,1753,1759,1777,1783,1787,1789,1801,1811,1823,
   1831,1847,1861,1867,1871,1873,1877,1879,1889,1901,1907,1913,1931,1933,
   1949,1951,1973,1979,1987,1993,1997,1999,2003,2011};
#define NPRIMES_SMALL (sizeof(primes_small)/sizeof(primes_small[0]))

static int _small_trial_factor(UV n, UV *factors, UV *newn, uint32_t *lastf)
{
  int nfactors = 0;
  uint32_t f = 7;

  if (n > 1) {
    while ( (n & 1) == 0 ) { factors[nfactors++] = 2; n /= 2; }
    while ( (n % 3) == 0 ) { factors[nfactors++] = 3; n /= 3; }
    while ( (n % 5) == 0 ) { factors[nfactors++] = 5; n /= 5; }
  }

  if (f*f <= n) {
    uint32_t const lastsp = 83;
    uint32_t sp = 4;
    /* Trial division from 7 to 421.  Use 32-bit if possible. */
    if (n <= 4294967295U) {
      uint32_t un = n;
      while (sp < lastsp) {
        while ( (un%f) == 0 ) {
          factors[nfactors++] = f;
          un /= f;
        }
        f = primes_small[++sp];
        if (f*f > un) break;
      }
      n = un;
    } else {
      while (sp < lastsp) {
        while ( (n%f) == 0 ) {
          factors[nfactors++] = f;
          n /= f;
        }
        f = primes_small[++sp];
        if (f*f > n) break;
      }
    }
    /* If n is small and still composite, finish it here */
    if (n < 2011*2011 && f*f <= n) {  /* Trial division from 431 to 2003 */
      uint32_t un = n;
      while (sp < NPRIMES_SMALL) {
        while ( (un%f) == 0 ) {
          factors[nfactors++] = f;
          un /= f;
        }
        f = primes_small[++sp];
        if (f*f > un) break;
      }
      n = un;
    }
  }
  if (f*f > n && n != 1) {
    factors[nfactors++] = n;
    n = 1;
  }
  if (newn) *newn = n;
  if (lastf) *lastf = f;
  return nfactors;
}

static int _power_factor(UV n, UV *factors)
{
  int nfactors, i, j, k = powerof(n);
  if (k > 1) {
    UV p = rootint(n, k);
    nfactors = factor(p, factors);
    for (i = nfactors; i >= 0; i--)
      for (j = 0; j < k; j++)
        factors[k*i+j] = factors[i];
    return k*nfactors;
  }
  factors[0] = n;
  return 1;
}

/* Find one factor of an input n. */
int factor_one(UV n, UV *factors, int primality, int trial)
{
  int nfactors;
  if (n < 4) {
    factors[0] = n;
    return (n == 1) ? 0 : 1;
  }
  /* TODO: deal with small n */
  if (trial) {
    uint32_t sp, f;
    if (!(n&1)) { factors[0] = 2; factors[1] = n >> 1; return 2; }
    if (!(n%3)) { factors[0] = 3; factors[1] = n / 3; return 2; }
    if (!(n%5)) { factors[0] = 5; factors[1] = n / 5; return 2; }
    for (sp = 4; (f = primes_small[sp]) < 2011; sp++) {
      if ( (n % f) == 0 )
        { factors[0] = f; factors[1] = n/f; return 2; }
    }
    if (n < f*f) { factors[0] = n; return 1; }
  }
  if (primality && is_prime(n)) {
    factors[0] = n;
    return 1;
  }
#if 0  /* Simple solution, just fine on x86_64 */
  nfactors = (n < 1073741824UL) ? holf32(n, factors, 1000000)
                                : pbrent_factor(n, factors, 500000, 1);
  if (nfactors < 2) croak("factor_one failed on %lu\n", n);
#endif
  {
    /* Adjust the number of rounds based on the number size and speed */
    UV const nbits = BITS_PER_WORD - clz(n);
#if USE_MONTMATH
    UV const br_rounds = 8000 + (9000 * ((nbits <= 45) ? 0 : (nbits-45)));
    UV const sq_rounds = 200000;
#elif MULMODS_ARE_FAST
    UV const br_rounds =  100 + ( 100 * ((nbits <= 45) ? 0 : (nbits-45)));
    UV const sq_rounds = 100000;
#else
    UV const br_rounds = (nbits >= 63) ? 120000 : (nbits >= 58) ? 500 : 0;
    UV const sq_rounds = 200000;
#endif

#if BITS_PER_WORD == 64
    /* For small semiprimes the fastest solution is HOLF under 32, then
     * Lehman (no trial) under 38.  However on random inputs, HOLF is
     * best only under 28-30 bits, and adding Lehman is always slower. */
    if (nbits <= 30) { /* This should always succeed */
      nfactors = holf32(n, factors, 1000000);
      if (nfactors > 1) return nfactors;
    }
#endif
    /* Almost all inputs are factored here */
    if (br_rounds > 0) {
      nfactors = pbrent_factor(n, factors, br_rounds, 1);
      if (nfactors > 1) return nfactors;
    }
#if USE_MONTMATH
    nfactors = pbrent_factor(n, factors, 2*br_rounds, 3);
    if (nfactors > 1) return nfactors;
#endif
    /* Random 64-bit inputs at this point:
     *   About 3.1% are small enough that we did with HOLF.
     *   montmath:  96.89% pbrent,  0.01% pbrent2
     *   fast:      73.43% pbrent, 21.97% squfof, 1.09% p-1, 0.49% prho, long
     *   slow:      75.34% squfof, 19.47% pbrent, 0.20% p-1, 0.06% prho
     */
    /* SQUFOF with these parameters gets 99.9% of everything left */
    if (nbits <= 62) {
      nfactors = squfof_factor(n, factors, sq_rounds);
      if (nfactors > 1) return nfactors;
    }
    /* At this point we should only have 16+ digit semiprimes. */
    nfactors = pminus1_factor(n, factors, 8000, 120000);
    if (nfactors > 1) return nfactors;
    /* Get the stragglers */
    nfactors = prho_factor(n, factors, 120000);
    if (nfactors > 1) return nfactors;
    nfactors = pbrent_factor(n, factors, 500000, 5);
    if (nfactors > 1) return nfactors;
    nfactors = prho_factor(n, factors, 120000);
    if (nfactors > 1) return nfactors;
    croak("factor_one failed on %lu\n", n);
  }
  return nfactors;
}

/******************************************************************************/
/*                              Main factor loop                              */
/*                                                                            */
/* Puts factors in factors[] and returns the number found.                    */
/******************************************************************************/
int factor(UV n, UV *factors)
{
  UV tofac_stack[MPU_MAX_FACTORS+1];
  int nsmallfactors, npowerfactors, nfactors, i, j, ntofac = 0;
  uint32_t f;

  nfactors = _small_trial_factor(n, factors, &n, &f);
  if (n == 1) return nfactors;

#if BITS_PER_WORD == 64
  /* For small values less than f^3, use simple factor to split semiprime */
  if (n < 100000000 && n < f*f*f) {
    if (MR32(n)) factors[nfactors++] = n;
    else         nfactors += holf32(n, factors+nfactors, 10000);
    return nfactors;
  }
#endif

  nsmallfactors = nfactors;

  /* Perfect powers.  Factor root only once. */
  npowerfactors = _power_factor(n, factors+nsmallfactors);
  if (npowerfactors > 1) return nsmallfactors + npowerfactors;

  /* loop over each remaining factor, until ntofac == 0 */
  do {
    while ( (n >= f*f) && (!is_def_prime(n)) ) {
      int split_success = factor_one(n, tofac_stack+ntofac, 0, 0) - 1;
      if (split_success != 1 || tofac_stack[ntofac] == 1 || tofac_stack[ntofac] == n)
        croak("internal: factor_one failed to factor %"UVuf"\n", n);
      ntofac++; /* Leave one on the to-be-factored stack */
      n = tofac_stack[ntofac];  /* Set n to the other one */
    }
    /* n is now prime (or 1), so add to already-factored stack */
    if (n != 1)  factors[nfactors++] = n;
    /* Pop the next number off the to-factor stack */
    if (ntofac > 0)  n = tofac_stack[ntofac-1];
  } while (ntofac-- > 0);

  /* Sort the non-small factors */
  for (i = nsmallfactors+1; i < nfactors; i++) {
    UV fi = factors[i];
    for (j = i; j > 0 && factors[j-1] > fi; j--)
      factors[j] = factors[j-1];
    factors[j] = fi;
  }
  return nfactors;
}


int factor_exp(UV n, UV *factors, UV* exponents)
{
  int i = 1, j = 1, nfactors;

  if (n == 1) return 0;
  nfactors = factor(n, factors);

  if (exponents == 0) {
    for (; i < nfactors; i++)
      if (factors[i] != factors[i-1])
        factors[j++] = factors[i];
  } else {
    exponents[0] = 1;
    for (; i < nfactors; i++) {
      if (factors[i] != factors[i-1]) {
        exponents[j] = 1;
        factors[j++] = factors[i];
      } else {
        exponents[j-1]++;
      }
    }
  }
  return j;
}

int prime_bigomega(UV n)
{
  UV factors[MPU_MAX_FACTORS+1];
  return factor(n, factors);
}
int prime_omega(UV n)
{
  UV factors[MPU_MAX_FACTORS+1];
  int i, j, nfactors;

  if (n <= 1) return (n==0);
  nfactors = factor(n, factors);

  for (i = 1, j = 1; i < nfactors; i++)
    if (factors[i] != factors[i-1])
      j++;
  return j;
}


int trial_factor(UV n, UV *factors, UV f, UV last)
{
  int sp, nfactors = 0;

  if (f < 2) f = 2;
  if (last == 0 || last*last > n) last = UV_MAX;

  if (n < 4 || last < f) {
    factors[0] = n;
    return (n == 1) ? 0 : 1;
  }

  /* possibly do uint32_t specific code here */

  if (f < primes_small[NPRIMES_SMALL-1]) {
    while ( (n & 1) == 0 ) { factors[nfactors++] = 2; n >>= 1; }
    if (3<=last) while ( (n % 3) == 0 ) { factors[nfactors++] = 3; n /= 3; }
    if (5<=last) while ( (n % 5) == 0 ) { factors[nfactors++] = 5; n /= 5; }
    for (sp = 4; sp < (int)NPRIMES_SMALL; sp++) {
      f = primes_small[sp];
      if (f*f > n || f > last) break;
      while ( (n%f) == 0 ) {
        factors[nfactors++] = f;
        n /= f;
      }
    }
  }
  /* Trial division using a mod-30 wheel for larger values */
  if (f*f <= n && f <= last) {
    UV m, newlimit, limit = isqrt(n);
    if (limit > last) limit = last;
    m = f % 30;
    while (f <= limit) {
      if ( (n%f) == 0 ) {
        do {
          factors[nfactors++] = f;
          n /= f;
        } while ( (n%f) == 0 );
        newlimit = isqrt(n);
        if (newlimit < limit)  limit = newlimit;
      }
      f += wheeladvance30[m];
      m = nextwheel30[m];
    }
  }
  /* All done! */
  if (n != 1)
    factors[nfactors++] = n;
  return nfactors;
}


static void _divisors_from_factors(UV nfactors, UV* fp, UV* fe, UV* res) {
  UV s, count = 1;

  res[0] = 1;
  for (s = 0; s < nfactors; s++) {
    UV i, j, scount = count, p = fp[s], e = fe[s], mult = 1;
    for (j = 0; j < e; j++) {
      mult *= p;
      for (i = 0; i < scount; i++)
        res[count++] = res[i] * mult;
    }
  }
}

UV* _divisor_list(UV n, UV *num_divisors)
{
  UV factors[MPU_MAX_FACTORS+1];
  UV exponents[MPU_MAX_FACTORS+1];
  UV* divs;
  int i, nfactors, ndivisors;

  if (n <= 1) {
    New(0, divs, 2, UV);
    if (n == 0) {  divs[0] = 0;  divs[1] = 1;  *num_divisors = 2;  }
    if (n == 1) {  divs[0] = 1;                *num_divisors = 1;  }
    return divs;
  }
  /* Factor and convert to factor/exponent pair */
  nfactors = factor_exp(n, factors, exponents);
  /* Calculate number of divisors, allocate space, fill with divisors */
  ndivisors = exponents[0] + 1;
  for (i = 1; i < nfactors; i++)
    ndivisors *= (exponents[i] + 1);
  New(0, divs, ndivisors, UV);
  _divisors_from_factors(nfactors, factors, exponents, divs);
  /* Sort divisors (numeric ascending) */
  qsort(divs, ndivisors, sizeof(UV), _numcmp);
  /* Return number of divisors and list */
  *num_divisors = ndivisors;
  return divs;
}


/* The usual method, on OEIS for instance, is:
 *    (p^(k*(e+1))-1) / (p^k-1)
 * but that overflows quicky.  Instead we rearrange as:
 *    1 + p^k + p^k^2 + ... p^k^e
 * Return 0 if the result overflowed.
 */
static const UV sigma_overflow[11] =
#if BITS_PER_WORD == 64
         {UVCONST(3000000000000000000),UVCONST(3000000000),2487240,64260,7026,
         1622, 566, 256, 139, 85, 57};
#else
         {UVCONST(845404560), 52560, 1548, 252, 84, 41, 24, 16, 12, 10, 8};
#endif
UV divisor_sum(UV n, UV k)
{
  UV factors[MPU_MAX_FACTORS+1];
  int nfac, i, j;
  UV product = 1;

  if (k > 11 || (k > 0 && n >= sigma_overflow[k-1])) return 0;
  if (n <= 1)                               /* n=0  divisors are [0,1] */
    return (n == 1) ? 1 : (k == 0) ? 2 : 1; /* n=1  divisors are [1]   */
  nfac = factor(n,factors);
  if (k == 0) {
    for (i = 0; i < nfac; i++) {
      UV e = 1,  f = factors[i];
      while (i+1 < nfac && f == factors[i+1]) { e++; i++; }
      product *= (e+1);
    }
  } else if (k == 1) {
    for (i = 0; i < nfac; i++) {
      UV f = factors[i];
      UV pke = f, fmult = 1 + f;
      while (i+1 < nfac && f == factors[i+1]) {
        pke *= f;
        fmult += pke;
        i++;
      }
      product *= fmult;
    }
  } else {
    for (i = 0; i < nfac; i++) {
      UV f = factors[i];
      UV fmult, pke, pk = f;
      for (j = 1; j < (int)k; j++)  pk *= f;
      fmult = 1 + pk;
      pke = pk;
      while (i+1 < nfac && f == factors[i+1]) {
        pke *= pk;
        fmult += pke;
        i++;
      }
      product *= fmult;
    }
  }
  return product;
}




static int found_factor(UV n, UV f, UV* factors)
{
  UV f2 = n/f;
  int i = f > f2;
  if (f == 1 || f2 == 1) {
    factors[0] = n;
    return 1;
  }
  factors[i] = f;
  factors[1-i] = f2;
  MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
  return 2;
}

/* Knuth volume 2, algorithm C.
 * Can't compete with HOLF, SQUFOF, pbrent, etc.
 */
int fermat_factor(UV n, UV *factors, UV rounds)
{
  IV sqn, x, y, r;
  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in fermat_factor");
  sqn = isqrt(n);
  x = 2 * sqn + 1;
  y = 1;
  r = (sqn*sqn) - n;

  while (r != 0) {
    if (rounds-- == 0) { factors[0] = n; return 1; }
    r += x;
    x += 2;
    do {
      r -= y;
      y += 2;
    } while (r > 0);
  }
  r = (x-y)/2;
  return found_factor(n, r, factors);
}

/* Hart's One Line Factorization. */
int holf_factor(UV n, UV *factors, UV rounds)
{
  UV i, s, m, f;

  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in holf_factor");

  /* We skip the perfect-square test for s in the loop, so we
   * will never succeed if n is a perfect square.  Test that now. */
  if (is_perfect_square(n))
    return found_factor(n, isqrt(n), factors);

  if (n <= (UV_MAX >> 6)) {    /* Try with premultiplier first */
    UV npre = n * ( (n <= (UV_MAX >> 13)) ? 720 :
                    (n <= (UV_MAX >> 11)) ? 480 :
                    (n <= (UV_MAX >> 10)) ? 360 :
                    (n <= (UV_MAX >>  8)) ?  60 : 30 );
    UV ni = npre;
#if 0                              /* Straightforward */
    while (rounds--) {
      s = isqrt(ni) + 1;
      m = (s*s) - ni;
      if (is_perfect_square(m)) {
        f = gcd_ui(n, s - isqrt(m));
        if (f > 1 && f < n)
          return found_factor(n, f, factors);
      }
      if (ni >= (ni+npre)) break;
      ni += npre;
    }
#else                              /* More optimized */
    while (rounds--) {
      s = 1 + (UV)sqrt((double)ni);
      m = (s*s) - ni;
      f = m & 127;
      if (!((f*0x8bc40d7d) & (f*0xa1e2f5d1) & 0x14020a)) {
        f = (UV)sqrt((double)m);
        if (m == f*f) {
          f = gcd_ui(n, s - f);
          if (f > 1 && f < n)
            return found_factor(n, f, factors);
        }
      }
      if (ni >= (ni+npre)) break;
      ni += npre;
    }
#endif
  }

  for (i = 1; i <= rounds; i++) {
    s = (UV) sqrt( (double)n * (double)i );
    /* Assume s^2 isn't a perfect square.  We're rapidly losing precision
     * so we won't be able to accurately detect it anyway. */
    s++;    /* s = ceil(sqrt(n*i)) */
    m = sqrmod(s, n);
    if (is_perfect_square(m)) {
      f = isqrt(m);
      f = gcd_ui( (s>f) ? s-f : f-s, n);
      /* This should always succeed, but with overflow concerns.... */
      return found_factor(n, f, factors);
    }
  }
  factors[0] = n;
  return 1;
}
static int holf32(uint32_t n, UV *factors, uint32_t rounds) {
  UV npre, ni;    /* These should be 64-bit */
  uint32_t s, m, f;

  if (n < 3) { factors[0] = n; return 1; }
  if (!(n&1)) { factors[0] = 2; factors[1] = n/2; return 2; }
  if (is_perfect_square(n)) { factors[0] = factors[1] = isqrt(n); return 2; }

  ni = npre = (UV) n * ((BITS_PER_WORD == 64) ? 5040 : 1);
  while (rounds--) {
    s = 1 + (uint32_t)sqrt((double)ni);
    m = ((UV)s*(UV)s) - ni;
    f = m & 127;
    if (!((f*0x8bc40d7d) & (f*0xa1e2f5d1) & 0x14020a)) {
      f = (uint32_t)sqrt((double)m);
      if (m == f*f) {
        f = gcd_ui(n, s - f);
        if (f > 1 && f < n)
          return found_factor(n, f, factors);
      }
    }
    if (ni >= (ni+npre)) break; /* We've overflowed */
    ni += npre;
  }
  factors[0] = n;
  return 1;
}


#define ABSDIFF(x,y)  (x>y) ? x-y : y-x
#if USE_MONTMATH
/* Pollard Rho with Brent's updates, using Montgomery reduction. */
int pbrent_factor(UV n, UV *factors, UV rounds, UV a)
{
  UV const nbits = BITS_PER_WORD - clz(n);
  const UV inner = (nbits <= 31) ? 32 : (nbits <= 35) ? 64 : (nbits <= 40) ? 160 : (nbits <= 52) ? 256 : 320;
  UV f, m, r, rleft, Xi, Xm, Xs;
  int irounds, fails = 6;
  const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);

  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in pbrent_factor");
  r = f = 1;
  Xi = Xm = Xs = mont1;
  a = mont_geta(a,n);

  while (rounds > 0) {
    rleft = (r > rounds) ? rounds : r;
    Xm = Xi;
    /* Do rleft rounds, inner at a time */
    while (rleft > 0) {
      irounds = (rleft > (UV)inner) ? inner : rleft;
      rleft -= irounds;
      rounds -= irounds;
      Xs = Xi;
      if (n < (1ULL << 63)) {
        Xi = mont_mulmod63(Xi,Xi+a,n);
        m = ABSDIFF(Xi,Xm);
        while (--irounds > 0) {
          Xi = mont_mulmod63(Xi,Xi+a,n);
          f = ABSDIFF(Xi,Xm);
          m = mont_mulmod63(m, f, n);
        }
      } else if (a == mont1) {
        Xi = mont_mulmod64(Xi,Xi+a,n);
        m = ABSDIFF(Xi,Xm);
        while (--irounds > 0) {
          Xi = mont_mulmod64(Xi,Xi+a,n);
          f = ABSDIFF(Xi,Xm);
          m = mont_mulmod64(m, f, n);
        }
      } else {
        Xi = addmod(mont_mulmod64(Xi,Xi,n), a, n);
        m = ABSDIFF(Xi,Xm);
        while (--irounds > 0) {
          Xi = addmod(mont_mulmod64(Xi,Xi,n), a, n);
          f = ABSDIFF(Xi,Xm);
          m = mont_mulmod64(m, f, n);
        }
      }
      f = gcd_ui(m, n);
      if (f != 1)
        break;
    }
    /* If f == 1, then we didn't find a factor.  Move on. */
    if (f == 1) {
      r *= 2;
      continue;
    }
    if (f == n) {  /* back up, with safety */
      Xi = Xs;
      do {
        if (n < (1ULL << 63) || a == mont1)
          Xi = mont_mulmod(Xi,Xi+a,n);
        else
          Xi = addmod(mont_mulmod(Xi,Xi,n),a,n);
        m = ABSDIFF(Xi,Xm);
        f = gcd_ui(m, n);
      } while (f == 1 && r-- != 0);
    }
    if (f == 0 || f == n) {
      if (fails-- <= 0) break;
      Xi = Xm = mont1;
      a = addmod(a, mont_geta(11,n), n);
      continue;
    }
    return found_factor(n, f, factors);
  }
  factors[0] = n;
  return 1;
}
#else
/* Pollard Rho with Brent's updates. */
int pbrent_factor(UV n, UV *factors, UV rounds, UV a)
{
  UV f, m, r, Xi, Xm;
  const UV inner = (n <= 4000000000UL) ? 32 : 160;
  int fails = 6;

  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in pbrent_factor");

  r = f = Xi = Xm = 1;
  while (rounds > 0) {
    UV rleft = (r > rounds) ? rounds : r;
    UV saveXi = Xi;
    /* Do rleft rounds, inner at a time */
    while (rleft > 0) {
      UV dorounds = (rleft > inner) ? inner : rleft;
      saveXi = Xi;
      rleft -= dorounds;
      rounds -= dorounds;
      Xi = sqraddmod(Xi, a, n);        /* First iteration, no mulmod needed */
      m = ABSDIFF(Xi,Xm);
      while (--dorounds > 0) {         /* Now do inner-1=63 more iterations */
        Xi = sqraddmod(Xi, a, n);
        f = ABSDIFF(Xi,Xm);
        m = mulmod(m, f, n);
      }
      f = gcd_ui(m, n);
      if (f != 1)
        break;
    }
    /* If f == 1, then we didn't find a factor.  Move on. */
    if (f == 1) {
      r *= 2;
      Xm = Xi;
      continue;
    }
    if (f == n) {  /* back up, with safety */
      Xi = saveXi;
      do {
        Xi = sqraddmod(Xi, a, n);
        f = gcd_ui( ABSDIFF(Xi,Xm), n);
      } while (f == 1 && r-- != 0);
    }
    if (f == 0 || f == n) {
      if (fails-- <= 0) break;
      Xm = addmod(Xm, 11, n);
      Xi = Xm;
      a++;
      continue;
    }
    return found_factor(n, f, factors);
  }
  factors[0] = n;
  return 1;
}
#endif

/* Pollard's Rho. */
int prho_factor(UV n, UV *factors, UV rounds)
{
  UV a, f, i, m, oldU, oldV;
  const UV inner = 64;
  UV U = 7;
  UV V = 7;
  int fails = 3;

  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in prho_factor");

  /* We could just as well say a = 1 */
  switch (n%8) {
    case 1:  a = 1; break;
    case 3:  a = 2; break;
    case 5:  a = 3; break;
    case 7:  a = 5; break;
    default: a = 7; break;
  }

  rounds = (rounds + inner - 1) / inner;

  while (rounds-- > 0) {
    m = 1; oldU = U; oldV = V;
    for (i = 0; i < inner; i++) {
      U = sqraddmod(U, a, n);
      V = sqraddmod(V, a, n);
      V = sqraddmod(V, a, n);
      f = (U > V) ? U-V : V-U;
      m = mulmod(m, f, n);
    }
    f = gcd_ui(m, n);
    if (f == 1)
      continue;
    if (f == n) {  /* back up to find a factor*/
      U = oldU; V = oldV;
      i = inner;
      do {
        U = sqraddmod(U, a, n);
        V = sqraddmod(V, a, n);
        V = sqraddmod(V, a, n);
        f = gcd_ui( (U > V) ? U-V : V-U, n);
      } while (f == 1 && i-- != 0);
    }
    if (f == 0 || f == n) {
      if (fails-- <= 0) break;
      U = addmod(U,2,n);
      V = U;
      a++;
      continue;
    }
    return found_factor(n, f, factors);
  }
  factors[0] = n;
  return 1;
}

/* Pollard's P-1 */
int pminus1_factor(UV n, UV *factors, UV B1, UV B2)
{
  UV f, k, kmin;
  UV a     = 2, q     = 2;
  UV savea = 2, saveq = 2;
  UV j = 1;
  UV sqrtB1 = isqrt(B1);
#if USE_MONTMATH
  const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);
  UV ma = mont_geta(a,n);
  #define PMINUS1_APPLY_POWER  ma = mont_powmod(ma, k, n)
  #define PMINUS1_RECOVER_A    a = mont_recover(ma,n)
#else
  #define PMINUS1_APPLY_POWER  a = powmod(a, k, n)
  #define PMINUS1_RECOVER_A
#endif
  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in pminus1_factor");

  if (B1 <= primes_small[NPRIMES_SMALL-2]) {
    UV i;
    for (i = 1; primes_small[i] <= B1; i++) {
      q = k = primes_small[i];
      if (q <= sqrtB1) {
        k = q*q;  kmin = B1/q;
        while (k <= kmin)  k *= q;
      }
      PMINUS1_APPLY_POWER;
      if ( (j++ % 32) == 0) {
        PMINUS1_RECOVER_A;
        if (a == 0 || gcd_ui(a-1, n) != 1)
          break;
        savea = a;  saveq = q;
      }
    }
    PMINUS1_RECOVER_A;
  } else {
    START_DO_FOR_EACH_PRIME(2, B1) {
      q = k = p;
      if (q <= sqrtB1) {
        k = q*q;  kmin = B1/q;
        while (k <= kmin)  k *= q;
      }
      PMINUS1_APPLY_POWER;
      if ( (j++ % 32) == 0) {
        PMINUS1_RECOVER_A;
        if (a == 0 || gcd_ui(a-1, n) != 1)
          break;
        savea = a;  saveq = q;
      }
    } END_DO_FOR_EACH_PRIME
    PMINUS1_RECOVER_A;
  }
  if (a == 0) { factors[0] = n; return 1; }
  f = gcd_ui(a-1, n);

  /* If we found more than one factor in stage 1, backup and single step */
  if (f == n) {
    a = savea;
    START_DO_FOR_EACH_PRIME(saveq, B1) {
      k = p;  kmin = B1/p;
      while (k <= kmin)  k *= p;
      a = powmod(a, k, n);
      f = gcd_ui(a-1, n);
      q = p;
      if (f != 1)
        break;
    } END_DO_FOR_EACH_PRIME
    /* If f == n again, we could do:
     * for (savea = 3; f == n && savea < 100; savea = next_prime(savea)) {
     *   a = savea;
     *   for (q = 2; q <= B1; q = next_prime(q)) {
     *     ...
     *   }
     * }
     * but this could be a huge time sink if B1 is large, so just fail.
     */
  }

  /* STAGE 2 */
  if (f == 1 && B2 > B1) {
    UV bm = a;
    UV b = 1;
    UV bmdiff;
    UV precomp_bm[111] = {0};    /* Enough for B2 = 189M */

    /* calculate (a^q)^2, (a^q)^4, etc. */
    bmdiff = sqrmod(bm, n);
    precomp_bm[0] = bmdiff;
    for (j = 1; j < 20; j++) {
      bmdiff = mulmod(bmdiff,bm,n);
      bmdiff = mulmod(bmdiff,bm,n);
      precomp_bm[j] = bmdiff;
    }

    a = powmod(a, q, n);
    j = 1;
    START_DO_FOR_EACH_PRIME( q+1, B2 ) {
      UV lastq = q;
      UV qdiff;
      q = p;
      /* compute a^q = a^lastq * a^(q-lastq) */
      qdiff = (q - lastq) / 2 - 1;
      if (qdiff >= 111) {
        bmdiff = powmod(bm, q-lastq, n);  /* Big gap */
      } else {
        bmdiff = precomp_bm[qdiff];
        if (bmdiff == 0) {
          if (precomp_bm[qdiff-1] != 0)
            bmdiff = mulmod(mulmod(precomp_bm[qdiff-1],bm,n),bm,n);
          else
            bmdiff = powmod(bm, q-lastq, n);
          precomp_bm[qdiff] = bmdiff;
        }
      }
      a = mulmod(a, bmdiff, n);
      if (a == 0) break;
      b = mulmod(b, a-1, n);   /* if b == 0, we found multiple factors */
      if ( (j++ % 64) == 0 ) {
        f = gcd_ui(b, n);
        if (f != 1)
          break;
      }
    } END_DO_FOR_EACH_PRIME
    f = gcd_ui(b, n);
  }
  return found_factor(n, f, factors);
}

/* Simple Williams p+1 */
static void pp1_pow(UV *cX, UV exp, UV n)
{
  UV X0 = *cX;
  UV X  = *cX;
  UV Y = mulsubmod(X, X, 2, n);
  UV bit = UVCONST(1) << (clz(exp)-1);
  while (bit) {
    UV T = mulsubmod(X, Y, X0, n);
    if ( exp & bit ) {
      X = T;
      Y = mulsubmod(Y, Y, 2, n);
    } else {
      Y = T;
      X = mulsubmod(X, X, 2, n);
    }
    bit >>= 1;
  }
  *cX = X;
}
int pplus1_factor(UV n, UV *factors, UV B1)
{
  UV X1, X2, f;
  UV sqrtB1 = isqrt(B1);
  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in pplus1_factor");

  X1 =  7 % n;
  X2 = 11 % n;
  f = 1;
  START_DO_FOR_EACH_PRIME(2, B1) {
    UV k = p;
    if (p < sqrtB1) {
      UV kmin = B1/p;
      while (k <= kmin)
        k *= p;
    }
    pp1_pow(&X1, k, n);
    if (X1 != 2) {
      f = gcd_ui( submod(X1, 2, n) , n);
      if (f != 1 && f != n) break;
    }
    pp1_pow(&X2, k, n);
    if (X2 != 2) {
      f = gcd_ui( submod(X2, 2, n) , n);
      if (f != 1 && f != n) break;
    }
  } END_DO_FOR_EACH_PRIME

  return found_factor(n, f, factors);
}


/* SQUFOF, based on Ben Buhrow's racing version. */
#if 1
  /* limit to 62-bit inputs, use 32-bit types, faster */
  #define SQUFOF_TYPE uint32_t
  #define SQUFOF_MAX  (UV_MAX >> 2)
#else
  /* All 64-bit inputs possible, though we severely limit multipliers */
  #define SQUFOF_TYPE UV
  #define SQUFOF_MAX  UV_MAX
#endif
typedef struct
{
  int valid;
  SQUFOF_TYPE P;
  SQUFOF_TYPE bn;
  SQUFOF_TYPE Qn;
  SQUFOF_TYPE Q0;
  SQUFOF_TYPE b0;
  SQUFOF_TYPE it;
  SQUFOF_TYPE imax;
  SQUFOF_TYPE mult;
} mult_t;

/* N < 2^63 (or 2^31).  Returns 0 or a factor */
static UV squfof_unit(UV n, mult_t* mult_save)
{
  SQUFOF_TYPE imax,i,Q0,Qn,bn,b0,P,bbn,Ro,S,So,t1,t2;

  P = mult_save->P;
  bn = mult_save->bn;
  Qn = mult_save->Qn;
  Q0 = mult_save->Q0;
  b0 = mult_save->b0;
  i  = mult_save->it;
  imax = i + mult_save->imax;

#define SQUARE_SEARCH_ITERATION \
      t1 = P; \
      P = bn*Qn - P; \
      t2 = Qn; \
      Qn = Q0 + bn*(t1-P); \
      Q0 = t2; \
      bn = (b0 + P) / Qn; \
      i++;

  while (1) {
    int j = 0;
    if (i & 0x1) {
      SQUARE_SEARCH_ITERATION;
    }
    /* i is now even */
    while (1) {
      /* We need to know P, bn, Qn, Q0, iteration count, i  from prev */
      if (i >= imax) {
        /* save state and try another multiplier. */
        mult_save->P = P;
        mult_save->bn = bn;
        mult_save->Qn = Qn;
        mult_save->Q0 = Q0;
        mult_save->it = i;
        return 0;
      }

      SQUARE_SEARCH_ITERATION;

      /* Even iteration.  Check for square: Qn = S*S */
      t2 = Qn & 127;
      if (!((t2*0x8bc40d7d) & (t2*0xa1e2f5d1) & 0x14020a)) {
        t1 = (uint32_t) sqrt(Qn);
        if (Qn == t1*t1)
          break;
      }

      /* Odd iteration. */
      SQUARE_SEARCH_ITERATION;
    }
    S = t1; /* isqrt(Qn); */
    mult_save->it = i;

    /* Reduce to G0 */
    Ro = P + S*((b0 - P)/S);
    So = (n - (UV)Ro*(UV)Ro)/(UV)S;
    bbn = (b0+Ro)/So;

    /* Search for symmetry point */
#define SYMMETRY_POINT_ITERATION \
      t1 = Ro; \
      Ro = bbn*So - Ro; \
      t2 = So; \
      So = S + bbn*(t1-Ro); \
      S = t2; \
      bbn = (b0+Ro)/So; \
      if (Ro == t1) break;

    j = 0;
    while (1) {
      SYMMETRY_POINT_ITERATION;
      SYMMETRY_POINT_ITERATION;
      SYMMETRY_POINT_ITERATION;
      SYMMETRY_POINT_ITERATION;
      if (j++ > 2000000) {
         mult_save->valid = 0;
         return 0;
      }
    }

    t1 = gcd_ui(Ro, n);
    if (t1 > 1)
      return t1;
  }
}

/* Gower and Wagstaff 2008:
 *    http://www.ams.org/journals/mcom/2008-77-261/S0025-5718-07-02010-8/
 * Section 5.3.  I've added some with 13,17,19.  Sorted by F(). */
static const UV squfof_multipliers[] =
  /* { 3*5*7*11, 3*5*7, 3*5*11, 3*5, 3*7*11, 3*7, 5*7*11, 5*7,
       3*11,     3,     5*11,   5,   7*11,   7,   11,     1   }; */
  { 3*5*7*11, 3*5*7,  3*5*7*11*13, 3*5*7*13, 3*5*7*11*17, 3*5*11,
    3*5*7*17, 3*5,    3*5*7*11*19, 3*5*11*13,3*5*7*19,    3*5*7*13*17,
    3*5*13,   3*7*11, 3*7,         5*7*11,   3*7*13,      5*7,
    3*5*17,   5*7*13, 3*5*19,      3*11,     3*7*17,      3,
    3*11*13,  5*11,   3*7*19,      3*13,     5,           5*11*13,
    5*7*19,   5*13,   7*11,        7,        3*17,        7*13,
    11,       1 };
#define NSQUFOF_MULT (sizeof(squfof_multipliers)/sizeof(squfof_multipliers[0]))

int squfof_factor(UV n, UV *factors, UV rounds)
{
  mult_t mult_save[NSQUFOF_MULT];
  UV i, nn64, sqrtnn64, mult, f64,rounds_done = 0;
  int mults_racing = NSQUFOF_MULT;

  /* Caller should have handled these trivial cases */
  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in squfof_factor");

  /* Too big */
  if (n > SQUFOF_MAX) {
    factors[0] = n;  return 1;
  }

  for (i = 0; i < NSQUFOF_MULT; i++) {
    mult_save[i].valid = -1;
    mult_save[i].it = 0;
  }

  /* Race each multiplier for a bit (20-20k rounds) */
  while (mults_racing > 0 && rounds_done < rounds) {
    for (i = 0; i < NSQUFOF_MULT && rounds_done < rounds; i++) {
      if (mult_save[i].valid == 0)  continue;
      mult = squfof_multipliers[i];
      nn64 = n * mult;
      if (mult_save[i].valid == -1) {
        if ((SQUFOF_MAX / mult) < n) {
          mult_save[i].valid = 0; /* This multiplier would overflow 64-bit */
          mults_racing--;
          continue;
        }
        sqrtnn64 = isqrt(nn64);
        mult_save[i].valid = 1;
        mult_save[i].Q0    = 1;
        mult_save[i].b0    = sqrtnn64;
        mult_save[i].P     = sqrtnn64;
        mult_save[i].Qn    = (SQUFOF_TYPE)(nn64 - sqrtnn64 * sqrtnn64);
        if (mult_save[i].Qn == 0)
          return found_factor(n, sqrtnn64, factors);
        mult_save[i].bn    = (2 * sqrtnn64) / (UV)mult_save[i].Qn;
        mult_save[i].it    = 0;
        mult_save[i].mult  = mult;
        mult_save[i].imax  = (UV) (sqrt(sqrtnn64) / 16);
        if (mult_save[i].imax < 20)     mult_save[i].imax = 20;
        if (mult_save[i].imax > rounds) mult_save[i].imax = rounds;
      }
      if (mults_racing == 1)  /* Do all rounds if only one multiplier left */
        mult_save[i].imax = (rounds - rounds_done);
      f64 = squfof_unit(nn64, &mult_save[i]);
      if (f64 > 1) {
        UV f64red = f64 / gcd_ui(f64, mult);
        if (f64red > 1) {
          /* unsigned long totiter = 0;
             {int K; for (K = 0; K < NSQUFOF_MULT; K++) totiter += mult_save[K].it; }
             printf("  n %lu mult %lu it %lu (%lu)\n",n,mult,totiter,(UV)mult_save[i].it); */
          return found_factor(n, f64red, factors);
        }
        /* Found trivial factor.  Quit working with this multiplier. */
        mult_save[i].valid = 0;
      }
      if (mult_save[i].valid == 0)
        mults_racing--;
      rounds_done += mult_save[i].imax;   /* Assume we did all rounds */
    }
  }

  /* No factors found */
  factors[0] = n;
  return 1;
}

#define SQR_TAB_SIZE 512
static int sqr_tab_init = 0;
static double sqr_tab[SQR_TAB_SIZE];
static void make_sqr_tab(void) {
  int i;
  for (i = 0; i < SQR_TAB_SIZE; i++)
    sqr_tab[i] = sqrt((double)i);
  sqr_tab_init = 1;
}

/* Lehman written and tuned by Warren D. Smith.
 * Revised by Ben Buhrow and Dana Jacobsen. */
int lehman_factor(UV n, UV *factors, int do_trial) {
  const double Tune = ((n >> 31) >> 5) ? 3.5 : 5.0;
  double x, sqrtn;
  UV a,c,kN,kN4,B2;
  uint32_t b,p,k,r,B,U,Bred,inc,ip=2;

  if (!(n&1)) return found_factor(n, 2, factors);

  B = Tune * (1+rootint(n,3));

  if (do_trial) {
    uint32_t FirstCut = 0.1 * B;
    if (FirstCut < 84) FirstCut = 84;
    if (FirstCut > 65535) FirstCut = 65535;
    for (ip = 2;  ip < NPRIMES_SMALL; ip++) {
      p = primes_small[ip];
      if (p >= FirstCut)
        break;
      if (n % p == 0)
        return found_factor(n, p, factors);
    }
  }

  if (n >= UVCONST(8796393022207)) { factors[0] = n; return 1; }
  Bred = B / (Tune * Tune * Tune);
  B2 = B*B;
  kN = 0;

  if (!sqr_tab_init) make_sqr_tab();
  sqrtn = sqrt(n);

  for (k = 1; k <= Bred; k++) {
    if (k&1) { inc = 4; r = (k+n) % 4; }
    else     { inc = 2; r = 1; }
    kN += n;
    if (kN >= UVCONST(1152921504606846976)) { factors[0] = n; return 1; }
    kN4 = kN*4;

    x = (k < SQR_TAB_SIZE) ? sqrtn * sqr_tab[k] : sqrt((double)kN);
    a = x;
    if ((UV)a * (UV)a == kN)
      return found_factor(n, gcd_ui(a,n), factors);
    x *= 2;
    a = x + 0.9999999665;  /* Magic constant */
    b = a % inc;
    b = a + (inc+r-b) % inc;
    c = (UV)b*(UV)b - kN4;
    U = x + B2/(2*x);
    for (a = b;  a <= U;  c += inc*(a+a+inc), a += inc) {
      /* Check for perfect square */
      b = c & 127;
      if (!((b*0x8bc40d7d) & (b*0xa1e2f5d1) & 0x14020a)) {
        b = (uint32_t) sqrt(c);
        if (c == b*b) {
          B2 = gcd_ui(a+b, n);
          return found_factor(n, B2, factors);
        }
      }
    }
  }
  if (do_trial) {
    if (B > 65535) B = 65535;
    /* trial divide from primes[ip] to B.  We could:
     *   1) use table of 6542 shorts for the primes.
     *   2) use a wheel
     *   3) let trial_factor handle it
     */
    if (ip >= NPRIMES_SMALL)  ip = NPRIMES_SMALL-1;
    return trial_factor(n, factors, primes_small[ip], B);
  }
  factors[0] = n;
  return 1;
}

static const uint32_t _fr_chunk = 8192;
static const uint32_t _fr_sieve_crossover = 10000000;  /* About 10^14 */

static void _vec_factor(UV lo, UV hi, UV *nfactors, UV *farray, UV noffset, int square_free)
{
  UV *N, j, n, sqrthi, sievelim;
  sqrthi = isqrt(hi);
  n = hi-lo+1;
  New(0, N, hi-lo+1, UV);
  for (j = 0; j < n; j++) {
    N[j] = 1;
    nfactors[j] = 0;
  }
  sievelim = (sqrthi < _fr_sieve_crossover) ? sqrthi : icbrt(hi);
  START_DO_FOR_EACH_PRIME(2, sievelim) {
    UV q, t, A;
    if (square_free == 0) {
      UV kmin = hi / p;
      for (q = p; q <= kmin; q *= p) {
        t = lo / q, A = t * q;
        if (A < lo) A += q;
        for (j = A-lo; j < n; j += q) {
          farray[ j*noffset + nfactors[j]++ ] = p;
          N[j] *= p;
        }
      }
    } else {
      q = p*p, t = lo / q, A = t * q;
      if (A < lo) A += q;
      for (j = A-lo; j < n; j += q) {
        N[j] = 0;
        nfactors[j] = 0;
      }
      q = p, t = lo / q, A = t * q;
      if (A < lo) A += q;
      for (j = A-lo; j < n; j += q) {
        if (N[j] > 0) {
          farray[ j*noffset + nfactors[j]++ ] = p;
          N[j] *= p;
        }
      }
    }
  } END_DO_FOR_EACH_PRIME

  if (sievelim == sqrthi) {
    /* Handle the unsieved results, which are prime */
    for (j = 0; j < n; j++) {
      if (N[j] == 1)
        farray[ j*noffset + nfactors[j]++ ] = j+lo;
      else if (N[j] > 0 && N[j] != j+lo)
        farray[ j*noffset + nfactors[j]++ ] = (j+lo) / N[j];
    }
  } else {
    /* Handle the unsieved results, which are prime or semi-prime */
    for (j = 0; j < n; j++) {
      UV rem = j+lo;
      if (N[j] > 0 && N[j] != rem) {
        if (N[j] != 1)
          rem /= N[j];
        if (square_free && is_perfect_square(rem)) {
          nfactors[j] = 0;
        } else {
          UV* f = farray + j*noffset + nfactors[j];
          nfactors[j] += factor_one(rem, f, 1, 0);
        }
      }
    }
  }
  Safefree(N);
}

factor_range_context_t factor_range_init(UV lo, UV hi, int square_free) {
  factor_range_context_t ctx;
  ctx.lo = lo;
  ctx.hi = hi;
  ctx.n = lo-1;
  ctx.is_square_free = square_free ? 1 : 0;
  if (hi-lo+1 > 100) {        /* Sieve in chunks */
    if (square_free) ctx._noffset = (hi <= 42949672965UL) ? 10 : 15;
    else             ctx._noffset = BITS_PER_WORD - clz(hi);
    ctx._coffset = _fr_chunk;
    New(0, ctx._nfactors, _fr_chunk, UV);
    New(0, ctx._farray, _fr_chunk * ctx._noffset, UV);
    { /* Prealloc all the sieving primes now. */
      UV t = isqrt(hi);
      if (t >= _fr_sieve_crossover) t = icbrt(hi);
      get_prime_cache(t, 0);
    }
  } else {                    /* factor each number */
    New(0, ctx.factors, square_free ? 15 : 63, UV);
    ctx._nfactors = 0;
    ctx._farray = ctx.factors;
    ctx._noffset = 0;
  }
  return ctx;
}

int factor_range_next(factor_range_context_t *ctx) {
  int j, nfactors;
  UV n;
  if (ctx->n >= ctx->hi)
    return -1;
  n = ++(ctx->n);
  if (ctx->_nfactors) {
    if (ctx->_coffset >= _fr_chunk) {
      UV clo = n;
      UV chi = n + _fr_chunk - 1;
      if (chi > ctx->hi) chi = ctx->hi;
      _vec_factor(clo, chi, ctx->_nfactors, ctx->_farray, ctx->_noffset, ctx->is_square_free);
      ctx->_coffset = 0;
    }
    nfactors = ctx->_nfactors[ctx->_coffset];
    ctx->factors = ctx->_farray + ctx->_coffset * ctx->_noffset;
    ctx->_coffset++;
  } else {
    if (ctx->is_square_free && n >= 49 && (!(n% 4) || !(n% 9) || !(n%25) || !(n%49)))
      return 0;
    nfactors = factor(n, ctx->factors);
    if (ctx->is_square_free) {
      for (j = 1; j < nfactors; j++)
        if (ctx->factors[j] == ctx->factors[j-1])
          break;
      if (j < nfactors) return 0;
    }
  }
  return nfactors;
}

void factor_range_destroy(factor_range_context_t *ctx) {
  if (ctx->_farray != 0) Safefree(ctx->_farray);
  if (ctx->_nfactors != 0) Safefree(ctx->_nfactors);
  ctx->_farray = ctx->_nfactors = ctx->factors = 0;
}

/******************************************************************************/
/* Find number of factors for all values in a range */
/******************************************************************************/

unsigned char* range_nfactor_sieve(UV lo, UV hi, int with_multiplicity) {
  unsigned char* nf;
  UV *N, i, range = hi-lo+1, sqrtn = isqrt(hi);

  Newz(0, nf, range, unsigned char);
  New(0, N, range, UV);

  /* We could set to 1 and sieve from 2, or do this initialization */
  for (i = lo; i <= hi && i >= lo; i++) {
    N[i-lo] = 1;
    if (!(i&1) && i >= 2) {
      UV k = i >> 1;
      int nz = 1;
      while (!(k&1)) { nz++; k >>= 1; }
      nf[i-lo] = (with_multiplicity) ? nz : 1;
      N[i-lo] = UVCONST(1) << nz;
    }
  }

  START_DO_FOR_EACH_PRIME(3, sqrtn) {
    UV pk, maxpk = UV_MAX/p; \
    for (i = P_GT_LO_0(p,p,lo); i < range; i += p)
      { N[i] *= p;  nf[i]++; }
    for (pk = p*p; pk <= hi; pk *= p) {
      for (i = P_GT_LO_0(pk,pk,lo); i < range; i += pk)
        { N[i] *= p;  if (with_multiplicity) nf[i]++; }
      if (pk >= maxpk) break;  /* Overflow protection */
    }
  } END_DO_FOR_EACH_PRIME

  for (i = 0; i < range; i++)
    if (N[i] < (lo+i))
      nf[i]++;
  Safefree(N);
  if (lo == 0) nf[0] = 1;
  return nf;
}


/******************************************************************************/
/* DLP */
/******************************************************************************/

static UV dlp_trial(UV a, UV g, UV p, UV maxrounds) {
  UV k, t;
  if (maxrounds > p) maxrounds = p;

#if USE_MONTMATH
  if (p&1) {
    const uint64_t npi = mont_inverse(p),  mont1 = mont_get1(p);
    g = mont_geta(g, p);
    a = mont_geta(a, p);
    for (t = g, k = 1; k < maxrounds; k++) {
      if (t == a)
        return k;
      t = mont_mulmod(t, g, p);
      if (t == g) break;   /* Stop at cycle */
    }
  } else
#endif
  {
    for (t = g, k = 1; k < maxrounds; k++) {
      if (t == a)
        return k;
      t = mulmod(t, g, p);
      if (t == g) break;   /* Stop at cycle */
    }
  }
  return 0;
}

/******************************************************************************/
/* DLP - Pollard Rho */
/******************************************************************************/

/* Compare with Pomerance paper (dartmouth dtalk4):
 * Type I/II/III = our case 1, 0, 2.
 * x_i = u, a_i = v, b_i = w
 *
 * Also see Bai/Brent 2008 for many ideas to speed this up.
 * https://maths-people.anu.edu.au/~brent/pd/rpb231.pdf
 * E.g. Teske adding-walk, Brent's cycle algo, Teske modified cycle
 */
#define pollard_rho_cycle(u,v,w,p,n,a,g) \
    switch (u % 3) { \
      case 0: u = mulmod(u,u,p);  v = mulmod(v,2,n);  w = mulmod(w,2,n); break;\
      case 1: u = mulmod(u,a,p);  v = addmod(v,1,n);                     break;\
      case 2: u = mulmod(u,g,p);                      w = addmod(w,1,n); break;\
    }

typedef struct prho_state_t {
  UV u;
  UV v;
  UV w;
  UV U;
  UV V;
  UV W;
  UV round;
  int failed;
  int verbose;
} prho_state_t;

static UV dlp_prho_uvw(UV a, UV g, UV p, UV n, UV rounds, prho_state_t *s) {
  UV i, k = 0;
  UV u=s->u, v=s->v, w=s->w;
  UV U=s->U, V=s->V, W=s->W;
  int const verbose = s->verbose;

  if (s->failed)  return 0;
  if (s->round + rounds > n)  rounds = n - s->round;

  for (i = 1; i <= rounds; i++) {
    pollard_rho_cycle(u,v,w,p,n,a,g);   /* xi, ai, bi */
    pollard_rho_cycle(U,V,W,p,n,a,g);
    pollard_rho_cycle(U,V,W,p,n,a,g);   /* x2i, a2i, b2i */
    if (verbose > 3) printf( "%3"UVuf"  %4"UVuf" %3"UVuf" %3"UVuf"  %4"UVuf" %3"UVuf" %3"UVuf"\n", i, u, v, w, U, V, W );
    if (u == U) {
      UV r1, r2, G, G2;
      r1 = submod(v, V, n);
      if (r1 == 0) {
        if (verbose) printf("DLP Rho failure, r=0\n");
        s->failed = 1;
        k = 0;
        break;
      }
      r2 = submod(W, w, n);

      G = gcd_ui(r1,n);
      G2 = gcd_ui(G,r2);
      k = divmod(r2/G2, r1/G2, n/G2);
      if (G > 1) {
        if (powmod(g,k,p) == a) {
          if (verbose > 2) printf("  common GCD %"UVuf"\n", G2);
        } else {
          UV m, l = divmod(r2, r1, n/G);
          for (m = 0; m < G; m++) {
            k = addmod(l, mulmod(m,(n/G),n), n);
            if (powmod(g,k,p) == a) break;
          }
          if (m<G && verbose > 2) printf("  GCD %"UVuf", found with m=%"UVuf"\n", G, m);
        }
      }

      if (powmod(g,k,p) != a) {
        if (verbose > 2) printf("r1 = %"UVuf"  r2 = %"UVuf" k = %"UVuf"\n", r1, r2, k);
        if (verbose) printf("Incorrect DLP Rho solution: %"UVuf"\n", k);
        s->failed = 1;
        k = 0;
      }
      break;
    }
  }
  s->round += i-1;
  if (verbose && k) printf("DLP Rho solution found after %"UVuf" steps\n", s->round + 1);
  s->u = u; s->v = v; s->w = w; s->U = U; s->V = V; s->W = W;
  return k;
}

#if 0
static UV dlp_prho(UV a, UV g, UV p, UV n, UV maxrounds) {
#ifdef DEBUG
  int const verbose = _XS_get_verbose()
#else
  int const verbose = 0;
#endif
  prho_state_t s = {1, 0, 0, 1, 0, 0,   0, 0, verbose};
  return dlp_prho_uvw(a, g, p, n, maxrounds, &s);
}
#endif


/******************************************************************************/
/* DLP - BSGS */
/******************************************************************************/

typedef struct bsgs_hash_t {
  UV M;    /* The baby step index */
  UV V;    /* The powmod value */
  struct bsgs_hash_t* next;
} bsgs_hash_t;

/****************************************/
/*  Simple and limited pool allocation  */
#define BSGS_ENTRIES_PER_PAGE 8000
typedef struct bsgs_page_top_t {
  struct bsgs_page_t* first;
  bsgs_hash_t** table;
  UV  size;
  int nused;
  int npages;
} bsgs_page_top_t;

typedef struct bsgs_page_t {
  bsgs_hash_t entries[BSGS_ENTRIES_PER_PAGE];
  struct bsgs_page_t* next;
} bsgs_page_t;

static bsgs_hash_t* get_entry(bsgs_page_top_t* top) {
  if (top->nused == 0 || top->nused >= BSGS_ENTRIES_PER_PAGE) {
    bsgs_page_t* newpage;
    Newz(0, newpage, 1, bsgs_page_t);
    newpage->next = top->first;
    top->first = newpage;
    top->nused = 0;
    top->npages++;
  }
  return top->first->entries + top->nused++;
}
static void destroy_pages(bsgs_page_top_t* top) {
  bsgs_page_t* head = top->first;
  while (head != 0) {
    bsgs_page_t* next = head->next;
    Safefree(head);
    head = next;
  }
  top->first = 0;
}
/****************************************/

static void bsgs_hash_put(bsgs_page_top_t* pagetop, UV v, UV i) {
  UV idx = v % pagetop->size;
  bsgs_hash_t** table = pagetop->table;
  bsgs_hash_t* entry = table[idx];

  while (entry && entry->V != v)
    entry = entry->next;

  if (!entry) {
    entry = get_entry(pagetop);
    entry->M = i;
    entry->V = v;
    entry->next = table[idx];
    table[idx] = entry;
  }
}

static UV bsgs_hash_get(bsgs_page_top_t* pagetop, UV v) {
  bsgs_hash_t* entry = pagetop->table[v % pagetop->size];
  while (entry && entry->V != v)
    entry = entry->next;
  return (entry) ? entry->M : 0;
}

static UV bsgs_hash_put_get(bsgs_page_top_t* pagetop, UV v, UV i) {
  UV idx = v % pagetop->size;
  bsgs_hash_t** table = pagetop->table;
  bsgs_hash_t* entry = table[idx];

  while (entry && entry->V != v)
    entry = entry->next;

  if (entry)
    return entry->M;

  entry = get_entry(pagetop);
  entry->M = i;
  entry->V = v;
  entry->next = table[idx];
  table[idx] = entry;
  return 0;
}

static UV dlp_bsgs(UV a, UV g, UV p, UV n, UV maxent, int race_rho) {
  bsgs_page_top_t PAGES;
  UV i, m, maxm, hashmap_count;
  UV aa, S, gm, T, gs_i, bs_i;
  UV result = 0;
#ifdef DEBUG
  int const verbose = _XS_get_verbose();
#else
  int const verbose = 0;
#endif
  prho_state_t rho_state = {1, 0, 0, 1, 0, 0,   0, 0, verbose};

  if (n <= 2) return 0;   /* Shouldn't be here with gorder this low */

  if (race_rho) {
    result = dlp_prho_uvw(a, g, p, n, 10000, &rho_state);
    if (result) {
      if (verbose) printf("rho found solution in BSGS step 0\n");
      return result;
    }
  }

  if (a == 0) return 0;  /* We don't handle this case */

  maxm = isqrt(n);
  m = (maxent > maxm) ? maxm : maxent;

  hashmap_count = (m < 65537) ? 65537 :
                  (m > 40000000) ? 40000003 :
                  next_prime(m);               /* Ave depth around 2 */

  /* Create table.  Size: 8*hashmap_count bytes. */
  PAGES.size = hashmap_count;
  PAGES.first = 0;
  PAGES.nused = 0;
  PAGES.npages = 0;
  Newz(0, PAGES.table, hashmap_count, bsgs_hash_t*);

  aa = mulmod(a,a,p);
  S = a;
  gm = powmod(g, m, p);
  T = gm;
  gs_i = 0;
  bs_i = 0;

  bsgs_hash_put(&PAGES, S, 0);   /* First baby step */
  S = mulmod(S, g, p);
  /* Interleaved Baby Step Giant Step */
  for (i = 1; i <= m; i++) {
    gs_i = bsgs_hash_put_get(&PAGES, S, i);
    if (gs_i) { bs_i = i; break; }
    S = mulmod(S, g, p);
    if (S == aa) {  /* We discovered the solution! */
      if (verbose) printf("  dlp bsgs: solution at BS step %"UVuf"\n", i+1);
      result = i+1;
      break;
    }
    bs_i = bsgs_hash_put_get(&PAGES, T, i);
    if (bs_i) { gs_i = i; break; }
    T = mulmod(T, gm, p);
    if (race_rho && (i % 2048) == 0) {
      result = dlp_prho_uvw(a, g, p, n, 100000, &rho_state);
      if (result) {
        if (verbose) printf("rho found solution in BSGS step %"UVuf"\n", i);
        break;
      }
    }
  }

  if (!result) {
    /* Extend Giant Step search */
    if (!(gs_i || bs_i)) {
      UV b = (p+m-1)/m;
      if (m < maxm && b > 8*m) b = 8*m;
      for (i = m+1; i < b; i++) {
        bs_i = bsgs_hash_get(&PAGES, T);
        if (bs_i) { gs_i = i; break; }
        T = mulmod(T, gm, p);
        if (race_rho && (i % 2048) == 0) {
          result = dlp_prho_uvw(a, g, p, n, 100000, &rho_state);
          if (result) {
            if (verbose) printf("rho found solution in BSGS step %"UVuf"\n", i);
            break;
          }
        }
      }
    }

    if (gs_i || bs_i) {
      result = submod(mulmod(gs_i, m, p), bs_i, p);
    }
  }
  if (verbose) printf("  dlp bsgs using %d pages (%.1fMB+%.1fMB) for hash\n", PAGES.npages, ((double)PAGES.npages * sizeof(bsgs_page_t)) / (1024*1024), ((double)hashmap_count * sizeof(bsgs_hash_t*)) / (1024*1024));

  destroy_pages(&PAGES);
  Safefree(PAGES.table);
  if (result != 0 && powmod(g,result,p) != a) {
    if (verbose) printf("Incorrect DLP BSGS solution: %"UVuf"\n", result);
    result = 0;
  }
  if (race_rho && result == 0) {
    result = dlp_prho_uvw(a, g, p, n, 2000000000U, &rho_state);
  }
  return result;
}

/* Find smallest k where a = g^k mod p */
#define DLP_TRIAL_NUM  10000
UV znlog_solve(UV a, UV g, UV p, UV n) {
  UV k, sqrtn;
  const int verbose = _XS_get_verbose();

  if (a >= p) a %= p;
  if (g >= p) g %= p;

  if (a == 1 || g == 0 || p <= 2)
    return 0;

  if (verbose > 1 && n != p-1) printf("  g=%"UVuf" p=%"UVuf", order %"UVuf"\n", g, p, n);

  /* printf(" solving znlog(%"UVuf",%"UVuf",%"UVuf") n=%"UVuf"\n", a, g, p, n); */

  if (n == 0 || n <= DLP_TRIAL_NUM) {
    k = dlp_trial(a, g, p, DLP_TRIAL_NUM);
    if (verbose) printf("  dlp trial 10k %s\n", (k!=0 || p <= DLP_TRIAL_NUM) ? "success" : "failure");
    if (k != 0 || (n > 0 && n <= DLP_TRIAL_NUM)) return k;
  }

  { /* Existence checks */
    UV aorder, gorder = n;
    if (gorder != 0 && powmod(a, gorder, p) != 1) return 0;
    aorder = znorder(a,p);
    if (aorder == 0 && gorder != 0) return 0;
    if (aorder != 0 && gorder % aorder != 0) return 0;
  }

  sqrtn = (n == 0) ? 0 : isqrt(n);
  if (n == 0) n = p-1;

  {
    UV maxent = (sqrtn > 0) ? sqrtn+1 : 100000;
    k = dlp_bsgs(a, g, p, n, maxent/2, /* race rho */ 1);
    if (verbose) printf("  dlp bsgs %"UVuf"k %s\n", maxent/1000, k!=0 ? "success" : "failure");
    if (k != 0) return k;
    if (sqrtn > 0 && sqrtn < maxent) return 0;
  }

  if (verbose) printf("  dlp doing exhaustive trial\n");
  k = dlp_trial(a, g, p, p);
  return k;
}

/* Silver-Pohlig-Hellman */
static UV znlog_ph(UV a, UV g, UV p, UV p1) {
  UV fac[MPU_MAX_FACTORS+1];
  UV exp[MPU_MAX_FACTORS+1];
  int i, nfactors;
  UV x, j;

  if (p1 == 0) return 0;   /* TODO: Should we plow on with p1=p-1? */
  nfactors = factor_exp(p1, fac, exp);
  if (nfactors == 1)
    return znlog_solve(a, g, p, p1);
  for (i = 0; i < nfactors; i++) {
    UV pi, delta, gamma;
    pi = fac[i];   for (j = 1; j < exp[i]; j++)  pi *= fac[i];
    delta = powmod(a,p1/pi,p);
    gamma = powmod(g,p1/pi,p);
    /* printf(" solving znlog(%"UVuf",%"UVuf",%"UVuf")\n", delta, gamma, p); */
    fac[i] = znlog_solve( delta, gamma, p, znorder(gamma,p) );
    exp[i] = pi;
  }
  if (chinese(&x, 0, fac, exp, nfactors) == 1 && powmod(g, x, p) == a)
    return x;
  return 0;
}

/* Find smallest k where a = g^k mod p */
UV znlog(UV a, UV g, UV p) {
  UV k, gorder, aorder;
  const int verbose = _XS_get_verbose();

  if (a >= p) a %= p;
  if (g >= p) g %= p;

  if (a == 1 || g == 0 || p <= 2)
    return 0;

  /* TODO: We call znorder with the same p many times.  We should have a
   * method for znorder given {phi,nfactors,fac,exp} */

  gorder = znorder(g,p);
  if (gorder != 0 && powmod(a, gorder, p) != 1) return 0;
  /* TODO: Can these tests every fail?  Do we need aorder? */
  aorder = znorder(a,p);
  if (aorder == 0 && gorder != 0) return 0;
  if (aorder != 0 && gorder % aorder != 0) return 0;

  /* TODO: Come up with a better solution for a=0 */
  if (a == 0 || p < DLP_TRIAL_NUM || (gorder > 0 && gorder < DLP_TRIAL_NUM)) {
    if (verbose > 1) printf("  dlp trial znlog(%"UVuf",%"UVuf",%"UVuf")\n",a,g,p);
    k = dlp_trial(a, g, p, p);
    return k;
  }

  if (!is_prob_prime(gorder)) {
    k = znlog_ph(a, g, p, gorder);
    if (verbose) printf("  dlp PH %s\n", k!=0 ? "success" : "failure");
    if (k != 0) return k;
  }

  return znlog_solve(a, g, p, gorder);
}


/* Compile with:
 *  gcc -O3 -fomit-frame-pointer -march=native -Wall -DSTANDALONE -DFACTOR_STANDALONE factor.c util.c primality.c cache.c sieve.c chacha.c csprng.c prime_counts.c prime_count_cache.c lmo.c legendre_phi.c real.c inverse_interpolate.c -lm
 */
#ifdef FACTOR_STANDALONE
#include <errno.h>
int main(int argc, char *argv[])
{
  UV n;
  UV factors[MPU_MAX_FACTORS+1];
  int nfactors, i, a;

  if (argc <= 1) {
    char line[1024];
    while (1) {
      if (!fgets(line,sizeof(line),stdin)) break;
      n = strtoull(line, 0, 10);
      nfactors = factor(n, factors);
      if (nfactors == 1) {
        printf("%"UVuf": %"UVuf"\n",n,n);
      } else if (nfactors == 2) {
        printf("%"UVuf": %"UVuf" %"UVuf"\n",n,factors[0],factors[1]);
      } else if (nfactors == 3) {
        printf("%"UVuf": %"UVuf" %"UVuf" %"UVuf"\n",n,factors[0],factors[1],factors[2]);
      } else {
        printf("%"UVuf": %"UVuf" %"UVuf" %"UVuf" %"UVuf"",n,factors[0],factors[1],factors[2],factors[3]);
        for (i = 4; i < nfactors; i++) printf(" %"UVuf"", factors[i]);
        printf("\n");
      }
    }
    exit(0);
  }

  for (a = 1; a < argc; a++) {
    n = strtoul(argv[a], 0, 10);
    if (n == ULONG_MAX && errno == ERANGE) { printf("Argument larger than ULONG_MAX\n"); return(-1); }
    nfactors = factor(n, factors);
    printf("%"UVuf":", n);
    for (i = 0; i < nfactors; i++)
      printf(" %"UVuf"", factors[i]);
    printf("\n");
  }

  return(0);
}
#endif
