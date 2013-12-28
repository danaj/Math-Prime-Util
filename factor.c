#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ptypes.h"
#include "factor.h"
#include "sieve.h"
#include "mulmod.h"
#include "cache.h"
#include "primality.h"
#define FUNC_isqrt  1
#define FUNC_gcd_ui 1
#include "util.h"

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

/* The main factoring loop */
/* Puts factors in factors[] and returns the number found. */
int factor(UV n, UV *factors)
{
  int nfactors = 0;           /* Number of factored in factors result */

  int const verbose = _XS_get_verbose();
  UV const tlim_lower = 401;  /* Trial division through this prime */
  UV const tlim = 409;        /* This means we've checked through here */
  UV tofac_stack[MPU_MAX_FACTORS+1];
  UV fac_stack[MPU_MAX_FACTORS+1];
  int ntofac = 0;             /* Number of items on tofac_stack */
  int nfac = 0;               /* Number of items on fac_stack */

  if (n < 10000000)
    return trial_factor(n, factors, 0);

  /* Trial division for all factors below tlim */
  nfactors = trial_factor(n, factors, tlim_lower);
  n = factors[--nfactors];

  /* loop over each remaining factor, until ntofac == 0 */
  do {
    while ( (n >= (tlim*tlim)) && (!_XS_is_prime(n)) ) {
      int split_success = 0;
      /* Adjust the number of rounds based on the number size */
      UV const br_rounds = ((n>>29) < 100000) ?  1500 :  1500;
      UV const sq_rounds = 80000; /* 20k 91%, 40k 98%, 80k 99.9%, 120k 99.99% */

      /* 99.7% of 32-bit, 94% of 64-bit random inputs factored here */
      if (!split_success) {
        split_success = pbrent_factor(n, tofac_stack+ntofac, br_rounds, 3)-1;
        if (verbose) { if (split_success) printf("pbrent 1:  %"UVuf" %"UVuf"\n", tofac_stack[ntofac], tofac_stack[ntofac+1]); else printf("pbrent 0\n"); }
      }
      /* SQUFOF with these parameters gets 99.9% of everything left */
      if (!split_success && n < (UV_MAX>>2)) {
        split_success = racing_squfof_factor(n,tofac_stack+ntofac, sq_rounds)-1;
        if (verbose) printf("rsqufof %d\n", split_success);
      }
      /* At this point we should only have 16+ digit semiprimes. */
      /* This p-1 gets about 2/3 of what makes it through the above */
      if (!split_success) {
        split_success = pminus1_factor(n, tofac_stack+ntofac, 5000, 80000)-1;
        if (verbose) printf("pminus1 %d\n", split_success);
      }
      /* Some rounds of HOLF, good for close to perfect squares which are
       * the worst case for the next step */
      if (!split_success) {
        split_success = holf_factor(n, tofac_stack+ntofac, 2000)-1;
        if (verbose) printf("holf %d\n", split_success);
      }
      /* The catch-all.  Should factor anything. */
      if (!split_success) {
        split_success = prho_factor(n, tofac_stack+ntofac, 256*1024)-1;
        if (verbose) printf("long prho %d\n", split_success);
      }

      if (split_success) {
        MPUassert( split_success == 1, "split factor returned more than 2 factors");
        ntofac++; /* Leave one on the to-be-factored stack */
        if ((tofac_stack[ntofac] == n) || (tofac_stack[ntofac] == 1))
          croak("bad factor\n");
        n = tofac_stack[ntofac];  /* Set n to the other one */
      } else {
        /* Factor via trial division.  Nothing should make it here. */
        UV f = tlim;
        UV m = tlim % 30;
        UV limit = isqrt(n);
        if (verbose) printf("doing trial on %"UVuf"\n", n);
        while (f <= limit) {
          if ( (n%f) == 0 ) {
            do {
              n /= f;
              fac_stack[nfac++] = f;
            } while ( (n%f) == 0 );
            limit = isqrt(n);
          }
          f += wheeladvance30[m];
          m =  nextwheel30[m];
        }
        break;  /* We just factored n via trial division.  Exit loop. */
      }
    }
    /* n is now prime (or 1), so add to already-factored stack */
    if (n != 1)  fac_stack[nfac++] = n;
    /* Pop the next number off the to-factor stack */
    if (ntofac > 0)  n = tofac_stack[ntofac-1];
  } while (ntofac-- > 0);

  /* Sort all the results from fac_stack and put into factors result */
  {
    int i, j;
    for (i = 0; i < nfac; i++) {
      int mini = i;
      for (j = i+1; j < nfac; j++)
        if (fac_stack[j] < fac_stack[mini])
          mini = j;
      if (mini != i) {
        UV t = fac_stack[mini];
        fac_stack[mini] = fac_stack[i];
        fac_stack[i] = t;
      }
      factors[nfactors++] = fac_stack[i];
    }
  }
  return nfactors;
}

int factor_exp(UV n, UV *factors, UV* exponents)
{
  int i, j, nfactors = factor(n, factors);

  if (n == 1) return 0;
  exponents[0] = 1;
  for (i = 1, j = 1; i < nfactors; i++) {
    if (factors[i] != factors[i-1]) {
      exponents[j] = 1;
      factors[j++] = factors[i];
    } else {
      exponents[j-1]++;
    }
  }
  return j;
}



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

int trial_factor(UV n, UV *factors, UV maxtrial)
{
  UV f, limit, newlimit;
  int nfactors = 0;

  if (maxtrial == 0)  maxtrial = UV_MAX;

  /* Cover the cases 0/1/2/3 now */
  if (n < 4 || maxtrial < 2) {
    factors[0] = n;
    return (n == 1) ? 0 : 1;
  }
  /* Trial division for 2, 3, 5, 7, and see if we're done */
  while ( (n & 1) == 0 ) { factors[nfactors++] = 2; n /= 2; }
  if (3<=maxtrial) while ( (n % 3) == 0 ) { factors[nfactors++] = 3; n /= 3; }
  if (5<=maxtrial) while ( (n % 5) == 0 ) { factors[nfactors++] = 5; n /= 5; }
  if (7<=maxtrial) while ( (n % 7) == 0 ) { factors[nfactors++] = 7; n /= 7; }
  f = 11;
  if ( (n < (f*f)) || (maxtrial < f) ) {
    if (n != 1)
      factors[nfactors++] = n;
    return nfactors;
  }

  /* Trial division to this number at most.  Reduced as we find factors. */
  limit = isqrt(n);
  if (limit > maxtrial)
    limit = maxtrial;

  /* Use the table of small primes to quickly do trial division. */
  {
    UV sp = 5;
    UV slimit = (limit < 2003) ? limit : 2003;
    f = primes_small[sp];
    while (f <= slimit) {
      if ( (n%f) == 0 ) {
        do {
          factors[nfactors++] = f;
          n /= f;
        } while ( (n%f) == 0 );
        newlimit = isqrt(n);
        if (newlimit < slimit)  slimit = newlimit;
        if (newlimit < limit)   limit = newlimit;
      }
      f = primes_small[++sp];
    }
  }

  /* Trial division using a mod-30 wheel for larger values */
  if (f <= limit) {
    UV m = f % 30;
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


/* Return 0 if n is not a perfect square.  Set sqrtn to int(sqrt(n)) if so.
 *
 * Some simple solutions:
 *
 *     return ( ((n&2)!= 0) || ((n&7)==5) || ((n&11) == 8) )  ?  0  :  1;
 *
 * or:
 *
 *     m = n & 31;
 *     if ( m==0 || m==1 || m==4 || m==9 || m==16 || m==17 || m==25 )
 *       ...test for perfect square...
 *
 * or:
 *
 *     if (  ((0x0202021202030213ULL >> (n & 63)) & 1) &&
 *           ((0x0402483012450293ULL >> (n % 63)) & 1) &&
 *           ((0x218a019866014613ULL >> ((n % 65) & 63)) & 1) &&
 *           ((0x23b                 >> (n % 11) & 1)) ) {
 *
 *
 * The following Bloom filter cascade works very well indeed.  Read all
 * about it here: http://mersenneforum.org/showpost.php?p=110896
 */
static int is_perfect_square(UV n, UV* sqrtn)
{
  UV m;
  m = n & 127;
  if ((m*0x8bc40d7d) & (m*0xa1e2f5d1) & 0x14020a)  return 0;
  /* 82% of non-squares rejected here */

#if 0
  /* The big deal with this technique is that you do two total operations,
   * one cheap (the & 127 above), one expensive (the modulo below) on n.
   * The rest of the operations are 32-bit operations.  This is a huge win
   * if n is multiprecision.
   * However, in this file we're doing native precision sqrt, so it just
   * isn't expensive enough to justify this second filter set.
   */
  lm = n % UVCONST(63*25*11*17*19*23*31);
  m = lm % 63;
  if ((m*0x3d491df7) & (m*0xc824a9f9) & 0x10f14008) return 0;
  m = lm % 25;
  if ((m*0x1929fc1b) & (m*0x4c9ea3b2) & 0x51001005) return 0;
  m = 0xd10d829a*(lm%31);
  if (m & (m+0x672a5354) & 0x21025115) return 0;
  m = lm % 23;
  if ((m*0x7bd28629) & (m*0xe7180889) & 0xf8300) return 0;
  m = lm % 19;
  if ((m*0x1b8bead3) & (m*0x4d75a124) & 0x4280082b) return 0;
  m = lm % 17;
  if ((m*0x6736f323) & (m*0x9b1d499) & 0xc0000300) return 0;
  m = lm % 11;
  if ((m*0xabf1a3a7) & (m*0x2612bf93) & 0x45854000) return 0;
  /* 99.92% of non-squares are rejected now */
#endif
#if 0
  /* This could save time on some platforms, but not on x86 */
  m = n % 63;
  if ((m*0x3d491df7) & (m*0xc824a9f9) & 0x10f14008) return 0;
#endif
  m = isqrt(n);
  if (n != (m*m))
    return 0;

  if (sqrtn != 0) *sqrtn = m;
  return 1;
}

static int _divisors_from_factors(UV v, UV npe, UV* fp, UV* fe, UV* res) {
  UV p, e, i;
  if (npe == 0) return 0;
  p = *fp++;
  e = *fe++;
  if (npe == 1) {
    for (i = 0; i <= e; i++) {
      *res++ = v;
      v *= p;
    }
    return e+1;
  } else {
    int nret = 0;;
    for (i = 0; i <= e; i++) {
      int nres = _divisors_from_factors(v, npe-1, fp, fe, res);
      v *= p;
      res += nres;
      nret += nres;
    }
    return nret;
  }
}

UV* _divisor_list(UV n, UV *num_divisors)
{
  UV factors[MPU_MAX_FACTORS+1];
  UV exponents[MPU_MAX_FACTORS+1];
  UV* divs;
  int i, j, nfactors, ndivisors;

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
  (void) _divisors_from_factors(1, nfactors, factors, exponents, divs);
  { /* Sort (Shell sort is easy and efficient) */
    static int gaps[] = {301, 132, 57, 23, 10, 4, 1, 0};
    int gap, gapi = 0;
    for (gap = gaps[gapi]; gap > 0; gap = gaps[++gapi]) {
      for (i = gap; i < ndivisors; i++) {
        UV v = divs[i];
        for (j = i; j >= gap && divs[j-gap] > v; j -= gap)
          divs[j] = divs[j-gap];
        divs[j] = v;
      }
    }
  }
  *num_divisors = ndivisors;
  return divs;
}


/*
 * The usual method, on OEIS for instance, is:
 *    (p^(k*(e+1))-1) / (p^k-1)
 * but that overflows quicky.  Instead we rearrange as:
 *    1 + p^k + p^k^2 + ... p^k^e
 */
UV _XS_divisor_sum(UV n, UV k)
{
  UV factors[MPU_MAX_FACTORS+1];
  int nfac, i, j;
  UV product = 1;

  if (n == 0) return (k == 0) ? 2 : 1;  /* divisors are [0,1] */
  if (n == 1) return 1;                 /* divisors are [1]   */
  nfac = factor(n, factors);
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




/* Knuth volume 2, algorithm C.
 * Very fast for small numbers, grows rapidly.
 * SQUFOF is better for numbers nearing the 64-bit limit.
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
    r += x;
    x += 2;
    do {
      r -= y;
      y += 2;
    } while (r > 0);
  }
  r = (x-y)/2;
  if ( (r != 1) && ((UV)r != n) ) {
    factors[0] = r;
    factors[1] = n/r;
    MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
    return 2;
  }
  factors[0] = n;
  return 1;
}

/* Hart's One Line Factorization.
 * Missing premult (hard to do in native precision without overflow)
 */
int holf_factor(UV n, UV *factors, UV rounds)
{
  UV i, s, m, f;

  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in holf_factor");

  for (i = 1; i <= rounds; i++) {
    s = (UV) sqrt( (double)n * (double)i );
    /* Assume s^2 isn't a perfect square.  We're rapidly losing precision
     * so we won't be able to accurately detect it anyway. */
    s++;    /* s = ceil(sqrt(n*i)) */
    m = sqrmod(s, n);
    if (is_perfect_square(m, &f)) {
      f = gcd_ui( (s>f) ? s-f : f-s, n);
      /* This should always succeed, but with overflow concerns.... */
      if ((f == 1) || (f == n))
        break;
      factors[0] = f;
      factors[1] = n/f;
      MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
      return 2;
    }
  }
  factors[0] = n;
  return 1;
}


/* Pollard / Brent.  Brent's modifications to Pollard's Rho.  Maybe faster. */
int pbrent_factor(UV n, UV *factors, UV rounds, UV a)
{
  UV f, i, r;
  UV Xi = 2;
  UV Xm = 2;
  const UV inner = 64;

  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in pbrent_factor");

  r = 1;
  while (rounds > 0) {
    UV rleft = (r > rounds) ? rounds : r;
    UV saveXi;
    /* Do rleft rounds, inner at a time */
    while (rleft > 0) {
      UV dorounds = (rleft > inner) ? inner : rleft;
      UV m = 1;
      saveXi = Xi;
      for (i = 0; i < dorounds; i++) {
        Xi = sqraddmod(Xi, a, n);
        f = (Xi>Xm) ? Xi-Xm : Xm-Xi;
        m = mulmod(m, f, n);
      }
      rleft -= dorounds;
      rounds -= dorounds;
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
        f = gcd_ui( (Xi>Xm) ? Xi-Xm : Xm-Xi, n);
      } while (f == 1 && r-- != 0);
      if ( (f == 1) || (f == n) ) break;
    }
    factors[0] = f;
    factors[1] = n/f;
    MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
    return 2;
  }
  factors[0] = n;
  return 1;
}

/* Pollard's Rho. */
int prho_factor(UV n, UV *factors, UV rounds)
{
  UV a, f, i, m, oldU, oldV;
  const UV inner = 64;
  UV U = 7;
  UV V = 7;

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
      if ( (f == 1) || (f == n) )
        break;
    }
    factors[0] = f;
    factors[1] = n/f;
    MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
    return 2;
  }
  factors[0] = n;
  return 1;
}

/* Pollard's P-1 */
int pminus1_factor(UV n, UV *factors, UV B1, UV B2)
{
  UV f;
  UV q = 2;
  UV a = 2;
  UV savea = 2;
  UV saveq = 2;
  UV j = 1;
  UV sqrtB1 = isqrt(B1);
  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in pminus1_factor");

  START_DO_FOR_EACH_PRIME(2, sqrtB1) {
    UV k = p*p;
    UV kmin = B1/p;
    while (k <= kmin)
      k *= p;
    a = powmod(a, k, n);
    q = p;
  } END_DO_FOR_EACH_PRIME
  if (a == 0) { factors[0] = n; return 1; }
  f = gcd_ui(a-1, n);
  if (f == 1) {
    savea = a;
    saveq = q;
    START_DO_FOR_EACH_PRIME(q+1, B1) {
      q = p;
      a = powmod(a, q, n);
      if ( (j++ % 32) == 0) {
        if (a == 0 || gcd_ui(a-1, n) != 1)
          break;
        savea = a;
        saveq = q;
      }
    } END_DO_FOR_EACH_PRIME
    if (a == 0) { factors[0] = n; return 1; }
    f = gcd_ui(a-1, n);
  }
  /* If we found more than one factor in stage 1, backup and single step */
  if (f == n) {
    a = savea;
    START_DO_FOR_EACH_PRIME(saveq, B1) {
      UV k = p;
      UV kmin = B1/p;
      while (k <= kmin)
        k *= p;
      a = powmod(a, k, n);
      f = gcd_ui(a-1, n);
      q = p;
      if (f != 1)
        break;
    } END_DO_FOR_EACH_PRIME
    /* If f == n again, we could do:
     * for (savea = 3; f == n && savea < 100; savea = _XS_next_prime(savea)) {
     *   a = savea;
     *   for (q = 2; q <= B1; q = _XS_next_prime(q)) {
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
  if ( (f != 1) && (f != n) ) {
    factors[0] = f;
    factors[1] = n/f;
    MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
    return 2;
  }
  factors[0] = n;
  return 1;
}

/* Simple Williams p+1 */
static void pp1_pow(UV *cX, unsigned long exp, UV n)
{
  UV X0 = *cX;
  UV X  = *cX;
  UV Y = mulsubmod(X, X, 2, n);
  unsigned long bit;
  {
    unsigned long v = exp;
    unsigned long b = 1;
    while (v >>= 1) b++;
    bit = 1UL << (b-2);
  }
  while (bit) {
    if ( exp & bit ) {
      X = mulsubmod(X, Y, X0, n);
      Y = mulsubmod(Y, Y, 2, n);
    } else {
      Y = mulsubmod(X, Y, X0, n);
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
  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in pminus1_factor");

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

  if ( (f != 1) && (f != n) ) {
    factors[0] = f;
    factors[1] = n/f;
    MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
    return 2;
  }
  factors[0] = n;
  return 1;
}


/* My modification of Ben Buhrow's modification of Bob Silverman's SQUFOF code.
 */
static IV qqueue[100+1];
static IV qpoint;
static void enqu(IV q, IV *iter) {
  qqueue[qpoint] = q;
  if (++qpoint >= 100) *iter = -1;
}

int squfof_factor(UV n, UV *factors, UV rounds)
{
  IV rounds2 = (IV) (rounds/16);
  UV temp;
  IV iq,ll,l2,p,pnext,q,qlast,r,s,t,i;
  IV jter, iter;
  int reloop;

  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in squfof_factor");

  /* TODO:  What value of n leads to overflow? */

  qlast = 1;
  s = isqrt(n);

  p = s;
  temp = n - (s*s);                 /* temp = n - floor(sqrt(n))^2   */
  if (temp == 0) {
    factors[0] = s;
    factors[1] = s;
    return 2;
  }

  q = temp;              /* q = excess of n over next smaller square */
  ll = 1 + 2*(IV)sqrt((double)(p+p));
  l2 = ll/2;
  qpoint = 0;

  /*  In the loop below, we need to check if q is a square right before   */
  /*  the end of the loop.  Is there a faster way? The current way is     */
  /*  EXPENSIVE! (many branches and double prec sqrt)                     */

  for (jter=0; (UV)jter < rounds; jter++) {
    iq = (s + p)/q;
    pnext = iq*q - p;
    if (q <= ll) {
      if ((q & 1) == 0)
        enqu(q/2,&jter);
      else if (q <= l2)
        enqu(q,&jter);
      if (jter < 0) {
        factors[0] = n;  return 1;
      }
    }

    t = qlast + iq*(p - pnext);
    qlast = q;
    q = t;
    p = pnext;                          /* check for square; even iter   */
    if (jter & 1) continue;             /* jter is odd:omit square test  */
    r = isqrt(q);                       /* r = floor(sqrt(q))      */
    if (q != r*r) continue;
    if (qpoint == 0) break;
    qqueue[qpoint] = 0;
    reloop = 0;
    for (i=0; i<qpoint-1; i+=2) {    /* treat queue as list for simplicity*/
      if (r == qqueue[i]) { reloop = 1; break; }
      if (r == qqueue[i+1]) { reloop = 1; break; }
    }
    if (reloop || (r == qqueue[qpoint-1])) continue;
    break;
  }   /* end of main loop */

  if ((UV)jter >= rounds) {
    factors[0] = n;  return 1;
  }

  qlast = r;
  p = p + r*((s - p)/r);
  q = (n - (p*p)) / qlast;               /* q = (n - p*p)/qlast (div is exact)*/
  for (iter=0; iter<rounds2; iter++) {   /* unrolled second main loop */
    iq = (s + p)/q;
    pnext = iq*q - p;
    if (p == pnext) break;
    t = qlast + iq*(p - pnext);
    qlast = q;
    q = t;
    p = pnext;
    iq = (s + p)/q;
    pnext = iq*q - p;
    if (p == pnext) break;
    t = qlast + iq*(p - pnext);
    qlast = q;
    q = t;
    p = pnext;
    iq = (s + p)/q;
    pnext = iq*q - p;
    if (p == pnext) break;
    t = qlast + iq*(p - pnext);
    qlast = q;
    q = t;
    p = pnext;
    iq = (s + p)/q;
    pnext = iq*q - p;
    if (p == pnext) break;
    t = qlast + iq*(p - pnext);
    qlast = q;
    q = t;
    p = pnext;
  }

  if (iter >= rounds2) {
    factors[0] = n;  return 1;
  }

  if ((q & 1) == 0) q/=2;      /* q was factor or 2*factor   */

  if ( (q == 1) || ((UV)q == n) ) {
    factors[0] = n;  return 1;
  }

  p = n/q;

  /* printf(" squfof found %lu = %lu * %lu in %ld/%ld rounds\n", n, p, q, jter, iter); */

  factors[0] = p;
  factors[1] = q;
  MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
  return 2;
}

/* Another version, based on Ben Buhrow's racing SQUFOF. */

typedef struct
{
  UV mult;
  UV valid;
  UV P;
  UV bn;
  UV Qn;
  UV Q0;
  UV b0;
  UV it;
  UV imax;
} mult_t;

/* N < 2^63 (or 2^31).  *f == 1 if no factor found */
static void squfof_unit(UV n, mult_t* mult_save, UV* f)
{
  UV imax,i,Q0,b0,Qn,bn,P,bbn,Ro,S,So,t1,t2;
  int j;

  *f = 0;

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
    j = 0;
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
        *f = 0;
        return;
      }

      SQUARE_SEARCH_ITERATION;

      /* Even iteration.  Check for square: Qn = S*S */
      if (is_perfect_square( Qn, &S ))
        break;

      /* Odd iteration. */
      SQUARE_SEARCH_ITERATION;
    }
    /* printf("found square %lu after %lu iterations with mult %d\n", Qn, i, mult_save->mult); */

    /*  Reduce to G0 */
    Ro = P + S*((b0 - P)/S);
    t1 = Ro;
    So = (n - t1*t1)/S;
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
         *f = 0;
         return;
      }
    }

    *f = gcd_ui(Ro, n);
    if (*f > 1)
      return;
  }
}

#define NSQUFOF_MULT (sizeof(multipliers)/sizeof(multipliers[0]))

int racing_squfof_factor(UV n, UV *factors, UV rounds)
{
  const UV multipliers[] = {
    3*5*7*11, 3*5*7, 3*5*11, 3*5, 3*7*11, 3*7, 5*7*11, 5*7,
    3*11,     3,     5*11,   5,   7*11,   7,   11,     1   };
  const UV big2 = UV_MAX >> 2;
  mult_t mult_save[NSQUFOF_MULT];
  int still_racing;
  UV i, nn64, mult, f64;
  UV rounds_done = 0;

  /* Caller should have handled these trivial cases */
  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in racing_squfof_factor");

  /* Too big */
  if (n > big2) {
    factors[0] = n;  return 1;
  }

  for (i = 0; i < NSQUFOF_MULT; i++) {
    mult = multipliers[i];
    nn64 = n * mult;
    mult_save[i].mult = mult;
    if ((big2 / mult) < n) {
      mult_save[i].valid = 0; /* This multiplier would overflow 64-bit */
      continue;
    }
    mult_save[i].valid = 1;

    mult_save[i].b0 = isqrt(nn64);
    mult_save[i].imax = (UV) (sqrt(mult_save[i].b0) / 3);
    if (mult_save[i].imax < 20)     mult_save[i].imax = 20;
    if (mult_save[i].imax > rounds) mult_save[i].imax = rounds;

    mult_save[i].Q0 = 1;
    mult_save[i].P  = mult_save[i].b0;
    mult_save[i].Qn = nn64 - (mult_save[i].b0 * mult_save[i].b0);
    if (mult_save[i].Qn == 0) {
      factors[0] = mult_save[i].b0;
      factors[1] = n / mult_save[i].b0;
      MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
      return 2;
    }
    mult_save[i].bn = (mult_save[i].b0 + mult_save[i].P) / mult_save[i].Qn;
    mult_save[i].it = 0;
  }

  /* Process the multipliers a little at a time: 0.33*(n*mult)^1/4: 20-20k */
  do {
    still_racing = 0;
    for (i = 0; i < NSQUFOF_MULT; i++) {
      if (!mult_save[i].valid)
        continue;
      nn64 = n * (UV)multipliers[i];
      squfof_unit(nn64, &mult_save[i], &f64);
      if (f64 > 1) {
        if (f64 != multipliers[i]) {
          f64 /= gcd_ui(f64, multipliers[i]);
          if (f64 != 1) {
            factors[0] = f64;
            factors[1] = n / f64;
            MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
            return 2;
          }
        }
        /* Found trivial factor.  Quit working with this multiplier. */
        mult_save[i].valid = 0;
      }
      if (mult_save[i].valid == 1)
        still_racing = 1;
      rounds_done += mult_save[i].imax;
      if (rounds_done >= rounds)
        break;
    }
  } while (still_racing && rounds_done < rounds);

  /* No factors found */
  factors[0] = n;
  return 1;
}
