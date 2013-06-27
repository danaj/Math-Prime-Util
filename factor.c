#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define FUNC_gcd_ui
#include "ptypes.h"
#include "factor.h"
#include "util.h"
#include "sieve.h"
#include "mulmod.h"
#include "cache.h"

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
  if ( (n < 4) || (maxtrial < 4) ) {
    factors[0] = n;
    return 1;
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
    f = primes_small[sp];
    while (f <= limit && f <= 2003) {
      if ( (n%f) == 0 ) {
        do {
          factors[nfactors++] = f;
          n /= f;
        } while ( (n%f) == 0 );
        newlimit = isqrt(n);
        if (newlimit < limit)  limit = newlimit;
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

static int jacobi_iu(IV in, UV m) {
  int j = 1;
  UV n = (in < 0) ? -in : in;

  if (m <= 0 || (m%2) == 0) return 0;
  if (in < 0 && (m%4) == 3) j = -j;
  while (n != 0) {
    while ((n % 2) == 0) {
      n >>= 1;
      if ( (m % 8) == 3 || (m % 8) == 5 )  j = -j;
    }
    { UV t = n; n = m; m = t; }
    if ( (n % 4) == 3 && (m % 4) == 3 )  j = -j;
    n = n % m;
  }
  return (m == 1) ? j : 0;
}


/* Fermat pseudoprime */
int _XS_is_pseudoprime(UV n, UV a)
{
  UV x;
  UV const nm1 = n-1;

  if (n == 2 || n == 3)  return 1;
  if (n < 5) return 0;
  if (a < 2) croak("Base %"UVuf" is invalid", a);
  if (a >= n) {
    a %= n;
    if ( a <= 1 || a == nm1 )
      return 1;
  }

  x = powmod(a, nm1, n);    /* x = a^(n-1) mod n */
  return (x == 1);
}


/* Miller-Rabin probabilistic primality test
 * Returns 1 if probably prime relative to the bases, 0 if composite.
 * Bases must be between 2 and n-2
 */
int _XS_miller_rabin(UV n, const UV *bases, int nbases)
{
  UV const nm1 = n-1;
  UV d = n-1;
  int b, r, s = 0;

  MPUassert(n > 3, "MR called with n <= 3");

  while ( (d&1) == 0 ) {
    s++;
    d >>= 1;
  }
  for (b = 0; b < nbases; b++) {
    UV x, a = bases[b];

    if (a < 2)
      croak("Base %"UVuf" is invalid", a);
    if (a >= n)
      a %= n;
    if ( (a <= 1) || (a == nm1) )
      continue;

    /* n is a strong pseudoprime to this base if either
     *   -  a^d = 1 mod n
     *   -  a^(d2^r) = -1 mod n for some r: 0 <= r <= s-1
     */

    x = powmod(a, d, n);
    if ( (x == 1) || (x == nm1) )  continue;

    /* cover r = 1 to s-1, r=0 was just done */
    for (r = 1; r < s; r++) {
      x = sqrmod(x, n);
      if ( x == nm1 )  break;
      if ( x == 1   )  return 0;
    }
    if (r >= s)
      return 0;
  }
  return 1;
}

/* M-R with a = 2 and some checks removed.  For internal use. */
int _SPRP2(UV n)
{
  UV const nm1 = n-1;
  UV d = n-1;
  UV x;
  int r, s = 0;

  MPUassert(n > 3, "S-PRP-2 called with n <= 3");
  if (!(n & 1)) return 0;
  while ( (d & 1) == 0 ) {  s++;  d >>= 1; }
  /* n is a strong pseudoprime to this base if either
   *   -  a^d = 1 mod n
   *   -  a^(d2^r) = -1 mod n for some r: 0 <= r <= s-1 */
  x = powmod(2, d, n);
  if (x == 1 || x == nm1)  return 1;

  /* just did r=0, now test r = 1 to s-1 */
  for (r = 1; r < s; r++) {
    x = sqrmod(x, n);
    if (x == nm1)  return 1;
  }
  return 0;
}

int _XS_is_prob_prime(UV n)
{
  UV bases[12];
  int nbases;
  int prob_prime;

  if (n == 2 || n == 3 || n == 5)           return 2;
  if (n < 2 || !(n%2) || !(n%3) || !(n%5))  return 0;
  if (n <   49) /* 7*7 */                   return 2;
  if (!(n% 7) || !(n%11) || !(n%13) || !(n%17) || !(n%19) || !(n%23)) return 0;
  if (!(n%29) || !(n%31) || !(n%37) || !(n%41) || !(n%43) || !(n%47)) return 0;
  if (n < 2809) /* 53*53 */                 return 2;

#if BITS_PER_WORD == 32

  /* These aren't ideal.  Could use 1 when n < 49191, 2 when n < 360018361 */
  if (n < UVCONST(9080191)) {
    bases[0] = 31; bases[1] = 73; nbases = 2;
  } else  {
    bases[0] = 2; bases[1] = 7; bases[2] = 61; nbases = 3;
  }
  prob_prime = _XS_miller_rabin(n, bases, nbases);
  return 2*prob_prime;

#else

  /* Verified with Feitsma database.  No counterexamples below 2^64.
   * This is faster than multiple M-R routines once we're over 32-bit */
  if (n >= UVCONST(4294967295)) {
    prob_prime = _SPRP2(n) && _XS_is_lucas_pseudoprime(n,2);
    return 2*prob_prime;
  }

  /* Better bases from http://miller-rabin.appspot.com/, 27 May 2013 */
  /* Verify with:
   *  cat /local/share/spsps-below-2-to-64.txt | perl -MMath::Prime::Util=:all
   *    -nE 'chomp; next unless is_strong_pseudoprime($_, @bases); say;'
   */
  if (n < UVCONST(341531)) {
    bases[0] = UVCONST(  9345883071009581737 );
    nbases = 1;
  } else if (n < UVCONST(885594169)) {
    bases[0] = UVCONST(   725270293939359937 );
    bases[1] = UVCONST(  3569819667048198375 );
    nbases = 2;
  } else if (n < UVCONST(350269456337)) {
    bases[0] = UVCONST(  4230279247111683200 );
    bases[1] = UVCONST( 14694767155120705706 );
    bases[2] = UVCONST( 16641139526367750375 );
    nbases = 3;
  } else if (n < UVCONST(55245642489451)) {      /* 38+ bits */
    bases[0] = 2;
    bases[1] = UVCONST(      141889084524735 );
    bases[2] = UVCONST(  1199124725622454117 );
    bases[3] = UVCONST( 11096072698276303650 );
    nbases = 4;
  } else if (n < UVCONST(7999252175582851)) {    /* 45+ bits */
    bases[0] = 2;
    bases[1] = UVCONST(        4130806001517 );
    bases[2] = UVCONST(   149795463772692060 );
    bases[3] = UVCONST(   186635894390467037 );
    bases[4] = UVCONST(  3967304179347715805 );
    nbases = 5;
  } else if (n < UVCONST(585226005592931977)) {  /* 52+ bits */
    bases[0] = 2;
    bases[1] = UVCONST(      123635709730000 );
    bases[2] = UVCONST(     9233062284813009 );
    bases[3] = UVCONST(    43835965440333360 );
    bases[4] = UVCONST(   761179012939631437 );
    bases[5] = UVCONST(  1263739024124850375 );
    nbases = 6;
  } else {                                       /* 59+ bits */
    bases[0] = 2;
    bases[1] = UVCONST( 325        );
    bases[2] = UVCONST( 9375       );
    bases[3] = UVCONST( 28178      );
    bases[4] = UVCONST( 450775     );
    bases[5] = UVCONST( 9780504    );
    bases[6] = UVCONST( 1795265022 );
    nbases = 7;
  }
  prob_prime = _XS_miller_rabin(n, bases, nbases);
  return 2*prob_prime;
#endif
}

/* Generic Lucas sequence for any appropriate P and Q */
void lucas_seq(UV* Uret, UV* Vret, UV* Qkret, UV n, IV P, IV Q, UV k)
{
  UV U, V, b, Dmod, Qmod, Pmod, Qk;
  IV D;

  if (k == 0) {
    *Uret = 0;
    *Vret = 2;
    *Qkret = Q;
    return;
  }

  Qmod = (Q < 0) ? Q + n : Q;
  Pmod = (P < 0) ? P + n : P;
  Dmod = submod( mulmod(Pmod, Pmod, n), mulmod(4, Qmod, n), n);
  MPUassert(Dmod != 0, "lucas_seq: D is 0");
  U = 1;
  V = Pmod;
  Qk = Qmod;
  { UV v = k; b = 1; while (v >>= 1) b++; }

  if (Q == 1) {
    while (b > 1) {
      U = mulmod(U, V, n);
      V = mulsubmod(V, V, 2, n);
      b--;
      if ( (k >> (b-1)) & UVCONST(1) ) {
        UV t2 = mulmod(U, Dmod, n);
        U = muladdmod(U, Pmod, V, n);
        if (U & 1) { U = (n>>1) + (U>>1) + 1; } else { U >>= 1; }
        V = muladdmod(V, P, t2, n);
        if (V & 1) { V = (n>>1) + (V>>1) + 1; } else { V >>= 1; }
      }
    }
  } else if (P == 1 && Q == -1) {
    /* This is about 30% faster than the generic code below.  Since 50% of
     * Lucas and strong Lucas tests come here, I think it's worth doing. */
    int sign = Q;
    while (b > 1) {
      U = mulmod(U, V, n);
      if (sign == 1) V = mulsubmod(V, V, 2, n);
      else           V = muladdmod(V, V, 2, n);
      sign = 1;   /* Qk *= Qk */
      b--;
      if ( (k >> (b-1)) & UVCONST(1) ) {
        UV t2 = mulmod(U, Dmod, n);
        U = addmod(U, V, n);
        if (U & 1) { U = (n>>1) + (U>>1) + 1; } else { U >>= 1; }
        V = addmod(V, t2, n);
        if (V & 1) { V = (n>>1) + (V>>1) + 1; } else { V >>= 1; }
        sign = -1;  /* Qk *= Q */
      }
    }
    if (sign == 1) Qk = 1;
  } else {
    while (b > 1) {
      U = mulmod(U, V, n);
      V = mulsubmod(V, V, addmod(Qk,Qk,n), n);
      Qk = sqrmod(Qk, n);
      b--;
      if ( (k >> (b-1)) & UVCONST(1) ) {
        UV t2 = mulmod(U, Dmod, n);
        U = muladdmod(U, Pmod, V, n);
        if (U & 1) { U = (n>>1) + (U>>1) + 1; } else { U >>= 1; }
        V = muladdmod(V, P, t2, n);
        if (V & 1) { V = (n>>1) + (V>>1) + 1; } else { V >>= 1; }
        Qk = mulmod(Qk, Qmod, n);
      }
    }
  }
  *Uret = U;
  *Vret = V;
  *Qkret = Qk;
}

/* Lucas tests:
 *  0: Standard
 *  1: Strong
 *  2: Extra Strong (Mo/Jones/Grantham)
 *
 * Goal:
 *       (1) no false results when combined with the SPRP-2 test.
 *       (2) fast enough to use SPRP-2 + this in place of 3+ M-R tests.
 *
 * For internal purposes, we typically want to use the extra strong test
 * because it is slightly faster (Q = 1 simplies things).  None of them have
 * any false positives for the BPSW test.
 *
 * This runs 4-7x faster than MPU::GMP, which means ~10x faster than most GMP
 * implementations.  It is about 2x slower than a single M-R test.
 */
int _XS_is_lucas_pseudoprime(UV n, int strength)
{
  IV P, Q, D;
  UV U, V, Qk, d, s;

  if (n == 2 || n == 3) return 1;
  if (n < 5 || (n%2) == 0) return 0;
  if (n == UV_MAX) return 0;

  if (strength < 2) {
    UV Du = 5;
    IV sign = 1;
    while (1) {
      D = Du * sign;
      if (gcd_ui(Du, n) > 1 && gcd_ui(Du, n) != n) return 0;
      if (jacobi_iu(D, n) == -1)
        break;
      if (Du == 21 && is_perfect_square(n, 0)) return 0;
      Du += 2;
      sign = -sign;
    }
    P = 1;
    Q = (1 - D) / 4;
  } else {
    P = 3;
    Q = 1;
    while (1) {
      D = P*P - 4;
      if (gcd_ui(D, n) > 1 && gcd_ui(D, n) != n) return 0;
      if (jacobi_iu(D, n) == -1)
        break;
      if (P == 21 && is_perfect_square(n, 0)) return 0;
      P++;
    }
  }
  MPUassert( D == (P*P - 4*Q) , "is_lucas_pseudoprime: incorrect DPQ");

  d = n+1;
  s = 0;
  if (strength > 0)
    while ( (d & 1) == 0 ) {  s++;  d >>= 1; }
  lucas_seq(&U, &V, &Qk, n, P, Q, d);

  if (strength == 0) {
    if (U == 0)
      return 1;
  } else if (strength == 1) {
    if (U == 0 || V == 0)
      return 1;
    while (s--) {
      V = muladdmod(V, V, n - addmod(Qk, Qk, n), n);
      if (V == 0)
        return 1;
      if (s)
        Qk = sqrmod(Qk, n);
    }
  } else {
    if ( (U == 0 && (V == 2 || V == (n-2))) || (V == 0) )
      return 1;
    while (s--) {
      V = mulsubmod(V, V, 2, n);
      if (V == 0)
        return 1;
    }
  }
  return 0;
}


UV _XS_divisor_sum(UV n)
{
  UV factors[MPU_MAX_FACTORS+1];
  int nfac, i;
  UV product = 1;

  if (n <= 1) return n;
  nfac = factor(n, factors);
  for (i = 0; i < nfac; i++) {
    if (i+1 < nfac && factors[i] == factors[i+1]) {
      UV fmult = factors[i]*factors[i];
      do {
        fmult *= factors[i++];
      } while (i+1 < nfac && factors[i] == factors[i+1]);
      product *= (fmult-1) / (factors[i]-1);
    } else {
      product *= factors[i]+1;
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
int pplus1_factor(UV n, UV *factors, UV B)
{
  UV i, f, b;
  UV P1, P2, P3, U1, U2, U3, V1, V2, V3;

  /* Calculate 3 sequences at once */
  P1 =  5 % n;
  P2 =  9 % n;
  P3 = 13 % n;
  U1 = P1;  U2 = P2;  U3 = P3;
  for (i = 1; i < B; i++) {
    { UV v = i; b = 1; while (v >>= 1) b++; }
    V1 = mulsubmod(P1, P1, 2, n);
    V2 = mulsubmod(P2, P2, 2, n);
    V3 = mulsubmod(P3, P3, 2, n);
    while (b > 1) {
      b--;
      if ( (i >> (b-1)) & UVCONST(1) ) {
        U1 = mulsubmod(U1, V1, P1, n);
        V1 = mulsubmod(V1, V1,  2, n);
        U2 = mulsubmod(U2, V2, P2, n);
        V2 = mulsubmod(V2, V2,  2, n);
        U3 = mulsubmod(U3, V3, P3, n);
        V3 = mulsubmod(V3, V3,  2, n);
      } else {
        V1 = mulsubmod(U1, V1, P1, n);
        U1 = mulsubmod(U1, U1,  2, n);
        V2 = mulsubmod(U2, V2, P2, n);
        U2 = mulsubmod(U2, U2,  2, n);
        V3 = mulsubmod(U3, V3, P3, n);
        U3 = mulsubmod(U3, U3,  2, n);
      }
    }
    P1 = U1;  P2 = U2;  P3 = U3;
    f = gcd_ui(n, submod(U1, 2, n));
    if (f == 1 || f == n)
      f = gcd_ui(n, submod(U2, 2, n));
    if (f == 1 || f == n)
      f = gcd_ui(n, submod(U3, 2, n));
    if (f > 1 && f != n) {
      factors[0] = f;
      factors[1] = n/f;
      MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
      return 2;
    }
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
  q = (n - (p*p)) / qlast;			/* q = (n - p*p)/qlast (div is exact)*/
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

// N < 2^63 (or 2^31).  *f == 1 if no factor found
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
    // i is now even
    while (1) {
      // We need to know P, bn, Qn, Q0, iteration count, i  from prev
      if (i >= imax) {
        // save state and try another multiplier.
        mult_save->P = P;
        mult_save->bn = bn;
        mult_save->Qn = Qn;
        mult_save->Q0 = Q0;
        mult_save->it = i;
        *f = 0;
        return;
      }

      SQUARE_SEARCH_ITERATION;

      // Even iteration.  Check for square: Qn = S*S
      if (is_perfect_square( Qn, &S ))
        break;

      // Odd iteration.
      SQUARE_SEARCH_ITERATION;
    }
    /* printf("found square %lu after %lu iterations with mult %d\n", Qn, i, mult_save->mult); */

    // Reduce to G0
    Ro = P + S*((b0 - P)/S);
    t1 = Ro;
    So = (n - t1*t1)/S;
    bbn = (b0+Ro)/So;

    // Search for symmetry point
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


/****************************************************************************/

/*
 *
 * The Frobenius-Underwood test has no known counterexamples below 10^13, but
 * has not been extensively tested above that.  This is the Minimal Lambda+2
 * test from section 9 of "Quadratic Composite Tests" by Paul Underwood.
 *
 * Given the script:
 *  time mpu 'forprimes { Math::Prime::Util::_XS_is_frobenius_underwood_pseudoprime($_); Math::Prime::Util::_XS_is_frobenius_underwood_pseudoprime($_+2); } 100_000_000'
 * and replacing the tests appropriately, I get these times:
 *
 *   0.57    $_ (cost of empty loop)
 *   6.89    _XS_is_pseudoprime($_,2)
 *   6.82    _XS_miller_rabin($_,2)
 *  17.42    _XS_is_lucas_pseudoprime($_,0) (standard)
 *  17.15    _XS_is_lucas_pseudoprime($_,1) (strong)
 *  11.81    _XS_is_lucas_pseudoprime($_,2) (extra strong)
 *  13.07    _XS_is_frobenius_underwood_pseudoprime($_)
 *   7.87    _XS_is_prob_prime($_)
 *   8.74    _XS_is_prime($_)
 *
 * At these sizes is_prob_prime is doing 1-2 M-R tests.  The input validation
 * is adding a noticeable overhead to is_prime.
 *
 * With a set of 10k 64-bit random primes; 'do { die unless ... } for 1..500'
 *
 *   0.36    empty loop
 *  12.38    _XS_is_pseudoprime($_,2)
 *  12.05    _XS_miller_rabin($_,2)
 *  24.95    _XS_is_extra_strong_lucas_pseudoprime($_)
 *  22.35    _XS_is_frobenius_underwood_pseudoprime($_)
 *  36.67    _XS_is_prob_prime($_)
 *  37.24    _XS_is_prime($_)
 *
 * At this point is_prob_prime has transitioned to BPSW.
 *
 * Calling a powmod a 'Selfridge' unit, then we see:
 *    1 Selfridge unit    M-R test
 *    2 Selfridge units   Lucas or Frobenius-Underwood
 *    3 Selfridge units   BPSW
 *
 * We try to structure the primality test like:
 *   1) simple divisibility    very fast       primes and ~10% of composites
 *   2) M-R with base 2        1 Selfridge     primes and .00000000002% comps
 *   3) Lucas test             2 Selfridge     only primes
 *
 * Hence given a random composite, about 90% of the time it costs us almost
 * nothing.  After spending 1 Selfridge on the first MR test, less than 32M
 * composites remain undecided out of 18 quintillion 64-bit composites.  The
 * final Lucas test has no false positives.
 * Replacing the Lucas test with the F-U test won't save any time.  Replacing
 * the whole thing with the F-U test (assuming it has no false results for
 * all 64-bit values, which hasn't been verified), doesn't help either.
 * It's 2/3 the cost for primes, but much more expensive for composites.  It
 * seems of interest for > 2^64 as a different test to do in addition to BPSW.
 */


int _XS_is_frobenius_underwood_pseudoprime(UV n)
{
  int bit;
  UV x, result, multiplier, a, b, np1, len, t1, t2, na;
  IV t;

  if (n < 2) return 0;
  if (n < 4) return 1;
  if ((n % 2) == 0) return 0;
  if (is_perfect_square(n,0)) return 0;
  if (n == UV_MAX) return 0;

  x = 0;
  t = -1;
  while ( jacobi_iu( t, n ) != -1 ) {
    x++;
    t = (IV)(x*x) - 4;
  }
  result = addmod( addmod(x, x, n), 5, n);
  multiplier = addmod(x, 2, n);

  a = 1;
  b = 2;
  np1 = n+1;
  { UV v = np1; len = 1;  while (v >>= 1) len++; }

  if (x == 0) {
    for (bit = len-2; bit >= 0; bit--) {
      t2 = addmod(b, b, n);
      na = mulmod(a, t2, n);
      t1 = addmod(b, a, n);
      t2 = submod(b, a, n);
      b = mulmod(t1, t2, n);
      a = na;
      if ( (np1 >> bit) & UVCONST(1) ) {
        t1 = mulmod(a, 2, n);
        na = addmod(t1, b, n);
        t1 = addmod(b, b, n);
        b = submod(t1, a, n);
        a = na;
      }
    }
  } else {
    for (bit = len-2; bit >= 0; bit--) {
      t1 = mulmod(a, x, n);
      t2 = addmod(b, b, n);
      t1 = addmod(t1, t2, n);
      na = mulmod(a, t1, n);
      t1 = addmod(b, a, n);
      t2 = submod(b, a, n);
      b = mulmod(t1, t2, n);
      a = na;
      if ( (np1 >> bit) & UVCONST(1) ) {
        t1 = mulmod(a, multiplier, n);
        na = addmod(t1, b, n);
        t1 = addmod(b, b, n);
        b = submod(t1, a, n);
        a = na;
      }
    }
  }
  if (_XS_get_verbose()>1) printf("%"UVuf" is %s with x = %"UVuf"\n", n, (a == 0 && b == result) ? "probably prime" : "composite", x);
  if (a == 0 && b == result)
    return 1;
  return 0;
}
