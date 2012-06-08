#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <math.h>

#include "factor.h"
#include "util.h"
#include "sieve.h"
#include "ptypes.h"

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


int trial_factor(UV n, UV *factors, UV maxtrial)
{
  UV f, m, limit;
  int nfactors = 0;

  if (maxtrial == 0)  maxtrial = UV_MAX;

  if ( (n < 2) || (maxtrial < 2) ) {
    factors[0] = n;
    return 1;
  }

  while ( (n & 1) == 0 ) {
    factors[nfactors++] = 2;
    n /= 2;
  }

  for (f = 3; (n > 1) && (f <= 7) && (f <= maxtrial); f += 2) {
    while ( (n % f) == 0 ) {
      factors[nfactors++] = f;
      n /= f;
    }
  }

  if ( (n < (7*7)) || (maxtrial < 11) ) {
    if (n != 1)
      factors[nfactors++] = n;
    return nfactors;
  }

  limit = sqrt((double) n);
  if (limit > maxtrial)
    limit = maxtrial;

  /* wheel 30 */
  f = 11;
  m = 11;
  while (f <= limit) {
    if ( (n%f) == 0 ) {
      UV newlimit;
      do {
        factors[nfactors++] = f;
        n /= f;
      } while ( (n%f) == 0 );
      newlimit = sqrt(n);
      if (newlimit < limit)  limit = newlimit;
    }
    f += wheeladvance30[m];
    m = nextwheel30[m];
  }
  if (n != 1)
    factors[nfactors++] = n;
  return nfactors;
}


#if UINT64_MAX > UV_MAX

  /* We have 64-bit available, but UV is 32-bit.  Do the math in 64-bit.
   * Even if it is emulated, it should be as fast or faster than us doing it.
   */
  #define addmod(n,a,m)  (UV)(((uint64_t)(n)+(uint64_t)(a)) % ((uint64_t)(m)))
  #define mulmod(a,b,m)  (UV)(((uint64_t)(a)*(uint64_t)(b)) % ((uint64_t)(m)))
  #define sqrmod(n,m)    (UV)(((uint64_t)(n)*(uint64_t)(n)) % ((uint64_t)(m)))

  static UV powmod(UV n, UV power, UV m) {
    UV t = 1;
    while (power) {
      if (power & 1)
        t = mulmod(t, n, m);
      n = sqrmod(n, m);
      power >>= 1;
    }
    return t;
  }

#else

  /* UV is the largest integral type available (that we know of). */

  /* if n is smaller than this, you can multiply without overflow */
  #define HALF_WORD (UVCONST(1) << (BITS_PER_WORD/2))

  static UV _mulmod(UV a, UV b, UV m) {
    UV r = 0;
    while (b > 0) {
      if (b & 1) {
        if (r == 0) {
          r = a;
        } else {
          r = m - r;
          r = (a >= r)  ?  a-r  :  m-r+a;
        }
      }
      a = (a > (m-a))  ?  (a-m)+a  :  a+a;
      b >>= 1;
    }
    return r;
  }

  #define addmod(n,a,m) ((((m)-(n)) > (a))  ?  ((n)+(a))  :  ((n)+(a)-(m)))
  #define mulmod(a,b,m) (((a)|(b)) < HALF_WORD) ? ((a)*(b))%(m) : _mulmod(a,b,m)
  #define sqrmod(n,m)   ((n) < HALF_WORD)       ? ((n)*(n))%(m) : _mulmod(n,n,m)

  /* n^power mod m */
  static UV powmod(UV n, UV power, UV m) {
    UV t = 1;
    if (m < HALF_WORD) {
      n %= m;
      while (power) {
        if (power & 1)
          t = (t*n)%m;
        n = (n*n)%m;
        power >>= 1;
      }
    } else {
      while (power) {
        if (power & 1)
          t = mulmod(t, n, m);
        n = sqrmod(n,m);
        power >>= 1;
      }
    }
    return t;
  }

#endif

/* n^power + a mod m */
#define powaddmod(n, p, a, m)  addmod(powmod(n,p,m),a,m)

/* n^2 + a mod m */
#define sqraddmod(n, a, m)     addmod(sqrmod(n,m),  a,m)


/* Miller-Rabin probabilistic primality test
 * Returns 1 if probably prime relative to the bases, 0 if composite.
 * Bases must be between 2 and n-2
 */
int miller_rabin(UV n, const UV *bases, UV nbases)
{
  int b;
  int s = 0;
  UV d = n-1;

  MPUassert(n > 3, "MR called with n <= 3");

  while ( (d&1) == 0 ) {
    s++;
    d >>= 1;
  }
  for (b = 0; b < nbases; b++) {
    int r;
    UV a = bases[b];
    UV x;

    /* Skip invalid bases */
    if ( (a < 2) || (a > (n-2)) )
      croak("Base %"UVuf" is invalid for input %"UVuf, a, n);

    x = powmod(a, d, n);
    if ( (x == 1) || (x == (n-1)) )  continue;

    for (r = 0; r < s; r++) {
      x = sqrmod(x, n);
      if (x == 1) {
        return 0;
      } else if (x == (n-1)) {
        a = 0;
        break;
      }
    }
    if (a != 0)
      return 0;
  }
  return 1;
}

int is_prob_prime(UV n)
{
  UV bases[12];
  int nbases;
  int prob_prime;

  if ( (n == 2) || (n == 3) || (n == 5) || (n == 7) )
    return 2;
  if ( (n<2) || ((n% 2)==0) || ((n% 3)==0) || ((n% 5)==0) || ((n% 7)==0) )
    return 0;
  if (n < 121)
    return 2;

#if BITS_PER_WORD == 32
  if (n < UVCONST(9080191)) {
    bases[0] = 31; bases[1] = 73; nbases = 2;
  } else  {
    bases[0] = 2; bases[1] = 7; bases[2] = 61; nbases = 3;
  }
#else
#if 1
  /* Better basis from:  http://miller-rabin.appspot.com/ */
  if (n < UVCONST(9080191)) {
    bases[0] = 31;
    bases[1] = 73;
    nbases = 2;
  } else if (n < UVCONST(4759123141)) {
    bases[0] = 2;
    bases[1] = 7;
    bases[2] = 61;
    nbases = 3;
  } else if (n < UVCONST(105936894253)) {
    bases[0] = 2;
    bases[1] = UVCONST( 1005905886 );
    bases[2] = UVCONST( 1340600841 );
    nbases = 3;
  } else if (n < UVCONST(31858317218647)) {
    bases[0] = 2;
    bases[1] = UVCONST( 642735     );
    bases[2] = UVCONST( 553174392  );
    bases[3] = UVCONST( 3046413974 );
    nbases = 4;
  } else if (n < UVCONST(3071837692357849)) {
    bases[0] = 2;
    bases[1] = UVCONST( 75088      );
    bases[2] = UVCONST( 642735     );
    bases[3] = UVCONST( 203659041  );
    bases[4] = UVCONST( 3613982119 );
    nbases = 5;
  } else {
    bases[0] = 2;
    bases[1] = UVCONST( 325        );
    bases[2] = UVCONST( 9375       );
    bases[3] = UVCONST( 28178      );
    bases[4] = UVCONST( 450775     );
    bases[5] = UVCONST( 9780504    );
    bases[6] = UVCONST( 1795265022 );
    nbases = 7;
  }
#else
  /* More standard bases */
  if (n < UVCONST(9080191)) {
    bases[0] = 31; bases[1] = 73; nbases = 2;
  } else if (n < UVCONST(4759123141)) {
    bases[0] = 2; bases[1] = 7; bases[2] = 61; nbases = 3;
  } else if (n < UVCONST(21652684502221)) {
    bases[0] = 2; bases[1] = 1215; bases[2] = 34862; bases[3] = 574237825;
    nbases = 4;
  } else if (n < UVCONST(341550071728321)) {
    bases[0] =  2; bases[1] =  3; bases[2] =  5; bases[3] =  7; bases[4] = 11;
    bases[5] = 13; bases[6] = 17; nbases = 7;
  } else if (n < UVCONST(3825123056546413051)) {
    bases[0] =  2; bases[1] =  3; bases[2] =  5; bases[3] =  7; bases[4] = 11;
    bases[5] = 13; bases[6] = 17; bases[7] = 19; bases[8] = 23; nbases = 9;
  } else {
    bases[0] =  2; bases[1] =  3; bases[2] =  5; bases[3] =  7; bases[4] = 11;
    bases[5] = 13; bases[6] = 17; bases[7] = 19; bases[8] = 23; bases[9] = 29;
    bases[10]= 31; bases[11]= 37;
    nbases = 12;
  }
#endif
#endif
  prob_prime = miller_rabin(n, bases, nbases);
  return 2*prob_prime;
}



/* Knuth volume 2, algorithm C.
 * Very fast for small numbers, grows rapidly.
 * SQUFOF is better for numbers nearing the 64-bit limit.
 */
int fermat_factor(UV n, UV *factors, UV rounds)
{
  IV sqn, x, y, r;

  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in fermat_factor");

  sqn = sqrt((double) n);
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
  if ( (r != 1) && (r != n) ) {
    factors[0] = r;
    factors[1] = n/r;
    MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
    return 2;
  }
  factors[0] = n;
  return 1;
}

/* Hart's One Line Factorization.
 * Translation from my GMP, missing perfect square calc and premult.
 */
int holf_factor(UV n, UV *factors, UV rounds)
{
  UV i, s, m, f;

  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in holf_factor");

  for (i = 1; i <= rounds; i++) {
    s = sqrt(n*i);                      /* TODO: overflow here */
    if ( (s*s) != (n*i) )  s++;
    m = sqrmod(s, n);
    /* Cheaper would be:
     *     if (m is probably not a perfect sequare)  continue;
     *     f = sqrt(m);
     *     if (f*f == m) { yay }
     */
    if ( ((m&2)!= 0) || ((m&7)==5) || ((m&11) == 8) )
      continue;
    f = sqrt(m);                        /* expensive */
    if ( (f*f) == m ) {
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


/* Pollard / Brent
 *
 * Probabilistic.  If you give this a prime number, it will loop
 * until it runs out of rounds.
 */
int pbrent_factor(UV n, UV *factors, UV rounds)
{
  UV a, f, i;
  UV Xi = 2;
  UV Xm = 2;

  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in pbrent_factor");

  switch (n%4) {
    case 0:  a =  1; break;
    case 1:  a =  3; break;
    case 2:  a =  5; break;
    case 3:  a =  7; break;
    default: a = 11; break;
  }

  for (i = 1; i <= rounds; i++) {
    Xi = sqraddmod(Xi, a, n);
    f = gcd_ui( (Xi>Xm) ? Xi-Xm : Xm-Xi, n);
    if ( (f != 1) && (f != n) ) {
      factors[0] = f;
      factors[1] = n/f;
      MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
      return 2;
    }
    if ( (i & (i-1)) == 0)   /* i is a power of 2 */
      Xm = Xi;
  }
  factors[0] = n;
  return 1;
}

/* Pollard's Rho
 *
 * Probabilistic.  If you give this a prime number, it will loop
 * until it runs out of rounds.
 */
int prho_factor(UV n, UV *factors, UV rounds)
{
  int in_loop = 0;
  UV a, f, i;
  UV U = 7;
  UV V = 7;

  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in prho_factor");

  switch (n%4) {
    case 0:  a =  5; break;
    case 1:  a =  7; break;
    case 2:  a = 11; break;
    case 3:  a =  1; break;
    default: a =  3; break;
  }

  for (i = 1; i <= rounds; i++) {
    U = sqraddmod(U, a, n);
    V = sqraddmod(V, a, n);
    V = sqraddmod(V, a, n);
    f = gcd_ui( (U > V) ? U-V : V-U, n);
    if (f == n) {
      if (in_loop++)     /* Mark that we've been here */
        break;           /* Exit now if we're cycling */
    } else if (f != 1) {
      factors[0] = f;
      factors[1] = n/f;
      MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
      return 2;
    }
  }
  factors[0] = n;
  return 1;
}

/* Pollard's P-1
 *
 * Probabilistic.  If you give this a prime number, it will loop
 * until it runs out of rounds.
 */
int pminus1_factor(UV n, UV *factors, UV rounds)
{
  UV f, i;
  UV kf = 13;

  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in pminus1_factor");

  for (i = 1; i <= rounds; i++) {
    kf = powmod(kf, i, n);
    if (kf == 0) kf = n;
    f = gcd_ui(kf-1, n);
    if ( (f != 1) && (f != n) ) {
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
  UV rounds2 = rounds/16;
  UV temp;
  IV iq,ll,l2,p,pnext,q,qlast,r,s,t,i;
  IV jter, iter;
  int reloop;

  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in squfof_factor");

  /* TODO:  What value of n leads to overflow? */

  qlast = 1;
  s = sqrt(n);

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

  for (jter=0; jter < rounds; jter++) {
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
    r = (int)sqrt((double)q);                 /* r = floor(sqrt(q))      */
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

  if (jter >= rounds) {
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

  if ( (q == 1) || (q == n) ) {
    factors[0] = n;  return 1;
  }

  p = n/q;

  /* printf(" squfof found %lu = %lu * %lu in %ld/%ld rounds\n", n, p, q, jter, iter); */

  factors[0] = p;
  factors[1] = q;
  MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
  return 2;
}


/* TODO: Add Jason Papadopoulos's racing SQUFOF */
