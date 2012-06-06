#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>

#include "factor.h"
#include "util.h"
#include "sieve.h"

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

static UV gcd_ui(UV x, UV y) {
  UV t;

  if (y < x) { t = x; x = y; y = t; }

  while (y > 0) {
    x = x % y;
    t = x; x = y; y = t;
  }
  return x;
}

/* n^power + a mod m */
/* This has serious overflow issues, making the programs that use it dubious */
static UV powmod(UV n, UV power, UV add, UV m) {
  UV t = 1;
  n = n % m;
  /* t and n will always be < m from now on */
  while (power) {
    if (power & 1)
      t = (t * n) % m;
    n = (n * n) % m;
    power >>= 1;
  }
  /* (t+a) % m, noting t is always < m */
  return ( ((m-t) > add) ? (t+add) : (t+add-m) );
}

/* Knuth volume 2, algorithm C.
 * Very fast for small numbers, grows rapidly.
 * SQUFOF is better for numbers nearing the 64-bit limit.
 */
int fermat_factor(UV n, UV *factors)
{
  int nfactors = 0;
  IV sqn, x, y, r;

  if (n < 2) {
    factors[0] = n;
    return 1;
  }

  while ((n & 1) == 0) {
    factors[nfactors++] = 2;
    n /= 2;
  }

  if (n == 1)
    return nfactors;

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
  if (r != 1)
    factors[nfactors++] = r;
  n /= r;
  if (n != 1)
    factors[nfactors++] = n;
  return nfactors;
}

/* Pollard / Brent
 *
 * Probabilistic.  If you give this a prime number, it will loop
 * until it runs out of rounds.
 */
int pbrent_factor(UV n, UV *factors, UV rounds)
{
  int nfactors = 0;
  UV a, f, Xi, Xm, i;

  if (n < 2) {
    factors[0] = n;
    return 1;
  }

  while ((n & 1) == 0) {
    factors[nfactors++] = 2;
    n /= 2;
  }

  if (n == 1)
    return nfactors;

  Xi = 2;
  Xm = 2;
  switch (n%4) {
    case 0:  a =  1; break;
    case 1:  a =  3; break;
    case 2:  a =  5; break;
    case 3:  a =  7; break;
    default: a = 11; break;
  }

  for (i = 1; i < rounds; i++) {
    Xi = powmod(Xi, 2, a, n);
    f = gcd_ui(Xi - Xm, n);
    if ( (f != 1) && (f != n) ) {
      factors[nfactors++] = f;
      factors[nfactors++] = n/f;
      return nfactors;
    }
    if ( (i & (i-1)) == 0)   /* i is a power of 2 */
      Xm = Xi;
  }
  factors[nfactors++] = n;
  return nfactors;
}

/* Pollard's Rho
 *
 * Probabilistic.  If you give this a prime number, it will loop
 * until it runs out of rounds.
 */
int prho_factor(UV n, UV *factors, UV rounds)
{
  int nfactors = 0;
  UV a, f, t, U, V, i;

  if (n < 2) {
    factors[0] = n;
    return 1;
  }

  while ((n & 1) == 0) {
    factors[nfactors++] = 2;
    n /= 2;
  }

  if (n == 1)
    return nfactors;

  switch (n%4) {
    case 0:  a =  5; break;
    case 1:  a =  7; break;
    case 2:  a = 11; break;
    case 3:  a =  1; break;
    default: a =  3; break;
  }

  U = 7;
  V = 7;

  for (i = 1; i < rounds; i++) {
    U = powmod(U, 2, a, n);
    V = powmod(V, 2, a, n);
    V = powmod(V, 2, a, n);

    f = gcd_ui( (U > V) ? U-V : V-U, n);
    if ( (f != 1) && (f != n) ) {
      factors[nfactors++] = f;
      factors[nfactors++] = n/f;
      return nfactors;
    }
  }
  factors[nfactors++] = n;
  return nfactors;
}

/* Pollard's P-1
 *
 * Probabilistic.  If you give this a prime number, it will loop
 * until it runs out of rounds.
 */
int pminus1_factor(UV n, UV *factors, UV rounds)
{
  int nfactors = 0;
  UV f, b, i;

  if (n < 2) {
    factors[0] = n;
    return 1;
  }

  while ((n & 1) == 0) {
    factors[nfactors++] = 2;
    n /= 2;
  }

  if (n == 1)
    return nfactors;

  b = 13;

  for (i = 1; i < rounds; i++) {
    b = powmod(b+1, i, 0, n);
    if (b == 0)  b = n;
    b--;
    f = gcd_ui(b, n);
    if (f == n) {
      factors[nfactors++] = n;
      return nfactors;
    } else if (f != 1) {
      factors[nfactors++] = f;
      factors[nfactors++] = n/f;
      return nfactors;
    }
  }
  factors[nfactors++] = n;
  return nfactors;
}

/* My modification of Ben Buhrow's modification of Bob Silverman's SQUFOF code.
 * I like Jason P's code a lot -- I should put it in. */
static IV qqueue[100+1];
static IV qpoint;
static void enqu(IV q, IV *iter) {
  qqueue[qpoint] = q;
  if (++qpoint >= 100) *iter = -1;
}

int squfof_factor(UV n, UV *factors, UV rounds, UV refactor_above)
{
  int nfactors = 0;
  UV rounds2 = rounds/16;
  UV temp;
  IV iq,ll,l2,p,pnext,q,qlast,r,s,t,i;
  IV jter, iter;
  int reloop;

  if ( (n < 2) ) {
    factors[0] = n;
    return 1;
  }

  while ((n & 1) == 0) {
    factors[nfactors++] = 2;
    n /= 2;
  }

  if (n == 1)
    return nfactors;

  /* TODO:  What value of n leads to overflow? */

  qlast = 1;
  s = sqrt(n);

  p = s;
  temp = n - (s*s);                 /* temp = n - floor(sqrt(n))^2   */
  if (temp == 0) {
    factors[nfactors++] = s;
    factors[nfactors++] = s;
    return nfactors;
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
        factors[nfactors++] = n;
        return nfactors;
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
    factors[nfactors++] = n;
    return nfactors;
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
    factors[nfactors++] = n;
    return nfactors;
  }

  if ((q & 1) == 0) q/=2;      /* q was factor or 2*factor   */

  if ( (q == 1) || (q == n) ) {
    factors[nfactors++] = n;
    return nfactors;
  }

  p = n/q;
  /* smallest factor first */
  if (q < p) {  t = p; p = q; q = t; }

  /* printf(" squfof found %lu = %lu * %lu in %ld/%ld rounds\n", n, p, q, jter, iter); */

#if 0
  factors[nfactors++] = p;
  factors[nfactors++] = q;
#else
  if (p >= refactor_above)
    nfactors += squfof_factor(p, factors+nfactors, rounds, refactor_above);
  else
    factors[nfactors++] = p;

  if (q >= refactor_above)
    nfactors += squfof_factor(q, factors+nfactors, rounds, refactor_above);
  else
    factors[nfactors++] = q;
#endif
  return nfactors;
}


/* TODO: Add Jason Papadopoulos's racing SQUFOF */
