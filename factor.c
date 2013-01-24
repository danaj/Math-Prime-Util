#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ptypes.h"
#include "factor.h"
#include "util.h"
#include "sieve.h"
#include "mulmod.h"

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
  while ( (n & 1) == 0 ) { factors[nfactors++] = 2; n /= 2; }

  if (3  <= maxtrial) while ( (n %  3) == 0 ) { factors[nfactors++] =  3; n /=  3; }
  if (5  <= maxtrial) while ( (n %  5) == 0 ) { factors[nfactors++] =  5; n /=  5; }
  if (7  <= maxtrial) while ( (n %  7) == 0 ) { factors[nfactors++] =  7; n /=  7; }
  f = 11;
  m = 11;

  if ( (n < (f*f)) || (maxtrial < f) ) {
    if (n != 1)
      factors[nfactors++] = n;
    return nfactors;
  }

  limit = (UV) (sqrt(n)+0.1);
  if (limit > maxtrial)
    limit = maxtrial;

  /* wheel 30 */
  while (f <= limit) {
    if ( (n%f) == 0 ) {
      UV newlimit;
      do {
        factors[nfactors++] = f;
        n /= f;
      } while ( (n%f) == 0 );
      newlimit = (UV) (sqrt(n)+0.1);
      if (newlimit < limit)  limit = newlimit;
    }
    f += wheeladvance30[m];
    m = nextwheel30[m];
  }
  if (n != 1)
    factors[nfactors++] = n;
  return nfactors;
}


/* Return 0 if n is not a perfect square.  Set sqrtn to int(sqrt(n)) if so.
 *
 * Some simple solutions:
 *
 *     return ( ((m&2)!= 0) || ((m&7)==5) || ((m&11) == 8) )  ?  0  :  1;
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
  UV m;  /* lm */
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
#else
  m = n % 63;
  if ((m*0x3d491df7) & (m*0xc824a9f9) & 0x10f14008) return 0;
#endif
  m = sqrt(n);
  if (n != (m*m))
    return 0;

  if (sqrtn != 0) *sqrtn = m;
  return 1;
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
#if 0
    if (a > (n-2))
      croak("Base %"UVuf" is invalid for input %"UVuf, a, n);
#else
    if (a >= n)
      a %= n;
    if ( (a <= 1) || (a == nm1) )
      continue;
#endif

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

int _XS_is_prob_prime(UV n)
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
  /* Better bases from http://miller-rabin.appspot.com/, 23 Jan 2013 */
  if (n < UVCONST(272161)) {
    bases[0] = UVCONST(62769592775616394);
    nbases = 1;
  } else if (n < UVCONST(466758181)) {
    bases[0] = UVCONST( 91869414    );
    bases[1] = UVCONST( 6346128598129234 );
    nbases = 2;
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
  /* Classic bases */
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
  prob_prime = _XS_miller_rabin(n, bases, nbases);
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

  sqn = (UV) (sqrt(n)+0.1);
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
    s = sqrt( (double)n * (double)i );
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
int pbrent_factor(UV n, UV *factors, UV rounds)
{
  UV f, i, r;
  UV Xi = 2;
  UV Xm = 2;
  UV a = 1;
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
  UV q, f;
  UV a = 2;
  UV savea = 2;
  UV saveq = 2;
  UV j = 1;
  UV sqrtB1 = sqrt(B1);
  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in pminus1_factor");

  for (q = 2; q <= sqrtB1; q = _XS_next_prime(q)) {
    UV k = q*q;
    UV kmin = B1/q;
    while (k <= kmin)
      k *= q;
    a = powmod(a, k, n);
  }
  if (a == 0) { factors[0] = n; return 1; }
  f = gcd_ui(a-1, n);
  if (f == 1) {
    savea = a;
    saveq = q;
    for (; q <= B1; q = _XS_next_prime(q)) {
      a = powmod(a, q, n);
      if ( (j++ % 32) == 0) {
        if (a == 0 || gcd_ui(a-1, n) != 1)
          break;
        savea = a;
        saveq = q;
      }
    }
    if (a == 0) { factors[0] = n; return 1; }
    f = gcd_ui(a-1, n);
  }
  /* If we found more than one factor in stage 1, backup and single step */
  if (f == n) {
    a = savea;
    for (q = saveq; q <= B1; q = _XS_next_prime(q)) {
      UV k = q;
      UV kmin = B1/q;
      while (k <= kmin)
        k *= q;
      a = powmod(a, k, n);
      f = gcd_ui(a-1, n);
      if (f != 1)
        break;
    }
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
    while (q <= B2) {
      UV lastq = q;
      UV qdiff;
      q = _XS_next_prime(q);
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
    }
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
  int i, still_racing;
  UV nn64, mult, f64;
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

    mult_save[i].b0 = sqrt( (double)nn64 );
    mult_save[i].imax = sqrt( (double)mult_save[i].b0 ) / 3;
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
