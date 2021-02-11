#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ptypes.h"
#define FUNC_is_strong_pseudoprime 1
#include "primality.h"
#include "mulmod.h"
#define FUNC_gcd_ui 1
#define FUNC_is_perfect_square
#include "util.h"
#include "montmath.h"  /* Fast Montgomery math */

/* Primality related functions */

/******************************************************************************/


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

static UV select_extra_strong_parameters(UV n, UV increment) {
  int j;
  UV D, P = 3;
  while (1) {
    D = P*P - 4;
    j = jacobi_iu(D, n);
    if (j == 0) { UV g = gcd_ui(D,n);  if (g != 1 && g != n) return 0; }
    if (j == -1) break;
    if (P == (3+20*increment) && is_perfect_square(n)) return 0;
    P += increment;
    if (P > 65535)
      croak("lucas_extrastrong_params: P exceeded 65535");
  }
  if (P >= n)  P %= n;   /* Never happens with increment < 4 */
  return P;
}

/* Fermat pseudoprime */
int is_pseudoprime(UV const n, UV a)
{
  if (n < 4) return (n == 2 || n == 3);
  if (!(n&1) && !(a&1)) return 0;
  if (a < 2) croak("Base %"UVuf" is invalid", a);
  if (a >= n) {
    a %= n;
    if (a <= 1)    return (a == 1);
    if (a == n-1)  return !(a & 1);
  }

#if USE_MONTMATH
  if (n & 1) {   /* The Montgomery code only works for odd n */
    const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);
    const uint64_t monta = (a == 2)  ?  mont_get2(n)  :  mont_geta(a, n);
    return mont_powmod(monta, n-1, n) == mont1;
  }
#endif
  return powmod(a, n-1, n) == 1;  /* a^(n-1) = 1 mod n */
}

/* Euler (aka Euler-Jacobi) pseudoprime:  a^((n-1)/2) = (a|n) mod n */
int is_euler_pseudoprime(UV const n, UV a)
{
  if (n < 5) return (n == 2 || n == 3);
  if (!(n&1)) return 0;
  if (a < 2) croak("Base %"UVuf" is invalid", a);
  if (a > 2) {
    if (a >= n) {
      a %= n;
      if (a <= 1)    return (a == 1);
      if (a == n-1)  return !(a & 1);
    }
    if ((n % a) == 0) return 0;
  }
  {
#if USE_MONTMATH
    const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);
    const uint64_t monta  = mont_geta(a, n);
    UV ap = mont_powmod(monta, (n-1)>>1, n);
    if (ap != mont1 && ap != n-mont1) return 0;
    if (a == 2) {
      uint32_t nmod8 = n & 0x7;
      return (nmod8 == 1 || nmod8 == 7)  ?  (ap == mont1)  :  (ap == n-mont1);
    } else {
      return (kronecker_uu(a,n) >= 0)    ?  (ap == mont1)  :  (ap == n-mont1);
    }
#else
    UV ap = powmod(a, (n-1)>>1, n);
    if (ap != 1 && ap != n-1) return 0;
    if (a == 2) {
      uint32_t nmod8 = n & 0x7;
      return (nmod8 == 1 || nmod8 == 7)  ?  (ap == 1)  :  (ap == n-1);
    } else {
      return (kronecker_uu(a,n) >= 0)    ?  (ap == 1)  :  (ap == n-1);
    }
#endif
  }
}

/* Colin Plumb's extended Euler Criterion test.
 * A tiny bit (~1 percent) faster than base 2 Fermat or M-R.
 * More stringent than base 2 Fermat, but a subset of base 2 M-R.
 */
int is_euler_plumb_pseudoprime(UV const n)
{
  UV ap;
  uint32_t nmod8 = n & 0x7;
  if (n < 5) return (n == 2 || n == 3);
  if (!(n&1)) return 0;
#if USE_MONTMATH
  {
    const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);
    const uint64_t mont2  = mont_get2(n);
    ap = mont_powmod(mont2, (n-1) >> (1 + (nmod8 == 1)), n);
    if (ap ==   mont1)  return (nmod8 == 1 || nmod8 == 7);
    if (ap == n-mont1)  return (nmod8 == 1 || nmod8 == 3 || nmod8 == 5);
  }
#else
  ap = powmod(2, (n-1) >> (1 + (nmod8 == 1)), n);
  if (ap ==   1)  return (nmod8 == 1 || nmod8 == 7);
  if (ap == n-1)  return (nmod8 == 1 || nmod8 == 3 || nmod8 == 5);
#endif
  return 0;
}

/* Miller-Rabin probabilistic primality test
 * Returns 1 if probably prime relative to the bases, 0 if composite.
 * Bases must be between 2 and n-2
 */
int miller_rabin(UV const n, const UV *bases, int nbases)
{
#if USE_MONTMATH
  MPUassert(n > 3, "MR called with n <= 3");
  if ((n & 1) == 0) return 0;
  {
    const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);
    uint64_t a, ma, md, u = n-1;
    int i, j, t = 0;

    while (!(u&1)) {  t++;  u >>= 1;  }
    for (j = 0; j < nbases; j++) {
      a = bases[j];
      if (a < 2)  croak("Base %"UVuf" is invalid", (UV)a);
      if (a >= n) {
        a %= n;
        if (a == 0 || (a == n-1 && a&1)) return 0;
      }
      ma = mont_geta(a,n);
      if (a == 1 || a == n-1 || !ma) continue;
      md = mont_powmod(ma, u, n);
      if (md != mont1 && md != n-mont1) {
        for (i=1; i<t; i++) {
          md = mont_sqrmod(md, n);
          if (md == mont1) return 0;
          if (md == n-mont1) break;
        }
        if (i == t)
          return 0;
      }
    }
  }
#else
  UV d = n-1;
  int b, r, s = 0;

  MPUassert(n > 3, "MR called with n <= 3");
  if ((n & 1) == 0) return 0;

  while (!(d&1)) {  s++;  d >>= 1;  }
  for (b = 0; b < nbases; b++) {
    UV x, a = bases[b];
    if (a < 2)  croak("Base %"UVuf" is invalid", a);
    if (a >= n) {
      a %= n;
      if (a == 0 || (a == n-1 && a&1)) return 0;
    }
    if (a == 1 || a == n-1) continue;
    /* n is a strong pseudoprime to base a if either:
     *    a^d = 1 mod n
     *    a^(d2^r) = -1 mod n for some r: 0 <= r <= s-1
     */
    x = powmod(a, d, n);
    if ( (x == 1) || (x == n-1) )  continue;
    for (r = 1; r < s; r++) {  /* r=0 was just done, test r = 1 to s-1 */
      x = sqrmod(x, n);
      if ( x == n-1 )  break;
      if ( x == 1   )  return 0;
    }
    if (r >= s)  return 0;
  }
#endif
  return 1;
}

int BPSW(UV const n)
{
  if (n < 7) return (n == 2 || n == 3 || n == 5);
  if ((n % 2) == 0 || n == UV_MAX) return 0;

#if !USE_MONTMATH
  return    is_strong_pseudoprime(n, 2)
         && is_almost_extra_strong_lucas_pseudoprime(n,1);
#else
  {
    const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);
    const uint64_t mont2 = mont_get2(n);
    uint64_t md, u = n-1;
    int i, t = 0;
    UV P, V, d, s;

    /* M-R with base 2 */
    while (!(u&1)) {  t++;  u >>= 1;  }
    md = mont_powmod(mont2, u, n);
    if (md != mont1 && md != n-mont1) {
      for (i=1; i<t; i++) {
        md = mont_sqrmod(md, n);
        if (md == mont1) return 0;
        if (md == n-mont1) break;
      }
      if (i == t)
        return 0;
    }
    /* AES Lucas test */
    P = select_extra_strong_parameters(n, 1);
    if (P == 0) return 0;

    d = n+1;
    s = 0;
    while ( (d & 1) == 0 ) {  s++;  d >>= 1; }

    {
      const uint64_t montP = mont_geta(P, n);
      UV W, b;
      W = submod(  mont_mulmod( montP, montP, n),  mont2, n);
      V = montP;
      { UV v = d; b = 1; while (v >>= 1) b++; }
      while (b-- > 1) {
        UV T = submod(  mont_mulmod(V, W, n),  montP, n);
        if ( (d >> (b-1)) & UVCONST(1) ) {
          V = T;
          W = submod(  mont_mulmod(W, W, n),  mont2, n);
        } else {
          W = T;
          V = submod(  mont_mulmod(V, V, n),  mont2, n);
        }
      }
    }

    if (V == mont2 || V == (n-mont2))
      return 1;
    while (s-- > 1) {
      if (V == 0)
        return 1;
      V = submod(  mont_mulmod(V, V, n),  mont2, n);
      if (V == mont2)
        return 0;
    }
  }
  return 0;
#endif
}

/* Alternate modular lucas sequence code.
 * A bit slower than the normal one, but works with even valued n. */
static void alt_lucas_seq(UV* Uret, UV* Vret, UV* Qkret, UV n, UV Pmod, UV Qmod, UV k)
{
  UV Uh, Vl, Vh, Ql, Qh;
  int j, s, m;

  Uh = 1;  Vl = 2;  Vh = Pmod;  Ql = 1;  Qh = 1;
  s = 0; m = 0;
  { UV v = k; while (!(v & 1)) { v >>= 1; s++; } }
  { UV v = k; while (v >>= 1) m++; }

  if (Pmod == 1 && Qmod == (n-1)) {
    int Sl = Ql, Sh = Qh;
    for (j = m; j > s; j--) {
      Sl *= Sh;
      Ql = (Sl==1) ? 1 : n-1;
      if ( (k >> j) & UVCONST(1) ) {
        Sh = -Sl;
        Uh = mulmod(Uh, Vh, n);
        Vl = submod(mulmod(Vh, Vl, n), Ql, n);
        Vh = submod(sqrmod(Vh, n), (Sh==1) ? 2 : n-2, n);
      } else {
        Sh = Sl;
        Uh = submod(mulmod(Uh, Vl, n), Ql, n);
        Vh = submod(mulmod(Vh, Vl, n), Ql, n);
        Vl = submod(sqrmod(Vl, n), (Sl==1) ? 2 : n-2, n);
      }
    }
    Sl *= Sh;
    Ql = (Sl==1) ? 1 : n-1;
    Uh = submod(mulmod(Uh, Vl, n), Ql, n);
    Vl = submod(mulmod(Vh, Vl, n), Ql, n);
    for (j = 0; j < s; j++) {
      Uh = mulmod(Uh, Vl, n);
      Vl = submod(sqrmod(Vl, n), (j>0) ? 2 : n-2, n);
    }
    *Uret = Uh;
    *Vret = Vl;
    *Qkret = (s>0)?1:n-1;
    return;
  }

  /* TODO: Measure speed:  X-Qh-Qh vs X-2*Qh */
  for (j = m; j > s; j--) {
    Ql = mulmod(Ql, Qh, n);
    if ( (k >> j) & UVCONST(1) ) {
      Qh = mulmod(Ql, Qmod, n);
      Uh = mulmod(Uh, Vh, n);
      Vl = submod(mulmod(Vh, Vl, n), mulmod(Pmod, Ql, n), n);
      Vh = submod(sqrmod(Vh, n), mulmod(2, Qh, n), n);
    } else {
      Qh = Ql;
      Uh = submod(mulmod(Uh, Vl, n), Ql, n);
      Vh = submod(mulmod(Vh, Vl, n), mulmod(Pmod, Ql, n), n);
      Vl = submod(sqrmod(Vl, n), mulmod(2, Ql, n), n);
    }
  }
  Ql = mulmod(Ql, Qh, n);
  Qh = mulmod(Ql, Qmod, n);
  Uh = submod(mulmod(Uh, Vl, n), Ql, n);
  Vl = submod(mulmod(Vh, Vl, n), mulmod(Pmod, Ql, n), n);
  Ql = mulmod(Ql, Qh, n);
  for (j = 0; j < s; j++) {
    Uh = mulmod(Uh, Vl, n);
    Vl = submod(sqrmod(Vl, n), mulmod(2, Ql, n), n);
    Ql = sqrmod(Ql, n);
  }
  *Uret = Uh;
  *Vret = Vl;
  *Qkret = Ql;
}

/* Generic Lucas sequence for any appropriate P and Q */
void lucas_seq(UV* Uret, UV* Vret, UV* Qkret, UV n, IV P, IV Q, UV k)
{
  UV U, V, b, Dmod, Qmod, Pmod, Qk;

  MPUassert(n > 0, "lucas_sequence:  modulus n must be > 0");
  if (n == 1) { *Uret = *Vret = *Qkret = 0; return; }

  Qmod = ivmod(Q, n);
  Pmod = ivmod(P, n);
  Dmod = submod( mulmod(Pmod, Pmod, n), mulmod(4, Qmod, n), n);

  if (k == 0) {
    *Uret = 0;
    *Vret = 2 % n;
    *Qkret = 1;
    return;
  }

  if (Dmod == 0 && (b = divmod(Pmod,2,n)) != 0) {
    *Uret = mulmod(k, powmod(b, k-1, n), n);
    *Vret = mulmod(2, powmod(b, k, n), n);
    *Qkret = powmod(Qmod, k, n);
    return;
  }
  if ((n % 2) == 0) {
    alt_lucas_seq(Uret, Vret, Qkret, n, Pmod, Qmod, k);
    return;
  }
  U = 1;
  V = Pmod;
  Qk = Qmod;
  { UV v = k; b = 0; while (v >>= 1) b++; }

  if (Q == 1) {
    while (b--) {
      U = mulmod(U, V, n);
      V = mulsubmod(V, V, 2, n);
      if ( (k >> b) & UVCONST(1) ) {
        UV t2 = mulmod(U, Dmod, n);
        U = muladdmod(U, Pmod, V, n);
        if (U & 1) { U = (n>>1) + (U>>1) + 1; } else { U >>= 1; }
        V = muladdmod(V, Pmod, t2, n);
        if (V & 1) { V = (n>>1) + (V>>1) + 1; } else { V >>= 1; }
      }
    }
  } else if (P == 1 && Q == -1) {
    /* This is about 30% faster than the generic code below.  Since 50% of
     * Lucas and strong Lucas tests come here, I think it's worth doing. */
    int sign = Q;
    while (b--) {
      U = mulmod(U, V, n);
      if (sign == 1) V = mulsubmod(V, V, 2, n);
      else           V = muladdmod(V, V, 2, n);
      sign = 1;   /* Qk *= Qk */
      if ( (k >> b) & UVCONST(1) ) {
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
    while (b--) {
      U = mulmod(U, V, n);
      V = mulsubmod(V, V, addmod(Qk,Qk,n), n);
      Qk = sqrmod(Qk, n);
      if ( (k >> b) & UVCONST(1) ) {
        UV t2 = mulmod(U, Dmod, n);
        U = muladdmod(U, Pmod, V, n);
        if (U & 1) { U = (n>>1) + (U>>1) + 1; } else { U >>= 1; }
        V = muladdmod(V, Pmod, t2, n);
        if (V & 1) { V = (n>>1) + (V>>1) + 1; } else { V >>= 1; }
        Qk = mulmod(Qk, Qmod, n);
      }
    }
  }
  *Uret = U;
  *Vret = V;
  *Qkret = Qk;
}

#define OVERHALF(v)  ( (UV)((v>=0)?v:-v) > (UVCONST(1) << (BITS_PER_WORD/2-1)) )
int lucasuv(IV* U, IV *V, IV P, IV Q, UV k)
{
  IV Uh, Vl, Vh, Ql, Qh;
  int j, s, n;

  if (k == 0) {
    if (U) *U = 0;
    if (V) *V = 2;
    return 1;
  }

  Uh = 1;  Vl = 2;  Vh = P;  Ql = 1;  Qh = 1;
  s = 0; n = 0;
  { UV v = k; while (!(v & 1)) { v >>= 1; s++; } }
  { UV v = k; while (v >>= 1) n++; }

  for (j = n; j > s; j--) {
    if (OVERHALF(Uh) || OVERHALF(Vh) || OVERHALF(Vl) || OVERHALF(Ql) || OVERHALF(Qh)) return 0;
    Ql *= Qh;
    if ( (k >> j) & UVCONST(1) ) {
      Qh = Ql * Q;
      Uh = Uh * Vh;
      Vl = Vh * Vl - P * Ql;
      Vh = Vh * Vh - 2 * Qh;
    } else {
      Qh = Ql;
      Uh = Uh * Vl - Ql;
      Vh = Vh * Vl - P * Ql;
      Vl = Vl * Vl - 2 * Ql;
    }
  }
  if (OVERHALF(Ql) || OVERHALF(Qh)) return 0;
  Ql = Ql * Qh;
  Qh = Ql * Q;
  if (OVERHALF(Uh) || OVERHALF(Vh) || OVERHALF(Vl) || OVERHALF(Ql) || OVERHALF(Qh)) return 0;
  Uh = Uh * Vl - Ql;
  Vl = Vh * Vl - P * Ql;
  Ql = Ql * Qh;
  for (j = 0; j < s; j++) {
    if (OVERHALF(Uh) || OVERHALF(Vl) || OVERHALF(Ql)) return 0;
    Uh *= Vl;
    Vl = Vl * Vl - 2 * Ql;
    Ql *= Ql;
  }
  if (U) *U = Uh;
  if (V) *V = Vl;
  return 1;
}


/******************************************************************************/
/******************************************************************************/

/* Lucas tests:
 *  0: Standard
 *  1: Strong
 *  2: Stronger (Strong + page 1401 extra tests)
 *  3: Extra Strong (Mo/Jones/Grantham)
 *
 * None of them have any false positives for the BPSW test.  Also see the
 * "almost extra strong" test.
 */
int is_lucas_pseudoprime(UV n, int strength)
{
  IV P, Q, D;
  UV U, V, Qk, d, s;

  if (n < 5) return (n == 2 || n == 3);
  if ((n % 2) == 0 || n == UV_MAX) return 0;

  if (strength < 3) {
    UV Du = 5;
    IV sign = 1;
    int j;
    while (1) {
      D = Du * sign;
      j = jacobi_iu(D, n);
      if (j != 1 && Du != n) break;
      if (Du == 21 && is_perfect_square(n)) return 0;
      Du += 2;
      sign = -sign;
    }
    if (j != -1) return 0;
    P = 1;
    Q = (1 - D) / 4;
    if (strength == 2 && Q == -1) P=Q=D=5;  /* Method A* */
    /* Check gcd(n,2QD). gcd(n,2D) already done. */
    Qk = (Q >= 0)  ?  Q % n  :  n-(((UV)(-Q)) % n);
    if (gcd_ui(Qk,n) != 1) return 0;
  } else {
    P = select_extra_strong_parameters(n, 1);
    if (P == 0) return 0;
    Q = 1;
    D = P*P - 4;
  }
  MPUassert( D == (P*P - 4*Q) , "is_lucas_pseudoprime: incorrect DPQ");

#if 0   /* Condition 2, V_n+1 = 2Q mod n */
{ UV us, vs, qs; lucas_seq(&us, &vs, &qs, n, P, Q, n+1); return (vs == addmod(Q,Q,n)); }
#endif
#if 0   /* Condition 3, n is a epsp(Q) */
return is_euler_pseudoprime(n,Qk);
#endif

  d = n+1;
  s = 0;
  if (strength > 0)
    while ( (d & 1) == 0 ) {  s++;  d >>= 1; }

#if USE_MONTMATH
  {
    const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);
    const uint64_t mont2 = mont_get2(n);
    const uint64_t montP = (P == 1) ? mont1
                         : (P >= 0) ? mont_geta(P, n)
                         : n - mont_geta(-P, n);
    const uint64_t montQ = (Q == 1) ? mont1
                         : (Q >= 0) ? mont_geta(Q, n)
                         : n - mont_geta(-Q, n);
    const uint64_t montD = (D >= 0) ? mont_geta(D, n)
                         : n - mont_geta(-D, n);
    UV b;
    { UV v = d; b = 0; while (v >>= 1) b++; }

    /* U, V, Qk, and mont* are in Montgomery space */
    U = mont1;
    V = montP;

    if (Q == 1 || Q == -1) {   /* Faster code for |Q|=1, also opt for P=1 */
      int sign = Q;
      while (b--) {
        U = mont_mulmod(U, V, n);
        if (sign == 1) V = submod( mont_sqrmod(V,n), mont2, n);
        else           V = addmod( mont_sqrmod(V,n), mont2, n);
        sign = 1;
        if ( (d >> b) & UVCONST(1) ) {
          UV t2 = mont_mulmod(U, montD, n);
          if (P == 1) {
            U = addmod(U, V, n);
            V = addmod(V, t2, n);
          } else {
            U = addmod( mont_mulmod(U, montP, n), V, n);
            V = addmod( mont_mulmod(V, montP, n), t2, n);
          }
          if (U & 1) { U = (n>>1) + (U>>1) + 1; } else { U >>= 1; }
          if (V & 1) { V = (n>>1) + (V>>1) + 1; } else { V >>= 1; }
          sign = Q;
        }
      }
      Qk = (sign == 1) ? mont1 : n-mont1;
    } else {
      Qk = montQ;
      while (b--) {
        U = mont_mulmod(U, V, n);
        V = submod( mont_sqrmod(V,n), addmod(Qk,Qk,n), n);
        Qk = mont_sqrmod(Qk,n);
        if ( (d >> b) & UVCONST(1) ) {
          UV t2 = mont_mulmod(U, montD, n);
          U = addmod( mont_mulmod(U, montP, n), V, n);
          if (U & 1) { U = (n>>1) + (U>>1) + 1; } else { U >>= 1; }
          V = addmod( mont_mulmod(V, montP, n), t2, n);
          if (V & 1) { V = (n>>1) + (V>>1) + 1; } else { V >>= 1; }
          Qk = mont_mulmod(Qk, montQ, n);
        }
      }
    }

    if (strength == 0) {
      if (U == 0)
        return 1;
    } else if (strength == 1) {
      if (U == 0)
        return 1;
      while (s--) {
        if (V == 0)
          return 1;
        if (s) {
          V = submod( mont_sqrmod(V,n), addmod(Qk,Qk,n), n);
          Qk = mont_sqrmod(Qk,n);
        }
      }
    } else if (strength == 2) {
      UV Ql = 0, Qj = 0;
      int qjacobi, is_slpsp = 0;
      if (U == 0)
        is_slpsp = 1;
      while (s--) {
        if (V == 0)
          is_slpsp = 1;
        Ql = Qk;
        V = submod( mont_sqrmod(V,n), addmod(Qk,Qk,n), n);
        Qk = mont_sqrmod(Qk,n);
      }
      if (!is_slpsp)                  return 0; /* slpsp */
      if (V != addmod(montQ,montQ,n)) return 0; /* V_{n+1} != 2Q mod n */
      qjacobi = jacobi_iu(Q,n);
      Qj = (qjacobi == 0) ? 0 : (qjacobi == 1) ? montQ : n-montQ;
      if (Ql != Qj)                   return 0; /* n is epsp base Q */
      return 1;
    } else {
      if ( U == 0 && (V == mont2 || V == (n-mont2)) )
        return 1;
      s--;
      while (s--) {
        if (V == 0)
          return 1;
        if (s)
          V = submod( mont_sqrmod(V,n), mont2, n);
      }
    }
    return 0;
  }
#else
  lucas_seq(&U, &V, &Qk, n, P, Q, d);

  if (strength == 0) {
    if (U == 0)
      return 1;
  } else if (strength == 1) {
    if (U == 0)
      return 1;
    /* Now check to see if V_{d*2^r} == 0 for any 0 <= r < s */
    while (s--) {
      if (V == 0)
        return 1;
      if (s) {
        V = mulsubmod(V, V, addmod(Qk,Qk,n), n);
        Qk = sqrmod(Qk, n);
      }
    }
  } else if (strength == 2) {
    UV Ql = 0, Qj = 0;
    UV Qu = (Q >= 0)  ?  Q % n  :  n-(((UV)(-Q)) % n);
    int qjacobi, is_slpsp = 0;
    if (U == 0)
      is_slpsp = 1;
    while (s--) {
      if (V == 0)
        is_slpsp = 1;
      Ql = Qk;
      V = mulsubmod(V, V, addmod(Qk,Qk,n), n);
      Qk = sqrmod(Qk, n);
    }
    if (!is_slpsp)                  return 0; /* slpsp */
    if (V != addmod(Qu,Qu,n))       return 0; /* V_{n+1} != 2Q mod n */
    qjacobi = jacobi_iu(Q,n);
    Qj = (qjacobi == 0) ? 0 : (qjacobi == 1) ? Qu : n-Qu;
    if (Ql != Qj)                   return 0; /* n is epsp base Q */
    return 1;
  } else {
    if ( U == 0 && (V == 2 || V == (n-2)) )
      return 1;
    /* Now check to see if V_{d*2^r} == 0 for any 0 <= r < s-1 */
    s--;
    while (s--) {
      if (V == 0)
        return 1;
      if (s)
        V = mulsubmod(V, V, 2, n);
    }
  }
  return 0;
#endif
}

/* A generalization of Pari's shortcut to the extra-strong Lucas test.
 *
 * This only calculates and tests V, which means less work, but it does result
 * in a few more pseudoprimes than the full extra-strong test.
 *
 * I've added a gcd check at the top, which needs to be done and also results
 * in fewer pseudoprimes.  Pari always does trial division to 100 first so
 * is unlikely to come up there.
 *
 * increment:  1 for Baillie OEIS, 2 for Pari.
 *
 * With increment = 1, these results will be a subset of the extra-strong
 * Lucas pseudoprimes.  With increment = 2, we produce Pari's results.
 */
int is_almost_extra_strong_lucas_pseudoprime(UV n, UV increment)
{
  UV P, V, W, d, s, b;

  if (n < 13) return (n == 2 || n == 3 || n == 5 || n == 7 || n == 11);
  if ((n % 2) == 0 || n == UV_MAX) return 0;
  if (increment < 1 || increment > 256)
    croak("Invalid lucas parameter increment: %"UVuf"\n", increment);

  /* Ensure small primes work with large increments. */
  if ( (increment >= 16 && n <= 331) || (increment > 148 && n <= 631) )
    return !!is_prob_prime(n);

  P = select_extra_strong_parameters(n, increment);
  if (P == 0) return 0;

  d = n+1;
  s = 0;
  while ( (d & 1) == 0 ) {  s++;  d >>= 1; }
  { UV v = d; b = 0; while (v >>= 1) b++; }

#if USE_MONTMATH
  {
    const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);
    const uint64_t mont2 = mont_get2(n);
    const uint64_t montP = mont_geta(P, n);
    W = submod(  mont_mulmod( montP, montP, n),  mont2, n);
    V = montP;
    while (b--) {
      UV T = submod(  mont_mulmod(V, W, n),  montP, n);
      if ( (d >> b) & UVCONST(1) ) {
        V = T;
        W = submod(  mont_mulmod(W, W, n),  mont2, n);
      } else {
        W = T;
        V = submod(  mont_mulmod(V, V, n),  mont2, n);
      }
    }

    if (V == mont2 || V == (n-mont2))
      return 1;
    s--;
    while (s--) {
      if (V == 0)
        return 1;
      if (s)
        V = submod(  mont_mulmod(V, V, n),  mont2, n);
    }
    return 0;
  }
#else
  W = mulsubmod(P, P, 2, n);
  V = P;
  while (b--) {
    UV T = mulsubmod(V, W, P, n);
    if ( (d >> b) & UVCONST(1) ) {
      V = T;
      W = mulsubmod(W, W, 2, n);
    } else {
      W = T;
      V = mulsubmod(V, V, 2, n);
    }
  }
  if (V == 2 || V == (n-2))
    return 1;
  while (s-- > 1) {
    if (V == 0)
      return 1;
    V = mulsubmod(V, V, 2, n);
    if (V == 2)
      return 0;
  }
  return 0;
#endif
}


typedef struct {
  unsigned short div;
  unsigned short period;
  unsigned short offset;
} _perrin;
#define NPERRINDIV 19
/* 1112 mask bytes */
static const uint32_t _perrinmask[] = {22,523,514,65890,8519810,130,4259842,0,526338,2147483904U,1644233728,1,8194,1073774592,1024,134221824,128,512,181250,2048,0,1,134217736,1049600,524545,2147500288U,0,524290,536870912,32768,33554432,2048,0,2,2,256,65536,64,536875010,32768,256,64,0,32,1073741824,0,1048576,1048832,371200000,0,0,536887552,32,2147487744U,2097152,32768,1024,0,1024,536870912,128,512,0,0,512,0,2147483650U,45312,128,0,8388640,0,8388608,8388608,0,2048,4096,92800000,262144,0,65536,4,0,4,4,4194304,8388608,1075838976,536870956,0,134217728,8192,0,8192,8192,0,2,0,268435458,134223392,1073741824,268435968,2097152,67108864,0,8192,1073741840,0,0,128,0,0,512,1450000,8,131136,536870928,0,4,2097152,4096,64,0,32768,0,0,131072,371200000,2048,33570816,4096,32,1024,536870912,1048576,16384,0,8388608,0,0,0,2,512,0,128,0,134217728,2,32,0,0,0,0,8192,0,1073742080,536870912,0,4096,16777216,526336,32,0,65536,33554448,708,67108864,2048,0,0,536870912,0,536870912,33554432,33554432,2147483648U,512,64,0,1074003968,512,0,524288,0,0,0,67108864,524288,1048576,0,131076,0,33554432,131072,0,2,8390656,16384,16777216,134217744,0,131104,0,2,32768,0,0,0,1450000,32768,0,0,0,0,0,16,0,1024,16400,1048576,32,1024,0,260,536870912,269484032,0,16384,0,524290,0,0,512,65536,0,0,0,134217732,0,67108880,536887296,0,0,32,0,65568,0,524288,2147483648U,0,4096,4096,134217984,268500992,0,33554432,131072,0,0,0,16777216,0,0,0,0,0,524288,0,0,67108864,0,0,2,0,2,32,1024,0};
static _perrin _perrindata[NPERRINDIV] = {
  {2, 7, 0},
  {3, 13, 1},
  {4, 14, 2},
  {5, 24, 3},
  {7, 48, 4},
  {9, 39, 6},
  {11, 120, 8},
  {13, 183, 12},
  {17, 288, 18},
  {19, 180, 27},
  {23, 22, 33},
  {25, 120, 34},
  {29, 871, 38},
  {31, 993, 66},
  {37, 1368, 98},
  {41, 1723, 141},
  {43, 231, 195},
  {47, 2257, 203},
  {223, 111, 274}
};

/* Calculate signature using the doubling rule from Adams and Shanks 1982 */
static void calc_perrin_sig(UV* S, UV n) {
#if USE_MONTMATH
  uint64_t npi = 0, mont1;
  int i;
#endif
  UV T[6], T01, T34, T45;
  int b;

  /* Signature for n = 1 */
  S[0] = 1; S[1] = n-1; S[2] = 3;   S[3] = 3; S[4] = 0; S[5] = 2;
  if (n <= 1) return;

#if USE_MONTMATH
  if ( (n&1) ) {
    npi = mont_inverse(n);
    mont1 = mont_get1(n);
    S[0] = mont1;  S[1] = n-mont1;  S[5] = addmod(mont1,mont1,n);
    S[2] = addmod(S[5],mont1,n);  S[3] = S[2];
  }
#endif

  /* Bits in n */
  { UV v = n; b = 1; while (v >>= 1) b++; }

  while (b-- > 1) {
    /* Double */
#if USE_MONTMATH
    if (n&1) {
      T[0] = submod(submod(mont_sqrmod(S[0],n), S[5],n), S[5],n);
      T[1] = submod(submod(mont_sqrmod(S[1],n), S[4],n), S[4],n);
      T[2] = submod(submod(mont_sqrmod(S[2],n), S[3],n), S[3],n);
      T[3] = submod(submod(mont_sqrmod(S[3],n), S[2],n), S[2],n);
      T[4] = submod(submod(mont_sqrmod(S[4],n), S[1],n), S[1],n);
      T[5] = submod(submod(mont_sqrmod(S[5],n), S[0],n), S[0],n);
    } else
#endif
    {
      T[0] = submod(submod(sqrmod(S[0],n), S[5],n), S[5],n);
      T[1] = submod(submod(sqrmod(S[1],n), S[4],n), S[4],n);
      T[2] = submod(submod(sqrmod(S[2],n), S[3],n), S[3],n);
      T[3] = submod(submod(sqrmod(S[3],n), S[2],n), S[2],n);
      T[4] = submod(submod(sqrmod(S[4],n), S[1],n), S[1],n);
      T[5] = submod(submod(sqrmod(S[5],n), S[0],n), S[0],n);
    }
    /* Move to S, filling in */
    T01 = submod(T[2], T[1], n);
    T34 = submod(T[5], T[4], n);
    T45 = addmod(T34, T[3], n);
    if ( (n >> (b-1)) & 1U ) {
      S[0] = T[0];   S[1] = T01;    S[2] = T[1];
      S[3] = T[4];   S[4] = T45;    S[5] = T[5];
   } else {
      S[0] = T01;    S[1] = T[1];   S[2] = addmod(T01,T[0],n);
      S[3] = T34;    S[4] = T[4];   S[5] = T45;
    }
  }
#if USE_MONTMATH
  if (n&1) { /* Recover result from Montgomery form */
    for (i = 0; i < 6; i++)
      S[i] = mont_recover(S[i],n);
  }
#endif
}

int is_perrin_pseudoprime(UV n, int restricted)
{
  int jacobi, i;
  UV S[6];

  if (n < 3) return (n >= 2);
  if (!(n&1) && restricted > 2) return 0;  /* Odds only for restrict > 2 */

  /* Hard code the initial tests.  60% of composites caught by 4 tests. */
  {
    uint32_t n32 = n % 10920;
    if (!(n32&1) && !((   22 >> (n32% 7)) & 1)) return 0;
    if (!(n32%3) && !((  523 >> (n32%13)) & 1)) return 0;
    if (!(n32%5) && !((65890 >> (n32%24)) & 1)) return 0;
    if (!(n32%4) && !((  514 >> (n32%14)) & 1)) return 0;
  }
  for (i = 4; i < NPERRINDIV; i++) {
    if ((n % _perrindata[i].div) == 0) {
      const uint32_t *mask = _perrinmask + _perrindata[i].offset;
      unsigned short mod = n % _perrindata[i].period;
      if (!((mask[mod/32] >> (mod%32)) & 1))
        return 0;
    }
  }
  /* Depending on which filters are used, 10-20% of composites are left. */

  calc_perrin_sig(S, n);

  if (S[4] != 0)        return 0;        /* P(n) = 0 mod n */
  if (restricted == 0)  return 1;

  if (S[1] != n-1)      return 0;        /* P(-n) = -1 mod n */
  if (restricted == 1)  return 1;

  /* Full restricted test looks for an acceptable signature.
   *
   * restrict = 2 is Adams/Shanks without quadratic form test
   *
   * restrict = 3 is Arno or Grantham: No qform, also reject mults of 2 and 23
   *
   * See:
   *     Adams/Shanks 1982 pages 257-261
   *     Arno 1991 pages 371-372
   *     Grantham 2000 pages 5-6
   */

  jacobi = kronecker_su(-23,n);

  if (jacobi == -1) { /* Q-type */

    UV B = S[2], B2 = sqrmod(B,n);
    UV A = submod(addmod(1,mulmod(B,3,n),n),B2,n);
    UV C = submod(mulmod(B2,3,n),2,n);
    if (S[0] == A  &&  S[2] == B  &&  S[3] == B  &&  S[5] == C  &&
        B != 3  &&  submod(mulmod(B2,B,n),B,n) == 1) {
      MPUverbose(2, "%"UVuf" Q-Type  %"UVuf" -1 %"UVuf"  %"UVuf" 0 %"UVuf"\n", n, A, B, B, C);
      return 1;
    }

  } else {            /* S-Type or I-Type */

    if (jacobi == 0 && n != 23 && restricted > 2) {
      MPUverbose(2, "%"UVuf" Jacobi %d\n",n,jacobi);
      return 0;  /* Adams/Shanks allows (-23|n) = 0 for S-Type */
    }

    if (S[0] == 1  &&  S[2] == 3  &&  S[3] == 3  &&  S[5] == 2) {
      MPUverbose(2, "%"UVuf" S-Type  1 -1 3  3 0 2\n",n);
      return 1;
    } else if (S[0] == 0  &&  S[5] == n-1  &&  S[2] != S[3]  &&
               addmod(S[2],S[3],n) == n-3  && sqrmod(submod(S[2],S[3],n),n) == n-(23%n)) {
      MPUverbose(2, "%"UVuf" I-Type  0 -1 %"UVuf"  %"UVuf" 0 -1\n",n, S[2], S[3]);
      return 1;
    }

  }
  MPUverbose(2, "%"UVuf" ? %2d ?  %"UVuf" -1 %"UVuf"  %"UVuf" 0 %"UVuf"\n", n, jacobi, S[0],S[2],S[3],S[5]);
  return 0;
}

int is_frobenius_pseudoprime(UV n, IV P, IV Q)
{
  UV U, V, Qk, Vcomp;
  int k = 0;
  IV D;
  UV Du, Pu, Qu;

  if (n < 7) return (n == 2 || n == 3 || n == 5);
  if ((n % 2) == 0 || n == UV_MAX) return 0;

  if (P == 0 && Q == 0) {
    P = -1; Q = 2;
    if (n == 7) P = 1;  /* So we don't test kronecker(-7,7) */
    do {
      P += 2;
      if (P == 3) P = 5;  /* P=3,Q=2 -> D=9-8=1 => k=1, so skip */
      D = P*P-4*Q;
      Du = D >= 0 ? D : -D;
      k = kronecker_su(D, n);
      if (P == 10001 && is_perfect_square(n)) return 0;
    } while (k == 1);
    if (k == 0) return 0;
    /* D=P^2-8 will not be a perfect square */
    MPUverbose(1, "%"UVuf" Frobenius (%"IVdf",%"IVdf") : x^2 - %"IVdf"x + %"IVdf"\n", n, P, Q, P, Q);
    Vcomp = 4;
  } else {
    D = P*P-4*Q;
    Du = D >= 0 ? D : -D;
    if (D != 5 && is_perfect_square(Du))
      croak("Frobenius invalid P,Q: (%"IVdf",%"IVdf")", P, Q);
  }
  Pu = (P >= 0 ? P : -P) % n;
  Qu = (Q >= 0 ? Q : -Q) % n;

  Qk = gcd_ui(n, Pu*Qu*Du);
  if (Qk != 1) {
    if (Qk == n) return !!is_prob_prime(n);
    return 0;
  }
  if (k == 0) {
    k = kronecker_su(D, n);
    if (k == 0) return 0;
    if (k == 1) {
      Vcomp = 2;
    } else {
      Qu = addmod(Qu,Qu,n);
      Vcomp = (Q >= 0)  ?  Qu  :  n-Qu;
    }
  }

  lucas_seq(&U, &V, &Qk, n, P, Q, n-k);
  /* MPUverbose(1, "%"UVuf" Frobenius U = %"UVuf" V = %"UVuf"\n", n, U, V); */
  if (U == 0 && V == Vcomp) return 1;
  return 0;
}

/*
 * Khashin, July 2018, https://arxiv.org/pdf/1807.07249.pdf
 * "Evaluation of the Effectiveness of the Frobenius Primality Test"
 *
 * See also the earlier https://arxiv.org/abs/1307.7920
 * "Counterexamples for Frobenius primality test"
 *
 * 1. select c as first in [-1,2,3,4,5,6,...] where (c|n)=-1
 * 2. Check this holds:
 *    (2+sqrt(c)^n = 2-sqrt(c)  mod n   for c = -1,2
 *    (1+sqrt(c)^n = 1-sqrt(c)  mod n   for c = 3,4,5,6,...
 *
 * The paper claims there are no 64-bit counterexamples.
 */
int is_frobenius_khashin_pseudoprime(UV n)
{
  int k = 2;
  UV ea, ra, rb, a, b, d = n-1, c = 1;

  if (n < 7) return (n == 2 || n == 3 || n == 5);
  if ((n % 2) == 0 || n == UV_MAX) return 0;
  if (is_perfect_square(n)) return 0;

  if      (n % 4 == 3) c = d;
  else if (n % 8 == 5) c = 2;
  else
  do { /* c = first odd prime where (c|n)=-1 */
    c += 2;
    if (c==9 || (c>=15 && (!(c%3) || !(c%5) || !(c%7) || !(c%11) || !(c%13))))
      continue;
    k = kronecker_uu(c, n);
  } while (k == 1);
  if (k == 0 || (k == 2 && n % 3 == 0)) return 0;

#if USE_MONTMATH
  {
    const uint64_t npi = mont_inverse(n);
    const uint64_t mont1 = mont_get1(n);
    const uint64_t montc = mont_geta(c, n);
    ra = a = ea = (k == 2) ? mont_get2(n) : mont1;
    rb = b = mont1;
    while (d) {
      if (d & 1) {
        UV ta=ra, tb=rb;
        ra = addmod( mont_mulmod(ta,a,n), mont_mulmod(mont_mulmod(tb,b,n),montc,n), n );
        rb = addmod( mont_mulmod(tb,a,n), mont_mulmod(ta,b,n), n);
      }
      d >>= 1;
      if (d) {
        UV t = mont_mulmod(mont_mulmod(b,b,n),montc,n);
        b = mont_mulmod(b,a,n);
        b = addmod(b,b,n);
        a = addmod(mont_mulmod(a,a,n),t,n);
      }
    }
    return (ra == ea && rb == n-mont1);
  }
#else
  ra = a = ea = (k == 2) ? 2 : 1;
  rb = b = 1;
  while (d) {
    if (d & 1) {
      /* This is faster than the 3-mulmod 5-addmod version */
      UV ta=ra, tb=rb;
      ra = addmod( mulmod(ta,a,n), mulmod(mulmod(tb,b,n),c,n), n );
      rb = addmod( mulmod(tb,a,n), mulmod(ta,b,n), n);
    }
    d >>= 1;
    if (d) {
      UV t = mulmod(sqrmod(b,n),c,n);
      b = mulmod(b,a,n);
      b = addmod(b,b,n);
      a = addmod(sqrmod(a,n),t,n);
    }
  }
  return (ra == ea && rb == n-1);
#endif
}

/*
 * The Frobenius-Underwood test has no known counterexamples below 2^50, but
 * has not been extensively tested above that.  This is the Minimal Lambda+2
 * test from section 9 of "Quadratic Composite Tests" by Paul Underwood.
 *
 * It is generally slower than the AES Lucas test, but for large values is
 * competitive with the BPSW test.  Since our BPSW is known to have no
 * counterexamples under 2^64, while the results of this test are unknown,
 * it is mainly useful for numbers larger than 2^64 as an additional
 * non-correlated test.
 */
int is_frobenius_underwood_pseudoprime(UV n)
{
  int j, bit;
  UV x, result, a, b, np1, len, t1;
  IV t;

  if (n < 7) return (n == 2 || n == 3 || n == 5);
  if ((n % 2) == 0 || n == UV_MAX) return 0;

  for (x = 0; x < 1000000; x++) {
    if (x==2 || x==4 || x==7 || x==8 || x==10 || x==14 || x==16 || x==18)
      continue;
    t = (IV)(x*x) - 4;
    j = jacobi_iu(t, n);
    if (j == -1) break;
    if (j == 0 || (x == 20 && is_perfect_square(n)))
      return 0;
  }
  if (x >= 1000000) croak("FU test failure, unable to find suitable a");
  t1 = gcd_ui(n, (x+4)*(2*x+5));
  if (t1 != 1 && t1 != n)
    return 0;
  np1 = n+1;
  { UV v = np1; len = 1;  while (v >>= 1) len++; }

#if USE_MONTMATH
  {
    const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);
    const uint64_t mont2 = mont_get2(n);
    const uint64_t mont5 = mont_geta(5, n);

    x = mont_geta(x, n);
    a = mont1;
    b = mont2;

    if (x == 0) {
      result = mont5;
      for (bit = len-2; bit >= 0; bit--) {
        t1 = addmod(b, b, n);
        b = mont_mulmod(submod(b, a, n), addmod(b, a, n), n);
        a = mont_mulmod(a, t1, n);
        if ( (np1 >> bit) & UVCONST(1) ) {
          t1 = b;
          b = submod( addmod(b, b, n), a, n);
          a = addmod( addmod(a, a, n), t1, n);
        }
      }
    } else {
      UV multiplier = addmod(x, mont2, n);
      result = addmod( addmod(x, x, n), mont5, n);
      for (bit = len-2; bit >= 0; bit--) {
        t1 = addmod( mont_mulmod(a, x, n), addmod(b, b, n), n);
        b = mont_mulmod(submod(b, a, n), addmod(b, a, n), n);
        a = mont_mulmod(a, t1, n);
        if ( (np1 >> bit) & UVCONST(1) ) {
          t1 = b;
          b = submod( addmod(b, b, n), a, n);
          a = addmod( mont_mulmod(a, multiplier, n), t1, n);
        }
      }
    }
    return (a == 0 && b == result);
  }
#else
  a = 1;
  b = 2;

  if (x == 0) {
    result = 5;
    for (bit = len-2; bit >= 0; bit--) {
      t1 = addmod(b, b, n);
      b = mulmod( submod(b, a, n), addmod(b, a, n), n);
      a = mulmod(a, t1, n);
      if ( (np1 >> bit) & UVCONST(1) ) {
        t1 = b;
        b = submod( addmod(b, b, n), a, n);
        a = addmod( addmod(a, a, n), t1, n);
      }
    }
  } else {
    UV multiplier = addmod(x, 2, n);
    result = addmod( addmod(x, x, n), 5, n);
    for (bit = len-2; bit >= 0; bit--) {
      t1 = addmod( mulmod(a, x, n), addmod(b, b, n), n);
      b = mulmod(submod(b, a, n), addmod(b, a, n), n);
      a = mulmod(a, t1, n);
      if ( (np1 >> bit) & UVCONST(1) ) {
        t1 = b;
        b = submod( addmod(b, b, n), a, n);
        a = addmod( mulmod(a, multiplier, n), t1, n);
      }
    }
  }

  MPUverbose(2, "%"UVuf" is %s with x = %"UVuf"\n", n, (a == 0 && b == result) ? "probably prime" : "composite", x);
  if (a == 0 && b == result)
    return 1;
  return 0;
#endif
}

/* We have a native-UV Lucas-Lehmer test with simple pretest.  If 2^p-1 is
 * prime but larger than a UV, we'll have to bail, and they'll run the nice
 * GMP version.  However, they're just asking if this is a Mersenne prime, and
 * there are millions of CPU years that have gone into enumerating them, so
 * instead we'll use a table. */
#define NUM_KNOWN_MERSENNE_PRIMES 51
static const uint32_t _mersenne_primes[NUM_KNOWN_MERSENNE_PRIMES] = {2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933};
#define LAST_CHECKED_MERSENNE 50845813
int is_mersenne_prime(UV p)
{
  int i;
  for (i = 0; i < NUM_KNOWN_MERSENNE_PRIMES; i++)
    if (p == _mersenne_primes[i])
      return 1;
  return (p < LAST_CHECKED_MERSENNE) ? 0 : -1;
}
int lucas_lehmer(UV p)
{
  UV k, V, mp;

  if (p == 2) return 1;
  if (!is_prob_prime(p))  return 0;
  if (p > BITS_PER_WORD) croak("lucas_lehmer with p > BITS_PER_WORD");
  V = 4;
  mp = UV_MAX >> (BITS_PER_WORD - p);
  for (k = 3; k <= p; k++) {
    V = mulsubmod(V, V, 2, mp);
  }
  return (V == 0);
}

/******************************************************************************/

/* Hashing similar to Forišek and Jančina 2015, trial with 2/3/5/7. */
static const uint16_t mr_bases_hash32[256] = {
29783,551,1651,1050,26680,311,635,9038,4165,4483,678,2337,5737,47570,1364,556,6066,16981,18532,3132,4309,2727,11084,895,1178,378,1598,7287,3162,3209,5733,3111,29966,572,3502,4646,3313,1473,274,6008,1706,11135,19247,430,2057,6193,609,2101,817,4717,139,1088,979,1347,321,3413,7525,2085,12629,477,9139,905,1411,683,19112,1136,7470,791,10777,1459,171,1988,3415,4269,326,2452,12800,14790,1077,1952,1879,1826,8240,211,8181,2665,1719,233,1097,3704,1768,37882,2703,1750,7906,5745,1313,509,723,3367,2928,6497,12208,312,2475,25734,8840,842,7051,230,6031,342,293,827,486,2683,331,2620,474,6707,3486,3780,1090,1132,9556,13156,14701,83,4689,71,18259,4422,5107,740,1778,542,2689,27289,488,7988,4343,6227,1266,1750,57,1836,2443,2759,4377,2404,2605,3084,194,10699,9557,13655,7045,2365,2495,13527,11969,1420,9144,296,37747,8627,2193,4511,270,5373,2248,1422,4838,131,4362,12629,2891,7770,20841,67,238,1464,1748,4152,6674,94,5067,13594,951,2378,1026,13625,380,4037,7245,1217,2244,11503,2076,667,3997,1889,2871,8201,11184,912,1460,2204,3730,103,7516,19769,971,2405,1077,3228,2082,846,4803,872,22946,101,1490,1223,8985,334,83,5587,280,18791,2067,10136,279,15162,7602,13298,4227,346,30,3472,2374,651,1948,304,5397,235,4086,2719,10042,145,572,1328,2204,4125,1762,5304
};

/* 1) This will be always correct for any 32-bit input.
 * 2) For performance, assume trial division has already been done.
 * 3) Make it as fast as possible.
 */
int MR32(uint32_t n) {
  uint32_t x = n;

  if (x < 11) return (x == 2 || x == 3 || x == 5 || x == 7);
  if (!(x&1) || !(x%3) || !(x%5) || !(x%7)) return 0;

  /* The {xorr,mul}+,xorr hash structure is pretty common.  Three rounds
   * gives wonderful dispersion.  Two (like the Mueller hash) gives plenty.
   * One is not enough in general, but it works for us here. */
  x = ((x >> 15) ^ x) * 0xda64aaab;
  x = ((x >> 16) ^ x) & 255;
  return is_strong_pseudoprime(n, mr_bases_hash32[x]);
}

/******************************************************************************/

int is_prob_prime(UV n)
{
#if BITS_PER_WORD == 64
  if (n > UVCONST(4294967295)) { /* input is >= 2^32, UV is 64-bit*/
    if (!(n%2) || !(n%3) || !(n%5) || !(n%7))       return 0;
    if (!(n%11) || !(n%13) || !(n%17) || !(n%19) ||
        !(n%23) || !(n%29) || !(n%31) || !(n%37) ||
        !(n%41) || !(n%43) || !(n%47) || !(n%53))   return 0;
    if (!(n%59) || !(n%61) || !(n%67) || !(n%71))   return 0;
    if (!(n%73) || !(n%79) || !(n%83) || !(n%89))   return 0;
    /* AESLSP test costs about 1.5 Selfridges, vs. ~2.2 for strong Lucas.
     * This makes the full BPSW test cost about 2.5x M-R tests for a prime. */
    return 2*BPSW(n);
  } else {
#else
  {
#endif
    uint32_t x = n;
    if (x < 11) return (x == 2 || x == 3 || x == 5 || x == 7) ? 2 : 0;
    if (!(x%2) || !(x%3) || !(x%5) || !(x%7))       return 0;
    if (x <  121) /* 11*11 */                       return 2;
    if (!(x%11) || !(x%13) || !(x%17) || !(x%19) ||
        !(x%23) || !(x%29) || !(x%31) || !(x%37) ||
        !(x%41) || !(x%43) || !(x%47) || !(x%53))   return 0;
    if (x < 3481) /* 59*59 */                       return 2;
    /* Trial division crossover point depends on platform */
    if (!USE_MONTMATH && n < 200000) {
      uint32_t f = 59;
      uint32_t limit = isqrt(n);
      while (f <= limit) {
        { if ((x%f) == 0)  return 0; }  f += 2;
        { if ((x%f) == 0)  return 0; }  f += 6;
        { if ((x%f) == 0)  return 0; }  f += 4;
        { if ((x%f) == 0)  return 0; }  f += 2;
        { if ((x%f) == 0)  return 0; }  f += 4;
        { if ((x%f) == 0)  return 0; }  f += 2;
        { if ((x%f) == 0)  return 0; }  f += 4;
        { if ((x%f) == 0)  return 0; }  f += 6;
      }
      return 2;
    }
    return 2*MR32(x);
  }
}
