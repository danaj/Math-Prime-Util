#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ptypes.h"
#include "primality.h"
#include "mulmod.h"
#define FUNC_gcd_ui 1
#define FUNC_is_perfect_square
#include "util.h"
#define FUNC_monty_mr64
#include "montmath.h"  /* Fast Montgomery math */

/* Primality related functions */

#if !USE_MONTMATH
static const UV mr_bases_const2[1] = {2};
#endif
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
  UV P = 3;
  while (1) {
    UV D = P*P - 4;
    /* TODO: We should not do this test. */
    if (gcd_ui(D, n) > 1 && gcd_ui(D, n) != n) return 0;
    if (jacobi_iu(D, n) == -1)
      break;
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
  return monty_mr64((uint64_t)n, bases, nbases);
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
  return 1;
#endif
}

int BPSW(UV const n)
{
  if (n < 7) return (n == 2 || n == 3 || n == 5);
  if ((n % 2) == 0 || n == UV_MAX) return 0;

#if !USE_MONTMATH
  return    miller_rabin(n, mr_bases_const2, 1)
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

  MPUassert(n > 1, "lucas_sequence:  modulus n must be > 1");
  if (k == 0) {
    *Uret = 0;
    *Vret = 2;
    *Qkret = Q;
    return;
  }

  Qmod = (Q < 0)  ?  (UV) (Q + (IV)(((-Q/n)+1)*n))  :  (UV)Q % n;
  Pmod = (P < 0)  ?  (UV) (P + (IV)(((-P/n)+1)*n))  :  (UV)P % n;
  Dmod = submod( mulmod(Pmod, Pmod, n), mulmod(4, Qmod, n), n);
  if (Dmod == 0) {
    b = Pmod >> 1;
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
int lucasu(IV* U, IV P, IV Q, UV k)
{
  IV Uh, Vl, Vh, Ql, Qh;
  int j, s, n;

  if (U == 0) return 0;
  if (k == 0) { *U = 0; return 1; }

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
  *U = Uh;
  return 1;
}
int lucasv(IV* V, IV P, IV Q, UV k)
{
  IV Vl, Vh, Ql, Qh;
  int j, s, n;

  if (V == 0) return 0;
  if (k == 0) { *V = 2; return 1; }

  Vl = 2;  Vh = P;  Ql = 1;  Qh = 1;
  s = 0; n = 0;
  { UV v = k; while (!(v & 1)) { v >>= 1; s++; } }
  { UV v = k; while (v >>= 1) n++; }

  for (j = n; j > s; j--) {
    if (OVERHALF(Vh) || OVERHALF(Vl) || OVERHALF(Ql) || OVERHALF(Qh)) return 0;
    Ql *= Qh;
    if ( (k >> j) & UVCONST(1) ) {
      Qh = Ql * Q;
      Vl = Vh * Vl - P * Ql;
      Vh = Vh * Vh - 2 * Qh;
    } else {
      Qh = Ql;
      Vh = Vh * Vl - P * Ql;
      Vl = Vl * Vl - 2 * Ql;
    }
  }
  if (OVERHALF(Ql) || OVERHALF(Qh)) return 0;
  Ql = Ql * Qh;
  Qh = Ql * Q;
  if (OVERHALF(Vh) || OVERHALF(Vl) || OVERHALF(Ql) || OVERHALF(Qh)) return 0;
  Vl = Vh * Vl - P * Ql;
  Ql = Ql * Qh;
  for (j = 0; j < s; j++) {
    if (OVERHALF(Vl) || OVERHALF(Ql)) return 0;
    Vl = Vl * Vl - 2 * Ql;
    Ql *= Ql;
  }
  *V = Vl;
  return 1;
}

/* Lucas tests:
 *  0: Standard
 *  1: Strong
 *  2: Extra Strong (Mo/Jones/Grantham)
 *
 * None of them have any false positives for the BPSW test.  Also see the
 * "almost extra strong" test.
 */
int is_lucas_pseudoprime(UV n, int strength)
{
  IV P, Q, D;
  UV U, V, Qk, d, s;

  if (n < 7) return (n == 2 || n == 3 || n == 5);
  if ((n % 2) == 0 || n == UV_MAX) return 0;

  if (strength < 2) {
    UV Du = 5;
    IV sign = 1;
    while (1) {
      D = Du * sign;
      if (gcd_ui(Du, n) > 1 && gcd_ui(Du, n) != n) return 0;
      if (jacobi_iu(D, n) == -1)
        break;
      if (Du == 21 && is_perfect_square(n)) return 0;
      Du += 2;
      sign = -sign;
    }
    P = 1;
    Q = (1 - D) / 4;
  } else {
    P = select_extra_strong_parameters(n, 1);
    if (P == 0) return 0;
    Q = 1;
    D = P*P - 4;
    if (gcd_ui(D, n) > 1 && gcd_ui(D, n) != n) return 0;
  }
  MPUassert( D == (P*P - 4*Q) , "is_lucas_pseudoprime: incorrect DPQ");

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
      const uint64_t montQ = (Q >= 0) ? mont_geta(Q, n)
                           : n - mont_geta(-Q, n);
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
      if (_XS_get_verbose()>1) printf("%"UVuf" Q-Type  %"UVuf" -1 %"UVuf"  %"UVuf" 0 %"UVuf"\n", n, A, B, B, C);
      return 1;
    }

  } else {            /* S-Type or I-Type */

    if (jacobi == 0 && n != 23 && restricted > 2) {
      if (_XS_get_verbose()>1) printf("%"UVuf" Jacobi %d\n",n,jacobi);
      return 0;  /* Adams/Shanks allows (-23|n) = 0 for S-Type */
    }

    if (S[0] == 1  &&  S[2] == 3  &&  S[3] == 3  &&  S[5] == 2) {
      if (_XS_get_verbose()>1) printf("%"UVuf" S-Type  1 -1 3  3 0 2\n",n);
      return 1;
    } else if (S[0] == 0  &&  S[5] == n-1  &&  S[2] != S[3]  &&
               addmod(S[2],S[3],n) == n-3  && sqrmod(submod(S[2],S[3],n),n) == n-(23%n)) {
      if (_XS_get_verbose()>1) printf("%"UVuf" I-Type  0 -1 %"UVuf"  %"UVuf" 0 -1\n",n, S[2], S[3]);
      return 1;
    }

  }
  if (_XS_get_verbose()>1) printf("%"UVuf" ? %2d ?  %"UVuf" -1 %"UVuf"  %"UVuf" 0 %"UVuf"\n", n, jacobi, S[0],S[2],S[3],S[5]);
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
    if (_XS_get_verbose()) printf("%"UVuf" Frobenius (%"IVdf",%"IVdf") : x^2 - %"IVdf"x + %"IVdf"\n", n, P, Q, P, Q);
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
  /* if (_XS_get_verbose()) printf("%"UVuf" Frobenius U = %"UVuf" V = %"UVuf"\n", n, U, V); */
  if (U == 0 && V == Vcomp) return 1;
  return 0;
}

/*
 * Khashin, 2013, "Counterexamples for Frobenius primality test"
 * http://arxiv.org/abs/1307.7920
 * 1. Select c as first odd prime where (c,n)=-1.
 * 2. Check (1 + sqrt(c))^n mod n  equiv  (1 - sqrt(c) mod n
 *
 * His Sep 2016 talk starts with c = -1,2 using checks:
 *    (2+sqrt(c)^n = 2-sqrt(c)  mod n   for c = -1,2
 *    (1+sqrt(c)^n = 1-sqrt(c)  mod n   for c = odd prime
 * There doesn't seem to be a big advantage for this change.
 */
int is_frobenius_khashin_pseudoprime(UV n)
{
  int k;
  UV ra, rb, a, b, d = n-1, c = 1;

  if (n < 7) return (n == 2 || n == 3 || n == 5);
  if ((n % 2) == 0 || n == UV_MAX) return 0;
  if (is_perfect_square(n)) return 0;

  /* c = first odd prime where (c|n)=-1 */
  do {
    c += 2;
    if (c==9 || (c>=15 && (!(c%3) || !(c%5) || !(c%7) || !(c%11) || !(c%13))))
      continue;
    k = kronecker_uu(c, n);
  } while (k == 1);
  if (k == 0) return 0;

#if USE_MONTMATH
  {
    const uint64_t npi = mont_inverse(n);
    const uint64_t mont1 = mont_get1(n);
    const uint64_t montc = mont_geta(c, n);
    ra = rb = a = b = mont1;
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
    return (ra == mont1 && rb == n-mont1);
  }
#else
  ra = rb = a = b = 1;
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
  return (ra == 1 && rb == n-1);
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

  if (_XS_get_verbose()>1) printf("%"UVuf" is %s with x = %"UVuf"\n", n, (a == 0 && b == result) ? "probably prime" : "composite", x);
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
#define NUM_KNOWN_MERSENNE_PRIMES 49
static const uint32_t _mersenne_primes[NUM_KNOWN_MERSENNE_PRIMES] = {2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281};
#define LAST_CHECKED_MERSENNE 40270751
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

/* Hashing similar to Forišek and Jančina 2015 */
static const uint16_t mr_bases_hash32[256] = {
157,1150,304,8758,362,15524,1743,212,1056,1607,140,3063,160,913,5842,2013,598,1929,696,1474,3006,524,155,705,694,1238,1851,1053,585,626,603,222,1109,1105,604,646,606,1249,1553,5609,515,548,1371,152,2824,532,3556,831,88,185,1355,501,1556,317,582,4739,4710,145,1045,2976,2674,318,1293,10934,1434,1178,3159,26,3526,1859,6467,602,699,5113,3152,2002,2361,101,464,68,813,446,1368,4637,368,1068,307,2820,6189,10457,569,1690,551,237,226,3235,405,3179,1101,610,56,14647,1687,247,8109,5172,1725,1248,536,2869,1047,899,12285,1026,250,1867,1432,336,5175,1632,5169,39,362,290,1372,11988,1329,2168,34,8781,495,399,34,29,4333,1669,166,6405,7357,694,579,746,1278,6347,7751,179,1085,11734,1615,3575,4253,7894,3097,591,1354,1676,151,702,7,5607,2565,440,566,112,3622,1241,1193,2324,1530,1423,548,3341,2012,6305,2410,39,106,3046,1507,1325,1807,2323,5645,1524,1301,1522,238,1226,2476,2126,1677,3288,1981,18481,287,1011,2877,563,7654,1231,776,3907,117,174,1124,199,16838,164,41,313,1692,1574,1021,2804,1093,1263,956,8508,1221,3743,1318,1304,1344,7628,10739,228,30,520,103,1621,6278,847,4537,272,2213,1989,1826,915,318,401,924,227,911,15505,1670,212,1391,700,3254,4931,3637,2822,1726,137,1843,1300
};


int is_prob_prime(UV n)
{
  int ret;

  if (n < 11) {
    if (n == 2 || n == 3 || n == 5 || n == 7)     return 2;
    else                                          return 0;
  }

#if BITS_PER_WORD == 64
  if (n <= UVCONST(4294967295)) {
#else
  {
#endif
    uint32_t x = n;
    UV base;
    if (!(x%2) || !(x%3) || !(x%5) || !(x%7))       return 0;
    if (x <  121) /* 11*11 */                       return 2;
    if (!(x%11) || !(x%13) || !(x%17) || !(x%19) ||
        !(x%23) || !(x%29) || !(x%31) || !(x%37) ||
        !(x%41) || !(x%43) || !(x%47) || !(x%53))   return 0;
    if (x < 3481) /* 59*59 */                       return 2;
    /* Trial division crossover point depends on platform */
    if (!USE_MONTMATH && n < 500000) {
      uint32_t f = 59;
      uint32_t limit = isqrt(n);
      while (f <= limit) {
        if ((x%f) == 0)  return 0;  f += 2;
        if ((x%f) == 0)  return 0;  f += 6;
        if ((x%f) == 0)  return 0;  f += 4;
        if ((x%f) == 0)  return 0;  f += 2;
        if ((x%f) == 0)  return 0;  f += 4;
        if ((x%f) == 0)  return 0;  f += 2;
        if ((x%f) == 0)  return 0;  f += 4;
        if ((x%f) == 0)  return 0;  f += 6;
      }
      return 2;
    }
    /* Use Mueller-like 32-bit hash to find single M-R base to use. */
    x = (((x >> 16) ^ x) * 0x45d9f3b) & 0xFFFFFFFFUL;
    x = ((x >> 16) ^ x) & 255;
    base = mr_bases_hash32[x];
    ret = miller_rabin(n, &base, 1);
#if BITS_PER_WORD == 64
  } else {  /* input is >= 2^32, UV is 64-bit*/
    if (!(n%2) || !(n%3) || !(n%5) || !(n%7))       return 0;
    if (!(n%11) || !(n%13) || !(n%17) || !(n%19) ||
        !(n%23) || !(n%29) || !(n%31) || !(n%37) ||
        !(n%41) || !(n%43) || !(n%47) || !(n%53))   return 0;
    if (!(n%59) || !(n%61) || !(n%67) || !(n%71))   return 0;
    if (!(n%73) || !(n%79) || !(n%83) || !(n%89))   return 0;
    /*
     * AESLSP test costs about 1.5 Selfridges, vs. ~2.2 for strong Lucas.
     * This makes the full BPSW test cost about 2.5x M-R tests for a prime.
     */
    ret = BPSW(n);
#endif
  }
  return 2*ret;
}
