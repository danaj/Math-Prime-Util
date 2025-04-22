#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ptypes.h"
#include "lucas_seq.h"
#include "mulmod.h"
#include "util.h"

/* TODO:  primality.c might be able to call these more often */

/* TODO:  montmath version of fastest, e.g. lucasvmod */



/******************************************************************************/

/* Alternate modular lucas sequence code.
 * A bit slower than the normal one, but works with even valued n. */
static void alt_lucas_seq(UV* Uret, UV* Vret, UV n, UV Pmod, UV Qmod, UV k)
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
    return;
  }

  for (j = m; j > s; j--) {
    Ql = mulmod(Ql, Qh, n);
    if ( (k >> j) & UVCONST(1) ) {
      Qh = mulmod(Ql, Qmod, n);
      Uh = mulmod(Uh, Vh, n);
      Vl = submod(mulmod(Vh, Vl, n), mulmod(Pmod, Ql, n), n);
      Vh = submod(sqrmod(Vh, n), addmod(Qh,Qh,n), n);
    } else {
      Qh = Ql;
      Uh = submod(mulmod(Uh, Vl, n), Ql, n);
      Vh = submod(mulmod(Vh, Vl, n), mulmod(Pmod, Ql, n), n);
      Vl = submod(sqrmod(Vl, n), addmod(Ql, Ql, n), n);
    }
  }
  Ql = mulmod(Ql, Qh, n);
  Qh = mulmod(Ql, Qmod, n);
  Uh = submod(mulmod(Uh, Vl, n), Ql, n);
  Vl = submod(mulmod(Vh, Vl, n), mulmod(Pmod, Ql, n), n);
  Ql = mulmod(Ql, Qh, n);
  for (j = 0; j < s; j++) {
    Uh = mulmod(Uh, Vl, n);
    Vl = submod(sqrmod(Vl, n), addmod(Ql, Ql, n), n);
    Ql = sqrmod(Ql, n);
  }
  *Uret = Uh;
  *Vret = Vl;
}

/* Generic Lucas sequence for any appropriate P and Q */
void lucas_seq(UV* Uret, UV* Vret, UV* Qkret, UV n, IV P, IV Q, UV k)
{
  MPUassert(n > 0, "lucas_sequence:  modulus n must be > 0");
  if (n == 1) { *Uret = *Vret = *Qkret = 0; return; }

  lucasuvmod(Uret, Vret, ivmod(P,n), ivmod(Q,n), k, n);
  *Qkret = powmod(ivmod(Q,n), k, n);
}


void lucasuvmod(UV* Uret, UV* Vret, UV P, UV Q, UV k, UV n)
{
  UV U, V, b, D, invD;

  MPUassert(n > 0, "lucasuvmod:  modulus n must be > 0");
  if (n == 1) { *Uret = *Vret = 0; return; }

  if (k == 0) {
    *Uret = 0;
    *Vret = 2 % n;
    return;
  }

  if (P >= n) P %= n;
  if (Q >= n) Q %= n;
  D = submod( mulmod(P, P, n), mulmod(4, Q, n), n);

  if (D == 0 && (b = divmod(P,2,n)) != 0) {
    *Uret = mulmod(k, powmod(b, k-1, n), n);
    *Vret = mulmod(2, powmod(b, k, n), n);
    return;
  }

  { UV v = k; b = 0; while (v >>= 1) b++; }
  U = 1;
  V = P;
  invD = modinverse(D, n);

  if (Q == 1 && invD != 0) { /* Inverting D: 2 mulmods/bit instead of 2-5 */
    U = mulsubmod(P,P,2,n);
    while (b--) {
      UV T = mulsubmod(U, V, P, n);
      if ( (k >> b) & UVCONST(1) ) {
        V = T;
        U = mulsubmod(U, U, 2, n);
      } else {
        U = T;
        V = mulsubmod(V, V, 2, n);
      }
    }
    U = addmod(U,U,n);
    U = submod(U, mulmod(V,P,n), n);
    U = mulmod(U, invD, n);
  } else if (P == 1 && Q == (n-1)) {  /* code for P=1 Q=-1 in here */
    alt_lucas_seq(&U, &V, n, P, Q, k);
  } else if ((n & 1) && (Q == 1 || (Q == (n-1)))) {
    int qs = (Q==1);
    while (b--) {
      U = mulmod(U, V, n);
      V = muladdmod(V, V, (qs) ? n-2 : 2, n);
      qs = 1;
      if ( (k >> b) & UVCONST(1) ) {
        UV t2 = mulmod(U, D, n);
        if (P != 1) U = mulmod(U, P, n);
        U = addmod(U, V, n);
        if (U & 1) { U = (n>>1) + (U>>1) + 1; } else { U >>= 1; }
        if (P != 1) V = mulmod(V, P, n);
        V = addmod(V, t2, n);
        if (V & 1) { V = (n>>1) + (V>>1) + 1; } else { V >>= 1; }
        qs = (Q==1);
      }
    }
  } else if (n & 1) {
    UV Qk = Q;
    while (b--) {
      U = mulmod(U, V, n);
      V = mulsubmod(V, V, addmod(Qk,Qk,n), n);
      Qk = sqrmod(Qk, n);
      if ( (k >> b) & UVCONST(1) ) {
        UV t2 = mulmod(U, D, n);
        U = muladdmod(U, P, V, n);
        if (U & 1) { U = (n>>1) + (U>>1) + 1; } else { U >>= 1; }
        V = muladdmod(V, P, t2, n);
        if (V & 1) { V = (n>>1) + (V>>1) + 1; } else { V >>= 1; }
        Qk = mulmod(Qk, Q, n);
      }
    }
  } else {
    /* This handles everything */
    alt_lucas_seq(&U, &V, n, P, Q, k);
  }
  *Uret = U;
  *Vret = V;
}

UV lucasvmod(UV P, UV Q, UV k, UV n)
{
  UV D, b, U, V, Qk;

  MPUassert(n > 0, "lucas_sequence:  modulus n must be > 0");
  if (n == 1) return 0;
  if (k == 0) return 2 % n;
  if (P >= n) P = P % n;
  if (Q >= n) Q = Q % n;

  D = submod(mulmod(P, P, n), mulmod(4, Q, n), n);
  if (D == 0 && (b = divmod(P,2,n)) != 0)
    return mulmod(2, powmod(b, k, n), n);

  { UV v = k; b = 0; while (v >>= 1) b++; }

  if (Q == 1) {

    V = P;
    U = mulsubmod(P, P, 2, n);
    while (b--) {
      UV T = mulsubmod(U, V, P, n);
      if ( (k >> b) & UVCONST(1) ) {
        V = T;
        U = mulsubmod(U, U, 2, n);
      } else {
        U = T;
        V = mulsubmod(V, V, 2, n);
      }
    }

  } else if ((n % 2) == 0) {

    alt_lucas_seq(&U, &V, n, P, Q, k);

  } else {

    U = 1;
    V = P;
    Qk = Q;
    while (b--) {
      U = mulmod(U, V, n);
      V = mulsubmod(V, V, addmod(Qk,Qk,n), n);
      Qk = sqrmod(Qk, n);
      if ( (k >> b) & UVCONST(1) ) {
        UV t2 = mulmod(U, D, n);
        U = muladdmod(U, P, V, n);
        if (U & 1) { U = (n>>1) + (U>>1) + 1; } else { U >>= 1; }
        V = muladdmod(V, P, t2, n);
        if (V & 1) { V = (n>>1) + (V>>1) + 1; } else { V >>= 1; }
        Qk = mulmod(Qk, Q, n);
      }
    }

  }
  return V;
}

UV lucasumod(UV P, UV Q, UV k, UV n)
{
  UV U, V;
  lucasuvmod(&U, &V, P, Q, k, n);
  return U;
}




#define OVERHALF(v)  ( (UV)((v>=0)?v:-v) > (UVCONST(1) << (BITS_PER_WORD/2-1)) )
bool lucasuv(IV* U, IV *V, IV P, IV Q, UV k)
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
