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

/* Primality related functions, including Montgomery math */

static const UV mr_bases_const2[1] = {2};

/******************************************************************************
  Code inside USE_MONT_PRIMALITY is Montgomery math and efficient M-R from
  Wojciech Izykowski.  See:  https://github.com/wizykowski/miller-rabin

Copyright (c) 2013-2014, Wojciech Izykowski
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * The name of the author may not be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/
#if USE_MONT_PRIMALITY
static INLINE uint64_t mont_prod64(uint64_t a, uint64_t b, uint64_t n, uint64_t npi)
{
  uint64_t t_hi, t_lo, m, mn_hi, mn_lo, u;
  /* t_hi * 2^64 + t_lo = a*b */
  asm("mulq %3" : "=a"(t_lo), "=d"(t_hi) : "a"(a), "rm"(b));
  m = t_lo * npi;
  /* mn_hi * 2^64 + mn_lo = m*n */
  asm("mulq %3" : "=a"(mn_lo), "=d"(mn_hi) : "a"(m), "rm"(n));
  u = t_hi + mn_hi;
  if (t_lo + mn_lo < t_lo) u++;
  return (u < t_hi || u >= n)  ?  u-n  :  u;
}
#define mont_square64(a, n, npi)  mont_prod64(a, a, n, npi)
static INLINE UV mont_powmod64(uint64_t a, uint64_t k, uint64_t one, uint64_t n, uint64_t npi)
{
  uint64_t t = one;
  while (k) {
    if (k & 1) t = mont_prod64(t, a, n, npi);
    k >>= 1;
    if (k)     a = mont_square64(a, n, npi);
  }
  return t;
}
/* Returns -a^-1 mod 2^64.  From B. Arazi "On Primality Testing Using Purely
 * Divisionless Operations", Computer Journal (1994) 37 (3): 219-222, Proc 5 */
static INLINE uint64_t modular_inverse64(const uint64_t a)
{
#if 0
  uint64_t S = 1, J = 0;
  int i;
  for (i = 0; i < 64; i++) {
    if (S & 1) {
      J |= (1ULL << i);
      S += a;
    }
    S >>= 1;
  }
#else
  /* 2 at a time. */
  uint64_t S = 1, J = 0;
  const uint64_t a1p = (a>>1)+1, a2 = a>>2, a2p = a2+1;
  int i;
  if (a & 2) {
    const uint64_t tab[4] = { 0,  a2p,  a1p,  a1p+a2p };
    for (i = 0; i < 32; i++) {
      uint64_t S2 = S >> 2;
      int index = S & 3;
      J |= (uint64_t)index << (2 * i);
      S = S2 + tab[index];
    }
  } else {
    const uint64_t tab[4] = { 0,  a1p+a2,  a1p,  a2p };
    for (i = 0; i < 32; i++) {
      uint64_t S2 = S >> 2;
      int index = S & 3;
      J |= (uint64_t)((4-index)&3) << (2 * i);
      S = S2 + tab[index];
    }
  }
#endif
  return J;
}
static INLINE uint64_t compute_modn64(const uint64_t n)
{

  if (n <= (1ULL << 63)) {
    uint64_t res = ((1ULL << 63) % n) << 1;
    return res < n ? res : res-n;
  } else
    return -n;
}
#define compute_a_times_2_64_mod_n(a, n, r)   mulmod(a, r, n)
static INLINE uint64_t compute_2_65_mod_n(const uint64_t n, const uint64_t modn)
{
  if (n <= (1ULL << 63)) {
    uint64_t res = modn << 1;
    return res < n ? res : res - n;
  } else {
    /* n can fit 2 or 3 times in 2^65 */
    if (n > UVCONST(12297829382473034410))
      return -n-n;    /* 2^65 mod n = 2^65 - 2*n */
    else
      return -n-n-n;  /* 2^65 mod n = 2^65 - 3*n */
  }
}
/* static INLINE int efficient_mr64(const uint64_t bases[], const int cnt, const uint64_t n) */
static int monty_mr64(const uint64_t n, const UV* bases, int cnt)
{
  int i, j, t;
  const uint64_t npi = modular_inverse64(n);
  const uint64_t r = compute_modn64(n);
  uint64_t u = n - 1;
  const uint64_t nr = n - r;

  t = 0;
  while (!(u&1)) {  t++;  u >>= 1;  }
  for (j = 0; j < cnt; j++) {
    const uint64_t a = bases[j];
    uint64_t       A = compute_a_times_2_64_mod_n(a, n, r);
    uint64_t       d;
    if (a < 2)  croak("Base %"UVuf" is invalid", (UV)a);
    if (!A) continue;  /* PRIME in subtest */
    d = mont_powmod64(A, u, r, n, npi);  /* compute a^u mod n */
    if (d == r || d == nr) continue;  /* PRIME in subtest */
    for (i=1; i<t; i++) {
      d = mont_square64(d, n, npi);
      if (d == r) return 0;
      if (d == nr) break;  /* PRIME in subtest */
    }
    if (i == t) return 0;
  }
  return 1;
}
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




/* Fermat pseudoprime */
int _XS_is_pseudoprime(UV const n, UV a)
{
  UV x;

  if (n < 5) return (n == 2 || n == 3);
  if (a < 2) croak("Base %"UVuf" is invalid", a);
  if (a >= n) {
    a %= n;
    if ( a <= 1 || a == n-1 )
      return 1;
  }
  x = powmod(a, n-1, n);    /* x = a^(n-1) mod n */
  return (x == 1);
}

/* Miller-Rabin probabilistic primality test
 * Returns 1 if probably prime relative to the bases, 0 if composite.
 * Bases must be between 2 and n-2
 */
int _XS_miller_rabin(UV const n, const UV *bases, int nbases)
{
  UV d = n-1;
  int b, r, s = 0;

  MPUassert(n > 3, "MR called with n <= 3");
  if ((n & 1) == 0) return 0;

#if USE_MONT_PRIMALITY
  return monty_mr64((uint64_t)n, bases, nbases);
#else
  while (!(d&1)) {  s++;  d >>= 1;  }

  for (b = 0; b < nbases; b++) {
    UV x, a = bases[b];
    if (a < 2)  croak("Base %"UVuf" is invalid", a);
    if (a >= n) a %= n;
    if ( (a <= 1) || (a == n-1) )
      continue;
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

int _XS_BPSW(UV const n)
{
  if (n < 7) return (n == 2 || n == 3 || n == 5);
  if ((n % 2) == 0 || n == UV_MAX) return 0;

#if !USE_MONT_PRIMALITY
  return    _XS_miller_rabin(n, mr_bases_const2, 1)
         && _XS_is_almost_extra_strong_lucas_pseudoprime(n,1);
#else
  {
    const uint64_t npi = modular_inverse64(n);
    const uint64_t montr = compute_modn64(n);
    const uint64_t mont2 = compute_2_65_mod_n(n, montr);
    uint64_t u = n-1;
    const uint64_t nr = n-montr;
    int i, t = 0;
    UV P, V, d, s;

    /* M-R with base 2 */
    while (!(u&1)) {  t++;  u >>= 1;  }
    {
      uint64_t A = mont2;
      if (A) {
        uint64_t d = mont_powmod64(A, u, montr, n, npi);
        if (d != montr && d != nr) {
          for (i=1; i<t; i++) {
            d = mont_square64(d, n, npi);
            if (d == montr) return 0;
            if (d == nr) break;
          }
          if (i == t)
            return 0;
        }
      }
    }
    /* AES Lucas test */
    P = 3;
    while (1) {
      UV D = P*P - 4;
      d = gcd_ui(D, n);
      if (d > 1 && d < n)
        return 0;
      if (jacobi_iu(D, n) == -1)
        break;
      if (P == (3+20) && is_perfect_square(n)) return 0;
      P++;
      if (P > 65535)
        croak("lucas_extrastrong_params: P exceeded 65535");
    }

    d = n+1;
    s = 0;
    while ( (d & 1) == 0 ) {  s++;  d >>= 1; }

    {
      const uint64_t montP = compute_a_times_2_64_mod_n(P, n, montr);
      UV W, b;
      W = submod(  mont_prod64( montP, montP, n, npi),  mont2, n);
      V = montP;
      { UV v = d; b = 1; while (v >>= 1) b++; }
      while (b-- > 1) {
        UV T = submod(  mont_prod64(V, W, n, npi),  montP, n);
        if ( (d >> (b-1)) & UVCONST(1) ) {
          V = T;
          W = submod(  mont_prod64(W, W, n, npi),  mont2, n);
        } else {
          W = T;
          V = submod(  mont_prod64(V, V, n, npi),  mont2, n);
        }
      }
    }

    if (V == mont2 || V == (n-mont2))
      return 1;
    while (s-- > 1) {
      if (V == 0)
        return 1;
      V = submod(  mont_prod64(V, V, n, npi),  mont2, n);
      if (V == mont2)
        return 0;
    }
  }
  return 0;
#endif
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

  Qmod = (Q < 0)  ?  (UV)(Q + (IV)n)  :  (UV)Q;
  Pmod = (P < 0)  ?  (UV)(P + (IV)n)  :  (UV)P;
  Dmod = submod( mulmod(Pmod, Pmod, n), mulmod(4, Qmod, n), n);
  if (Dmod == 0) {
    b = Pmod >> 1;
    *Uret = mulmod(k, powmod(b, k-1, n), n);
    *Vret = mulmod(2, powmod(b, k, n), n);
    *Qkret = powmod(Qmod, k, n);
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

/* Lucas tests:
 *  0: Standard
 *  1: Strong
 *  2: Extra Strong (Mo/Jones/Grantham)
 *
 * None of them have any false positives for the BPSW test.  Also see the
 * "almost extra strong" test.
 */
int _XS_is_lucas_pseudoprime(UV n, int strength)
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
    P = 3;
    Q = 1;
    while (1) {
      D = P*P - 4;
      if (gcd_ui(D, n) > 1 && gcd_ui(D, n) != n) return 0;
      if (jacobi_iu(D, n) == -1)
        break;
      if (P == 21 && is_perfect_square(n)) return 0;
      P++;
    }
  }
  MPUassert( D == (P*P - 4*Q) , "is_lucas_pseudoprime: incorrect DPQ");

  d = n+1;
  s = 0;
  if (strength > 0)
    while ( (d & 1) == 0 ) {  s++;  d >>= 1; }

#if USE_MONT_PRIMALITY
  {
    const uint64_t npi = modular_inverse64(n);
    const uint64_t mont1 = compute_modn64(n);
    const uint64_t mont2 = compute_2_65_mod_n(n, mont1);
    const uint64_t montP = (P == 1) ? mont1
                         : (P >= 0) ? compute_a_times_2_64_mod_n(P, n, mont1)
                         : n - compute_a_times_2_64_mod_n(-P, n, mont1);
    const uint64_t montD = (D >= 0) ? compute_a_times_2_64_mod_n(D, n, mont1)
                         : n - compute_a_times_2_64_mod_n(-D, n, mont1);
    UV b;
    { UV v = d; b = 0; while (v >>= 1) b++; }

    /* U, V, Qk, and mont* are in Montgomery space */
    U = mont1;
    V = montP;

    if (Q == 1 || Q == -1) {   /* Faster code for |Q|=1, also opt for P=1 */
      int sign = Q;
      while (b--) {
        U = mont_prod64(U, V, n, npi);
        if (sign == 1) V = submod( mont_square64(V,n,npi), mont2, n);
        else           V = addmod( mont_square64(V,n,npi), mont2, n);
        sign = 1;
        if ( (d >> b) & UVCONST(1) ) {
          UV t2 = mont_prod64(U, montD, n, npi);
          if (P == 1) {
            U = addmod(U, V, n);
            V = addmod(V, t2, n);
          } else {
            U = addmod( mont_prod64(U, montP, n, npi), V, n);
            V = addmod( mont_prod64(V, montP, n, npi), t2, n);
          }
          if (U & 1) { U = (n>>1) + (U>>1) + 1; } else { U >>= 1; }
          if (V & 1) { V = (n>>1) + (V>>1) + 1; } else { V >>= 1; }
          sign = Q;
        }
      }
      Qk = (sign == 1) ? mont1 : n-mont1;
    } else {
      const uint64_t montQ = (Q >= 0) ? compute_a_times_2_64_mod_n(Q, n, mont1)
                           : n - compute_a_times_2_64_mod_n(-Q, n, mont1);
      Qk = montQ;
      while (b--) {
        U = mont_prod64(U, V, n, npi);
        V = submod( mont_square64(V,n,npi), addmod(Qk,Qk,n), n);
        Qk = mont_square64(Qk,n,npi);
        if ( (d >> b) & UVCONST(1) ) {
          UV t2 = mont_prod64(U, montD, n, npi);
          U = addmod( mont_prod64(U, montP, n, npi), V, n);
          if (U & 1) { U = (n>>1) + (U>>1) + 1; } else { U >>= 1; }
          V = addmod( mont_prod64(V, montP, n, npi), t2, n);
          if (V & 1) { V = (n>>1) + (V>>1) + 1; } else { V >>= 1; }
          Qk = mont_prod64(Qk, montQ, n, npi);
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
          V = submod( mont_square64(V,n,npi), addmod(Qk,Qk,n), n);
          Qk = mont_square64(Qk,n,npi);
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
          V = submod( mont_square64(V,n,npi), mont2, n);
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
int _XS_is_almost_extra_strong_lucas_pseudoprime(UV n, UV increment)
{
  UV P, V, W, d, s, b;

  if (n < 7) return (n == 2 || n == 3 || n == 5);
  if ((n % 2) == 0 || n == UV_MAX) return 0;
  if (increment < 1 || increment > 256)
    croak("Invalid lucas parameter increment: %"UVuf"\n", increment);

  P = 3;
  while (1) {
    UV D = P*P - 4;
    d = gcd_ui(D, n);
    if (d > 1 && d < n)
      return 0;
    if (jacobi_iu(D, n) == -1)
      break;
    if (P == (3+20*increment) && is_perfect_square(n)) return 0;
    P += increment;
    if (P > 65535)
      croak("lucas_extrastrong_params: P exceeded 65535");
  }
  if (P >= n)  P %= n;   /* Never happens with increment < 4 */

  d = n+1;
  s = 0;
  while ( (d & 1) == 0 ) {  s++;  d >>= 1; }
  { UV v = d; b = 0; while (v >>= 1) b++; }

#if USE_MONT_PRIMALITY
  {
    const uint64_t npi = modular_inverse64(n);
    const uint64_t montr = compute_modn64(n);
    const uint64_t mont2 = compute_2_65_mod_n(n, montr);
    const uint64_t montP = compute_a_times_2_64_mod_n(P, n, montr);
    W = submod(  mont_prod64( montP, montP, n, npi),  mont2, n);
    V = montP;
    while (b--) {
      UV T = submod(  mont_prod64(V, W, n, npi),  montP, n);
      if ( (d >> b) & UVCONST(1) ) {
        V = T;
        W = submod(  mont_prod64(W, W, n, npi),  mont2, n);
      } else {
        W = T;
        V = submod(  mont_prod64(V, V, n, npi),  mont2, n);
      }
    }

    if (V == mont2 || V == (n-mont2))
      return 1;
    s--;
    while (s--) {
      if (V == 0)
        return 1;
      if (s)
        V = submod(  mont_prod64(V, V, n, npi),  mont2, n);
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

/*
 * The Frobenius-Underwood test has no known counterexamples below 10^13, but
 * has not been extensively tested above that.  This is the Minimal Lambda+2
 * test from section 9 of "Quadratic Composite Tests" by Paul Underwood.
 *
 * It is generally slower than the AES Lucas test, but for large values is
 * competitive with the BPSW test.  Since our BPSW is known to have no
 * counterexamples under 2^64, while the results of this test are unknown,
 * it is mainly useful for numbers larger than 2^64 as an additional
 * non-correlated test.
 */
int _XS_is_frobenius_underwood_pseudoprime(UV n)
{
  int bit;
  UV x, result, a, b, np1, len, t1;
  IV t;

  if (n < 7) return (n == 2 || n == 3 || n == 5);
  if ((n % 2) == 0 || n == UV_MAX) return 0;
  if (is_perfect_square(n)) return 0;

  x = 0;
  t = -1;
  while ( jacobi_iu( t, n ) != -1 ) {
    x++;
    t = (IV)(x*x) - 4;
  }
  np1 = n+1;
  { UV v = np1; len = 1;  while (v >>= 1) len++; }

#if USE_MONT_PRIMALITY
  {
    const uint64_t npi = modular_inverse64(n);
    const uint64_t mont1 = compute_modn64(n);
    const uint64_t mont2 = compute_2_65_mod_n(n, mont1);
    const uint64_t mont5 = compute_a_times_2_64_mod_n(5, n, mont1);

    x = compute_a_times_2_64_mod_n(x, n, mont1);
    a = mont1;
    b = mont2;

    if (x == 0) {
      result = mont5;
      for (bit = len-2; bit >= 0; bit--) {
        t1 = addmod(b, b, n);
        b = mont_prod64(submod(b, a, n), addmod(b, a, n), n, npi);
        a = mont_prod64(a, t1, n, npi);
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
        t1 = addmod( mont_prod64(a, x, n, npi), addmod(b, b, n), n);
        b = mont_prod64(submod(b, a, n), addmod(b, a, n), n, npi);
        a = mont_prod64(a, t1, n, npi);
        if ( (np1 >> bit) & UVCONST(1) ) {
          t1 = b;
          b = submod( addmod(b, b, n), a, n);
          a = addmod( mont_prod64(a, multiplier, n, npi), t1, n);
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


/******************************************************************************/


/* Select M-R bases from http://miller-rabin.appspot.com/, 26 July 2013 */
#if BITS_PER_WORD == 32
static const UV mr_bases_small_2[2] = {31, 73};
static const UV mr_bases_small_3[3] = {2, 7, 61};
#else
static const UV mr_bases_large_1[1] = { UVCONST(  9345883071009581737 ) };
static const UV mr_bases_large_2[2] = { UVCONST(         336781006125 ),
                                        UVCONST(     9639812373923155 ) };
static const UV mr_bases_large_3[3] = { UVCONST(  4230279247111683200 ),
                                        UVCONST( 14694767155120705706 ),
                                        UVCONST( 16641139526367750375 ) };
static const UV mr_bases_large_7[7] = { 2, 325, 9375, 28178, 450775, 9780504, 1795265022 };
#endif

int is_prob_prime(UV n)
{
  int ret;

  if (n < 11) {
    if (n == 2 || n == 3 || n == 5 || n == 7)     return 2;
    else                                          return 0;
  }
  if (!(n%2) || !(n%3) || !(n%5) || !(n%7))       return 0;
  if (n <  121) /* 11*11 */                       return 2;
  if (!(n%11) || !(n%13) || !(n%17) || !(n%19) ||
      !(n%23) || !(n%29) || !(n%31) || !(n%37) ||
      !(n%41) || !(n%43) || !(n%47) || !(n%53))   return 0;
  if (n < 3481) /* 59*59 */                       return 2;

#if BITS_PER_WORD == 32
  /* We could use one base when n < 49191, two when n < 360018361. */
  if (n < UVCONST(9080191))
    ret = _XS_miller_rabin(n, mr_bases_small_2, 2);
  else
    ret = _XS_miller_rabin(n, mr_bases_small_3, 3);
#else
  /* AESLSP test costs about 1.5 Selfridges, vs. ~2.2 for strong Lucas.
   * So it works out to be faster to do AES-BPSW vs. 3 M-R tests. */
  if (n < UVCONST(341531))
    ret = _XS_miller_rabin(n, mr_bases_large_1, 1);
  else if (n < UVCONST(1050535501))
    ret = _XS_miller_rabin(n, mr_bases_large_2, 2);
  else
    ret = _XS_BPSW(n);
  /*
    ret = efficient_mr64(mr_bases_large_7, 7, n);
    ret = _XS_miller_rabin(n, mr_bases_large_7, 7);
  */
#endif
  return 2*ret;
}
