#ifndef MPU_MONTMATH_H
#define MPU_MONTMATH_H

#include "ptypes.h"
#include "mulmod.h"

#if BITS_PER_WORD == 64 && HAVE_STD_U64 && defined(__GNUC__) && defined(__x86_64__)
#define USE_MONTMATH 1
#else
#define USE_MONTMATH 0
#endif

#if USE_MONTMATH
/******************************************************************************/
/*                   This is the interface we'll use                          */
/******************************************************************************/
#define mont_inverse(n)           modular_inverse64(n)
#define mont_get1(n)              compute_modn64(n)
/* Must have npi = mont_inverse(n), mont1 = mont_get1(n) */
#define mont_get2(n)              addmod(mont1,mont1,n)
#define mont_geta(a,n)            mulmod(a,mont1,n)
#define mont_mulmod(a,b,n)        mont_prod64(a,b,n,npi)
#define mont_sqrmod(a,n)          mont_prod64(a,a,n,npi)
#define mont_powmod(a,k,n)        mont_powmod64(a,k,mont1,n,npi)
#define mont_recover(a,n)         mont_mulmod(a,1,n)
/******************************************************************************/
/******************************************************************************/



/******************************************************************************
  Montgomery math and efficient M-R from
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

#if defined(__GNUC__) && (__GNUC__ >= 3)
 #define MPU_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
 #define MPU_UNLIKELY(x) (x)
#endif

static INLINE uint64_t mont_prod64(uint64_t a, uint64_t b, uint64_t n, uint64_t npi)
{
  uint64_t t_hi, t_lo, m, mn_hi, mn_lo, u;
  /* t_hi * 2^64 + t_lo = a*b */
  asm("mulq %3" : "=a"(t_lo), "=d"(t_hi) : "a"(a), "rm"(b));
  if (MPU_UNLIKELY(t_lo == 0)) return t_hi;  /* Allows carry=1 below */
  m = t_lo * npi;
  /* mn_hi * 2^64 + mn_lo = m*n */
  asm("mulq %3" : "=a"(mn_lo), "=d"(mn_hi) : "a"(m), "rm"(n));
  u = t_hi + mn_hi + 1;
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
 * Divisionless Operations", Computer Journal (1994) 37 (3): 219-222, Proc 5.
 * Hensel lifting idea from Robert Gerbicz. */
static INLINE uint64_t modular_inverse64(const uint64_t a)
{
  static const char mask[128] = {255,85,51,73,199,93,59,17,15,229,195,89,215,237,203,33,31,117,83,105,231,125,91,49,47,5,227,121,247,13,235,65,63,149,115,137,7,157,123,81,79,37,3,153,23,45,11,97,95,181,147,169,39,189,155,113,111,69,35,185,55,77,43,129,127,213,179,201,71,221,187,145,143,101,67,217,87,109,75,161,159,245,211,233,103,253,219,177,175,133,99,249,119,141,107,193,191,21,243,9,135,29,251,209,207,165,131,25,151,173,139,225,223,53,19,41,167,61,27,241,239,197,163,57,183,205,171,1};
  uint64_t ret = mask[(a >> 1) & 127];
  ret *= 2 + a * ret;
  ret *= 2 + a * ret;
  ret *= 2 + a * ret;
  return ret;
}
static INLINE uint64_t compute_modn64(const uint64_t n)
{

  if (n <= (1ULL << 63)) {
    uint64_t res = ((1ULL << 63) % n) << 1;
    return res < n ? res : res-n;
  } else
    return -n;
}
#ifdef FUNC_monty_mr64
#define compute_a_times_2_64_mod_n(a, n, r)   mulmod(a, r, n)
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
    uint64_t d, a = bases[j];
    uint64_t A = compute_a_times_2_64_mod_n(a, n, r);
    if (a < 2)  croak("Base %"UVuf" is invalid", (UV)a);
    if (a >= n) {
      a %= n;
      if (a == 0 || (a == n-1 && a&1)) return 0;
    }
    if (a == 1 || a == n-1) continue;
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
#endif /* use_montmath */

#endif
