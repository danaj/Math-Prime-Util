#ifndef MPU_MULMOD_H
#define MPU_MULMOD_H

#include "ptypes.h"

/* if n is smaller than this, you can multiply without overflow */
#define HALF_WORD (UVCONST(1) << (BITS_PER_WORD/2))
/* This will be true if we think mulmods are fast */
#define MULMODS_ARE_FAST 1

#if (BITS_PER_WORD == 32) && HAVE_STD_U64

  /* We have 64-bit available, but UV is 32-bit.  Do the math in 64-bit.
   * Even if it is emulated, it should be as fast or faster than us doing it.
   */
  #define addmod(a,b,n)  (UV)( ((uint64_t)(a) + (b)) % (n) )
  #define mulmod(a,b,n)  (UV)( ((uint64_t)(a) * (b)) % (n) )
  #define sqrmod(a,n)    (UV)( ((uint64_t)(a) * (a)) % (n) )

#elif defined(__GNUC__) && defined(__x86_64__)

  /* GCC on a 64-bit Intel x86, help from WraithX and Wojciech Izykowski */
  /* Beware: if (a*b)/c > 2^64, there will be an FP exception */
  static INLINE UV _mulmod(UV a, UV b, UV n) {
    UV d, dummy;                    /* d will get a*b mod c */
    asm ("mulq %3\n\t"              /* mul a*b -> rdx:rax */
         "divq %4\n\t"              /* (a*b)/c -> quot in rax remainder in rdx */
         :"=a"(dummy), "=&d"(d)     /* output */
         :"a"(a), "rm"(b), "rm"(n)  /* input */
         :"cc"                      /* mulq and divq can set conditions */
        );
    return d;
  }
  #define mulmod(a,b,n) _mulmod(a,b,n)
  #define sqrmod(a,n)   _mulmod(a,a,n)
  /* A version for _MSC_VER:
   *    __asm { mov rax, qword ptr a
   *            mul qword ptr b
   *            div qword ptr c
   *            mov qword ptr d, rdx }   */

  /* addmod from Kruppa 2010 page 67 */
  static INLINE UV _addmod(UV a, UV b, UV n) {
    UV r = a+b;
    UV t = a-n;
    asm ("add %2, %1\n\t"    /* t := t + b */
         "cmovc %1, %0\n\t"  /* if (carry) r := t */
         :"+r" (r), "+&r" (t)
         :"rm" (b)
         :"cc"
        );
    return r;
  }
  #define addmod(a,b,n) _addmod(a,b,n)

#elif BITS_PER_WORD == 64 && HAVE_UINT128

  /* We're 64-bit, using a modern gcc, and the target has some 128-bit type.
   * The actual number of targets that have this implemented are limited. */

  #define mulmod(a,b,n) (UV)( ((uint128_t)(a) * (b)) % (n) )
  #define sqrmod(a,n)   (UV)( ((uint128_t)(a) * (a)) % (n) )

#else

  /* UV is the largest integral type available (that we know of). */
#undef MULMODS_ARE_FAST
#define MULMODS_ARE_FAST 0

  /* Do it by hand */
  static INLINE UV _mulmod(UV a, UV b, UV n) {
    UV r = 0;
    if (a >= n) a %= n;   /* Careful attention from the caller should make */
    if (b >= n) b %= n;   /* these unnecessary.                            */
    if ((a|b) < HALF_WORD) return (a*b) % n;
    if (a < b) { UV t = a; a = b; b = t; }
    if (n <= (UV_MAX>>1)) {
      while (b > 0) {
        if (b & 1)  { r += a;  if (r >= n) r -= n; }
        b >>= 1;
        if (b)      { a += a;  if (a >= n) a -= n; }
      }
    } else {
      while (b > 0) {
        if (b & 1)  r = ((n-r) > a) ? r+a : r+a-n;    /* r = (r + a) % n */
        b >>= 1;
        if (b)      a = ((n-a) > a) ? a+a : a+a-n;    /* a = (a + a) % n */
      }
    }
    return r;
  }

  #define mulmod(a,b,n) _mulmod(a,b,n)
  #define sqrmod(a,n)   _mulmod(a,a,n)

#endif

#ifndef addmod
  static INLINE UV addmod(UV a, UV b, UV n) {
    return ((n-a) > b) ?  a+b  :  a+b-n;
  }
#endif

static INLINE UV submod(UV a, UV b, UV n) {
  UV t = n-b;  /* Evaluate as UV, then hand to addmod */
  return addmod(a, t, n);
}

/* a^2 + c mod n */
#define sqraddmod(a, c, n)     addmod(sqrmod(a,n),   c, n)
/* a*b + c mod n */
#define muladdmod(a, b, c, n)  addmod(mulmod(a,b,n), c, n)
/* a*b - c mod n */
#define mulsubmod(a, b, c, n)  submod(mulmod(a,b,n), c, n)

/* a^k mod n */
#ifndef HALF_WORD
  static INLINE UV powmod(UV a, UV k, UV n) {
    UV t = 1;
    if (a >= n) a %= n;
    while (k) {
      if (k & 1) t = mulmod(t, a, n);
      k >>= 1;
      if (k)     a = sqrmod(a, n);
    }
    return t;
  }
#else
  static INLINE UV powmod(UV a, UV k, UV n) {
    UV t = 1;
    if (a >= n) a %= n;
    if (n < HALF_WORD) {
      while (k) {
        if (k & 1) t = (t*a)%n;
        k >>= 1;
        if (k)     a = (a*a)%n;
      }
    } else {
      while (k) {
        if (k & 1) t = mulmod(t, a, n);
        k >>= 1;
        if (k)     a = sqrmod(a, n);
      }
    }
    return t;
  }
#endif

/* a^k + c mod n */
#define powaddmod(a, k, c, n)  addmod(powmod(a,k,n),c,n)

#endif
