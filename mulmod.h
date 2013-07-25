#ifndef MPU_MULMOD_H
#define MPU_MULMOD_H

#include "ptypes.h"

#if defined(__GNUC__)
  #define INLINE inline
#elif defined(_MSC_VER)
  #define INLINE __inline
#else
  #define INLINE
#endif

/* if n is smaller than this, you can multiply without overflow */
#define HALF_WORD (UVCONST(1) << (BITS_PER_WORD/2))

#if (BITS_PER_WORD == 32) && HAVE_STD_U64

  /* We have 64-bit available, but UV is 32-bit.  Do the math in 64-bit.
   * Even if it is emulated, it should be as fast or faster than us doing it.
   */
  #define addmod(n,a,m)  (UV)(((uint64_t)(n)+(uint64_t)(a)) % ((uint64_t)(m)))
  #define mulmod(a,b,m)  (UV)(((uint64_t)(a)*(uint64_t)(b)) % ((uint64_t)(m)))
  #define sqrmod(n,m)    (UV)(((uint64_t)(n)*(uint64_t)(n)) % ((uint64_t)(m)))

#elif defined(__GNUC__) && defined(__x86_64__)

  /* GCC on a 64-bit Intel x86 */
  static INLINE UV _mulmod(UV a, UV b, UV c) {
    UV d; /* to hold the result of a*b mod c */
    /* calculates a*b mod c, stores result in d */
    asm ("mov %1, %%rax;"        /* put a into rax */
         "mul %2;"               /* mul a*b -> rdx:rax */
         "div %3;"               /* (a*b)/c -> quot in rax remainder in rdx */
         "mov %%rdx, %0;"        /* store result in d */
         :"=r"(d)                /* output */
         :"r"(a), "r"(b), "r"(c) /* input */
         :"%rax", "%rdx"         /* clobbered registers */
        );
    /* A version for _MSC_VER:
     *
     *    __asm { mov rax, qword ptr a
     *            mul qword ptr b
     *            div qword ptr c
     *            mov qword ptr d, rdx }
     */
    return d;
  }
  #define mulmod(a,b,m) _mulmod(a,b,m)
  #define sqrmod(n,m)   _mulmod(n,n,m)

#elif BITS_PER_WORD == 64 && __GNUC__ == 4 && __GNUC_MINOR__ >= 4 && (defined(__LP64__) || defined(__x86_64__) || defined(__powerpc64__) || defined(_M_X64) || defined(_M_IX86))

  /* We're 64-bit, using a modern gcc, and the target has some 128-bit type */

  #if __GNUC__ == 4 && __GNUC_MINOR__ >= 4 && __GNUC_MINOR__ < 6
    typedef unsigned int uint128_t __attribute__ ((__mode__ (TI)));
  #else
    typedef unsigned __int128 uint128_t;
  #endif

  #define mulmod(a,b,m) (UV)(((uint128_t)(a)*(uint128_t)(b)) % ((uint128_t)(m)))
  #define sqrmod(n,m)   (UV)(((uint128_t)(n)*(uint128_t)(n)) % ((uint128_t)(m)))

#else

  /* UV is the largest integral type available (that we know of). */

  /* Do it by hand */
  static INLINE UV _mulmod(UV a, UV b, UV m) {
    UV r = 0;
    if (a >= m) a %= m;   /* Careful attention from the caller should make */
    if (b >= m) b %= m;   /* these unnecessary.                            */
    if (a < b) { UV t = a; a = b; b = t; }
    while (b > 0) {
      if (b & 1)  r = ((m-r) > a) ? r+a : r+a-m;    /* r = (r + a) % m */
      b >>= 1;
      if (b)      a = ((m-a) > a) ? a+a : a+a-m;    /* a = (a + a) % m */
    }
    return r;
  }

  #define mulmod(a,b,m) (((a)|(b)) < HALF_WORD) ? ((a)*(b))%(m):_mulmod(a,b,m)
  #define sqrmod(n,m)   ((n) < HALF_WORD)       ? ((n)*(n))%(m):_mulmod(n,n,m)

#endif

#ifndef addmod
  #define addmod(n,a,m) ((((m)-(n)) > (a))  ?  ((n)+(a))  :  ((n)+(a)-(m)))
#endif

/* We need to make sure a and b get evaluated into UVs, then do the
 * subtract into a UV before the addmod. */
static INLINE UV submod(UV a, UV b, UV m) {
  UV t1 = m - b;
  return addmod(a, t1, m);
}

/* n^2 + a mod m */
#define sqraddmod(n, a, m)     addmod(sqrmod(n,m),  a,m)
/* i*j + a mod m */
#define muladdmod(i, j, a, m)  addmod(mulmod(i,j,m), a, m)
/* i*j - a mod m */
#define mulsubmod(i, j, a, m)  submod(mulmod(i,j,m), a, m)

/* n^power mod m */
#ifndef HALF_WORD
  static INLINE UV powmod(UV n, UV power, UV m) {
    UV t = 1;
    n %= m;
    while (power) {
      if (power & 1) t = mulmod(t, n, m);
      power >>= 1;
      if (power)     n = sqrmod(n, m);
    }
    return t;
  }
#else
  static INLINE UV powmod(UV n, UV power, UV m) {
    UV t = 1;
    n %= m;
    if (m < HALF_WORD) {
      while (power) {
        if (power & 1) t = (t*n)%m;
        power >>= 1;
        if (power)     n = (n*n)%m;
      }
    } else {
      while (power) {
        if (power & 1) t = mulmod(t, n, m);
        power >>= 1;
        if (power)     n = sqrmod(n,m);
      }
    }
    return t;
  }
#endif

/* n^power + a mod m */
#define powaddmod(n, p, a, m)  addmod(powmod(n,p,m),a,m)

#endif
