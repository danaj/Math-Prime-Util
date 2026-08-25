#ifndef MPU_UTIL_BITS_H
#define MPU_UTIL_BITS_H

#include "ptypes.h"

/* I think uint32_t is a better return type, but we follow GCC's prototype.
 *
 * log2floor(n) gives the location of the first set bit (starting from left)
 * ctz(n)       gives the number of times n is divisible by 2
 * clz(n)       gives the number of zeros on the left
 */
#if defined(__GNUC__) && 100*__GNUC__ + __GNUC_MINOR >= 304
 #if BITS_PER_WORD == 64
  #define ctz(n)        ((n) ?    __builtin_ctzll(n) : 64)
  #define clz(n)        ((n) ?    __builtin_clzll(n) : 64)
  #define log2floor(n)  ((n) ? 63-__builtin_clzll(n) : 0)
 #else
  #define ctz(n)        ((n) ?    __builtin_ctzl(n) : 32)
  #define clz(n)        ((n) ?    __builtin_clzl(n) : 32)
  #define log2floor(n)  ((n) ? 31-__builtin_clzl(n) : 0)
 #endif

 /* For MSC, we need to use _BitScanForward and _BitScanReverse.  The way to
  * get to them has changed, so we're going to only use them on new systems.
  * The 64-bit intrinsics are unavailable on x86 even with a 64-bit UV.
  * The performance of these functions are not super critical.
  * What is:  popcnt, mulmod, and muladd.
  */
#elif defined(_MSC_VER) && _MSC_VER >= 1400 && !defined(__clang__) && \
      !defined(_WIN32_WCE) && (BITS_PER_WORD == 32 || defined(_M_X64))
 #include <intrin.h>
 static MAYBE_UNUSED int ctz(UV n) {
   unsigned long tz = 0;
   #if BITS_PER_WORD == 64
    if (_BitScanForward64(&tz, n)) return (int)tz; else return 64;
   #else
    if (_BitScanForward(&tz, n))   return (int)tz; else return 32;
   #endif
 }
 static MAYBE_UNUSED int log2floor(UV n) {
   unsigned long lz = 0;
   #if BITS_PER_WORD == 64
    if (_BitScanReverse64(&lz, n)) return (int)lz; else return 0;
   #else
    if (_BitScanReverse(&lz, n))   return (int)lz; else return 0;
   #endif
 }
#elif BITS_PER_WORD == 64
 static const unsigned char _debruijn64[64] = {
     63, 0,58, 1,59,47,53, 2, 60,39,48,27,54,33,42, 3, 61,51,37,40,49,18,28,20,
     55,30,34,11,43,14,22, 4, 62,57,46,52,38,26,32,41, 50,36,17,19,29,10,13,21,
     56,45,25,31,35,16, 9,12, 44,24,15, 8,23, 7, 6, 5 };
 static MAYBE_UNUSED int ctz(UV n) {
   return n ? _debruijn64[((n & -n)*UVCONST(0x07EDD5E59A4E28C2)) >> 58] : 64;
 }
 static MAYBE_UNUSED int log2floor(UV n) {
   if (n == 0) return 0;
   n |= n >> 1;   n |= n >> 2;   n |= n >> 4;
   n |= n >> 8;   n |= n >> 16;  n |= n >> 32;
   return _debruijn64[((n-(n>>1))*UVCONST(0x07EDD5E59A4E28C2)) >> 58];
 }
#else
 static const unsigned char _trail_debruijn32[32] = {
      0, 1,28, 2,29,14,24, 3,30,22,20,15,25,17, 4, 8,
     31,27,13,23,21,19,16, 7,26,12,18, 6,11, 5,10, 9 };
 static MAYBE_UNUSED int ctz(UV n) {
   return n ? _trail_debruijn32[((n & -n) * UVCONST(0x077CB531)) >> 27] : 32;
 }
 static const unsigned char _lead_debruijn32[32] = {
      0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
      8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31 };
 static MAYBE_UNUSED int log2floor(UV n) {
   if (n == 0) return 0;
   n |= n >> 1;   n |= n >> 2;   n |= n >> 4;   n |= n >> 8;   n |= n >> 16;
   return _lead_debruijn32[(n * UVCONST(0x07C4ACDD)) >> 27];
 }
#endif

#if !defined(clz)
 #define clz(n)  ( (n) ? BITS_PER_WORD-1-log2floor(n) : BITS_PER_WORD )
#endif

/* GCC 3.4 - 4.1 has broken 64-bit popcount.
 * GCC 4.2+ can generate awful code when it doesn't have asm (GCC bug 36041).
 * When the asm is present (e.g. compile with -march=native on a platform that
 * has them, like Nahelem+), then it is almost as fast as manually written asm.
 */
#if BITS_PER_WORD == 64
 #if defined(__POPCNT__) && defined(__GNUC__) && 100*__GNUC__ + __GNUC_MINOR >= 402
   #define popcnt(b)  __builtin_popcountll(b)
 #else
   static MAYBE_UNUSED int popcnt(UV b) {
     b -= (b >> 1) & 0x5555555555555555;
     b = (b & 0x3333333333333333) + ((b >> 2) & 0x3333333333333333);
     b = (b + (b >> 4)) & 0x0f0f0f0f0f0f0f0f;
     return (b * 0x0101010101010101) >> 56;
   }
 #endif
#else
 static MAYBE_UNUSED int popcnt(UV b) {
   b -= (b >> 1) & 0x55555555;
   b = (b & 0x33333333) + ((b >> 2) & 0x33333333);
   b = (b + (b >> 4)) & 0x0f0f0f0f;
   return (b * 0x01010101) >> 24;
 }
#endif

#endif
