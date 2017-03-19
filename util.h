#ifndef MPU_UTIL_H
#define MPU_UTIL_H

#include "ptypes.h"

extern int  _XS_get_verbose(void);
extern void _XS_set_verbose(int v);
extern int  _XS_get_callgmp(void);
extern void _XS_set_callgmp(int v);

extern int _XS_is_prime(UV x);
extern UV  next_prime(UV x);
extern UV  prev_prime(UV x);

extern UV  _XS_prime_count(UV low, UV high);
extern UV  nth_prime(UV x);
extern UV  nth_prime_upper(UV x);
extern UV  nth_prime_lower(UV x);
extern UV  nth_prime_approx(UV x);
extern UV  prime_count_upper(UV x);
extern UV  prime_count_lower(UV x);
extern UV  prime_count_approx(UV x);
extern UV  twin_prime_count(UV low, UV high);
extern UV  twin_prime_count_approx(UV n);
extern UV  nth_twin_prime(UV n);
extern UV  nth_twin_prime_approx(UV n);
extern UV* n_ramanujan_primes(UV n);
extern UV* n_range_ramanujan_primes(UV nlo, UV nhi);
extern UV* ramanujan_primes(UV* first, UV* last, UV low, UV high);
extern int is_ramanujan_prime(UV n);
extern UV  ramanujan_prime_count(UV lo, UV hi);
extern UV  ramanujan_prime_count_upper(UV n);
extern UV  ramanujan_prime_count_lower(UV n);
extern UV  ramanujan_prime_count_approx(UV n);
extern UV  nth_ramanujan_prime(UV n);
extern UV  nth_ramanujan_prime_upper(UV n);
extern UV  nth_ramanujan_prime_lower(UV n);
extern int sum_primes(UV low, UV high, UV *sum);
extern int sum_primes128(UV n, UV *hisum, UV *losum);
extern void print_primes(UV low, UV high, int fd);

extern int powerof(UV n);
extern int is_power(UV n, UV a);
extern UV rootof(UV n, UV k);
extern int primepower(UV n, UV* prime);
extern UV valuation(UV n, UV k);
extern UV logint(UV n, UV b);
extern UV mpu_popcount(UV n);
extern UV mpu_popcount_string(const char* ptr, int len);

extern signed char* _moebius_range(UV low, UV high);
extern UV*    _totient_range(UV low, UV high);
extern IV     mertens(UV n);
extern long double chebyshev_function(UV n, int which); /* 0 = theta, 1 = psi */

extern long double _XS_ExponentialIntegral(long double x);
extern long double _XS_LogarithmicIntegral(long double x);
extern long double ld_riemann_zeta(long double x);
extern long double _XS_RiemannR(long double x);
extern long double lambertw(long double k);
extern UV inverse_li(UV x);

extern int kronecker_uu(UV a, UV b);
extern int kronecker_su(IV a, UV b);
extern int kronecker_ss(IV a, IV b);

extern UV factorial(UV n);
extern UV binomial(UV n, UV k);
extern IV gcdext(IV a, IV b, IV* u, IV* v, IV* s, IV* t); /* Ext Euclidean */
extern UV modinverse(UV a, UV p);              /* Returns 1/a mod p */
extern UV divmod(UV a, UV b, UV n);            /* Returns a/b mod n */
extern int sqrtmod(UV* s, UV a, UV p);         /* sqrt(a) mod p */
extern int sqrtmod_composite(UV* s, UV a,UV n);/* sqrt(a) mod n */
extern UV chinese(UV* a, UV* n, UV num, int *status); /* Chinese Remainder */

extern UV totient(UV n);
extern int moebius(UV n);
extern UV exp_mangoldt(UV n);
extern UV carmichael_lambda(UV n);
extern UV jordan_totient(UV k, UV n);
extern UV znprimroot(UV n);
extern UV znorder(UV a, UV n);
extern int is_primitive_root(UV a, UV n, int nprime);
#define is_square_free(n)  (moebius(n) != 0)
extern int is_semiprime(UV n);
extern int is_carmichael(UV n);
extern UV  is_quasi_carmichael(UV n);  /* Returns number of bases */
extern UV  pillai_v(UV n);             /* v: v! % n == n-1 && n % v != 1 */

extern UV stirling3(UV n, UV m);
extern IV stirling2(UV n, UV m);
extern IV stirling1(UV n, UV m);

extern IV hclassno(UV n);
extern IV ramanujan_tau(UV n);

extern int strnum_minmax(int min, char* a, STRLEN alen, char* b, STRLEN blen);

extern int from_digit_string(UV* n, const char* s, int base);
extern int from_digit_to_UV(UV* rn, UV* r, int len, int base);
extern int from_digit_to_str(char** rstr, UV* r, int len, int base);
extern int to_digit_array(int* bits, UV n, int base, int length);
extern int to_digit_string(char *s, UV n, int base, int length);
extern int to_string_128(char s[40], IV hi, UV lo);

extern int is_catalan_pseudoprime(UV n);

extern UV gcdz(UV x, UV y);

#if defined(FUNC_isqrt) || defined(FUNC_is_perfect_square)
static UV isqrt(UV n) {
  UV root;
#if BITS_PER_WORD == 32
  if (n >= UVCONST(4294836225)) return UVCONST(65535);
#else
  if (n >= UVCONST(18446744065119617025)) return UVCONST(4294967295);
#endif
  root = (UV) sqrt((double)n);
  while (root*root > n)  root--;
  while ((root+1)*(root+1) <= n)  root++;
  return root;
}
#endif

#if defined(FUNC_icbrt) || defined(FUNC_is_perfect_cube)
static UV icbrt(UV n) {
  UV b, root = 0;
#if BITS_PER_WORD == 32
  int s = 30;
  if (n >= UVCONST(4291015625)) return UVCONST(1625);
#else
  int s = 63;
  if (n >= UVCONST(18446724184312856125)) return UVCONST(2642245);
#endif
  for ( ; s >= 0; s -= 3) {
    root += root;
    b = 3*root*(root+1)+1;
    if ((n >> s) >= b) {
      n -= b << s;
      root++;
    }
  }
  return root;
}
#endif

#if defined(FUNC_ipow)
static UV ipow(UV n, UV k) {
  UV p = 1;
  while (k) {
    if (k & 1) p *= n;
    k >>= 1;
    if (k)     n *= n;
  }
  return p;
}
#endif

#if defined(FUNC_gcd_ui) || defined(FUNC_lcm_ui)
/* If we have a very fast ctz, then use the fast FLINT version of gcd */
#if defined(__GNUC__) && (__GNUC__ >= 4 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4))
#define gcd_ui(x,y) gcdz(x,y)
#else
static UV gcd_ui(UV x, UV y) {
  UV t;
  if (y < x) { t = x; x = y; y = t; }
  while (y > 0) {
    t = y;  y = x % y;  x = t;  /* y1 <- x0 % y0 ; x1 <- y0 */
  }
  return x;
}
#endif
#endif

#ifdef FUNC_lcm_ui
static UV lcm_ui(UV x, UV y) {
  /* Can overflow if lcm(x,y) > 2^64 (e.g. two primes each > 2^32) */
  return x * (y / gcd_ui(x,y));
}
#endif

#ifdef FUNC_is_perfect_square
/* See:  http://mersenneforum.org/showpost.php?p=110896 */
static int is_perfect_square(UV n)
{
  UV m;
  m = n & 127;
  if ((m*0x8bc40d7d) & (m*0xa1e2f5d1) & 0x14020a)  return 0;
  /* This cuts out another 80%: */
  m = n % 63; if ((m*0x3d491df7) & (m*0xc824a9f9) & 0x10f14008) return 0;
  /* m = n % 25; if ((m*0x1929fc1b) & (m*0x4c9ea3b2) & 0x51001005) return 0; */
  m = isqrt(n);
  return m*m == n;
}
#endif

#ifdef FUNC_is_perfect_cube
static int is_perfect_cube(UV n)
{
  UV m;
  if ((n & 3) == 2) return 0;
  /* m = n & 511; if ((m*5016427) & (m*95638165) & 438)  return 0; */
  m = n % 117; if ((m*833230740) & (m*120676722) & 813764715)  return 0;
  m = n % 133; if ((m*76846229) & (m*305817297) & 306336544)  return 0;
  m = icbrt(n);
  return m*m*m == n;
}
#endif

#ifdef FUNC_is_perfect_fifth
static int is_perfect_fifth(UV n)
{
  UV m;
  if ((n & 3) == 2) return 0;
  m = n %  88; if ((m*85413603) & (m*76260301) & 26476550)  return 0;
  m = n %  31; if ((m*80682551) & (m*73523539) & 45414528)  return 0;
  m = n %  41; if ((m*92806493) & (m*130690042) & 35668129)  return 0;
  /* m = n %  25; if ((m*109794298) & (m*105535723) & 16097553)  return 0; */
  m = rootof(n, 5);
  return m*m*m*m*m == n;
}
#endif

#ifdef FUNC_is_perfect_seventh
static int is_perfect_seventh(UV n)
{
  UV m;
  /* if ((n & 3) == 2) return 0; */
  m = n & 511; if ((m*97259473) & (m*51311663) & 894)  return 0;
  m = n %  49; if ((m*109645301) & (m*76482737) & 593520192)  return 0;
  m = n %  71; if ((m*71818386) & (m*38821587) & 35299393)  return 0;
  /* m = n %  43; if ((m*101368253) & (m*814158665) & 142131408)  return 0; */
  /* m = n %  29; if ((m*81935611) & (m*84736134) & 37831965)  return 0; */
  /* m = n % 116; if ((m*348163737) & (m*1539055705) & 2735997248)  return 0; */
  m = rootof(n, 7);
  return m*m*m*m*m*m*m == n;
}
#endif

#if defined(FUNC_clz) || defined(FUNC_ctz) || defined(FUNC_log2floor)
/* log2floor(n) gives the location of the first set bit (starting from left)
 * ctz(n)       gives the number of times n is divisible by 2
 * clz(n)       gives the number of zeros on the left                       */
#if defined(__GNUC__) && (__GNUC__ >= 4 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4))
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
  * The performance of these functions are not super critical.
  * What is:  popcnt, mulmod, and muladd.
  */
#elif defined (_MSC_VER) && _MSC_VER >= 1400 && !defined(__clang__) && !defined(_WIN32_WCE)
 #include <intrin.h>
 #ifdef FUNC_ctz
  static int ctz(UV n) {
    UV tz = 0;
   #if BITS_PER_WORD == 64
    if (_BitScanForward64(&tz, n)) return tz; else return 64;
   #else
    if (_BitScanForward(&tz, n))   return tz; else return 32;
   #endif
  }
 #endif
 #if defined(FUNC_clz) || defined(FUNC_log2floor)
  static int log2floor(UV n) {
    UV lz = 0;
   #if BITS_PER_WORD == 64
    if (_BitScanReverse64(&lz, n)) return lz; else return 0;
   #else
    if (_BitScanReverse(&lz, n))   return lz; else return 0;
   #endif
  }
 #endif
#elif BITS_PER_WORD == 64
 static const unsigned char _debruijn64[64] = {
     63, 0,58, 1,59,47,53, 2, 60,39,48,27,54,33,42, 3, 61,51,37,40,49,18,28,20,
     55,30,34,11,43,14,22, 4, 62,57,46,52,38,26,32,41, 50,36,17,19,29,10,13,21,
     56,45,25,31,35,16, 9,12, 44,24,15, 8,23, 7, 6, 5 };
 #ifdef FUNC_ctz
   static unsigned int ctz(UV n) {
     return n ? _debruijn64[((n & -n)*UVCONST(0x07EDD5E59A4E28C2)) >> 58] : 64;
   }
 #endif
 #if defined(FUNC_clz) || defined(FUNC_log2floor)
   static unsigned int log2floor(UV n) {
     if (n == 0) return 0;
     n |= n >> 1;   n |= n >> 2;   n |= n >> 4;
     n |= n >> 8;   n |= n >> 16;  n |= n >> 32;
     return _debruijn64[((n-(n>>1))*UVCONST(0x07EDD5E59A4E28C2)) >> 58];
   }
 #endif
#else
 #ifdef FUNC_ctz
   static const unsigned char _trail_debruijn32[32] = {
      0, 1,28, 2,29,14,24, 3,30,22,20,15,25,17, 4, 8,
     31,27,13,23,21,19,16, 7,26,12,18, 6,11, 5,10, 9 };
   static unsigned int ctz(UV n) {
     return n ? _trail_debruijn32[((n & -n) * UVCONST(0x077CB531)) >> 27] : 32;
   }
 #endif
 #if defined(FUNC_clz) || defined(FUNC_log2floor)
   static const unsigned char _lead_debruijn32[32] = {
      0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
      8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31 };
   static unsigned int log2floor(UV n) {
     if (n == 0) return 0;
     n |= n >> 1;   n |= n >> 2;   n |= n >> 4;   n |= n >> 8;   n |= n >> 16;
     return _lead_debruijn32[(n * UVCONST(0x07C4ACDD)) >> 27];
   }
 #endif
#endif
#if defined(FUNC_clz) && !defined(clz)
 #define clz(n)  ( (n) ? BITS_PER_WORD-1-log2floor(n) : BITS_PER_WORD )
#endif
#endif  /* End of log2floor, clz, and ctz */

#endif
