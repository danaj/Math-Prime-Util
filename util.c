#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <float.h>

/* Use long double to get a little more precision when we're calculating the
 * math functions -- especially those calculated with a series.  Long double
 * is defined in C89 (ISO C).  Note that 'long double' on many platforms is
 * identical to 'double so it may buy us nothing.  But it's worth trying.
 *
 * While the type was in C89, math functions using it are in C99.  Some
 * systems didn't really get it right (e.g. NetBSD which left out some
 * functions for 13 years).
 */
#include <math.h>
#if _MSC_VER || defined(__IBMC__) | defined(__IBMCPP__) || (defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L)
  /* math.h should give us these as functions or macros.
   *
   *  extern long double fabsl(long double);
   *  extern long double floorl(long double);
   *  extern long double ceill(long double);
   *  extern long double sqrtl(long double);
   *  extern long double powl(long double, long double);
   *  extern long double expl(long double);
   *  extern long double logl(long double);
   *  extern long double log1pl(long double);
   */
#else
  #define fabsl(x)    (long double) fabs( (double) (x) )
  #define floorl(x)   (long double) floor( (double) (x) )
  #define ceill(x)    (long double) ceil( (double) (x) )
  #define sqrtl(x)    (long double) sqrt( (double) (x) )
  #define powl(x, y)  (long double) pow( (double) (x), (double) (y) )
  #define expl(x)     (long double) exp( (double) (x) )
  #define logl(x)     (long double) log( (double) (x) )
  #define log1pl(x)   (long double) log1p( (double) (x) )
#endif

#ifdef LDBL_INFINITY
  #undef INFINITY
  #define INFINITY LDBL_INFINITY
#elif !defined(INFINITY)
  #define INFINITY (DBL_MAX + DBL_MAX)
#endif
#ifndef LDBL_EPSILON
  #define LDBL_EPSILON 1e-16
#endif
#ifndef LDBL_MAX
  #define LDBL_MAX DBL_MAX
#endif

#include "ptypes.h"
#define FUNC_isqrt 1
#define FUNC_icbrt 1
#define FUNC_lcm_ui 1
#define FUNC_ctz 1
#define FUNC_log2floor 1
#define FUNC_is_perfect_square
#define FUNC_is_perfect_cube
#define FUNC_is_perfect_fifth
#define FUNC_is_perfect_seventh
#define FUNC_next_prime_in_sieve 1
#define FUNC_prev_prime_in_sieve 1
#define FUNC_ipow 1
#include "util.h"
#include "sieve.h"
#include "primality.h"
#include "cache.h"
#include "lmo.h"
#include "factor.h"
#include "mulmod.h"
#include "constants.h"
#include "montmath.h"
#include "csprng.h"
#include "keyval.h"
#include "inverse_interpolate.h"

#define KAHAN_INIT(s) \
  LNV s ## _y, s ## _t; \
  LNV s ## _c = 0.0; \
  LNV s = 0.0;

#define KAHAN_SUM(s, term) \
  do { \
    s ## _y = (term) - s ## _c; \
    s ## _t = s + s ## _y; \
    s ## _c = (s ## _t - s) - s ## _y; \
    s = s ## _t; \
  } while (0)

int _numcmp(const void *a, const void *b) {
  const UV *x = a, *y = b;
  return (*x > *y) ? 1 : (*x < *y) ? -1 : 0;
}

static int _verbose = 0;
void _XS_set_verbose(int v) { _verbose = v; }
int _XS_get_verbose(void) { return _verbose; }

static int _call_gmp = 0;
void _XS_set_callgmp(int v) { _call_gmp = v; }
int  _XS_get_callgmp(void) { return _call_gmp; }

static int _secure = 0;
void _XS_set_secure(void) { _secure = 1; }
int  _XS_get_secure(void) { return _secure; }

/* We'll use this little static sieve to quickly answer small values of
 *   is_prime, next_prime, prev_prime, prime_count
 * for non-threaded Perl it's basically the same as getting the primary
 * cache.  It guarantees we'll have an answer with no waiting on any version.
 */
static const unsigned char prime_sieve30[] =
  {0x01,0x20,0x10,0x81,0x49,0x24,0xc2,0x06,0x2a,0xb0,0xe1,0x0c,0x15,0x59,0x12,
   0x61,0x19,0xf3,0x2c,0x2c,0xc4,0x22,0xa6,0x5a,0x95,0x98,0x6d,0x42,0x87,0xe1,
   0x59,0xa9,0xa9,0x1c,0x52,0xd2,0x21,0xd5,0xb3,0xaa,0x26,0x5c,0x0f,0x60,0xfc,
   0xab,0x5e,0x07,0xd1,0x02,0xbb,0x16,0x99,0x09,0xec,0xc5,0x47,0xb3,0xd4,0xc5,
   0xba,0xee,0x40,0xab,0x73,0x3e,0x85,0x4c,0x37,0x43,0x73,0xb0,0xde,0xa7,0x8e,
   0x8e,0x64,0x3e,0xe8,0x10,0xab,0x69,0xe5,0xf7,0x1a,0x7c,0x73,0xb9,0x8d,0x04,
   0x51,0x9a,0x6d,0x70,0xa7,0x78,0x2d,0x6d,0x27,0x7e,0x9a,0xd9,0x1c,0x5f,0xee,
   0xc7,0x38,0xd9,0xc3,0x7e,0x14,0x66,0x72,0xae,0x77,0xc1,0xdb,0x0c,0xcc,0xb2,
   0xa5,0x74,0xe3,0x58,0xd5,0x4b,0xa7,0xb3,0xb1,0xd9,0x09,0xe6,0x7d,0x23,0x7c,
   0x3c,0xd3,0x0e,0xc7,0xfd,0x4a,0x32,0x32,0xfd,0x4d,0xb5,0x6b,0xf3,0xa8,0xb3,
   0x85,0xcf,0xbc,0xf4,0x0e,0x34,0xbb,0x93,0xdb,0x07,0xe6,0xfe,0x6a,0x57,0xa3,
   0x8c,0x15,0x72,0xdb,0x69,0xd4,0xaf,0x59,0xdd,0xe1,0x3b,0x2e,0xb7,0xf9,0x2b,
   0xc5,0xd0,0x8b,0x63,0xf8,0x95,0xfa,0x77,0x40,0x97,0xea,0xd1,0x9f,0xaa,0x1c,
   0x48,0xae,0x67,0xf7,0xeb,0x79,0xa5,0x55,0xba,0xb2,0xb6,0x8f,0xd8,0x2d,0x6c,
   0x2a,0x35,0x54,0xfd,0x7c,0x9e,0xfa,0xdb,0x31,0x78,0xdd,0x3d,0x56,0x52,0xe7,
   0x73,0xb2,0x87,0x2e,0x76,0xe9,0x4f,0xa8,0x38,0x9d,0x5d,0x3f,0xcb,0xdb,0xad,
   0x51,0xa5,0xbf,0xcd,0x72,0xde,0xf7,0xbc,0xcb,0x49,0x2d,0x49,0x26,0xe6,0x1e,
   0x9f,0x98,0xe5,0xc6,0x9f,0x2f,0xbb,0x85,0x6b,0x65,0xf6,0x77,0x7c,0x57,0x8b,
   0xaa,0xef,0xd8,0x5e,0xa2,0x97,0xe1,0xdc,0x37,0xcd,0x1f,0xe6,0xfc,0xbb,0x8c,
   0xb7,0x4e,0xc7,0x3c,0x19,0xd5,0xa8,0x9e,0x67,0x4a,0xe3,0xf5,0x97,0x3a,0x7e,
   0x70,0x53,0xfd,0xd6,0xe5,0xb8,0x1c,0x6b,0xee,0xb1,0x9b,0xd1,0xeb,0x34,0xc2,
   0x23,0xeb,0x3a,0xf9,0xef,0x16,0xd6,0x4e,0x7d,0x16,0xcf,0xb8,0x1c,0xcb,0xe6,
   0x3c,0xda,0xf5,0xcf};
#define NPRIME_SIEVE30 (sizeof(prime_sieve30)/sizeof(prime_sieve30[0]))

static const unsigned short primes_tiny[] =
  {0,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
   101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,
   193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,
   293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,
   409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503};
#define NPRIMES_TINY (sizeof(primes_tiny)/sizeof(primes_tiny[0]))

/* Return of 2 if n is prime, 0 if not.  Do it fast. */
int is_prime(UV n)
{
  if (n <= 10)
    return (n == 2 || n == 3 || n == 5 || n == 7) ? 2 : 0;

  if (n < UVCONST(200000000)) {
    UV d = n/30;
    UV m = n - d*30;
    unsigned char mtab = masktab30[m];  /* Bitmask in mod30 wheel */
    const unsigned char* sieve;

    /* Return 0 if a multiple of 2, 3, or 5 */
    if (mtab == 0)
      return 0;

    /* Check static tiny sieve */
    if (d < NPRIME_SIEVE30)
      return (prime_sieve30[d] & mtab) ? 0 : 2;

    if (!(n%7) || !(n%11) || !(n%13)) return 0;

    /* Check primary cache */
    if (n <= get_prime_cache(0,0)) {
      int isprime = -1;
      if (n <= get_prime_cache(0, &sieve))
        isprime = 2*((sieve[d] & mtab) == 0);
      release_prime_cache(sieve);
      if (isprime >= 0)
        return isprime;
    }
  }
  return is_prob_prime(n);
}


UV next_prime(UV n)
{
  UV m, next;

  if (n < 30*NPRIME_SIEVE30) {
    next = next_prime_in_sieve(prime_sieve30, n, 30*NPRIME_SIEVE30);
    if (next != 0) return next;
  }

  if (n >= MPU_MAX_PRIME) return 0; /* Overflow */

  if (n < get_prime_cache(0,0)) {
    const unsigned char* sieve;
    UV sieve_size = get_prime_cache(0, &sieve);
    next = (n < sieve_size)  ?  next_prime_in_sieve(sieve, n, sieve_size)  :  0;
    release_prime_cache(sieve);
    if (next != 0) return next;
  }

  m = n % 30;
  do { /* Move forward one. */
    n += wheeladvance30[m];
    m = nextwheel30[m];
  } while (!is_prob_prime(n));
  return n;
}


UV prev_prime(UV n)
{
  UV m, prev;

  if (n < 30*NPRIME_SIEVE30)
    return prev_prime_in_sieve(prime_sieve30, n);

  if (n < get_prime_cache(0,0)) {
    const unsigned char* sieve;
    UV sieve_size = get_prime_cache(0, &sieve);
    prev = (n < sieve_size)  ?  prev_prime_in_sieve(sieve, n)  :  0;
    release_prime_cache(sieve);
    if (prev != 0) return prev;
  }

  m = n % 30;
  do { /* Move back one. */
    n -= wheelretreat30[m];
    m = prevwheel30[m];
  } while (!is_prob_prime(n));
  return n;
}

/******************************************************************************/
/*                                 PRINTING                                   */
/******************************************************************************/

static int my_sprint(char* ptr, UV val) {
  int nchars;
  UV t;
  char *s = ptr;
  do {
    t = val / 10; *s++ = (char) ('0' + val - 10 * t);
  } while ((val = t));
  nchars = s - ptr + 1;  *s = '\n';
  while (--s > ptr) { char c = *s; *s = *ptr; *ptr++ = c; }
  return nchars;
}
static char* write_buf(int fd, char* buf, char* bend) {
  int res = (int) write(fd, buf, bend-buf);
  if (res == -1) croak("print_primes write error");
  return buf;
}
void print_primes(UV low, UV high, int fd) {
  char buf[8000+25];
  char* bend = buf;
  if ((low <= 2) && (high >= 2)) bend += my_sprint(bend,2);
  if ((low <= 3) && (high >= 3)) bend += my_sprint(bend,3);
  if ((low <= 5) && (high >= 5)) bend += my_sprint(bend,5);
  if (low < 7) low = 7;

  if (low <= high) {
    unsigned char* segment;
    UV seg_base, seg_low, seg_high;
    void* ctx = start_segment_primes(low, high, &segment);
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
        bend += my_sprint(bend,p);
        if (bend-buf > 8000) { bend = write_buf(fd, buf, bend); }
      END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
  }
  if (bend > buf) { bend = write_buf(fd, buf, bend); }
}

/******************************************************************************/
/*                     TOTIENT, MOEBIUS, MERTENS                              */
/******************************************************************************/

/* Return a char array with lo-hi+1 elements. mu[k-lo] = µ(k) for k = lo .. hi.
 * It is the callers responsibility to call Safefree on the result. */
signed char* range_moebius(UV lo, UV hi)
{
  signed char* mu;
  UV i, sqrtn = isqrt(hi), count = hi-lo+1;

  if (hi < lo) croak("range_mobius error hi %"UVuf" < lo %"UVuf"\n", hi, lo);

  /* Kuznetsov indicates that the Deléglise & Rivat (1996) method can be
   * modified to work on logs, which allows us to operate with no
   * intermediate memory at all.  Same time as the D&R method, less memory. */
  unsigned char logp;
  UV nextlog, nextlogi;

  Newz(0, mu, count, signed char);
  if (sqrtn*sqrtn != hi && sqrtn < (UVCONST(1)<<(BITS_PER_WORD/2))-1) sqrtn++;

  /* For small ranges, do it by hand */
  if (hi < 100 || count <= 10 || (hi > (1UL<<25) && count < icbrt(hi)/4)) {
    for (i = 0; i < count; i++)
      mu[i] = moebius(lo+i);
    return mu;
  }

  logp = 1; nextlog = 3; /* 2+1 */
  START_DO_FOR_EACH_PRIME(2, sqrtn) {
    UV p2 = p*p;
    if (p > nextlog) {
      logp += 2;   /* logp is 1 | ceil(log(p)/log(2)) */
      nextlog = ((nextlog-1)*4)+1;
    }
    for (i = P_GT_LO(p, p, lo); i >= lo && i <= hi; i += p)
      mu[i-lo] += logp;
    for (i = P_GT_LO(p2, p2, lo); i >= lo && i <= hi; i += p2)
      mu[i-lo] = 0x80;
  } END_DO_FOR_EACH_PRIME

  logp = log2floor(lo);
  nextlogi = (UVCONST(2) << logp) - lo;
  for (i = 0; i < count; i++) {
    unsigned char a = mu[i];
    if (i >= nextlogi) nextlogi = (UVCONST(2) << ++logp) - lo;
    if (a & 0x80)       { a = 0; }
    else if (a >= logp) { a =  1 - 2*(a&1); }
    else                { a = -1 + 2*(a&1); }
    mu[i] = a;
  }
  if (lo == 0)  mu[0] = 0;

  return mu;
}

static short* mertens_array(UV hi)
{
  signed char* mu;
  short* M;
  UV i;

  /* We could blend this with range_moebius but it seems not worth it. */
  mu = range_moebius(0, hi);
  New(0, M, hi+1, short);
  M[0] = 0;
  for (i = 1; i <= hi; i++)
    M[i] = M[i-1] + mu[i];
  Safefree(mu);

  return M;
}

UV* range_totient(UV lo, UV hi) {
  UV* totients;
  UV i, seg_base, seg_low, seg_high, count = hi-lo+1;
  unsigned char* segment;
  void* ctx;

  if (hi < lo) croak("range_totient error hi %"UVuf" < lo %"UVuf"\n", hi, lo);
  New(0, totients, count, UV);

  /* Do via factoring if very small or if we have a small range */
  if (hi < 100 || count <= 10 || hi/count > 1000) {
    for (i = 0; i < count; i++)
      totients[i] = totient(lo+i);
    return totients;
  }

  if (hi == UV_MAX) {
    totients[--count] = totient(UV_MAX);
    hi--;
  }

  /* If doing a full sieve, do it monolithic.  Faster. */
  if (lo == 0) {
    UV* prime;
    double loghi = log(hi);
    UV max_index = (hi < 67)     ? 18
                 : (hi < 355991) ? 15+(hi/(loghi-1.09))
                 : (hi/loghi) * (1.0+1.0/loghi+2.51/(loghi*loghi));
    UV j, index, nprimes = 0;

    New(0, prime, max_index, UV);  /* could use prime_count_upper(hi) */
    memset(totients, 0, count * sizeof(UV));
    for (i = 2; i <= hi/2; i++) {
      index = 2*i;
      if ( !(i&1) ) {
        if (i == 2) { totients[2] = 1; prime[nprimes++] = 2; }
        totients[index] = totients[i]*2;
      } else {
        if (totients[i] == 0) {
          totients[i] = i-1;
          prime[nprimes++] = i;
        }
        for (j=0; j < nprimes && index <= hi; index = i*prime[++j]) {
          if (i % prime[j] == 0) {
            totients[index] = totients[i]*prime[j];
            break;
          } else {
            totients[index] = totients[i]*(prime[j]-1);
          }
        }
      }
    }
    Safefree(prime);
    /* All totient values have been filled in except the primes.  Mark them. */
    for (i = ((hi/2) + 1) | 1; i <= hi; i += 2)
      if (totients[i] == 0)
        totients[i] = i-1;
    totients[1] = 1;
    totients[0] = 0;
    return totients;
  }

  for (i = 0; i < count; i++) {
    UV v = lo+i, nv = v;
    if (v % 2 == 0)  nv -= nv/2;
    if (v % 3 == 0)  nv -= nv/3;
    if (v % 5 == 0)  nv -= nv/5;
    totients[i] = nv;
  }

  ctx = start_segment_primes(7, hi/2, &segment);
  while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
    START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high ) {
      for (i = P_GT_LO(2*p,p,lo); i >= lo && i <= hi; i += p)
        totients[i-lo] -= totients[i-lo]/p;
    } END_DO_FOR_EACH_SIEVE_PRIME
  }
  end_segment_primes(ctx);

  /* Fill in all primes */
  for (i = (lo | 1) - lo; i < count; i += 2)
    if (totients[i] == i+lo)
      totients[i]--;
  if (lo <= 1) totients[1-lo] = 1;

  return totients;
}

#if 0
IV mertens(UV n) {
  /* See Deléglise and Rivat (1996) for O(n^2/3 log(log(n))^1/3) algorithm.
   * This implementation uses their lemma 2.1 directly, so is ~ O(n).
   * In serial it is quite a bit faster than segmented summation of mu
   * ranges, though the latter seems to be a favored method for GPUs.
   */
  UV u, j, m, nmk, maxmu;
  signed char* mu;
  short* M;   /* 16 bits is enough range for all 32-bit M => 64-bit n */
  IV sum;

  if (n <= 1)  return n;
  u = isqrt(n);
  maxmu = (n/(u+1));              /* maxmu lets us handle u < sqrt(n) */
  if (maxmu < u) maxmu = u;
  mu = range_moebius(0, maxmu);
  New(0, M, maxmu+1, short);      /* Works up to maxmu < 7613644886 */
  M[0] = 0;
  for (j = 1; j <= maxmu; j++)
    M[j] = M[j-1] + mu[j];
  sum = M[u];
  for (m = 1; m <= u; m++) {
    if (mu[m] != 0) {
      IV inner_sum = 0;
      UV lower = (u/m) + 1;
      UV last_nmk = n/(m*lower);
      UV this_k = 0;
      UV next_k = n/(m*1);
      UV nmkm = m * 2;
      for (nmk = 1; nmk <= last_nmk; nmk++, nmkm += m) {
        this_k = next_k;
        next_k = n/nmkm;
        inner_sum += M[nmk] * (this_k - next_k);
      }
      sum += (mu[m] > 0) ? -inner_sum : inner_sum;
    }
  }
  Safefree(M);
  Safefree(mu);
  return sum;
}
#endif

typedef struct {
  UV n;
  IV sum;
} mertens_value_t;
static void _insert_mert_hash(mertens_value_t *H, UV hsize, UV n, IV sum) {
  UV idx = n % hsize;
  H[idx].n = n;
  H[idx].sum = sum;
}
static int _get_mert_hash(mertens_value_t *H, UV hsize, UV n, IV *sum) {
  UV idx = n % hsize;
  if (H[idx].n == n) {
    *sum = H[idx].sum;
    return 1;
  }
  return 0;
}

/* Thanks to Trizen for this algorithm. */
static IV _rmertens(UV n, UV maxmu, short *M, mertens_value_t *H, UV hsize) {
  UV s, k, ns, nk, nk1, mk, mnk;
  IV sum;

  if (n <= maxmu)
    return M[n];

  if (_get_mert_hash(H, hsize, n, &sum))
    return sum;

  s = isqrt(n);
  ns = n / (s+1);
  sum = 1;

#if 0
  for (k = 2; k <= ns; k++)
    sum -= _rmertens(n/k, maxmu, M, H, hsize);
  for (k = 1; k <= s; k++)
    sum -= M[k] * (n/k - n/(k+1));
#else
  /* Take the above: merge the loops and iterate the divides. */
  if (s != ns && s != ns+1) croak("mertens  s / ns");
  nk  = n;
  nk1 = n/2;
  sum -= (nk - nk1);
  for (k = 2; k <= ns; k++) {
    nk = nk1;
    nk1 = n/(k+1);
    mnk = (nk <= maxmu)  ?  M[nk]  :  _rmertens(nk, maxmu, M, H, hsize);
    mk  = (k  <= maxmu)  ?  M[k]   :  _rmertens(k,  maxmu, M, H, hsize);
    sum -= mnk + mk * (nk-nk1);
  }
  if (s > ns)
    sum -= _rmertens(s, maxmu, M, H, hsize) * (n/s - n/(s+1));
#endif

  _insert_mert_hash(H, hsize, n, sum);
  return sum;
}

static short* _prep_rmertens(UV n, UV* pmaxmu, UV* phsize) {
  UV j = icbrt(n);
  UV maxmu = 1 * j * j;
  UV hsize = next_prime(100 + 8*j);

  /* At large sizes, start clamping memory use. */
  if (maxmu > 100000000UL) {
    /* Exponential decay, reduce by factor of 1 to 8 */
    float rfactor = 1.0 + 7.0 * (1.0 - exp(-(float)maxmu/15000000000.0));
    maxmu /= rfactor;
    hsize = next_prime(hsize * 16);  /* Increase the result cache size */
  }

#if BITS_PER_WORD == 64
  /* A 16-bit signed short will overflow at maxmu > 7613644883 */
  if (maxmu > UVCONST(7613644883))  maxmu = UVCONST(7613644883);
#endif

  *pmaxmu = maxmu;
  *phsize = hsize;
  return mertens_array(maxmu);
}

IV mertens(UV n) {
  UV j, maxmu, hsize;
  short* M;   /* 16 bits is enough range for all 32-bit M => 64-bit n */
  mertens_value_t *H;  /* Cache of calculated values */
  IV sum;

  if (n <= 512) {
    static signed char MV16[33] = {0,-1,-4,-3,-1,-4,2,-4,-2,-1,0,-4,-5,-3,3,-1,-1,-3,-7,-2,-4,2,1,-1,-2,1,1,-3,-6,-6,-6,-5,-4};
    j = n/16;
    sum = MV16[j];
    for (j = j*16 + 1; j <= n; j++)
      sum += moebius(j);
    return sum;
  }

  M = _prep_rmertens(n, &maxmu, &hsize);
  Newz(0, H, hsize, mertens_value_t);

  sum = _rmertens(n, maxmu, M, H, hsize);

  Safefree(H);
  Safefree(M);
  return sum;
}

static const signed char _small_liouville[16] = {-1,1,-1,-1,1,-1,1,-1,-1,1,1,-1,-1,-1,1,1};

static signed char* liouville_array(UV hi)
{
  signed char* l;
  UV a, b, k;

  if (hi < 16) hi = 15;
  New(0, l, hi+1, signed char);
  memcpy(l, _small_liouville, 16);
  if (hi >= 16) memset(l+16, -1, hi-16+1);

  for (a = 16; a <= hi; a = b+1) {
    /* TODO: 2*a >= UV_MAX */
    b = (2*a-1 <= hi)  ?  2*a-1  :  hi;
    START_DO_FOR_EACH_PRIME(2, isqrt(b)) {
      for (k = 2*p; k <= b; k += p) {
        if (k >= a)
          l[k] = -1 * l[k/p];
      }
    } END_DO_FOR_EACH_PRIME
  }

  return l;
}

int liouville(UV n) {
  if (n < 16)
    return _small_liouville[n];
  else
    return( (prime_bigomega(n) & 1) ? -1 : 1 );
}

IV sumliouville(UV n) {
  short* M;
  mertens_value_t *H;
  UV j, maxmu, hsize, k, nk, sqrtn;
  IV sum;

  if (n <= 96) {
    signed char* l = liouville_array(n);
    for (sum = 0, j = 1; j <= n; j++)
      sum += l[j];
    Safefree(l);
    return sum;
  }

  M = _prep_rmertens(n, &maxmu, &hsize);
  Newz(0, H, hsize, mertens_value_t);

  sqrtn = isqrt(n);
  sum = _rmertens(n, maxmu, M, H, hsize);
  for (k = 2; k <= sqrtn; k++) {
    nk = n / (k*k);
    if (nk == 1) break;
    sum += (nk <= maxmu) ? M[nk] : _rmertens(nk, maxmu, M, H, hsize);
  }
  sum += (sqrtn + 1 - k);  /* all k where n/(k*k) == 1 */
  /* TODO: find method to get exact number of n/(k*k)==1 .. 4.  Halves k */
  /*       Ends up with method like Lehmer's g. */

  Safefree(H);
  Safefree(M);
  return sum;
}

/* This paper shows an algorithm for sieving an interval:
 *https://www.ams.org/journals/mcom/2008-77-263/S0025-5718-08-02036-X/S0025-5718-08-02036-X.pdf */
signed char* range_liouville(UV lo, UV hi)
{
  UV i;
  signed char *l;
  if (hi < lo) croak("range_liouvillle error hi %"UVuf" < lo %"UVuf"\n",hi,lo);
  unsigned char *nf = range_nfactor_sieve(lo, hi, 1);
  New(0, l, hi-lo+1, signed char);
  for (i = 0; i < hi-lo+1; i++)
    l[i] = (nf[i] & 1) ? -1 : 1;
  Safefree(nf);
  return l;
}

/******************************************************************************/
/*                             POWERS and ROOTS                               */
/******************************************************************************/

/* There are at least 4 ways to do this, plus hybrids.
 * 1) use a table.  Great for 32-bit, too big for 64-bit.
 * 2) Use pow() to check.  Relatively slow and FP is always dangerous.
 * 3) factor or trial factor.  Slow for 64-bit.
 * 4) Dietzfelbinger algorithm 2.3.5.  Quite slow.
 * This currently uses a hybrid of 1 and 2.
 */
uint32_t powerof(UV n) {
  UV t;
  if ((n <= 3) || (n == UV_MAX)) return 1;
  if ((n & (n-1)) == 0)          return ctz(n);  /* powers of 2    */
  if (is_perfect_square(n))      return 2 * powerof(isqrt(n));
  if (is_perfect_cube(n))        return 3 * powerof(icbrt(n));

  /* Simple rejection filter for non-powers of 5-37.  Rejects 47.85%. */
  t = n & 511; if ((t*77855451) & (t*4598053) & 862)  return 1;

  if (is_perfect_fifth(n))       return 5 * powerof(rootint(n,5));
  if (is_perfect_seventh(n))     return 7 * powerof(rootint(n,7));

  if (n > 177146 && n <= UVCONST(1977326743)) {
    switch (n) { /* Check for powers of 11, 13, 17, 19 within 32 bits */
      case 177147: case 48828125: case 362797056: case 1977326743: return 11;
      case 1594323: case 1220703125: return 13;
      case 129140163: return 17;
      case 1162261467: return 19;
      default:  break;
    }
  }
#if BITS_PER_WORD == 64
  if (n >= UVCONST(8589934592)) {
    /* The Bloom filters reject about 90% of inputs each, about 99% for two.
     * Bach/Sorenson type sieves do about as well, but are much slower due
     * to using a powmod. */
    if ( (t = n %121, !((t*19706187) & (t*61524433) & 876897796)) &&
         (t = n % 89, !((t*28913398) & (t*69888189) & 2705511937U)) ) {
      /* (t = n % 67, !((t*117621317) & (t*48719734) & 537242019)) ) { */
      UV root = rootint(n,11);
      if (n == ipow(root,11)) return 11;
    }
    if ( (t = n %131, !((t*1545928325) & (t*1355660813) & 2771533888U)) &&
         (t = n % 79, !((t*48902028) & (t*48589927) & 404082779)) ) {
      /* (t = n % 53, !((t*79918293) & (t*236846524) & 694943819)) ) { */
      UV root = rootint(n,13);
      if (n == ipow(root,13)) return 13;
    }
    switch (n) {
      case UVCONST(762939453125):
      case UVCONST(16926659444736):
      case UVCONST(232630513987207):
      case UVCONST(100000000000000000):
      case UVCONST(505447028499293771):
      case UVCONST(2218611106740436992):
      case UVCONST(8650415919381337933):  return 17;
      case UVCONST(19073486328125):
      case UVCONST(609359740010496):
      case UVCONST(11398895185373143):
      case UVCONST(10000000000000000000): return 19;
      case UVCONST(94143178827):
      case UVCONST(11920928955078125):
      case UVCONST(789730223053602816):   return 23;
      case UVCONST(68630377364883):       return 29;
      case UVCONST(617673396283947):      return 31;
      case UVCONST(450283905890997363):   return 37;
      default:  break;
    }
  }
#endif
  return 1;
}
int is_power(UV n, UV a)
{
  int ret;
  if (a > 0) {
    if (a == 1 || n <= 1) return 1;
    if ((a % 2) == 0)
      return !is_perfect_square(n) ? 0 : (a == 2) ? 1 : is_power(isqrt(n),a>>1);
    if ((a % 3) == 0)
      return !is_perfect_cube(n) ? 0 : (a == 3) ? 1 : is_power(icbrt(n),a/3);
    if ((a % 5) == 0)
      return !is_perfect_fifth(n) ? 0 : (a == 5) ? 1 :is_power(rootint(n,5),a/5);
  }
  ret = powerof(n);
  if (a != 0) return !(ret % a);  /* Is the max power divisible by a? */
  return (ret == 1) ? 0 : ret;
}

#if BITS_PER_WORD == 64
#define ROOT_MAX_3 41
static const uint32_t root_max[ROOT_MAX_3] = {0,0,4294967295U,2642245,65535,7131,1625,565,255,138,84,56,40,30,23,19,15,13,11,10,9,8,7,6,6,5,5,5,4,4,4,4,3,3,3,3,3,3,3,3,3};
#else
#define ROOT_MAX_3 21
static const uint32_t root_max[ROOT_MAX_3] = {0,0,65535,1625,255,84,40,23,15,11,9,7,6,5,4,4,3,3,3,3,3};
#endif

UV rootint(UV n, UV k) {
  UV lo, hi, max;
  if (k == 0) return 0;
  if (k == 1) return n;
  if (k == 2) return isqrt(n);
  if (k == 3) return icbrt(n);

  /* Bracket between powers of 2, but never exceed max power so ipow works */
  max = 1 + ((k >= ROOT_MAX_3) ? 2 : root_max[k]);
  lo = UVCONST(1) << (log2floor(n)/k);
  hi = ((lo*2) < max) ? lo*2 : max;

  /* Binary search */
  while (lo < hi) {
    UV mid = lo + (hi-lo)/2;
    if (ipow(mid,k) <= n) lo = mid+1;
    else                  hi = mid;
  }
  return lo-1;
}

/* Like ipow but returns UV_MAX if overflow */
UV ipowsafe(UV n, UV k) {
  UV p = 1;

  if (k < ROOT_MAX_3)
    return (n <= root_max[k]) ? ipow(n,k) : UV_MAX;

  while (k) {
    if (k & 1) { if (UV_MAX/n < p) return UV_MAX;  p *= n; }
    k >>= 1;
    if (k)     { if (UV_MAX/n < n) return UV_MAX;  n *= n; }
  }
  return p;
}

int primepower(UV n, UV* prime)
{
  int power = 0;
  if (n < 2) return 0;
  /* Check for small divisors */
  if (!(n&1)) {
    if (n & (n-1)) return 0;
    *prime = 2;
    return ctz(n);
  }
  if ((n%3) == 0) {
    /* if (UVCONST(12157665459056928801) % n) return 0; */
    do { n /= 3; power++; } while (n > 1 && (n%3) == 0);
    if (n != 1) return 0;
    *prime = 3;
    return power;
  }
  if ((n%5) == 0) {
    do { n /= 5; power++; } while (n > 1 && (n%5) == 0);
    if (n != 1) return 0;
    *prime = 5;
    return power;
  }
  if ((n%7) == 0) {
    do { n /= 7; power++; } while (n > 1 && (n%7) == 0);
    if (n != 1) return 0;
    *prime = 7;
    return power;
  }
  if (is_prob_prime(n))
    { *prime = n; return 1; }
  /* Composite.  Test for perfect power with prime root. */
  power = powerof(n);
  if (power == 1) power = 0;
  if (power) {
    UV root = rootint(n, (UV)power);
    if (is_prob_prime(root))
      *prime = root;
    else
      power = 0;
  }
  return power;
}

UV valuation(UV n, UV k)
{
  UV v = 0;
  UV kpower = k;
  if (k < 2 || n < 2) return 0;
  if (k == 2) return ctz(n);
  while ( !(n % kpower) ) {
    kpower *= k;
    v++;
  }
  return v;
}

UV logint(UV n, UV b)
{
  /* UV e;  for (e=0; n; n /= b) e++;  return e-1; */
  UV v, e = 0;
  if (b == 2)
    return log2floor(n);
  if (n > UV_MAX/b) {
    n /= b;
    e = 1;
  }
  for (v = b; v <= n; v *= b)
    e++;
  return e;
}

UV mpu_popcount_string(const char* ptr, uint32_t len)
{
  uint32_t count = 0, i, j, d, v, power, slen, *s, *sptr;

  while (len > 0 && (*ptr == '0' || *ptr == '+' || *ptr == '-'))
    {  ptr++;  len--;  }

  /* Create s as array of base 10^8 numbers */
  slen = (len + 7) / 8;
  Newz(0, s, slen, uint32_t);
  for (i = 0; i < slen; i++) {  /* Chunks of 8 digits */
    for (j = 0, d = 0, power = 1;  j < 8 && len > 0;  j++, power *= 10) {
      v = ptr[--len] - '0';
      if (v > 9) croak("Parameter '%s' must be a positive integer",ptr);
      d += power * v;
    }
    s[slen - 1 - i] = d;
  }
  /* Repeatedly count and divide by 2 across s */
  while (slen > 1) {
    if (s[slen-1] & 1)  count++;
    sptr = s;
    if (s[0] == 1) {
      if (--slen == 0) break;
      *++sptr += 100000000;
    }
    for (i = 0; i < slen; i++) {
      if ( (i+1) < slen  &&  sptr[i] & 1 ) sptr[i+1] += 100000000;
      s[i] = sptr[i] >> 1;
    }
  }
  /* For final base 10^8 number just do naive popcnt */
  for (d = s[0]; d > 0; d >>= 1)
    if (d & 1)
      count++;
  Safefree(s);
  return count;
}


/* How many times does 2 divide n? */
#define padic2(n)  ctz(n)
#define IS_MOD8_3OR5(x)  (((x)&7)==3 || ((x)&7)==5)

static int kronecker_uu_sign(UV a, UV b, int s) {
  while (a) {
    int r = padic2(a);
    if (r) {
      if ((r&1)  &&  IS_MOD8_3OR5(b))  s = -s;
      a >>= r;
    }
    if (a & b & 2)  s = -s;
    { UV t = b % a;  b = a;  a = t; }
  }
  return (b == 1) ? s : 0;
}

int kronecker_uu(UV a, UV b) {
  int r, s;
  if (b & 1)   return kronecker_uu_sign(a, b, 1);
  if (!(a&1))  return 0;
  s = 1;
  r = padic2(b);
  if (r) {
    if ((r&1) && IS_MOD8_3OR5(a))  s = -s;
    b >>= r;
  }
  return kronecker_uu_sign(a, b, s);
}

int kronecker_su(IV a, UV b) {
  int r, s;
  if (a >= 0)  return kronecker_uu(a, b);
  if (b == 0)  return (a == 1 || a == -1) ? 1 : 0;
  s = 1;
  r = padic2(b);
  if (r) {
    if (!(a&1))  return 0;
    if ((r&1) && IS_MOD8_3OR5(a))  s = -s;
    b >>= r;
  }
  a %= (IV) b;
  if (a < 0)  a += b;
  return kronecker_uu_sign(a, b, s);
}

int kronecker_ss(IV a, IV b) {
  if (a >= 0 && b >= 0)
    return (b & 1)  ?  kronecker_uu_sign(a, b, 1)  :  kronecker_uu(a,b);
  if (b >= 0)
    return kronecker_su(a, b);
  return kronecker_su(a, -b) * ((a < 0) ? -1 : 1);
}

UV primorial(UV n) {
  UV pi, prim = 1, max = (sizeof(UV) <= 4) ? 28 : 52;
  if (n > max) return 0;
  for (pi = 1; pi <= max; pi++) {
    if (n < primes_tiny[pi]) break;
    prim *= primes_tiny[pi];
  }
  return prim;
}
UV factorial(UV n) {
  UV i, r = 1;
  if ( (n > 12 && sizeof(UV) <= 4) || (n > 20 && sizeof(UV) <= 8) ) return 0;
  for (i = 2; i <= n; i++)
    r *= i;
  return r;
}

UV binomial(UV n, UV k) {    /* Thanks to MJD and RosettaCode for ideas */
  UV d, g, r = 1;
  if (k == 0) return 1;
  if (k == 1) return n;
  if (k >= n) return (k == n);
  if (k > n/2) k = n-k;
  for (d = 1; d <= k; d++) {
    if (r >= UV_MAX/n) {  /* Possible overflow */
      UV nr, dr;  /* reduced numerator / denominator */
      g = gcd_ui(n, d);  nr = n/g;  dr = d/g;
      g = gcd_ui(r, dr);  r = r/g;  dr = dr/g;
      if (r >= UV_MAX/nr) return 0;  /* Unavoidable overflow */
      r *= nr;
      r /= dr;
      n--;
    } else {
      r *= n--;
      r /= d;
    }
  }
  return r;
}

UV stirling3(UV n, UV m) {   /* Lah numbers */
  UV f1, f2;

  if (m == n) return 1;
  if (n == 0 || m == 0 || m > n) return 0;
  if (m == 1) return factorial(n);

  f1 = binomial(n, m);
  if (f1 == 0)  return 0;
  f2 = binomial(n-1, m-1);
  if (f2 == 0 || f1 >= UV_MAX/f2)  return 0;
  f1 *= f2;
  f2 = factorial(n-m);
  if (f2 == 0 || f1 >= UV_MAX/f2)  return 0;
  return f1 * f2;
}

IV stirling2(UV n, UV m) {
  UV f;
  IV j, k, t, s = 0;

  if (m == n) return 1;
  if (n == 0 || m == 0 || m > n) return 0;
  if (m == 1) return 1;

  if ((f = factorial(m)) == 0) return 0;
  for (j = 1; j <= (IV)m; j++) {
    t = binomial(m, j);
    for (k = 1; k <= (IV)n; k++) {
      if (t == 0 || j >= IV_MAX/t) return 0;
      t *= j;
    }
    if ((m-j) & 1) t *= -1;
    s += t;
  }
  return s/f;
}

IV stirling1(UV n, UV m) {
  IV k, t, b1, b2, s2, s = 0;

  if (m == n) return 1;
  if (n == 0 || m == 0 || m > n) return 0;
  if (m == 1) {
    UV f = factorial(n-1);
    if (f>(UV)IV_MAX) return 0;
    return (n&1) ? ((IV)f) : -((IV)f);
  }

  for (k = 1; k <= (IV)(n-m); k++) {
    b1 = binomial(k + n - 1, n - m + k);
    b2 = binomial(2 * n - m, n - m - k);
    s2 = stirling2(n - m + k, k);
    if (b1 == 0 || b2 == 0 || s2 == 0 || b1 > IV_MAX/b2) return 0;
    t = b1 * b2;
    if (s2 > IV_MAX/t) return 0;
    t *= s2;
    s += (k & 1) ? -t : t;
  }
  return s;
}

UV totient(UV n) {
  UV i, nfacs, totient, lastf, facs[MPU_MAX_FACTORS+1];
  if (n <= 1) return n;
  totient = 1;
  /* phi(2m) = 2phi(m) if m even, phi(m) if m odd */
  while ((n & 0x3) == 0) { n >>= 1; totient <<= 1; }
  if ((n & 0x1) == 0) { n >>= 1; }
  /* factor and calculate totient */
  nfacs = factor(n, facs);
  lastf = 0;
  for (i = 0; i < nfacs; i++) {
    UV f = facs[i];
    if (f == lastf) { totient *= f;               }
    else            { totient *= f-1;  lastf = f; }
  }
  return totient;
}

static const UV jordan_overflow[5] =
#if BITS_PER_WORD == 64
  {UVCONST(4294967311), 2642249, 65537, 7133, 1627};
#else
  {UVCONST(     65537),    1627,   257,   85,   41};
#endif
UV jordan_totient(UV k, UV n) {
  UV factors[MPU_MAX_FACTORS+1];
  int nfac, i;
  UV totient;
  if (k == 0 || n <= 1) return (n == 1);
  if (k > 6 || (k > 1 && n >= jordan_overflow[k-2])) return 0;

  totient = 1;
  /* Similar to Euler totient, shortcut even inputs */
  while ((n & 0x3) == 0) { n >>= 1; totient *= (1<<k); }
  if ((n & 0x1) == 0) { n >>= 1; totient *= ((1<<k)-1); }
  nfac = factor(n,factors);
  for (i = 0; i < nfac; i++) {
    UV p = factors[i];
    UV pk = ipow(p,k);
    totient *= (pk-1);
    while (i+1 < nfac && p == factors[i+1]) {
      i++;
      totient *= pk;
    }
  }
  return totient;
}

UV carmichael_lambda(UV n) {
  UV fac[MPU_MAX_FACTORS+1];
  int i, nfactors;
  UV lambda = 1;

  if (n < 8) return totient(n);
  if ((n & (n-1)) == 0) return n >> 2;

  i = ctz(n);
  if (i > 0) {
    n >>= i;
    lambda <<= (i>2) ? i-2 : i-1;
  }
  nfactors = factor(n, fac);
  for (i = 0; i < nfactors; i++) {
    UV p = fac[i], pk = p-1;
    while (i+1 < nfactors && p == fac[i+1]) {
      i++;
      pk *= p;
    }
    lambda = lcm_ui(lambda, pk);
  }
  return lambda;
}

int is_carmichael(UV n) {
  UV fac[MPU_MAX_FACTORS+1];
  UV exp[MPU_MAX_FACTORS+1];
  int i, nfactors;

  /* Small or even is not a Carmichael number */
  if (n < 561 || !(n&1)) return 0;

  /* Simple pre-test for square free (odds only) */
  if (!(n% 9) || !(n%25) || !(n%49) || !(n%121) || !(n%169))
    return 0;

  /* Check Korselt's criterion for small divisors */
  if (!(n% 5) && ((n-1) %  4 != 0)) return 0;
  if (!(n% 7) && ((n-1) %  6 != 0)) return 0;
  if (!(n%11) && ((n-1) % 10 != 0)) return 0;
  if (!(n%13) && ((n-1) % 12 != 0)) return 0;
  if (!(n%17) && ((n-1) % 16 != 0)) return 0;
  if (!(n%19) && ((n-1) % 18 != 0)) return 0;
  if (!(n%23) && ((n-1) % 22 != 0)) return 0;

  /* Fast check without having to factor */
  if (n > 5000000) {
    if (!(n%29) && ((n-1) % 28 != 0)) return 0;
    if (!(n%31) && ((n-1) % 30 != 0)) return 0;
    if (!(n%37) && ((n-1) % 36 != 0)) return 0;
    if (!(n%41) && ((n-1) % 40 != 0)) return 0;
    if (!(n%43) && ((n-1) % 42 != 0)) return 0;
    if (!is_pseudoprime(n,2)) return 0;
  }

  nfactors = factor_exp(n, fac, exp);
  if (nfactors < 3)
    return 0;
  for (i = 0; i < nfactors; i++) {
    if (exp[i] > 1  ||  ((n-1) % (fac[i]-1)) != 0)
      return 0;
  }
  return 1;
}

static int is_quasi_base(int nfactors, UV *fac, UV p, UV b) {
  int i;
  for (i = 0; i < nfactors; i++) {
    UV d = fac[i] - b;
    if (d == 0 || (p % d) != 0)
      return 0;
  }
  return 1;
}

/* Returns number of bases that pass */
UV is_quasi_carmichael(UV n) {
  UV nbases, fac[MPU_MAX_FACTORS+1], exp[MPU_MAX_FACTORS+1];
  UV spf, lpf, ndivisors, *divs;
  int i, nfactors;

  if (n < 35) return 0;

  /* Simple pre-test for square free */
  if (!(n% 4) || !(n% 9) || !(n%25) || !(n%49) || !(n%121) || !(n%169))
    return 0;

  nfactors = factor_exp(n, fac, exp);
  /* Must be composite */
  if (nfactors < 2)
    return 0;
  /* Must be square free */
  for (i = 0; i < nfactors; i++)
    if (exp[i] > 1)
      return 0;

  nbases = 0;
  spf = fac[0];
  lpf = fac[nfactors-1];

  /* Algorithm from Hiroaki Yamanouchi, 2015 */
  if (nfactors == 2) {
    divs = _divisor_list(n / spf - 1, &ndivisors);
    for (i = 0; i < (int)ndivisors; i++) {
      UV d = divs[i];
      UV k = spf - d;
      if (d >= spf) break;
      if (is_quasi_base(nfactors, fac, n-k, k))
        nbases++;
    }
  } else {
    divs = _divisor_list(lpf * (n / lpf - 1), &ndivisors);
    for (i = 0; i < (int)ndivisors; i++) {
      UV d = divs[i];
      UV k = lpf - d;
      if (lpf > d && k >= spf) continue;
      if (k != 0 && is_quasi_base(nfactors, fac, n-k, k))
        nbases++;
    }
  }
  Safefree(divs);
  return nbases;
}

int is_semiprime(UV n) {
  UV sp, p, n3, factors[2];

  if (n < 6) return (n == 4);
  if (!(n&1)) return !!is_prob_prime(n>>1);
  if (!(n%3)) return !!is_prob_prime(n/3);
  if (!(n%5)) return !!is_prob_prime(n/5);
  /* 27% of random inputs left */
  n3 = icbrt(n);
  for (sp = 4; sp < 60; sp++) {
    p = primes_tiny[sp];
    if (p > n3)
      break;
    if ((n % p) == 0)
      return !!is_prob_prime(n/p);
  }
  /* 9.8% of random inputs left */
  if (is_def_prime(n)) return 0;
  if (p > n3) return 1; /* past this, n is a composite and larger than p^3 */
  /* 4-8% of random inputs left */
  if (factor_one(n, factors, 0, 0) != 2) return 0;
  return (is_def_prime(factors[0]) && is_def_prime(factors[1]));
}
int is_almost_prime(UV k, UV n) {
  UV p, sp;

  if (k == 0) return (n == 1);
  if (k == 1) return !!is_prob_prime(n);
  if (k == 2) return is_semiprime(n);

  if ((n >> k) == 0) return 0;  /* The smallest k-almost prime is 2^k */

  while (k > 0 && !(n& 1)) { k--; n >>= 1; }
  while (k > 0 && !(n% 3)) { k--; n /=  3; }
  while (k > 0 && !(n% 5)) { k--; n /=  5; }
  while (k > 0 && !(n% 7)) { k--; n /=  7; }
  p = 11;
  if (k >= 5) {
    for (sp = 5; k > 1 && n > 1 && sp < NPRIMES_TINY-1; sp++) {
      p = primes_tiny[sp];
      if (n < ipowsafe(p,k))
        return 0;
      while ((n % p) == 0 && k > 0)
        { k--; n /= p; }
    }
    p = primes_tiny[sp];
  }
  if (k == 0) return (n == 1);
  if (k == 1) return !!is_prob_prime(n);
  if (k == 2) return is_semiprime(n);
  if (n < ipowsafe(p,k)) return 0;

  return ((UV)prime_bigomega(n) == k);
}

int is_fundamental(UV n, int neg) {
  UV r = n & 15;
  if (r) {
    if (!neg) {
      switch (r & 3) {
        case 0:  return (r ==  4) ? 0 : is_square_free(n >> 2);
        case 1:  return is_square_free(n);
        default: break;
      }
    } else {
      switch (r & 3) {
        case 0:  return (r == 12) ? 0 : is_square_free(n >> 2);
        case 3:  return is_square_free(n);
        default: break;
      }
    }
  }
  return 0;
}

static int _totpred(UV n, UV maxd) {
  UV i, ndivisors, *divs;
  int res;

  if (n & 1) return 0;
  if ((n & (n-1)) == 0) return 1;
  n >>= 1;
  if (n == 1) return 1;
  if (n < maxd && is_prime(2*n+1)) return 1;

  divs = _divisor_list(n, &ndivisors);
  for (i = 0, res = 0; i < ndivisors && divs[i] < maxd && res == 0; i++) {
    UV r, d = divs[i], p = 2*d+1;
    if (!is_prime(p)) continue;
    r = n/d;
    while (1) {
      if (r == p || _totpred(r, d)) { res = 1; break; }
      if (r % p) break;
      r /= p;
    }
  }
  Safefree(divs);
  return res;
}

int is_totient(UV n) {
  return (n == 0 || (n & 1))  ?  (n==1)  :  _totpred(n,n);
}


UV inverse_totient_count(UV n) {
  set_t set, sumset;
  keyval_t keyval;
  UV res, i, ndivisors, *divs;

  if (n == 1) return 2;
  if (n < 1 || n & 1) return 0;
  if (is_prime(n >> 1)) { /* Coleman Remark 3.3 (Thm 3.1) and Prop 6.2 */
    if (!is_prime(n+1)) return 0;
    if (n >= 10)        return 2;
  }

  divs = _divisor_list(n, &ndivisors);

  init_set(&set, 2*ndivisors);
  keyval.key = 1;  keyval.val = 1;
  set_addsum(&set, keyval);

  for (i = 0; i < ndivisors; i++) {
    UV d = divs[i],  p = d+1;
    if (is_prime(p)) {
      UV j,  np = d,  v = valuation(n, p);
      init_set(&sumset, ndivisors/2);
      for (j = 0; j <= v; j++) {
        UV k, ndiv = n/np;  /* Loop over divisors of n/np */
        if (np == 1) {
          keyval_t kv;  kv.key = 1; kv.val = 1;
          set_addsum(&sumset, kv);
        } else {
          for (k = 0; k < ndivisors && divs[k] <= ndiv; k++) {
            UV val, d2 = divs[k];
            if ((ndiv % d2) != 0) continue;
            val = set_getval(set, d2);
            if (val > 0) {
              keyval_t kv;  kv.key = d2*np; kv.val = val;
              set_addsum(&sumset, kv);
            }
          }
        }
        /* if (j < v && np > UV_MAX/p) croak("overflow np  d %lu", d); */
        np *= p;
      }
      set_merge(&set, sumset);
      free_set(&sumset);
    }
  }
  Safefree(divs);
  res = set_getval(set, n);
  free_set(&set);
  return res;
}

UV* inverse_totient_list(UV *ntotients, UV n) {
  set_list_t setlist, divlist;
  UV i, ndivisors, *divs, *tlist;
  UV *totlist = 0;

  MPUassert(n <= UV_MAX/7.5, "inverse_totient_list n too large");

  if (n == 1) {
    New(0, totlist, 2, UV);
    totlist[0] = 1;  totlist[1] = 2;
    *ntotients = 2;
    return totlist;
  }
  if (n < 1 || n & 1) {
    *ntotients = 0;
    return totlist;
  }
  if (is_prime(n >> 1)) { /* Coleman Remark 3.3 (Thm 3.1) and Prop 6.2 */
    if (!is_prime(n+1)) {
      *ntotients = 0;
      return totlist;
    }
    if (n >= 10) {
      New(0, totlist, 2, UV);
      totlist[0] = n+1;  totlist[1] = 2*n+2;
      *ntotients = 2;
      return totlist;
    }
  }

  divs = _divisor_list(n, &ndivisors);

  init_setlist(&setlist, 2*ndivisors);
  setlist_addval(&setlist, 1, 1);   /* Add 1 => [1] */

  for (i = 0; i < ndivisors; i++) {
    UV d = divs[i],  p = d+1;
    if (is_prime(p)) {
      UV j,  dp = d,  pp = p,  v = valuation(n, p);
      init_setlist(&divlist, ndivisors/2);
      for (j = 0; j <= v; j++) {
        UV k, ndiv = n/dp;  /* Loop over divisors of n/dp */
        if (dp == 1) {
          setlist_addval(&divlist, 1, 2);   /* Add 1 => [2] */
        } else {
          for (k = 0; k < ndivisors && divs[k] <= ndiv; k++) {
            UV nvals, *vals, d2 = divs[k];
            if ((ndiv % d2) != 0) continue;
            vals = setlist_getlist(&nvals, setlist, d2);
            if (vals != 0)
              setlist_addlist(&divlist, d2 * dp, nvals, vals, pp);
          }
        }
        dp *= p;
        pp *= p;
      }
      setlist_merge(&setlist, divlist);
      free_setlist(&divlist);
    }
  }
  Safefree(divs);
  tlist = setlist_getlist(ntotients, setlist, n);
  if (tlist != 0 && *ntotients > 0) {
    New(0, totlist, *ntotients, UV);
    memcpy(totlist, tlist, *ntotients * sizeof(UV));
    qsort(totlist, *ntotients, sizeof(UV), _numcmp);
  }
  free_setlist(&setlist);
  return totlist;
}


UV pillai_v(UV n) {
  UV v, fac;
  if (n == 0) return 0;
  for (v = 8, fac = 5040 % n; v < n-1 && fac != 0; v++) {
    fac = (n < HALF_WORD) ? (fac*v) % n : mulmod(fac,v,n);
    if (fac == n-1 && (n % v) != 1)
      return v;
  }
  return 0;
}


int moebius(UV n) {
  UV factors[MPU_MAX_FACTORS+1];
  int i, nfactors;

  if (n <= 5) return (n == 1) ? 1 : (n % 4) ? -1 : 0;
  if (n >=  49 && (!(n %   4) || !(n %   9) || !(n %  25) || !(n %  49)))
    return 0;
  if (n >= 361 && (!(n % 121) || !(n % 169) || !(n % 289) || !(n % 361)))
    return 0;
  if (n >= 961 && (!(n % 529) || !(n % 841) || !(n % 961)))
    return 0;

  nfactors = factor(n, factors);
  for (i = 1; i < nfactors; i++)
    if (factors[i] == factors[i-1])
      return 0;
  return (nfactors % 2) ? -1 : 1;
}

UV exp_mangoldt(UV n) {
  UV p;
  if (!primepower(n,&p)) return 1;     /* Not a prime power */
  return p;
}

/* least quadratic non-residue mod p (p may be composite) */
/* The returned result will always be 0 or a prime */
UV qnr(UV n) {
  UV a;

  if (n <= 2) return n;

  /* If n is not a prime, this may or may not succeed */
  if (kronecker_uu(2,n) == -1) return 2;

  if (is_prime(n)) {
    for (a = 3; a < n; a += 2)
      if (kronecker_uu(a,n) == -1)
        return a;
  } else {
#if 0 /* Not terrible, but does more work than we need. */
    for (a = 2; a < n; a = next_prime(a))
      if (!sqrtmod_composite(0, a, n))
        return a;
#endif
    UV fac[MPU_MAX_FACTORS+1];
    int i, nfactors;
    if (!(n&1)) { /* Check and remove all multiples of 2 */
      int e = ctz(n);
      n >>= e;
      if (e >= 2 || n == 1) return 2;
    }
    if (!(n % 3) || !(n % 5) || !(n % 11) || !(n % 13) || !(n % 19)) return 2;
    nfactors = factor_exp(n, fac, 0);
    for (a = 2; a < n; a = next_prime(a)) {
      for (i = 0; i < nfactors; i++)
        if (a < fac[i] && kronecker_uu(a,fac[i]) == -1)
          return a;
    }
  }
  return 0;
}


UV znorder(UV a, UV n) {
  UV fac[MPU_MAX_FACTORS+1];
  UV exp[MPU_MAX_FACTORS+1];
  int i, nfactors;
  UV k, phi;

  if (n <= 1) return n;   /* znorder(x,0) = 0, znorder(x,1) = 1          */
  if (a <= 1) return a;   /* znorder(0,x) = 0, znorder(1,x) = 1  (x > 1) */
  if (gcd_ui(a,n) > 1)  return 0;

  /* Cohen 1.4.3 using Carmichael Lambda */
  phi = carmichael_lambda(n);
  nfactors = factor_exp(phi, fac, exp);
  k = phi;
#if USE_MONTMATH
  if (n & 1) {
    const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);
    UV ma = mont_geta(a, n);
    for (i = 0; i < nfactors; i++) {
      UV b, a1, ek, pi = fac[i], ei = exp[i];
      b = ipow(pi,ei);
      k /= b;
      a1 = mont_powmod(ma, k, n);
      for (ek = 0; a1 != mont1 && ek++ <= ei; a1 = mont_powmod(a1, pi, n))
        k *= pi;
      if (ek > ei) return 0;
    }
  } else
#endif
  for (i = 0; i < nfactors; i++) {
    UV b, a1, ek, pi = fac[i], ei = exp[i];
    b = ipow(pi,ei);
    k /= b;
    a1 = powmod(a, k, n);
    for (ek = 0; a1 != 1 && ek++ <= ei; a1 = powmod(a1, pi, n))
      k *= pi;
    if (ek > ei) return 0;
  }
  return k;
}

UV znprimroot(UV n) {
  UV fac[MPU_MAX_FACTORS+1];
  UV phi_div_fac[MPU_MAX_FACTORS+1];
  UV a, phi, on, r;
  int i, nfactors;

  if (n <= 4) return (n == 0) ? 0 : n-1;
  if (n % 4 == 0)  return 0;

  on = (n&1) ? n : (n>>1);
  a = powerof(on);
  r = rootint(on, a);
  if (!is_prob_prime(r)) return 0;        /* c^a or 2c^a */
  phi = (r-1) * (on/r);                   /* p^a or 2p^a */

  nfactors = factor_exp(phi, fac, 0);
  for (i = 0; i < nfactors; i++)
    phi_div_fac[i] = phi / fac[i];

#if USE_MONTMATH
  if (n & 1) {
    const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);
    for (a = 2; a < n; a++) {
      if (a == 4 || a == 8 || a == 9) continue;  /* Skip some perfect powers */
      /* Skip values we know can't be right: (a|n) = 0 (or 1 for odd primes) */
      if (phi == n-1) { if (kronecker_uu(a, n) != -1) continue; }
      else            { if (gcd_ui(a,n) != 1) continue; }
      r = mont_geta(a, n);
      for (i = 0; i < nfactors; i++)
        if (mont_powmod(r, phi_div_fac[i], n) == mont1)
          break;
      if (i == nfactors) return a;
    }
  } else
#endif
  for (a = 2; a < n; a++) {
    if (a == 4 || a == 8 || a == 9) continue;  /* Skip some perfect powers */
    /* Skip values we know can't be right: (a|n) = 0 (or 1 for odd primes) */
    if (phi == n-1) { if (kronecker_uu(a, n) != -1) continue; }
    else            { if (gcd_ui(a,n) != 1) continue; }
    for (i = 0; i < nfactors; i++)
      if (powmod(a, phi_div_fac[i], n) == 1)
        break;
    if (i == nfactors) return a;
  }
  return 0;
}

int is_primitive_root(UV a, UV n, int nprime) {
  UV s, fac[MPU_MAX_FACTORS+1];
  int i, nfacs;

  if (n <= 1) return n;
  if (a >= n) a %= n;
  if (n <= 4) return a == n-1;
  if (n % 4 == 0)  return 0;

  /* Very simple, but not fast:
   *     s = nprime ? n-1 : totient(n);
   *     return s == znorder(a, n);
   */

  if (gcd_ui(a,n) != 1) return 0;
  if (nprime) {
    s = n-1;
  } else {
    UV on = (n&1) ? n : (n>>1);
    UV k = powerof(on);
    UV r = rootint(on, k);
    if (!is_prob_prime(r)) return 0;        /* c^a or 2c^a */
    s = (r-1) * (on/r);                     /* p^a or 2p^a */
  }
  if (s == n-1 && kronecker_uu(a,n) != -1) return 0;

  /* a^x can be a primitive root only if gcd(x,s) = 1 */
  i = is_power(a,0);
  if (i > 1 && gcd_ui(i, s) != 1) return 0;

#if USE_MONTMATH
  if (n & 1) {
    const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);
    a = mont_geta(a, n);
    /* Quick check for small factors before full factor */
    if ((s % 2) == 0 && mont_powmod(a, s/2, n) == mont1) return 0;
    if ((s % 3) == 0 && mont_powmod(a, s/3, n) == mont1) return 0;
    if ((s % 5) == 0 && mont_powmod(a, s/5, n) == mont1) return 0;
    nfacs = factor_exp(s, fac, 0);
    for (i = 0; i < nfacs; i++)
      if (fac[i] > 5 && mont_powmod(a, s/fac[i], n) == mont1)
        return 0;
  } else
#endif
  {
    /* Quick check for small factors before full factor */
    if ((s % 2) == 0 && powmod(a, s/2, n) == 1) return 0;
    if ((s % 3) == 0 && powmod(a, s/3, n) == 1) return 0;
    if ((s % 5) == 0 && powmod(a, s/5, n) == 1) return 0;
    /* Complete factor and check each one not found above. */
    nfacs = factor_exp(s, fac, 0);
    for (i = 0; i < nfacs; i++)
      if (fac[i] > 5 && powmod(a, s/fac[i], n) == 1)
        return 0;
  }
  return 1;
}

IV gcdext(IV a, IV b, IV* u, IV* v, IV* cs, IV* ct) {
  IV s = 0;  IV os = 1;
  IV t = 1;  IV ot = 0;
  IV r = b;  IV or = a;
  if (a == 0 && b == 0) { os = 0; t = 0; }
  while (r != 0) {
    IV quot = or / r;
    { IV tmp = r; r = or - quot * r;  or = tmp; }
    { IV tmp = s; s = os - quot * s;  os = tmp; }
    { IV tmp = t; t = ot - quot * t;  ot = tmp; }
  }
  if (or < 0) /* correct sign */
    { or = -or; os = -os; ot = -ot; }
  if (u  != 0) *u = os;
  if (v  != 0) *v = ot;
  if (cs != 0) *cs = s;
  if (ct != 0) *ct = t;
  return or;
}

/* Calculate 1/a mod n. */
UV modinverse(UV a, UV n) {
  IV t = 0;  UV nt = 1;
  UV r = n;  UV nr = a;
  while (nr != 0) {
    UV quot = r / nr;
    { UV tmp = nt;  nt = t - quot*nt;  t = tmp; }
    { UV tmp = nr;  nr = r - quot*nr;  r = tmp; }
  }
  if (r > 1) return 0;  /* No inverse */
  if (t < 0) t += n;
  return t;
}

UV divmod(UV a, UV b, UV n) {   /* a / b  mod n */
  UV binv = modinverse(b, n);
  if (binv == 0)  return 0;
  return mulmod(a, binv, n);
}

static UV _powfactor(UV p, UV d, UV m) {
  UV e = 0;
  do { d /= p; e += d; } while (d > 0);
  return powmod(p, e, m);
}

UV factorialmod(UV n, UV m) {  /*  n! mod m */
  UV i, d = n, res = 1;

  if (n >= m || m == 1) return 0;

  if (n <= 10) { /* Keep things simple for small n */
    for (i = 2; i <= n && res != 0; i++)
      res = (res * i) % m;
    return res;
  }

  if (n > m/2 && is_prime(m))    /* Check if we can go backwards */
    d = m-n-1;
  if (d < 2)
    return (d == 0) ? m-1 : 1;   /* Wilson's Theorem: n = m-1 and n = m-2 */

  if (d == n && d > 5000000) {   /* Check for composite m that leads to 0 */
    UV fac[MPU_MAX_FACTORS+1], exp[MPU_MAX_FACTORS+1];
    int j, k, nfacs = factor_exp(m, fac, exp);
    for (j = 0; j < nfacs; j++) {
      UV t = fac[j];
      for (k = 1; (UV)k < exp[j]; k++)
        t *= fac[j];
      if (n >= t) return 0;
    }
  }

#if USE_MONTMATH
  if (m & 1 && d < 40000) {
    const uint64_t npi = mont_inverse(m),  mont1 = mont_get1(m);
    uint64_t monti = mont1;
    res = mont1;
    for (i = 2; i <= d && res != 0; i++) {
      monti = addmod(monti,mont1,m);
      res = mont_mulmod(res,monti,m);
    }
    res = mont_recover(res, m);
  } else
#endif
  if (d < 10000) {
    for (i = 2; i <= d && res != 0; i++)
      res = mulmod(res,i,m);
  } else {
#if 0    /* Monolithic prime walk */
    START_DO_FOR_EACH_PRIME(2, d) {
      UV k = (p > (d>>1))  ?  p  :  _powfactor(p, d, m);
      res = mulmod(res, k, m);
      if (res == 0) break;
    } END_DO_FOR_EACH_PRIME;
#else    /* Segmented prime walk */
    unsigned char* segment;
    UV seg_base, seg_low, seg_high;
    void* ctx = start_segment_primes(7, d, &segment);
    for (i = 1; i <= 3; i++)    /* Handle 2,3,5 assume d>10*/
      res = mulmod(res, _powfactor(2*i - (i>1), d, m), m);
    while (res != 0 && next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
        UV k = (p > (d>>1))  ?  p  :  _powfactor(p, d, m);
        res = mulmod(res, k, m);
        if (res == 0) break;
      END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
#endif
  }

  if (d != n && res != 0) {      /* Handle backwards case */
    if (!(d&1)) res = submod(m,res,m);
    res = modinverse(res,m);
  }

  return res;
}

static int verify_sqrtmod(UV s, UV *rs, UV a, UV p) {
  if (p-s < s)  s = p-s;
  if (mulmod(s, s, p) != a) return 0;
  if (rs != 0) *rs = s;
  return 1;
}
#if !USE_MONTMATH
UV _sqrtmod_prime(UV a, UV p) {
  if ((p % 4) == 3) {
    return powmod(a, (p+1)>>2, p);
  }
  if ((p % 8) == 5) { /* Atkin's algorithm.  Faster than Legendre. */
    UV a2, alpha, beta, b;
    a2 = addmod(a,a,p);
    alpha = powmod(a2,(p-5)>>3,p);
    beta  = mulmod(a2,sqrmod(alpha,p),p);
    b     = mulmod(alpha, mulmod(a, (beta ? beta-1 : p-1), p), p);
    return b;
  }
  if ((p % 16) == 9) { /* Müller's algorithm extending Atkin */
    UV a2, alpha, beta, b, d = 1;
    a2 = addmod(a,a,p);
    alpha = powmod(a2, (p-9)>>4, p);
    beta  = mulmod(a2, sqrmod(alpha,p), p);
    if (sqrmod(beta,p) != p-1) {
      do { d += 2; } while (kronecker_uu(d,p) != -1 && d < p);
      alpha = mulmod(alpha, powmod(d,(p-9)>>3,p), p);
      beta  = mulmod(a2, mulmod(sqrmod(d,p),sqrmod(alpha,p),p), p);
    }
    b = mulmod(alpha, mulmod(a, mulmod(d,(beta ? beta-1 : p-1),p),p),p);
    return b;
  }

  /* Verify Euler condition for odd p */
  if ((p & 1) && powmod(a,(p-1)>>1,p) != 1) return 0;

  /* Algorithm 1.5.1 from Cohen.  Tonelli/Shanks. */
  {
    UV x, q, e, t, z, r, m, b;
    q = p-1;
    e = valuation(q, 2);
    q >>= e;
    t = 3;
    while (kronecker_uu(t, p) != -1) {
      t += 2;
      if (t == 201) {           /* exit if p looks like a composite */
        if ((p % 2) == 0 || powmod(2, p-1, p) != 1 || powmod(3, p-1, p) != 1)
          return 0;
      } else if (t >= 20000) {  /* should never happen */
        return 0;
      }
    }
    z = powmod(t, q, p);
    b = powmod(a, q, p);
    r = e;
    q = (q+1) >> 1;
    x = powmod(a, q, p);
    while (b != 1) {
      t = b;
      for (m = 0; m < r && t != 1; m++)
        t = sqrmod(t, p);
      if (m >= r) break;
      t = powmod(z, UVCONST(1) << (r-m-1), p);
      x = mulmod(x, t, p);
      z = mulmod(t, t, p);
      b = mulmod(b, z, p);
      r = m;
    }
    return x;
  }
  return 0;
}
#else
UV _sqrtmod_prime(UV a, UV p) {
  const uint64_t npi = mont_inverse(p),  mont1 = mont_get1(p);
  a = mont_geta(a,p);

  if ((p % 4) == 3) {
    UV b = mont_powmod(a, (p+1)>>2, p);
    return mont_recover(b, p);
  }

  if ((p % 8) == 5) { /* Atkin's algorithm.  Faster than Legendre. */
    UV a2, alpha, beta, b;
    a2 = addmod(a,a,p);
    alpha = mont_powmod(a2,(p-5)>>3,p);
    beta  = mont_mulmod(a2,mont_sqrmod(alpha,p),p);
    beta  = submod(beta, mont1, p);
    b     = mont_mulmod(alpha, mont_mulmod(a, beta, p), p);
    return mont_recover(b, p);
  }
  if ((p % 16) == 9) { /* Müller's algorithm extending Atkin */
    UV a2, alpha, beta, b, d = 1;
    a2 = addmod(a,a,p);
    alpha = mont_powmod(a2, (p-9)>>4, p);
    beta  = mont_mulmod(a2, mont_sqrmod(alpha,p), p);
    if (mont_sqrmod(beta,p) != submod(0,mont1,p)) {
      do { d += 2; } while (kronecker_uu(d,p) != -1 && d < p);
      d = mont_geta(d,p);
      alpha = mont_mulmod(alpha, mont_powmod(d,(p-9)>>3,p), p);
      beta  = mont_mulmod(a2, mont_mulmod(mont_sqrmod(d,p),mont_sqrmod(alpha,p),p), p);
      beta  = mont_mulmod(submod(beta,mont1,p), d, p);
    } else {
      beta  = submod(beta, mont1, p);
    }
    b = mont_mulmod(alpha, mont_mulmod(a, beta, p), p);
    return mont_recover(b, p);
  }

  /* Verify Euler condition for odd p */
  if ((p & 1) && mont_powmod(a,(p-1)>>1,p) != mont1) return 0;

  /* Algorithm 1.5.1 from Cohen.  Tonelli/Shanks. */
  {
    UV x, q, e, t, z, r, m, b;
    q = p-1;
    e = valuation(q, 2);
    q >>= e;
    t = 3;
    while (kronecker_uu(t, p) != -1) {
      t += 2;
      if (t == 201) {           /* exit if p looks like a composite */
        if ((p % 2) == 0 || powmod(2, p-1, p) != 1 || powmod(3, p-1, p) != 1)
          return 0;
      } else if (t >= 20000) {  /* should never happen */
        return 0;
      }
    }
    t = mont_geta(t, p);
    z = mont_powmod(t, q, p);
    b = mont_powmod(a, q, p);
    r = e;
    q = (q+1) >> 1;
    x = mont_powmod(a, q, p);
    while (b != mont1) {
      t = b;
      for (m = 0; m < r && t != mont1; m++)
        t = mont_sqrmod(t, p);
      if (m >= r) break;
      t = mont_powmod(z, UVCONST(1) << (r-m-1), p);
      x = mont_mulmod(x, t, p);
      z = mont_mulmod(t, t, p);
      b = mont_mulmod(b, z, p);
      r = m;
    }
    return mont_recover(x, p);
  }
  return 0;
}
#endif

int sqrtmod(UV *s, UV a, UV p) {
  if (p == 0) return 0;
  if (a >= p) a %= p;
  if (p <= 2 || a <= 1) return verify_sqrtmod(a, s,a,p);

  return verify_sqrtmod(_sqrtmod_prime(a,p), s,a,p);
}

int sqrtmod_composite(UV *s, UV a, UV n) {
  UV fac[MPU_MAX_FACTORS+1];
  UV exp[MPU_MAX_FACTORS+1];
  UV sqr[MPU_MAX_FACTORS+1];
  UV p, j, k, gcdan;
  int i, nfactors;

  if (n == 0) return 0;
  if (a >= n) a %= n;
  if (n <= 2 || a <= 1) return verify_sqrtmod(a, s,a,n);

  /* Simple existence check.  It's still possible no solution exists.*/
  if (kronecker_uu(a, ((n%4) == 2) ? n/2 : n) == -1) return 0;

  /* if 8|n 'a' must = 1 mod 8, else if 4|n 'a' must = 1 mod 4 */
  if ((n % 4) == 0) {
    if ((n % 8) == 0) {
      if ((a % 8) != 1) return 0;
    } else {
      if ((a % 4) != 1) return 0;
    }
  }

  /* More detailed existence check before factoring.  Still possible. */
  gcdan = gcd_ui(a, n);
  if (gcdan == 1) {
    if ((n % 3) == 0 && kronecker_uu(a, 3) != 1) return 0;
    if ((n % 5) == 0 && kronecker_uu(a, 5) != 1) return 0;
    if ((n % 7) == 0 && kronecker_uu(a, 7) != 1) return 0;
  }

  /* Factor n */
  nfactors = factor_exp(n, fac, exp);

  /* If gcd(a,n)==1, this answers conclusively if a solution exists. */
  if (gcdan == 1) {
    for (i = 0; i < nfactors; i++)
      if (fac[i] > 7 && kronecker_uu(a, fac[i]) != 1) return 0;
    if (s == 0)  /* They only care about existence */
      return 1;
  }

  for (i = 0; i < nfactors; i++) {

    /* Powers of 2 */
    if (fac[i] == 2) {
      if (exp[i] == 1) {
        sqr[i] = a & 1;
      } else if (exp[i] == 2) {
        sqr[i] = 1;  /* and 3 */
      } else {
        UV this_roots[256], next_roots[256];
        UV nthis = 0, nnext = 0;
        this_roots[nthis++] = 1;
        this_roots[nthis++] = 3;
        for (j = 2; j < exp[i]; j++) {
          p = UVCONST(1) << (j+1);
          nnext = 0;
          for (k = 0; k < nthis && nnext < 254; k++) {
            UV r = this_roots[k];
            if (sqrmod(r,p) == (a % p))
              next_roots[nnext++] = r;
            if (sqrmod(p-r,p) == (a % p))
              next_roots[nnext++] = p-r;
          }
          if (nnext == 0) return 0;
          /* copy next exponent's found roots to this one */
          nthis = nnext;
          for (k = 0; k < nnext; k++)
            this_roots[k] = next_roots[k];
        }
        sqr[i] = this_roots[0];
      }
      continue;
    }

    /* p is an odd prime */
    p = fac[i];
    if (!sqrtmod(&(sqr[i]), a, p))
      return 0;

    /* Lift solution of x^2 = a mod p  to  x^2 = a mod p^e */
    for (j = 1; j < exp[i]; j++) {
      UV xk2, yk, expected, sol;
      xk2 = addmod(sqr[i],sqr[i],p);
      yk = modinverse(xk2, p);
      expected = mulmod(xk2,yk,p);
      p *= fac[i];
      sol = submod(sqr[i], mulmod(submod(sqrmod(sqr[i],p), a % p, p), yk, p), p);
      if (expected != 1 || sqrmod(sol,p) != (a % p)) {
        /* printf("a %lu failure to lift to %lu^%d\n", a, fac[i], j+1); */
        return 0;
      }
      sqr[i] = sol;
    }
  }

  /* raise fac[i] */
  for (i = 0; i < nfactors; i++)
    fac[i] = ipow(fac[i], exp[i]);

  p = chinese(sqr, fac, nfactors, &i);
  return (i == 1) ? verify_sqrtmod(p, s, a, n) : 0;
}

UV rootmodp(UV n, UV k, UV p) {
  UV i;

  /* TODO k > 3 */

  n = (p <= 1)  ?  0  :  (n % p);
  if (n <= 2) return n;
  if (k == 0) return 1;
  if (k == 1) return n;
  if (!is_prime(p)) croak("rootmodp: modulus must be prime");
  if (k == 2) { if (sqrtmod(&i, n, p)) return i; return 0; }

  if (k == 3) {
    /* https://www.sciencedirect.com/science/article/pii/S0893965902000319 */
    if ( (p % 3) == 2)
      return powmod(n, 2*((p-2)/3)+1, p);  /* (2p-1)/3 */
    if (powmod(n, (p-1)/3, p) != 1)
      return 0;                        /* no root exists */
    /* It would be nice to return the smallest root */
    if ((p % 9) == 4)
      return powmod(n, 2*((p-4)/9)+1, p);  /* (2p+1)/9 */
    if ((p % 9) == 7)
      return powmod(n, (p+2)/9, p);
    { /* Tonelli-Shanks as shown in Padró and Sáez (2002) */
      UV e, q, h, g, s, y, r, b, x, m, B, t;
      for (e = 0, q = p-1; !(q % 3); q /= 3)
        e++;
      /* MPUassert(e >= 1, "Error in rootmod T-S cube root: e = 0"); */
      for (i = 1, s = 1; i < 10000 && s == 1; i++) {
        h = i % p;  /* "random" */
        if (h == 0) break;
        s = powmod(h, (p-1)/3, p);
      }
      if (s == 1) return 0;
      g = powmod(h, q, p);
      y = g;
      r = e;
      x = ((q % 3) == 2) ? powmod(n, (q-2)/3, p) : powmod(n, (2*q-2)/3, p);
      b = mulmod(sqrmod(n,p),powmod(x,3,p),p);
      x = mulmod(n,x,p);
      while (1) {
        if ((b % p) == 1) {   /* return smallest root x,x*s,x*s*s */
          UV x2 = mulmod(x, s, p), x3 = mulmod(x2, s, p);
          return (x <= x2 && x <= x3) ? x : (x2 <= x3) ? x2 : x3;
        }
        for (m = 1, B = b; m < r; m++) {
          B = powmod(B, 3, p);
          if ((B % p) == 1) break;
        }
        if (m == r) return 0;
        if (s == powmod(b, ipow(3,m-1), p)) {
          t = sqrmod(y,p);
          s = sqrmod(s,p);
        } else {
          t = y;
        }
        t = powmod(t, ipow(3, r-m-1), p);
        y = powmod(t, 3, p);
        r = m;
        x = mulmod(x, t, p);
        b = mulmod(b, y, p);
      }
    }
  }

  /* Easy case: always exists */
  if (gcd_ui(k, p-1) == 1)
    return powmod(n, modinverse(k, p-1), p);

  /* General Euler Criterion for odd p */
  if (powmod(n, (p-1)/gcd_ui(k,p-1), p) != 1)
    return 0;

  /* TODO: generic k-th roots */

  /* .... */
  for (i = 2; i < p; i++) {
    if (powmod(i, k, p) == n)
      return i;
  }
  return 0;
}

UV rootmod(UV n, UV k, UV p) {
  UV i;

  /* TODO:
   *  1) prime powers
   *  2) composites (Factor + CRT + Hensel)
   */

  n = (p <= 1)  ?  0  :  (n % p);
  if (n <= 2) return n;
  if (k == 0) return 1;
  if (k == 1) return n;
  if (k == 2) { if (sqrtmod_composite(&i, n, p)) return i; return 0; }

  /* TODO: combine calls to rootmodp */

  /* .... */
  for (i = 2; i < p; i++) {
    if (powmod(i, k, p) == n)
      return i;
  }
  return 0;
}

/* works only for co-prime inputs and also slower than the algorithm below,
 * but handles the case where IV_MAX < lcm <= UV_MAX.
 */
static UV _simple_chinese(UV* a, UV* n, UV num, int* status) {
  UV i, lcm = 1, res = 0;
  *status = 0;
  if (num == 0) return 0;

  for (i = 0; i < num; i++) {
    UV ni = n[i];
    UV gcd = gcd_ui(lcm, ni);
    if (gcd != 1) return 0;           /* not coprime */
    ni /= gcd;
    if (ni > (UV_MAX/lcm)) return 0;  /* lcm overflow */
    lcm *= ni;
  }
  for (i = 0; i < num; i++) {
    UV p, inverse, term;
    p = lcm / n[i];
    inverse = modinverse(p, n[i]);
    if (inverse == 0) return 0;       /* n's coprime so should never happen */
    term = mulmod(p, mulmod(a[i], inverse, lcm), lcm);
    res = addmod(res, term, lcm);
  }
  *status = 1;
  return res;
}

/* status: 1 ok, -1 no inverse, 0 overflow */
UV chinese(UV* a, UV* n, UV num, int* status) {
  static unsigned short sgaps[] = {7983,3548,1577,701,301,132,57,23,10,4,1,0};
  UV gcd, i, j, lcm, sum, gi, gap;
  *status = 1;
  if (num == 0) return 0;

  /* Sort modulii, largest first */
  for (gi = 0, gap = sgaps[gi]; gap >= 1; gap = sgaps[++gi]) {
    for (i = gap; i < num; i++) {
      UV tn = n[i], ta = a[i];
      for (j = i; j >= gap && n[j-gap] < tn; j -= gap)
        {  n[j] = n[j-gap];  a[j] = a[j-gap];  }
      n[j] = tn;  a[j] = ta;
    }
  }

  if (n[0] > IV_MAX) return _simple_chinese(a,n,num,status);
  lcm = n[0]; sum = a[0] % n[0];
  for (i = 1; i < num; i++) {
    IV u, v, t, s;
    UV vs, ut;
    gcd = gcdext(lcm, n[i], &u, &v, &s, &t);
    if (gcd != 1 && ((sum % gcd) != (a[i] % gcd))) { *status = -1; return 0; }
    if (s < 0) s = -s;
    if (t < 0) t = -t;
    if (s > (IV)(IV_MAX/lcm)) return _simple_chinese(a,n,num,status);
    lcm *= s;
    if (u < 0) u += lcm;
    if (v < 0) v += lcm;
    vs = mulmod((UV)v, (UV)s, lcm);
    ut = mulmod((UV)u, (UV)t, lcm);
    sum = addmod(  mulmod(vs, sum, lcm),  mulmod(ut, a[i], lcm),  lcm  );
  }
  return sum;
}

NV chebyshev_psi(UV n)
{
  UV k;
  KAHAN_INIT(sum);

  for (k = log2floor(n); k > 0; k--) {
    KAHAN_SUM(sum, chebyshev_theta(rootint(n,k)));
  }
  return sum;
}

#if BITS_PER_WORD == 64
typedef struct {
  UV n;
  LNV theta;
} cheby_theta_t;
static const cheby_theta_t _cheby_theta[] = { /* >= quad math precision */
  { UVCONST(      67108864),LNVCONST(    67100507.6357700963903836828562472350035880) },
  { UVCONST(     100000000),LNVCONST(    99987730.0180220043832124342600487053812729) },
  { UVCONST(     134217728),LNVCONST(   134204014.5735572091791081610859055728165544) },
  { UVCONST(     268435456),LNVCONST(   268419741.6134308193112682817754501071404173) },
  { UVCONST(     536870912),LNVCONST(   536842885.8045763840625719515011160692495056) },
  { UVCONST(    1000000000),LNVCONST(   999968978.5775661447991262386023331863364793) },
  { UVCONST(    1073741824),LNVCONST(  1073716064.8860663337617909073555831842945484) },
  { UVCONST(    2147483648),LNVCONST(  2147432200.2475857676814950053003448716360822) },
  { UVCONST(    4294967296),LNVCONST(  4294889489.1735446386752045191908417183337361) },
  { UVCONST(    8589934592),LNVCONST(  8589863179.5654263491545135406516173629373070) },
  { UVCONST(   10000000000),LNVCONST(  9999939830.6577573841592219954033850595228736) },
  { UVCONST(   12884901888),LNVCONST( 12884796620.4324254952601520445848183460347362) },
  { UVCONST(   17179869184),LNVCONST( 17179757715.9924077567777285147574707468995695) },
  { UVCONST(   21474836480),LNVCONST( 21474693322.0998273969188369449626287713082943) },
  { UVCONST(   25769803776),LNVCONST( 25769579799.3751535467593954636665656772211515) },
  { UVCONST(   30064771072),LNVCONST( 30064545001.2305211029215168703433831598544454) },
  { UVCONST(   34359738368),LNVCONST( 34359499180.0126643918259085362039638823175054) },
  { UVCONST(   51539607552),LNVCONST( 51539356394.9531019037592855639826469993402730) },
  { UVCONST(   68719476736),LNVCONST( 68719165213.6369838785284711480925219076501720) },
  { UVCONST(   85899345920),LNVCONST( 85899083852.3471545629838432726841470626910905) },
  { UVCONST(  100000000000),LNVCONST( 99999737653.1074446948519125729820679772770146) },
  { UVCONST(  103079215104),LNVCONST(103079022007.113299711630969211422868856259124) },
  { UVCONST(  120259084288),LNVCONST(120258614516.787336970535750737470005730125261) },
  { UVCONST(  137438953472),LNVCONST(137438579206.444595884982301543904849253294539) },
  { UVCONST(  171798691840),LNVCONST(171798276885.585945657918751085729734540334501) },
  { UVCONST(  206158430208),LNVCONST(206158003808.160276853604927822609009916573462) },
  { UVCONST(  240518168576),LNVCONST(240517893445.995868018331936763125264759516048) },
  { UVCONST(  274877906944),LNVCONST(274877354651.045354829956619821889825596300686) },
  { UVCONST(  309237645312),LNVCONST(309237050379.850690561796126460858271984023198) },
  { UVCONST(  343597383680),LNVCONST(343596855806.595496630500062749631211394707114) },
  { UVCONST(  377957122048),LNVCONST(377956498560.227794386327526022452943941537993) },
  { UVCONST(  412316860416),LNVCONST(412316008796.349553568121442261222464590518293) },
  { UVCONST(  446676598784),LNVCONST(446675972485.936512329625489223180824947531484) },
  { UVCONST(  481036337152),LNVCONST(481035608287.572961376833237046440177624505864) },
  { UVCONST(  515396075520),LNVCONST(515395302740.633513931333424447688399032397200) },
  { UVCONST(  549755813888),LNVCONST(549755185085.539613556787409928561107952681488) },
  { UVCONST(  584115552256),LNVCONST(584115015741.698143680148976236958207248900725) },
  { UVCONST(  618475290624),LNVCONST(618474400071.621528348965919774195984612254220) },
  { UVCONST(  652835028992),LNVCONST(652834230470.583317059774197550110194348469358) },
  { UVCONST(  687194767360),LNVCONST(687193697328.927006867624832386534836384752774) },
  { UVCONST(  721554505728),LNVCONST(721553211683.605313067593521060195071837766347) },
  { UVCONST(  755914244096),LNVCONST(755913502349.878525212441903698096011352015192) },
  { UVCONST(  790273982464),LNVCONST(790273042590.053075430445971969285969445183076) },
  { UVCONST(  824633720832),LNVCONST(824633080997.428352876758261549475609957696369) },
  { UVCONST(  858993459200),LNVCONST(858992716288.318498931165663742671579465316192) },
  { UVCONST(  893353197568),LNVCONST(893352235882.851072417721659027263613727927680) },
  { UVCONST(  927712935936),LNVCONST(927711881043.628817668337317445143018372892386) },
  { UVCONST(  962072674304),LNVCONST(962071726126.508938539006575212272731584070786) },
  { UVCONST(  996432412672),LNVCONST(996431411588.361462717402562171913706963939018) },
  { UVCONST( 1099511627776),LNVCONST(1099510565082.05800550569923209414874779035972) },
  { UVCONST( 1168231104512),LNVCONST(1168230478726.83399452743801182220790107593115) },
  { UVCONST( 1236950581248),LNVCONST(1236949680081.02610603189530371762093291521116) },
  { UVCONST( 1305670057984),LNVCONST(1305668780900.04255251887970870257110498423202) },
  { UVCONST( 1374389534720),LNVCONST(1374388383792.63751003694755359184583212193880) },
  { UVCONST( 1443109011456),LNVCONST(1443107961091.80955496949174183091839841371227) },
  { UVCONST( 1511828488192),LNVCONST(1511827317611.91227277802426032456922797572429) },
  { UVCONST( 1580547964928),LNVCONST(1580546753969.30607547506449941085747942395437) },
  { UVCONST( 1649267441664),LNVCONST(1649265973878.75361554498682516738256005501353) },
  { UVCONST( 1717986918400),LNVCONST(1717985403764.24562741452793071287954107946922) },
  { UVCONST( 1786706395136),LNVCONST(1786704769212.04241689416220650800274263053933) },
  { UVCONST( 1855425871872),LNVCONST(1855425013030.54920163513184322741954734357404) },
  { UVCONST( 1924145348608),LNVCONST(1924143701943.02957992419280264060220278182021) },
  { UVCONST( 1992864825344),LNVCONST(1992863373568.84039296068619447120308124302086) },
  { UVCONST( 2061584302080),LNVCONST(2061583632335.91985095534685076604018573279204) },
  { UVCONST( 2130303778816),LNVCONST(2113122935598.01727180199783433992649406589029) },
  { UVCONST( 2199023255552),LNVCONST(2199021399611.18488312543276191461914978761981) },
  { UVCONST( 2267742732288),LNVCONST(2267740947106.05038218811506263712808318234921) },
  { UVCONST( 2336462209024),LNVCONST(2336460081480.34962633829077377680844065198307) },
  { UVCONST( 2405181685760),LNVCONST(2405179969505.38642629423585641169740223940265) },
  { UVCONST( 2473901162496),LNVCONST(2473899311193.37872375168104562948639924654178) },
  { UVCONST( 2542620639232),LNVCONST(2542619362554.88893589220737167756411653816418) },
  { UVCONST( 2611340115968),LNVCONST(2611338370515.94936514022501267847930999670553) },
  { UVCONST( 2680059592704),LNVCONST(2680057722824.52981820001574883706268873541107) },
  { UVCONST( 2748779069440),LNVCONST(2748777610452.18903407570165081726781627254885) },
  { UVCONST( 2817498546176),LNVCONST(2817497017165.31924616507392971415494161401775) },
  { UVCONST( 2886218022912),LNVCONST(2886216579432.32232322707222172612181994322081) },
  { UVCONST( 2954937499648),LNVCONST(2954936100812.97301730406598982753121204977388) },
  { UVCONST( 3023656976384),LNVCONST(3023654789503.82041452274471455184651411931920) },
  { UVCONST( 3298534883328),LNVCONST(3298533215621.76606493931157388037915263658637) },
  { UVCONST( 3573412790272),LNVCONST(3573411344351.74163523704886736624674718378131) },
  { UVCONST( 3848290697216),LNVCONST(3848288415701.82534219216958446478503907262807) },
  { UVCONST( 4123168604160),LNVCONST(4123166102085.86116301709394219323327831487542) },
  { UVCONST( 4398046511104),LNVCONST(4398044965678.05143041707871320554940671182665) },
  { UVCONST( 4672924418048),LNVCONST(4672922414672.04998927945349278916525727295687) },
  { UVCONST( 4947802324992),LNVCONST(4947800056419.04384937181159608905993450182729) },
  { UVCONST( 5222680231936),LNVCONST(5222678728087.69487334278665824384732845008859) },
  { UVCONST( 5497558138880),LNVCONST(5497555766573.55159115560501595606332808978878) },
  { UVCONST( 5772436045824),LNVCONST(5772433560746.27053256770924553245647027548204) },
  { UVCONST( 6047313952768),LNVCONST(6047310750621.24497633828761530843255989494448) },
  { UVCONST( 6322191859712),LNVCONST(6322189275338.39747421237532473168802646234745) },
  { UVCONST( 6597069766656),LNVCONST(6579887620000.56226807898107616294821989189226) },
  { UVCONST( 6871947673600),LNVCONST(6871945430474.61791600096091374271286154432006) },
  { UVCONST( 7146825580544),LNVCONST(7146823258390.34361980709600216319269118247416) },
  { UVCONST( 7421703487488),LNVCONST(7421700443390.35536080251964387835425662360121) },
  { UVCONST( 7696581394432),LNVCONST(7696578975137.73249441643024336954233783264803) },
  { UVCONST( 7971459301376),LNVCONST(7971457197928.90863708984184849978605273042512) },
  { UVCONST( 8246337208320),LNVCONST(8246333982863.77146812177727648999195989358960) },
  { UVCONST( 8521215115264),LNVCONST(8529802085075.55635100929751669785228592926043) },
  { UVCONST( 8796093022208),LNVCONST(8796089836425.34909684634625258535266362465034) },
  { UVCONST( 9345848836096),LNVCONST(9345845828116.77456046925508587313) },
  { UVCONST( 9895604649984),LNVCONST(9895601077915.26821447819584407150) },
  { UVCONST(10000000000000),LNVCONST(9999996988293.03419965318214160284) },
  { UVCONST(15000000000000),LNVCONST(14999996482301.7098815115045166858) },
  { UVCONST(20000000000000),LNVCONST(19999995126082.2286880312461318496) },
  { UVCONST(25000000000000),LNVCONST(24999994219058.4086216020475916538) },
  { UVCONST(30000000000000),LNVCONST(29999995531389.8454274046657200568) },
  { UVCONST(35000000000000),LNVCONST(34999992921190.8049427456456479005) },
  { UVCONST(40000000000000),LNVCONST(39999993533724.3168289589273168844) },
  { UVCONST(45000000000000),LNVCONST(44999993567606.9795798378256194424) },
  { UVCONST(50000000000000),LNVCONST(49999992543194.2636545758235373677) },
  { UVCONST(55000000000000),LNVCONST(54999990847877.2435105757522625171) },
  { UVCONST(60000000000000),LNVCONST(59999990297033.6261976055811111726) },
  { UVCONST(65000000000000),LNVCONST(64999990861395.5522142429859245014) },
  { UVCONST(70000000000000),LNVCONST(69999994316409.8717306521862685981) },
  { UVCONST(75000000000000),LNVCONST(74999990126219.8344899338374090165) },
  { UVCONST(80000000000000),LNVCONST(79999990160858.3042387288372250950) },
  { UVCONST(85000000000000),LNVCONST(84999987096970.5915212896832780715) },
  { UVCONST(90000000000000),LNVCONST(89999989501395.0738966599857919767) },
  { UVCONST(95000000000000),LNVCONST(94999990785908.6672552042792168144) },
  { UVCONST(100000000000000),LNVCONST(99999990573246.9785384070303475639) },
};
#define NCHEBY_VALS (sizeof(_cheby_theta)/sizeof(_cheby_theta[0]))
#endif

NV chebyshev_theta(UV n)
{
  uint16_t i = 0;
  UV tp, startn, seg_base, seg_low, seg_high;
  unsigned char* segment;
  void* ctx;
  LNV initial_sum, prod = LNV_ONE;
  KAHAN_INIT(sum);

  if (n < 500) {
    for (i = 1;  (tp = primes_tiny[i]) <= n; i++) {
      KAHAN_SUM(sum, loglnv(tp));
    }
    return sum;
  }

#if defined NCHEBY_VALS
  if (n >= _cheby_theta[0].n) {
    for (i = 1; i < NCHEBY_VALS; i++)
      if (n < _cheby_theta[i].n)
        break;
    startn = _cheby_theta[i-1].n;
    initial_sum = _cheby_theta[i-1].theta;
  } else
#endif
  {
    KAHAN_SUM(sum, loglnv(2*3*5*7*11*13));
    startn = 17;
    initial_sum = 0;
  }

  ctx = start_segment_primes(startn, n, &segment);
#if 0
  while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
    START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high ) {
      KAHAN_SUM(sum, loglnv(p));
    } END_DO_FOR_EACH_SIEVE_PRIME
  }
#else
  while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
    START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high ) {
      prod *= (LNV) p;
      if (++i >= (LNV_IS_QUAD ? 64 : 8)) {
        KAHAN_SUM(sum, loglnv(prod));
        prod = LNV_ONE;
        i = 0;
      }
    } END_DO_FOR_EACH_SIEVE_PRIME
  }
  if (prod > 1.0) { KAHAN_SUM(sum, loglnv(prod));  prod = LNV_ONE; }
#endif
  end_segment_primes(ctx);

  if (initial_sum > 0) KAHAN_SUM(sum, initial_sum);
  return sum;
}



/*
 * See:
 *  "Multiple-Precision Exponential Integral and Related Functions"
 *      by David M. Smith
 *  "On the Evaluation of the Complex-Valued Exponential Integral"
 *      by Vincent Pegoraro and Philipp Slusallek
 *  "Numerical Recipes" 3rd edition
 *      by William H. Press et al.
 *  "Rational Chevyshev Approximations for the Exponential Integral E_1(x)"
 *      by W. J. Cody and Henry C. Thacher, Jr.
 *  "High-precision Computation of Uniform Asymptotic Expansions for Special Functions"
 *      by Guillermo Navas-Palencia (2019)
 *
 * Any mistakes here are mine.  This code has not been rigorously verified.
 * Alternates: Navas-Palencia, Boost, MPFR, Pari/GP, Arb.
 *
 * We are trying to get close to maximum precision for all x with double, long
 * double, and quadmath.  Hence the rational Chebyshev approximations should
 * not be used with quadmath (unless they are are modified).
 *
 * Performance, i7-6700HQ, 2.6GHz, 1e-9 to 1000 step 0.001
 * range x > 0:
 *   0.22 microseconds, NV = double         max rel error 1.4e-14
 *   0.19 microseconds, NV = long double    max rel error 4.3e-17
 *  18.97 microseconds, NV = quad           max rel error 4.4e-32
 * range x < 0:
 *   0.18 microseconds, NV = double         max rel error 1.4e-14
 *   0.15 microseconds, NV = long double    max rel error 1.2e-17
 *   9.31 microseconds, NV = quad           max rel error 1.7e-32
 *
 * max error is at 0.372.  The relative error profile for double precision
 * is essentially identical to the Navas-Palencia expintei(x) function.
 * Using long double on x86 greatly improves the results with no time penalty.
 * Using quadmath gives improved results at a substantial time penalty.
 */

static LNV const euler_mascheroni = LNVCONST(0.57721566490153286060651209008240243104215933593992);
static LNV const li2 = LNVCONST(1.045163780117492784844588889194613136522615578151);

/* Rational Chebyshev approximation (Cody, Thacher), good for -1 < x < 0 */
static LNV _ei_chebyshev_neg(const LNV x) {
  static const LNV C6p[7] = { LNVCONST(-148151.02102575750838086),
                              LNVCONST( 150260.59476436982420737),
                              LNVCONST(  89904.972007457256553251),
                              LNVCONST(  15924.175980637303639884),
                              LNVCONST(   2150.0672908092918123209),
                              LNVCONST(    116.69552669734461083368),
                              LNVCONST(      5.0196785185439843791020) };
  static const LNV C6q[7] = { LNVCONST( 256664.93484897117319268),
                              LNVCONST( 184340.70063353677359298),
                              LNVCONST(  52440.529172056355429883),
                              LNVCONST(   8125.8035174768735759866),
                              LNVCONST(    750.43163907103936624165),
                              LNVCONST(     40.205465640027706061433),
                              LNVCONST(      1.0000000000000000000000) };
  LNV sumn = C6p[0]-x*(C6p[1]-x*(C6p[2]-x*(C6p[3]-x*(C6p[4]-x*(C6p[5]-x*C6p[6])))));
  LNV sumd = C6q[0]-x*(C6q[1]-x*(C6q[2]-x*(C6q[3]-x*(C6q[4]-x*(C6q[5]-x*C6q[6])))));
  return loglnv(-x) - sumn/sumd;
}
/* Cody / Thacher rational Chebyshev for x > 24 */
static LNV _ei_chebyshev_pos24(const LNV x) {
  static const LNV P2[10] = {
      LNVCONST( 1.75338801265465972390E02),
      LNVCONST(-2.23127670777632409550E02),
      LNVCONST(-1.81949664929868906455E01),
      LNVCONST(-2.79798528624305389340E01),
      LNVCONST(-7.63147701620253630855E00),
      LNVCONST(-1.52856623636929636839E01),
      LNVCONST(-7.06810977895029358836E00),
      LNVCONST(-5.00006640413131002475E00),
      LNVCONST(-3.00000000320981265753E00),
      LNVCONST( 1.00000000000000485503E00) };
  static const LNV Q2[9] = {
      LNVCONST( 3.97845977167414720840E04),
      LNVCONST( 3.97277109100414518365E00),
      LNVCONST( 1.37790390235747998793E02),
      LNVCONST( 1.17179220502086455287E02),
      LNVCONST( 7.04831847180424675988E01),
      LNVCONST(-1.20187763547154743238E01),
      LNVCONST(-7.99243595776339741065E00),
      LNVCONST(-2.99999894040324959612E00),
      LNVCONST( 1.99999999999048104167E00) };
  LNV invx = LNV_ONE / x, frac = 0.0;
  uint32_t n;
  for (n = 0; n <= 8; n++)
    frac = Q2[n] / (P2[n] + x + frac);
  frac += P2[9];
  return explnv(x) * (invx + invx*invx*frac);
}
#if 0
/* Continued fraction, good for x < -1 */
static LNV _ei_cfrac_neg(const LNV x) {
  LNV lc = 0, ld = LNV_ONE / (LNV_ONE - x);
  LNV val = ld * (-explnv(x));
  uint32_t n;
  for (n = 1; n <= 20000; n++) {
    LNV old, t, n2 = n * n;
    t = (LNV)(2*n + 1) - x;
    lc = LNV_ONE / (t - n2 * lc);
    ld = LNV_ONE / (t - n2 * ld);
    old = val;
    val *= ld/lc;
    if ( fabslnv(val-old) <= LNV_EPSILON*fabslnv(val) )
      break;
  }
  return val;
}
#endif
/* eint_v using Laguerre series, Navas-Palencia (2019). */
static LNV _eintv_laguerre_series(const LNV v, const LNV x) {
  LNV L_k = 1.0, L_k1 = x + v;
  LNV q, r, u = LNV_ONE, d = LNV_ONE;
  uint32_t k;
  KAHAN_INIT(sum);
  KAHAN_SUM(sum, (LNV_ONE/L_k1));
  for (k = 1; k < 500; k++) {
    u *= v + k - 1;
    d *= 1 + k;
    q = L_k1 * (x + 2*k + v) / (k + 1)  -  L_k * (k + v - 1) / (k + 1);
    r = u / (d * (q * L_k1));
    KAHAN_SUM(sum, r);
    L_k = L_k1;
    L_k1 = q;
    if (fabslnv(r) < 0.1 * LNV_EPSILON)
      break;
  }
  return sum * explnv(-x);
}
/* Convergent series for small negative x through medium positive x */
static LNV _ei_series_convergent(LNV const x) {
  LNV val, term, fact_n = x;
  uint32_t n;
  KAHAN_INIT(sum);
  for (n = 2; n <= 400; n++) {
    LNV invn = LNV_ONE / n;
    fact_n *= (LNV)x * invn;
    term = fact_n * invn;
    KAHAN_SUM(sum, term);
    /* printf("C  after adding %.20Lf, val = %.20Lf\n", term, sum); */
    if (fabslnv(term) < LNV_EPSILON*fabslnv(sum)) break;
  }
  KAHAN_SUM(sum, euler_mascheroni);
  KAHAN_SUM(sum, loglnv(fabslnv(x)));
  KAHAN_SUM(sum, x);
  return sum;
}
/* Asymptotic divergent series, for large positive x */
static LNV _ei_series_divergent(LNV const x) {
  LNV val, invx = LNV_ONE / x, term = invx;
  unsigned int n;
  KAHAN_INIT(sum);
  for (n = 2; n <= 400; n++) {
    LNV last_term = term;
    term = term * ( (LNV)n * invx );
    if (term < LNV_EPSILON*sum) break;
    if (term < last_term) {
      KAHAN_SUM(sum, term);
      /* printf("A  after adding %.20llf, sum = %.20llf\n", term, sum); */
    } else {
      KAHAN_SUM(sum, (-last_term/1.07) );
      /* printf("A  after adding %.20llf, sum = %.20llf\n", -last_term/1.07, sum); */
      break;
    }
  }
  KAHAN_SUM(sum, invx);
  KAHAN_SUM(sum, LNV_ONE);
  return explnv(x) * sum * invx;
}

NV Ei(NV x) {
  if (x == 0) croak("Invalid input to ExponentialIntegral:  x must be != 0");
  /* Protect against messed up rounding modes */
  if (x >=  12000) return INFINITY;
  if (x <= -12000) return 0;

  if (x < 0) {
    if (x >= -1.0 && !LNV_IS_QUAD) return _ei_chebyshev_neg(x);
    else if (x < -0.80)            return -_eintv_laguerre_series(1, -x);
    else                           return _ei_series_convergent(x);
  } else {
    if (x < (-2 * loglnv(LNV_EPSILON)))         return _ei_series_convergent(x);
    if (x >= 24 && (!LNV_IS_QUAD || x <= 43.2)) return _ei_chebyshev_pos24(x);
    else                                        return _ei_series_divergent(x);
  }
}

NV Li(NV x) {
  if (x == 0) return 0;
  if (x == 1) return -INFINITY;
  if (x == 2) return li2;
  if (x < 0) croak("Invalid input to LogarithmicIntegral:  x must be >= 0");
  if (x >= NV_MAX) return INFINITY;

  /* Calculate directly using Ramanujan's series. */
  if (x > 1) {
    const LNV logx = loglnv(x);
    LNV sum = 0, inner_sum = 0, old_sum, factorial = 1, power2 = 1;
    LNV q, p = -1;
    int k = 0, n = 0;

    for (n = 1, k = 0; n < 200; n++) {
      factorial *= n;
      p *= -logx;
      q = factorial * power2;
      power2 *= 2;
      for (; k <= (n - 1) / 2; k++)
        inner_sum += LNV_ONE / (2 * k + 1);
      old_sum = sum;
      sum += (p / q) * inner_sum;
      if (fabslnv(sum - old_sum) <= LNV_EPSILON) break;
    }
    return euler_mascheroni + loglnv(logx) + sqrtlnv(x) * sum;
  }

  return Ei(loglnv(x));
}

static long double ld_inverse_li(long double lx) {
  int i;
  long double t, term, old_term = 0;
  /* Iterate Halley's method until error grows. */
  t = (lx <= 2)  ?  2  :  lx * logl(lx);
  for (i = 0; i < 4; i++) {
    long double dn = Li(t) - lx;
    term = dn*logl(t) / (1.0L + dn/(2*t));
    if (i > 0 && fabsl(term) >= fabsl(old_term)) { t -= term/4; break; }
    old_term = term;
    t -= term;
  }
  return t;
}

UV inverse_li(UV x) {
  UV r, i;
  long double lx = (long double) x;

  if (x <= 2) return x + (x > 0);
  r = (UV) ceill( ld_inverse_li(lx) );
  /* Meet our more stringent goal of an exact answer. */
  i = (x > 4e16) ? 2048 : 128;
  if (Li(r-1) >= lx) {
    while (Li(r-i) >= lx) r -= i;
    for (i = i/2; i > 0; i /= 2)
      if (Li(r-i) >= lx) r -= i;
  } else {
    while (Li(r+i-1) < lx) r += i;
    for (i = i/2; i > 0; i /= 2)
      if (Li(r+i-1) < lx) r += i;
  }
  return r;
}

static long double ld_inverse_R(long double lx) {
  int i;
  long double t, dn, term, old_term = 0;

  /* Rough estimate */
  if (lx <= 3.5) {
    t = lx + 2.24*(lx-1)/2;
  } else {
    t = lx * logl(lx);
    if      (lx <   50) { t *= 1.2; }
    else if (lx < 1000) { t *= 1.15; }
    else {   /* use inverse Li (one iteration) for first inverse R approx */
      dn = Li(t) - lx;
      term = dn * logl(t) / (1.0L + dn/(2*t));
      t -= term;
    }
  }
  /* Iterate 1-n rounds of Halley, usually only 3 needed. */
  for (i = 0; i < 100; i++) {
    dn = RiemannR(t, 1e-12) - lx;
#if 1  /* Use f(t) = li(t) for derivatives */
    term = dn * logl(t) / (1.0L + dn/(2*t));
#else  /* Use f(t) = li(t) - li(sqrt(t))/2 for derivatives */
    long double logt = logl(t);
    long double sqrtt = sqrtl(t);
    long double FA = 2 * sqrtt * logt;
    long double FB = 2 * sqrtt - 1;
    long double ifz = FA / FB;
    long double iffz = (logt - 2*FB) / (2 * sqrtt * FA * FA * FA * FA);
    term = dn * ifz * (1.0L - dn * iffz);
#endif
    if (i > 0 && fabsl(term) >= fabsl(old_term)) { t -= term/4; break; }
    old_term = term;
    t -= term;
  }
  return t;
}

UV inverse_R(UV x) {
  if (x < 2) return x + (x > 0);
  return (UV) ceill( ld_inverse_R( (long double) x) );
}


/*
 * Storing the first 10-20 Zeta values makes sense.  Past that it is purely
 * to avoid making the call to generate them ourselves.  We could cache the
 * calculated values. These all have 1 subtracted from them.  */
static const long double riemann_zeta_table[] = {
  0.6449340668482264364724151666460251892L,  /* zeta(2) */
  0.2020569031595942853997381615114499908L,
  0.0823232337111381915160036965411679028L,
  0.0369277551433699263313654864570341681L,
  0.0173430619844491397145179297909205279L,
  0.0083492773819228268397975498497967596L,
  0.0040773561979443393786852385086524653L,
  0.0020083928260822144178527692324120605L,
  0.0009945751278180853371459589003190170L,
  0.0004941886041194645587022825264699365L,
  0.0002460865533080482986379980477396710L,
  0.0001227133475784891467518365263573957L,
  0.0000612481350587048292585451051353337L,
  0.0000305882363070204935517285106450626L,
  0.0000152822594086518717325714876367220L,
  0.0000076371976378997622736002935630292L,  /* zeta(17)  Past here all we're */
  0.0000038172932649998398564616446219397L,  /* zeta(18)  getting is speed.   */
  0.0000019082127165539389256569577951013L,
  0.0000009539620338727961131520386834493L,
  0.0000004769329867878064631167196043730L,
  0.0000002384505027277329900036481867530L,
  0.0000001192199259653110730677887188823L,
  0.0000000596081890512594796124402079358L,
  0.0000000298035035146522801860637050694L,
  0.0000000149015548283650412346585066307L,
  0.0000000074507117898354294919810041706L,
  0.0000000037253340247884570548192040184L,
  0.0000000018626597235130490064039099454L,
  0.0000000009313274324196681828717647350L,
  0.0000000004656629065033784072989233251L,
  0.0000000002328311833676505492001455976L,
  0.0000000001164155017270051977592973835L,
  0.0000000000582077208790270088924368599L,
  0.0000000000291038504449709968692942523L,
  0.0000000000145519218910419842359296322L,
  0.0000000000072759598350574810145208690L,
  0.0000000000036379795473786511902372363L,
  0.0000000000018189896503070659475848321L,
  0.0000000000009094947840263889282533118L,
  0.0000000000004547473783042154026799112L,
  0.0000000000002273736845824652515226821L,
  0.0000000000001136868407680227849349105L,
  0.0000000000000568434198762758560927718L,
  0.0000000000000284217097688930185545507L,
  0.0000000000000142108548280316067698343L,
  0.00000000000000710542739521085271287735L,
  0.00000000000000355271369133711367329847L,
  0.00000000000000177635684357912032747335L,
  0.000000000000000888178421093081590309609L,
  0.000000000000000444089210314381336419777L,
  0.000000000000000222044605079804198399932L,
  0.000000000000000111022302514106613372055L,
  0.0000000000000000555111512484548124372374L,
  0.0000000000000000277555756213612417258163L,
  0.0000000000000000138777878097252327628391L,
};
#define NPRECALC_ZETA (sizeof(riemann_zeta_table)/sizeof(riemann_zeta_table[0]))

/* Riemann Zeta on the real line, with 1 subtracted.
 * Compare to Math::Cephes zetac.  Also zeta with q=1 and subtracting 1.
 *
 * The Cephes zeta function uses a series (2k)!/B_2k which converges rapidly
 * and has a very wide range of values.  We use it here for some values.
 *
 * Note: Calculations here are done on long doubles and we try to generate as
 *       much accuracy as possible.  They will get returned to Perl as an NV,
 *       which is typically a 64-bit double with 15 digits.
 *
 * For values 0.5 to 5, this code uses the rational Chebyshev approximation
 * from Cody and Thacher.  This method is extraordinarily fast and very
 * accurate over its range (slightly better than Cephes for most values).  If
 * we had quad floats, we could use the 9-term polynomial.
 */
long double ld_riemann_zeta(long double x) {
  int i;

  if (x < 0)  croak("Invalid input to RiemannZeta:  x must be >= 0");
  if (x == 1) return INFINITY;

  if (x == (unsigned int)x) {
    int k = x - 2;
    if ((k >= 0) && (k < (int)NPRECALC_ZETA))
      return riemann_zeta_table[k];
  }

  /* Cody / Thacher rational Chebyshev approximation for small values */
  if (x >= 0.5 && x <= 5.0) {
    static const long double C8p[9] = { 1.287168121482446392809e10L,
                                        1.375396932037025111825e10L,
                                        5.106655918364406103683e09L,
                                        8.561471002433314862469e08L,
                                        7.483618124380232984824e07L,
                                        4.860106585461882511535e06L,
                                        2.739574990221406087728e05L,
                                        4.631710843183427123061e03L,
                                        5.787581004096660659109e01L };
    static const long double C8q[9] = { 2.574336242964846244667e10L,
                                        5.938165648679590160003e09L,
                                        9.006330373261233439089e08L,
                                        8.042536634283289888587e07L,
                                        5.609711759541920062814e06L,
                                        2.247431202899137523543e05L,
                                        7.574578909341537560115e03L,
                                       -2.373835781373772623086e01L,
                                        1.000000000000000000000L    };
    long double sumn = C8p[0]+x*(C8p[1]+x*(C8p[2]+x*(C8p[3]+x*(C8p[4]+x*(C8p[5]+x*(C8p[6]+x*(C8p[7]+x*C8p[8])))))));
    long double sumd = C8q[0]+x*(C8q[1]+x*(C8q[2]+x*(C8q[3]+x*(C8q[4]+x*(C8q[5]+x*(C8q[6]+x*(C8q[7]+x*C8q[8])))))));
    long double sum = (sumn - (x-1)*sumd) / ((x-1)*sumd);
    return sum;
  }

  if (x > 17000.0)
    return 0.0;

#if 0
  {
    KAHAN_INIT(sum);
    /* Simple defining series, works well. */
    for (i = 5; i <= 1000000; i++) {
      long double term = powl(i, -x);
      KAHAN_SUM(sum, term);
      if (term < LDBL_EPSILON*sum) break;
    }
    KAHAN_SUM(sum, powl(4, -x) );
    KAHAN_SUM(sum, powl(3, -x) );
    KAHAN_SUM(sum, powl(2, -x) );
    return sum;
  }
#endif

  /* The 2n!/B_2k series used by the Cephes library. */
  {
    /* gp/pari:
     *   for(i=1,13,printf("%.38g\n",(2*i)!/bernreal(2*i)))
     * MPU:
     *   use bignum;
     *   say +(factorial(2*$_)/bernreal(2*$_))->bround(38) for 1..13;
     */
    static const long double A[] = {
      12.0L,
     -720.0L,
      30240.0L,
     -1209600.0L,
      47900160.0L,
     -1892437580.3183791606367583212735166425L,
      74724249600.0L,
     -2950130727918.1642244954382084600497650L,
      116467828143500.67248729113000661089201L,
     -4597978722407472.6105457273596737891656L,
      181521054019435467.73425331153534235290L,
     -7166165256175667011.3346447367083352775L,
      282908877253042996618.18640556532523927L,
    };
    long double a, b, s, t;
    const long double w = 10.0;
    s = 0.0;
    b = 0.0;
    for (i = 2; i < 11; i++) {
      b = powl( i, -x );
      s += b;
      if (fabsl(b) < fabsl(LDBL_EPSILON * s))
        return s;
    }
    s = s + b*w/(x-1.0) - 0.5 * b;
    a = 1.0;
    for (i = 0; i < 13; i++) {
      long double k = 2*i;
      a *= x + k;
      b /= w;
      t = a*b/A[i];
      s = s + t;
      if (fabsl(t) < fabsl(LDBL_EPSILON * s))
        break;
      a *= x + k + 1.0;
      b /= w;
    }
    return s;
  }
}

long double RiemannR(long double x, long double eps) {
  long double part_term, term, flogx, ki, old_sum;
  unsigned int k;
  KAHAN_INIT(sum);

  if (x <= 0) croak("Invalid input to RiemannR:  x must be > 0");
  if (eps < LDBL_EPSILON) eps = LDBL_EPSILON;

  if (x > 1e19) {
    const signed char* amob = range_moebius(0, 100);
    KAHAN_SUM(sum, Li(x));
    for (k = 2; k <= 100; k++) {
      if (amob[k] == 0) continue;
      ki = 1.0L / (long double) k;
      part_term = powl(x,ki);
      if (part_term > LDBL_MAX) return INFINITY;
      term = amob[k] * ki * Li(part_term);
      old_sum = sum;
      KAHAN_SUM(sum, term);
      if (fabsl(sum - old_sum) <= eps) break;
    }
    Safefree(amob);
    return sum;
  }

  KAHAN_SUM(sum, 1.0);
  flogx = logl(x);
  part_term = 1;

  for (k = 1; k <= 10000; k++) {
    ki = (k-1 < NPRECALC_ZETA) ? riemann_zeta_table[k-1] : ld_riemann_zeta(k+1);
    part_term *= flogx / k;
    term = part_term / (k + k * ki);
    old_sum = sum;
    KAHAN_SUM(sum, term);
    /* printf("R %5d after adding %.18Lg, sum = %.19Lg (%Lg)\n", k, term, sum, fabsl(sum-old_sum)); */
    if (fabsl(sum - old_sum) <= eps) break;
  }

  return sum;
}

/* Options for LambertW initial approximation:
 *
 * - Four regions, we used before:
 *      Pade(3,2), Winitzki 2003, Vargas 2013, Corless 1993
 *   Has issues near -1/e but ok around zero.
 *
 * - Iacono and Boyd (2017).  Very simple function over whole range.
 *   Doesn't work right very near -1/e and around zero.
 *
 * - Vazquez-Leal et al. (2019).  Divides into four regions, power
 *   series for each.  Great results.  Also has issues near -1/e and zero.
 *
 * We use known solutions for near -1/e and around zero.  See Fukushima (2013)
 * and Johannson (2017,2020) for lots of discussion and solutions.
 * Use Vazquez-Leal (PSEM Approximations) for the rest.
 */
static long double _lambertw_approx(long double x) {
  long double w, k1, k2, k3;

  if (x < -0.312) {
    /* Use Puiseux series, e.g. Verberic 2009, Boost, Johannson (2020). */
    /* Near the branch point.  See Fukushima (2013) section 2.5. */
    k2 = 2.0L * (1.0L + 2.7182818284590452353603L * x);
    if (k2 <= 0) return -1.0L + 1*LDBL_EPSILON;
    k1 = sqrtl(k2);
    w = -1.0L + (1.0L + (-1.0L/3.0L + (11.0L/72.0L + (-43.0L/540.0L + (769.0L/17280.0L + (-221.0L/8505.0L + (680863.0L/43545600.0L + (-1963.0L/204120.0L + 226287557.0L/37623398400.0L
    * k1) * k1) * k1) * k1) * k1) * k1) * k1) * k1) * k1;

  } else if (x > -0.14 && x < 0.085) {
    /* Around zero.  See Fukushima (2013) section 2.6. */
    w = (1.0L + (-1.0L + (3.0L/2.0L + (-8.0L/3.0L + (125.0L/24.0L + (-54.0L/5.0L + (16807.0L/720.0L + (-16384.0L/315.0L + 531441.0L/4480.0L
        * x) * x) * x) * x) * x) * x) * x) * x) * x;

  } else if (x < 1) {
    /* This and the rest from Vazquez-Leal et al. (2019). */
    k1 = sqrtl(1.0L + 2.7182818284590452353603L * x);
    k2 = 0.33333333333333333333333L + 0.70710678118654752440084L / k1 - 0.058925565098878960366737L * k1 +
         (x + 0.36787944117144L) * (0.050248489761611L + (0.11138904851051 + 0.040744556245195L * x) * x)
         /
         (1.0L + (2.7090878606183L + (1.5510922597820L + 0.095477712183841L * x) * x) * x);
    w = -(k2-1)/k2;

  } else if (x < 40) {
    k1 = 1.0L + (5.950065500550155L + (13.96586471370701L + (10.52192021050505L + (3.065294254265870L + 0.1204576876518760L * x) * x) * x) * x) * x;
    w = 0.1600049638651493L * logl(k1);

  } else if (x < 20000) {
    k1 = 1.0L + (-3.16866642511229e11L + (3.420439800038598e10L +
         (-1.501433652432257e9L + (3.44887729947585e7L + (-4.453783741137856e5L +
         (3257.926478908996L + (-10.82545259305382L + (0.6898058947898353e-1L +
         0.4703653406071575e-4L * x) * x) * x) * x) * x) * x) * x) * x) * x;
    w = 0.9898045358731312e-1L * logl(k1);

  } else {
    k1 = 1.0L / (1.0L + logl(1.0L + x));
    k2 = 1.0L / k1;
    k3 = logl(k2);
    w = k2-1-k3+(1+k3+(-1/2+(1/2)*k3*k3 +(-1/6+(-1+(-1/2+
        (1/3) * k3) * k3) * k3) * k1) * k1) * k1;
  }
  return w;
}

NV lambertw(NV x) {
  long double w;
  int i;

  if (x < -0.36787944117145L)
    croak("Invalid input to LambertW:  x must be >= -1/e");
  if (x == 0.0L) return 0.0L;

  /* Estimate initial value */
  w = _lambertw_approx(x);

  /* TODO: this section might not be best for quad precision */
  /* If input is too small, return .99999.... */
  /* if (w <= -1.0L) return -1.0L + LDBL_EPSILON; */
  /* For very small inputs, don't iterate, return approx directly. */
  if (x < -0.36768) return w;

#if 0  /* Halley */
  long double lastw = w;
  for (i = 0; i < 100; i++) {
    long double ew = expl(w);
    long double wew = w * ew;
    long double wewx = wew - x;
    long double w1 = w + 1;
    w = w - wewx / (ew * w1 - (w+2) * wewx/(2*w1));
    if (w != 0.0L && fabsl((w-lastw)/w) <= 8*LDBL_EPSILON) break;
    lastw = w;
  }
#else  /* Fritsch, see Veberic 2009.  1-2 iterations are enough. */
  for (i = 0; i < 6 && w != 0.0L; i++) {
    long double w1 = 1 + w;
    long double zn = logl((long double)x/w) - w;
    long double qn = 2 * w1 * (w1+(2.0L/3.0L)*zn);
    long double en = (zn/w1) * (qn-zn)/(qn-2.0L*zn);
    /* w *= 1.0L + en;  if (fabsl(en) <= 16*LDBL_EPSILON) break; */
    long double wen = w * en;
    if (isnan(wen)) return 0;
    w += wen;
    if (fabsl(wen) <= 64*LDBL_EPSILON) break;
  }
#endif

#if LNV_IS_QUAD /* For quadmath, one high precision correction */
  if (w != LNV_ZERO) {
    LNV lw = w;
    LNV w1 = LNV_ONE + lw;
    LNV zn = loglnv((LNV)x/lw) - lw;
    LNV qn = LNVCONST(2.0) * w1 * (w1+(LNVCONST(2.0)/LNVCONST(3.0))*zn);
    LNV en = (zn/w1) * (qn-zn)/(qn-LNVCONST(2.0)*zn);
    return lw + lw * en;
  }
#endif

  /* With long double = 64-bit double, we have 15 digits precision
   * near the branch point, and 16 over the rest of the range.
   * With long double = x86 extended precision, we have over 17 digits
   * over the entire range.
   * Correcting to the exact LDBL_EPSILON does not improve this. */

  return w;
}

#if HAVE_STD_U64
  #define U64T uint64_t
#else
  #define U64T UV
#endif

/* Spigot from Arndt, Haenel, Winter, and Flammenkamp. */
/* Modified for larger digits and rounding by Dana Jacobsen */
char* pidigits(int digits)
{
  char* out;
  uint32_t *a, b, c, d, e, g, i, d4, d3, d2, d1;
  uint32_t const f = 10000;
  U64T d64;  /* 64-bit intermediate for 2*2*10000*b > 2^32 (~30k digits) */

  if (digits <= 0) return 0;
  if (digits >= 1 && digits <= DBL_DIG && digits <= 18) {
    Newz(0, out, 19, char);
    (void)sprintf(out, "%.*lf", (digits-1), 3.141592653589793238);
    return out;
  }
  digits++;   /* For rounding */
  c = 14*(digits/4 + 2);
  New(0, out, digits+5+1, char);
  *out++ = '3';  /* We'll turn "31415..." into "3.1415..." */
  New(0, a, c, uint32_t);
  for (b = 0; b < c; b++)  a[b] = 2000;

  d = i = 0;
  while ((b = c -= 14) > 0 && i < (uint32_t)digits) {
    d = e = d % f;
    if (b > 107000) {  /* Use 64-bit intermediate while necessary. */
      for (d64 = d; --b > 107000; ) {
        g = (b << 1) - 1;
        d64 = d64 * b  +  f * (U64T)a[b];
        a[b] = d64 % g;
        d64 /= g;
      }
      d = d64;
      b++;
    }
    while (--b > 0) {
      g = (b << 1) - 1;
      d = d * b  +  f * a[b];
      a[b] = d % g;
      d /= g;
    }
    /* sprintf(out+i, "%04d", e+d/f);   i += 4; */
    d4 = e + d/f;
    if (d4 > 9999) {
      d4 -= 10000;
      out[i-1]++;
      for (b=i-1; out[b] == '0'+1; b--) { out[b]='0'; out[b-1]++; }
    }
    d3 = d4/10;  d2 = d3/10;  d1 = d2/10;
    out[i++] = '0' + d1;
    out[i++] = '0' + d2-d1*10;
    out[i++] = '0' + d3-d2*10;
    out[i++] = '0' + d4-d3*10;
  }
  Safefree(a);
  if (out[digits-1] >= '5') out[digits-2]++;  /* Round */
  for (i = digits-2; out[i] == '9'+1; i--)    /* Keep rounding */
    { out[i] = '0';  out[i-1]++; }
  digits--;  /* Undo the extra digit we used for rounding */
  out[digits] = '\0';
  *out-- = '.';
  return out;
}

/* 1. Perform signed integer validation on b/blen.
 * 2. Compare to a/alen using min or max based on first arg.
 * 3. Return 0 to select a, 1 to select b.
 */
int strnum_minmax(int min, char* a, STRLEN alen, char* b, STRLEN blen)
{
  int aneg, bneg;
  STRLEN i;
  /* a is checked, process b */
  if (b == 0 || blen == 0) croak("Parameter must be a positive integer");
  bneg = (b[0] == '-');
  if (b[0] == '-' || b[0] == '+') { b++; blen--; }
  while (blen > 0 && *b == '0') { b++; blen--; }
  for (i = 0; i < blen; i++)
    if (!isDIGIT(b[i]))
      break;
  if (blen == 0 || i < blen)
    croak("Parameter must be a positive integer");

  if (a == 0) return 1;

  aneg = (a[0] == '-');
  if (a[0] == '-' || a[0] == '+') { a++; alen--; }
  while (alen > 0 && *a == '0') { a++; alen--; }

  if (aneg != bneg)  return  min  ?  (bneg == 1)  :  (aneg == 1);
  if (aneg == 1)  min = !min;
  if (alen != blen)  return  min  ?  (alen > blen) :  (blen > alen);

  for (i = 0; i < blen; i++)
    if (a[i] != b[i])
      return  min  ?  (a[i] > b[i])  :  (b[i] > a[i]);
  return 0; /* equal */
}

int from_digit_string(UV* rn, const char* s, int base)
{
  UV max, n = 0;
  int i, len;

  /* Skip leading -/+ and zeros */
  if (s[0] == '-' || s[0] == '+') s++;
  while (s[0] == '0') s++;

  len = strlen(s);
  max = (UV_MAX-base+1)/base;

  for (i = 0; i < len; i++) {
    const char c = s[i];
    int d = !isalnum(c) ? 255 : (c <= '9') ? c-'0' : (c <= 'Z') ? c-'A'+10 : c-'a'+10;
    if (d >= base) croak("Invalid digit for base %d", base);
    if (n > max) return 0;   /* Overflow */
    n = n * base + d;
  }
  *rn = n;
  return 1;
}

int from_digit_to_UV(UV* rn, UV* r, int len, int base)
{
  UV d, n = 0;
  int i;
  if (len < 0 || len > BITS_PER_WORD)
    return 0;
  for (i = 0; i < len; i++) {
    d = r[i];
    if (n > (UV_MAX-d)/base) break;  /* overflow */
    n = n * base + d;
  }
  *rn = n;
  return (i >= len);
}


int from_digit_to_str(char** rstr, UV* r, int len, int base)
{
  char *so, *s;
  int i;

  if (len < 0 || !(base == 2 || base == 10 || base == 16)) return 0;

  if (r[0] >= (UV) base) return 0;  /* TODO: We don't apply extended carry */

  New(0, so, len + 3, char);
  s = so;
  if (base == 2 || base == 16) {
    *s++ = '0';
    *s++ = (base == 2) ? 'b' : 'x';
  }
  for (i = 0; i < len; i++) {
    UV d = r[i];
    s[i] = (d < 10) ? '0'+d : 'a'+d-10;
  }
  s[len] = '\0';
  *rstr = so;
  return 1;
}

int to_digit_array(int* bits, UV n, int base, int length)
{
  int d;

  if (base < 2 || length > 128) return -1;

  if (base == 2) {
    for (d = 0; n; n >>= 1)
      bits[d++] = n & 1;
  } else {
    for (d = 0; n; n /= base)
      bits[d++] = n % base;
  }
  if (length < 0) length = d;
  while (d < length)
    bits[d++] = 0;
  return length;
}

int to_digit_string(char* s, UV n, int base, int length)
{
  int digits[128];
  int i, len = to_digit_array(digits, n, base, length);

  if (len < 0) return -1;
  if (base > 36) croak("invalid base for string: %d", base);

  for (i = 0; i < len; i++) {
    int dig = digits[len-i-1];
    s[i] = (dig < 10) ? '0'+dig : 'a'+dig-10;
  }
  s[len] = '\0';
  return len;
}

int to_string_128(char str[40], IV hi, UV lo)
{
  int i, slen = 0, isneg = 0;

  if (hi < 0) {
    isneg = 1;
    hi = -(hi+1);
    lo = UV_MAX - lo + 1;
  }
#if BITS_PER_WORD == 64 && HAVE_UINT128
  {
    uint128_t dd, sum = (((uint128_t) hi) << 64) + lo;
    do {
      dd = sum / 10;
      str[slen++] = '0' + (sum - dd*10);
      sum = dd;
    } while (sum);
  }
#else
  {
    UV d, r;
    uint32_t a[4];
    a[0] = hi >> (BITS_PER_WORD/2);
    a[1] = hi & (UV_MAX >> (BITS_PER_WORD/2));
    a[2] = lo >> (BITS_PER_WORD/2);
    a[3] = lo & (UV_MAX >> (BITS_PER_WORD/2));
    do {
      r = a[0];
      d = r/10;  r = ((r-d*10) << (BITS_PER_WORD/2)) + a[1];  a[0] = d;
      d = r/10;  r = ((r-d*10) << (BITS_PER_WORD/2)) + a[2];  a[1] = d;
      d = r/10;  r = ((r-d*10) << (BITS_PER_WORD/2)) + a[3];  a[2] = d;
      d = r/10;  r = r-d*10;  a[3] = d;
      str[slen++] = '0'+(r%10);
    } while (a[0] || a[1] || a[2] || a[3]);
  }
#endif
  /* Reverse the order */
  for (i=0; i < slen/2; i++) {
    char t=str[i];
    str[i]=str[slen-i-1];
    str[slen-i-1] = t;
  }
  /* Prepend a negative sign if needed */
  if (isneg) {
    for (i = slen; i > 0; i--)
      str[i] = str[i-1];
    str[0] = '-';
    slen++;
  }
  /* Add terminator */
  str[slen] = '\0';
  return slen;
}

/* Oddball primality test.
 * In this file rather than primality.c because it uses factoring (!).
 * Algorithm from Charles R Greathouse IV, 2015 */
static INLINE uint32_t _catalan_v32(uint32_t n, uint32_t p) {
  uint32_t s = 0;
  while (n /= p)  s += n % 2;
  return s;
}
static INLINE uint32_t _catalan_v(UV n, UV p) {
  uint32_t s = 0;
  while (n /= p)  s += n % 2;
  return s;
}
static UV _catalan_mult(UV m, UV p, UV n, UV a) {
  if (p > a) {
    m = mulmod(m, p, n);
  } else {
    UV pow = (n <= 4294967295UL) ? _catalan_v32(a<<1,p) : _catalan_v(a<<1,p);
    m = (pow == 0) ? m
      : (pow == 1) ? mulmod(m,p,n)
                   : mulmod(m,powmod(p,pow,n),n);
  }
  return m;
}
static int _catalan_vtest(UV n, UV p) {
  while (n /= p)
    if (n % 2)
      return 1;
  return 0;
}
int is_catalan_pseudoprime(UV n) {
  UV m, a;
  int i;

  if (n < 2 || ((n % 2) == 0 && n != 2)) return 0;
  if (is_prob_prime(n)) return 1;

  m = 1;
  a = n >> 1;
  /*
   * Ideally we could use some of the requirements for a mod 4/8/64 here:
   * http://www.combinatorics.net/conf/Z60/sp/sp/Shu-Chung%20Liu.pdf
   * But, how do we make +/-2 = X mod n into a solution for x = X mod 8?
   *
   * We could also exploit the exhaustive testing that shows there only
   * exist three below 1e10:  5907, 1194649, and 12327121.
   */
  {
    UV factors[MPU_MAX_FACTORS+1];
    int nfactors = factor_exp(n, factors, 0);
#if BITS_PER_WORD == 32
    if (nfactors == 2) return 0;  /* Page 9, all 32-bit semiprimes */
#else
    if (nfactors == 2) {   /* Conditions from Aebi and Cairns (2008) */
      if (n < UVCONST(10000000000)) return 0;     /* Page 9 */
      if (2*factors[0]+1 >= factors[1]) return 0; /* Corollary 2 and 3 */
    }
#endif
    /* Test every factor */
    for (i = 0; i < nfactors; i++) {
      if (_catalan_vtest(a << 1, factors[i]))
        return 0;
    }
  }
  {
    UV seg_base, seg_low, seg_high;
    unsigned char* segment;
    void* ctx;
    m = _catalan_mult(m, 2, n, a);
    m = _catalan_mult(m, 3, n, a);
    m = _catalan_mult(m, 5, n, a);
    ctx = start_segment_primes(7, n, &segment);
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high ) {
        m = _catalan_mult(m, p, n, a);
      } END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
  }
  return (a & 1) ? (m==(n-1)) : (m==1);
}

/* If we have fast CTZ, use this GCD.  See Brent Alg V and FLINT Abhinav Baid */
UV gcdz(UV x, UV y) {
  UV f, x2, y2;

  if (x == 0) return y;

  if (y & 1) {  /* Optimize y odd */
    x >>= ctz(x);
    while (x != y) {
      if (x < y) { y -= x; y >>= ctz(y); }
      else       { x -= y; x >>= ctz(x); }
    }
    return x;
  }

  if (y == 0) return x;

  /* Alternately:  f = ctz(x|y); x >>= ctz(x); y >>= ctz(y); */
  x2 = ctz(x);
  y2 = ctz(y);
  f = (x2 <= y2) ? x2 : y2;
  x >>= x2;
  y >>= y2;

  while (x != y) {
    if (x < y) {
      y -= x;
      y >>= ctz(y);
    } else {
      x -= y;
      x >>= ctz(x);
    }
  }
  return x << f;
}

/* The intermediate values are so large that we can only stay in 64-bit
 * up to 53 or so using the divisor_sum calculations.  So just use a table.
 * Save space by just storing the 32-bit values. */
static const int32_t tau_table[] = {
  0,1,-24,252,-1472,4830,-6048,-16744,84480,-113643,-115920,534612,-370944,-577738,401856,1217160,987136,-6905934,2727432,10661420,-7109760,-4219488,-12830688,18643272,21288960,-25499225,13865712,-73279080,24647168,128406630,-29211840,-52843168,-196706304,134722224,165742416,-80873520,167282496,-182213314,-255874080,-145589976,408038400,308120442,101267712,-17125708,-786948864,-548895690,-447438528
};
#define NTAU (sizeof(tau_table)/sizeof(tau_table[0]))
IV ramanujan_tau(UV n) {
  return (n < NTAU)  ?  tau_table[n]  :  0;
}

static UV _count_class_div(UV s, UV b2) {
  UV h = 0, i, ndivisors, *divs, lim;

  lim = isqrt(b2);
  if (lim*lim == b2) lim--;
  if (s > lim) return 0;

  if ((lim-s) < 70) {  /* Iterate looking for divisors */
    for (i = s; i <= lim; i++)
      if (b2 % i == 0)
        h++;
  } else {             /* Walk through all the divisors */
    divs = _divisor_list(b2, &ndivisors);
    for (i = 0; i < ndivisors && divs[i] <= lim; i++)
      if (divs[i] >= s)
        h++;
    Safefree(divs);
  }
  return h;
}

/* Returns 12 * H(n).  See Cohen 5.3.5 or Pari/GP.
 * Pari/GP uses a different method for n > 500000, which is quite a bit
 * faster, but assumes the GRH. */
IV hclassno(UV n) {
  UV nmod4 = n % 4, b2, b, h;
  int square;

  if (n == 0) return -1;
  if (nmod4 == 1 || nmod4 == 2) return 0;
  if (n == 3) return 4;

  b = n & 1;
  b2 = (n+1) >> 2;
  square = is_perfect_square(b2);

  h = divisor_sum(b2,0) >> 1;
  if (b == 1)
    h = 1 + square + ((h - 1) << 1);
  b += 2;

  for (;  b2 = (n + b*b) >> 2, 3*b2 < n;  b += 2) {
    h += (b2 % b == 0)
      +  is_perfect_square(b2)
      +  (_count_class_div(b+1, b2) << 1);
  }
  return 12*h + ((b2*3 == n) ? 4 : square && !(n&1) ? 6 : 0);
}

UV polygonal_root(UV n, UV k, int* overflow) {
  UV D, R;
  MPUassert(k >= 3, "is_polygonal root < 3");
  *overflow = 0;
  if (n <= 1) return n;
  if (k == 4) return is_perfect_square(n) ? isqrt(n) : 0;
  if (k == 3) {
    if (n >= UV_MAX/8) *overflow = 1;
    D = n << 3;
    R = 1;
  } else {
    if (k > UV_MAX/k || n > UV_MAX/(8*k-16)) *overflow = 1;
    D = (8*k-16) * n;
    R = (k-4) * (k-4);
  }
  if (D+R <= D) *overflow = 1;
  D += R;
  if (*overflow || !is_perfect_square(D)) return 0;
  D = isqrt(D) + (k-4);
  R = 2*k - 4;
  if ((D % R) != 0) return 0;
  return D/R;
}

/* These rank/unrank are O(n^2) algorithms using O(n) in-place space.
 * Bonet 2008 gives O(n log n) algorithms using a bit more space.
 */

int num_to_perm(UV k, int n, int *vec) {
  int i, j, t, si = 0;
  UV f = factorial(n-1);

  while (f == 0) /* We can handle n! overflow if we have a valid k */
    f = factorial(n - 1 - ++si);

  if (k/f >= (UV)n)
    k %= f*n;

  for (i = 0; i < n; i++)
    vec[i] = i;
  for (i = si; i < n-1; i++) {
    UV p = k/f;
    k -= p*f;
    f /= n-i-1;
    if (p > 0) {
      for (j = i+p, t = vec[j]; j > i; j--)
        vec[j] = vec[j-1];
      vec[i] = t;
    }
  }
  return 1;
}

int perm_to_num(int n, int *vec, UV *rank) {
  int i, j, k;
  UV f, num = 0;
  f = factorial(n-1);
  if (f == 0) return 0;
  for (i = 0; i < n-1; i++) {
    for (j = i+1, k = 0; j < n; j++)
      if (vec[j] < vec[i])
        k++;
    if ((UV)k > (UV_MAX-num)/f) return 0;  /* overflow */
    num += k*f;
    f /= n-i-1;
  }
  *rank = num;
  return 1;
}

/*
 * For k<n, an O(k) time and space method is shown on page 39 of
 *    https://www.math.upenn.edu/~wilf/website/CombinatorialAlgorithms.pdf
 * Note it requires an O(k) complete shuffle as the results are sorted.
 *
 * This seems to be 4-100x faster than NumPy's random.{permutation,choice}
 * for n under 100k or so.  It's even faster with larger n.  For example
 *   from numpy.random import choice;  choice(100000000, 4, replace=False)
 * uses 774MB and takes 55 seconds.  We take less than 1 microsecond.
 */
void randperm(void* ctx, UV n, UV k, UV *S) {
  UV i, j;

  if (k > n)  k = n;

  if        (k == 0) {                  /* 0 of n */
  } else if (k == 1) {                  /* 1 of n.  Pick one at random */
    S[0] = urandomm64(ctx,n);
  } else if (k == 2 && n == 2) {        /* 2 of 2.  Flip a coin */
    S[0] = urandomb(ctx,1);
    S[1] = 1-S[0];
  } else if (k == 2) {                  /* 2 of n.  Pick 2 skipping dup */
    S[0] = urandomm64(ctx,n);
    S[1] = urandomm64(ctx,n-1);
    if (S[1] >= S[0]) S[1]++;
  } else if (k < n/100 && k < 30) {     /* k of n.  Pick k with loop */
    for (i = 0; i < k; i++) {
      do {
        S[i] = urandomm64(ctx,n);
        for (j = 0; j < i; j++)
          if (S[j] == S[i])
            break;
      } while (j < i);
    }
  } else if (k < n/100 && n > 1000000) {/* k of n.  Pick k with dedup retry */
    for (j = 0; j < k; ) {
      for (i = j; i < k; i++) /* Fill S[j .. k-1] then sort S */
        S[i] = urandomm64(ctx,n);
      qsort(S, k, sizeof(UV), _numcmp);
      for (j = 0, i = 1; i < k; i++)  /* Find and remove dups.  O(n). */
        if (S[j] != S[i])
          S[++j] = S[i];
      j++;
    }
    /* S is sorted unique k-selection of 0 to n-1.  Shuffle. */
    for (i = 0; i < k; i++) {
      j = urandomm64(ctx,k-i);
      { UV t = S[i]; S[i] = S[i+j]; S[i+j] = t; }
    }
  } else if (k < n/4) {                 /* k of n.  Pick k with mask */
    uint32_t *mask, smask[8] = {0};
    if (n <= 32*8) mask = smask;
    else           Newz(0, mask, n/32 + ((n%32)?1:0), uint32_t);
    for (i = 0; i < k; i++) {
      do {
        j = urandomm64(ctx,n);
      } while ( mask[j>>5] & (1U << (j&0x1F)) );
      S[i] = j;
      mask[j>>5] |= (1U << (j&0x1F));
    }
    if (mask != smask) Safefree(mask);
  } else if (k < n) {                   /* k of n.  FYK shuffle n, pick k */
    UV *T;
    New(0, T, n, UV);
    for (i = 0; i < n; i++)
      T[i] = i;
    for (i = 0; i < k && i <= n-2; i++) {
      j = urandomm64(ctx,n-i);
      S[i] = T[i+j];
      T[i+j] = T[i];
    }
    Safefree(T);
  } else {                              /* n of n.  FYK shuffle. */
    for (i = 0; i < n; i++)
      S[i] = i;
    for (i = 0; i < k && i <= n-2; i++) {
      j = urandomm64(ctx,n-i);
      { UV t = S[i]; S[i] = S[i+j]; S[i+j] = t; }
    }
  }
}


int is_smooth(UV n, UV k) {
  UV fac[MPU_MAX_FACTORS+1];
  int nfac;

  /* True if no prime factors of n are larger than k. */
  if (n <= 1) return 1;   /* (0,k) = 1, (1,k) = 1 */
  if (k <= 1) return 0;   /* (n,0) = (n,1) = 0 if n > 1 */
  if (n <= k) return 1;
  /* k >= 2, n >= 2 */
  if (k == 2) return ((n & (n-1)) == 0);
  while (n > 1 && !(n&1)) n >>= 1;
  if (n <= k) return 1;
  /* k >= 3, n >= 3 */

  /* Could do: gcd with primorial, p-1/pbrent/ecm to find small factors */

  if (k <= 9000) {
    nfac = trial_factor(n, fac, 2, k);
    return (fac[nfac-1] <= k);
  }

  /* Get rid of small factors */
  nfac = trial_factor(n, fac, 2, 500);
  n = fac[nfac-1];
  if (n <= k) return 1;

  /* Complete factoring including primality test */
  nfac = factor_exp(n, fac, 0);
  return (fac[nfac-1] <= k);
}
int is_rough(UV n, UV k) {
  UV fac[MPU_MAX_FACTORS+1];
  int nfac;

  /* True if no prime factors of n are smaller than k. */

  if (n == 0) return (k == 0);
  if (n == 1) return 1;
  /* n >= 2 */
  if (k <= 1) return 1;
  if (k == 2) return (n >= 1);
  if (k == 3) return (n > 1 && (n&1));
  /* k >= 4 */

  if (!(n&1)) return 0;
  if (!(n%3)) return 0;
  if (k <= 5) return 1;
  if (!(n%5)) return 0;

  if (k <= 2500) {
    nfac = trial_factor(n, fac, 7, k);
    return (fac[0] >= k);
  }

  /* TODO: look into factor_one. */
  /* But it doesn't guarantee returning a prime factor or the smallest. */

  nfac = trial_factor(n, fac, 7, 200);
  if (nfac > 1 && fac[nfac-2] <= k) return 0;
  n = fac[nfac-1];

  if ( (n >> 30) >= 64) {  /* Arbitrarily chose 2^36 for more tests */
    if (is_prime(n)) return 1;
    nfac = pminus1_factor(n, fac, 500, 500);
    if (nfac > 1) {  /* 2 factors, but they could be composites */
      UV f1 = fac[0], f2 = fac[1];
      nfac = factor_exp(f1, fac, 0);
      if (fac[0] < k) return 0;
      nfac = factor_exp(f2, fac, 0);
      if (fac[0] < k) return 0;
      return 1;
    }
  }

  nfac = factor_exp(n, fac, 0);
  return (fac[0] >= k);
}

int is_powerful(UV n, UV k) {
  UV pk;
  int res;

  if (n <= 1 || k <= 1) return 1;

  if (!(n&1)) { /* Check and remove all multiples of 2 */
    if (n & ((UVCONST(1) << k)-1)) return 0;
    n >>= ctz(n);
    if (n == 1) return 1;
  }

  /* Quick checks */
  if (k == 2) {
    if (   (!(n %  3) && (n %     9))
        || (!(n %  5) && (n %    25))
        || (!(n %  7) && (n %    49))
        || (!(n % 11) && (n %   121))
        || (!(n % 13) && (n %   169)) ) return 0;
  } else if (k == 3) {
    if (   (!(n %  3) && (n %    27))
        || (!(n %  5) && (n %   125))
        || (!(n %  7) && (n %   343))
        || (!(n % 11) && (n %  1331)) ) return 0;
  } else {
    if (   (!(n %  3) && (n %    81))
        || (!(n %  5) && (n %   625))
        || (!(n %  7) && (n %  2401))
        || (!(n % 11) && (n % 14641)) ) return 0;
  }

#if 0   /* Full factoring.  Very simple and reasonably efficient. */
  {
    UV fac[MPU_MAX_FACTORS+1];
    UV exp[MPU_MAX_FACTORS+1];
    int i, nfactors;
    nfactors = factor_exp(n, fac, exp);
    for (i = 0; i < nfactors; i++)
      if (exp[i] < k)
        return 0;
    return 1;
  }
#endif

  /* Rather than full factoring, we'll use trial division.  For k=2, we
   * only need to check up to the fourth root of n, and k=3 to the sixth.
   * Even for k=2 this is faster than full factoring on average. */

  /* At every checkpoint (new prime p) for k=2 either:
   *   1) N < p^4  and N=1, p^2, q^2, p^3, or q^3. (q>p).  Return 1.
   *   2) N < p^4  otherwise, N cannot be powerful.  Return 0;
   *   3) N = p^4  is_square caught this and returned 1.
   * So the next possibility is p^2 * q^2 where q = next_prime(p).
   * Check n < p^4 before each new prime, and condition 1 after modifying n.
   */

  /* We could save 5-10% by special casing k=2 and pulling 3 from the loop */

  if (n == 1 || powerof(n) >= k) return 1;
  res = 0;
  START_DO_FOR_EACH_PRIME(3, rootint(n, 2*k)) {
    pk = ipow(p,k);
    if (n < pk*pk) break;
    if (!(n%p)) {
      if (n%pk) break;
      for (n /= pk; (n%p) == 0; n /= p) ;
      if (n == 1 || powerof(n) >= k) { res = 1; break; }
    }
  } END_DO_FOR_EACH_PRIME
  return res;
}

static UV _divsum1(UV prod, UV f, uint32_t e) {
  UV pke, fmult;
  for (pke = f, fmult = 1+f; e > 1; e--) {
    pke *= f;
    fmult += pke;
  }
  return prod * fmult;
}

int is_practical(UV n) {
  UV fac[MPU_MAX_FACTORS+1];
  UV exp[MPU_MAX_FACTORS+1];
  UV prod;
  int i, nfactors;

  if (n == 0 || (n & 1)) return (n == 1);
  if ((n & (n-1)) == 0) return 1;  /* All powers of 2 are practical */
  /* Allowable prefixes: {6,4} => {6,20,28,8} => {6,20,28,88,104,16} */
  if ((n % 6) && (n % 20) && (n % 28) && (n % 88) && (n % 104) && (n % 16))
    return 0;

  nfactors = factor_exp(n, fac, exp);
  /* fac[0] must be 2 */
  prod = 1;  /* running divisor sum of product of all previous factors */
  for (i = 1; i < nfactors; i++) {
    /* prod *= ipow(fac[i-1],exp[i-1]);  sum = 1 + divisor_sum(prod,1); */
    prod = _divsum1(prod, fac[i-1], exp[i-1]);
    if (fac[i] > (1 + prod))
      return 0;
  }
  return 1;
}

static unsigned char* _squarefree_range(UV lo, UV hi) {
  unsigned char* isf;
  UV i, i2, j, range = hi-lo+1, sqrthi = isqrt(hi);

  if (hi < lo) return 0;

  New(0, isf, range, unsigned char);
  memset(isf, 1, range);
  if (lo == 0) isf[0] = 0;
  for (i = 2; i <= sqrthi; i++)
    for (i2 = i*i, j = P_GT_LO(i2, i2, lo); j <= hi; j += i2)
      isf[j] = 0;
  return isf;
}

static UV _pcr(UV n, UV k, unsigned char* isf, UV m, UV r) {
  UV i, sum = 0, lim = rootint(n/m, r);

  if (r <= k) return lim;

  if (r-1 == k) {
    for (i = 1; i <= lim; i++)
      if (isf[i] && gcd_ui(m,i) == 1)
        sum += rootint(n/(m*ipow(i,r)),k);
  } else {
    for (i = 1; i <= lim; i++)
      if (isf[i] && gcd_ui(m,i) == 1)
        sum += _pcr(n, k, isf,  m * ipow(i,r), r-1);
  }
  return sum;
}

UV powerful_count(UV n, UV k) {
  UV i, r, lim, sum = 0;
  unsigned char *isf;

  if (k == 0) return 0;
  if (k == 1 || n <= 1) return n;
  if (k >= BITS_PER_WORD) return 1;

  lim = rootint(n, k+1);
  isf = _squarefree_range(0, lim);

  if (k == 2) {
    for (i = 1; i <= lim; i++)
      if (isf[i])
        sum += isqrt(n/(i*i*i));
  } else {
    /* sum = _pcr(n, k, isf,  1, 2*k-1); */
    r = 2*k-1;
    lim = rootint(n, r);
    for (i = 1; i <= lim; i++)
      if (isf[i])
        sum += _pcr(n, k, isf,  ipow(i,r), r-1);
  }

  Safefree(isf);
  return sum;
}

/* We want:
 *    k=0 turned into k=2 in XS (0 here ok)
 *    n=0 undef in XS (0 here ok)
 *    k=1 => n
 *    n=1 => 1
 *    n=2 => 1<<k
 *    overflow here should return 0
 */
UV nth_powerful(UV n, UV k) {
  static UV const maxpow[11] = {0,UV_MAX,9330124695,11938035,526402,85014,25017,10251,5137,2903,1796};
  static unsigned char const mink[20+1] = {0,0,1,2,4,6,7,9,11,12,14,16,18,19,21,23,24,26,28,30,31};
  UV lo, hi, mid;

  if (k == 0 || k >= BITS_PER_WORD) return 0;
  if (k == 1 || n <= 1) return n;

  if (k <= 10 && n > maxpow[k]) return 0;
  if (k > 10 && (n > maxpow[10] || n > powerful_count(UV_MAX,k))) return 0;

  if (n <= 20 && k >= mink[n]) return UVCONST(1) << (k+(n-2));
  /* Now k >= 2, n >= 4 */

  if (k == 2) { /* From Mincu and Panaitopol 2009 */
    double n53 = pow(n, 5.0/3.0);
    double nc = pow(n, 2) / pow(2.1732543125195541, 2);
    double dlo = nc + 0.3 * n53;
    double dhi = nc + 0.5 * n53;
    lo = (UV) dlo;
    hi = (n < 170) ? 8575 : (dhi >= UV_MAX) ? UV_MAX : 1 + (UV) dhi;
  } else { /* min/max so not good at all */
    lo = (UVCONST(1) << (k+1))+1;
    hi = UV_MAX;
  }

  return inverse_interpolate(lo, hi, n, k, &powerful_count, 0);
}

/* n A069623; 10^n A070428 */
UV perfect_power_count(UV n) {
  uint32_t k, log2n;
  UV sum = 1;

  if (n <= 1) return n;

  log2n = log2floor(n);
  for (k = 2; k <= log2n; k++) {
    IV m = moebius(k);
    if (m != 0) sum -= m * (rootint(n, k)-1);
  }
  return sum;
}

#if 0
/* This is the de Bruijn approximation, not exact! */
static long double dickman_rho(long double u) {
  int i;
  long double zeta;

  if (u <= 1) return 1;
  if (u <= 2) return 1-logl(u);

  /* Also see Granville (2008) Smooth numbers */

  /* Calculate zeta.  See Bach and Sorenson (2013) page 10 */
  zeta = 2*(u-1);
  for (i = 0; i < 7; i++) {
    long double uz1 = 1 + u*zeta;
    zeta = zeta - ( (zeta-logl(uz1))*uz1 ) / (uz1-u);
  }
  /* Alternately:  zeta = -1/u - LambertW1(-exp(-1/u)/u) */

  return expl(-u*zeta+Ei(zeta)) / (zeta * sqrtl(2*3.1415926535*u));
}
#endif

static const unsigned char _psi_cache_v__7[128] = {8,9,10,10,11,11,12,13,14,14,15,15,16,17,17,17,18,19,19,20,21,21,22,22,23,23,23,24,25,25,25,25,26,26,27,27,27,28,28,28,29,30,31,31,31,31,32,32,33,33,33,33,34,34,34,35,36,36,36,36,36,36,37,37,38,38,38,39,39,39,39,39,40,41,41,41,42,42,42,42,42,42,43,43,43,43,43,43,44,44,45,45,46,46,46,46,46,47,47,47,48,48,48,48,49,49,49,49,49,49,49,49,50,50,50,50,50,51,52,52,53,53,53,53,53,53,53,54};
static const unsigned char _psi_cache_v_11[96] = {12,12,13,14,15,15,16,16,17,18,19,19,20,21,21,22,23,23,24,24,25,26,26,27,28,28,28,28,29,29,30,30,31,32,32,32,33,34,35,35,35,35,36,37,38,38,38,38,39,39,39,40,41,41,42,42,42,42,43,43,44,44,44,45,45,46,46,46,47,48,48,48,49,49,49,49,50,50,51,51,51,51,51,51,52,52,53,54,55,55,55,55,55,56,56,56};
static const unsigned char _psi_cache_v_13[64] = {14,15,16,16,17,17,18,19,20,20,21,22,23,24,25,25,26,26,27,28,28,29,30,30,30,31,32,32,33,33,34,35,35,35,36,37,38,38,39,39,40,41,42,42,42,42,43,43,43,44,45,46,47,47,47,47,48,48,49,49,49,50,50,51};

UV debruijn_psi(UV x, UV y) {
  UV sum, x3, x5;
  if (x < 1) return 0;
  if (y <= 1) return 1;
  if (y >= x) return x;
  if (y == 2) return 1 + log2floor(x);
  if (!(y&1)) y--; /* Make y odd for simplicity */

  /* Caches etc. to avoid recursion - about 1.6x speedup for big inputs */
  if (y ==  7 && x- 7 <=128) return _psi_cache_v__7[x-1- 7];
  if (y == 11 && x-11 <= 96) return _psi_cache_v_11[x-1-11];
  if (y == 13 && x-13 <= 64) return _psi_cache_v_13[x-1-13];
  if (y >= 17 && x <= 128) {
    /* mpu 'for (7..128) { $f=(factor($_))[-1]; push(@$X,$_),push(@$Y,$f) if $f > 17; }  say scalar(@$X); say join(",",@$_) for ($X,$Y);' */
    static const unsigned char xt[48] = {19,23,29,31,37,38,41,43,46,47,53,57,58,59,61,62,67,69,71,73,74,76,79,82,83,86,87,89,92,93,94,95,97,101,103,106,107,109,111,113,114,115,116,118,122,123,124,127};
    static const unsigned char yt[48] = {19,23,29,31,37,19,41,43,23,47,53,19,29,59,61,31,67,23,71,73,37,19,79,41,83,43,29,89,23,31,47,19,97,101,103,53,107,109,37,113,19,23,29,59,61,41,31,127};
    unsigned char i;
    for (i = 0, sum = x; i < 48 && x >= xt[i]; i++)
      if (y < yt[i])
        sum--;
    return sum;
  }

  /*  given z < y < x,  (e.g. z=2 or z=19)
   *  psi(x,y) = psi(x,z) + sum[z+1..y] psi(x/p,p) */

  sum = 1 + log2floor(x);  /* debruijn_psi(x,2) */
  /* if (y >= 3)  sum += debruijn_psi(x/3, 3); */
  /* if (y >= 5)  sum += debruijn_psi(x/5, 5); */
  if (y >= 3) {
    for (x3 = x/3; x3 > 3; x3 /= 3)
      sum += 1+log2floor(x3);
    sum += x3;
  }
  if (y >= 5) {
    for (x5 = x/5; x5 > 5; x5 /= 5) {
      sum += 1+log2floor(x5);
      for (x3 = x5/3; x3 > 3; x3 /= 3)
        sum += 1+log2floor(x3);
      sum += x3;
    }
    sum += x5;
  }
  if (y >=  7) sum += debruijn_psi(x/ 7, 7);
  if (y >= 11) sum += debruijn_psi(x/11,11);
  if (y >= 13) sum += debruijn_psi(x/13,13);
  if (y >= 17) sum += debruijn_psi(x/17,17);
  if (y >= 19) sum += debruijn_psi(x/19,19);
  if (y >= 23) sum += debruijn_psi(x/23,23);
  if (y >= 29) {
    START_DO_FOR_EACH_PRIME(29, y) {
      UV xp = x/p;
      sum += (p >= xp)  ?  xp  :  debruijn_psi(xp, p);
    } END_DO_FOR_EACH_PRIME
  }

  return sum;
}

UV buchstab_phi(UV x, UV y) {
  if (y <= 2) return x;
  if (y <= 3) return x - x/2;
  if (y <= 5) return x - x/2 - x/3 + x/6;
  /* We'll use the legendre_phi function we already have. */
  return legendre_phi(x, LMO_prime_count(y-1));
}


UV random_factored_integer(void* ctx, UV n, int *nf, UV *factors) {
  UV r, s, nfac;
  if (n < 1)
    return 0;
#if BITS_PER_WORD == 64 && (USE_MONTMATH || MULMODS_ARE_FAST)
  if (1)   /* Our factoring is very fast, just use it */
#elif BITS_PER_WORD == 64
  if (n < UVCONST(1000000000000))
#endif
  {
    r = 1 + urandomm64(ctx, n);
    *nf = factor(r, factors);
    return r;
  }
  do {  /* Kalai's algorithm */
    for (s = n, r = 1, nfac = 0;  s > 1;  ) {
      s = 1 + urandomm64(ctx, s);
      if (!is_prime(s)) continue;
      if (s > n / r) { r = 0; break; } /* overflow */
      r *= s;
      factors[nfac++] = s;
    }
  } while (r == 0 || r > n || (1 + urandomm64(ctx,n)) > r);
  *nf = nfac;
  return r;
}

/******************************************************************************/
/*                             LUCKY NUMBERS                                  */
/******************************************************************************/

static const char _lmask63[63+2] = {1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0,0,1,1,0,1,1,0,0,0,0,1,1,0,1,1,0,1,1,0,0,0,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0,0,0,0,0,1,1};
static const unsigned char _small_lucky[48] = {1,3,7,9,13,15,21,25,31,33,37,43,49,51,63,67,69,73,75,79,87,93,99,105,111,115,127,129,133,135,141,151,159,163,169,171,189,193,195,201,205,211,219,223,231,235,237,241};

/* Lucky Number sieve for 32-bit inputs.
 * Pre-sieve for first 4-7 levels, then in-place deletion using memmove,
 * plus an optimization for a single pass for all single skips.
 * On x86 at least, faster than the other sieves, but uses more memory.
 */
uint32_t* lucky_sieve32(UV *size, UV n) {
  UV i, m, l15, l21, lsize = 0, level, init_level;
  uint32_t *lucky;

  if (n == 0) { *size = 0; return 0; }
  MPUassert(n < 0xFFFFFFFFUL, "lucky_sieve32 n too large");

  /* make initial list, with more culling if it seems worth doing */
  if (n < 300000) {
    const UV fsize = 96*(n+377)/378;  /* 378 = 2*3*7*9 */
    New(0, lucky, 1 + fsize, uint32_t);
    /* Cut 2,3 using mod 6 wheel, 7,9 using a mod 63 mask */
    for (i = 1, m = 1; i <= n; i += 6) {
      if (_lmask63[m  ])  lucky[lsize++] = i;
      if (_lmask63[m+2])  lucky[lsize++] = i+2;
      if ((m += 6) >= 63) m -= 63;
    }
    init_level = 4;
  } else {
    char mask819[819+2];
    const uint16_t v13[16] = {45,96,147,198,252,300,351,402,453,507,555,609,660,708,762,810};
    const UV fsize = (n < 1000000) ? 16128*(n+73709)/73710
                                   : 322560*(n+1547909)/1547910;
    New(0, lucky, 1 + fsize, uint32_t);
    /* Create the mod 819 mask from the smaller one */
    for (i = 0; i < 13; i++) memcpy(mask819+63*i,_lmask63,65);
    for (i = 0; i < 16; i++) mask819[v13[i]] = mask819[v13[i]+1] = 0;
    /* Use the mask and additionally two counters for another two levels */
    for (i = 1, m = 1, l15 = 0, l21 = 0; i <= n; i += 6) {
      if (mask819[m  ] && ++l15 != 15 && ++l21 != 21)  lucky[lsize++] = i;
      if (mask819[m+2] && ++l15 != 15 && ++l21 != 21)  lucky[lsize++] = i+2;
      if ((m += 6) >= 819) m -= 819;
      if (l15 >= 15) l15 -= 15;
      if (l21 >= 21) l21 -= 21;
    }
    init_level = 7;
  }
  if (lucky[lsize-1] > n) lsize--;   /* The +2 could have gone past N */

  /* After the fill-in, we'll start deleting at 13 or 25 */
  for (level = init_level; level < lsize && lucky[level]-1 < lsize; level++) {
    UV skip = lucky[level]-1, nlsize = skip;
    if (2*(skip+1) > lsize) break;  /* Only single skips left */
    for (i = skip+1; i < lsize; i += skip+1) {
      UV ncopy = (skip <= (lsize-i)) ? skip : (lsize-i);
      memmove( lucky + nlsize, lucky + i, ncopy * sizeof(uint32_t) );
      nlsize += ncopy;
    }
    lsize = nlsize;
  }
  /* Now we just have single skips.  Process them all in one pass. */
  if (level < lsize && lucky[level]-1 < lsize) {
    UV skip = lucky[level], nlsize = skip-1;
    while (skip < lsize) {
      UV ncopy = lucky[level+1] - lucky[level];
      if (ncopy > lsize-skip)  ncopy = lsize - skip;
      memmove(lucky + nlsize, lucky + skip, ncopy * sizeof(uint32_t));
      nlsize += ncopy;
      skip += ncopy + 1;
      level++;
    }
    lsize = nlsize;
  }
  *size = lsize;
  return lucky;
}
/* Lucky Number sieve for 64-bit inputs.
 * Uses running counters to skip entries while we add them.
 * Based substantially on Hugo van der Sanden's cgen_lucky.c.
 */
UV* lucky_sieve(UV *size, UV n) {
  UV i, j, c3, lsize, lmax, lindex, *lucky, *count;

  if (n == 0) { *size = 0; return 0; }

  /* Init */
  lmax = (n < 1000) ? 153 : 100 + n/log(n);
  New(0, lucky, lmax, UV);
  New(0, count, lmax, UV);
  lucky[0] = 1;
  lucky[1] = 3;
  lucky[2] = 7;
  lindex = 2;
  lsize = 1;
  c3 = 2;

  for (i = 3; i <= n; i += 2) {
    if (!--c3) { c3 = 3; continue; }  /* Shortcut count[1] */
    for (j = 2; j < lindex; j++) {
      if (--count[j] == 0) {
        count[j] = lucky[j];
        break;
      }
    }
    if (j < lindex) continue;

    if (lsize >= lmax) {  /* Given the estimate, we probably never do this. */
      lmax = 1 + lsize * 1.2;
      Renew(lucky, lmax, UV);
      Renew(count, lmax, UV);
    }
    lucky[lsize] = count[lsize] = i;
    lsize++;

    if (lucky[lindex] == lsize) {
      lindex++;  lsize--;  /* Discard immediately */
    }
  }
  Safefree(count);
  *size = lsize;
  return lucky;
}

UV nth_lucky(UV n) {
  UV i, k, nlucky;

  if (n <= 48)  return (n == 0) ? 0 : _small_lucky[n-1];

  /* Apply the backward sieve, ala Wilson, for entry n */
  if (n <= UVCONST(2000000000)) {
    uint32_t *lucky32 = lucky_sieve32(&nlucky, n);
    for (i = nlucky-1, k = n-1; i >= 1; i--)
      k += k/(lucky32[i]-1);
    Safefree(lucky32);
  } else {
    UV *lucky64 = lucky_sieve(&nlucky, n);
    for (i = nlucky-1, k = n-1; i >= 1; i--)
      k += k/(lucky64[i]-1);
    Safefree(lucky64);
  }
  return (2 * k + 1);
}

#if 0 /* Simple versions */
static UV lucky_count_approx(UV n) { return 0.5 + 0.970 * n / log(n); }
static UV lucky_count_upper(UV n) { return 200 + lucky_count_approx(n) * 1.025; }
#else /* Hacky but better */
static UV lucky_count_approx(UV n) {
  return   (n <        7)  ?  (n > 0) + (n > 2)
         : (n <= 1000000)  ?  0.9957 * n/log(n)
                           : (1.03670 - log(n)/299) * n/log(n);
}
static UV lucky_count_upper(UV n) {   /* Holds under 1e9 */
  return (n <= 10000) ?  10 + lucky_count_approx(n) * 1.1
                      : 140 + lucky_count_approx(n) * 1.004;
}
#endif

int is_lucky(UV n) {
  uint32_t *lucky32;
  UV i, l, quo, pos, nlucky, lsize;

  /* Simple pre-tests */
  if ( !(n & 1) || (n%6) == 5 || !_lmask63[n % 63]) return 0;
  if (n < 45) return 1;

  /* Check valid position using the static list */
  pos = (n+1) >> 1;  /* Initial position in odds */
  for (i = 1; i < 48; i++) {
    l = _small_lucky[i];
    if (pos < l) return 1;
    quo = pos / l;
    if (pos == quo*l) return 0;
    pos -= quo;
  }

  /* Check more small values */
  if (n >= 1000000U) {
    lucky32 = lucky_sieve32(&nlucky, lsize = lucky_count_upper(n)/25);
    while (i < nlucky) {
      l = lucky32[i++];
      if (pos < l) break;
      quo = pos / l;
      if (pos == quo*l) { Safefree(lucky32); return 0; }
      pos -= quo;
    }
    Safefree(lucky32);
    if (pos < l) return 1;
  }

  /* Generate all needed values and continue checking from where we left off. */
  lucky32 = lucky_sieve32(&nlucky, lsize = lucky_count_upper(n));
  while (1) {
    if (i >= nlucky) { Safefree(lucky32); lucky32 = lucky_sieve32(&nlucky, lsize *= 1.02); }
    l = lucky32[i++];
    if (pos < l)  break;
    quo = pos / l;
    if (pos == quo*l) { Safefree(lucky32); return 0; }
    pos -= quo;
  }
  Safefree(lucky32);
  return 1;
}
