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
#if _MSC_VER || defined(__IBMC__) | defined(__IBMCPP__) || (defined(__STDC_VERSION__) && __STDC_VERSION >= 199901L)
  /* math.h should give us these as functions or macros.
   *
   *  extern long double fabsl(long double);
   *  extern long double floorl(long double);
   *  extern long double ceill(long double);
   *  extern long double sqrtl(long double);
   *  extern long double powl(long double, long double);
   *  extern long double expl(long double);
   *  extern long double logl(long double);
   */
#else
  #define fabsl(x)    (long double) fabs( (double) (x) )
  #define floorl(x)   (long double) floor( (double) (x) )
  #define ceill(x)    (long double) ceil( (double) (x) )
  #define sqrtl(x)    (long double) sqrt( (double) (x) )
  #define powl(x, y)  (long double) pow( (double) (x), (double) (y) )
  #define expl(x)     (long double) exp( (double) (x) )
  #define logl(x)     (long double) log( (double) (x) )
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

#define KAHAN_INIT(s) \
  long double s ## _y, s ## _t; \
  long double s ## _c = 0.0; \
  long double s = 0.0;

#define KAHAN_SUM(s, term) \
  do { \
    s ## _y = (term) - s ## _c; \
    s ## _t = s + s ## _y; \
    s ## _c = (s ## _t - s) - s ## _y; \
    s = s ## _t; \
  } while (0)


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
#include "prime_nth_count.h"
#include "factor.h"
#include "mulmod.h"
#include "constants.h"
#include "montmath.h"
#include "csprng.h"

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


/* Return a char array with lo-hi+1 elements. mu[k-lo] = µ(k) for k = lo .. hi.
 * It is the callers responsibility to call Safefree on the result. */
#define PGTLO(p,lo)  ((p) >= lo) ? (p) : ((p)*(lo/(p)) + ((lo%(p))?(p):0))
#define P2GTLO(pinit, p, lo) \
   ((pinit) >= lo) ? (pinit) : ((p)*(lo/(p)) + ((lo%(p))?(p):0))
signed char* _moebius_range(UV lo, UV hi)
{
  signed char* mu;
  UV i;
  UV sqrtn = isqrt(hi);

  /* Kuznetsov indicates that the Deléglise & Rivat (1996) method can be
   * modified to work on logs, which allows us to operate with no
   * intermediate memory at all.  Same time as the D&R method, less memory. */
  unsigned char logp;
  UV nextlog;

  Newz(0, mu, hi-lo+1, signed char);
  if (sqrtn*sqrtn != hi) sqrtn++;  /* ceil sqrtn */

  logp = 1; nextlog = 3; /* 2+1 */
  START_DO_FOR_EACH_PRIME(2, sqrtn) {
    UV p2 = p*p;
    if (p > nextlog) {
      logp += 2;   /* logp is 1 | ceil(log(p)/log(2)) */
      nextlog = ((nextlog-1)*4)+1;
    }
    for (i = PGTLO(p, lo); i <= hi; i += p)
      mu[i-lo] += logp;
    for (i = PGTLO(p2, lo); i <= hi; i += p2)
      mu[i-lo] = 0x80;
  } END_DO_FOR_EACH_PRIME

  logp = log2floor(lo);
  nextlog = UVCONST(2) << logp;
  for (i = lo; i <= hi; i++) {
    unsigned char a = mu[i-lo];
    if (i >= nextlog) {  logp++;  nextlog *= 2;  } /* logp is log(p)/log(2) */
    if (a & 0x80)       { a = 0; }
    else if (a >= logp) { a =  1 - 2*(a&1); }
    else                { a = -1 + 2*(a&1); }
    mu[i-lo] = a;
  }
  if (lo == 0)  mu[0] = 0;

  return mu;
}

UV* _totient_range(UV lo, UV hi) {
  UV* totients;
  UV i, seg_base, seg_low, seg_high;
  unsigned char* segment;
  void* ctx;

  if (hi < lo) croak("_totient_range error hi %"UVuf" < lo %"UVuf"\n", hi, lo);
  New(0, totients, hi-lo+1, UV);

  /* Do via factoring if very small or if we have a small range */
  if (hi < 100 || (hi-lo) < 10 || hi/(hi-lo+1) > 1000) {
    for (i = lo; i <= hi; i++)
      totients[i-lo] = totient(i);
    return totients;
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
    memset(totients, 0, (hi-lo+1) * sizeof(UV));
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

  for (i = lo; i <= hi; i++) {
    UV v = i;
    if (i % 2 == 0)  v -= v/2;
    if (i % 3 == 0)  v -= v/3;
    if (i % 5 == 0)  v -= v/5;
    totients[i-lo] = v;
  }

  ctx = start_segment_primes(7, hi/2, &segment);
  while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
    START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high ) {
      for (i = P2GTLO(2*p,p,lo); i <= hi; i += p)
        totients[i-lo] -= totients[i-lo]/p;
    } END_DO_FOR_EACH_SIEVE_PRIME
  }
  end_segment_primes(ctx);

  /* Fill in all primes */
  for (i = lo | 1; i <= hi; i += 2)
    if (totients[i-lo] == i)
      totients[i-lo]--;
  if (lo <= 1) totients[1-lo] = 1;

  return totients;
}

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
  mu = _moebius_range(0, maxmu);
  New(0, M, maxmu+1, short);
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

/* There are at least 4 ways to do this, plus hybrids.
 * 1) use a table.  Great for 32-bit, too big for 64-bit.
 * 2) Use pow() to check.  Relatively slow and FP is always dangerous.
 * 3) factor or trial factor.  Slow for 64-bit.
 * 4) Dietzfelbinger algorithm 2.3.5.  Quite slow.
 * This currently uses a hybrid of 1 and 2.
 */
int powerof(UV n) {
  UV t;
  if ((n <= 3) || (n == UV_MAX)) return 1;
  if ((n & (n-1)) == 0)          return ctz(n);  /* powers of 2    */
  if (is_perfect_square(n))      return 2 * powerof(isqrt(n));
  if (is_perfect_cube(n))        return 3 * powerof(icbrt(n));

  /* Simple rejection filter for non-powers of 5-37.  Rejects 47.85%. */
  t = n & 511; if ((t*77855451) & (t*4598053) & 862)  return 1;

  if (is_perfect_fifth(n))       return 5 * powerof(rootof(n,5));
  if (is_perfect_seventh(n))     return 7 * powerof(rootof(n,7));

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
      UV root = rootof(n,11);
      if (n == ipow(root,11)) return 11;
    }
    if ( (t = n %131, !((t*1545928325) & (t*1355660813) & 2771533888U)) &&
         (t = n % 79, !((t*48902028) & (t*48589927) & 404082779)) ) {
      /* (t = n % 53, !((t*79918293) & (t*236846524) & 694943819)) ) { */
      UV root = rootof(n,13);
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
      return !is_perfect_fifth(n) ? 0 : (a == 5) ? 1 :is_power(rootof(n,5),a/5);
  }
  ret = powerof(n);
  if (a != 0) return !(ret % a);  /* Is the max power divisible by a? */
  return (ret == 1) ? 0 : ret;
}

#if BITS_PER_WORD == 64
#define ROOT_MAX_3 41
static const uint32_t root_max[ROOT_MAX_3] = {0,0,0,2642245,65535,7131,1625,565,255,138,84,56,40,30,23,19,15,13,11,10,9,8,7,6,6,5,5,5,4,4,4,4,3,3,3,3,3,3,3,3,3};
#else
#define ROOT_MAX_3 21
static const uint32_t root_max[ROOT_MAX_3] = {0,0,0,1625,255,84,40,23,15,11,9,7,6,5,4,4,3,3,3,3,3};
#endif

UV rootof(UV n, UV k) {
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
    UV root = rootof(n, (UV)power);
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
    p = nth_prime(sp);
    if (p > n3)
      break;
    if ((n % p) == 0)
      return !!is_prob_prime(n/p);
  }
  /* 9.8% of random inputs left */
  if (is_prob_prime(n)) return 0;
  if (p > n3) return 1;
  /* 4-8% of random inputs left */
  /* n is a composite and larger than p^3 */
  if (   pbrent_factor(n, factors, 70000, 1) == 2
         /* The only things we normally see by now are 16+ digit semiprimes */
      || pminus1_factor(n, factors, 4000, 4000) == 2
         /* 0.09% of random 64-bit inputs left */
      || pbrent_factor(n, factors, 180000, 7) == 2 )
    return (is_prob_prime(factors[0]) && is_prob_prime(factors[1]));
  /* 0.002% of random 64-bit inputs left */
  {
    UV facs[MPU_MAX_FACTORS+1];
    return (factor(n,facs) == 2);
  }
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

UV pillai_v(UV n) {
  UV v, fac = 5040 % n;
  if (n == 0) return 0;
  for (v = 8; v < n-1 && fac != 0; v++) {
    fac = (n < HALF_WORD) ? (fac*v) % n : mulmod(fac,v,n);
    if (fac == n-1 && (n % v) != 1)
      return v;
  }
  return 0;
}


int moebius(UV n) {
  UV factors[MPU_MAX_FACTORS+1];
  UV i, nfactors;

  if (n <= 1) return (int)n;
  if ( n >= 49 && (!(n% 4) || !(n% 9) || !(n%25) || !(n%49)) )
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
  UV exp[MPU_MAX_FACTORS+1];
  UV a, phi, on, r;
  int i, nfactors;

  if (n <= 4) return (n == 0) ? 0 : n-1;
  if (n % 4 == 0)  return 0;

  on = (n&1) ? n : (n>>1);
  a = powerof(on);
  r = rootof(on, a);
  if (!is_prob_prime(r)) return 0;        /* c^a or 2c^a */
  phi = (r-1) * (on/r);                   /* p^a or 2p^a */

  nfactors = factor_exp(phi, fac, exp);
  for (i = 0; i < nfactors; i++)
    exp[i] = phi / fac[i];  /* exp[i] = phi(n) / i-th-factor-of-phi(n) */
  for (a = 2; a < n; a++) {
    /* Skip first few perfect powers */
    if (a == 4 || a == 8 || a == 9) continue;
    /* Skip values we know can't be right: (a|n) = 0 (or 1 for odd primes) */
    if (phi == n-1) {
      if (kronecker_uu(a, n) != -1)  continue;
    } else {
      if (kronecker_uu(a, n) == 0)  continue;
    }
    for (i = 0; i < nfactors; i++)
      if (powmod(a, exp[i], n) == 1)
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
  if (gcd_ui(a,n) != 1) return 0;
  s = nprime ? n-1 : totient(n);

  /* a^x can be a primitive root only if gcd(x,s) = 1 */
  i = is_power(a,0);
  if (i > 1 && gcd_ui(i, s) != 1) return 0;

  /* Quick check for small factors before full factor */
  if ((s % 2) == 0 && powmod(a, s/2, n) == 1) return 0;

#if USE_MONTMATH
  if (n & 1) {
    const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);
    a = mont_geta(a, n);
    if ((s % 3) == 0 && mont_powmod(a, s/3, n) == mont1) return 0;
    if ((s % 5) == 0 && mont_powmod(a, s/5, n) == mont1) return 0;
    nfacs = factor_exp(s, fac, 0);
    for (i = 0; i < nfacs; i++) {
      if (fac[i] > 5 && mont_powmod(a, s/fac[i], n) == mont1) return 0;
    }
  } else
#endif
  {
    if ((s % 3) == 0 && powmod(a, s/3, n) == 1) return 0;
    if ((s % 5) == 0 && powmod(a, s/5, n) == 1) return 0;
    /* Complete factor and check each one not found above. */
    nfacs = factor_exp(s, fac, 0);
    for (i = 0; i < nfacs; i++) {
      if (fac[i] > 5 && powmod(a, s/fac[i], n) == 1) return 0;
    }
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
      for (k = 1; k < exp[j]; k++)
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
  *rs = s;
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

  /* If gcd(a,n)==1, this answers comclusively if a solution exists. */
  if (gcdan == 1) {
    for (i = 0; i < nfactors; i++)
      if (fac[i] > 7 && kronecker_uu(a, fac[i]) != 1) return 0;
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
      UV xk2, yk, expect, sol;
      xk2 = addmod(sqr[i],sqr[i],p);
      yk = modinverse(xk2, p);
      expect = mulmod(xk2,yk,p);
      p *= fac[i];
      sol = submod(sqr[i], mulmod(submod(sqrmod(sqr[i],p), a % p, p), yk, p), p);
      if (expect != 1 || sqrmod(sol,p) != (a % p)) {
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

long double chebyshev_function(UV n, int which)
{
  long double logp, logn = logl(n);
  UV sqrtn = which ? isqrt(n) : 0;  /* for theta, p <= sqrtn always false */
  KAHAN_INIT(sum);

  if (n < 500) {
    UV p, pi;
    for (pi = 1;  (p = nth_prime(pi)) <= n; pi++) {
      logp = logl(p);
      if (p <= sqrtn) logp *= floorl(logn/logp+1e-15);
      KAHAN_SUM(sum, logp);
    }
  } else {
    UV seg_base, seg_low, seg_high;
    unsigned char* segment;
    void* ctx;
    long double logl2 = logl(2);
    long double logl3 = logl(3);
    long double logl5 = logl(5);
    if (!which) {
      KAHAN_SUM(sum,logl2); KAHAN_SUM(sum,logl3); KAHAN_SUM(sum,logl5);
    } else {
      KAHAN_SUM(sum, logl2 * floorl(logn/logl2 + 1e-15));
      KAHAN_SUM(sum, logl3 * floorl(logn/logl3 + 1e-15));
      KAHAN_SUM(sum, logl5 * floorl(logn/logl5 + 1e-15));
    }
    ctx = start_segment_primes(7, n, &segment);
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high ) {
        logp = logl(p);
        if (p <= sqrtn) logp *= floorl(logn/logp+1e-15);
        KAHAN_SUM(sum, logp);
      } END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
  }
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
 *
 * Any mistakes here are completely my fault.  This code has not been
 * verified for anything serious.  For better results, see:
 *    http://www.trnicely.net/pi/pix_0000.htm
 * which although the author claims are demonstration programs, will
 * undoubtedly produce more reliable results than this code does (I don't
 * know of any obvious issues with this code, but it just hasn't been used
 * by many people).
 */

static long double const euler_mascheroni = 0.57721566490153286060651209008240243104215933593992L;
static long double const li2 = 1.045163780117492784844588889194613136522615578151L;

long double Ei(long double x) {
  long double val, term;
  unsigned int n;
  KAHAN_INIT(sum);

  if (x == 0) croak("Invalid input to ExponentialIntegral:  x must be != 0");
  /* Protect against messed up rounding modes */
  if (x >=  12000) return INFINITY;
  if (x <= -12000) return 0;

  if (x < -1) {
    /* Continued fraction, good for x < -1 */
    long double lc = 0;
    long double ld = 1.0L / (1.0L - (long double)x);
    val = ld * (-expl(x));
    for (n = 1; n <= 100000; n++) {
      long double old, t, n2;
      t = (long double)(2*n + 1) - (long double) x;
      n2 = n * n;
      lc = 1.0L / (t - n2 * lc);
      ld = 1.0L / (t - n2 * ld);
      old = val;
      val *= ld/lc;
      if ( fabsl(val-old) <= LDBL_EPSILON*fabsl(val) )
        break;
    }
  } else if (x < 0) {
    /* Rational Chebyshev approximation (Cody, Thacher), good for -1 < x < 0 */
    static const long double C6p[7] = { -148151.02102575750838086L,
                                    150260.59476436982420737L,
                                     89904.972007457256553251L,
                                     15924.175980637303639884L,
                                      2150.0672908092918123209L,
                                       116.69552669734461083368L,
                                         5.0196785185439843791020L };
    static const long double C6q[7] = {  256664.93484897117319268L,
                                    184340.70063353677359298L,
                                     52440.529172056355429883L,
                                      8125.8035174768735759866L,
                                       750.43163907103936624165L,
                                        40.205465640027706061433L,
                                         1.0000000000000000000000L };
    long double sumn = C6p[0]-x*(C6p[1]-x*(C6p[2]-x*(C6p[3]-x*(C6p[4]-x*(C6p[5]-x*C6p[6])))));
    long double sumd = C6q[0]-x*(C6q[1]-x*(C6q[2]-x*(C6q[3]-x*(C6q[4]-x*(C6q[5]-x*C6q[6])))));
    val = logl(-x) - sumn/sumd;
  } else if (x < (-2 * logl(LDBL_EPSILON))) {
    /* Convergent series.  Accurate but slow especially with large x. */
    long double fact_n = x;
    for (n = 2; n <= 200; n++) {
      long double invn = 1.0L / n;
      fact_n *= (long double)x * invn;
      term = fact_n * invn;
      KAHAN_SUM(sum, term);
      /* printf("C  after adding %.20Lf, val = %.20Lf\n", term, sum); */
      if (term < LDBL_EPSILON*sum) break;
    }
    KAHAN_SUM(sum, euler_mascheroni);
    KAHAN_SUM(sum, logl(x));
    KAHAN_SUM(sum, x);
    val = sum;
  } else if (x >= 24) {
    /* Cody / Thacher rational Chebyshev */
    static const long double P2[10] = {
        1.75338801265465972390E02L,-2.23127670777632409550E02L,
        -1.81949664929868906455E01L,-2.79798528624305389340E01L,
        -7.63147701620253630855E00L,-1.52856623636929636839E01L,
        -7.06810977895029358836E00L,-5.00006640413131002475E00L,
        -3.00000000320981265753E00L, 1.00000000000000485503E00L };
    static const long double Q2[9] = {
        3.97845977167414720840E04L, 3.97277109100414518365E00L,
        1.37790390235747998793E02L, 1.17179220502086455287E02L,
        7.04831847180424675988E01L,-1.20187763547154743238E01L,
        -7.99243595776339741065E00L,-2.99999894040324959612E00L,
        1.99999999999048104167E00L };
    long double invx = 1.0L / x;
    long double frac = 0.0;
    for (n = 0; n <= 8; n++)
      frac = Q2[n] / (P2[n] + x + frac);
    frac += P2[9];
    val = expl(x) * (invx + invx*invx*frac);
  } else {
    /* Asymptotic divergent series */
    long double invx = 1.0L / x;
    term = 1.0;
    for (n = 1; n <= 200; n++) {
      long double last_term = term;
      term = term * ( (long double)n * invx );
      if (term < LDBL_EPSILON*sum) break;
      if (term < last_term) {
        KAHAN_SUM(sum, term);
        /* printf("A  after adding %.20llf, sum = %.20llf\n", term, sum); */
      } else {
        KAHAN_SUM(sum, (-last_term/3) );
        /* printf("A  after adding %.20llf, sum = %.20llf\n", -last_term/3, sum); */
        break;
      }
    }
    KAHAN_SUM(sum, 1.0L);
    val = expl(x) * sum * invx;
  }

  return val;
}

long double Li(long double x) {
  if (x == 0) return 0;
  if (x == 1) return -INFINITY;
  if (x == 2) return li2;
  if (x < 0) croak("Invalid input to LogarithmicIntegral:  x must be >= 0");
  if (x >= LDBL_MAX) return INFINITY;

  /* Calculate directly using Ramanujan's series. */
  if (x > 1) {
    const long double logx = logl(x);
    long double sum = 0, inner_sum = 0, old_sum, factorial = 1, power2 = 1;
    long double q, p = -1;
    int k = 0, n = 0;

    for (n = 1, k = 0; n < 200; n++) {
      factorial *= n;
      p *= -logx;
      q = factorial * power2;
      power2 *= 2;
      for (; k <= (n - 1) / 2; k++)
        inner_sum += 1.0L / (2 * k + 1);
      old_sum = sum;
      sum += (p / q) * inner_sum;
      if (fabsl(sum - old_sum) <= LDBL_EPSILON) break;
    }
    return euler_mascheroni + logl(logx) + sqrtl(x) * sum;
  }

  return Ei(logl(x));
}

UV inverse_li(UV x) {
  UV r;
  int i;
  long double t, lx = (long double)x, term, old_term = 0;
  if (x <= 2) return x + (x > 0);
  /* Iterate Halley's method until error grows. */
  for (i = 0, t = lx*logl(x); i < 4; i++) {
    long double dn = Li(t) - lx;
    term = dn*logl(t) / (1.0L + dn/(2*t));
    if (i > 0 && fabsl(term) >= fabsl(old_term)) { t -= term/4; break; }
    old_term = term;
    t -= term;
  }
  r = (UV)ceill(t);

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

UV inverse_R(UV x) {
  int i;
  long double t, dn, lx = (long double) x, term, old_term = 0;
  if (x <= 2) return x + (x > 0);

  /* Rough estimate */
  t = lx * logl(x);
  /* Improve: approx inverse li with one round of Halley */
  dn = Li(t) - lx;
  t = t - dn * logl(t) / (1.0L + dn/(2*t));
  /* Iterate 1-4 rounds of Halley */
  for (i = 0; i < 4; i++) {
    dn = RiemannR(t) - lx;
    term = dn * logl(t) / (1.0L + dn/(2*t));
    if (i > 0 && fabsl(term) >= fabsl(old_term)) { t -= term/4; break; }
    old_term = term;
    t -= term;
  }
  return (UV)ceill(t);
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

long double RiemannR(long double x) {
  long double part_term, term, flogx, ki, old_sum;
  unsigned int k;
  KAHAN_INIT(sum);

  if (x <= 0) croak("Invalid input to ReimannR:  x must be > 0");

  if (x > 1e19) {
    const signed char* amob = _moebius_range(0, 100);
    KAHAN_SUM(sum, Li(x));
    for (k = 2; k <= 100; k++) {
      if (amob[k] == 0) continue;
      ki = 1.0L / (long double) k;
      part_term = powl(x,ki);
      if (part_term > LDBL_MAX) return INFINITY;
      term = amob[k] * ki * Li(part_term);
      old_sum = sum;
      KAHAN_SUM(sum, term);
      if (fabsl(sum - old_sum) <= LDBL_EPSILON) break;
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
    if (fabsl(sum - old_sum) <= LDBL_EPSILON) break;
  }

  return sum;
}

static long double _lambertw_approx(long double x) {
  /* See Veberic 2009 for other approximations */
  if (x < -0.060) {  /* Pade(3,2) */
    long double ti = 5.4365636569180904707205749L * x + 2.0L;
    long double t  = (ti <= 0.0L) ? 0.0L : sqrtl(ti);
    long double t2 = t*t;
    long double t3 = t*t2;
    return (-1.0L + (1.0L/6.0L)*t + (257.0L/720.0L)*t2 + (13.0L/720.0L)*t3) / (1.0L + (5.0L/6.0L)*t + (103.0L/720.0L)*t2);
  } else if (x < 1.363) {  /* Winitzki 2003 section 3.5 */
    long double l1 = logl(1.0L+x);
    return l1 * (1.0L - logl(1.0L+l1) / (2.0L+l1));
  } else if (x < 3.7) {    /* Modification of Vargas 2013 */
    long double l1 = logl(x);
    long double l2 = logl(l1);
    return l1 - l2 - logl(1.0L - l2/l1)/2.0L;
  } else {                 /* Corless et al. 1993, page 22 */
    long double l1 = logl(x);
    long double l2 = logl(l1);
    long double d1 = 2.0L*l1*l1;
    long double d2 = 3.0L*l1*d1;
    long double d3 = 2.0L*l1*d2;
    long double d4 = 5.0L*l1*d3;
    long double w = l1 - l2 + l2/l1 + l2*(l2-2.0L)/d1;
    w += l2*(6.0L+l2*(-9.0L+2.0L*l2))/d2;
    w += l2*(-12.0L+l2*(36.0L+l2*(-22.0L+3.0L*l2)))/d3;
    w += l2*(60.0L+l2*(-300.0L+l2*(350.0L+l2*(-125.0L+12.0L*l2))))/d4;
    return w;
  }
}

long double lambertw(long double x) {
  long double w;
  int i;

  if (x < -0.36787944117145L)
    croak("Invalid input to LambertW:  x must be >= -1/e");
  if (x == 0.0L) return 0.0L;

  /* Estimate initial value */
  w = _lambertw_approx(x);
  /* If input is too small, return .99999.... */
  if (w <= -1.0L) return -1.0L + 8*LDBL_EPSILON;
  /* For very small inputs, don't iterate, return approx directly. */
  if (x < -0.36783) return w;

#if 0  /* Halley */
  lastw = w;
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
    long double zn = logl(x/w) - w;
    long double qn = 2 * w1 * (w1+(2.0L/3.0L)*zn);
    long double en = (zn/w1) * (qn-zn)/(qn-2.0L*zn);
    /* w *= 1.0L + en;  if (fabsl(en) <= 16*LDBL_EPSILON) break; */
    long double wen = w * en;
    w += wen;
    if (fabsl(wen) <= 64*LDBL_EPSILON) break;
  }
#endif

  return w;
}

char* pidigits(int digits)
{
  char* out;
  IV *a, b, c, d, e, f, g, i,  d4, d3, d2, d1;
  if (digits <= 0) return 0;
  if (digits <= DBL_DIG && digits <= 18) {
    Newz(0, out, 19, char);
    (void)sprintf(out, "%.*lf", (digits-1), 3.141592653589793238);
    return out;
  }
  digits++;   /* For rounding */
  b = d = e = g = i = 0;  f = 10000;
  c = 14*(digits/4 + 2);
  New(0, a, c, IV);
  New(0, out, digits+5+1, char);
  *out++ = '3';  /* We'll turn "31415..." into "3.1415..." */
  for (b = 0; b < c; b++)  a[b] = 20000000;

  while ((b = c -= 14) > 0 && i < digits) {
    d = e = d % f;
    while (--b > 0) {
      d = d * b + a[b];
      g = (b << 1) - 1;
      a[b] = (d % g) * f;
      d /= g;
    }
    /* sprintf(out+i, "%04d", e+d/f);   i += 4; */
    d4 = e+d/f;
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

  for (i = 0, len = strlen(s); i < len; i++) {
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
  {
    UV factors[MPU_MAX_FACTORS+1];
    int nfactors = factor_exp(n, factors, 0);
    /* Aebi and Cairns 2008, page 9 */
#if BITS_PER_WORD == 32
    if (nfactors == 2)
#else
    if (nfactors == 2 && (n < UVCONST(10000000000)))
#endif
      return 0;
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

static int numcmp(const void *a, const void *b)
  { const UV *x = a, *y = b; return (*x > *y) ? 1 : (*x < *y) ? -1 : 0; }

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
      qsort(S, k, sizeof(UV), numcmp);
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
