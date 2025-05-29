#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <float.h>

#include "ptypes.h"
#define FUNC_isqrt 1
#define FUNC_lcm_ui 1
#define FUNC_ctz 1
#define FUNC_log2floor 1
#define FUNC_is_perfect_square
#define FUNC_next_prime_in_sieve 1
#define FUNC_prev_prime_in_sieve 1
#define FUNC_ipow 1
#include "util.h"
#include "sieve.h"
#include "primality.h"
#include "cache.h"
#include "legendre_phi.h"
#include "prime_counts.h"
#include "prime_powers.h"
#include "factor.h"
#include "mulmod.h"
#include "constants.h"
#include "montmath.h"
#include "csprng.h"
#include "inverse_interpolate.h"
#include "rootmod.h"
#include "lucas_seq.h"
#include "sort.h"

static int _verbose = 0;
void _XS_set_verbose(int v) { _verbose = v; }
int _XS_get_verbose(void) { return _verbose; }

static int _call_gmp = 0;
void _XS_set_callgmp(int v) { _call_gmp = v; }
int  _XS_get_callgmp(void) { return _call_gmp; }

static bool _secure = 0;
void _XS_set_secure(void) { _secure = 1; }
bool  _XS_get_secure(void) { return _secure; }

/******************************************************************************/

/* Returns 0 if not found, index+1 if found (returns leftmost if dups) */
unsigned long index_in_sorted_uv_array(UV v, UV* L, unsigned long len)
{
  unsigned long lo, hi;
  if (len == 0 || v < L[0] || v > L[len-1])
    return 0;
  lo = 0;
  hi = len-1;
  while (lo < hi) {
    unsigned long mid = lo + ((hi-lo) >> 1);
    if (L[mid] < v)  lo = mid + 1;
    else             hi = mid;
  }
  return (L[lo] == v)  ?  lo+1  :  0;
}
unsigned long index_in_sorted_iv_array(IV v, IV* L, unsigned long len)
{
  unsigned long lo, hi;
  if (len == 0 || v < L[0] || v > L[len-1])
    return 0;
  lo = 0;
  hi = len-1;
  while (lo < hi) {
    unsigned long mid = lo + ((hi-lo) >> 1);
    if (L[mid] < v)  lo = mid + 1;
    else             hi = mid;
  }
  return (L[lo] == v)  ?  lo+1  :  0;
}

/******************************************************************************/

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

/* Return true if n is prime, false if not.  Do it fast. */
bool is_prime(UV n)
{
  if (n < UVCONST(500000000)) {

    if (n < 11) return 0xAC >> n & 1;
    if (is_divis_2_3_5_7(n)) return 0;

    /* Check static tiny sieve */
    if (n < 30*NPRIME_SIEVE30) {
      UV d = n/30,  m = n - d*30;
      return ((prime_sieve30[d] & masktab30[m]) == 0);
    }

    /* Check primary cache */
    if (n <= get_prime_cache(0,0)) {
      const unsigned char* sieve;
      int isprime = -1;
      if (!(n%11) || !(n%13)) return 0;
      if (n <= get_prime_cache(0, &sieve)) {
        UV d = n/30,  m = n - d*30;
        isprime = ((sieve[d] & masktab30[m]) == 0);
      }
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

/* We're trying to quickly give a reasonable monotonic upper prime count */
UV max_nprimes(UV n)
{
  /* 2-bit error term of the 1..726 func so 0-143 gives exact results */
  static const uint32_t _cor[9] = {0x415556af,0x01400001,0x00014140,0x01150100,0x14001515,0xa5515014,0x01555696,0xbea95501,0xeaabfaba};
  double r;

  if (n < 727)
    return (13 + n - 7*n*n/16384)/4
           -  (n < 144  ?  _cor[n/16] >> n%16*2 & 3  :  0);

  r = 1/log(n);

  if (n <     59471) /* Special */
    return (UV)(n*r * (1 + r*(1 + 2.47687*r))) + 1;

  if (n <   1333894) /* Dusart 2018   x > 1 */
    return n*r * (1 + r*(1 + 2.53816*r));

  if (n < 883495117) /* Dusart 2022   x > 1 */
    return n*r * (1 + r*(1 + r*(2 + 7.59*r)));

  /* We could use better bounds with Li(n) but that is MUCH slower. */
  /* Use prime_count_upper(n) if you want tighter bounds. */

  /* Axler 2022 x > 1    Prp 4.6 */
  return n*r * (1 + r*(1 + r*(2 + r*(6.024334 + r*(24.024334 + r*(120.12167 + r*(720.73002 + 6098*r)))))));
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

  /* Kuznetsov indicates that the Deléglise & Rivat (1996) method can be
   * modified to work on logs, which allows us to operate with no
   * intermediate memory at all.  Same time as the D&R method, less memory. */
  unsigned char logp;
  UV nextlog, nextlogi;

  if (hi < lo) croak("range_mobius error hi %"UVuf" < lo %"UVuf"\n", hi, lo);

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
    float rfactor = 1.0 + 7.0 * (1.0 - exp(-(float)maxmu/8000000000.0));
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
  unsigned char *nf;

  if (hi < lo) croak("range_liouvillle error hi %"UVuf" < lo %"UVuf"\n",hi,lo);
  nf = range_nfactor_sieve(lo, hi, 1);
  New(0, l, hi-lo+1, signed char);
  for (i = 0; i < hi-lo+1; i++)
    l[i] = (nf[i] & 1) ? -1 : 1;
  Safefree(nf);
  return l;
}

UV carmichael_lambda(UV n) {
  unsigned char _totient[8] = {0,1,1,2,2,4,2,6};
  UV fac[MPU_MAX_FACTORS+1];
  int i, nfactors;
  UV lambda = 1;

  if (n < 8) return _totient[n];
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


/******************************************************************************/
/*                             POWERS and ROOTS                               */
/******************************************************************************/

static float _cbrtf(float x)
{
  float t, r;
  union { float f; uint32_t i; } xx = { x };
  xx.i = (xx.i + 2129874493)/3;
  t = xx.f;
  /* One round of Halley's method gets to 15.53 bits */
  r = t * t * t;
  t *= (x + (x + r)) / ((x + r) + r);
#if BITS_PER_WORD > 45
  /* A second round gets us the 21.5 bits we need. */
  r = t * t * t;
  t += t * (x - r) / (x + (r + r));
#endif
  return t;
}
uint32_t icbrt(UV n) {
  if (n > 0) {
    uint32_t root = (float)(_cbrtf((float)n) + 0.375f);
    UV rem = n - (UV)root * root * root;
    return root - ((IV)rem < 0);
  }
  return 0;
}

/******************************************************************************/

static UV _ipow(unsigned b, unsigned e, unsigned bit)
{
  UV r = b;
  while (bit >>= 1) {
    r *= r;
    if (e & bit)
      r *= b;
  }
  return r;
}
/* Estimate the kth root of n.
 *
 * Returns exact root if n is a perfect power, otherwise either root or root+1.
 * Requires k >= 3 so a float can exactly represent the kth root.
 *
 * This version is heavily trimmed for internal use with rootint's prefilters.
 *
 *      n > 1
 *      n>>k != 0   <=>   n < 1<<k
 *      4  <  k  <=  MAX_IROOTN (32-bit: 10  64-bit: 15)
 */
static uint32_t _est_root(UV n, unsigned k, unsigned msbit)
{
  const float y = n;
  union { float f; uint32_t i; } both32 = { y };
  const uint32_t float_one = (uint32_t)127 << 23;
  float x, err, xk;

  if (k == 4) return (int)(sqrtf(sqrtf(y)) + 0.5f);

  /* The standard floating-point trick for an initial estimate,
   * but using two constants for variable k.  The constants
   * are chosen to be perfect for k=5 and very close to ideal
   * for k=6.  As k increases, the relative accuracy needed
   * decreases, so higher k can tolerate a lot of slop.
   *
   * One problem is that n==1 underflows.  We could fix this
   * (add k<<20 before division and subtract 1<<10 after), but
   * it's simpelr just to special case n==1. */

  both32.i = (both32.i - float_one - 89788) / k + float_one - 282298;
  x = both32.f;

  /* Improve it with one round of Halley's method.
   *
   * Newton's (quadratic) method for a root of x^k - y == 0 is
   *     x += (y - x^k) / (k * x^(k-1))
   * which simplifies a lot for fixed k, but for variable k,
   * Halley's (cubic) method is not much more complex:
   *     x += 2*x*(y-x^k) / (k*(y+x^k) - (y - x^k))
   * For all k >= 5, one round suffices.
   * Since k < 5 is handled already, this works for us. */

  xk = x;
  while (msbit >>= 1) {
    xk *= xk;
    if (k & msbit)
      xk *= x;
  }

  err = y - xk;
  x += 2.0f*x*err / ((float)(int)k*(y+xk) - err);
  return (int)(x + 0.5f);
}

/* Trimmed for internal use.  k MUST be between 4 and 15, n > 1 */
#define MAX_IROOTN ((BITS_PER_WORD == 64) ? 15 : 10)
static uint32_t _irootn(UV n, uint32_t k)
{
  uint32_t const msb = 4 << (k >= 8);
  uint32_t const r   = _est_root(n,k,msb);
  return r - ((IV)(n - _ipow(r,k,msb)) < 0);
}

/******************************************************************************/

#if BITS_PER_WORD == 64
static const uint32_t root_max[1+MPU_MAX_POW3] = {0,0,4294967295U,2642245,65535,7131,1625,565,255,138,84,56,40,30,23,19,15,13,11,10,9,8,7,6,6,5,5,5,4,4,4,4,3,3,3,3,3,3,3,3,3};
#else
static const uint32_t root_max[1+MPU_MAX_POW3] = {0,0,65535,1625,255,84,40,23,15,11,9,7,6,5,4,4,3,3,3,3,3};
#endif

UV rootint(UV n, uint32_t k)
{
  if (n <= 1) return (k != 0 && n != 0);

  switch (k) {
    case 0:  return 0;
    case 1:  return n;
    case 2:  return isqrt(n);
    case 3:  return icbrt(n);
    case 4:  return _irootn(n,4);
    case 5:  return _irootn(n,5);
    default: break;
  }

  /*            MAX_IROOTN  <  BITS_PER_WORD/2  <  MPU_MAX_POW3  */
  /*  32-bit:       10               16                 20       */
  /*  64-bit:       15               32                 40       */

  if (n >> k == 0)           return 1;

  if (k <= MAX_IROOTN)       return _irootn(n,k);

  if (k > MPU_MAX_POW3)      return 1 + (k < BITS_PER_WORD);
  if (k >= BITS_PER_WORD/2)  return 2 + (n >= ipow(3,k));

  /* k is now in range 11-15 (32-bit), 16-31 (64-bit).  Binary search. */
  {
    uint32_t lo = 1U << (log2floor(n)/k);
    uint32_t hi = root_max[k];
    if (hi >= lo*2) hi = lo*2 - 1;

    while (lo < hi) {
      uint32_t mid = lo + (hi-lo+1)/2;
      if (ipow(mid,k) > n) hi = mid-1;
      else                 lo = mid;
    }
    return lo;
  }
}

/* Like ipow but returns UV_MAX if overflow */
UV ipowsafe(UV n, UV k) {
  UV p = 1;

  if (k <= MPU_MAX_POW3) {
    if (k == 0) return 1;
    if (k == 1) return n;
    return (n <= root_max[k]) ? ipow(n,k) : UV_MAX;
  }

  while (k) {
    if (k & 1) { if (UV_MAX/n < p) return UV_MAX;  p *= n; }
    k >>= 1;
    if (k)     { if (UV_MAX/n < n) return UV_MAX;  n *= n; }
  }
  return p;
}


/******************************************************************************/

/* Mod 32 filters for allowable k-th root */
static const uint32_t _rootmask32[41] = {
  0x00000000,0x00000000,0xfdfcfdec,0x54555454,0xfffcfffc,           /* 0-4   */
  0x55555554,0xfdfdfdfc,0x55555554,0xfffffffc,0x55555554,0xfdfdfdfc,/* 5-10  */
  0x55555554,0xfffdfffc,0xd5555556,0xfdfdfdfc,0xf57d57d6,0xfffffffc,/* 11-16 */
  0xffffd556,0xfdfdfdfe,0xd57ffffe,0xfffdfffc,0xffd7ff7e,0xfdfdfdfe,/* 17-22 */
  0xffffd7fe,0xfffffffc,0xffffffd6,0xfdfffdfe,0xd7fffffe,0xfffdfffe,/* 23-28 */
  0xfff7fffe,0xfdfffffe,0xfffff7fe,0xfffffffc,0xfffffff6,0xfffffdfe,/* 29-34 */
  0xf7fffffe,0xfffdfffe,0xfff7fffe,0xfdfffffe,0xfffff7fe,0xfffffffc /* 35-40 */
};

bool is_power_ret(UV n, uint32_t k, uint32_t *root)
{
  uint32_t r, msbit;

  /* Simple edge cases */
  if (n < 2 || k == 1) {
    if (root) *root = n;
    return 1;
  }
  if (k == 0)
    return 0;
  if (k > MPU_MAX_POW3) {
    if (root) *root = 2;
    return (k < BITS_PER_WORD && n == (UV)1 << k);
  }

  if (k == 2) return is_perfect_square_ret(n,root);

  /* Filter out many numbers which cannot be k-th roots */
  if ((1U << (n&31)) & _rootmask32[k]) return 0;

  if (k == 3) {
    r = n % 117; if ((r*833230740) & (r*120676722) & 813764715) return 0;
    r = icbrt(n);
    if (root) *root = r;
    return (UV)r*r*r == n;
  }

  for (msbit = 8 /* k >= 4 */; k >= msbit; msbit <<= 1)  ;
  msbit >>= 1;
  r = _est_root(n, k, msbit);
  if (root) *root = r;
  return _ipow(r, k, msbit) == n;
}

#define PORET(base,exp)  do { \
  uint32_t n_ = base;  /* In case base uses k or exp uses n */ \
  k *= exp; \
  n = n_; \
  goto poreturn; \
} while (0)

/* max power for 64-bit inputs */
static const uint8_t _maxpow128[128] = {31,7,0,11,2,7,0,17,3,17,0,13,0,11,0,11,2,11,0,29,0,13,0,11,3,11,0,7,0,7,0,7,5,7,0,13,2,7,0,11,3,7,0,31,0,7,0,11,0,11,0,11,0,11,0,13,3,13,0,13,0,19,0,7,3,7,0,17,2,17,0,11,3,11,0,23,0,17,0,13,0,13,0,13,0,7,0,19,3,19,0,19,0,11,0,11,5,11,0,7,2,13,0,13,3,13,0,7,0,23,0,7,0,7,0,37,0,7,0,11,3,11,0,11,0,13,0,7};

/* Returns maximal k for c^k = n for k > 1, n > 1.  0 otherwise. */
uint32_t powerof_ret(UV n, uint32_t *root) {
  uint32_t r, t, k = 1;

  /* SPECIAL: For n = 0 and n = 1, return k=1 with root n. */
  /* This matches SAGE's .perfect_power(n) method (FLINT chooses k=2). */
  if (n <= 1)                         PORET(n,1);

  if ((n <= 3) || (n == UV_MAX))      return 0;
  if ((n & (n-1)) == 0)               PORET(2,ctz(n));

  while (is_perfect_square_ret(n,&r)) { n = r; k *= 2; }
  while (is_power_ret(n, 3, &r))      { n = r; k *= 3; }
  while (is_power_ret(n, 5, &r))      { n = r; k *= 5; }

  if (is_power_ret(n, 7, &r))         PORET(r,7);

  if ( !(((n%121)*0x8dd6295a) & 0x2088081) &&
       is_power_ret(n, 11, &r) )      PORET(r,11);

  /* Reject 78% of inputs as not powers of 13,17,19,... */
  if (_maxpow128[n % 128] < 13)       goto poreturn;

  if (is_power_ret(n, 13, &r))        PORET(r,13);
  if (is_power_ret(n, 17, &r))        PORET(r,17);

  if (n >= 1162261467) {
    r = t = 0;
    switch (n) {
      case UVCONST(1162261467):           t=19; r=3; break;
#if BITS_PER_WORD == 64
      case UVCONST(19073486328125):       t=19; r=5; break;
      case UVCONST(609359740010496):      t=19; r=6; break;
      case UVCONST(11398895185373143):    t=19; r=7; break;
      case UVCONST(10000000000000000000): t=19; r=10;break;
      case UVCONST(94143178827):          t=23; r=3; break;
      case UVCONST(11920928955078125):    t=23; r=5; break;
      case UVCONST(789730223053602816):   t=23; r=6; break;
      case UVCONST(68630377364883):       t=29; r=3; break;
      case UVCONST(617673396283947):      t=31; r=3; break;
      case UVCONST(450283905890997363):   t=37; r=3; break;
#endif
      default:  break;
    }
    if (t != 0) { n = r; k *= t; }
  }

 poreturn:
  if (k <= 1) return 0;
  if (root) *root = n;
  return k;
}


/******************************************************************************/


/* Like lcm_ui, but returns 0 if overflow */
UV lcmsafe(UV x, UV y) {
  y /= gcd_ui(x,y);
  if (UV_MAX/x < y) return 0;
  return x*y;
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
/* N => k^s * t   =>   s = valuation_remainder(N, k, &t); */
UV valuation_remainder(UV n, UV k, UV *r) {
  UV v;
  if      (k <= 1) { v = 0; }
  else if (k == 2) { v = ctz(n); n >>= v; }
  else {
    for (v=0;  !(n % k);  v++)
      n /= k;
  }
  *r = n;
  return v;
}

UV logint(UV n, UV b)
{
  /* UV e;  for (e=0; n; n /= b) e++;  return e-1; */
  UV v, e = 0;
  if (b == 2)
    return log2floor(n);
  if (b > n)
    return 0;
  if (n > UV_MAX/b) {
    n /= b;
    e = 1;
  }
  for (v = b; v <= n; v *= b)
    e++;
  return e;
}

unsigned char* range_issquarefree(UV lo, UV hi) {
  unsigned char* isf;
  UV i, p2, range = hi-lo+1, sqrthi = isqrt(hi);
  if (hi < lo) return 0;
  New(0, isf, range, unsigned char);
  memset(isf, 1, range);
  if (lo == 0) isf[0] = 0;

  { /* Sieve multiples of 2^2,3^2,5^2 */
    UV p = 2;
    while (p < 7 && p <= sqrthi) {
      for (p2=p*p, i = P_GT_LO(p2, p2, lo); i >= lo && i <= hi; i += p2)
        isf[i-lo] = 0;
      p += 1 + (p > 2);
    }
  }
  if (sqrthi >= 7) { /* Sieve multiples of higher prime squares */
    unsigned char* segment;
    UV seg_base, seg_low, seg_high;
    void* ctx = start_segment_primes(7, sqrthi, &segment);
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
        for (p2=p*p, i = P_GT_LO(p2, p2, lo); i >= lo && i <= hi; i += p2)
          isf[i-lo] = 0;
      END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
  }
  return isf;
}


#if BITS_PER_WORD == 32
static const uint32_t _maxpowersumn[32] = {0,92681,2343,361,116,53,30,20,14,11,8,7,6,5,4,4,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2};
#else
static const UV _maxpowersumn[64] = {0,6074000999,3810777,92681,9839,2190,745,331,175,105,69,49,36,28,22,18,15,13,11,10,9,8,7,6,6,5,5,5,4,4,4,4,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
#endif

UV powersum(UV n, UV k)
{
  UV a, a2, i, sum;

  if (n <= 1 || k == 0) return n;
  if (k >= BITS_PER_WORD || n > _maxpowersumn[k]) return 0;
  if (n == 2) return 1 + (UVCONST(1) << k);

  a = (n & 1)  ?  n*((n+1)>>1)  :  (n>>1)*(n+1);
  a2 = a*a;
  if (k == 1) return a;
  if (k == 3) return a2;

#if BITS_PER_WORD == 64
  if (k == 2 && n <=2642245) return a * (2*n+1) / 3;
  if (k == 4 && n <=   5724) return a * (2*n+1) * (3*n*(n+1)-1) / 15;
  if (k == 5 && n <=   1824) return a2 * (4*a - 1) / 3;
  if (k == 6 && n <=    482) return a * (2*n+1) * (n*((n*(n*(3*n+6)))-3)+1) /21;
  if (k == 7 && n <=    288) return a2 * (6*a2 - 4*a + 1) / 3;
  if (k == 8 && n <=    115) return a * (2*n+1) * (n*(n*(n*(n*(n*(5*n+15)+5)-15)-1)+9)-3) / 45;
#else
  /* TODO:  find the 32-bit limits */
#endif

  if (k <= 8 && k < n) {
    UV r, fac = 1;
    for (sum = 0, r = 1; r <= k; r++) {
      /* sum += factorial(r) * stirling2(k,r) * binomial(n+1,r+1); */
      sum += fac * stirling2(k,r) * binomial(n+1,r+1);;
      fac *= (r+1);
    }
    return sum;
  }

  sum = 1 + (UVCONST(1)<<k);
  for (i = 3; i <= n; i++)
    sum += ipow(i, k);
  return sum;
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
      if (v > 9) croak("Parameter '%s' must be a single decimal number",ptr);
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
  UV rem;
  if (a >= 0)  return kronecker_uu(a, b);
  if (b == 0)  return (a == 1 || a == -1) ? 1 : 0;
  s = 1;
  r = padic2(b);
  if (r) {
    if (!(a&1))  return 0;
    if ((r&1) && IS_MOD8_3OR5(a))  s = -s;
    b >>= r;
  }
  rem = (-a) % b;
  a = (rem == 0) ? 0 : b-rem;
  return kronecker_uu_sign(a, b, s);
}

int kronecker_ss(IV a, IV b) {
  if (a >= 0 && b >= 0)
    return (b & 1)  ?  kronecker_uu_sign(a, b, 1)  :  kronecker_uu(a,b);
  if (b >= 0)
    return kronecker_su(a, b);
  return kronecker_su(a, -b) * ((a < 0) ? -1 : 1);
}

#define MAX_PNPRIM ( (BITS_PER_WORD == 64) ? 15 : 9 )
#define MAX_PRIM   ( (BITS_PER_WORD == 64) ? 52 : 28 )
#if BITS_PER_WORD == 64
  static const UV _pn_prim[MAX_PNPRIM+1] =
    {1,2,6,30,210,2310,30030,510510,9699690,223092870,
     UVCONST(6469693230),UVCONST(200560490130),UVCONST(7420738134810),UVCONST(304250263527210),UVCONST(13082761331670030),UVCONST(614889782588491410)};
  static const unsigned char _prim_map[MAX_PRIM+1] =
    {0,0,1,2,2,3,3,4,4,4,4,5,5,6,6,6,6,7,7,8,8,8,8,9,9,9,9,9,9,10,10,11,11,11,11,11,11,12,12,12,12,13,13,14,14,14,14,15,15,15,15,15,15};
#else
  static const UV _pn_prim[MAX_PNPRIM+1] =
    {1,2,6,30,210,2310,30030,510510,9699690,223092870};
  static const unsigned char _prim_map[MAX_PRIM+1] =
    {0,0,1,2,2,3,3,4,4,4,4,5,5,6,6,6,6,7,7,8,8,8,8,9,9,9,9,9,9};
#endif

UV pn_primorial(UV n) {
  return (n > MAX_PNPRIM)  ?  0  :  _pn_prim[n];
}
UV primorial(UV n) {
  return (n > MAX_PRIM)  ?  0  :  _pn_prim[_prim_map[n]];
}
UV factorial(UV n) {
  UV i, r = 1;
  if ( (n > 12 && sizeof(UV) <= 4) || (n > 20 && sizeof(UV) <= 8) ) return 0;
  for (i = 2; i <= n; i++)
    r *= i;
  return r;
}
UV subfactorial(UV n) {
  if (n <= 3)  return (n ? n-1 : 1);
  if (n >= ((BITS_PER_WORD == 64) ? 21 : 14))  return 0;
  return (n * subfactorial(n-1) + ((n & 1) ? -1 : 1));
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

UV fubini(UV n) {
  UV k, sum;
  if (n == 0) return 1;
  if (n >= ((BITS_PER_WORD == 64) ? 16 : 10))  return 0;
  for (sum = 1, k = 2; k <= n; k++)
    sum += factorial(k) * stirling2(n, k);
  return sum;
}

UV falling_factorial(UV n, UV m)
{
  UV i, r = n;
  if (m == 0) return 1;
  if (m > n) return 0;
  for (i = 1; i < m; i++) {
    if (UV_MAX/(n-1) < r) return UV_MAX;  /* Overflow */
    r *= (n-i);
  }
  return r;
}
UV rising_factorial(UV n, UV m)
{
  if (m == 0) return 1;
  if ((m-1) > (UV_MAX-n)) return UV_MAX;  /* Overflow */
  return falling_factorial(n+m-1, m);
}

IV falling_factorial_s(IV n, UV m)
{
  UV r = (n>=0) ? falling_factorial(n,m) : rising_factorial(-n,m);
  if (r >= IV_MAX) return IV_MAX;  /* Overflow */
  return (n < 0 && (m&1)) ? -(IV)r : (IV)r;
}
IV rising_factorial_s(IV n, UV m)
{
  UV r = (n>=0) ? rising_factorial(n,m) : falling_factorial(-n,m);
  if (r >= IV_MAX) return IV_MAX;  /* Overflow */
  return (n < 0 && (m&1)) ? -(IV)r : (IV)r;
}

bool is_cyclic(UV n) {
  UV phi, facs[MPU_MAX_FACTORS+1];
  int i, nfacs;

  if (n < 4) return (n != 0);

  /* Fast filters for necessary conditions */
  if (   !(n & 1)                                              /* 2 only even */
      || !(n% 9) || !(n%25) || !(n%49)                         /* not sq free */
      || !(n%21) || !(n%39) || !(n%55) || !(n%57) || !(n%93)   /* q = 1 mod p */
      || !(n%121) || !(n%169)                                  /* not sq free */
      || !(n%111) || !(n%129) || !(n%155) || !(n%183))         /* q = 1 mod p */
    return 0;

  if (n <= 200) return 1;    /* Filters above were sufficient for tiny inputs */

  /* return gcd_ui(n, totient(n)) == 1; */

  nfacs = factor(n, facs);
  if (nfacs == 1)
    return 1;                        /* prime => cyclic */
  for (i = 1; i < nfacs; i++)
    if (facs[i] == facs[i-1])
      return 0;                      /* repeated factor => not cyclic */
  for (phi = 1, i = 0; i < nfacs; i++)
    phi *= facs[i]-1;
  return gcd_ui(n, phi) == 1;        /* cyclic <=> coprime with totient */
}

bool is_carmichael(UV n) {
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

static bool is_quasi_base(int nfactors, UV *fac, UV p, UV b) {
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
    divs = divisor_list(n / spf - 1, &ndivisors, UV_MAX);
    for (i = 0; i < (int)ndivisors; i++) {
      UV d = divs[i];
      UV k = spf - d;
      if (d >= spf) break;
      if (is_quasi_base(nfactors, fac, n-k, k))
        nbases++;
    }
  } else {
    divs = divisor_list(lpf * (n / lpf - 1), &ndivisors, UV_MAX);
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

bool is_semiprime(UV n) {
  UV sp, p, factors[2];
  uint32_t n2, n3;

  if (n < 6) return (n == 4);
  if (!(n&1)) return is_prob_prime(n>>1);
  if (!(n%3)) return is_prob_prime(n/3);
  if (!(n%5)) return is_prob_prime(n/5);
  /* 27% of random inputs left */
  n3 = icbrt(n);
  for (sp = 4; sp < 60; sp++) {
    p = primes_tiny[sp];
    if (p > n3)
      break;
    if ((n % p) == 0)
      return is_prob_prime(n/p);
  }
  /* 9.8% of random inputs left */
  if (is_def_prime(n)) return 0;
  if (p > n3) return 1; /* past this, n is a composite and larger than p^3 */
  /* 4-8% of random inputs left */

  if (is_perfect_square_ret(n,&n2)) /* Fast square check */
    return is_def_prime(n2);

  /* Find one factor, check primality of factor and co-factor */
  if (factor_one(n, factors, 0, 0) != 2) return 0;
  return (is_def_prime(factors[0]) && is_def_prime(factors[1]));
}
bool is_almost_prime(UV k, UV n) {
  UV p, sp;

  if (k == 0) return (n == 1);
  if (k == 1) return is_prob_prime(n);
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
  if (k == 1) return is_prob_prime(n);
  if (k == 2) return is_semiprime(n);
  if (n < ipowsafe(p,k)) return 0;

  return ((UV)prime_bigomega(n) == k);
}

bool is_fundamental(UV n, bool neg) {
  uint32_t r = n & 15;
  if (r) {
    if (neg) r = 16-r;
    if ((r & 3) == 0 && r != 4) return is_square_free(n >> 2);
    if ((r & 3) == 1)           return is_square_free(n);
  }
  return 0;
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


#define MOB_TESTP(p) \
  { uint32_t psq = p*p;  if (n >= psq && (n % psq) == 0) return 0; }

/* mpu 'for (0..255) { $x=moebius($_)+1; $b[$_ >> 4] |= ($x << (2*($_%16))); } say join ",",@b;' */
static const uint32_t _smoebius[16] = {2703565065,23406865,620863913,1630114197,157354249,2844895525,2166423889,363177345,2835441929,2709852521,1095049497,92897577,1772687649,162113833,160497957,689538385};
int moebius(UV n) {
  if (n < 256)  return (int)((_smoebius[n >> 4] >> (2*(n % 16))) & 3) - 1;

  if (!(n % 4) || !(n % 9) || !(n % 25) || !(n % 49) || !(n %121) || !(n %169))
    return 0;

  MOB_TESTP(17); MOB_TESTP(19); MOB_TESTP(23);
  MOB_TESTP(29); MOB_TESTP(31); MOB_TESTP(37);

  {
    UV factors[MPU_MAX_FACTORS+1];
    int i, nfactors;
    nfactors = factor(n, factors);
    for (i = 1; i < nfactors; i++)
      if (factors[i] == factors[i-1])
        return 0;
    return (nfactors % 2) ? -1 : 1;
  }
}

#define ISF_TESTP(p) \
  { uint32_t psq = p*p;  if (psq > n) return 1;  if ((n % psq) == 0) return 0; }

static const uint32_t _isf[8] = {3840601326,1856556782,3941394158,2362371810,3970362990,3471729898,4008603310,3938642668};
bool is_square_free(UV n) {
  if (n < 256)  return (_isf[n >> 5] & (1U << (n % 32))) != 0;

  if (!(n % 4) || !(n % 9) || !(n % 25) || !(n % 49) || !(n %121) || !(n %169))
    return 0;

  ISF_TESTP(17); ISF_TESTP(19); ISF_TESTP(23);
  ISF_TESTP(29); ISF_TESTP(31); ISF_TESTP(37);

  /* return (moebius(n) != 0); */
  {
    UV factors[MPU_MAX_FACTORS+1];
    int i, nfactors;
    nfactors = factor(n, factors);
    for (i = 1; i < nfactors; i++)
      if (factors[i] == factors[i-1])
        return 0;
    return 1;
  }
}

bool is_perfect_number(UV n) {
  UV v, m;
  if (n == 0 || (n & 1)) return 0;

  v = valuation(n,2);
  m = n >> v;
  if (m & (m+1)) return 0;
  if ((m >> v) != 1) return 0;
  return is_mersenne_prime(v+1);
}

UV exp_mangoldt(UV n) {
  UV p;
  if (!prime_power(n,&p)) return 1;     /* Not a prime power */
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
      if (!sqrtmod(0, a, n))
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

bool is_qr(UV a, UV n) {
  bool res;
  if (n == 0) return (a == 1);    /* Should return undef */
  if (n <= 2) return 1;
  if (a >= n) a %= n;
  if (a <= 1) return 1;

  if (is_prob_prime(n)) {
    res = (kronecker_uu(a,n) == 1);
  } else {
    UV fac[MPU_MAX_FACTORS+1], exp[MPU_MAX_FACTORS+1];
    int i, nfactors;

    nfactors = factor_exp(n, fac, exp);
    for (i = 0, res = 1;  res && i < nfactors;  i++) {
      if (exp[i] == 1 && (fac[i] == 2 || gcd_ui(a,fac[i]) != 1))
        res = 1;
      else if (exp[i] == 1 || (fac[i] != 2 && gcd_ui(a,fac[i]) == 1))
        res = (kronecker_uu(a,fac[i]) == 1);
      else {
        res = sqrtmod(0, a, ipow(fac[i],exp[i]));
      }
    }
  }
  return res;
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
  UV p, phi, a, psquared;
  int i, nfactors, isneven, ispow;
  uint32_t root;

  if (n <= 4) return (n == 0) ? 0 : n-1;
  if (n % 4 == 0)  return 0;

  isneven = !(n & 1);
  if (isneven)  n >>= 1;

  ispow = powerof_ret(n,&root);
  p = ispow ? root : n;
  if (p == 3 && isneven) return 5;
  if (!is_prob_prime(p)) return 0;

  phi = p-1;  /* p an odd prime */
  psquared = ispow ? p*p : 0;

  nfactors = factor_exp(phi, fac, 0);
  for (i = 1; i < nfactors; i++)
    phi_div_fac[i] = phi / fac[i];

#if USE_MONTMATH
 {
  UV r;
  const uint64_t npi = mont_inverse(p),  mont1 = mont_get1(p);
  for (a = 2; a < p; a++) {
    if (isneven && !(a&1)) continue;
    if (a == 4 || a == 8 || a == 9) continue;  /* Skip some perfect powers */
    if (kronecker_uu(a, p) != -1) continue;
    r = mont_geta(a, p);
    for (i = 1; i < nfactors; i++)
      if (mont_powmod(r, phi_div_fac[i], p) == mont1)
        break;
    if (i == nfactors)
      if (!ispow || powmod(a, phi, psquared) != 1)
        return a;
  }
 }
#else
  for (a = 2; a < p; a++) {
    if (isneven && !(a&1)) continue;
    if (a == 4 || a == 8 || a == 9) continue;  /* Skip some perfect powers */
    if (kronecker_uu(a, p) != -1) continue;
    for (i = 1; i < nfactors; i++)
      if (powmod(a, phi_div_fac[i], p) == 1)
        break;
    if (i == nfactors)
      if (!ispow || powmod(a, phi, psquared) != 1)
        return a;
  }
#endif
  return 0;
}

bool is_primitive_root(UV a, UV n, bool nprime) {
  UV p, phi, fac[MPU_MAX_FACTORS+1];
  int i, nfacs;

  /* Trivial but very slow:  return totient(n) == znorder(a,n) */

  if (n <= 1) return n;
  if (a >= n) a %= n;
  if (a == 0) return (n == 1);
  if (a == 1) return (n <= 2);
  if (n <= 4) return a == n-1;
  if (n % 4 == 0)  return 0;

  if (!(n&1)) {             /* If n is even, */
    if (!(a&1)) return 0;   /* 'a' cannot also be even */
    n >>= 1;                /* since 'a' is odd, it is also a root of p^k */
  }

  if (is_perfect_square(a)) return 0;
  if (gcd_ui(a,n) != 1) return 0;

  if (!nprime) {
    UV k = prime_power(n, &p);
    if (!k) return 0;  /* Not a prime power */
    n = p;
    /* Check if a isn't a root for a power, only two known <= 10^16 */
    if (k > 1 && powmod(a, p-1, p*p) == 1) return 0;
  }
  if (kronecker_uu(a,n) != -1) return 0;
  phi = n-1;

  /* a^x can be a primitive root only if gcd(x,phi) = 1. */
  /* Checking powerof(a) will typically take more time than it saves. */
  /* We already checked 'a' not a perfect square */
  if (is_power(a,3) && gcd_ui(3,phi) != 1) return 0;

#if USE_MONTMATH
  if (n & 1) {
    const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);
    a = mont_geta(a, n);
    /* Quick check for small factors before full factor */
    if ((phi % 2) == 0 && mont_powmod(a, phi/2, n) == mont1) return 0;
    if ((phi % 3) == 0 && mont_powmod(a, phi/3, n) == mont1) return 0;
    if ((phi % 5) == 0 && mont_powmod(a, phi/5, n) == mont1) return 0;
    nfacs = factor_exp(phi, fac, 0);
    for (i = 0; i < nfacs; i++)
      if (fac[i] > 5 && mont_powmod(a, phi/fac[i], n) == mont1)
        return 0;
  } else
#endif
  {
    /* Quick check for small factors before full factor */
    if ((phi % 2) == 0 && powmod(a, phi/2, n) == 1) return 0;
    if ((phi % 3) == 0 && powmod(a, phi/3, n) == 1) return 0;
    if ((phi % 5) == 0 && powmod(a, phi/5, n) == 1) return 0;
    /* Complete factor and check each one not found above. */
    nfacs = factor_exp(phi, fac, 0);
    for (i = 0; i < nfacs; i++)
      if (fac[i] > 5 && powmod(a, phi/fac[i], n) == 1)
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
UV gcddivmod(UV a, UV b, UV n) {
  UV g = gcd_ui(a,b);
  UV binv = modinverse(b/g, n);
  if (binv == 0)  return 0;
  return mulmod(a/g, binv, n);
}

/* In C89, the division and modulo operators are implementation-defined
 * for negative inputs.  C99 fixed this. */
#if __STDC_VERSION__ >= 199901L
  #define _tdivrem(q,r, D,d)   q = D/d, r = D % d
#else
  #define _tdivrem(q,r, D,d) \
    q = ((D>=0) ? ( (d>=0) ? D/d : -(D/-d) ) \
                : ( (d>=0) ? -(-D/d) : (-D/-d) ) ), \
    r = D - d*q
#endif

IV tdivrem(IV *Q, IV *R, IV D, IV d) {
  IV q,r;
  _tdivrem(q,r,D,d);
  if (Q) *Q=q;
  if (R) *R=r;
  return r;
}
IV fdivrem(IV *Q, IV *R, IV D, IV d) {
  IV q,r;
  _tdivrem(q,r,D,d);
  if ((r > 0 && d < 0) || (r < 0 && d > 0)) { q--; r += d; }
  if (Q) *Q=q;
  if (R) *R=r;
  return r;
}
IV cdivrem(IV *Q, IV *R, IV D, IV d) {
  IV q,r;
  _tdivrem(q,r,D,d);
  if (r != 0 && ((D >= 0) == (d >= 0))) { q++; r -= d; }
  if (Q) *Q=q;
  if (R) *R=r;
  return r;
}
IV edivrem(IV *Q, IV *R, IV D, IV d) {
  IV q,r;
  _tdivrem(q,r,D,d);
  if (r < 0) {
    if (d > 0) { q--; r += d; }
    else       { q++; r -= d; }
  }
  if (Q) *Q=q;
  if (R) *R=r;
  return r;
}

UV ivmod(IV a, UV n) {   /* a mod n with signed a (0 <= r < n) */
  if (a >= 0) {
    return (UV)(a) % n;
  } else {
    UV r = (UV)(-a) % n;
    return (r == 0)  ?  0  :  n-r;
  }
}

#if 0
int is_regular(UV a, UV n) {  /* there exists an x s.t. a^2*x = a mod n */
  UV d;
  if (a == 0) return 1;
  d = gcd_ui(a, n);
  return ( (d % n) == 0 && gcd_ui(d, n/d) == 1);
}
#endif


/******************************************************************************/
/*                                  N! MOD M                                  */
/******************************************************************************/

static UV _powersin(UV p, UV d) {
  UV td = d/p, e = td;
  do { td/=p; e += td; } while (td > 0);
  return e;
}

static UV _facmod(UV n, UV m) {
  UV i, res = 1;

  if (n < 1000) {

    for (i = 2; i <= n && res != 0; i++)
      res = mulmod(res,i,m);

  } else {

    unsigned char* segment;
    UV seg_base, seg_low, seg_high;
    UV sqn = isqrt(n), nsqn = n/sqn, j = sqn, nlo = 0, nhi = 0, s1 = 1;
    void* ctx = start_segment_primes(7, n, &segment);

    for (i = 1; i <= 3; i++) {  /* Handle 2,3,5 assume n>=25*/
      UV p = primes_tiny[i];
      res = mulmod(res, powmod(p,_powersin(p, n),m), m);
    }
    while (res!=0 && next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
        if (p <= nsqn) {
          res = mulmod(res, powmod(p,_powersin(p,n),m), m);
        } else {
          while (p > nhi) {
            res = mulmod(res, powmod(s1,j,m), m);
            s1 = 1;
            j--;
            nlo = n/(j+1)+1;
            nhi = n/j;
          }
          if (p >= nlo)
            s1 = mulmod(s1, p, m);
        }
        if (res == 0) break;
      END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
    res = mulmod(res, s1, m);

  }

  return res;
}
#if USE_MONTMATH
static UV _facmod_mont(UV n, UV m) {
  const uint64_t npi = mont_inverse(m),  mont1 = mont_get1(m);
  uint64_t monti = mont1;
  UV i, res = mont1;

  if (n < 1000) {

    for (i = 2; i <= n && res != 0; i++) {
      monti = addmod(monti,mont1,m);
      res = mont_mulmod(res,monti,m);
    }

  } else {

    unsigned char* segment;
    UV seg_base, seg_low, seg_high;
    UV sqn = isqrt(n), nsqn = n/sqn, j = sqn, nlo = 0, nhi = 0;
    UV s1 = mont1;
    void* ctx = start_segment_primes(7, n, &segment);

    for (i = 1; i <= 3; i++) {  /* Handle 2,3,5 assume n>=25*/
      UV p = primes_tiny[i];
      UV mp = mont_geta(p,m);
      res = mont_mulmod(res, mont_powmod(mp,_powersin(p,n),m), m);
    }
    while (res!=0 && next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
        UV mp = mont_geta(p,m);
        if (p <= nsqn) {
          res = mont_mulmod(res, mont_powmod(mp,_powersin(p,n),m), m);
        } else {
          while (p > nhi) {
            res = mont_mulmod(res, mont_powmod(s1,j,m), m);
            s1 = mont1;
            j--;
            nlo = n/(j+1)+1;
            nhi = n/j;
          }
          if (p >= nlo)
            s1 = mont_mulmod(s1, mp, m);
        }
        if (res == 0) break;
      END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
    res = mont_mulmod(res, s1, m);

  }

  res = mont_recover(res, m);
  return res;
}
#endif

UV factorialmod(UV n, UV m) {  /*  n! mod m */
  UV i, d = n, res = 1;
  bool m_prime;

  if (n >= m || m == 1) return 0;
  if (n <= 1 || m == 2) return (n <= 1);

  if (n <= 10) { /* Keep things simple for small n */
    for (i = 2; i <= n && res != 0; i++)
      res = (res * i) % m;
    return res;
  }

  m_prime = is_prime(m);
  if (n > m/2 && m_prime)    /* Check if we can go backwards */
    d = m-n-1;
  if (d < 2)
    return (d == 0) ? m-1 : 1;   /* Wilson's Theorem: n = m-1 and n = m-2 */

  if (d > 100 && !m_prime) {   /* Check for composite m that leads to 0 */
    UV maxpk = 0, fac[MPU_MAX_FACTORS+1], exp[MPU_MAX_FACTORS+1];
    int j, nfacs = factor_exp(m, fac, exp);
    for (j = 0; j < nfacs; j++) {
      fac[j] = fac[j] * exp[j];   /* Possibly too high if exp[j] > fac[j] */
      if (fac[j] > maxpk)
        maxpk = fac[j];
    }
    /* Maxpk is >= S(m), the Kempner number A002034 */
    if (n >= maxpk)
      return 0;
  }

#if USE_MONTMATH
  if (m & 1) {
    res = _facmod_mont(d, m);
  } else
#endif
  {
    res = _facmod(d, m);
  }

  if (d != n && res != 0) {      /* Handle backwards case */
    if (!(d&1)) res = submod(m,res,m);
    res = modinverse(res,m);
  }

  return res;
}


/******************************************************************************/
/*                          BINOMIAL(N,K) MOD M                               */
/******************************************************************************/

static UV _factorial_valuation(UV n, UV p) {
  UV k = 0;
  while (n >= p) {
    n /= p;
    k += n;
  }
  return k;
}
static int _binoval(UV n, UV k, UV m) {
  return _factorial_valuation(n,m) - _factorial_valuation(k,m) - _factorial_valuation(n-k,m);
}
static UV _factorialmod_without_prime(UV n, UV p, UV m) {
  UV i, pmod, r = 1;
  MPUassert(p >= 2 && m >= p && (m % p) == 0, "_factorialmod called with wrong args");
  if (n <= 1) return 1;

  if (n >= m) {
    /* Note with p=2 the behaviour is different */
    if ( ((n/m) & 1) && (p > 2 || m == 4) )  r = m-1;
    n %= m;
  }

#if USE_MONTMATH
  if (m & 1) {
    const uint64_t npi = mont_inverse(m),  mont1 = mont_get1(m);
    uint64_t mi = mont1;
    r = mont_geta(r, m);
    for (i = pmod = 2; i <= n; i++) {
      mi = addmod(mi, mont1, m);
      if (pmod++ == p) pmod = 1;
      else             r = mont_mulmod(r, mi, m);
    }
    r = mont_recover(r, m);
  } else
#endif
  {
    for (i = pmod = 2; i <= n; i++) {
      if (pmod++ == p) pmod = 1;
      else             r = mulmod(r, i, m);
    }
  }
  return r;
}
static UV _factorialmod_without_prime_powers(UV n, UV p, UV m) {
  UV ip, r = 1;

  for (ip = n; ip > 1; ip /= p)
    r = mulmod(r, _factorialmod_without_prime(ip, p, m), m);

  return r;
}
static UV _binomial_mod_prime_power(UV n, UV k, UV p, UV e) {
  UV r, b, m, i, num, den, ip, ires;

  if (k > n) return 0;
  if (k == 0 || k == n) return 1;
  if (k > n/2)  k = n-k;

  b = _binoval(n,k,p);
  if (e <= b) return 0;
  m = ipow(p,e);

  if (k == 1) return n % m;

  /* Both methods work fine -- choose based on performance. */
  den  = _factorialmod_without_prime_powers(k, p, m);
  if (k >= m) {
    num  = _factorialmod_without_prime_powers(n, p, m);
    ip   = _factorialmod_without_prime_powers(n-k, p, m);
    den = mulmod(den, ip, m);
  } else {
#if USE_MONTMATH
    if (m & 1) {
      const uint64_t npi = mont_inverse(m),  mont1 = mont_get1(m);
      num = mont1;
      for (i = n-k+1, ires = (i-1)%p; i <= n; i++) {
        ip = i;
        if (++ires == p) { ires = 0; do { ip /= p; } while ((ip % p) == 0); }
        num = mont_mulmod(num, mont_geta(ip, m), m);
      }
      num = mont_recover(num, m);
    } else
#endif
    {
      num = 1;
      for (i = n-k+1, ires = (i-1) % p; i <= n; i++) {
        ip = i;
        if (++ires == p) { ires = 0; do { ip /= p; } while ((ip % p) == 0); }
        num = mulmod(num, ip, m);
      }
    }
  }

  r = divmod(num, den, m);
  if (b > 0) r = mulmod(r, ipow(p,b), m);
  return r;
}

static UV _binomial_lucas_mod_prime(UV n, UV k, UV p) {
  UV res, t, vn[BITS_PER_WORD], vk[BITS_PER_WORD];
  int i, ln, lk;

  if (p < 2) return 0;
  if (p == 2) return !(~n & k);

  for (t = n, ln = 0; t > 0; t /= p)
    vn[ln++] = t % p;
  for (t = k, lk = 0; t > 0; t /= p)
    vk[lk++] = t % p;

  res = 1;
  for (i = ln-1; i >= 0; i--) {
    UV ni = vn[i];
    UV ki = (i < lk) ? vk[i] : 0;
    res = mulmod(res, _binomial_mod_prime_power(ni, ki, p, 1), p);
  }
  return res;
}

/* Based on Granville's paper on the generalization of Lucas's theorem to
 * prime powers: https://www.dms.umontreal.ca/~andrew/Binomial/genlucas.html
 * and Max Alekseyev's binomod.gp program. */
static UV _binomial_lucas_mod_prime_power(UV n, UV k, UV p, UV q) {
  UV N[BITS_PER_WORD], K[BITS_PER_WORD], R[BITS_PER_WORD], e[BITS_PER_WORD];
  UV i, d, m, n1, k1, r1, m1, res;

  MPUassert(q < BITS_PER_WORD, "bad exponent in binomialmod generalized lucas");
  m = ipow(p, q);

  /* Construct the digits for N, K, and N-K (R). */
  n1 = n;   k1 = k;  r1 = n-k;
  for (d = 0; n1 > 0; d++) {
    N[d] = n1 % p;  n1 /= p;
    K[d] = k1 % p;  k1 /= p;
    R[d] = r1 % p;  r1 /= p;
  }
  /* Compute the number of carries. */
  for (i = 0; i < d; i++)
    e[i] = (N[i] < (K[i] + ((i > 0) ? e[i-1] : 0)));
  /* Turn the carries into a cumulative count. */
  for (i = d-1; i >= 1; i--)
    e[i-1] += e[i];

  if (e[0] >= q) return 0;
  q -= e[0];
  m1 = ipow(p, q);

  /* Now make the digits for the reduced N, K, N-K */
  n1 = n;   k1 = k;  r1 = n-k;
  for (d = 0; n1 > 0; d++) {
    N[d] = n1 % m1;  n1 /= p;
    K[d] = k1 % m1;  k1 /= p;
    R[d] = r1 % m1;  r1 /= p;
  }

  /* Theorem 1 from Granville indicates the +/- 1.  */
  res = ((p > 2 || q < 3) && q < d && e[q-1] % 2)  ?  m-1  :  1;
  res = mulmod(res, powmod(p, e[0], m), m);

  /* Compute the individual binomials (again, theorem 1) */
  for (i = 0; i < d; i++) {
    UV ni = _factorialmod_without_prime(N[i], p, m);
    UV ki = _factorialmod_without_prime(K[i], p, m);
    UV ri = _factorialmod_without_prime(R[i], p, m);
    UV r = divmod(ni, mulmod(ki, ri, m), m);
    res = mulmod(res, r, m);
  }
  return res;
}

bool binomialmod(UV *res, UV n, UV k, UV m) {

  if (m <= 1)           { *res = 0; return 1; }
  if (k == 0 || k >= n) { *res = (k == 0 || k == n); return 1; }

  if (m == 2) { *res = !(~n & k); return 1; }

#if 0
    if ( (*res = binomial(n,k)) )
      { *res %= m; return 1; }
#endif

  if (is_prime(m)) {
    *res = _binomial_lucas_mod_prime(n, k, m);
    return 1;
  }
  {
    UV bin[MPU_MAX_FACTORS+1];
    UV fac[MPU_MAX_FACTORS+1];
    UV exp[MPU_MAX_FACTORS+1];
    int i, cret, nfactors = factor_exp(m, fac, exp);
    for (i = 0; i < nfactors; i++) {
      if (exp[i] == 1) {
        bin[i] = _binomial_lucas_mod_prime(n, k, fac[i]);
      } else {
        /* bin[i] = _binomial_mod_prime_power(n, k, fac[i], exp[i]); */
        /* Use generalized Lucas */
        bin[i] = _binomial_lucas_mod_prime_power(n, k, fac[i], exp[i]);
        fac[i] = ipow(fac[i], exp[i]);
      }
    }
    /* chinese with p^e as modulos, so should never get -1 back */
    cret = chinese(res, 0, bin, fac, nfactors);
    return (cret == 1);
  }
}

/* Pisano period.  */
/* Thanks to Trizen & Charles R Greathouse IV for ideas and working examples. */
/* Algorithm from Charles R Greathouse IV, https://oeis.org/A001175 */
static UV _pisano_prime_power(UV p, UV e)
{
  UV k;
  if (e == 0) return 1;
  if (p == 2) return 3UL << (e-1);
  if      (p == 3) k = 8;
  else if (p == 5) k = 20;
  else if (p == 7) k = 16;
  else if (p < 300) {          /* Simple search */
    UV a = 1,b = 1, t;
    k = 1;
    while (!(a == 0 && b == 1)) {
      k++;
      t = b; b = addmod(a,b,p); a = t;
    }
  } else {                     /* Look through divisors of p-(5|p) */
    int i, nfactors;
    UV j, fac[MPU_MAX_FACTORS+1], exp[MPU_MAX_FACTORS+1];

    k = p - kronecker_uu(5,p);
    nfactors = factor_exp(k, fac, exp);
    for (i = 0; i < nfactors; i++) {
      for (j = 0; j < exp[i]; j++) {
        if (lucasumod(1, p-1, k/fac[i], p) != 0) break;
        k /= fac[i];
      }
    }
  }
  return (e == 1)  ?  k  :  k * ipow(p, e-1);
}
UV pisano_period(UV n)
{
  UV r, lim, k, fac[MPU_MAX_FACTORS+1], exp[MPU_MAX_FACTORS+1];
  int i, nfactors;

  if (n <= 1) return (n == 1);

  nfactors = factor_exp(n, fac, exp);
  for (i = 0, k = 1; i < nfactors; i++) {
    k = lcmsafe(k, _pisano_prime_power(fac[i], exp[i]));
    if (k == 0) return 0;
  }

  /* Do this carefully to avoid overflow */
  r = 0;
  lim = (UV_MAX/6 < n) ? UV_MAX : 6*n;
  do {
    r += k;
    if (lucasumod(1, n-1, r-1, n) == 1)
      return r;
  } while (r <= (lim-k));

  return 0;
}

/******************************************************************************/
/*                                   HAPPY                                    */
/******************************************************************************/

static UV sum_of_digits(UV n, uint32_t base, uint32_t k) {
  UV t, r, sum = 0;
  while (n) {
    t = n / base;
    r = n - base * t;
    switch (k) {
      case 0:  sum += 1;         break;
      case 1:  sum += r;         break;
      case 2:  sum += r*r;       break;
      default: sum += ipow(r,k); break;
    }
    n = t;
  }
  return sum;
}
static UV sum_of_squared_digits(UV n) {
  UV t, r, sum = 0;
  while (n) {
    t = n / 10;
    r = n - 10 * t;
    sum += r*r;
    n = t;
  }
  return sum;
}

int happy_height(UV n, uint32_t base, uint32_t exponent) {
  int h;

  if (base == 10 && exponent == 2) {
    static const char sh[101] = {0,1,0,0,0,0,0,6,0,0,2,0,0,3,0,0,0,0,0,5,0,0,0,4,0,0,0,0,4,0,0,3,4,0,0,0,0,0,0,0,0,0,0,0,5,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,6,0,0,0,0,0,0,0,0,4,0,0,4,0,0,0,3,0,0,0,0,5,0,0,5,0,0,4,0,0,2};
    for (h = 0;  n > 100;  h++)
      n = sum_of_squared_digits(n);
    return (sh[n] == 0) ? 0 : h+sh[n];
  } else {
    UV ncheck = 0;
    for (h = 1;  n > 1 && n != ncheck;  h++) {
      if ((h & (h-1)) == 0) ncheck = n;         /* Brent cycle finding */
      n = sum_of_digits(n, base, exponent);
    }
  }
  return (n == 1) ? h : 0;
}


/******************************************************************************/
/*                                    CRT                                     */
/******************************************************************************/

/* works only for co-prime inputs and also slower than the algorithm below,
 * but handles the case where IV_MAX < lcm <= UV_MAX.
 * status = 1 means good result, 0 means try another method.
 */
static bool _simple_chinese(UV *r, UV *mod, const UV* a, const UV* n, UV num) {
  UV i, lcm = 1, res = 0;
  if (num == 0) { *r = 0; if (mod) *mod = 0; return 1; }  /* Dubious return */

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
  *r = res;
  if (mod) *mod = lcm;
  return 1;
}

/* status: 1 ok, -1 no inverse, 0 overflow */
int chinese(UV *r, UV *mod, UV* a, UV* n, UV num) {
  static unsigned short sgaps[] = {7983,3548,1577,701,301,132,57,23,10,4,1,0};
  UV gcd, i, j, lcm, sum, gi, gap;
  if (num == 0) { *r = 0; if (mod) *mod = 0; return 1; }  /* Dubious return */

  /* Sort modulii, largest first */
  for (gi = 0, gap = sgaps[gi]; gap >= 1; gap = sgaps[++gi]) {
    for (i = gap; i < num; i++) {
      UV tn = n[i], ta = a[i];
      for (j = i; j >= gap && n[j-gap] < tn; j -= gap)
        {  n[j] = n[j-gap];  a[j] = a[j-gap];  }
      n[j] = tn;  a[j] = ta;
    }
  }

  if (n[num-1] == 0) return -1;  /* mod 0 */
  if (n[0] > IV_MAX) return _simple_chinese(r,mod,a,n,num);
  lcm = n[0]; sum = a[0] % n[0];
  for (i = 1; i < num; i++) {
    IV u, v, t, s;
    UV vs, ut;
    gcd = gcdext(lcm, n[i], &u, &v, &s, &t);
    if (gcd != 1 && ((sum % gcd) != (a[i] % gcd))) return -1;
    if (s < 0) s = -s;
    if (t < 0) t = -t;
    if (s > (IV)(IV_MAX/lcm)) return _simple_chinese(r,mod,a,n,num);
    lcm *= s;
    if (u < 0) u += lcm;
    if (v < 0) v += lcm;
    vs = mulmod((UV)v, (UV)s, lcm);
    ut = mulmod((UV)u, (UV)t, lcm);
    sum = addmod(  mulmod(vs, sum, lcm),  mulmod(ut, a[i], lcm),  lcm  );
  }
  *r = sum;
  if (mod) *mod = lcm;
  return 1;
}

bool prep_pow_inv(UV *a, UV *k, int kstatus, UV n) {
  if (n == 0) return 0;
  if (kstatus < 0) {
    if (*a != 0) *a = modinverse(*a, n);
    if (*a == 0) return 0;
    *k = -(IV)*k;
  }
  return 1;
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

static int strnum_parse(const char **sp, STRLEN *slen)
{
  const char* s = *sp;
  STRLEN i = 0, len = *slen;
  int neg = 0;

  if (s != 0 && len > 0) {
    neg = (s[0] == '-');
    if (s[0] == '-' || s[0] == '+') { s++; len--; }
    while (len > 0 && *s == '0') { s++; len--; }
    for (i = 0; i < len; i++)
      if (!isDIGIT(s[i]))
        break;
  }
  if (s == 0 || len == 0 || i < len) croak("Parameter must be an integer");
  *sp = s;
  *slen = len;
  return neg;
}
int strnum_cmp(const char* a, STRLEN alen, const char* b, STRLEN blen) {
  STRLEN i;
  int aneg = strnum_parse(&a, &alen);
  int bneg = strnum_parse(&b, &blen);
  if (aneg != bneg)  return (bneg) ? 1 : -1;
  if (aneg) { /* swap a and b if both negative */
    const char* t = a;  STRLEN tlen = alen;
    a = b; b = t;  alen = blen;  blen = tlen;
  }
  if (alen != blen)  return (alen > blen) ? 1 : -1;
  for (i = 0; i < blen; i++)
    if (a[i] != b[i])
      return  (a[i] > b[i]) ? 1 : -1;
  return 0;
}

/* 1. Perform signed integer validation on b/blen.
 * 2. Compare to a/alen using min or max based on first arg.
 * 3. Return 0 to select a, 1 to select b.
 */
bool strnum_minmax(bool min, const char* a, STRLEN alen, const char* b, STRLEN blen)
{
  int aneg, bneg;
  STRLEN i;

  /* a is checked, process b */
  bneg = strnum_parse(&b, &blen);

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

bool from_digit_string(UV* rn, const char* s, int base)
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

bool from_digit_to_UV(UV* rn, const UV* r, int len, int base)
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


bool from_digit_to_str(char** rstr, const UV* r, int len, int base)
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

#if BITS_PER_WORD == 64
  #define MAX_FIB_LEN 92
  #define MAX_FIB_STR "10100101000100000101000100010010001001000000001001000100100010101000100000101000101000001010"
#else
  #define MAX_FIB_LEN 46
  #define MAX_FIB_STR "1010001000010101000101000100000001000100100100"
#endif
#define MAX_FIB_VAL (MAX_FIB_LEN+1)

/* 0 = bad,   -1 = not canonical,   1 = good,   2 = ok but out of UV range */
int validate_zeckendorf(const char* str)
{
  int i;
  if (str == 0)
    return 0;
  if (str[0] != '1')
    return (str[0] == '0' && str[1] == '\0');
  /* str[0] = 1 */
  for (i = 1; str[i] != '\0'; i++) {
    if (str[i] == '1') {
      if (str[i-1] == '1')
        return -1;
    } else if (str[i] != '0') {
      return 0;
    }
  }
  /* Valid number.  Check if in range. */
  if (i > MAX_FIB_LEN || (i == MAX_FIB_LEN && strcmp(str, MAX_FIB_STR) > 0))
    return 2;
  return 1;
}

UV from_zeckendorf(const char* str)
{
  int i, len;
  UV n, fa = 0, fb = 1, fc = 1;  /* fc = fib(2) */

  if (str == 0) return 0;
  for (len = 0; len+1 <= MAX_FIB_LEN && str[len] != '\0'; len++)
    if (str[len] != '0' && str[len] != '1')
      return 0;
  if (len == 0 || len > MAX_FIB_LEN) return 0;
  n = (str[len-1] == '1');
  for (i = len-2; i >= 0; i--) {
    fa = fb; fb = fc; fc = fa+fb;  /* Advance */
    if (str[i] == '1') n += fc;
  }
  return n;
}

char* to_zeckendorf(UV n)
{
  char *str;
  int i, k, spos = 0;
  UV fa = 0, fb = 1, fc = 1;  /* fc = fib(2) */

  New(0, str, MAX_FIB_LEN+1, char);
  if (n == 0) {
    str[spos++] = '0';
  } else {
    UV rn = n;
    for (k = 2; k <= MAX_FIB_VAL && fc <= rn; k++) {
      fa = fb; fb = fc; fc = fa+fb;  /* Advance: fc = fib(k) */
    }
    for (i = k-1; i >= 2; i--) {
      fc = fb; fb = fa; fa = fc-fb;  /* Reverse: fc = fib(i) */
      str[spos++] = '0' + (fc <= rn);
      if (fc <= rn) rn -= fc;
    }
  }
  str[spos++] = '\0';
#if 0
  if (validate_zeckendorf(str) != 1) croak("to_zeckendorf bad for %lu\n",n);
  if (from_zeckendorf(str) != n) croak("to_zeckendorf wrong for %lu\n",n);
#endif
  return str;
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
bool is_catalan_pseudoprime(UV n) {
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
    divs = divisor_list(b2, &ndivisors, lim);
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

UV polygonal_root(UV n, UV k, bool* overflow) {
  UV D, R;
  MPUassert(k >= 3, "is_polygonal root < 3");
  *overflow = 0;
  if (n <= 1) return n;
  if (k == 4) {
    uint32_t root;
    return is_perfect_square_ret(n,&root)  ?  root  :  0;
  }
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

/*
 # On Mac M1.  The combinatorial solution that we use is both slower and
 # has *much* worse growth than the Rademacher implementation that uses high
 # precision floating point (e.g. Pari, MPFR, Arb).
 #
 #               10^5      10^6      10^7    10^8    10^9    10^10
 #   Perl-comb   78        ----
 #   GMP-comb     0.32     44        ----
 #   Sympy 1.7.1  0.0045    0.018    0.091    0.62     5.3     51
 #   Pari 2.14    0.00043   0.0018   0.013    0.19     4.5     54
 #   Bober 0.6    0.00010   0.00085  0.062    0.91    10.9     15
 #   Arb 2.19     0.00018   0.00044  0.004    0.011   0.031     0.086
 #
 #   Arb 2.19 takes only 62 seconds for 10^14.
*/

UV npartitions(UV n) {
  UV *part, *pent, i, j, k, d, npart;

  if (n <= 3)  return (n == 0) ? 1 : n;
  if (n > ((BITS_PER_WORD == 32) ? 127 : 416)) return 0;  /* Overflow */

  d = isqrt(n+1);
  New(0, pent, 2*d+2, UV);
  pent[0] = 0;
  pent[1] = 1;
  for (i = 1; i <= d; i++) {
    pent[2*i  ] = ( i   *(3*i+1)) / 2;
    pent[2*i+1] = ((i+1)*(3*i+2)) / 2;
  }
  New(0, part, n+1, UV);
  part[0] = 1;
  for (j = 1; j <= n; j++) {
    UV psum = 0;
    for (k = 1; pent[k] <= j; k++) {
      if ((k+1) & 2) psum += part[ j - pent[k] ];
      else           psum -= part[ j - pent[k] ];
    }
    part[j] = psum;
  }
  npart = part[n];
  Safefree(part);
  Safefree(pent);
  return npart;
}

UV consecutive_integer_lcm(UV n)
{
  UV i, ilcm, sqrtn;

  if (n <= 2)  return (n == 0) ? 1 : n;

  ilcm = 1;
  sqrtn = isqrt(n);
  for (i = 1; i < NPRIMES_TINY; i++) {
    uint32_t p = primes_tiny[i];
    if (p > n) break;
    if (p <= sqrtn) p = ipow(p, logint(n,p));
    if (ilcm > UV_MAX/p) return 0;
    ilcm *= p;
  }
  return ilcm;
}

UV frobenius_number(UV* A, uint32_t alen)
{
  UV g, i, j, max, *N, nlen;

  if (alen <= 1) return 0;
  sort_uv_array(A, alen);
  if (A[0] <= 1) return 0;

  for (g = A[0], i = 1; i < alen; i++)
    g = gcd_ui(g, A[i]);
  if (g != 1) croak("Frobenius number set must be coprime");

  if (UV_MAX/A[0] < A[1]) return UV_MAX;   /* Overflow */

  if (alen == 2)
    return A[0] * A[1] - A[0] - A[1];

  /* Algorithm "Round Robin" by Böcker and Lipták
   *
   * https://bio.informatik.uni-jena.de/wp/wp-content/uploads/2024/01/BoeckerLiptak_FastSimpleAlgorithm_reprint_2007.pdf
   *
   * This is the basic version, not the optimized one.  It's quite fast
   * in general, but the time is more or less O(A[0] * alen) and uses
   * A[0] * sizeof(UV) memory.  This means it's not going to work with very
   * large inputs, where something like DQQDU would work much better.
   *
   * See https://www.combinatorics.org/ojs/index.php/eljc/article/view/v12i1r27/pdf
   */

  nlen = A[0];
  /* if (nlen > 1000000000U) croak("overflow in frobenius number"); */
  New(0, N, nlen+1, UV);
  N[0] = 0;
  for (j = 1; j < nlen; j++)
    N[j] = UV_MAX;

  for (i = 1; i < alen; i++) {
    UV r, d, np, ai = A[i];
    np = N[ai % nlen];
    if (np != UV_MAX && np <= ai) continue;  /* Redundant basis (opt 3) */
    d = gcd_ui(nlen, ai);
    for (r = 0; r < d; r++) {
      UV p, q, n = 0;
      if (r > 0) {
        for (n = UV_MAX, q = r; q < nlen; q += d)
          if (N[q] < n)
            n = N[q];
      }
      if (n != UV_MAX) {
        for (j = 0; j < (nlen / d); j++) {
          n += ai;
          p = n % nlen;
          if (N[p] >= n)  N[p] = n;
          else            n = N[p];
        }
      }
    }
  }
  max = 0;
  for (i = 0; i < nlen; i++)
    if (N[i] == UV_MAX || (N[i] != UV_MAX && N[i] > max))
      max = N[i];
  Safefree(N);
  if (max == UV_MAX)  return UV_MAX;
  return max - nlen;
}


/* These rank/unrank are O(n^2) algorithms using O(n) in-place space.
 * Bonet 2008 gives O(n log n) algorithms using a bit more space.
 */

bool num_to_perm(UV k, int n, int *vec) {
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

bool perm_to_num(int n, int *vec, UV *rank) {
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
      sort_uv_array(S, k);
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

#define SMOOTH_TEST(n,k,p,nextprime) \
  if (n < p*p) return (n <= k);  /* p*p > n means n is prime */ \
  if ((n%p) == 0) { \
    do { n /= p; } while ((n%p) == 0); \
    if (n < nextprime) return 1; \
  } \
  if (k < nextprime) return (n <= k);

bool is_smooth(UV n, UV k) {
  UV fac[MPU_MAX_FACTORS+1];
  uint32_t i, p, pn;
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

  SMOOTH_TEST(n, k, 3,  5);  /* after this, k >=  5, n > 3*3 */
  SMOOTH_TEST(n, k, 5,  7);  /* after this, k >=  7, n > 5*5 */
  SMOOTH_TEST(n, k, 7, 11);  /* after this, k >= 11, n > 7*7 */

  /* Remove tiny factors.  Tests to 499. */
  for (i = 5, pn = primes_tiny[i]; i < NPRIMES_TINY-1; i++) {
    p = pn;  pn = primes_tiny[i+1];
    SMOOTH_TEST(n, k, p, pn);
  }
  if (k < pn || n < pn*pn) return (n <= k);   /* k >= 503 and n >= 503*503. */

  if (is_prime(n)) return 0;
  if (k <= 290000) {
    nfac = trial_factor(n, fac, pn, k);
    return (fac[nfac-1] <= k);
  }

  nfac = trial_factor(n, fac, pn, 4999);
  n = fac[nfac-1];
  pn = 5003;
  if (k < pn || n < pn*pn) return (n <= k);  /* k > 290k, n > 25M */

  /* Complete factoring including primality test */
  nfac = factor_exp(n, fac, 0);
  return (fac[nfac-1] <= k);
}
bool is_rough(UV n, UV k) {
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


static UV _divsum1(UV prod, UV f, uint32_t e) {
  UV pke, fmult;
  for (pke = f, fmult = 1+f; e > 1; e--) {
    pke *= f;
    fmult += pke;
  }
  return prod * fmult;
}

bool is_practical(UV n) {
  UV fac[MPU_MAX_FACTORS+1];
  UV exp[MPU_MAX_FACTORS+1];
  UV prod;
  int i, nfactors;

  if (n == 0 || (n & 1)) return (n == 1);
  if ((n & (n-1)) == 0) return 1;  /* All powers of 2 are practical */
  /* Allowable prefixes: {6,4} => {6,20,28,8} => {6,20,28,88,104,16} */
  if ((n % 6) && (n % 20) && (n % 28) && (n % 88) && (n % 104) && (n % 16))
    return 0;

  /* In theory for better performance we should test with small primes
   * before fully factoring.  On average it doesn't seem to help. */
  nfactors = factor_exp(n, fac, exp);
  MPUassert(fac[0] == 2, "is_practical first factor must be 2");
  prod = _divsum1(1, 2, exp[0]);
  for (i = 1; i < nfactors; i++) {
    if (fac[i] > (1 + prod))
      return 0;
    prod = _divsum1(prod, fac[i], exp[i]);
  }
  return 1;
}

int is_delicate_prime(UV n, uint32_t b) {

  if (b < 2) croak("is_delicate_prime base must be >= 2");
  if (b == 10 && n < 100)  return 0;  /* All 1,2,3,4 digit inputs are false */
  if (b ==  3 && n ==  2)  return 1;

  if (!is_prime(n)) return 0;

  if (b == 10) {

    UV d, dold, dnew, digpow, maxd = (BITS_PER_WORD == 32) ? 9 : 19;

    if (n >= ipow(10,maxd)) return -1;  /* We can't check all values */

    /* Check the last digit, given a > 1 digit prime, must be one of these. */
    dold = n % 10;
    if ( (dold != 1 && is_prime(n - dold + 1)) ||
         (dold != 3 && is_prime(n - dold + 3)) ||
         (dold != 7 && is_prime(n - dold + 7)) ||
         (dold != 9 && is_prime(n - dold + 9)) )
      return 0;

    /* Check the rest of the digits. */
    for (d = 1, digpow = 10;  d <= maxd && n >= digpow;  digpow *= 10, d++) {
      dold = (n / digpow) % 10;
      for (dnew = 0; dnew < 10; dnew++)
        if (dnew != dold && is_prime(n - dold*digpow + dnew*digpow))
          return 0;
    }

  } else if (b == 2) {

    UV bit;
    if (n < 127)  return 0;
    for (bit = log2floor(n); bit > 0; bit--)
      if (is_prime(n ^ (UVCONST(1) << bit)))
        return 0;

  } else {

#if 0   /* Our simpler method, but must add proper overflow check. */
    UV dold, dnew, digpow, N;
    for (digpow = 1;  n >= digpow;  digpow *= b) {
      dold = (n / digpow) % b;
      if ( (UV_MAX-(b-1)*digpow) < (n-dold*digpow) ) return -1;
      for (dnew = 0, N = n-dold*digpow;  dnew < b;  dnew++, N += digpow)
        if (dnew != dold && is_prime(N))
          return 0;
    }
#endif

    /* Algorithm isWeakly from Emily Stamm, 2020. */
    UV current, bm;
    for (bm = 1;  n >= bm;  bm *= b) {
      uint32_t j, counter;
      UV bmb = bm * b;
      if ( ((UV_MAX/b) < bm) || ((UV_MAX-bmb) < n) ) return -1; /* overflow */
      /* Check all n + j * b^m are composite */
      for (counter = 0, current = n+bm;
           (n % bm) != (current % bmb);
           counter++,  current += bm) {
        if (counter >= b-1) croak("is_delicate_prime overflow failure\n");
        if (is_prime(current))
          return 0;
      }
      /* Check all n - j * b^m are composite */
      for (j = 1, current = n-bm;  j < b-counter;  j++, current -= bm) {
        if (is_prime(current))
          return 0;
      }
    }

  }
  return 1;
}


bool is_sum_of_two_squares(UV n) {
  UV i, nfacs, fac[MPU_MAX_FACTORS+1], exp[MPU_MAX_FACTORS+1];
  if (n < 3) return 1;

  while (!(n&1)) n >>= 1;  /* Remove all factors of two */

  if ((n % 4) == 3) return 0;

  /* if (is_prime(n)) return ((n % 4) == 1); */

  /* TODO: a factor iterator should handle this reasonably */
  for (i = 0;  !(n %  3);  n /=  3) { i++; }    if ((i & 1) == 1)   return 0;
  for (i = 0;  !(n %  7);  n /=  7) { i++; }    if ((i & 1) == 1)   return 0;
  for (i = 0;  !(n % 11);  n /= 11) { i++; }    if ((i & 1) == 1)   return 0;
  for (i = 0;  !(n % 19);  n /= 19) { i++; }    if ((i & 1) == 1)   return 0;
  for (i = 0;  !(n % 23);  n /= 23) { i++; }    if ((i & 1) == 1)   return 0;
  for (i = 0;  !(n % 31);  n /= 31) { i++; }    if ((i & 1) == 1)   return 0;

  nfacs = factor_exp(n, fac, exp);
  for (i = 0; i < nfacs; i++)
    if ( (fac[i] % 4) == 3 && (exp[i] & 1) == 1 )
      return 0;

  return 1;
}

bool is_sum_of_three_squares(UV n) {
  UV tz = valuation(n,2);
  return ((tz & 1) == 1) || (((n>>tz) % 8) != 7);
}

#if 0  /* https://eprint.iacr.org/2023/807.pdf */
static UV halfgcd(UV m, UV u) {
  UV l = isqrt(m);
  UV a = m, b = u;
  while (a > l) {
    UV r = a % b;
    a = b;
    b = r;
  }
  return a;
}

/* Given an initial root, solve */
static bool corn_one(UV *x, UV *y, UV u, UV d, UV p) {
  UV rk = halfgcd(p, u);
  u = negmod(sqrmod(rk,p),p);
  u = (u % d == 0)  ?  u/d  :  0;
  if (u && is_perfect_square(u)) {
    *x = rk;
    *y = isqrt(u);
    return 1;
  }
  return 0;
}
#else
/* Given an initial root, solve.  Algorithm 2.3.12 of C&P */
static bool corn_one(UV *x, UV *y, UV u, UV d, UV p) {
  UV a = p;
  UV b = (2*u < p) ? p-u : u;
  uint32_t c = isqrt(p);
  while (b > c) {  UV t = a % b; a = b; b = t;  }
  u = p - b*b;
  u = (u % d == 0)  ?  u/d  :  0;
  if (u && is_perfect_square_ret(u,&c)) {
    *x = b;
    *y = c;
    return 1;
  }
  return 0;
}
#endif

  /* Cornacchia-Smith run over each root. */
static bool corn_all(UV *x, UV *y, UV d, UV p) {
  UV negd = negmod(d,p),  i, nroots, *roots;
  bool success = 0;
  roots = allsqrtmod(&nroots, negd, p);
  if (roots) {
    for (i = 0; i < nroots/2 && !success; i++)
      success = corn_one(x, y, roots[i], d, p);
    Safefree(roots);
  }
  return success;
}

bool cornacchia(UV *x, UV *y, UV d, UV p) {
  UV u, negd, limu;
  uint32_t root;

  if (p == 0) {
    *x = *y = 0;
    return 1;
  }

  if (d == 0) {
    if (!is_perfect_square_ret(p,&root))  return 0;
    *x = root;  *y = 0;
    return 1;
  }

  negd = negmod(d, p);

  if (is_prime(p)) {
    if (kronecker_uu(negd,p) == -1) return 0;
    if (!sqrtmodp(&u, negd, p))     return 0;
    return corn_one(x, y, u, d, p);
  }

  if (((p >> 31) >> 22) && kronecker_uu(negd,p) != -1 && corn_all(x, y, d, p))
    return 1;

  /* Loop through all valid integers until one is found.
   * Until p is quite large, this is faster than using allsqrtmod.
   * It also finds many solutions for composites.
   */
  for (u = 0, limu = isqrt(p/d);  u <= limu;  u++) {
    UV t = p - d*u*u;
    if (is_perfect_square_ret(t,&root)) {
      *x = root;  *y = u;  return 1;
    }
  }

  return 0;
}




/* TODO: */
/* See https://arxiv.org/pdf/2208.01725.pdf for new info on smooth count
 * estimate methods.  Consider adding an estimate function. */

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
  return legendre_phi(x, prime_count(y-1));
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
