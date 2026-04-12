#include <string.h>

#define FUNC_isqrt 1
#define FUNC_log2floor 1
#include "ptypes.h"
#include "cache.h"
#include "sieve.h"
#include "util.h"
#include "factor.h"
#include "moebius.h"


#define MOB_TESTP(p) \
  { uint32_t psq = p*p;  if (n >= psq && (n % psq) == 0) return 0; }

/* mpu 'for (0..255) { $x=moebius($_)+1; $b[$_ >> 4] |= ($x << (2*($_%16))); } say join ",",@b;' */
static const uint32_t _smoebius[16] = {2703565065U,23406865,620863913,1630114197,157354249,2844895525U,2166423889U,363177345,2835441929U,2709852521U,1095049497,92897577,1772687649,162113833,160497957,689538385};
int moebius(UV n) {
  if (n < 256)  return (int)((_smoebius[n >> 4] >> (2*(n % 16))) & 3) - 1;

  if (!(n % 4) || !(n % 9) || !(n % 25) || !(n % 49) || !(n %121) || !(n %169))
    return 0;

  MOB_TESTP(17); MOB_TESTP(19); MOB_TESTP(23);
  MOB_TESTP(29); MOB_TESTP(31); MOB_TESTP(37);

  return factored_moebius(factorint(n));
}

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

  if (hi < lo) croak("range_moebius error hi %"UVuf" < lo %"UVuf"\n", hi, lo);

  Newz(0, mu, count, signed char);
  if (sqrtn*sqrtn != hi && sqrtn < (UVCONST(1)<<(BITS_PER_WORD/2))-1) sqrtn++;

  /* For small ranges, do it by hand */
  if (hi < 100 || count <= 10 || (hi > (1UL<<25) && count < icbrt(hi)/4)) {
    for (i = 0; i < count; i++)
      mu[i] = (signed char)moebius(lo+i);
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
      mu[i-lo] = (signed char)0x80;
  } END_DO_FOR_EACH_PRIME

  logp = (unsigned char)log2floor(lo);
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
    double rfactor = 1.0 + 7.0 * (1.0 - exp(-(double)maxmu/8000000000.0));
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
    UV twoa1 = a > (UV_MAX >> 1)  ?  UV_MAX  :  2*a-1;
    b = twoa1 < hi  ?  twoa1  :  hi;
    START_DO_FOR_EACH_PRIME(2, isqrt(b)) {
      for (k = 2*p; k <= b; k += p) {
        if (k >= a)
          l[k] = -1 * l[k/p];
      }
    } END_DO_FOR_EACH_PRIME
    if (b == hi) break;
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

  if (hi < lo) croak("range_liouville error hi %"UVuf" < lo %"UVuf"\n",hi,lo);
  nf = range_nfactor_sieve(lo, hi, 1);
  New(0, l, hi-lo+1, signed char);
  for (i = 0; i < hi-lo+1; i++)
    l[i] = (nf[i] & 1) ? -1 : 1;
  Safefree(nf);
  return l;
}
