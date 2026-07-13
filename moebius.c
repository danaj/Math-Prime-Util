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

static UV _moebius_sieve_limit(UV hi)
{
  UV sqrtn = isqrt(hi);
  if (sqrtn*sqrtn != hi && sqrtn < (UVCONST(1)<<(BITS_PER_WORD/2))-1)
    sqrtn++;
  return sqrtn;
}

static INLINE void _sieve_moebius_prime(signed char* mu, UV lo, UV hi, UV p,
                                        unsigned char* logp, UV* nextlog)
{
  UV i, p2 = p*p;
  if (p > *nextlog) {
    *logp += 2;   /* logp is 1 | ceil(log(p)/log(2)) */
    *nextlog = ((*nextlog-1)*4)+1;
  }
  for (i = P_GT_LO(p, p, lo); i >= lo && i <= hi; i += p)
    mu[i-lo] += *logp;
  for (i = P_GT_LO(p2, p2, lo); i >= lo && i <= hi; i += p2)
    mu[i-lo] = (signed char)0x80;
}

static void _sieve_moebius(const unsigned char* sieve, signed char* mu,
                           UV lo, UV hi)
{
  UV i, sqrtn = _moebius_sieve_limit(hi), count = hi-lo+1;
  unsigned char logp;
  UV nextlog, nextlogi;

  memset(mu, 0, count * sizeof(signed char));

  logp = 1; nextlog = 3; /* 2+1 */
  if (sqrtn >= 2) _sieve_moebius_prime(mu,lo,hi,2,&logp,&nextlog);
  if (sqrtn >= 3) _sieve_moebius_prime(mu,lo,hi,3,&logp,&nextlog);
  if (sqrtn >= 5) _sieve_moebius_prime(mu,lo,hi,5,&logp,&nextlog);
  if (sqrtn >= 7) {
    START_DO_FOR_EACH_SIEVE_PRIME(sieve, 0, 7, sqrtn) {
      _sieve_moebius_prime(mu,lo,hi,p,&logp,&nextlog);
    } END_DO_FOR_EACH_SIEVE_PRIME
  }

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
}

static bool _moebius_by_hand(UV lo, UV hi)
{
  UV width = hi-lo;
  return hi < 100 || width < 10 ||
         (hi > (UVCONST(1)<<25) && width < icbrt(hi)/4);
}

/* Return a char array with lo-hi+1 elements. mu[k-lo] = µ(k) for k = lo .. hi.
 * It is the callers responsibility to call Safefree on the result. */
signed char* range_moebius(UV lo, UV hi)
{
  signed char* mu;
  UV i, count;

  if (hi < lo) croak("range_moebius error hi %"UVuf" < lo %"UVuf"\n", hi, lo);
  if (hi-lo == UV_MAX) croak("range_moebius range too large\n");
  count = hi-lo+1;
  New(0, mu, count, signed char);

  /* Small high-offset ranges are faster by individual factorization. */
  if (_moebius_by_hand(lo, hi)) {
    for (i = 0; i < count; i++)
      mu[i] = (signed char)moebius(lo+i);
  } else {
    const unsigned char* sieve;
    get_prime_cache(_moebius_sieve_limit(hi), &sieve);
    _sieve_moebius(sieve, mu, lo, hi);
    release_prime_cache(sieve);
  }

  return mu;
}

#define MOEBIUS_SEGMENT_SIZE UVCONST(32*1024)

typedef struct {
  UV lo;
  UV hi;
  UV segment_size;
  bool done;
  bool by_hand;
  signed char* segment;
  const unsigned char* sieve;
} moebius_segment_context_t;

void* start_segment_moebius(UV lo, UV hi, signed char** segmentmem)
{
  moebius_segment_context_t* ctx;
  UV width;

  New(0, ctx, 1, moebius_segment_context_t);
  ctx->lo = lo;
  ctx->hi = hi;
  ctx->done = hi < lo;
  ctx->by_hand = ctx->done || _moebius_by_hand(lo, hi);
  ctx->sieve = 0;

  width = ctx->done ? 0 : hi-lo;
  ctx->segment_size = (width < MOEBIUS_SEGMENT_SIZE)
                    ? width+1 : MOEBIUS_SEGMENT_SIZE;
  if (ctx->segment_size == 0) ctx->segment_size = 1;
  New(0, ctx->segment, ctx->segment_size, signed char);

  if (!ctx->by_hand)
    get_prime_cache(_moebius_sieve_limit(hi), &ctx->sieve);

  *segmentmem = ctx->segment;
  return ctx;
}

bool next_segment_moebius(void* vctx, UV* lo, UV* hi)
{
  moebius_segment_context_t* ctx = vctx;
  UV seglo, seghi, i, count;

  if (ctx->done) return 0;
  seglo = ctx->lo;
  seghi = (ctx->hi-seglo >= ctx->segment_size)
        ? seglo + ctx->segment_size - 1 : ctx->hi;
  count = seghi-seglo+1;

  if (ctx->by_hand) {
    for (i = 0; i < count; i++)
      ctx->segment[i] = (signed char)moebius(seglo+i);
  } else {
    _sieve_moebius(ctx->sieve, ctx->segment, seglo, seghi);
  }

  *lo = seglo;
  *hi = seghi;
  if (seghi == ctx->hi) ctx->done = 1;
  else                  ctx->lo = seghi+1;
  return 1;
}

void end_segment_moebius(void* vctx)
{
  moebius_segment_context_t* ctx = vctx;
  if (ctx->sieve != 0) release_prime_cache(ctx->sieve);
  Safefree(ctx->segment);
  Safefree(ctx);
}

static short* mertens_array(UV hi)
{
  signed char* mu;
  short* M;
  void* mctx;
  UV i, lo, seghi;
  IV sum = 0;

  New(0, M, hi+1, short);
  mctx = start_segment_moebius(0, hi, &mu);
  while (next_segment_moebius(mctx, &lo, &seghi)) {
    for (i = lo; i <= seghi; i++) {
      sum += mu[i-lo];
      M[i] = (short)sum;
    }
  }
  end_segment_moebius(mctx);

  return M;
}

typedef struct {
  UV n;
  IV sum;
} mertens_value_t;
#define MERTENS_HASH_PROBES 8
static void _insert_mert_hash(mertens_value_t *H, UV hsize, UV n, IV sum) {
  UV i, idx = n % hsize, base = idx;
  for (i = 0; i < MERTENS_HASH_PROBES && H[idx].n != 0 && H[idx].n != n; i++)
    if (++idx == hsize) idx = 0;
  if (i == MERTENS_HASH_PROBES) idx = base;
  H[idx].n = n;
  H[idx].sum = sum;
}
static int _get_mert_hash(mertens_value_t *H, UV hsize, UV n, IV *sum) {
  UV i, idx = n % hsize;
  for (i = 0; i < MERTENS_HASH_PROBES && H[idx].n != 0; i++) {
    if (H[idx].n == n) {
      *sum = H[idx].sum;
      return 1;
    }
    if (++idx == hsize) idx = 0;
  }
  return 0;
}

typedef struct {
  UV maxmu;
  short* M;             /* 16 bits covers M(n) through n = 7,613,644,883 */
  UV hsize;
  mertens_value_t* H;
} mertens_context_t;

static IV _mertens_uv_to_iv(UV n)
{
  return (n <= (UV)IV_MAX) ? (IV)n : -1 - (IV)(UV_MAX-n);
}

/* Thanks to Trizen for this algorithm. */
static IV _rmertens(mertens_context_t* ctx, UV n) {
  const UV maxmu = ctx->maxmu;
  const short* M = ctx->M;
  UV s, k, ns, nk, nk1;
  UV sum;
  IV cached, mnk;

  if (n <= maxmu)
    return M[n];

  if (_get_mert_hash(ctx->H, ctx->hsize, n, &cached))
    return cached;

  s = isqrt(n);
  ns = n / (s+1);

  /* Merge the two defining sums and reuse each quotient division. */
  if (s != ns && s != ns+1) croak("mertens  s / ns");
  nk  = n;
  nk1 = n/2;
  sum = 1 - nk1 - (nk & 1);
  for (k = 2; k <= ns; k++) {
    nk = nk1;
    nk1 = n/(k+1);
    mnk = (nk <= maxmu) ? M[nk] : _rmertens(ctx, nk);
    sum -= (UV)mnk + (UV)(IV)M[k] * (nk-nk1);
  }
  if (s > ns)
    sum -= (UV)(IV)M[s] * (n/s - n/(s+1));

  cached = _mertens_uv_to_iv(sum);
  _insert_mert_hash(ctx->H, ctx->hsize, n, cached);
  return cached;
}

static mertens_context_t* _mertens_context_create(UV n) {
  mertens_context_t* ctx;
  UV j = icbrt(n);
  UV maxmu = 1 * j * j;
  UV hsize = next_prime(100 + 8*j);
  UV sqrtn = isqrt(n);

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

  /* _rmertens uses direct M[k] and M[s] lookups through sqrt(n). */
  if (maxmu < sqrtn) maxmu = sqrtn;

  New(0, ctx, 1, mertens_context_t);
  ctx->maxmu = maxmu;
  ctx->M = mertens_array(maxmu);
  ctx->hsize = hsize;
  Newz(0, ctx->H, hsize, mertens_value_t);
  return ctx;
}

static void _mertens_context_destroy(mertens_context_t* ctx)
{
  Safefree(ctx->H);
  Safefree(ctx->M);
  Safefree(ctx);
}

static IV _small_mertens(UV n)
{
  static const signed char MV16[33] = {0,-1,-4,-3,-1,-4,2,-4,-2,-1,0,-4,-5,-3,3,-1,-1,-3,-7,-2,-4,2,1,-1,-2,1,1,-3,-6,-6,-6,-5,-4};
  UV j = n/16;
  IV sum = MV16[j];
  for (j = j*16 + 1; j <= n; j++)
    sum += moebius(j);
  return sum;
}

static IV _segmented_mertens(UV lo, UV hi)
{
  signed char* mu;
  void* mctx;
  UV seglo, seghi, i, count;
  IV sum = 0;

  mctx = start_segment_moebius(lo, hi, &mu);
  while (next_segment_moebius(mctx, &seglo, &seghi)) {
    count = seghi-seglo+1;
    for (i = 0; i < count; i++)
      sum += mu[i];
  }
  end_segment_moebius(mctx);
  return sum;
}

IV mertens_range(UV lo, UV hi) {
  mertens_context_t* ctx;
  UV j, width;
  IV sum = 0;

  if (lo < 1) lo = 1;
  if (hi < lo) return 0;
  if (hi <= 512)
    return _small_mertens(hi) - ((lo <= 1) ? 0 : _small_mertens(lo-1));

  /* Endpoint Mertens builds about hi^(2/3) prefix values.  Use half that
   * as a conservative crossover after accounting for segmented marking. */
  width = hi-lo+1;
  j = icbrt(hi);
  if (width <= (j*j)/2)
    return _segmented_mertens(lo, hi);

  ctx = _mertens_context_create(hi);
  sum = _rmertens(ctx, hi);
  if (lo > 1) sum -= _rmertens(ctx, lo-1);
  _mertens_context_destroy(ctx);
  return sum;
}

IV mertens(UV n) {
  return mertens_range(1, n);
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
  mertens_context_t* ctx;
  UV j, k, nk, sqrtn;
  IV sum;

  if (n <= 96) {
    signed char* l = liouville_array(n);
    for (sum = 0, j = 1; j <= n; j++)
      sum += l[j];
    Safefree(l);
    return sum;
  }

  ctx = _mertens_context_create(n);

  sqrtn = isqrt(n);
  sum = _rmertens(ctx, n);
  for (k = 2; k <= sqrtn; k++) {
    nk = n / (k*k);
    if (nk == 1) break;
    sum += _rmertens(ctx, nk);
  }
  sum += (sqrtn + 1 - k);  /* all k where n/(k*k) == 1 */
  /* TODO: find method to get exact number of n/(k*k)==1 .. 4.  Halves k */
  /*       Ends up with method like Lehmer's g. */

  _mertens_context_destroy(ctx);
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
