#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "constants.h"
#define FUNC_isqrt 1
#define FUNC_icbrt 1
#define FUNC_ctz 1
#include "util.h"
#include "cache.h"
#include "sieve.h"
#include "lmo.h"
#include "inverse_interpolate.h"
#include "semi_primes.h"

#define SP_SIEVE_THRESH 100  /* When to sieve vs. iterate */

/******************************************************************************/
/*                                SEMI PRIMES                                 */
/******************************************************************************/

#if 0
static const unsigned char _semiprimelist[] =
  {0,4,6,9,10,14,15,21,22,25,26,33,34,35,38,39,46,49,51,55,57,58,62,65,69,74,
   77,82,85,86,87,91,93,94,95,106,111,115,118,119,121,122,123,129,133,134,141,
   142,143,145,146,155,158,159,161,166,169,177,178,183,185,187,194,201,202,
   203,205,206,209,213,214,215,217,218,219,221,226,235,237,247,249,253,254};
#else
static const unsigned short _semiprimelist[] =
  {0,4,6,9,10,14,15,21,22,25,26,33,34,35,38,39,46,49,51,55,57,58,62,65,69,74,
   77,82,85,86,87,91,93,94,95,106,111,115,118,119,121,122,123,129,133,134,141,
   142,143,145,146,155,158,159,161,166,169,177,178,183,185,187,194,201,202,
   203,205,206,209,213,214,215,217,218,219,221,226,235,237,247,249,253,254,
   259,262,265,267,274,278,287,289,291,295,298,299,301,302,303,305,309,314,
   319,321,323,326,327,329,334,335,339,341,346,355,358,361,362,365,371,377,
   381,382,386,391,393,394,395,398,403,407,411,413,415,417,422,427,437,445,
   446,447,451,453,454,458,466,469,471,473,478,481,482,485,489,493,497,501,
   502,505,511,514,515,517,519,526,527,529,533,535,537,538,542,543,545,551,
   553,554,559,562,565,566,573,579,581,583,586,589,591,597,611,614,622,623};
#endif
#define NSEMIPRIMELIST (sizeof(_semiprimelist)/sizeof(_semiprimelist[0]))

#if 1
static UV _bs_count(UV n, UV const* const primes, UV lastidx)
{
  UV i = 0, j = lastidx;   /* primes may not start at 0 */
  MPUassert(n >= primes[0] && n < primes[lastidx], "prime count via binary search out of range");
  while (i < j) {
    UV mid = i + (j-i)/2;
    if (primes[mid] <= n)  i = mid+1;
    else                   j = mid;
  }
  return i-1;
}

static UV _semiprime_count(UV n)
{
  UV pc = 0, sum = 0, sqrtn = prev_prime(isqrt(n)+1);
  UV xbeg = 0, xend = 0, xlim = 0, xoff = 0, xsize, *xarr = 0;
  UV const xmax = 200000000UL;

  if (n > 1000000) { /* Upfront work to speed up the many small calls */
    UV nprecalc = (UV) pow(n, .75);
    if (nprecalc > _MPU_LMO_CROSSOVER)  nprecalc = _MPU_LMO_CROSSOVER;
    prime_precalc(nprecalc);
    /* Make small calls even faster using binary search on a list */
    xlim = (UV) pow(n, 0.70);
  }

  if (sqrtn >= 2)  sum += LMO_prime_count(n/2) - pc++;
  if (sqrtn >= 3)  sum += LMO_prime_count(n/3) - pc++;
  if (sqrtn >= 5)  sum += LMO_prime_count(n/5) - pc++;
  if (sqrtn >= 7) {
    unsigned char* segment;
    UV seg_base, seg_low, seg_high, np, cnt;
    void* ctx = start_segment_primes(7, sqrtn, &segment);
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
        np = n/p;
        if (np < xlim) {
          if (xarr == 0 || np < xbeg) {
            if (xarr != 0) { Safefree(xarr); xarr = 0; }
            xend = np;
            xbeg = n/sqrtn;
            if (xend - xbeg > xmax) xbeg = xend - xmax;
            xbeg = prev_prime(xbeg);
            xend = next_prime(xend);
            xoff = LMO_prime_count(xbeg);
            xsize = range_prime_sieve(&xarr, xbeg, xend);
            xend = xarr[xsize-1];
          }
          cnt = xoff + _bs_count(np, xarr, xsize-1);
        } else {
          cnt = LMO_prime_count(np);
        }
        sum += cnt - pc++;
      END_DO_FOR_EACH_SIEVE_PRIME
    }
    if (xarr != 0) { Safefree(xarr); xarr = 0; }
    end_segment_primes(ctx);
  }
  return sum;
}
#else

/* This is much cleaner, but ends up being a little slower. */

#include "prime_count_cache.h"
#define CACHED_PC(cache,n) prime_count_cache_lookup(cache,n)

static UV _semiprime_count(UV n)
{
  UV sum = 0, sqrtn = prev_prime(isqrt(n)+1), pc_sqrtn;
  void *cache = prime_count_cache_create( (UV)pow(n,0.70) );

  if (sqrtn >= 2)  sum += CACHED_PC(cache,n/2);
  if (sqrtn >= 3)  sum += CACHED_PC(cache,n/3);
  if (sqrtn >= 5)  sum += CACHED_PC(cache,n/5);
  if (sqrtn >= 7) {
    unsigned char* segment;
    UV seg_base, seg_low, seg_high;
    void* ctx = start_segment_primes(7, sqrtn, &segment);
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
        sum += CACHED_PC(cache, n/p);
      END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
  }
  pc_sqrtn = CACHED_PC(cache, sqrtn);
  sum -= (pc_sqrtn * pc_sqrtn - pc_sqrtn) / 2;
  prime_count_cache_destroy(cache);
  return sum;
}
#endif

/* TODO: This overflows, see p=3037000507,lo=10739422018595509581.
 *       p2 = 9223372079518257049 => 9223372079518257049 + 9223372079518257049
 *       Also with lo=18446744073709551215,hi=18446744073709551515.
 *       Using P_GT_LO_0 might help, but the biggest issue is 2*p*p overflows.
 */
#define MARKSEMI(p,arr,lo,hi) \
    do {  UV i, p2=(p)*(p); \
      for (i = P_GT_LO(p2, p, lo); i >= lo && i <= hi; i += p) arr[i-lo]++; \
      for (i = P_GT_LO(2*p2, p2, lo); i >= lo && i <= hi; i += p2) arr[i-lo]++; \
    } while (0);

UV range_semiprime_sieve(UV** semis, UV lo, UV hi)
{
  UV *S, i, count = 0;

  if (lo < 4) lo = 4;
  if (hi > MPU_MAX_SEMI_PRIME) hi = MPU_MAX_SEMI_PRIME;

  if (hi <= _semiprimelist[NSEMIPRIMELIST-1]) {
    if (semis == 0) {
      for (i = 1; i < NSEMIPRIMELIST && _semiprimelist[i] <= hi; i++)
        if (_semiprimelist[i] >= lo)
          count++;
    } else {
      Newz(0, S, NSEMIPRIMELIST+1, UV);
      for (i = 1; i < NSEMIPRIMELIST && _semiprimelist[i] <= hi; i++)
        if (_semiprimelist[i] >= lo)
          S[count++] = _semiprimelist[i];
      *semis = S;
    }
  } else {
    unsigned char* nfacs;
    UV cutn, sqrtn = isqrt(hi);
    Newz(0, nfacs, hi-lo+1, unsigned char);
    if (sqrtn*sqrtn < hi && sqrtn < (UVCONST(1)<<(BITS_PER_WORD/2))-1) sqrtn++;
    cutn = (sqrtn > 30000) ? 30000 : sqrtn;
    START_DO_FOR_EACH_PRIME(2, cutn) {
      MARKSEMI(p,nfacs,lo,hi);
    } END_DO_FOR_EACH_PRIME
    if (cutn < sqrtn) {
      unsigned char* segment;
      UV seg_base, seg_low, seg_high;
      void* ctx = start_segment_primes(cutn, sqrtn, &segment);
      while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
        START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
          MARKSEMI(p,nfacs,lo,hi);
        END_DO_FOR_EACH_SIEVE_PRIME
      }
      end_segment_primes(ctx);
    }
    if (semis == 0) {
      for (i = lo; i <= hi; i++)
        if (nfacs[i-lo] == 1)
          count++;
    } else {
      UV cn = 50 + 1.01 * (semiprime_count_approx(hi) - semiprime_count_approx(lo));
      New(0, S, cn, UV);
      for (i = lo; i <= hi; i++) {
        if (nfacs[i-lo] == 1) {
          if (count >= cn)
            Renew(S, cn += 4000, UV);
          S[count++] = i;
        }
      }
      *semis = S;
    }
    Safefree(nfacs);
  }
  return count;
}

static UV _range_semiprime_count_iterate(UV lo, UV hi)
{
  UV sum = 0;
  for (; lo < hi; lo++)     /* TODO: We should walk composites */
    if (is_semiprime(lo))
      sum++;
  if (is_semiprime(hi))
    sum++;
  return sum;
}

#if 0
static UV _range_semiprime_selection(UV** semis, UV lo, UV hi)
{
  UV *S = 0, *pr, cn = 0, count = 0;
  UV i, xsize, lim = hi/2 + 1000, sqrtn = isqrt(hi);

  if (lo < 4) lo = 4;
  if (hi > MPU_MAX_SEMI_PRIME) hi = MPU_MAX_SEMI_PRIME;

  if (semis != 0) {
    cn = 50 + 1.01 * (semiprime_count_approx(hi) - semiprime_count_approx(lo));
    New(0, S, cn, UV);
  }

  xsize = range_prime_sieve(&pr, 0, lim);

  for (i = 0; pr[i] <= sqrtn; i++) {
    UV const pi = pr[i], jlo = (lo+pi-1)/pi, jhi = hi/pi;
    UV skip, j = i;
    if (pr[j] < jlo)
     for (skip = 2048; skip > 0; skip >>= 1)
       while (j+skip-1 < xsize && pr[j+skip-1] < jlo)
         j += skip;
    if (semis == 0) {
      while (pr[j++] <= jhi)
        count++;
    } else {
      for (; pr[j] <= jhi; j++) {
        if (count >= cn)
          Renew(S, cn += 4000, UV);
        S[count++] = pi * pr[j];
      }
    }
  }
  Safefree(pr);
  if (semis != 0) {
    qsort(S, count, sizeof(UV), _numcmp);
    *semis = S;
  }
  return count;
}
#endif



UV semiprime_count(UV lo, UV hi)
{
  if (lo > hi || hi < 4)
    return 0;

  /* tiny sizes fastest with the sieving code */
  if (hi <= 400) return range_semiprime_sieve(0, lo, hi);
  /* Large sizes best with the prime count method */
  if (lo <= 4) return _semiprime_count(hi);

  /* Now it gets interesting.  lo > 4, hi > 400. */

  if ((hi-lo+1) < hi / (isqrt(hi)*200)) {
    MPUverbose(2, "semiprimes %"UVuf"-%"UVuf" via iteration\n", lo, hi);
    return _range_semiprime_count_iterate(lo,hi);
  }
  /* TODO: Determine when _range_semiprime_selection(0,lo,hi) is better */
  if ((hi-lo+1) < hi / (isqrt(hi)/4)) {
    MPUverbose(2, "semiprimes %"UVuf"-%"UVuf" via sieving\n", lo, hi);
    return range_semiprime_sieve(0, lo, hi);
  }
  MPUverbose(2, "semiprimes %"UVuf"-%"UVuf" via prime count\n", lo, hi);
  return _semiprime_count(hi) - _semiprime_count(lo-1);
}

static UV _cb_sca(UV mid, UV k) { return semiprime_count_approx(mid); }

UV semiprime_count_approx(UV n) {
  UV i;

  if (n <= _semiprimelist[NSEMIPRIMELIST-1]) {

    for (i = 0;  i < NSEMIPRIMELIST-1 && n >= _semiprimelist[i+1];  i++)
      ;
    return i;

  } else {

    /* Crisan and Erban (2020)  https://arxiv.org/abs/2006.16491 */
    UV L, res;
    double logn = log(n), loglogn = log(logn);
    double series = 0, den = 1, mc;
    static const double C[19] = {
         0.26149721284764278375L,
        -2.0710850628855780875L,
        -7.6972777412176108802L,
        -35.345660320564161516L,
        -206.71503925406509339L,
        -1511.1997871316530251L,
        -13546.323682845914021L,
        -146229.10675883565523L,
        -1867579.6280076650637L,
        -27733045.258413542557L,
        -470983423.57703294361L,
 /*
  * Values for C_11+ are not exact, but that's ok here.
  * \p 80
  * zetald(n) = { zeta'(n) / zeta(n) }
  * zetalim(n) = { derivnum(s = 1-1e-40, zetald(s) + 1/(s-1), n-1) }
  * B(n,x=100) = { if(n==0,return(0.2614972128476427837554268386086958590516)); (-1)^n * (sum(i=2, x, moebius(i) * i^(n-1) * derivnum(X=i,zetald(X),n-1)) + zetalim(n)) }
  * BN = vector(20,n,B(n-1,500));
  * C(n) = { n!*(sum(i=0,n,BN[i+1]/i!) - sum(i=1,n,1/i)) }
  */
        -9011500984.0L,
        -191744069149.0L,
        -4487573459711.0L,
        -114472069580580.0L,
        -3158610502077136.0L,
        -93682567786528912.0L,
        -2970838770257639696.0L,
        -100274471240063911725.0L };   /* ~ C_18 */
    static const double CROSS[16] =
     {     632, 9385, 136411, 4610076, 66358000, 440590000, 2557200000, 53032001000.0, 1151076796431.0, 20416501389724.0, 165815501587300.0, 2316269995165180.0, 12857302245660782.0, 82444629251749248.0, 516965078559580064.0, 1875312882821895168.0 };
     /*                                                          ^^ last value is not right */
    static const double mincount[17] =
     { 82, 195, 2485,  31446,  906319, 11741185,  72840337,  398702652,    7538564737,  150382042176.0,  2482510001499.0,  19204997230933.0,  253305676016061.0,  1357320342479135.0,  8390553919247938.0,  50820734885440616.0,  180076215506140384.0 };

    /* Pick truncation point */
    for (L = 3; L <= 18 && (double)n >= CROSS[L-3]; L++)  ;
    /* TODO: Tune crossover for 18 */

    /* Calculate truncated asymptotic value */
    for (i = 1; i <= L; i++) {
      series += factorial(i-1) * (loglogn / den);
      series += C[i-1] / den;
      den *= logn;
    }
    res = (UV) ( (n / logn) * series + 0.5L );

    /* Check for overflow */
    if (res >= MPU_MAX_SEMI_PRIME_IDX)  return MPU_MAX_SEMI_PRIME_IDX;

    /* Ensure monotonic using simple clamping */
    mc = mincount[L-3];
    /* mc = (L == 3) ? 82 : semiprime_count_approx(CROSS[L-4]-1); */
    if ((double)res < mc)  return mc;

    return res;

  }
}

UV nth_semiprime_approx(UV n) {
  double logn,log2n,log3n,log4n, err_lo, err_md, err_hi, err_factor, est;
  UV lo, hi;

  if (n < NSEMIPRIMELIST)
    return _semiprimelist[n];
  if (n > MPU_MAX_SEMI_PRIME_IDX)
    return 0;

  /* Piecewise with blending.  Hacky and maybe overkill.  It makes a good
   * estimator by itself, but our count approximation is even better, so we
   * use this as an excellent initial estimate, then use inverse binary
   * search to lower the error another order of magnitude.
   *
   *     Interp Range    Crossover to next
   * lo   2^8  - 2^28      2^26 - 2^27
   * md   2^25 - 2^48      2^46 - 2^47
   * hi   2^45 - 2^64
   */

  logn = log(n);   log2n = log(logn);   log3n = log(log2n);   log4n=log(log3n);
  err_lo = 1.000 - 0.00018216088*logn + 0.18099609886*log2n - 0.51962474356*log3n - 0.01136143381*log4n;
  err_md = 0.968 - 0.00073297945*logn + 0.09731690314*log2n - 0.25212500749*log3n - 0.01366795346*log4n;
  err_hi = 0.968 - 0.00008034109*logn + 0.01522628393*log2n - 0.04020257367*log3n - 0.01266447175*log4n;

  if        (n <= (1UL<<26)) {
    err_factor = err_lo;
  } else if (n < (1UL<<27)) { /* Linear interpolate the two in the blend area */
    double x = (n - 67108864.0L) / 67108864.0L;
    err_factor = ((1.0L-x) * err_lo) + (x * err_md);
  } else if (logn <= 31.88477030575) {
    err_factor = err_md;
  } else if (logn <  32.57791748632) {
    double x = (n - 70368744177664.0L) / 70368744177664.0L;
    err_factor = ((1.0L-x) * err_md) + (x * err_hi);
  } else {
    err_factor = err_hi;
  }
  est = err_factor * n * logn / log2n;
  if (est >= MPU_MAX_SEMI_PRIME) return MPU_MAX_SEMI_PRIME;

  /* Use inverse interpolation to improve the result. */
  lo = 0.98 * est - 5;
  hi = 1.03 * est;
  return inverse_interpolate(lo, hi, n, 0, &_cb_sca, 0);
}

static UV _next_semiprime(UV n) {
  while (!is_semiprime(++n))
    ;
  return n;
}
static UV _prev_semiprime(UV n) {
  while (!is_semiprime(--n))
    ;
  return n;
}
UV nth_semiprime(UV n)
{
  UV guess, spcnt, sptol, gn, ming = 0, maxg = UV_MAX;

  if (n < NSEMIPRIMELIST)
    return _semiprimelist[n];

  guess = nth_semiprime_approx(n);    /* Initial guess */
  sptol = 16*icbrt(n);                /* Guess until within this many SPs */
  MPUverbose(2, "  using exact counts until within %"UVuf"\n",sptol);

  /* Make successive interpolations until small enough difference */
  for (gn = 2; gn < 20; gn++) {
    IV adjust;
    while (!is_semiprime(guess)) guess++;  /* Guess is a semiprime */
    MPUverbose(2, "  %"UVuf"-th semiprime is around %"UVuf" ... ", n, guess);
    /* Compute exact count at our nth-semiprime guess */
    spcnt = _semiprime_count(guess);
    MPUverbose(2, "(%"IVdf")\n", (IV)(n-spcnt));
    /* Stop guessing if within our tolerance */
    if (n==spcnt || (n>spcnt && n-spcnt < sptol) || (n<spcnt && spcnt-n < sptol)) break;
    /* Determine how far off we think we are */
    adjust = (IV) (nth_semiprime_approx(n) - nth_semiprime_approx(spcnt));
    /* When computing new guess, ensure we don't overshoot.  Rarely used. */
    if (spcnt <= n && guess > ming) ming = guess;   /* Previous guesses */
    if (spcnt >= n && guess < maxg) maxg = guess;
    guess += adjust;
    if (guess <= ming || guess >= maxg) MPUverbose(2, "  fix min/max for %"UVuf"\n",n);
    if (guess <= ming) guess = ming + sptol - 1;
    if (guess >= maxg) guess = maxg - sptol + 1;
  }

  /* If we have far enough to go, sieve for semiprimes */
  if (n > spcnt && (n-spcnt) > SP_SIEVE_THRESH) {         /* sieve forwards */
    UV *S, count, i, range;
    while (n > spcnt) {
      range = nth_semiprime_approx(n) - nth_semiprime_approx(spcnt);
      range = 1.10 * range + 100;
      if (range > guess) range = guess;          /* just in case */
      if (range > 125000000) range = 125000000;  /* Not too many at a time */
      /* Get a bunch of semiprimes */
      MPUverbose(2, "  sieving forward %"UVuf"\n", range);
      count = range_semiprime_sieve(&S, guess+1, guess+range);
      if (spcnt+count <= n) {
        guess = S[count-1];
        spcnt += count;
      } else { /* Walk forwards */
        for (i = 0; i < count && spcnt < n; i++) {
          guess = S[i];
          spcnt++;
        }
      }
      Safefree(S);
    }
  } else if (n < spcnt && (spcnt-n) > SP_SIEVE_THRESH) {  /* sieve backwards */
    UV *S, count, range;
    while (n < spcnt) {
      range = nth_semiprime_approx(spcnt) - nth_semiprime_approx(n);
      range = 1.10 * range + 100;
      if (range > guess) range = guess;          /* just in case */
      if (range > 125000000) range = 125000000;  /* Not too many at a time */
      /* Get a bunch of semiprimes */
      MPUverbose(2, "  sieving backward %"UVuf"\n", range);
      count = range_semiprime_sieve(&S, guess-range, guess-1);
      if (spcnt-count >= n) {
        guess = S[0];
        spcnt -= count;
      } else { /* Walk backwards */
        while (count > 0 && n < spcnt) {
          guess = S[--count];
          spcnt--;
        }
      }
      Safefree(S);
    }
  }

  /* Finally, iterate over semiprimes until we hit the exact spot */
  for (; spcnt > n; spcnt--)
    guess = _prev_semiprime(guess);
  for (; spcnt < n; spcnt++)
    guess = _next_semiprime(guess);
  return guess;
}
