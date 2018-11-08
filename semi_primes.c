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
#include "semi_primes.h"

#define SP_SIEVE_THRESH 100  /* When to sieve vs. iterate */

/******************************************************************************/
/*                                SEMI PRIMES                                 */
/******************************************************************************/

static const unsigned char _semiprimelist[] =
  {0,4,6,9,10,14,15,21,22,25,26,33,34,35,38,39,46,49,51,55,57,58,62,65,69,74,
   77,82,85,86,87,91,93,94,95,106,111,115,118,119,121,122,123,129,133,134,141,
   142,143,145,146,155,158,159,161,166,169,177,178,183,185,187,194,201,202,
   203,205,206,209,213,214,215,217,218,219,221,226,235,237,247,249,253,254};
#define NSEMIPRIMELIST (sizeof(_semiprimelist)/sizeof(_semiprimelist[0]))

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
            xarr = array_of_primes_in_range(&xsize, xbeg, xend);
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

/* TODO: This overflows, see p=3037000507,lo=10739422018595509581.
 *       p2 = 9223372079518257049 => 9223372079518257049 + 9223372079518257049
 *       Also with lo=18446744073709551215,hi=18446744073709551515.
 */
#define PGTLO(ip,p,lo)  ((ip)>=(lo)) ? (ip) : ((p)*((lo)/(p)) + (((lo)%(p))?(p):0))
#define MARKSEMI(p,arr,lo,hi) \
    do {  UV i, p2=(p)*(p); \
      for (i = PGTLO(p2, p, lo); i >= lo && i <= hi; i += p) arr[i-lo]++; \
      for (i = PGTLO(2*p2, p2, lo); i >= lo && i <= hi; i += p2) arr[i-lo]++; \
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
static int _numcmp(const void *a, const void *b)
  { const UV *x = a, *y = b; return (*x > *y) ? 1 : (*x < *y) ? -1 : 0; }

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

  pr = array_of_primes_in_range(&xsize, 0, lim);

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

UV semiprime_count_approx(UV n) {
  if (n <= _semiprimelist[NSEMIPRIMELIST-1]) {
    UV i = 0;
    while (i < NSEMIPRIMELIST-1 && n >=  _semiprimelist[i+1])
      i++;
    return i;
  } else {
    UV lo, hi;
    double init, logn = log(n), loglogn = log(logn);
    /* init = n * loglogn / logn; */
    /* init = (n/logn) * (0.11147910114 + 0.00223801350*logn + 0.44233207922*loglogn + 1.65236647896*log(loglogn)); */
    init = n * (loglogn + 0.302) / logn;
    if (1.05*init >= (double)UV_MAX)
      return init;

    lo = 0.90 * init - 5,   hi = 1.05 * init;
    while (lo < hi) {
      UV mid = lo + (hi-lo)/2;
      if (nth_semiprime_approx(mid) < n) lo = mid+1;
      else                               hi = mid;
    }
    return lo;
  }
}

UV nth_semiprime_approx(UV n) {
  double logn,log2n,log3n,log4n, err_lo, err_md, err_hi, err_factor, est;

  if (n < NSEMIPRIMELIST)
    return _semiprimelist[n];

  /* Piecewise with blending.  Hacky and maybe overkill, but it makes
   * a big performance difference, especially at the high end.
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
  est = 0.5 + err_factor * n * logn / log2n;
  if (est >= UV_MAX) return 0;
  return (UV)est;
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
