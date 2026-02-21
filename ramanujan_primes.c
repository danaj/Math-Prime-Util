#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#define FUNC_log2floor 1
#include "util.h"
#define FUNC_is_prime_in_sieve 1
#include "prime_counts.h"
#include "sieve.h"
#include "inverse_interpolate.h"
#include "ramanujan_primes.h"

/******************************************************************************/
/*                             RAMANUJAN PRIMES                               */
/******************************************************************************/

static const uint8_t small_ram_primes[] = {
2,11,17,29,41,47,59,67,71,97,101,107,127,149,151,167,179,181,227,229,233,239,241
};
#define NSMALL_RAM (sizeof(small_ram_primes)/sizeof(small_ram_primes[0]))

#define FAST_SMALL_NTH(n) \
  if (n <= NSMALL_RAM) \
    { return (n == 0) ? 0 : small_ram_primes[n-1]; }
#define FAST_SMALL_COUNT(n) \
 if (n <= small_ram_primes[NSMALL_RAM-1]) { \
   UV i; \
   for (i = 0; i < NSMALL_RAM; i++) \
     if (n < small_ram_primes[i]) break; \
   return i; \
 }


/*******************************     Bounds     *******************************/

/* Upper and lower bounds done using Axler 2017:
 *    https://arxiv.org/pdf/1711.04588.pdf
 * The parameter values have been computed using exact nth_prime,
 * so does not depend on the nth_prime_upper / nth_prime_lower method.
 */

UV nth_ramanujan_prime_upper(UV n) {
  long double R = 0, D = 0.565;

  FAST_SMALL_NTH(n);

  if (n < 12581) {
    if      (n <       168) R = ramanujan_axler(n, -4.7691, -6.2682);
    else if (n <      2290) R = ramanujan_axler(n, -0.9315, -0.5635);
    else if (n <      5225) R = ramanujan_axler(n, -0.5318, -0.0710);
    else if (n <     12581) R = ramanujan_axler(n,  0.1212,  0.7973);
    return nth_prime_upper(R);
  }
  if      (n <     18175) D =  0.3548;
  else if (n <     82883) D = -0.2450;
  else if (n <    316314) D = -0.6384;
  else if (n <   1000001) D = -0.9353;
  else if (n <   4000001) D = -1.1271;
  else if (n <  16000001) D = -1.4152;
  else if (n <  64000001) D = -1.6671;
  else if (n < 128000001) D = -1.8855;
  else if (n < 256000001) D = -1.9325;
  else if (n < 384000001) D = -2.0190;
  else if (n < 512000001) D = -2.0310;
  else {
    D = -2.0884;
    if (n > UVCONST(     3999654659)) D = -2.235;
#if BITS_PER_WORD == 64
    if (n > UVCONST(    84086679236)) D = -2.435;
    if (n > UVCONST(   514808375201)) D = -2.535;
    if (n > UVCONST(  3594243587299)) D = -2.635;
    if (n > UVCONST( 28330126673435)) D = -2.735;
    if (n > UVCONST(117462814829787)) D = -2.8;
#endif
  }

  return nth_prime_upper(ramanujan_axler(n, 0.0, D));
}

UV nth_ramanujan_prime_lower(UV n) {
  double R = 0, D = 0;

  FAST_SMALL_NTH(n);

  if (n < 34816) {
    if      (n <       189) R = ramanujan_axler(n,  4.2720,  0.340);
    else if (n <      1245) R = ramanujan_axler(n, -0.2179, -6.179);
    else if (n <      2984) R = ramanujan_axler(n,  0.1446, -4.8693);
    else if (n <     14303) R = ramanujan_axler(n, -0.3570, -5.1154);
    else if (n <     34816) R = ramanujan_axler(n, -1.5770, -7.5332);
    return nth_prime_lower(R);
  }

  if      (n <     76400) D = 0.0126;
  else if (n <    280816) D = 0.5132;
  else if (n <    915887) D = 0.9967;
  else if (n <   4000001) D = 1.5004;
  else if (n <  16000001) D = 1.7184;
  else if (n <  64000001) D = 1.9860;
  else if (n < 128000001) D = 2.1352;
  else if (n < 256000001) D = 2.1658;
  else if (n < 384000001) D = 2.1999;
  else if (n < 512000001) D = 2.2047;
  else if (n < 640000001) D = 2.2324;
  else {
    D = 2.2385;
#if BITS_PER_WORD == 64
    if (n > UVCONST(    14888378285)) D = 2.29;
    if (n > UVCONST(   467037926604)) D = 2.31;
    if (n > UVCONST(  2778491401197)) D = 2.315;
    if (n > UVCONST( 10656144781918)) D = 2.317;
    if (n > UVCONST( 63698770351741)) D = 2.319;
#endif
  }

  return nth_prime_lower(ramanujan_axler(n, 1.472, D));
}

/* For Ramanujan prime count bounds, use binary searches on the inverses. */

UV ramanujan_prime_count_lower(UV n) {
  UV lo, hi;
  FAST_SMALL_COUNT(n);
  /* We know we're between p_2n and p_3n, probably close to the former. */
  lo = prime_count_lower(n)/3;
  hi = prime_count_upper(n) >> 1;
  return inverse_interpolate(lo, hi, n, &nth_ramanujan_prime_upper, 0);
}
UV ramanujan_prime_count_upper(UV n) {
  UV lo, hi;
  FAST_SMALL_COUNT(n);
  /* We know we're between p_2n and p_3n, probably close to the former. */
  lo = prime_count_lower(n)/3;
  hi = prime_count_upper(n) >> 1;
  return inverse_interpolate(lo, hi, n, &nth_ramanujan_prime_lower, 0);
}

/****************************     Approximate     ****************************/

UV ramanujan_prime_count_approx(UV n)
{
  FAST_SMALL_COUNT(n);
  /* Extremely accurate but a bit slow */
  return prime_count_approx(n) - prime_count_approx(n >> 1);
}

UV nth_ramanujan_prime_approx(UV n)
{
  UV lo, hi;
  FAST_SMALL_NTH(n);
  /* Interpolating using ramanujan prime count approximation */
  lo = nth_ramanujan_prime_lower(n) - 1;
  hi = nth_ramanujan_prime_upper(n);
  return inverse_interpolate(lo, hi, n, &ramanujan_prime_count_approx, 0);
}

/*******************************     Arrays     *******************************/

/* Return array of first n ramanujan primes.  Use Noe's algorithm. */
UV* n_ramanujan_primes(UV n) {
  UV max, k, s, *L;
  unsigned char* sieve;

  if (n <= NSMALL_RAM) {
    New(0, L, n, UV);
    for (k = 0; k < n; k++)
      L[k] = small_ram_primes[k];
    return L;
  }

  max = nth_ramanujan_prime_upper(n); /* Rn <= max, so we can sieve to there */
  MPUverbose(2, "sieving to %"UVuf" for first %"UVuf" Ramanujan primes\n", max, n);
  Newz(0, L, n, UV);
  L[0] = 2;
  sieve = sieve_erat30(max);
  for (s = 0, k = 7; k <= max; k += 2) {
    if (is_prime_in_sieve(sieve, k)) s++;
    if (s < n) L[s] = k+1;
    if ((k & 3) == 1 && is_prime_in_sieve(sieve, (k+1)>>1)) s--;
    if (s < n) L[s] = k+2;
  }
  Safefree(sieve);
  return L;
}

UV* n_range_ramanujan_primes(UV nlo, UV nhi) {
  UV mink, maxk, k, s, *L;

  if (nlo == 0) nlo = 1;
  if (nhi == 0) nhi = 1;

  /* If we're starting from 1, just do single monolithic sieve */
  if (nlo == 1)  return n_ramanujan_primes(nhi);

  Newz(0, L, nhi-nlo+1, UV);
  if (nlo <= 1 && nhi >= 1) L[1-nlo] =  2;
  if (nlo <= 2 && nhi >= 2) L[2-nlo] = 11;
  if (nhi < 3) return L;

  mink = nth_ramanujan_prime_lower(nlo) - 1;
  maxk = nth_ramanujan_prime_upper(nhi) + 1;

  if (mink < 15) mink = 15;
  if (mink % 2 == 0) mink--;
  MPUverbose(2, "Rn[%"UVuf"] to Rn[%"UVuf"]     Noe's: %"UVuf" to %"UVuf"\n", nlo, nhi, mink, maxk);

  s = 1 + prime_count(mink-2) - prime_count((mink-1)>>1);
  {
    unsigned char *segment, *seg2 = 0;
    void* ctx = start_segment_primes(mink, maxk, &segment);
    UV seg_base, seg_low, seg_high, new_size, seg2beg, seg2end, seg2size = 0;
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      seg2beg = 30 * (((seg_low+1)>>1)/30);
      seg2end = 30 * ((((seg_high+1)>>1)+29)/30);
      new_size = (seg2end - seg2beg)/30 + 1;
      if (new_size > seg2size) {
        if (seg2size > 0) Safefree(seg2);
        New(0, seg2, new_size, unsigned char);
        seg2size = new_size;
      }
      (void) sieve_segment(seg2, seg2beg/30, seg2end/30);
      for (k = seg_low; k <= seg_high; k += 2) {
        if (is_prime_in_sieve(segment, k-seg_base)) s++;
        if (s >= nlo && s <= nhi) L[s-nlo] = k+1;
        if ((k & 3) == 1 && is_prime_in_sieve(seg2, ((k+1)>>1)-seg2beg)) s--;
        if (s >= nlo && s <= nhi) L[s-nlo] = k+2;
      }
    }
    end_segment_primes(ctx);
    Safefree(seg2);
  }
  MPUverbose(2, "Generated %"UVuf" Ramanujan primes from %"UVuf" to %"UVuf"\n", nhi-nlo+1, L[0], L[nhi-nlo]);
  return L;
}

/* Returns array of Ram primes between low and high, results from first->last */
UV* ramanujan_primes(UV* first, UV* last, UV low, UV high)
{
  UV nlo, nhi, *L, lo, hi, mid;

  if (high < 2 || high < low) return 0;
  if (low < 2) low = 2;

  nlo = ramanujan_prime_count_lower(low);
  nhi = ramanujan_prime_count_upper(high);
  L = n_range_ramanujan_primes(nlo, nhi);

  /* Search for first entry in range */
  for (lo = 0, hi = nhi-nlo+1;  lo < hi;  ) {
    mid = lo + (hi-lo)/2;
    if (L[mid]  <  low)  lo = mid+1;
    else                 hi = mid;
  }
  *first = lo;
  /* Search for last entry in range */
  for (hi = nhi-nlo+1;  lo < hi;  ) {
    mid = lo + (hi-lo)/2;
    if (L[mid] <= high)  lo = mid+1;
    else                 hi = mid;
  }
  *last = lo-1;
  return L;
}

UV range_ramanujan_prime_sieve(UV** list, UV lo, UV hi)
{
  UV first, last, *L;
  L = ramanujan_primes(&first, &last, lo, hi);
  if (L == 0 || first > last) { *list = 0; return 0; }
  if (first > 0)
    memmove( L + 0,  L + first,  (last-first+1) * sizeof(UV) );
  *list = L;
  return last-first+1;
}

/* Generate a small window of Rp's around n */
static UV* _ramanujan_prime_window(UV n, UV* winsize, UV* npos) {
  UV i, v, *L, window, swin, ewin, wlen, winmult = 1;

  MPUverbose(1, "ramanujan_prime_count calculating Pi(%"UVuf")\n",n);
  v = prime_count(n) - prime_count(n >> 1);

  /* For large enough n make a slightly bigger window */
  if (n > 1000000000U) winmult = 16;

  while (1) {
    window = 20 * winmult;
    swin = (v <= window) ? 1 : v-window;
    ewin = v+window;
    wlen = ewin-swin+1;
    L = n_range_ramanujan_primes(swin, ewin);
    if (L[0] < n && L[wlen-1] > n) {
      /* Naive linear search from the start. */
      for (i = 1; i < wlen; i++)
        if (L[i] > n && L[i-1] <= n)
          break;
      if (i < wlen) break;
    }
    winmult *= 2;
    MPUverbose(1, "  %s increasing window\n", "ramanujan_prime_count");
  }
  *winsize = swin;
  *npos = i-1;
  return L;
}

/*******************************     Exact     *******************************/

UV nth_ramanujan_prime(UV n) {
  UV rn, *L;
  FAST_SMALL_NTH(n);
  L = n_range_ramanujan_primes(n, n);
  rn = L[0];
  Safefree(L);
  return rn;
}

bool is_ramanujan_prime(UV n) {
  UV i, d, *L, swin, rn;
  bool res;

  if (!is_prime(n))  return 0;
  if (n < 17)        return (n == 2 || n == 11);

  /* Pre-test: Check if Pi(n/2) increases before Pi(n) does. */
  if (is_prime(n/2+1)) return 0;
  d = (next_prime(n) - n)/2;
  for (i = 2; i <= d; i++)
    if (is_prime(n/2+i)) return 0;

  /* Very straightforward, but not the fastest method:
   *   return nth_ramanujan_prime(ramanujan_prime_count(n)) == n;
   *
   * Slower than below for most input sizes:
   *   L = ramanujan_primes(&beg, &end, n, n);
   *   Safefree(L);
   *   return (beg <= end);
   */

  L = _ramanujan_prime_window(n, &swin, &rn);
  res = (L[rn] == n);
  Safefree(L);
  return res;
}

#if BITS_PER_WORD == 64
#define RAMPC2 56
static const UV ramanujan_counts_pow2[RAMPC2+1] = {
   0, 1, 1, 1, 2, 4, 7, 13, 23, 42, 75, 137, 255, 463, 872, 1612,
   3030, 5706, 10749, 20387, 38635, 73584, 140336, 268216, 513705,
   985818, 1894120, 3645744, 7027290, 13561906, 26207278, 50697533,
   98182656, 190335585, 369323301, 717267167,
   UVCONST(     1394192236), UVCONST(     2712103833), UVCONST(     5279763823),
   UVCONST(    10285641777), UVCONST(    20051180846), UVCONST(    39113482639),
   UVCONST(    76344462797), UVCONST(   149100679004), UVCONST(   291354668495),
   UVCONST(   569630404447), UVCONST(  1114251967767), UVCONST(  2180634225768),
   UVCONST(  4269555883751), UVCONST(  8363243713305), UVCONST( 16388947026629),
   UVCONST( 32129520311897), UVCONST( 63012603695171), UVCONST(123627200537929),
   UVCONST(242637500756376), UVCONST(476379740340417), UVCONST(935609435783647) };
#else
#define RAMPC2 31  /* input limited */
static const UV ramanujan_counts_pow2[RAMPC2+1] = {
   0, 1, 1, 1, 2, 4, 7, 13, 23, 42, 75, 137, 255, 463, 872, 1612,
   3030, 5706, 10749, 20387, 38635, 73584, 140336, 268216, 513705,
   985818, 1894120, 3645744, 7027290, 13561906, 26207278, 50697533 };
#endif

UV ramanujan_prime_count(UV n) {
  UV swin, rn, *L, log2 = log2floor(n);

  if ((n & (n-1)) == 0 && log2 <= RAMPC2) /* Powers of two from table */
    return ramanujan_counts_pow2[log2];
  FAST_SMALL_COUNT(n);

  L = _ramanujan_prime_window(n, &swin, &rn);
  Safefree(L);
  return swin+rn;
}

UV ramanujan_prime_count_range(UV lo, UV hi)
{
  if (hi < 2 || hi < lo) return 0;
  return ramanujan_prime_count(hi)  -  ((lo <= 2) ? 0 : ramanujan_prime_count(lo-1));
}
