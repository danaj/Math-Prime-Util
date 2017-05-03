#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#define FUNC_log2floor 1
#include "util.h"
#define FUNC_is_prime_in_sieve 1
#include "prime_nth_count.h"
#include "sieve.h"
#include "ramanujan_primes.h"

/******************************************************************************/
/*                             RAMANUJAN PRIMES                               */
/******************************************************************************/

/* For Ramanujan prime estimates:
 *  - counts are done via inverse nth, so only one thing to tune.
 *  - For nth tables, upper values ok if too high, lower values ok if too low.
 *  - both upper & lower empirically tested to 175e9 (175 thousand million),
 *    with a return value of over 10^13.
 */

/* These are playing loose with Sondow/Nicholson/Noe 2011 theorem 4.
 * The last value should be rigorously checked using actual R_n values. */
static const uint32_t small_ram_upper_idx[] = {
  3971,3980,5219,5222,5225,5261,5264,5270,5276,5278,5324,5326,5554,7430,
  7448,7451,8580,8584,8607,12589,12603,12620,12729,18119,18134,18174,18289,
  18300,18401,18419,25799,27247,27267,28663,39635,40061,40366,45338,51320,
  64439,65566,65829,84761,89055,104959,107852,146968,151755,186499,217258,
  223956,270700,332195,347223,440804,508096,565039,768276,828377,1090285,
  1277320,1568165,1896508,2375799,3300765,4162908,5124977,6522443,9298256,
  11406250, 15528199, 20637716, 28239295, 39711166,
  55623925, 80161468, 117683515, 174200145, 261514813,
  404761877, 633278258, 1024431762, 1683645810, UVCONST(2868095461),
  UVCONST(4294967295)
};
#define SMALL_NRAM_UPPER_MULT 2852
#define SMALL_NRAM_UPPER (sizeof(small_ram_upper_idx)/sizeof(small_ram_upper_idx[0]))

#if BITS_PER_WORD == 64
static const UV large_ram_upper_idx[] = {
  UVCONST(     2209850625), UVCONST(     2505047527), UVCONST(     2868095461),
  /* 11071, 11070, 11069,   11068, 11067, 11066,   11065, 11064, 11063 */
  UVCONST(     3300302120), UVCONST(     3790618144), UVCONST(     4350498841),
  UVCONST(     5046044184), UVCONST(     5822787413), UVCONST(     6745591770),
  UVCONST(     7843974826), UVCONST(     9136430799), UVCONST(    10622435139),
  /* 11062, 11061, 11060,   11059, 11058, 11057,   11056, 11055, 11054 */
  UVCONST(    12458964476), UVCONST(    14547850676), UVCONST(    17105209669),
  UVCONST(    20137461647), UVCONST(    23754200841), UVCONST(    28113266506),
  UVCONST(    33244053524), UVCONST(    39635847319), UVCONST(    47216781912),
  /* 11053, 11052, 11051,   11050, 11049, 11048,   11047, 11046, 11045 */
  UVCONST(    56400647995), UVCONST(    67708204893), UVCONST(    81253508908),
  UVCONST(    98059403896), UVCONST(   118606509360), UVCONST(   143852101796),
  UVCONST(   175186371899), UVCONST(   213694229671), UVCONST(   261738209720),
  /* 11044, 11043, 11042,   11041, 11040, 11039,   11038, 11037, 11036 */
  UVCONST(   321608703183), UVCONST(   396839992097), UVCONST(   490698410355),
  UVCONST(   609608518202), UVCONST(   760145301247), UVCONST(   951394830802),
  UVCONST(  1195260877070), UVCONST(  1507711305081), UVCONST(  1910188609050),
  /* 11035, 11034, 11033,   11032, 11031, 11030,   11029, 11028, 11027 */
  UVCONST(  2429318332761), UVCONST(  3104355953395), UVCONST(  3983636468774),
  UVCONST(  5136852322734), UVCONST(  6654037637976), UVCONST(  8663089898701),
  UVCONST( 11333068937712), UVCONST( 14903446657241), UVCONST( 19702880416316),
  /* 11026, 11025, 11024,   11023, 11022, 11021,   11020, 11019, 11018 */
  UVCONST( 26192198741217), UVCONST( 35007464587700), UVCONST( 47072080175048),
  UVCONST( 63662369320719), UVCONST( 86631548367209), UVCONST(118640390543125),
  UVCONST(163548371026622), UVCONST(226944845224518), UVCONST(317123850443715),
  /* 11017, 11016, 11015,   11014, 11013, 11012,   11011, 11010, 11009 */
  UVCONST(446332665856000), /* estimates here and further */
  UVCONST(632810110373375),
  UVCONST(906250000000000),
  1.50*UVCONST(906250000000000),
  2.25*UVCONST(906250000000000),
  3.40*UVCONST(906250000000000),
  5.10*UVCONST(906250000000000),
  7.60*UVCONST(906250000000000),
};
#define LARGE_NRAM_UPPER_MULT 11075
#define LARGE_NRAM_UPPER (sizeof(large_ram_upper_idx)/sizeof(large_ram_upper_idx[0]))
#endif

UV nth_ramanujan_prime_upper(UV n) {
  UV i, mult, res;

  if (n <= 2) return (n==0) ? 0 : (n==1) ? 2 : 11;
  res = nth_prime_upper(3*n);

  if (n < UVCONST(2209850625) || BITS_PER_WORD < 64) {
    /* While p_3n is a complete upper bound, Rp_n tends to p_2n, and
     * SNN(2011) theorem 4 shows how we can find (m,c) values where m < 1,
     * Rn < m*p_3n for all n > c.  Here we use various quantized m values
     * and the table gives us c values where it applies. */
    if      (n < 20) mult = 3580;
    else if (n < 98) mult = 3340;
    else if (n < 1580) mult = 3040;
    else if (n < 3242) mult = 2885;
    else {
      for (i = 0; i < SMALL_NRAM_UPPER; i++)
        if (small_ram_upper_idx[i] > n)
          break;
      mult = SMALL_NRAM_UPPER_MULT-i;
    }
    if (res > (UV_MAX/mult)) res = (UV) (((long double) mult / 4096.0L) * res);
    else                     res = (res * mult) >> 12;
#if BITS_PER_WORD == 64
  } else {
    for (i = 0; i < LARGE_NRAM_UPPER; i++)
      if (large_ram_upper_idx[i] > n)
        break;
    mult = (LARGE_NRAM_UPPER_MULT-i);
    if (res > (UV_MAX/mult)) res = (UV) (((long double) mult / 16384.0L) * res);
    else                     res = (res * mult) >> 14;
#endif
  }
  return res;
}
static const uint32_t small_ram_lower_idx[] = {
  2786, 2801, 4275, 5935, 6107, 8797, 9556, 13314, 13641, 20457, 23745,
  34432, 50564, 69194, 97434, 149399, 224590, 337116, 514260, 804041,
  1367781, 2448102, 4513600, 8581572, 17555268, 39295429, 96043141,
  255707435, 771748027, UVCONST(2691176863), UVCONST(4294967295)
};
#define SMALL_NRAM_LOWER_MULT 557
#define SMALL_NRAM_LOWER (sizeof(small_ram_lower_idx)/sizeof(small_ram_lower_idx[0]))

#if BITS_PER_WORD == 64
static const UV large_ram_lower_idx[] = {
  UVCONST(    2287326511), UVCONST(    2691176863), UVCONST(    3186734059),
  UVCONST(    3773549435), UVCONST(    4490874244), UVCONST(    5351635422),
  UVCONST(    6379440261), UVCONST(    7669615016), UVCONST(    9222420352),
  UVCONST(   11116948392), UVCONST(   13468245052), UVCONST(   16342357141),
  /* 4213, 4212, 4211,   4210, 4209, 4208,   4207, 4206, 4205 */
  UVCONST(   19922979418), UVCONST(   24348237409), UVCONST(   29901283153),
  UVCONST(   36755324577), UVCONST(   45489957420), UVCONST(   56422464192),
  UVCONST(   70346667919), UVCONST(   87960867006), UVCONST(  110516423444),
  /* 4204, 4203, 4202,   4201, 4200, 4199,   4198, 4197, 4196 */
  UVCONST(  139376328492), UVCONST(  176769212370), UVCONST(  224987344210),
  UVCONST(  287637947625), UVCONST(  369783640662), UVCONST(  477404340380),
  UVCONST(  619830096553), UVCONST(  808881647349), UVCONST( 1061320360069),
  /* 4195, 4194, 4193,   4192, 4191, 4190,   4189, 4188, 4187 */
  UVCONST(  1400537557830), UVCONST(  1858188046854), UVCONST(  2480691884624),
  UVCONST(  3332007446416), UVCONST(  4503079799802), UVCONST(  6126600987695),
  UVCONST(  8391511555810), UVCONST( 11573178646156), UVCONST( 16075623921843),
  /* 4186, 4185, 4184,   4183, 4182, 4181,   4180, 4179, 4178 */
  UVCONST( 22497437657452), UVCONST( 31726512746004), UVCONST( 45096957879914),
  UVCONST( 64630470145902), UVCONST( 93416220154664), UVCONST(136212405571722),
  UVCONST(200433884291811),
  UVCONST(297735302623021),
  UVCONST(446609489920000), /* estimates here and further */
  /* 4177, 4176, 4175,   4174, 4173, 4172,   4171, 4170, 4169 */
  UVCONST(676524345920000),
  UVCONST(1025000000000000),
  1.50*UVCONST(1025000000000000),
  2.25*UVCONST(1025000000000000),
  3.37*UVCONST(1025000000000000),
  5.06*UVCONST(1025000000000000),
  7.59*UVCONST(1025000000000000),
};
#define LARGE_NRAM_LOWER_MULT 4225
#define LARGE_NRAM_LOWER (sizeof(large_ram_lower_idx)/sizeof(large_ram_lower_idx[0]))
#endif

UV nth_ramanujan_prime_lower(UV n) {
  UV res, i, mult;
  if (n <= 2) return (n==0) ? 0 : (n==1) ? 2 : 11;

  res = nth_prime_lower(2*n);

  if (n < UVCONST(2287326511) || BITS_PER_WORD < 64) {
    for (i = 0; i < SMALL_NRAM_LOWER; i++)
      if (small_ram_lower_idx[i] > n)
        break;
    mult = (SMALL_NRAM_LOWER_MULT-i);
    if (res > (UV_MAX/mult)) res = (UV) (((long double) mult / 512.0L) * res);
    else                     res = (res * mult) >> 9;
#if BITS_PER_WORD == 64
  } else {
    if (n < large_ram_lower_idx[LARGE_NRAM_LOWER-1]) {
      for (i = 0; i < LARGE_NRAM_LOWER; i++)
        if (large_ram_lower_idx[i] > n)
          break;
      mult = (LARGE_NRAM_LOWER_MULT-i);
      if (res > (UV_MAX/mult)) res = (UV) (((long double) mult / 4096.0L) * res);
      else                     res = (res * mult) >> 12;
    }
#endif
  }
  return res;
}

/* An advantage of making these binary searches on the inverse is that we
 * don't have to tune them separately, and nothing changes if the prime
 * count bounds are modified.  We do need to keep up to date with any
 * changes to nth_prime_{lower,upper} however. */

UV ramanujan_prime_count_lower(UV n) {
  UV lo, hi;
  if (n < 29) return (n < 2) ? 0 : (n < 11) ? 1 : (n < 17) ? 2 : 3;
  /* Binary search on nth_ramanujan_prime_upper */
  /* We know we're between p_2n and p_3n, probably close to the former. */
  lo = prime_count_lower(n)/3;
  hi = prime_count_upper(n) >> 1;
  while (lo < hi) {
    UV mid = lo + (hi-lo)/2;
    if (nth_ramanujan_prime_upper(mid) < n) lo = mid+1;
    else                                    hi = mid;
  }
  return lo-1;
}
UV ramanujan_prime_count_upper(UV n) {
  /* return prime_count_upper(n) >> 1; */       /* Simple bound */
  UV lo, hi;
  if (n < 29) return (n < 2) ? 0 : (n < 11) ? 1 : (n < 17) ? 2 : 3;
  /* Binary search on nth_ramanujan_prime_upper */
  /* We know we're between p_2n and p_3n, probably close to the former. */
  lo = prime_count_lower(n)/3;
  hi = prime_count_upper(n) >> 1;
  while (lo < hi) {
    UV mid = lo + (hi-lo)/2;
    if (nth_ramanujan_prime_lower(mid) < n) lo = mid+1;
    else                                    hi = mid;
  }
  return lo-1;
}

/* Return array of first n ramanujan primes.  Use Noe's algorithm. */
UV* n_ramanujan_primes(UV n) {
  UV max, k, s, *L;
  unsigned char* sieve;
  max = nth_ramanujan_prime_upper(n); /* Rn <= max, so we can sieve to there */
  if (_XS_get_verbose() >= 2) { printf("sieving to %"UVuf" for first %"UVuf" Ramanujan primes\n", max, n); fflush(stdout); }
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
  int verbose = _XS_get_verbose();

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
  if (verbose >= 2) { printf("Rn[%"UVuf"] to Rn[%"UVuf"]     Noe's: %"UVuf" to %"UVuf"\n", nlo, nhi, mink, maxk); fflush(stdout); }

  s = 1 + prime_count(2,mink-2) - prime_count(2,(mink-1)>>1);
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
  if (verbose >= 2) { printf("Generated %"UVuf" Ramanujan primes from %"UVuf" to %"UVuf"\n", nhi-nlo+1, L[0], L[nhi-nlo]); fflush(stdout); }
  return L;
}

UV nth_ramanujan_prime(UV n) {
  UV rn, *L;
  if (n <= 2) return (n == 0) ? 0 : (n == 1) ? 2 : 11;
  L = n_range_ramanujan_primes(n, n);
  rn = L[0];
  Safefree(L);
  return rn;
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

int is_ramanujan_prime(UV n) {
  UV beg, end, *L;

  if (!is_prime(n))  return 0;
  if (n < 17)        return (n == 2 || n == 11);

  /* Generate Ramanujan primes and see if we're in the list.  Slow. */
  L = ramanujan_primes(&beg, &end, n, n);
  Safefree(L);
  return (beg <= end);
}

UV ramanujan_prime_count_approx(UV n)
{
  /* Binary search on nth_ramanujan_prime_approx */
  UV lo, hi;
  if (n < 29) return (n < 2) ? 0 : (n < 11) ? 1 : (n < 17) ? 2 : 3;
  lo = ramanujan_prime_count_lower(n);
  hi = ramanujan_prime_count_upper(n);
  while (lo < hi) {
    UV mid = lo + (hi-lo)/2;
    if (nth_ramanujan_prime_approx(mid) < n) lo = mid+1;
    else                                     hi = mid;
  }
  return lo-1;
}

UV nth_ramanujan_prime_approx(UV n)
{
  UV lo = nth_ramanujan_prime_lower(n),  hi = nth_ramanujan_prime_upper(n);
  /* Our upper bounds come out much closer, so weight toward them. */
  double weight = (n <= UVCONST(4294967295))  ?  1.62  :  1.51;
  return lo + weight * ((hi-lo) >> 1);
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

static UV _ramanujan_prime_count(UV n) {
  UV i, v, rn, *L, window, swin, ewin, wlen, log2 = log2floor(n), winmult = 1;

  if (n <= 10) return (n < 2) ? 0 : 1;

  /* We have some perfect powers of 2 in our table */
  if ((n & (n-1)) == 0 && log2 <= RAMPC2)
    return ramanujan_counts_pow2[log2];

  if (_XS_get_verbose()) { printf("ramanujan_prime_count calculating Pi(%lu)\n",n); fflush(stdout); }
  v = prime_count(2,n) - prime_count(2,n >> 1);

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
    if (_XS_get_verbose()) { printf("  ramanujan_prime_count increasing window\n"); fflush(stdout); }
  }
  rn = swin + i - 1;
  Safefree(L);
  return rn;
}

UV ramanujan_prime_count(UV lo, UV hi)
{
  UV count;

  if (hi < 2 || hi < lo) return 0;

#if 1
  count = _ramanujan_prime_count(hi);
  if (lo > 2)
    count -= _ramanujan_prime_count(lo-1);
#else
  {
    UV beg, end, *L;
    /* Generate all Rp from lo to hi */
    L = ramanujan_primes(&beg, &end, lo, hi);
    count = (L && end >= beg) ? end-beg+1 : 0;
    Safefree(L);
  }
#endif
  return count;
}
