#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FUNC_isqrt 1
#define FUNC_icbrt 1
#include "ptypes.h"
#include "constants.h"
#include "util.h"
#include "cache.h"
#include "sieve.h"
#include "prime_sums.h"

/******************************************************************************/
/*                                   SUMS                                     */
/******************************************************************************/

/* As an aside, good information about bounds and approximations can be
 * found in Axler (2019) "On the sum of the first n prime numbers"
 * https://jtnb.centre-mersenne.org/item/10.5802/jtnb.1081.pdf
 */

static const unsigned char byte_zeros[256] =
  {8,7,7,6,7,6,6,5,7,6,6,5,6,5,5,4,7,6,6,5,6,5,5,4,6,5,5,4,5,4,4,3,
   7,6,6,5,6,5,5,4,6,5,5,4,5,4,4,3,6,5,5,4,5,4,4,3,5,4,4,3,4,3,3,2,
   7,6,6,5,6,5,5,4,6,5,5,4,5,4,4,3,6,5,5,4,5,4,4,3,5,4,4,3,4,3,3,2,
   6,5,5,4,5,4,4,3,5,4,4,3,4,3,3,2,5,4,4,3,4,3,3,2,4,3,3,2,3,2,2,1,
   7,6,6,5,6,5,5,4,6,5,5,4,5,4,4,3,6,5,5,4,5,4,4,3,5,4,4,3,4,3,3,2,
   6,5,5,4,5,4,4,3,5,4,4,3,4,3,3,2,5,4,4,3,4,3,3,2,4,3,3,2,3,2,2,1,
   6,5,5,4,5,4,4,3,5,4,4,3,4,3,3,2,5,4,4,3,4,3,3,2,4,3,3,2,3,2,2,1,
   5,4,4,3,4,3,3,2,4,3,3,2,3,2,2,1,4,3,3,2,3,2,2,1,3,2,2,1,2,1,1,0};

/* The fastest way to compute the sum of primes is using a combinatorial
 * algorithm such as Deleglise-Rivat or Gourdon.  This is what Kim Walisch's
 * primesum program does.  Note that one quickly needs 128-bit or larger
 * storage, as the sums grow rapidly.
 *
 * We are using much simpler methods.  Performance at small sizes is also a
 * consideration.  Using tables combined with summing over sieved primes can
 * work well with small input sizes.
 */

/* Simplified Legendre method giving pisum(n) for n <= 65535 or 4294967295. */

UV sum_primes64(UV n) {
  uint32_t *V, j, k, r, r2, p;
  UV *S, sum;

  if (n < 2 || (n >> (BITS_PER_WORD/2)) > 0)  /* S[] will overflow */
    return 0;
  r = isqrt(n);
  r2 = r + n/(r+1);

  New(0, V, r2+1, uint32_t);
  New(0, S, r2+1, UV);
  for (k = 1; k <= r2; k++) {
    UV v = (k <= r)  ?  k  :  n/(r2-k+1);
    V[k] = v;
    S[k] = ((v*(v-1))>>1) + (v-1);
  }

  for (p = 2; p <= r; p++) {
    if (S[p] > S[p-1]) { /* For each prime p from 2 to r */
      UV sp = S[p-1], p2 = p*p;
      for (j = r2; j > 1 && V[j] >= p2; j--) {
        uint32_t a = V[j], b = a/p;
        if (a > r) a = r2 - n/a + 1;
        if (b > r) b = r2 - n/b + 1;
        S[a] -= (UV)p * (S[b] - sp);   /* sp = sum of primes less than p */
      }
    }
  }
  sum = S[r2];
  Safefree(V);
  Safefree(S);
  return sum;
}

/* Simplified Legendre method giving pisum(n) for any 64-bit input n,
 * assuming the uint128_t type is available.  The result is returned as
 * two 64-bit results. */

int sum_primes128(UV n, UV *hi_sum, UV *lo_sum) {
#if HAVE_SUM_PRIMES128
  uint128_t *S;
  UV *V, j, k, r, r2, p;

  /* pisum(2^64-1) < 2^128-1, so no overflow issues */
  r = isqrt(n);
  r2 = r + n/(r+1);

  New(0, V, r2+1, UV);
  New(0, S, r2+1, uint128_t);
  for (k = 0; k <= r2; k++) {
    uint128_t v = (k <= r)  ?  k  :  n/(r2-k+1);
    V[k] = v;
    S[k] = (v*(v+1))/2 - 1;
  }

  for (p = 2; p <= r; p++) {
    if (S[p] > S[p-1]) { /* For each prime p from 2 to r */
      uint128_t sp = S[p-1], p2 = ((uint128_t)p) * p;
      for (j = r2; j > 1 && V[j] >= p2; j--) {
        UV a = V[j], b = a/p;
        if (a > r) a = r2 - n/a + 1;
        if (b > r) b = r2 - n/b + 1;
        S[a] -= p * (S[b] - sp);   /* sp = sum of primes less than p */
      }
    }
  }
  *hi_sum = (S[r2] >> 64) & UV_MAX;
  *lo_sum = (S[r2]      ) & UV_MAX;
  Safefree(V);
  Safefree(S);
  return 1;
#else
  return 0;
#endif
}


/* sum primes in a 64-bit range using a sieving with table acceleration */

static const unsigned char byte_sum[256] =
  {120,119,113,112,109,108,102,101,107,106,100,99,96,95,89,88,103,102,96,95,92,
   91,85,84,90,89,83,82,79,78,72,71,101,100,94,93,90,89,83,82,88,87,81,80,77,
   76,70,69,84,83,77,76,73,72,66,65,71,70,64,63,60,59,53,52,97,96,90,89,86,85,
   79,78,84,83,77,76,73,72,66,65,80,79,73,72,69,68,62,61,67,66,60,59,56,55,49,
   48,78,77,71,70,67,66,60,59,65,64,58,57,54,53,47,46,61,60,54,53,50,49,43,42,
   48,47,41,40,37,36,30,29,91,90,84,83,80,79,73,72,78,77,71,70,67,66,60,59,74,
   73,67,66,63,62,56,55,61,60,54,53,50,49,43,42,72,71,65,64,61,60,54,53,59,58,
   52,51,48,47,41,40,55,54,48,47,44,43,37,36,42,41,35,34,31,30,24,23,68,67,61,
   60,57,56,50,49,55,54,48,47,44,43,37,36,51,50,44,43,40,39,33,32,38,37,31,30,
   27,26,20,19,49,48,42,41,38,37,31,30,36,35,29,28,25,24,18,17,32,31,25,24,21,
   20,14,13,19,18,12,11,8,7,1,0};

#if BITS_PER_WORD == 64
/* We have a much more limited range, so use a fixed interval.  We should be
 * able to get any 64-bit sum in under a half-second. */
static const UV sum_table_2e8[] =
  {1075207199997324,3071230303170813,4990865886639877,6872723092050268,8729485610396243,10566436676784677,12388862798895708,14198556341669206,15997206121881531,17783028661796383,19566685687136351,21339485298848693,23108856419719148,
   24861364231151903,26619321031799321,28368484289421890,30110050320271201,31856321671656548,33592089385327108,35316546074029522,37040262208390735,38774260466286299,40490125006181147,42207686658844380,43915802985817228,45635106002281013,
   47337822860157465,49047713696453759,50750666660265584,52449748364487290,54152689180758005,55832433395290183,57540651847418233,59224867245128289,60907462954737468,62597192477315868,64283665223856098,65961576139329367,67641982565760928,
   69339211720915217,71006044680007261,72690896543747616,74358564592509127,76016548794894677,77694517638354266,79351385193517953,81053240048141953,82698120948724835,84380724263091726,86028655116421543,87679091888973563,89348007111430334,
   90995902774878695,92678527127292212,94313220293410120,95988730932107432,97603162494502485,99310622699836698,100935243057337310,102572075478649557,104236362884241550,105885045921116836,107546170993472638,109163445284201278,
   110835950755374921,112461991135144669,114116351921245042,115740770232532531,117408250788520189,119007914428335965,120652479429703269,122317415246500401,123951466213858688,125596789655927842,127204379051939418,128867944265073217,
   130480037123800711,132121840147764197,133752985360747726,135365954823762234,137014594650995101,138614165689305879,140269121741383097,141915099618762647,143529289083557618,145146413750649432,146751434858695468,148397902396643807,
   149990139346918801,151661665434334577,153236861034424304,154885985064643097,156500983286383741,158120868946747299,159735201435796748,161399264792716319,162999489977602579,164566400448130092,166219688860475191,167836981098849796,
   169447127305804401,171078187147848898,172678849082290997,174284436375728242,175918609754056455,177525046501311788,179125593738290153,180765176633753371,182338473848291683,183966529541155489,185585792988238475,187131988176321434,
   188797837140841381,190397649440649965,191981841583560122,193609739194967419,195166830650558070,196865965063113041,198400070713177440,200057161591648721,201621899486413406,203238279253414934,204790684829891896,206407676204061001,
   208061050481364659,209641606658938873,211192088300183855,212855420483750498,214394145510853736,216036806225784861,217628995137940563,219277567478725189,220833877268454872,222430818525363309,224007307616922530,225640739533952807,
   227213096159236934,228853318075566255,230401824696558125,231961445347821085,233593317860593895,235124654760954338,236777716068869769,238431514923528303,239965003913481640,241515977959535845,243129874530821395};
#define N_SUM_TABLE  (sizeof(sum_table_2e8)/sizeof(sum_table_2e8[0]))
#endif

int sum_primes(UV low, UV high, UV *return_sum) {
  UV sum = 0;
  int overflow = 0;

  if (low <= 2 && high >= 100000) {
    *return_sum = sum_primes64(high);
    if (*return_sum != 0)
      return 1;
  }
  /* TODO: performance: more cases where using sum_primes64 is faster. */

  if ((low <= 2) && (high >= 2)) sum += 2;
  if ((low <= 3) && (high >= 3)) sum += 3;
  if ((low <= 5) && (high >= 5)) sum += 5;
  if (low < 7) low = 7;

  /* If we know the range will overflow, return now */
#if BITS_PER_WORD == 64
  if (low == 7 && high >= 29505444491)  return 0;
  if (low >= 1e10 && (high-low) >= 32e9) return 0;
  if (low >= 1e13 && (high-low) >=  5e7) return 0;
#else
  if (low == 7 && high >= 323381)  return 0;
#endif

#if 1 && BITS_PER_WORD == 64    /* Tables */
  if (low == 7 && high >= 2e8) {
    UV step;
    for (step = 1; high >= (step * 2e8) && step < N_SUM_TABLE; step++) {
      sum += sum_table_2e8[step-1];
      low = step * 2e8;
    }
  }
#endif

  if (low <= high) {
    unsigned char* segment;
    UV seg_base, seg_low, seg_high;
    void* ctx = start_segment_primes(low, high, &segment);
    while (!overflow && next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      UV bytes = seg_high/30 - seg_low/30 + 1;
      unsigned char s;
      unsigned char* sp = segment;
      unsigned char* const spend = segment + bytes - 1;
      UV i, p, pbase = 30*(seg_low/30);

      /* Clear primes before and after our range */
      p = pbase;
      for (i = 0; i < 8 && p+wheel30[i] < low; i++)
        if ( (*sp & (1<<i)) == 0 )
          *sp |= (1 << i);

      p = 30*(seg_high/30);
      for (i = 0; i < 8;  i++)
        if ( (*spend & (1<<i)) == 0 && p+wheel30[i] > high )
          *spend |= (1 << i);

      while (sp <= spend) {
        s = *sp++;
        if (sum < (UV_MAX >> 3) && pbase < (UV_MAX >> 5)) {
          /* sum block of 8 all at once */
          sum += pbase * byte_zeros[s] + byte_sum[s];
        } else {
          /* sum block of 8, checking for overflow at each step */
          for (i = 0; i < byte_zeros[s]; i++) {
            if (sum+pbase < sum) overflow = 1;
            sum += pbase;
          }
          if (sum+byte_sum[s] < sum) overflow = 1;
          sum += byte_sum[s];
          if (overflow) break;
        }
        pbase += 30;
      }
    }
    end_segment_primes(ctx);
  }
  if (!overflow && return_sum != 0)  *return_sum = sum;
  return !overflow;
}
