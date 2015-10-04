#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*****************************************************************************
 *
 * Prime counts using the extended Lagarias-Miller-Odlyzko combinatorial method.
 *
 * Copyright (c) 2013-2014 Dana Jacobsen (dana@acm.org)
 * This is free software; you can redistribute it and/or modify it under
 * the same terms as the Perl 5 programming language system itself.
 *
 * This file is part of the Math::Prime::Util Perl module, but it should
 * not be difficult to turn it into standalone code.
 *
 * The structure of the main routine is based on Christian Bau's earlier work.
 *
 * References:
 *  - Christian Bau's paper and example implementation, 2003, Christian Bau
 *    This was of immense help.  References to "step #" refer to this preprint.
 *  - "Computing Pi(x): the combinatorial method", 2006, Tomás Oliveira e Silva
 *  - "Computing Pi(x): The Meissel, Lehmer, Lagarias, Miller, Odlyzko Method"
 *    1996, Deléglise and Rivat.
 *
 * Comparisons to the other prime counting implementations in this package:
 *
 * Sieve:   Segmented, single threaded, thread-safe.  Small table enhanced,
 *          fastest for n < 60M.  Bad growth rate (like all sieves will have).
 * Legendre:Simple.  Recursive caching phi.
 * Meissel: Simple.  Non-recursive phi, lots of memory.
 * Lehmer:  Non-recursive phi, tries to restrict memory.
 * LMOS:    Simple.  Non-recursive phi, less memory than Lehmer above.
 * LMO:     Sieve phi.  Much faster and less memory than the others.
 *
 * Timing below is single core Haswell 4770K using Math::Prime::Util.
 *
 *  |   n   | Legendre |  Meissel |  Lehmer  |   LMOS   |    LMO    |
 *  +-------+----------+----------+----------+----------+-----------+
 *  | 10^19 |          |          |          |          | 2493.4    |
 *  | 10^18 |          |          |          |          |  498.16   |
 *  | 10^17 |10459.3   | 4348.3   | 6109.7   | 3478.0   |  103.03   |
 *  | 10^16 | 1354.6   |  510.8   |  758.6   |  458.4   |   21.64   |
 *  | 10^15 |  171.2   |   97.1   |  106.4   |   68.11  |    4.707  |
 *  | 10^14 |   23.56  |   18.59  |   16.51  |   10.44  |    1.032  |
 *  | 10^13 |    3.783 |    3.552 |    2.803 |    1.845 |    0.237  |
 *  | 10^12 |    0.755 |    0.697 |    0.505 |    0.378 |    54.9ms |
 *  | 10^11 |    0.165 |    0.144 |    93.7ms|    81.6ms|    13.80ms|
 *  | 10^10 |    35.9ms|    29.9ms|    19.9ms|    17.8ms|     3.64ms|
 *
 * Run with high memory limits: Meissel uses 1GB for 10^16, ~3GB for 10^17.
 * Lehmer is limited at high n values by sieving speed.  It is much faster
 * using parallel primesieve, though cannot come close to LMO.
 */

/* Below this size, just sieve (with table speedup). */
#define SIEVE_LIMIT  60000000
/* Adjust to get best performance.  Alpha from TOS paper. */
#define M_FACTOR(n)     (UV) ((double)n * (log(n)/log(5.2)) * (log(log(n))-1.4))
/* Size of segment used for previous primes, must be >= 21 */
#define PREV_SIEVE_SIZE 512
/* Phi sieve multiplier, adjust for best performance and memory use. */
#define PHI_SIEVE_MULT 13

#define FUNC_isqrt 1
#define FUNC_icbrt 1
#include "lmo.h"
#include "util.h"
#include "cache.h"
#include "sieve.h"

#ifdef _MSC_VER
  typedef unsigned __int8   uint8;
  typedef unsigned __int16  uint16;
  typedef unsigned __int32  uint32;
#else
  typedef unsigned char  uint8;
  typedef unsigned short uint16;
  typedef uint32_t       uint32;
#endif

/* UV is either uint32 or uint64 depending on Perl.  We use this native size
 * for the basic unit of the phi sieve.  It can be easily overridden here. */
typedef UV sword_t;
#define SWORD_BITS  BITS_PER_WORD
#define SWORD_ONES  UV_MAX
#define SWORD_MASKBIT(bits)  (UVCONST(1) << ((bits) % SWORD_BITS))
#define SWORD_CLEAR(s,bits)  s[bits/SWORD_BITS] &= ~SWORD_MASKBIT(bits)

/* GCC 3.4 - 4.1 has broken 64-bit popcount.
 * GCC 4.2+ can generate awful code when it doesn't have asm (GCC bug 36041).
 * When the asm is present (e.g. compile with -march=native on a platform that
 * has them, like Nahelem+), then it is almost as fast as the direct asm. */
#if SWORD_BITS == 64
 #if defined(__POPCNT__) && defined(__GNUC__) && (__GNUC__> 4 || (__GNUC__== 4 && __GNUC_MINOR__> 1))
   #define bitcount(b)  __builtin_popcountll(b)
 #else
   static sword_t bitcount(sword_t b) {
     b -= (b >> 1) & 0x5555555555555555;
     b = (b & 0x3333333333333333) + ((b >> 2) & 0x3333333333333333);
     b = (b + (b >> 4)) & 0x0f0f0f0f0f0f0f0f;
     return (b * 0x0101010101010101) >> 56;
   }
 #endif
#else
   /* An 8-bit table version is usually a little faster, but this is simpler. */
   static sword_t bitcount(sword_t b) {
     b -= (b >> 1) & 0x55555555;
     b = (b & 0x33333333) + ((b >> 2) & 0x33333333);
     b = (b + (b >> 4)) & 0x0f0f0f0f;
     return (b * 0x01010101) >> 24;
   }
#endif


/* Create array of small primes: 0,2,3,5,...,prev_prime(n+1) */
static uint32_t* make_primelist(uint32 n, uint32* number_of_primes)
{
  uint32 i = 0;
  uint32_t* plist;
  double logn = log(n);
  uint32 max_index = (n < 67)     ? 18
                   : (n < 355991) ? 15+(n/(logn-1.09))
                   : (n/logn) * (1.0+1.0/logn+2.51/(logn*logn));
  *number_of_primes = 0;
  New(0, plist, max_index+1, uint32_t);
  plist[0] = 0;
  /* We could do a simple SoE here.  This is not time critical. */
  START_DO_FOR_EACH_PRIME(2, n) {
    plist[++i] = p;
  } END_DO_FOR_EACH_PRIME;
  *number_of_primes = i;
  return plist;
}
#if 0  /* primesieve 5.0 example */
#include <primesieve.h>
static uint32_t* make_primelist(uint32 n, uint32* number_of_primes) {
  uint32_t plist;
  uint32_t* psprimes = generate_primes(2, n, number_of_primes, UINT_PRIMES);
  New(0, plist, *number_of_primes + 1, uint32_t);
  plist[0] = 0;
  memcpy(plist+1, psprimes, *number_of_primes * sizeof(uint32_t));
  primesieve_free(psprimes);
  return plist;
}
#endif

/* Given a max prime in small prime list, return max prev prime input */
static uint32 prev_sieve_max(UV maxprime) {
  UV limit = maxprime*maxprime - (maxprime*maxprime % (16*PREV_SIEVE_SIZE)) - 1;
  return (limit > U32_CONST(4294967295)) ? U32_CONST(4294967295) : limit;
}

/* Simple SoE filling a segment */
static void _prev_sieve_fill(UV start, uint8* sieve, const uint32_t* primes) {
  UV i, j, p;
  memset( sieve, 0xFF, PREV_SIEVE_SIZE );
  for (i = 2, p = 3; p*p < start + (16*PREV_SIEVE_SIZE); p = primes[++i])
    for (j = (start == 0) ? p*p/2 : (p-1) - ((start+(p-1))/2) % p;
         j < (8*PREV_SIEVE_SIZE); j += p)
      sieve[j/8] &= ~(1U << (j%8));
}

/* Calculate previous prime using small segment */
static uint32 prev_sieve_prime(uint32 n, uint8* sieve, uint32* segment_start, uint32 sieve_max, const uint32_t* primes)
{
  uint32 sieve_start, bit_offset;
  if (n <= 3) return (n == 3) ? 2 : 0;
  if (n > sieve_max) croak("ps overflow\n");

  /* If n > 3 && n <= sieve_max, then there is an odd prime we can find. */
  n -= 2;
  bit_offset = n % (16*PREV_SIEVE_SIZE);
  sieve_start = n - bit_offset;
  bit_offset >>= 1;

  while (1) {
    if (sieve_start != *segment_start) {  /* Fill sieve if necessary */
      _prev_sieve_fill(sieve_start, sieve, primes);
      *segment_start = sieve_start;
    }
    do {                                  /* Look for a set bit in sieve */
      if (sieve[bit_offset / 8] & (1u << (bit_offset % 8)))
        return sieve_start + 2*bit_offset + 1;
    } while (bit_offset-- > 0);
    sieve_start -= (16 * PREV_SIEVE_SIZE);
    bit_offset = ((16 * PREV_SIEVE_SIZE) - 1) / 2;
  }
}

/* Create factor table.
 * In lehmer.c we create mu and lpf arrays.  Here we use Christian Bau's
 * method, which is slightly more memory efficient and also a bit faster than
 * the code there (which does not use our fast ranged moebius).  It makes
 * very little difference -- mainly using this table is more convenient.
 *
 * In a uint16 we have stored:
 *    0     moebius(n) = 0
 *    even  moebius(n) = 1
 *    odd   moebius(n) = -1   (last bit indicates even/odd number of factors)
 *    v     smallest odd prime factor of n is v&1
 *    65535 large prime
 */
static uint16* ft_create(uint32 max)
{
  uint16* factor_table;
  uint32 i;
  uint32 tableLimit = max + 338 + 1;  /* At least one more prime */
  uint32 tableSize = tableLimit/2;
  uint32 max_prime = (tableLimit - 1) / 3 + 1;

  New(0, factor_table, tableSize, uint16);

  /* Set all values to 65535 (a large prime), set 0 to 65534. */
  factor_table[0] = 65534;
  for (i = 1; i < tableSize; ++i)
    factor_table[i] = 65535;

  /* Process each odd. */
  for (i = 1; i < tableSize; ++i) {
    uint32 factor, max_factor;
    uint32 p = i*2+1;
    if (factor_table[i] != 65535)  /* Already marked. */
      continue;
    if (p < 65535)                 /* p is a small prime, so set the number. */
      factor_table[i] = p;
    if (p >= max_prime)            /* No multiples will be in the table */
      continue;

    max_factor = (tableLimit - 1) / p + 1;
    /* Look for odd multiples of the prime p. */
    for (factor = 3; factor < max_factor; factor += 2) {
      uint32 index = (p*factor)/2;
      if (factor_table[index] == 65535)  /* p is smallest factor */
        factor_table[index] = p;
      else if (factor_table[index] > 0)  /* Change number of factors */
        factor_table[index] ^= 0x01;
    }

    /* Change all odd multiples of p*p to 0 to indicate non-square-free. */
    for (factor = p; factor < max_factor; factor += 2*p)
      factor_table[ (p*factor) / 2] = 0;
  }
  return factor_table;
}

#define PHIC 6

/* static const uint8_t _s0[ 1] = {0};
  static const uint8_t _s1[ 2] = {0,1};
  static const uint8_t _s2[ 6] = {0,1,1,1,1,2}; */
static const uint8_t _s3[30] = {0,1,1,1,1,1,1,2,2,2,2,3,3,4,4,4,4,5,5,6,6,6,6,7,7,7,7,7,7,8};
static const uint8_t _s4[210]= {0,1,1,1,1,1,1,1,1,1,1,2,2,3,3,3,3,4,4,5,5,5,5,6,6,6,6,6,6,7,7,8,8,8,8,8,8,9,9,9,9,10,10,11,11,11,11,12,12,12,12,12,12,13,13,13,13,13,13,14,14,15,15,15,15,15,15,16,16,16,16,17,17,18,18,18,18,18,18,19,19,19,19,20,20,20,20,20,20,21,21,21,21,21,21,21,21,22,22,22,22,23,23,24,24,24,24,25,25,26,26,26,26,27,27,27,27,27,27,27,27,28,28,28,28,28,28,29,29,29,29,30,30,30,30,30,30,31,31,32,32,32,32,33,33,33,33,33,33,34,34,35,35,35,35,35,35,36,36,36,36,36,36,37,37,37,37,38,38,39,39,39,39,40,40,40,40,40,40,41,41,42,42,42,42,42,42,43,43,43,43,44,44,45,45,45,45,46,46,47,47,47,47,47,47,47,47,47,47,48};
static UV tablephi(UV x, uint32 a)
{
  switch (a) {
    case 0: return x;
    case 1: return x-x/2;
    case 2: return x-x/2-x/3+x/6;
    case 3: return (x/    30U) *     8U + _s3[x %     30U];
    case 4: return (x/   210U) *    48U + _s4[x %    210U];
    case 5: {
              UV xp  = x / 11U;
              return ((x /210) * 48 + _s4[x  % 210]) -
                     ((xp/210) * 48 + _s4[xp % 210]);
             }
    case 6:
    default:{
              UV xp  = x / 11U;
              UV x2  = x / 13U;
              UV x2p = x2 / 11U;
              return ((x  /210) * 48 + _s4[x  % 210]) -
                     ((xp /210) * 48 + _s4[xp % 210]) -
                     ((x2 /210) * 48 + _s4[x2 % 210]) +
                     ((x2p/210) * 48 + _s4[x2p% 210]);
            }
    /* case 7: return tablephi(x,a-1)-tablephi(x/17,a-1); */  /* Hack hack */
  }
}

/****************************************************************************/
/*              Legendre Phi.  Not used by LMO, but exported.               */
/****************************************************************************/

/*
 * Choices include:
 *   1) recursive, memory-less.  We use this for small values.
 *   2) recursive, caching.  We use a this for larger values w/ 32MB cache.
 *   3) a-walker sorted list.  lehmer.c has this implementation.  It is
 *      faster for some values, but big and memory intensive.
 */
static UV _phi_recurse(UV x, UV a) {
  UV i, c = (a > PHIC) ? PHIC : a;
  UV sum = tablephi(x, c);
  if (a > c) {
    UV p  = nth_prime(c);
    UV pa = nth_prime(a);
    for (i = c+1; i <= a; i++) {
      UV xp;
      p = next_prime(p);
      xp = x/p;
      if (xp < p) {
        while (x < pa) {
          a--;
          pa = prev_prime(pa);
        }
        return (sum - a + i - 1);
      }
      sum -= legendre_phi(xp, i-1);
    }
  }
  return sum;
}

#define PHICACHEA 256
#define PHICACHEX 65536
#define PHICACHE_EXISTS(x,a) \
  ((x < PHICACHEX && a < PHICACHEA) ? cache[a*PHICACHEX+x] : 0)
static IV _phi(UV x, UV a, int sign, const uint32_t* const primes, const uint32_t lastidx, uint16_t* cache)
{
  IV sum;
  if      (PHICACHE_EXISTS(x,a))  return sign * cache[a*PHICACHEX+x];
  else if (a <= PHIC)             return sign * tablephi(x, a);
  else if (x < primes[a+1])       sum = sign;
  else {
    /* sum = _phi(x, a-1, sign, primes, lastidx, cache) +              */
    /*       _phi(x/primes[a], a-1, -sign, primes, lastidx, cache);    */
    UV a2, iters = (a*a > x)  ?  _XS_prime_count(2,isqrt(x))  :  a;
    UV c = (iters > PHIC) ? PHIC : iters;
    IV phixc = PHICACHE_EXISTS(x,c) ? cache[a*PHICACHEX+x] : tablephi(x,c);
    sum = sign * (iters - a + phixc);
    for (a2 = c+1; a2 <= iters; a2++)
      sum += _phi(x/primes[a2], a2-1, -sign, primes, lastidx, cache);
  }
  if (x < PHICACHEX && a < PHICACHEA && sign*sum <= SHRT_MAX)
    cache[a*PHICACHEX+x] = sign * sum;
  return sum;
}
UV legendre_phi(UV x, UV a)
{
  /* If 'x' is very small, give a quick answer with any 'a' */
  if (x <= PHIC)
    return tablephi(x, (a > PHIC) ? PHIC : a);

  /* Shortcuts for large values, from R. Andrew Ohana */
  if (a > (x >> 1))  return 1;
  /* If a > prime_count(2^32), then we need not be concerned with composite
   * x values with all factors > 2^32, as x is limited to 64-bit. */
  if (a > 203280221) {  /* prime_count(2**32) */
    UV pc = _XS_LMO_pi(x);
    return (a > pc)  ?  1  :  pc - a + 1;
  }
  /* If a is large enough, check the ratios */
  if (a > 1000000 && x < a*21) {  /* x always less than 2^32 */
    if ( _XS_LMO_pi(x) < a)  return 1;
  }

  /* TODO:  R. Andrew Ohana's 2011 SAGE code is faster as the a value
   * increases.  It uses a primelist as in the caching code below, as
   * well as a binary search prime count on it (like in our lehmer). */

  if ( a > 254 || (x > 1000000000 && a > 30) ) {
    uint16_t* cache;
    uint32_t* primes;
    uint32_t lastidx;
    UV res, max_cache_a = (a >= PHICACHEA) ? PHICACHEA : a+1;
    Newz(0, cache, PHICACHEX * max_cache_a, uint16_t);
    primes = make_primelist(nth_prime(a+1), &lastidx);
    res = (UV) _phi(x, a, 1, primes, lastidx, cache);
    Safefree(primes);
    Safefree(cache);
    return res;
  }

  return _phi_recurse(x, a);
}
/****************************************************************************/


typedef struct {
  sword_t  *sieve;                 /* segment bit mask */
  uint8    *word_count;            /* bit count in each 64-bit word */
  uint32   *word_count_sum;        /* cumulative sum of word_count */
  UV       *totals;                /* total bit count for all phis at index */
  uint32   *prime_index;           /* index of prime where phi(n/p/p(k+1))=1 */
  uint32   *first_bit_index;       /* offset relative to start for this prime */
  uint8    *multiplier;            /* mod-30 wheel of each prime */
  UV        start;                 /* x value of first bit of segment */
  UV        phi_total;             /* cumulative bit count before removal */
  uint32    size;                  /* segment size in bits */
  uint32    first_prime;           /* index of first prime in segment */
  uint32    last_prime;            /* index of last prime in segment */
  uint32    last_prime_to_remove;  /* index of last prime p, p^2 in segment */
} sieve_t;

/* Size of phi sieve in words.  Multiple of 3*5*7*11 words. */
#define PHI_SIEVE_WORDS (1155 * PHI_SIEVE_MULT)

/* Bit counting using cumulative sums.  A bit slower than using a running sum,
 * but a little simpler and can be run in parallel. */
static uint32 make_sieve_sums(uint32 sieve_size, const uint8* sieve_word_count, uint32* sieve_word_count_sum) {
  uint32 i, bc, words = (sieve_size + 2*SWORD_BITS-1) / (2*SWORD_BITS);
  sieve_word_count_sum[0] = 0;
  for (i = 0, bc = 0; i+7 < words; i += 8) {
    const uint8* cntptr = sieve_word_count + i;
    uint32* sumptr = sieve_word_count_sum + i;
    sumptr[1] = bc += cntptr[0];
    sumptr[2] = bc += cntptr[1];
    sumptr[3] = bc += cntptr[2];
    sumptr[4] = bc += cntptr[3];
    sumptr[5] = bc += cntptr[4];
    sumptr[6] = bc += cntptr[5];
    sumptr[7] = bc += cntptr[6];
    sumptr[8] = bc += cntptr[7];
  }
  for (; i < words; i++)
    sieve_word_count_sum[i+1] = sieve_word_count_sum[i] + sieve_word_count[i];
  return sieve_word_count_sum[words];
}

static UV _sieve_phi(UV segment_x, const sword_t* sieve, const uint32* sieve_word_count_sum) {
  uint32 bits = (segment_x + 1) / 2;
  uint32 words = bits / SWORD_BITS;
  uint32 sieve_sum = sieve_word_count_sum[words];
  sieve_sum += bitcount( sieve[words] & ~(SWORD_ONES << (bits % SWORD_BITS)) );
  return sieve_sum;
}

/* Erasing primes from the sieve is done using Christian Bau's
 * case statement walker.  It's not pretty, but it is short, fast,
 * clever, and does the job. */

#define sieve_zero(sieve, si, wordcount) \
  { uint32  index_ = si/SWORD_BITS; \
    sword_t mask_  = SWORD_MASKBIT(si); \
    if (sieve[index_] & mask_) { \
      sieve[index_] &= ~mask_; \
      wordcount[index_]--; \
    }  }

#define sieve_case_zero(casenum, skip, si, p, size, mult, sieve, wordcount) \
  case casenum: sieve_zero(sieve, si, wordcount); \
                si += skip * p; \
                mult = (casenum+1) % 8; \
                if (si >= size) break;

static void remove_primes(uint32 index, uint32 last_index, sieve_t* s, const uint32_t* primes)
{
  uint32    size = (s->size + 1) / 2;
  sword_t  *sieve = s->sieve;
  uint8    *word_count = s->word_count;

  s->phi_total = s->totals[last_index];
  for ( ;index <= last_index; index++) {
    if (index >= s->first_prime && index <= s->last_prime) {
      uint32 b = (primes[index] - (uint32) s->start - 1) / 2;
      sieve_zero(sieve, b, word_count);
    }
    if (index <= s->last_prime_to_remove) {
      uint32 b = s->first_bit_index[index];
      if (b < size) {
        uint32 p    = primes[index];
        uint32 mult = s->multiplier[index];
        switch (mult) {
          reloop: ;
            sieve_case_zero(0, 3, b, p, size, mult, sieve, word_count);
            sieve_case_zero(1, 2, b, p, size, mult, sieve, word_count);
            sieve_case_zero(2, 1, b, p, size, mult, sieve, word_count);
            sieve_case_zero(3, 2, b, p, size, mult, sieve, word_count);
            sieve_case_zero(4, 1, b, p, size, mult, sieve, word_count);
            sieve_case_zero(5, 2, b, p, size, mult, sieve, word_count);
            sieve_case_zero(6, 3, b, p, size, mult, sieve, word_count);
            sieve_case_zero(7, 1, b, p, size, mult, sieve, word_count);
          goto reloop;
        }
        s->multiplier[index] = mult;
      }
      s->first_bit_index[index] = b - size;
    }
  }
  s->totals[last_index] += make_sieve_sums(s->size, s->word_count, s->word_count_sum);
}

static void word_tile (sword_t* source, uint32 from, uint32 to) {
  while (from < to) {
    uint32 words = (2*from > to) ? to-from : from;
    memcpy(source+from, source, sizeof(sword_t)*words);
    from += words;
  }
}

static void init_segment(sieve_t* s, UV segment_start, uint32 size, uint32 start_prime_index, uint32 sieve_last, const uint32_t* primes)
{
  uint32    i, words;
  sword_t*  sieve = s->sieve;
  uint8*    word_count = s->word_count;

  s->start = segment_start;
  s->size  = size;

  if (segment_start == 0) {
    s->last_prime = 0;
    s->last_prime_to_remove = 0;
  }
  s->first_prime = s->last_prime + 1;
  while (s->last_prime < sieve_last) {
    uint32 p = primes[s->last_prime + 1];
    if (p >= segment_start + size)
      break;
    s->last_prime++;
  }
  while (s->last_prime_to_remove < sieve_last) {
    UV p = primes[s->last_prime_to_remove + 1];
    UV p2 = p*p;
    if (p2 >= segment_start + size)
      break;
    s->last_prime_to_remove++;
    s->first_bit_index[s->last_prime_to_remove] = (p2 - segment_start - 1) / 2;
    s->multiplier[s->last_prime_to_remove] = (uint8) ((p % 30) * 8 / 30);
  }

  memset(sieve, 0xFF, 3*sizeof(sword_t));  /* Set first 3 words to all 1 bits */
  if (start_prime_index >= 3)      /* Remove multiples of 3. */
    for (i = 3/2; i < 3 * SWORD_BITS; i += 3)
      SWORD_CLEAR(sieve, i);

  word_tile(sieve, 3, 15);         /* Copy to first 15 = 3*5 words */
  if (start_prime_index >= 3)      /* Remove multiples of 5. */
    for (i = 5/2; i < 15 * SWORD_BITS; i += 5)
      SWORD_CLEAR(sieve, i);

  word_tile(sieve, 15, 105);       /* Copy to first 105 = 3*5*7 words */
  if (start_prime_index >= 4)      /* Remove multiples of 7. */
    for (i = 7/2; i < 105 * SWORD_BITS; i += 7)
      SWORD_CLEAR(sieve, i);

  word_tile(sieve, 105, 1155);     /* Copy to first 1155 = 3*5*7*11 words */
  if (start_prime_index >= 5)      /* Remove multiples of 11. */
    for (i = 11/2; i < 1155 * SWORD_BITS; i += 11)
      SWORD_CLEAR(sieve, i);

  size = (size+1) / 2;             /* size to odds */
  words = (size + SWORD_BITS-1) / SWORD_BITS;   /* sieve size in words */
  word_tile(sieve, 1155, words);   /* Copy first 1155 words to rest */
  /* Zero all unused bits and words */
  if (size % SWORD_BITS)
    sieve[words-1] &= ~(SWORD_ONES << (size % SWORD_BITS));
  memset(sieve + words, 0x00, sizeof(sword_t)*(PHI_SIEVE_WORDS+2 - words));

  /* Create counts, remove primes (updating counts and sums). */
  for (i = 0; i < words; i++)
    word_count[i] = (uint8) bitcount(sieve[i]);
  remove_primes(6, start_prime_index, s, primes);
}

/* However we want to handle reduced prime counts */
#define simple_pi(n)  _XS_LMO_pi(n)
/* Macros to hide all the variables being passed */
#define prev_sieve_prime(n) \
  prev_sieve_prime(n, &prev_sieve[0], &ps_start, ps_max, primes)
#define sieve_phi(x) \
  ss.phi_total + _sieve_phi((x) - ss.start, ss.sieve, ss.word_count_sum)


UV _XS_LMO_pi(UV n)
{
  UV        N2, N3, K2, K3, M, sum1, sum2, phi_value;
  UV        sieve_start, sieve_end, least_divisor, step7_max, last_phi_sieve;
  uint32    j, k, piM, KM, end, prime, prime_index;
  uint32    ps_start, ps_max, smallest_divisor, nprimes;
  uint8     prev_sieve[PREV_SIEVE_SIZE];
  uint32_t *primes;
  uint16   *factor_table;
  sieve_t   ss;

  const uint32 c = PHIC;  /* We can use our fast function for this */

  /* For "small" n, use our table+segment sieve. */
  if (n < SIEVE_LIMIT || n < 10000)  return _XS_prime_count(2, n);
  /* n should now be reasonably sized (not tiny). */

  N2 = isqrt(n);             /* floor(N^1/2) */
  N3 = icbrt(n);             /* floor(N^1/3) */
  K2 = simple_pi(N2);        /* Pi(N2) */
  K3 = simple_pi(N3);        /* Pi(N3) */

  /* M is N^1/3 times a tunable performance factor. */
  M = (N3 > 500) ? M_FACTOR(N3) : N3+N3/2;
  if (M >= N2) M = N2 - 1;         /* M must be smaller than N^1/2 */
  if (M < N3) M = N3;              /* M must be at least N^1/3 */

  /* Create the array of small primes, and least-prime-factor/moebius table */
  primes = make_primelist( M + 500, &nprimes );
  factor_table = ft_create( M );

  /* Create other arrays */
  New(0, ss.sieve,           PHI_SIEVE_WORDS   + 2, sword_t);
  New(0, ss.word_count,      PHI_SIEVE_WORDS   + 2, uint8);
  New(0, ss.word_count_sum,  PHI_SIEVE_WORDS   + 2, uint32);
  New(0, ss.totals,          K3+2, UV);
  New(0, ss.prime_index,     K3+2, uint32);
  New(0, ss.first_bit_index, K3+2, uint32);
  New(0, ss.multiplier,      K3+2, uint8);

  if (ss.sieve == 0 || ss.word_count == 0 || ss.word_count_sum == 0 ||
      ss.totals == 0 || ss.prime_index == 0 || ss.first_bit_index == 0 ||
      ss.multiplier == 0)
    croak("Allocation failure in LMO Pi\n");

  /* Variables for fast prev_prime using small segment sieves (up to M^2) */
  ps_max   = prev_sieve_max( primes[nprimes] );
  ps_start = U32_CONST(0xFFFFFFFF);

  /* Look for the smallest divisor: the smallest number > M which is
   * square-free and not divisible by any prime covered by our Mapes
   * small-phi case.  The largest value we will look up in the phi
   * sieve is n/smallest_divisor. */
  for (j = (M+1)/2; factor_table[j] <= primes[c]; j++) /* */;
  smallest_divisor = 2*j+1;
  /* largest_divisor = (N2 > (UV)M * (UV)M)  ?  N2  :  (UV)M * (UV)M; */

  M = smallest_divisor - 1;   /* Increase M if possible */
  piM = simple_pi(M);
  if (piM < c)  croak("N too small for LMO\n");
  last_phi_sieve = n / smallest_divisor + 1;

  /* KM = smallest k, c <= k <= piM, s.t. primes[k+1] * primes[k+2] > M. */
  for (KM = c; primes[KM+1] * primes[KM+2] <= M && KM < piM; KM++) /* */;
  if (K3 < KM)  K3 = KM;  /* Ensure K3 >= KM */

  /* Start calculating Pi(n).  Steps 4-10 from Bau. */
  sum1 = (K2 - 1) + (UV) (piM - K3 - 1) * (UV) (piM - K3) / 2;
  sum2 = 0;
  end = (M+1)/2;

  /* Start at index K2, which is the prime preceeding N^1/2 */
  prime = prev_sieve_prime( (N2 >= ps_start) ? ps_start : N2+1 );
  prime_index = K2 - 1;
  step7_max = K3;

  /* Step 4:  For 1 <= x <= M where x is square-free and has no
   * factor <= primes[c], sum phi(n / x, c). */
  for (j = 0; j < end; j++) {
    uint32 lpf = factor_table[j];
    if (lpf > primes[c]) {
      phi_value = tablephi(n / (2*j+1), c);   /* x = 2j+1 */
      if (lpf & 0x01) sum2 += phi_value; else sum1 += phi_value;
    }
  }

  /* Step 5:  For 1+M/primes[c+1] <= x <= M, x square-free and
   * has no factor <= primes[c+1], sum phi(n / (x*primes[c+1]), c). */
  if (c < piM) {
    UV pc_1 = primes[c+1];
    for (j = (1+M/pc_1)/2; j < end; j++) {
      uint32 lpf = factor_table[j];
      if (lpf > pc_1) {
        phi_value = tablephi(n / (pc_1 * (2*j+1)), c);   /* x = 2j+1 */
        if (lpf & 0x01) sum1 += phi_value; else sum2 += phi_value;
      }
    }
  }

  for (k = 0; k <= K3; k++)     ss.totals[k] = 0;
  for (k = 0; k < KM; k++)      ss.prime_index[k] = end;

  /* Instead of dividing by all primes up to pi(M), once a divisor is large
   * enough then phi(n / (p*primes[k+1]), k) = 1. */
  {
    uint32 last_prime = piM;
    for (k = KM; k < K3; k++) {
      UV pk = primes[k+1];
      while (last_prime > k+1 && pk * pk * primes[last_prime] > n)
        last_prime--;
      ss.prime_index[k] = last_prime;
      sum1 += piM - last_prime;
    }
  }

  for (sieve_start = 0; sieve_start < last_phi_sieve; sieve_start = sieve_end) {
    /* This phi segment goes from sieve_start to sieve_end. */
    sieve_end = ((sieve_start + 2*SWORD_BITS*PHI_SIEVE_WORDS) <  last_phi_sieve)
              ?   sieve_start + 2*SWORD_BITS*PHI_SIEVE_WORDS  :  last_phi_sieve;
    /* Only divisors s.t. sieve_start <= N / divisor < sieve_end considered. */
    least_divisor = n / sieve_end;
    /* Initialize the sieve segment and all associated variables. */
    init_segment(&ss, sieve_start, sieve_end - sieve_start, c, K3, primes);

    /* Step 6:  For c < k < KM:  For 1+M/primes[k+1] <= x <= M, x square-free
     * and has no factor <= primes[k+1], sum phi(n / (x*primes[k+1]), k). */
    for (k = c+1; k < KM; k++) {
      UV pk = primes[k+1];
      uint32 start = (least_divisor >= pk * U32_CONST(0xFFFFFFFE))
                   ? U32_CONST(0xFFFFFFFF)
                   : (least_divisor / pk + 1)/2;
      remove_primes(k, k, &ss, primes);
      for (j = ss.prime_index[k] - 1; j >= start; j--) {
        uint32 lpf = factor_table[j];
        if (lpf > pk) {
          phi_value = sieve_phi(n / (pk * (2*j+1)));
          if (lpf & 0x01) sum1 += phi_value; else sum2 += phi_value;
        }
      }
      if (start < ss.prime_index[k])
        ss.prime_index[k] = start;
    }
    /* Step 7:  For KM <= K < Pi_M:  For primes[k+2] <= x <= M, sum
     * phi(n / (x*primes[k+1]), k).  The inner for loop can be parallelized. */
    for (; k < step7_max; k++) {
      remove_primes(k, k, &ss, primes);
      j = ss.prime_index[k];
      if (j >= k+2) {
        UV pk = primes[k+1];
        UV endj = j;
        while (endj > 7 && endj-7 >= k+2 && pk*primes[endj-7] > least_divisor) endj -= 8;
        while (            endj   >= k+2 && pk*primes[endj  ] > least_divisor) endj--;
        /* Now that we know how far to go, do the summations */
        for ( ; j > endj; j--)
          sum1 += sieve_phi(n / (pk*primes[j]));
        ss.prime_index[k] = endj;
      }
    }
    /* Restrict work for the above loop when we know it will be empty. */
    while (step7_max > KM && ss.prime_index[step7_max-1] < (step7_max-1)+2)
      step7_max--;

    /* Step 8:  For KM <= K < K3, sum -phi(n / primes[k+1], k) */
    remove_primes(k, K3, &ss, primes);
    /* Step 9:  For K3 <= k < K2, sum -phi(n / primes[k+1], k) + (k-K3). */
    while (prime > least_divisor && prime_index >= piM) {
      sum1 += prime_index - K3;
      sum2 += sieve_phi(n / prime);
      prime_index--;
      prime = prev_sieve_prime(prime);
    }
  }

  Safefree(ss.sieve);
  Safefree(ss.word_count);
  Safefree(ss.word_count_sum);
  Safefree(ss.totals);
  Safefree(ss.prime_index);
  Safefree(ss.first_bit_index);
  Safefree(ss.multiplier);
  Safefree(factor_table);
  Safefree(primes);

  return sum1 - sum2;
}
