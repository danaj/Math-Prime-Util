#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*****************************************************************************
 *
 * Prime counts using the extended Lagarias-Miller-Odlyzko combinatorial method.
 *
 * Copyright (c) 2013 Dana Jacobsen (dana@acm.org)
 * This is free software; you can redistribute it and/or modify it under
 * the same terms as the Perl 5 programming language system itself.
 *
 * This file is part of the Math::Prime::Util Perl module, but it should
 * not be difficult to turn it into standalone code.
 *
 * The algorithm details are based on Christian Bau's earlier work.
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
 * Legendre:Combinatorial phi.  Simple implementation.
 * Meissel: Combinatorial phi.  Simple implementation.
 * Lehmer:  Combinatorial phi.  Memory use grows rapidly.
 * LMOS:    Combinatorial phi.  Basic LMO implementation.
 * LMO:     Sieve phi.  10-50x faster than LMOS, better growth rate,
 *          Much, much better memory use than the others.
 *
 * Performance is slightly worse than Christian Bau's implementation, but
 * in theory this implementation has more parallelization opportunities.
 * Timing below is single core using MPU.
 *
 *  |   n   |  Meissel |  Lehmer  |   LMOS   |    LMO    |
 *  +-------+----------+----------+----------+-----------+
 *  | 10^19 |          |          |          | 3765.02   |
 *  | 10^18 |          |          |          |  778.21   |
 *  | 10^17 |          |          | 8844.2   |  163.77   |
 *  | 10^16 | 1410.2   | 1043.6   | 1058.9   |   34.81   |
 *  | 10^15 |  137.1   |  137.3   |  142.7   |    7.905  |
 *  | 10^14 |   26.18  |   21.74  |   21.29  |    1.726  |
 *  | 10^13 |    5.155 |    3.947 |    3.353 |    0.405  |
 *  | 10^12 |    1.059 |    0.700 |    0.626 |    0.0936 |
 *  | 10^11 |    0.227 |    0.138 |    0.124 |    0.0227 |
 *  | 10^10 |    0.0509|    0.0309|    0.0286|    0.00589|
 */

/* Below this size, just sieve (with table speedup). */
#define SIEVE_LIMIT  60000000
/* Adjust to get best performance.  Alpha from TOS paper. */
#define M_FACTOR(n)     (UV) ((double)n * (log(n)/log(5.2)) * (log(log(n))-1.4))
/* Size of segment used for previous primes, must be >= 21 */
#define PREV_SIEVE_SIZE 512
/* Phi sieve multiplier, adjust for best performance and memory use. */
#define PHI_SIEVE_MULT 13

#define FUNC_icbrt 1
#include "lmo.h"
#include "util.h"
#include "cache.h"
#include "sieve.h"

typedef unsigned char  uint8;
typedef unsigned short uint16;
typedef uint32_t       uint32;

static const unsigned char byte_ones[256] =
  {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};

static uint32_t bitcount(uint32_t b) {
  /* The simple table method is faster than __builtin_popcount or
   * using a 11-bit word table.  It is slower than the Nahalem asm. */
#if 1
  return byte_ones[b&0xFF] + byte_ones[(b>>8)&0xFF] + byte_ones[(b>>16)&0xFF] + byte_ones[b>>24];
#else
  uint32_t ret;
  __asm__("popcnt %1, %0" : "=r" (ret) : "r" (b));
  return ret;
#endif
}

/* static uint32 count_one_bits(const unsigned char* m, uint32 nbytes) {
  uint32 count = 0;
  while (nbytes--)  count += byte_ones[*m++];
  return count;
} */

/* Create array of small primes: 0,2,3,5,...,prev_prime(n+1) */
static uint32_t* make_primelist(uint32 n, uint32* number_of_primes)
{
  uint32 i = 0;
  uint32_t* plist;
  double logn = log(n);
  uint32 max_index = (n < 67)     ? 18
                   : (n < 355991) ? 15+(n/(logn-1.09))
                   : (n/log(n)) * (1.0+1.0/logn+2.51/(logn*logn));
  *number_of_primes = 0;
  New(0, plist, max_index+1, uint32_t);
  if (plist == 0)
    croak("Can not allocate small primes\n");
  plist[0] = 0;
  START_DO_FOR_EACH_PRIME(2, n) {
    plist[++i] = p;
  } END_DO_FOR_EACH_PRIME;
  *number_of_primes = i;
  return plist;
}

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
  if (factor_table == 0)  return 0;

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


static UV mapes(UV x, uint32 a)
{
  UV val;
  if (a == 0)  return x;
  if (a == 1)  return x-x/2;
  val = x-x/2-x/3+x/6;
  if (a >= 3) val += 0-x/5+x/10+x/15-x/30;
  if (a >= 4) val += 0-x/7+x/14+x/21-x/42+x/35-x/70-x/105+x/210;
  if (a >= 5) val += 0-x/11+x/22+x/33-x/66+x/55-x/110-x/165+x/330+x/77-x/154-x/231+x/462-x/385+x/770+x/1155-x/2310;
  if (a >= 6) val += 0-x/13+x/26+x/39-x/78+x/65-x/130-x/195+x/390+x/91-x/182-x/273+x/546-x/455+x/910+x/1365-x/2730+x/143-x/286-x/429+x/858-x/715+x/1430+x/2145-x/4290-x/1001+x/2002+x/3003-x/6006+x/5005-x/10010-x/15015+x/30030;
  if (a >= 7) val += 0-x/17+x/34+x/51-x/102+x/85-x/170-x/255+x/510+x/119-x/238-x/357+x/714-x/595+x/1190+x/1785-x/3570+x/187-x/374-x/561+x/1122-x/935+x/1870+x/2805-x/5610-x/1309+x/2618+x/3927-x/7854+x/6545-x/13090-x/19635+x/39270+x/221-x/442-x/663+x/1326-x/1105+x/2210+x/3315-x/6630-x/1547+x/3094+x/4641-x/9282+x/7735-x/15470-x/23205+x/46410-x/2431+x/4862+x/7293-x/14586+x/12155-x/24310-x/36465+x/72930+x/17017-x/34034-x/51051+x/102102-x/85085+x/170170+x/255255-x/510510;
  return (UV) val;
}

typedef struct {
  uint32_t *sieve;                 /* segment bit mask */
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

/* Size of phi sieve in 32-bit words.  Multiple of 2*2*3*5*7*11 bytes. */
#define PHI_SIEVE_WORDS (1155 * PHI_SIEVE_MULT)

/* Bit counting using cumulative sums.  A bit slower than using a running sum,
 * but a little simpler and can be run in parallel. */
static void make_sieve_counts(const uint32_t* sieve, uint32 sieve_size, uint8* sieve_word_count) {
  uint32 lwords = (sieve_size + 127) / 128;
  while (lwords--) {
    *sieve_word_count++ = bitcount(sieve[0]) + bitcount(sieve[1]);
    sieve += 2;
  }
}
static uint32 make_sieve_sums(uint32 sieve_size, const uint8* sieve_word_count, uint32* sieve_word_count_sum) {
  uint32 i, bc, lwords = (sieve_size + 127) / 128;
  sieve_word_count_sum[0] = 0;
  for (i = 0, bc = 0; i+7 < lwords; i += 8) {
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
  for (; i < lwords; i++)
    sieve_word_count_sum[i+1] = sieve_word_count_sum[i] + sieve_word_count[i];
  return sieve_word_count_sum[lwords];
}

static UV _sieve_phi(UV segment_x, const uint32_t* sieve, const uint32* sieve_word_count_sum) {
  uint32 bits = (segment_x + 1) / 2;
  uint32 words = bits / 32;
  uint32 sieve_sum = sieve_word_count_sum[words/2];
  if (words & 1) sieve_sum += bitcount( sieve[words-1] );
  sieve_sum += bitcount( sieve[words] & ~(0xfffffffful << (bits % 32)) );
  return sieve_sum;
}

/* Erasing primes from the sieve is done using Christian Bau's
 * case statement walker.  It's not pretty, but it is short, fast,
 * clever, and does the job. */

#define sieve_zero(sieve, si, wordcount) \
  { uint32 index = si/32; \
    uint32 mask = 1UL << (si % 32); \
    if (sieve[index] & mask) { \
      sieve[index] &= ~mask; \
      wordcount[index/2]--; \
    }  }

#define sieve_case_zero(casenum, skip, si, p, size, mult, sieve, wordcount) \
  case casenum: sieve_zero(sieve, si, wordcount); \
                si += skip * p; \
                mult = (casenum+1) % 8; \
                if (si >= size) break;

static void remove_primes(uint32 index, uint32 last_index, sieve_t* s, const uint32_t* primes)
{
  uint32    size = (s->size + 1) / 2;
  uint32_t *sieve = s->sieve;
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

static void word_tile (uint32_t* source, uint32 from, uint32 to) {
  while (from < to) {
    uint32 words = (2*from > to) ? to-from : from;
    memcpy(source+from, source, 4*words);
    from += words;
  }
}

static void init_segment(sieve_t* s, UV segment_start, uint32 size, uint32 start_prime_index, uint32 sieve_last, const uint32_t* primes)
{
  uint32    i, words;
  uint32_t* sieve = s->sieve;
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
    s->multiplier[s->last_prime_to_remove] = (p % 30) * 8 / 30;
  }

  memset(sieve, 0xFF, 3*4);        /* Set first 3 words to all 1 bits */
  if (start_prime_index >= 3)      /* Remove multiples of 3. */
    for (i = 3/2; i < 3 * 32; i += 3)
      sieve[i / 32] &= ~(1ul << (i % 32));

  word_tile(sieve, 3, 15);         /* Copy to first 15 = 3*5 words */
  if (start_prime_index >= 3)      /* Remove multiples of 5. */
    for (i = 5/2; i < 15 * 32; i += 5)
      sieve[i / 32] &= ~(1ul << (i % 32));

  word_tile(sieve, 15, 105);       /* Copy to first 105 = 3*5*7 words */
  if (start_prime_index >= 4)      /* Remove multiples of 7. */
    for (i = 7/2; i < 105 * 32; i += 7)
      sieve[i / 32] &= ~(1ul << (i % 32));

  word_tile(sieve, 105, 1155);     /* Copy to first 1155 = 3*5*7*11 words */
  if (start_prime_index >= 5)      /* Remove multiples of 11. */
    for (i = 11/2; i < 1155 * 32; i += 11)
      sieve[i / 32] &= ~(1ul << (i % 32));

  size = (size+1) / 2;             /* size to odds */
  words = (size + 31) / 32;        /* sieve size in 32-bit words */
  word_tile(sieve, 1155, words);   /* Copy first 1155 words to rest */
  /* Zero all unused bits and words */
  if (size % 32)
    sieve[words-1] &= ~ (0xffffffffUL << (size % 32));
  memset(sieve + words, 0x00, 4*(PHI_SIEVE_WORDS+2-words));

  /* Create counts, remove primes (updating counts and sums). */
  make_sieve_counts(sieve, s->size, word_count);
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

  const uint32 c = 7;  /* We can use out Mapes function for c <= 7 */

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
  if (primes == 0 || factor_table == 0)
    croak("Allocation failure in LMO Pi\n");

  /* Create other arrays */
  New(0, ss.sieve,           PHI_SIEVE_WORDS   + 2, uint32_t);
  New(0, ss.word_count,      PHI_SIEVE_WORDS/2 + 2, uint8);
  New(0, ss.word_count_sum,  PHI_SIEVE_WORDS/2 + 2, uint32);
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
      phi_value = mapes(n / (2*j+1), c);   /* x = 2j+1 */
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
        phi_value = mapes(n / (pc_1 * (2*j+1)), c);   /* x = 2j+1 */
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
    sieve_end = ((sieve_start + 64*PHI_SIEVE_WORDS) <  last_phi_sieve)
              ?   sieve_start + 64*PHI_SIEVE_WORDS  :  last_phi_sieve;
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
      if (ss.prime_index[k] >= k+2) {
        UV pk = primes[k+1];
        UV endj = ss.prime_index[k];
        while (endj > 7 && endj-7 >= k+2 && pk*primes[endj-7] > least_divisor) endj -= 8;
        while (            endj   >= k+2 && pk*primes[endj  ] > least_divisor) endj--;
        /* Now that we know how far to go, do the summations */
        for (j = ss.prime_index[k]; j > endj; j--)
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
