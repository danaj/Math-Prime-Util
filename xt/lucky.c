/*
 * Lucky number sieves and utilities.  Standalone version 0.10, March 2023.
 * by Dana Jacobsen
 *
 * Includes small implementations of David Wilson's and Hugo van der Sanden's
 * algorithms from OEIS for comparison purposes.  The bitmask implementation
 * here is faster and uses less memory.
 *
 * Originally part of Math::Prime::Util available here:
 *     https://github.com/danaj/Math-Prime-Util
 *
 * The full pagelist implementation is not included here, as the bitmask is
 * superior in both memory and CPU time.  The small monolithic version is
 * included and is quite fast for very small inputs if counting microseconds,
 * but is inappropriate for large sieves.
 *
 * Timings for different sizes and algorithms on 2020 M1 Macbook.
 * Generate lucky numbers to    100k     1M       10M       100M      1000M
 *                              ----     -----  ------  --------   --------
 * bitmask126                   0.003    0.018    0.100     0.781     9.708
 * pagelist (monolithic)        0.002    0.137    6.824  1395.2
 * cgen (using count array)     0.005    0.133    5.741   344.04
 * wilson (iterated nth lucky)  0.024    0.879   54.069  3457.5
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
typedef unsigned long UV;
#define UVCONST(x) ((unsigned long)x##UL)
#define BMTYPE UV

#define croak(fmt,...)            { printf(fmt,##__VA_ARGS__); exit(3); }
#define New(id, mem, size, type)  mem = (type*) malloc((size)*sizeof(type))
#define Newz(id, mem, size, type) mem = (type*) calloc(size, sizeof(type))
#define Renew(mem, size, type)    mem =(type*)realloc(mem,(size)*sizeof(type))
#define Safefree(mem)             free((void*)mem)


/******************************************************************************/

#define MPU_MAX_LUCKY          UVCONST(18446744073709551517) /* Maybe */
static int _verbose = 0;

/******************************************************************************/
/*                         BITMASK126 DATA STRUCTURE                          */
/******************************************************************************/

/*
 * This is a bitmask for lucky numbers, using a 32-bit word for 126 integers.
 * Crucially, we use a tree of counts so we can skip to a given index in a
 * reasonable amount of time.
 *
 * The amount of memory used is about n/25.  This is about 20x smaller than
 * the 64-bit pagelist or cgen method, and 10x smaller than Wilson's list,
 * in addition to being much faster than cgen or Wilson.
 */

#define SSHIFT     4
#define TSHIFT     3

#define ADDSIZE(bm, wi, n) \
  { int _i; \
    bm->size[wi] += n; \
    bm->bsize[(wi) >> 3] += n; \
    bm->sbsize[(wi) >> (3+SSHIFT)] += n; \
    for (_i = 0; _i < bm->nilevels; _i++) \
      bm->tbsize[_i][(wi) >> (3+SSHIFT+(_i+1)*TSHIFT)] += n; \
  }

static unsigned char _bm_offset[32] = {1,3,7,9,13,15,21,25,31,33,37,43,45,49,51,55,63,67,69,73,75,79,85,87,93,97,99,105,109,111,115,117};
static unsigned char _bm_bit[63] = {0,1,1,2,3,3,4,5,5,5,6,6,7,7,7,8,9,9,10,10,10,11,12,12,13,14,14,15,15,15,15,16,16,17,18,18,19,20,20,21,21,21,22,23,23,23,24,24,25,26,26,26,27,27,28,29,29,30,31,31,31,31,31};

#define BM_WORD(n) (((n)>>1) / 63)
#define BM_BITN(n)  _bm_bit[(((n)>>1) % 63)]
#define BM_BITM(n)  (1U << _bm_bit[(((n)>>1) % 63)])

/* From Stanford Bit Twiddling Hacks, via "Nominal Animal" */
static uint32_t _nth_bit_set(uint32_t n, uint32_t word) {
  const uint32_t  pop2  = (word  & 0x55555555u) + ((word  >> 1) & 0x55555555u);
  const uint32_t  pop4  = (pop2  & 0x33333333u) + ((pop2  >> 2) & 0x33333333u);
  const uint32_t  pop8  = (pop4  & 0x0f0f0f0fu) + ((pop4  >> 4) & 0x0f0f0f0fu);
  const uint32_t  pop16 = (pop8  & 0x00ff00ffu) + ((pop8  >> 8) & 0x00ff00ffu);
  const uint32_t  pop32 = (pop16 & 0x000000ffu) + ((pop16 >>16) & 0x000000ffu);
  uint32_t        temp, rank = 0;

  if (n++ >= pop32)  return 32;

  temp = pop16 & 0xffu;
  if (n > temp) { n -= temp; rank += 16; }
  temp = (pop8 >> rank) & 0xffu;
  if (n > temp) { n -= temp; rank += 8; }
  temp = (pop4 >> rank) & 0x0fu;
  if (n > temp) { n -= temp; rank += 4; }
  temp = (pop2 >> rank) & 0x03u;
  if (n > temp) { n -= temp; rank += 2; }
  temp = (word >> rank) & 0x01u;
  if (n > temp) rank += 1;

  return rank;
}

typedef struct bitmask126_t {
  BMTYPE     n;
  BMTYPE     nelems;
  BMTYPE     nwords;
  int        nilevels;
  uint32_t*  data;
  uint8_t*   size;
  uint8_t*   bsize;
  uint16_t*  sbsize;
  BMTYPE*    tbsize[12];  /* Further index levels */
} bitmask126_t;

static bitmask126_t* bitmask126_create(BMTYPE n) {
  BMTYPE nblocks, nlevels;
  bitmask126_t *bm;
  New(0, bm, 1, bitmask126_t);

  bm->n = n;
  bm->nelems = 0;
  bm->nwords = (n+125)/126;
  nblocks = (bm->nwords + 7) / 8;
  Newz(0, bm->data, bm->nwords, uint32_t);
  Newz(0, bm->size, bm->nwords, uint8_t);
  Newz(0, bm->bsize, nblocks, uint8_t);
  nblocks = (nblocks + (1U << SSHIFT) - 1) >> SSHIFT;
  Newz(0, bm->sbsize, nblocks, uint16_t);

  for (nlevels=0;  nlevels < 12 && nblocks > 2*(1U<<TSHIFT);  nlevels++) {
    nblocks = (nblocks + (1U << TSHIFT) - 1) >> TSHIFT;
    Newz(0, bm->tbsize[nlevels], nblocks, BMTYPE);
  }
  bm->nilevels = nlevels;

  return bm;
}

static void bitmask126_destroy(bitmask126_t *bm) {
  int i;

  Safefree(bm->data);
  Safefree(bm->size);
  Safefree(bm->bsize);
  Safefree(bm->sbsize);
  for (i = 0; i < bm->nilevels; i++)
    Safefree(bm->tbsize[i]);
  bm->nelems = 0;
  bm->n = 0;
  Safefree(bm);
}

static void bitmask126_append(bitmask126_t *bm, BMTYPE n) {
  BMTYPE w = BM_WORD(n);
  bm->data[w] |= BM_BITM(n);
  ADDSIZE(bm, w, 1);
  bm->nelems++;
}

static BMTYPE* bitmask126_to_array(UV *size, bitmask126_t *bm) {
  BMTYPE nelem, wi, nwords, *arr;

  New(0, arr, bm->nelems, BMTYPE);
  nwords = bm->nwords;
  nelem = 0;

  for (wi = 0; wi < nwords; wi++) {
    uint32_t bit, w = bm->data[wi];
    for (bit = 0;  bit < 32;  bit++, w >>= 1)
      if (w & 1)
        arr[nelem++] = wi*126 + _bm_offset[bit];
  }
  if (nelem != bm->nelems) croak("bitmask126: bad number of elements in array");
  *size = nelem;
  return arr;
}
static uint32_t* bitmask126_to_array32(UV *size, bitmask126_t *bm) {
  uint32_t nelem, wi, nwords, *arr;

  New(0, arr, bm->nelems, uint32_t);
  nwords = bm->nwords;
  nelem = 0;

  for (wi = 0; wi < nwords; wi++) {
    uint32_t bit, w = bm->data[wi];
    for (bit = 0;  bit < 32;  bit++, w >>= 1)
      if (w & 1)
        arr[nelem++] = wi*126 + _bm_offset[bit];
  }
  if (nelem != bm->nelems) croak("bitmask126: bad number of elements in array");
  *size = nelem;
  return arr;
}


/* We want to find the e.g. 101'st set value, returns the array index wn. */
static BMTYPE _bitmask126_find_index(bitmask126_t *bm, BMTYPE *idx) {
  int lev;
  BMTYPE i = *idx, tbi, sbi, bi, wi;

  if (i > bm->nelems) croak("index higher than number of elements");

  /* Skip though superblock tree (128,2048,32768,524288,... words) */
  for (lev = bm->nilevels-1, tbi = 0;  lev >= 0;  lev--) {
    BMTYPE *tbsizei = bm->tbsize[lev];
    for (tbi = tbi << TSHIFT;  i >= tbsizei[tbi];  tbi++)
      i -= tbsizei[tbi];
  }
  for (sbi = tbi << TSHIFT;  i >= bm->sbsize[sbi];  sbi++)/* Skip superblocks */
    i -= bm->sbsize[sbi];
  for (bi = sbi << SSHIFT;  i >= bm->bsize[bi];  bi++)    /* Skip 8w blocks */
    i -= bm->bsize[bi];
  for (wi = bi << 3;  i >= bm->size[wi];  wi++)           /* Skip words */
    i -= bm->size[wi];

  *idx = i;
  return wi;
}

#if 0
static BMTYPE bitmask126_val(bitmask126_t *bm, BMTYPE idx) {
  BMTYPE wi;
  uint32_t bit;

  wi = _bitmask126_find_index(bm, &idx);
  bit = _nth_bit_set(idx, bm->data[wi]);
  return wi * 126 + _bm_offset[bit];
}
#endif

static void bitmask126_delete(bitmask126_t *bm, BMTYPE idx) { /* idx 0,1,... */
  BMTYPE wi = _bitmask126_find_index(bm, &idx);

  if (bm->size[wi] == 1) {   /* Only 1 value, zero the word. */
    bm->data[wi] = 0;
  } else {                   /* Find the index bit and zero it */
    uint32_t bit = _nth_bit_set(idx, bm->data[wi]);
    bm->data[wi] &= ~(1U << bit);
  }

  ADDSIZE(bm, wi, -1);
  bm->nelems--;
}



typedef struct bitmask126_iter_t {
  bitmask126_t *bm;
  uint32_t     *data;
  BMTYPE        wi;
  uint32_t      bit;
} bitmask126_iter_t;

static bitmask126_iter_t bitmask126_iterator_create(bitmask126_t *bm, BMTYPE idx) {
  bitmask126_iter_t iter;
  if (idx >= bm->nelems) croak("bitmask126: invalid iterator initial position\n");
  iter.bm = bm;
  iter.data = bm->data;
  iter.wi = _bitmask126_find_index(bm, &idx);
  iter.bit = _nth_bit_set(idx, bm->data[iter.wi]);
  return iter;
}

static BMTYPE bitmask126_iterator_next(bitmask126_iter_t *iter) {
  BMTYPE   v, wi  = iter->wi;
  uint32_t    bit = iter->bit;
  uint32_t    w   = iter->data[wi] >> bit;

  while (w == 0) {   /* skip any empty words */
    w = iter->data[++wi];
    bit = 0;
  }

#if defined(__GNUC__) && (__GNUC__ >= 4 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4))
  bit += __builtin_ctzl(w);
#else
  for ( ; bit < 32; bit++, w >>= 1)    /* Find next set bit */
    if (w & 1)
      break;
#endif

  v = wi * 126 + _bm_offset[bit];

  if (++bit > 31) {   /* Skip to next set bit */
    bit = 0;
    wi++;
  }
  iter->bit = bit;
  iter->wi  = wi;
  return v;
}

static BMTYPE bitmask126_iterator_prev(bitmask126_iter_t *iter) {
  BMTYPE   v, wi  = iter->wi;
  int         bit = iter->bit;
  uint32_t    w   = iter->data[wi];

  do {
    if (bit < 0) {
      if (wi == 0) croak("bitmask126: iterator underflow");
      w = iter->data[--wi];
      bit = 31;
    }
    for ( ; bit >= 0; bit--) {      /* Find prev set bit */
      if (w & (1U << bit))
        break;
    }
  } while (bit < 0);

  v = wi * 126 + _bm_offset[bit];

  if (bit > 0) {
    iter->bit = bit-1;
    iter->wi = wi;
  } else {
    iter->bit = 31;
    iter->wi = wi-1;
  }
  return v;
}

#undef BMTYPE
#undef SSHIFT
#undef TSHIFT
#undef ADDSIZE
#undef BM_WORD
#undef BM_BITN
#undef BM_BITM

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/






/******************************************************************************/
/*                             LUCKY NUMBERS                                  */
/******************************************************************************/

static const unsigned char _small_lucky[48] = {1,3,7,9,13,15,21,25,31,33,37,43,49,51,63,67,69,73,75,79,87,93,99,105,111,115,127,129,133,135,141,151,159,163,169,171,189,193,195,201,205,211,219,223,231,235,237,241};
static const unsigned char _small_lucky_count[48] = {0,1,1,2,2,2,2,3,3,4,4,4,4,5,5,6,6,6,6,6,6,7,7,7,7,8,8,8,8,8,8,9,9,10,10,10,10,11,11,11,11,11,11,12,12,12,12,12};
/* True for any position where (n % 7*9) could be a lucky number */
static const char _lmask63[63+2] = {1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0,0,1,1,0,1,1,0,0,0,0,1,1,0,1,1,0,1,1,0,0,0,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0,0,0,0,0,1,1};
/* mpufile '$n++; chomp; $v=$_; next unless $v > 10000; $m[ ($v>>1) % 4095 ]++; END { for (0..4094) { next unless $m[$_]; $b[$_ >> 5] |= (1 << ($_%32)); } say join ",",@b; }' ~/misc/ntheory/lucky_1e8.txt */
/* A large bitmask for ((n>>1) % 3*7*3*13) (819).  Covers 2,3,7,9,13. */
static const uint32_t _lmask5[26] = {2334495963,2261929142,1169344621,2204739155,2727961910,1639207725,3513561243,2430232978,1754683725,3630970059,3025873062,1278646881,3658323539,3055177010,1830209833,3406669457,3054200212,1837519692,1531293898,650340770,757258597,2606838995,2530306226,1169218145,3408442969,11572};


/* Lucky Number sieves.
 *
 * Mask presieving for the first 5 levels, followed by pre-sieving with a small
 * number of initial values.
 *
 * For fairly small sieves, less than 250k or so, we use a simplied pagelist.
 * Unlike the full pagelist method, this does not use an index tree.
 *
 * For sieving of non-small sizes, a bitmask (32 bits per 126 integers) is
 * used, with an index tree allowing log(n) time index lookups.  This is much
 * faster and uses substantially less memory than the other methods.  Memory
 * use grows linearly with the sieve size n.
 *
 * Generate first 10M lucky numbers (from 1 to 196502733) on 2020 M1 Mac:
 *           1.8s  bitmask126        memory:  n/25                (  8MB)
 *           3.1s  pagelist_sieve32  memory:  4 * count * ~2.5    (100MB)
 *           4.2s  pagelist_sieve64  memory:  8 * count * ~2.3    (190MB)
 *        1356s    lucky_cgen        memory:  8 * count * 2       (160MB)
 *        8950s    Wilson            memory:  8 * count * 1       ( 80MB)
 *
 * pagelist:
 * nth_lucky(1<<31):   55291335127    47 sec using lucky_sieve32  930MB
 * nth_lucky(1<<32):  113924214621   140 sec using lucky_sieve64  3.2GB
 * nth_lucky(1<<33):  234516370291   312 sec using lucky_sieve64  6.3GB
 * nth_lucky(1<<34):  482339741617   733 sec using lucky_sieve64 12.1GB
 *
 * bitmask:
 * nth_lucky(1<<31):   55291335127    23 sec using lucky_sieve32   89MB
 * nth_lucky(1<<32):  113924214621    50 sec using lucky_sieve64  173MB
 * nth_lucky(1<<33):  234516370291   107 sec using lucky_sieve64  341MB
 * nth_lucky(1<<34):  482339741617   224 sec using lucky_sieve64  675MB
 * nth_lucky(1<<35):  991238156013   469 sec using lucky_sieve64  1.3GB
 * nth_lucky(1<<36): 2035487409679   987 sec using lucky_sieve64  2.6GB
 * nth_lucky(1<<37): 4176793875529  2063 sec using lucky_sieve64  5.3GB
 *
 * A Graviton3 r7g takes about 1.6x more CPU time.
 * nth_lucky(1<<39)  17551419620869 in 258min on Graviton3 r7g, 21GB.
 * nth_lucky(1<<40)  35944896074391 in 523min on Graviton3 r7g, 42GB.
 * nth_lucky(1<<41)  73571139180453 in 1112min on Graviton3 r7g, 84GB.
 * nth_lucky(1<<42) 150499648533909 in 2303min on Graviton3 r7g, 168GB.
 * nth_lucky(1<<43) 307703784778627 in 3691min on Graviton3 r7g, 334GB.
 */


/* Simple 32-bit pagelist: fast for small (less than 10M or so) inputs.
 * Simple filtering, then sieve a big block using memmove.
 * This is memory intensive and has poor performance with large n.
 */
static uint32_t* _small_lucky_sieve32(UV *size, uint32_t n) {
  uint32_t i, m, c13, level, init_level, fsize, lsize, *lucky;

  if (n < 259) {
    if (n == 0) { *size = 0; return 0; }
    New(0, lucky, 5+n/5, uint32_t);
    for (lsize = 0; lsize < 48 && _small_lucky[lsize] <= n; lsize++)
      lucky[lsize] = _small_lucky[lsize];
    *size = lsize;
    return lucky;
  }

  /* @l=(2,3,7,9,13); $n=vecprod(@l); $n -= divint($n,$_) for @l; say $n */
  fsize = (uint64_t)1152*(n+4913)/4914;
  if (_verbose) { printf(" fsize %u\n", fsize); fflush(stdout); }
  New(0, lucky, 1 + fsize, uint32_t);
  lsize = c13 = 0;

  /* Create initial list, filtering out 3,7,9,13 */
  for (i = 1, m = 1; i <= n; i += 6) {
    if (_lmask63[m  ]) {
      if (++c13 == 13) c13 = 0;   else lucky[lsize++] = i;
    }
    if (_lmask63[m+2] && (i+2) <= n) {
      if (++c13 == 13) c13 = 0;   else lucky[lsize++] = i+2;
    }
    if ((m += 6) >= 63) m -= 63;
  }
  init_level = 5;
  if (_verbose) { printf(" finished initial list, lsize %u\n", lsize); fflush(stdout); }

  /* After the fill-in, we'll start deleting at 15 */
  for (level = init_level; level < lsize && lucky[level]-1 < lsize; level++) {
    uint32_t skip = lucky[level]-1, nlsize = skip;
    if (2*(skip+1) > lsize) break;  /* Only single skips left */
    for (i = skip+1; i < lsize; i += skip+1) {
      uint32_t ncopy = (skip <= (lsize-i)) ? skip : (lsize-i);
      memmove( lucky + nlsize, lucky + i, ncopy * sizeof(uint32_t) );
      nlsize += ncopy;
    }
    lsize = nlsize;
  }
  if (_verbose) { printf(" finished multi-skip deletes, level %u, lsize %u\n", level, lsize); fflush(stdout); }
  /* Now we just have single skips.  Process them all in one pass. */
  if (level < lsize && lucky[level]-1 < lsize) {
    uint32_t skip = lucky[level], nlsize = skip-1;
    while (skip < lsize) {
      uint32_t ncopy = lucky[level+1] - lucky[level];
      if (ncopy > lsize-skip)  ncopy = lsize - skip;
      memmove(lucky + nlsize, lucky + skip, ncopy * sizeof(uint32_t));
      nlsize += ncopy;
      skip += ncopy + 1;
      level++;
    }
    lsize = nlsize;
  }
  if (_verbose) { printf(" finished all deletes, lsize %u\n", lsize); fflush(stdout); }
  *size = lsize;
  return lucky;
}


static bitmask126_t* _bitmask126_sieve(UV* size, UV n) {
  UV i, lsize, level, init_level;
  bitmask126_t *pl;

  pl = bitmask126_create(n);

  {
    uint8_t count[48] = {0};
    uint32_t m, sln, ln, lbeg, lend;

    /* Decide how much additional filtering we'll do. */
    sln =  (n <=  200000000)  ?   21  :
           (n <= 0xFFFFFFFF)  ?   25  :   87;
    for (lbeg = lend = 5; lend < 48; lend++)
      if (_small_lucky[lend] >= sln)
        break;

    if (_verbose) { printf("bitmask lucky pre-sieve using %u lucky numbers up to %u\n", lend, _small_lucky[lend]); fflush(stdout); }

    /* Construct the initial list */
    for (i = 1, m = 0; i <= n; i += 2, m += 1) {
      if (m >= 819) m -= 819;  /* m = (i>>1) % 819 */
      if (_lmask5[m >> 5] & (1U << (m & 0x1F))) {
        for (ln = lbeg; ln <= lend; ln++) {
          if (++count[ln] == _small_lucky[ln]) {
            count[ln] = 0;
            break;
          }
        }
        if (ln > lend)
          bitmask126_append(pl,i);
      }
    }
    init_level = lend+1;
  }

  lsize = pl->nelems;
  if (_verbose) { printf("bitmask lucky done inserting.  values:  %lu\n",lsize); fflush(stdout); }

  if (init_level < lsize) {
    bitmask126_iter_t iter = bitmask126_iterator_create(pl, init_level);
    for (level = init_level; level < lsize; level++) {
      UV skip = bitmask126_iterator_next(&iter) - 1;
      if (skip >= lsize) break;
      for (i = skip; i < lsize; i += skip) {
        bitmask126_delete(pl, i);
        lsize--;
      }
    }
    if (_verbose) { printf("bitmask lucky done sieving.    values:  %lu\n",lsize); fflush(stdout); }
  }
  *size = lsize;
  return pl;
}

uint32_t* lucky_sieve32(UV *size, uint32_t n) {
  uint32_t *lucky;
  bitmask126_t *pl;

  if (n == 0) { *size = 0; return 0; }
  if (n > 4294967275U)  n = 4294967275U;  /* Max 32-bit lucky number */

  if (n <=   240000U) return _small_lucky_sieve32(size, n);

  pl = _bitmask126_sieve(size, n);

  lucky = bitmask126_to_array32(size, pl);
  if (_verbose) { printf("lucky_sieve32 done copying.\n"); fflush(stdout); }
  bitmask126_destroy(pl);
  return lucky;
}

UV* lucky_sieve64(UV *size, UV n) {
  UV *lucky;
  bitmask126_t *pl;

  if (n == 0) { *size = 0; return 0; }

  pl = _bitmask126_sieve(size, n);

  lucky = bitmask126_to_array(size, pl);
  if (_verbose) { printf("lucky_sieve64 done copying.\n"); fflush(stdout); }
  bitmask126_destroy(pl);
  return lucky;
}

UV* lucky_sieve_range(UV *size, UV beg, UV end) {
  UV i, nlucky, startcount, *lucky;
  bitmask126_t *pl;
  bitmask126_iter_t iter;

  if (end == 0 || beg > end) { *size = 0; return 0; }

  if (beg <= 1)  return lucky_sieve64(size, end);

  startcount = 0;  /* lucky_count_lower(beg) - 1; */
  pl = _bitmask126_sieve(size, end);
  New(0, lucky, *size - startcount, UV);
  iter = bitmask126_iterator_create(pl, startcount);
  for (i = startcount, nlucky = 0;  i < *size;  i++) {
    UV l = bitmask126_iterator_next(&iter);
    if (l >= beg)
      lucky[nlucky++] = l;
  }
  bitmask126_destroy(pl);
  *size = nlucky;
  return lucky;
}


/* Lucky Number sieve for 64-bit inputs.
 * Uses running counters to skip entries while we add them.
 * Based substantially on Hugo van der Sanden's cgen_lucky.c.
 */
UV* lucky_sieve_cgen(UV *size, UV n) {
  UV i, j, c3, lsize, lmax, lindex, *lucky, *count;

  if (n == 0) { *size = 0; return 0; }

  /* Init */
  lmax = (n < 1000) ? 153 : 100 + n/log(n);
  New(0, lucky, lmax, UV);
  New(0, count, lmax, UV);
  lucky[0] = 1;
  lucky[1] = 3;
  lucky[2] = 7;
  lindex = 2;
  lsize = 1;
  c3 = 2;

  for (i = 3; i <= n; i += 2) {
    if (!--c3) { c3 = 3; continue; }  /* Shortcut count[1] */
    for (j = 2; j < lindex; j++) {
      if (--count[j] == 0) {
        count[j] = lucky[j];
        break;
      }
    }
    if (j < lindex) continue;

    if (lsize >= lmax) {  /* Given the estimate, we probably never do this. */
      lmax = 1 + lsize * 1.2;
      Renew(lucky, lmax, UV);
      Renew(count, lmax, UV);
    }
    lucky[lsize] = count[lsize] = i;
    lsize++;

    if (lucky[lindex] == lsize) {
      lindex++;  lsize--;  /* Discard immediately */
    }
  }
  Safefree(count);
  *size = lsize;
  return lucky;
}

/* Lucky Number sieve based on David W. Wilson's generator */
UV* lucky_sieve_wilson(UV *size, UV n) {
  UV lmax, nlucky, g, k, i, *lucky;

  if (n == 0) { *size = 0; return 0; }
  lmax = (n <= 6)  ?  1+(n>2)  :  5 + 1.039 * n/log(n);
  New(0, lucky, lmax, UV);
  nlucky = 0;
  if (n >= 1)  lucky[nlucky++] = 1;
  if (n >= 3)  lucky[nlucky++] = 3;
  for (g = 0;  nlucky < lmax; nlucky++) {
    if (lucky[g+1] <= nlucky+1) g++;
    for (k = nlucky, i = g;  i >= 1;  i--)
      k = k * lucky[i] / (lucky[i]-1);
    k = 2*k;
    if (k+1 > n) break;
    lucky[nlucky] = k+1;
  }
  *size = nlucky;
  return lucky;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static UV _simple_lucky_count_upper(UV n) {
  double a, logn = log(n);
  if (n <=    6) return (n > 0) + (n > 2);
  if (n <= 7000) return 5 + 1.039 * n/logn;
  /* Don't make discontinities */
  a = (n < 10017000) ?   0.58003 - 3.00e-9 * (n-7000)   : 0.55;
  return n/(1.065*logn - a - 3.1/logn - 2.85/(logn*logn));
}
UV lucky_count_upper(UV n) {   /* Holds under 1e9 */
  if (n <       48) return _small_lucky_count[n];
  return _simple_lucky_count_upper(n);
}

UV lucky_count_range(UV lo, UV hi) {
  UV nlucky, lsize;

  if (hi < lo)
    return 0;
  if (hi < 48)
    return _small_lucky_count[hi] - (lo == 0 ? 0 : _small_lucky_count[lo-1]);

  /*
   * Analogous to how nth_lucky works, we sieve enough lucky numbers to
   * ensure we cover everything up to 'hi'.  We can then get an exact
   * count by determining exactly how many values will be removed.
   */

  if ((lo & 1)) lo--;    /* Both lo and hi will be even */
  if ((hi & 1)) hi++;
  lsize = 1+lucky_count_upper(hi);

  if (hi <= UVCONST(2000000000)) {
    uint32_t i, hicount = hi/2, locount = lo/2;
    uint32_t *lucky32 = lucky_sieve32(&nlucky, lsize);
    for (i = 1; i < nlucky && lucky32[i] <= lo; i++) {
      locount -= locount/lucky32[i];
      hicount -= hicount/lucky32[i];
    }
    for ( ; i < nlucky && lucky32[i] <= hicount; i++)
      hicount -= hicount/lucky32[i];
    Safefree(lucky32);
    return hicount - locount;
  } else {
    /* We use the iterator here to cut down on memory use. */
    UV i, hicount = hi/2, locount = lo/2;
    bitmask126_t* pl = _bitmask126_sieve(&nlucky, lsize);
    bitmask126_iter_t iter = bitmask126_iterator_create(pl, 1);
    for (i = 1; i < nlucky; i++) {
      UV l = bitmask126_iterator_next(&iter);
      if (l <= lo)  locount -= locount/l;
      if (l > hicount)  break;
      hicount -= hicount/l;
    }
    bitmask126_destroy(pl);
    return hicount - locount;
  }
}
UV lucky_count(UV n) {
  return lucky_count_range(0,n);
}

UV nth_lucky(UV n) {
  UV i, k, nlucky;

  if (n <= 48)  return (n == 0) ? 0 : _small_lucky[n-1];

  /* Apply the backward sieve, ala Wilson, for entry n */
  if (n <= UVCONST(100000000)) {
    uint32_t *lucky32 = lucky_sieve32(&nlucky, n);
    for (i = nlucky-1, k = n-1; i >= 1; i--)
      k += k/(lucky32[i]-1);
    Safefree(lucky32);
  } else { /* Iterate backwards through the sieve directly to save memory. */
    bitmask126_t* pl = _bitmask126_sieve(&nlucky, n);
    bitmask126_iter_t iter = bitmask126_iterator_create(pl, nlucky-1);
    for (i = nlucky-1, k = n-1; i >= 1; i--)
      k += k / (bitmask126_iterator_prev(&iter) - 1);
    bitmask126_destroy(pl);
  }
  return (2 * k + 1);
}


static int _test_lucky_to(UV lsize, UV *beg, UV *end) {
  UV i = *beg, pos = *end, l, quo, nlucky;
  int ret = -1;

  if (lsize <= 700000000U) {
    uint32_t *lucky32 = lucky_sieve32(&nlucky, lsize);
    while (i < nlucky) {
      l = lucky32[i++];
      if (pos < l)      { ret = 1; break; }
      quo = pos / l;
      if (pos == quo*l) { ret = 0; break; }
      pos -= quo;
    }
    Safefree(lucky32);
  } else {
    /* For 64-bit, iterate directly through the bit-mask to save memory. */
    bitmask126_t* pl = _bitmask126_sieve(&nlucky, lsize);
    if (i < nlucky) {
      bitmask126_iter_t iter = bitmask126_iterator_create(pl, i);
      while (i < nlucky) {
        l = bitmask126_iterator_next(&iter);
        i++;
        if (pos < l)      { ret = 1; break; }
        quo = pos / l;
        if (pos == quo*l) { ret = 0; break; }
        pos -= quo;
      }
    }
    bitmask126_destroy(pl);
  }
  /* printf("tested lsize = %lu  from %lu to %lu\n", lsize, *beg, i-1); */
  *beg = i;
  *end = pos;
  return ret;
}

int is_lucky(UV n) {
  UV i, l, quo, pos, lsize;
  int res;

  /* Simple pre-tests */
  if ( !(n & 1) || (n%6) == 5 || !_lmask63[n % 63]) return 0;
  if (n < 45) return 1;
  if (n > MPU_MAX_LUCKY) return 0;

  /* Check valid position using the static list */
  pos = (n+1) >> 1;  /* Initial position in odds */

  for (i = 1; i < 48; i++) {
    l = _small_lucky[i];
    if (pos < l) return 1;
    quo = pos / l;
    if (pos == quo*l) return 0;
    pos -= quo;
  }

  lsize = 1+lucky_count_upper(n);

  { /* Check more small values */
    UV psize = 600, gfac = 6;
    while (psize < lsize/3) {
      res = _test_lucky_to(psize, &i, &pos);
      if (res != -1) return res;
      psize *= gfac;
      gfac += 1;
    }
  }
  res = _test_lucky_to(lsize, &i, &pos);
  return (res == 0) ? 0 : 1;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static void _dieusage(const char* pname) {
  printf("lucky number sieve version 0.10.  Dana Jacobsen, 2023.\n\n");
  printf("Usage: %s [options] <n>\n\n", pname);
  printf("With no options, print lucky numbers up to <n> (inclusive)\n\n");
  printf("Options:\n");
  printf("  -v     set verbose (a few extra lines of information)\n");
  printf("  -q     set quiet (do work but don't print anything)\n");
  printf("  -help  this message\n");
  printf("\n");
  printf("  -t     test if <n> is a lucky number\n");
  printf("  -c     count of lucky numbers up to <n> (inclusive)\n");
  printf("  -n     the <n>-th lucky number\n");
  printf("\n");
  printf("  -s     sieve with counts (Hugo van der Sanden's method)\n");
  printf("  -w     compute using nth-lucky method (David Wilson's method)\n");
  printf("  -p     sieve with monolithic pagelist\n");
  printf("  -b     sieve with bitmask (the default)\n");
  exit(0);
}

#include <errno.h>
int main(int argc, char *argv[])
{
  UV i, n, nlucky;

  int a, optdone = 0;
  int flag_verbose = 0,  flag_quiet = 0;
  int flag_test = 0,  flag_count = 0,  flag_nth = 0;
  int flag_cgen = 0,  flag_wilson = 0,  flag_pagelist = 0,  flag_bitmask = 1;

  if (argc < 2) _dieusage(argv[0]);

  for (a = 1; a < argc; a++) {
    if (!optdone && argv[a][0] == '-') {
      if (strcmp(argv[a], "--") == 0) {
        optdone = 1;
      } else if (strcmp(argv[a], "-v") == 0) {
        flag_verbose = 1;
      } else if (strcmp(argv[a], "-q") == 0) {
        flag_quiet = 1;
      } else if (strcmp(argv[a], "-help") == 0 || strcmp(argv[a], "--help") == 0) {
        _dieusage(argv[0]);
      } else if (strcmp(argv[a], "-t") == 0) {
        flag_test = 1;  flag_count = 0;  flag_nth = 0;
      } else if (strcmp(argv[a], "-c") == 0) {
        flag_test = 0;  flag_count = 1;  flag_nth = 0;
      } else if (strcmp(argv[a], "-n") == 0) {
        flag_test = 0;  flag_count = 0;  flag_nth = 1;
      } else if (strcmp(argv[a], "-s") == 0) {
        flag_cgen = 1;  flag_wilson = 0;  flag_pagelist = 0;  flag_bitmask = 0;
      } else if (strcmp(argv[a], "-w") == 0) {
        flag_cgen = 0;  flag_wilson = 1;  flag_pagelist = 0;  flag_bitmask = 0;
      } else if (strcmp(argv[a], "-p") == 0) {
        flag_cgen = 0;  flag_wilson = 0;  flag_pagelist = 1;  flag_bitmask = 0;
      } else if (strcmp(argv[a], "-b") == 0) {
        flag_cgen = 0;  flag_wilson = 0;  flag_pagelist = 0;  flag_bitmask = 1;
      } else {
        printf("Unknown option: %s\n\n", argv[a]);
        _dieusage(argv[0]);
      }
      continue;
    }
    break;
  }
  n = strtoul(argv[a], 0, 10);
  if (n == ULONG_MAX && errno == ERANGE) { printf("Argument larger than ULONG_MAX\n"); return(-1); }

  if (flag_verbose) _verbose = 1;

  /************/

  if (flag_test) {
    int is = is_lucky(n);
    if (!flag_quiet) printf("%lu %s a lucky number\n", n, is ? "is" : "is not");
    return is;
  }
  if (flag_nth) {
    UV nth = nth_lucky(n);
    if (!flag_quiet) printf("%lu\n", nth);
    return 0;
  }
  if (flag_count) {
    UV count = lucky_count(n);
    if (!flag_quiet) printf("%lu\n", count);
    return 0;
  }

  /************/

  if (flag_wilson) {
    UV* lucky64 = lucky_sieve_wilson(&nlucky, n);
    if (!flag_quiet)
      for (i = 0; i < nlucky; i++)
        printf("%lu\n", lucky64[i]);
    Safefree(lucky64);
    return 0;
  }

  if (flag_cgen) {
    UV* lucky64 = lucky_sieve_cgen(&nlucky, n);
    if (!flag_quiet)
      for (i = 0; i < nlucky; i++)
        printf("%lu\n", lucky64[i]);
    Safefree(lucky64);
    return 0;
  }

  if (flag_pagelist) {
    uint32_t *lucky32;
    if (n > 0xFFFFFFFF) { printf("Small monolithic pagelist only supports 32-bit inputs\n"); exit(0); }
    lucky32 = _small_lucky_sieve32(&nlucky, n);
    if (!flag_quiet)
      for (i = 0; i < nlucky; i++)
        printf("%u\n", lucky32[i]);
    Safefree(lucky32);
    return 0;
  }

#if 0    /* Sieve to array.  Convenient. */
  if (n <= 0xFFFFFFFF) {
    uint32_t* lucky32 = lucky_sieve32(&nlucky, n);
    if (!flag_quiet)
      for (i = 0; i < nlucky; i++)
        printf("%u\n", lucky32[i]);
    Safefree(lucky32);
  } else {
    UV* lucky64 = lucky_sieve_range(&nlucky, 0, n);
    if (!flag_quiet)
      for (i = 0; i < nlucky; i++)
        printf("%lu\n", lucky64[i]);
    Safefree(lucky64);
  }
#else     /* Use the iterator */
  if (n > 0) {
    bitmask126_t* pl = _bitmask126_sieve(&nlucky, n);
    if (!flag_quiet) {
      bitmask126_iter_t iter = bitmask126_iterator_create(pl, 0);
      for (i = 0; i < nlucky; i++)
        printf("%lu\n", bitmask126_iterator_next(&iter));
    }
    bitmask126_destroy(pl);
  }
#endif

  return(0);
}

