#ifndef MPU_DS_BITMASK126_H
#define MPU_DS_BITMASK126_H

#include "ptypes.h"

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
 * in addition to being orders of magnitude faster than cgen or Wilson.
 */

#ifndef BMTYPE
#define BMTYPE UV
#endif

#define BMDEBUG    0
/*  Not clear if SSHIFT/TSHIFT should be 3/3, 3/4, 4/3, or 4/4  */
#define SSHIFT     4
#define TSHIFT     3

static unsigned char _bm_offset[32] = {
   1, 3, 7, 9,13,15,21,25,31,33,37, 43, 45, 49, 51, 55,
  63,67,69,73,75,79,85,87,93,97,99,105,109,111,115,117};
static unsigned char _bm_bit[63] = {
   0, 1, 1, 2, 3, 3, 4, 5, 5, 5, 6, 6, 7, 7, 7, 8,
   9, 9,10,10,10,11,12,12,13,14,14,15,15,15,15,16,
  16,17,18,18,19,20,20,21,21,21,22,23,23,23,24,24,
  25,26,26,26,27,27,28,29,29,30,31,31,31,31,31   };

#define BM_WORD(n)  ((n) / 2 / 63)
#define BM_BITN(n)  _bm_bit[(n) / 2 % 63]
#define BM_BITM(n)  (1U << BM_BITN(n))

/* Modified from Stanford Bit Twiddling Hacks, via "Nominal Animal" */
static uint32_t _nth_bit_set(uint32_t n, uint32_t word) {
  const uint32_t  pop2  =  word                - (word >> 1 & 0x55555555u);
  const uint32_t  pop4  = (pop2 & 0x33333333u) + (pop2 >> 2 & 0x33333333u);
  const uint32_t  pop8  =  pop4  + (pop4  >>  4) & 0x0f0f0f0fu;
  const uint32_t  pop16 =  pop8  + (pop8  >>  8) & 0x00ff00ffu;
  const uint32_t  pop32 =  pop16 + (pop16 >> 16) & 0x000000ffu;
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

/* We are trying to balance space and performance. */
/* Also note that we do not support a bitmask with 256 consecutive bits set,
 * as that would overflow bsize.
 * A "block" could be changed to 16 or 32 words * with a uint16_t bsize. */
typedef struct bitmask126_t {
  BMTYPE     n;           /* The upper limit on the sieve */
  BMTYPE     nelems;      /* The total number of bits set */
  BMTYPE     nwords;      /* The number of words in data[] */
  int        nilevels;    /* The number of tbsize[] arrays actually used */
  uint32_t*  data;        /* The bitmap itself */
  uint8_t*   size;        /* The number of bits set in each data[] word */
  uint8_t*   bsize;       /* Sums over 8-word blocks */
  uint16_t*  sbsize;      /* Sums over 1<<SSHIFT block superblocks */
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

  for (nlevels=0;  nlevels < 12;  nlevels++) {
    nblocks = (nblocks + (1U << TSHIFT) - 1) >> TSHIFT;
    if (nblocks < 3) break;
#if BMDEBUG
    printf("    level %lu blocks = %lu\n", nlevels, nblocks);
#endif
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

/* Update all index levels for adding (subtracting) n bits to word wi. */
#define ADDSIZE(bm, wi, n) \
  do { int _i; \
       BMTYPE _j = wi; \
       bm->size[_j] += n; \
       bm->bsize[_j >>= 3] += n; \
       bm->sbsize[_j >>= SSHIFT] += n; \
       for (_i = 0; _i < bm->nilevels; _i++) \
         bm->tbsize[_i][_j >>= TSHIFT] += n; \
       bm->nelems += n; \
  } while (0)

static void bitmask126_append(bitmask126_t *bm, BMTYPE n) {
  BMTYPE w = BM_WORD(n);
#if BMDEBUG
  if (n >= bm->n) croak("bitmask126: bad n in append");
#endif
  bm->data[w] |= BM_BITM(n);
  ADDSIZE(bm, w, 1);
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
static uint32_t* bitmask126_to_array32(UV *size, const bitmask126_t *bm) {
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


/* We want to find the e.g. 101'st set value, returns the array word index,
 * and set *idx to the number of bits to skip within that word.            */
static BMTYPE _bitmask126_find_index(const bitmask126_t *bm, BMTYPE *idx) {
  int lev;
  BMTYPE i = *idx, j = 0;

  if (i > bm->nelems) croak("index higher than number of elements");

  /* Skip though superblock tree (128,2048,32768,524288,... words) */
  for (lev = bm->nilevels-1;  lev >= 0;  lev--) {
    const BMTYPE *tbsizei = bm->tbsize[lev];
    for (j <<= TSHIFT;  i >= tbsizei[j];  j++)
      i -= tbsizei[j];
  }
  for (j <<= TSHIFT;  i >= bm->sbsize[j];  j++) /* Skip superblocks */
    i -= bm->sbsize[j];
  for (j <<= SSHIFT;  i >= bm->bsize[j];   j++) /* Skip 8w blocks */
    i -= bm->bsize[j];
  for (j <<= 3;       i >= bm->size[j];    j++) /* Skip words */
    i -= bm->size[j];

  *idx = i;
  return j;
}

static INLINE BMTYPE bitmask126_val(const bitmask126_t *bm, BMTYPE idx) {
  BMTYPE wi;
  uint32_t bit;

  wi = _bitmask126_find_index(bm, &idx);
  bit = _nth_bit_set(idx, bm->data[wi]);
  return wi * 126 + _bm_offset[bit];
}

static void bitmask126_delete(bitmask126_t *bm, BMTYPE idx) { /* idx 0,1,... */
  BMTYPE wi;

#if BMDEBUG
  if (idx >= bm->nelems) croak("bitmask126: bad index in delete");
#endif

  wi = _bitmask126_find_index(bm, &idx);

  if (bm->size[wi] == 1) {   /* Only 1 value, zero the word. */
    bm->data[wi] = 0;
  } else {                   /* Find the index bit and zero it */
    uint32_t bit = _nth_bit_set(idx, bm->data[wi]);
    bm->data[wi] &= ~(1U << bit);
  }

  ADDSIZE(bm, wi, -1);
}



typedef struct bitmask126_iter_t {
  const bitmask126_t *bm;
  const uint32_t     *data;
  BMTYPE              wi;
  uint32_t            bit;
} bitmask126_iter_t;

static bitmask126_iter_t bitmask126_iterator_create(const bitmask126_t *bm, BMTYPE idx) {
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

#if defined(__GNUC__) && 100*__GNUC__ + __GNUC_MINOR >= 304
  bit += __builtin_ctzl(w);
#else
  for ( ; bit < 32; bit++, w >>= 1)    /* Find next set bit */
    if (w & 1)
      break;
#endif

  v = wi * 126 + _bm_offset[bit];

  iter->bit = ++bit & 31;
  iter->wi  = wi + (bit>>5);
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
      if (w & 1U << bit)
        break;
    }
  } while (bit < 0);

  v = wi * 126 + _bm_offset[bit];

  iter->bit = --bit & 31;
  iter->wi = wi - (bit >> 5 & 1);
  return v;
}

#undef BMTYPE
#undef BMDEBUG
#undef SSHIFT
#undef TSHIFT
#undef ADDSIZE
#undef BM_WORD
#undef BM_BITN
#undef BM_BITM

#endif
