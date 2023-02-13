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
#define SSHIFT     3
#define TSHIFT     4

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

static uint32_t _nth_bit_set(uint32_t n, uint32_t word) {
  uint32_t bit;
  for (bit = 0;  bit < 32;  bit++, word >>= 1)
    if (word & 1)
      if (n-- == 0)
        break;
  if (bit > 31) croak("bitmask126: bad nth_bit");
  return bit;
}

typedef struct bitmask126_t {
  BMTYPE     n;
  BMTYPE     nelems;
  BMTYPE     nwords;
  BMTYPE     nilevels;
  uint32_t*  data;
  uint8_t*   size;
  uint8_t*   bsize;
  uint16_t*  sbsize;
  BMTYPE*    tbsize[8];  /* Further index levels */
} bitmask126_t;

static bitmask126_t* bitmask126_create(BMTYPE n) {
  BMTYPE i, nblocks, nlevels;
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

  for (i=0, nlevels=0;  nlevels < 8 && nblocks > 2*(1U<<TSHIFT);  nlevels++) {
    nblocks = (nblocks + (1U << TSHIFT) - 1) >> TSHIFT;
#if BMDEBUG
    printf("    level %lu blocks = %lu\n", nlevels, nblocks);
#endif
    Newz(0, bm->tbsize[nlevels], nblocks, BMTYPE);
  }
  bm->nilevels = nlevels;

  return bm;
}

static void bitmask126_destroy(bitmask126_t *bm) {
  BMTYPE i;

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
#if BMDEBUG
  if (n >= bm->n) croak("bitmask126: bad n in append");
#endif
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

static BMTYPE bitmask126_val(bitmask126_t *bm, BMTYPE idx) {
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
  bm->nelems--;
}



/* TODO:  The iterator might be faster using a function that gets all values
 * from a word and inserts into an array.  The iterator could hand these off
 * and refill as needed.
 */

typedef struct bitmask126_iter_t {
  bitmask126_t *bm;
  uint32_t     *data;
  BMTYPE        wi;
  uint32_t      bit;
} bitmask126_iter_t;

static bitmask126_iter_t bitmask126_iterator_create(bitmask126_t *bm, BMTYPE idx) {
  bitmask126_iter_t iter;
  if (idx >= bm->nelems) croak("bitmask 126: invalid iterator initial position\n");
  iter.bm = bm;
  iter.data = bm->data;
  iter.wi = _bitmask126_find_index(bm, &idx);
  iter.bit = _nth_bit_set(idx, bm->data[iter.wi]);
  return iter;
}

static BMTYPE bitmask126_iterator_next(bitmask126_iter_t *iter) {
  BMTYPE v;
  uint32_t w = iter->data[iter->wi] >> iter->bit;

  while (w == 0) {   /* skip any empty words */
    w = iter->bm->data[++iter->wi];
    iter->bit = 0;
  }

  for ( ; iter->bit < 32; iter->bit++, w >>= 1)    /* Find next set bit */
    if (w & 1)
      break;
  if (iter->bit > 31) croak("bitmask126: iterator bad nth_bit");

  v = iter->wi * 126 + _bm_offset[iter->bit];

  iter->bit++;              /* Skip to next set bit */
  if (iter->bit > 31) {
    iter->bit = 0;
    iter->wi++;
  }
  return v;
}

static BMTYPE bitmask126_iterator_prev(bitmask126_iter_t *iter) {
  BMTYPE v;
  uint32_t w;
  int bit = iter->bit;

  if (iter->wi >= iter->bm->nwords) croak("bitmask126: iterator underflow");
  w = iter->data[iter->wi];

  do {
    if (bit < 0) {
      if (iter->wi == 0) croak("bitmask126: iterator underflow");
      w = iter->data[--iter->wi];
      bit = 31;
    }
    for ( ; bit >= 0; bit--) {      /* Find prev set bit */
      if (w & (1U << bit))
        break;
    }
  } while (bit < 0);

  v = iter->wi * 126 + _bm_offset[bit];

  if (bit > 0) {
    iter->bit = bit-1;
  } else {
    iter->bit = 31;
    iter->wi--;
  }
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
