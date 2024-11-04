#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "ds_iset.h"
#include "util.h"

#define FILL_RATIO 0.50

#if BITS_PER_WORD == 32
/* 16 0x45d9f3b  16 0x45d9f3b  16 */
/* 16 0x21f0aaad 15 0x735a2d97 15 */
static UV _hash(UV x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}
#else
/* 30 0xbf58476d1ce4e5b9 27 0x94d049bb133111eb 31     SplitMix64/Meuller */
/* 32 0xd6e8feb86659fd93 32 0xd6e8feb86659fd93 32     degski */
/* 33 0xff51afd7ed558ccd 33 0xc4ceb9fe1a85ec53 33     Murmur64 */
static UV _hash(UV x) {
    x = (x ^ (x >> 30)) * UVCONST(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UVCONST(0x94d049bb133111eb);
    x = x ^ (x >> 31);
    return x;
}
#endif

#define HVAL(x,mask) (_hash(x) & mask)

/******************************************************************************/

iset_t iset_create(unsigned long init_size) {
  iset_t set;
  int bits = 4;

  set.size = set.contains_zero = set.type = 0;

  while (bits < BITS_PER_WORD-1 && ((UVCONST(1) << bits) * FILL_RATIO + 1) < init_size)
    bits++;
  set.maxsize = UVCONST(1) << bits;
  set.mask = set.maxsize - 1;
  Newz(0, set.arr, set.maxsize, UV);
  return set;
}

void iset_destroy(iset_t *set) {
  set->maxsize = set->size = set->contains_zero = set->type = 0;
  Safefree(set->arr);
  set->arr = 0;
}

static unsigned long _iset_pos(const UV* arr, UV mask, UV val) {
  UV h = HVAL(val,mask);
  while (arr[h] != 0 && arr[h] != val)
    h = (h+1) & mask;
  return h;
}

int iset_contains(const iset_t set, UV val) {
  if (val == 0) return set.contains_zero;
  return set.arr[_iset_pos(set.arr, set.mask, val)] == val;
}

static void _iset_resize(iset_t *set) {
  UV v, newmask, *narr;
  unsigned long i, oldsize, newsize;

  oldsize = set->maxsize;
  newsize = oldsize << 1;
  if (newsize < oldsize) croak("iset: max set size overflow");
  newmask = newsize - 1;

  Newz(0, narr, newsize, UV);
  for (i = 0; i < oldsize; i++)
    if (v = set->arr[i], v != 0)
      narr[ _iset_pos(narr,newmask,v) ] = v;
  Safefree(set->arr);
  set->arr = narr;
  set->maxsize = newsize;
  set->mask = newmask;
}

int iset_add(iset_t *set, UV val, int sign) {
  if (sign == 0)
    set->type = ISET_TYPE_INVALID;
  else if (val > (UV)IV_MAX)
    set->type |= ((sign > 0) ? ISET_TYPE_UV : ISET_TYPE_IV);
  if (val == 0) {
    if (set->contains_zero)
      return 0;
    set->contains_zero = 1;
    set->size++;
  } else {
    UV h = _iset_pos(set->arr, set->mask, val);
    if (set->arr[h] == val)
      return 0;
    set->arr[h] = val;
    if (++set->size > FILL_RATIO * (double)set->maxsize)
      _iset_resize(set);
  }
  return 1;
}

iset_t iset_create_from_array(UV* d, unsigned long dlen, int dsign) {
  unsigned long i;
  unsigned char typemask = ((dsign > 0) ? ISET_TYPE_UV : ISET_TYPE_IV);
  iset_t s = iset_create(dlen);

  if (dsign != 0) {
    for (i = 0; i < dlen; i++) {
      UV h, val = d[i];
      if (val == 0) {
        if (!s.contains_zero) { s.contains_zero = 1;  s.size++; }
      } else {
        if (val > (UV)IV_MAX)
          s.type |= typemask;
        h = _iset_pos(s.arr, s.mask, val);
        if (s.arr[h] != val) {
          s.arr[h] = val;
          s.size++;
        }
      }
    }
  }
  return s;
}

void iset_allvals(const iset_t set, UV* array) {
  unsigned long j, i = 0;
  if (set.contains_zero)
    array[i++] = 0;
  for (j = 0; j < set.maxsize; j++)
    if (set.arr[j] != 0)
      array[i++] = set.arr[j];
  if (i != set.size) croak("iset_allvals bad size");
  if (set.type == ISET_TYPE_IV) sort_iv_array((IV*)array, i);
  else                          sort_uv_array(array, i);
}

#if 0
void iset_minmax(const iset_t set, UV *min, UV *max) {
  unsigned long i;
  UV v;

  if (set.type == ISET_TYPE_INVALID || set.size == 0) { *min=*max=0; return; }

  if (set.type != ISET_TYPE_IV) {
    *min = UV_MAX;
    *max = 0;
    for (i = 0; i < set.maxsize; i++)
      if (v = A.arr[i], v != 0) {
        if (v < *min) *min = v;
        if (v > *max) *max = v;
      }
    }
  } else {
    *min = IV_MAX;
    *max = IV_MIN;
    for (i = 0; i < set.maxsize; i++)
      if (v = A.arr[i], v != 0) {
        if ((IV)v < *min) *min = v;
        if ((IV)v > *max) *max = v;
      }
    }
  }
}
#endif


/******************************************************************************/

void iset_union_with(iset_t *set, const iset_t L) {
  unsigned long i, lsize;
  UV v, *larr;
  int lsign = iset_sign(L);;

  lsize = L.maxsize;
  larr = L.arr;
  for (i = 0; i < lsize; i++)
    if (v = larr[i], v != 0)
      iset_add(set, v, lsign);
  if (L.contains_zero && !set->contains_zero)  iset_add(set,0,1);
}

void iset_intersect_with(iset_t *set, const iset_t L) {
  iset_t s = iset_intersection_of(*set, L);
  iset_destroy(set);
  *set = s;
}

void iset_difference_with(iset_t *set, const iset_t L) {
  iset_t s = iset_difference_of(*set, L);
  iset_destroy(set);
  *set = s;
}

void iset_symdiff_with(iset_t *set, const iset_t L) {
  iset_t s = iset_symdiff_of(*set, L);
  iset_destroy(set);
  *set = s;
}

/******************************************************************************/

iset_t iset_union_of(const iset_t A, const iset_t B) {
  unsigned long i;
  UV v;
  int asign = iset_sign(A), bsign = iset_sign(B);
  iset_t s = iset_create(A.size + B.size);

  for (i = 0; i < A.maxsize; i++)
    if (v = A.arr[i], v != 0)
      iset_add(&s, v, asign);
  for (i = 0; i < B.maxsize; i++)
    if (v = B.arr[i], v != 0)
      iset_add(&s, v, bsign);
  if (A.contains_zero || B.contains_zero)  iset_add(&s,0,1);
  return s;
}

iset_t iset_intersection_of(const iset_t A, const iset_t B) {
  int asign = iset_sign(A), bsign = iset_sign(B);
  int samesign = (asign == bsign);
  unsigned long i;
  UV v;
  iset_t s;

  if (A.size > B.maxsize)               /* Swap for performance. */
    return iset_intersection_of(B,A);

  s = iset_create((A.size > B.size) ? A.size : B.size);

  for (i = 0; i < A.maxsize; i++)
    if (v = A.arr[i], v != 0)
      if ( !((v > (UV)IV_MAX) && !samesign) && iset_contains(B, v))
        iset_add(&s, v, asign);
  if (A.contains_zero && B.contains_zero)  iset_add(&s,0,1);
  return s;
}
iset_t iset_difference_of(const iset_t A, const iset_t B) {
  int asign = iset_sign(A), bsign = iset_sign(B);
  int samesign = (asign == bsign);
  unsigned long i;
  UV v;
  iset_t s = iset_create((A.size > B.size) ? A.size : B.size);

  for (i = 0; i < A.maxsize; i++)
    if (v = A.arr[i], v != 0)
      if ( ((v > (UV)IV_MAX) && !samesign) || !iset_contains(B, v) )
        iset_add(&s, v, asign);
  if (A.contains_zero && !B.contains_zero)  iset_add(&s,0,1);
  return s;
}
iset_t iset_symdiff_of(const iset_t A, const iset_t B) {
  int asign = iset_sign(A), bsign = iset_sign(B);
  int samesign = (asign == bsign);
  unsigned long i;
  UV v;
  iset_t s = iset_create((A.size > B.size) ? A.size : B.size);

  for (i = 0; i < A.maxsize; i++)
    if (v = A.arr[i], v != 0)
      if ( ((v > (UV)IV_MAX) && !samesign) || !iset_contains(B, v) )
        iset_add(&s, v, asign);
  for (i = 0; i < B.maxsize; i++)
    if (v = B.arr[i], v != 0)
      if ( ((v > (UV)IV_MAX) && !samesign) || !iset_contains(A, v) )
        iset_add(&s, v, bsign);
  if ((A.contains_zero + B.contains_zero) == 1)  iset_add(&s,0,1);
  return s;
}
int iset_is_subset_of(const iset_t A, const iset_t B) {
  int samesign = (iset_sign(A) == iset_sign(B));
  unsigned long i;
  UV v;

  if (A.size > B.size)
    return 0;
  if (A.contains_zero && !B.contains_zero)
    return 0;
  for (i = 0; i < A.maxsize; i++)
    if (v = A.arr[i], v != 0)
      if ( ((v > (UV)IV_MAX) && !samesign) || !iset_contains(B, v) )
        return 0;
  return 1;
}

/******************************************************************************/

void iset_test(void) {
  iset_t s;
  UV i, *S;
  const unsigned long ts = 30000000;

  printf("create .. "); fflush(stdout);
  s = iset_create(0);
  printf("done\n"); fflush(stdout);
  for (i = ts/2; i < ts; i++) {
    iset_add(&s, i, 1);
  }
  printf("done adding.  size is %lu\n", iset_size(s)); fflush(stdout);
  if (iset_contains(s,0) != 0) croak("fail 0");
  for (i = 0; i < ts; i++) {
    iset_add(&s, i, 1);
  }
  printf("done adding.  size is %lu\n", iset_size(s)); fflush(stdout);

  if (iset_contains(s,1) != 1) croak("fail 1");
  if (iset_contains(s,ts-1) != 1) croak("fail 999");
  if (iset_contains(s,ts) != 0) croak("fail 1000");
  if (iset_contains(s,0) != 1) croak("fail 0");
  if (iset_sign(s) != 1) croak("fail sign");
  if (iset_is_invalid(s) != 0) croak("fail invalid");
  if (iset_size(s) != ts) croak("fail size");

  New(0,S,iset_size(s),UV);
  iset_allvals(s,S);
  for (i = 0; i < ts; i++)
    if (S[i] != i)
      croak("fail element %lu expected %lu got %lu\n", i, i, S[i]);
  iset_destroy(&s);
}
