#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "ds_iset.h"
#include "util.h"

#define FILL_RATIO 0.50

#if BITS_PER_WORD == 32
static UV _hash(UV x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}
#else
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

  set.size = set.contains_zero = set.seen_uv = 0;
  set.sign = 1;

  while (bits < BITS_PER_WORD-1 && ((UVCONST(1) << bits) * FILL_RATIO + 1) < init_size)
    bits++;
  set.maxsize = UVCONST(1) << bits;
  set.mask = set.maxsize - 1;
  Newz(0, set.arr, set.maxsize, UV);
  return set;
}

void iset_destroy(iset_t *set) {
  set->maxsize = set->size = 0;
  set->contains_zero = set->seen_uv = 0;
  set->sign = 1;
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
  UV h;
  if (sign != 1 && sign != -1)
    set->sign = 0;
  if (val > (UV)IV_MAX) {
    if (sign == 1) set->seen_uv = 1;
    if (set->seen_uv && (sign != set->sign)) set->sign = 0;
    else if (sign == -1 && set->sign ==  1)  set->sign = -1;
  }
  if (val == 0) {
    if (set->contains_zero)  return 0;
    set->contains_zero = 1;
    set->size++;
    return 1;
  }
  h = _iset_pos(set->arr, set->mask, val);
  if (set->arr[h] == val)
    return 0;
  set->arr[h] = val;
  if (++set->size > FILL_RATIO * (double)set->maxsize)
    _iset_resize(set);
  return 1;
}

void iset_allvals(const iset_t set, UV* array) {
  unsigned long j, i = 0;
  if (set.contains_zero)
    array[i++] = 0;
  for (j = 0; j < set.maxsize; j++)
    if (set.arr[j] != 0)
      array[i++] = set.arr[j];
  if (i != set.size) croak("iset_allvals bad size");
  if (set.sign == -1) sort_iv_array((IV*)array, i);
  else                sort_uv_array(array, i);
}


/******************************************************************************/

void iset_union_with(iset_t *set, const iset_t L) {
  unsigned long i, lsize;
  UV v, *larr;
  int lsign;

  lsize = L.maxsize;
  larr = L.arr;
  lsign = L.sign;
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
  iset_t s = iset_create(A.size + B.size);

  for (i = 0; i < A.maxsize; i++)
    if (v = A.arr[i], v != 0)
      iset_add(&s, v, A.sign);
  for (i = 0; i < B.maxsize; i++)
    if (v = B.arr[i], v != 0)
      iset_add(&s, v, B.sign);
  if (A.contains_zero || B.contains_zero)  iset_add(&s,0,1);
  return s;
}

iset_t iset_intersection_of(const iset_t A, const iset_t B) {
  int samesign = (A.sign == B.sign);
  unsigned long i;
  UV v;
  iset_t s = iset_create((A.size > B.size) ? A.size : B.size);

  for (i = 0; i < A.maxsize; i++)
    if (v = A.arr[i], v != 0)
      if ( !((v > (UV)IV_MAX) && !samesign) && iset_contains(B, v))
        iset_add(&s, v, A.sign);
  if (A.contains_zero && B.contains_zero)  iset_add(&s,0,1);
  return s;
}
iset_t iset_difference_of(const iset_t A, const iset_t B) {
  int samesign = (A.sign == B.sign);
  unsigned long i;
  UV v;
  iset_t s = iset_create((A.size > B.size) ? A.size : B.size);

  for (i = 0; i < A.maxsize; i++)
    if (v = A.arr[i], v != 0)
      if ( ((v > (UV)IV_MAX) && !samesign) || !iset_contains(B, v) )
        iset_add(&s, v, A.sign);
  if (A.contains_zero && !B.contains_zero)  iset_add(&s,0,1);
  return s;
}
iset_t iset_symdiff_of(const iset_t A, const iset_t B) {
  int samesign = (A.sign == B.sign);
  unsigned long i;
  UV v;
  iset_t s = iset_create((A.size > B.size) ? A.size : B.size);

  for (i = 0; i < A.maxsize; i++)
    if (v = A.arr[i], v != 0)
      if ( ((v > (UV)IV_MAX) && !samesign) || !iset_contains(B, v) )
        iset_add(&s, v, A.sign);
  for (i = 0; i < B.maxsize; i++)
    if (v = B.arr[i], v != 0)
      if ( ((v > (UV)IV_MAX) && !samesign) || !iset_contains(A, v) )
        iset_add(&s, v, B.sign);
  if ((A.contains_zero + B.contains_zero) == 1)  iset_add(&s,0,1);
  return s;
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
  if (iset_size(s) != ts) croak("fail size");

  New(0,S,iset_size(s),UV);
  iset_allvals(s,S);
  for (i = 0; i < ts; i++)
    if (S[i] != i)
      croak("fail element %lu expected %lu got %lu\n", i, i, S[i]);
  iset_destroy(&s);
}
