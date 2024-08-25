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

  set.size = set.contains_zero = 0;

  while (bits < BITS_PER_WORD-1 && ((UVCONST(1) << bits) * FILL_RATIO + 1) < init_size)
    bits++;
  set.maxsize = UVCONST(1) << bits;
  set.mask = set.maxsize - 1;
  Newz(0, set.arr, set.maxsize, UV);
  return set;
}

void iset_destroy(iset_t *set) {
  set->maxsize = set->size = set->contains_zero = 0;
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

int iset_add(iset_t *set, UV val) {
  UV h;
  if (val == 0) {
    if (set->contains_zero)  return 0;
    set->contains_zero = 1;
    set->size++;
    return 1;
  }
  h = _iset_pos(set->arr, set->mask, val);
  if (set->arr[h] == val) return 0;
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
  sort_uv_array(array, i);
}

void iset_test(void) {
  iset_t s;
  UV i, sz, *S;
  const unsigned long ts = 30000000;

  printf("create .. "); fflush(stdout);
  s = iset_create(0);
  printf("done\n"); fflush(stdout);
  for (i = ts/2; i < ts; i++) {
    iset_add(&s, i);
  }
  printf("done adding.  size is %lu\n", iset_size(s)); fflush(stdout);
  if (iset_contains(s,0) != 0) croak("fail 0");
  for (i = 0; i < ts; i++) {
    iset_add(&s, i);
  }
  printf("done adding.  size is %lu\n", iset_size(s)); fflush(stdout);

  if (iset_contains(s,1) != 1) croak("fail 1");
  if (iset_contains(s,ts-1) != 1) croak("fail 999");
  if (iset_contains(s,ts) != 0) croak("fail 1000");
  if (iset_contains(s,0) != 1) croak("fail 0");
  if (iset_size(s) != ts) croak("fail size");

  New(0,S,iset_size(s),UV);
  iset_allvals(s,S);
  for (i = 0; i < ts; i++)
    if (S[i] != i)
      croak("fail element %lu expected %lu got %lu\n", i, i, S[i]);
  iset_destroy(&s);
}
