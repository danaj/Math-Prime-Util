#ifndef MPU_KEYVAL_H
#define MPU_KEYVAL_H

/* This includes:
 *   keyval_t   simple key/val type, both UV
 *   set        a key value set, with "add" function for new=old+new
 *   setlist    key (UV) plus dynamic array of UVs.  "append" functionality
 *
 * Key=0 is not allowed.
 */


#include "ptypes.h"

typedef struct {
  UV key;
  UV val;
} keyval_t;

typedef struct {
  keyval_t *keyval;
  Size_t mask;
  Size_t maxsize;
  Size_t size;
} set_t;




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


static Size_t _keyval_grow_size(Size_t size, Size_t elemsize) {
  Size_t maxelems = MAX_SIZET / elemsize;
  if (size > maxelems/2)
    croak("keyval allocation is too large");
  return 2 * size;
}

static Size_t _keyval_initial_size(Size_t isize, Size_t elemsize) {
  Size_t maxsize = 8;
  if (maxsize > MAX_SIZET/elemsize)
    croak("keyval allocation is too large");
  while (maxsize <= isize)
    maxsize = _keyval_grow_size(maxsize, elemsize);
  return maxsize;
}

static bool _keyval_too_full(Size_t size, Size_t maxsize) {
  Size_t limit = (maxsize/20) * 13 + ((maxsize%20) * 13) / 20;
  return size > limit;
}


/******************************************************************************/

static void init_set(set_t *S, Size_t isize) {
  S->size = 0;
  S->maxsize = _keyval_initial_size(isize, sizeof(keyval_t));
  S->mask = S->maxsize - 1;
  Newz(0,S->keyval,S->maxsize,keyval_t);
}

static void free_set(set_t *S) {
  Safefree(S->keyval);
  S->keyval = 0;
  S->mask = S->size = S->maxsize = 0;
}

static void _set_expand(set_t *S) {
  Size_t i, max = S->maxsize, newsize = 0;
  Size_t newmax = _keyval_grow_size(max, sizeof(keyval_t));
  Size_t newmask = newmax-1;
  keyval_t *nkv;
  Newz(0, nkv, newmax, keyval_t);
  for (i = 0; i < max; i++) {
    UV key = S->keyval[i].key;
    if (key != 0) {
      Size_t h = (Size_t)_hash(key) & newmask;
      while (nkv[h].key > 0 && nkv[h].key != key)
        h = (h+1) & newmask;
      nkv[h] = S->keyval[i];
      newsize++;
    }
  }
  Safefree(S->keyval);
  S->keyval = nkv;
  S->maxsize = newmax;
  S->mask = newmax-1;
  MPUassert(newsize == S->size, "keyval set size mismatch");
}

static SSize_t set_search(set_t S, UV key) {
  Size_t h;
  MPUassert(key != 0, "keyval key must not be zero");
  h = (Size_t)_hash(key) & S.mask;
  while (S.keyval[h].key > 0 && S.keyval[h].key != key)
    h = (h+1) & S.mask;   /* Linear probe */
  return (S.keyval[h].key == key) ? (SSize_t)h : -1;
}

static UV set_getval(set_t S, UV key) {
  SSize_t i = set_search(S, key);
  return (i == -1) ? 0 : S.keyval[i].val;
}

static void set_addsum(set_t *S, keyval_t kv) {
  Size_t h;
  MPUassert(kv.key != 0, "keyval key must not be zero");
  h = (Size_t)_hash(kv.key) & S->mask;
  while (S->keyval[h].key > 0 && S->keyval[h].key != kv.key)
    h = (h+1) & S->mask;
  if (S->keyval[h].key == kv.key) {
    if (kv.val > UV_MAX - S->keyval[h].val)
      croak("keyval set value overflow");
    S->keyval[h].val += kv.val;
  } else {
    S->keyval[h] = kv;
    S->size++;
    if (_keyval_too_full(S->size, S->maxsize))
      _set_expand(S);
  }
}

static void set_merge(set_t *S, set_t T) {
  Size_t j;
  for (j = 0; j < T.maxsize; j++)
    if (T.keyval[j].key > 0)
      set_addsum(S, T.keyval[j]);
}

/******************************************************************************/

typedef struct {
  UV key;
  UV *vals;
  Size_t size;
  Size_t maxsize;
} keylist_t;

typedef struct {
  keylist_t *keylist;
  Size_t mask;
  Size_t maxsize;
  Size_t size;
} set_list_t;

static void init_setlist(set_list_t *L, Size_t isize) {
  L->size = 0;
  L->maxsize = _keyval_initial_size(isize, sizeof(keylist_t));
  L->mask = L->maxsize - 1;
  Newz(0, L->keylist, L->maxsize, keylist_t);
}

static void free_setlist(set_list_t *L) {
  Size_t i;
  for (i = 0; i < L->maxsize; i++)
    if (L->keylist[i].vals != 0)
      Safefree(L->keylist[i].vals);
  Safefree(L->keylist);
  L->keylist = 0;
  L->mask = L->size = L->maxsize = 0;
}

static void _setlist_expand(set_list_t *L) {
  Size_t i, max = L->maxsize, newsize = 0;
  Size_t newmax = _keyval_grow_size(max, sizeof(keylist_t));
  Size_t newmask = newmax-1;
  keylist_t *nlist;
  Newz(0, nlist, newmax, keylist_t);
  for (i = 0; i < max; i++) {
    UV key = L->keylist[i].key;
    if (key != 0) {
      Size_t h = (Size_t)_hash(key) & newmask;
      while (nlist[h].key > 0 && nlist[h].key != key)
        h = (h+1) & newmask;
      nlist[h] = L->keylist[i];
      newsize++;
    }
  }
  Safefree(L->keylist);
  L->keylist = nlist;
  L->maxsize = newmax;
  L->mask = newmax-1;
  MPUassert(newsize == L->size, "setlist size mismatch");
}

static SSize_t setlist_search(set_list_t L, UV key) {
  Size_t h;
  MPUassert(key != 0, "keyval key must not be zero");
  h = (Size_t)_hash(key) & L.mask;
  while (L.keylist[h].key > 0 && L.keylist[h].key != key)
    h = (h+1) & L.mask;   /* Linear probe */
  return (L.keylist[h].key == key) ? (SSize_t)h : -1;
}

static void setlist_addlist(set_list_t *L, UV key, Size_t nvals,
                            const UV *list, UV mult) {
  UV *vptr;
  UV maxval;
  Size_t j, h;

  MPUassert(key != 0, "keyval key must not be zero");
  if (nvals == 0) return;
  MPUassert(list != 0, "keyval list must not be null");
  maxval = (mult == 0) ? UV_MAX : UV_MAX/mult;
  h = (Size_t)_hash(key) & L->mask;
  while (L->keylist[h].key > 0 && L->keylist[h].key != key)
    h = (h+1) & L->mask;
  if (L->keylist[h].key == key) {
    Size_t size = L->keylist[h].size;
    Size_t maxsize = L->keylist[h].maxsize;
    Size_t newsize;
    bool selflist = (list == L->keylist[h].vals);
    if (nvals > MAX_SIZET-size)
      croak("keyval list size overflow");
    newsize = size + nvals;
    if (newsize > maxsize) {
      maxsize = _keyval_grow_size(newsize, sizeof(UV));
      Renew(L->keylist[h].vals, maxsize, UV);
      L->keylist[h].maxsize = maxsize;
      if (selflist) list = L->keylist[h].vals;
    }
    vptr = L->keylist[h].vals + size;
    for (j = 0; j < nvals; j++) {
      if (list[j] > maxval)
        croak("keyval list value overflow");
      vptr[j] = list[j] * mult;
    }
    L->keylist[h].size = newsize;
  } else {
    Size_t maxsize;
    if (nvals < 5) {
      maxsize = 12;
    } else {
      if (nvals == MAX_SIZET)
        croak("keyval list size overflow");
      maxsize = _keyval_grow_size(nvals+1, sizeof(UV));
    }
    New(0, L->keylist[h].vals, maxsize, UV);
    L->keylist[h].maxsize = maxsize;
    vptr = L->keylist[h].vals;
    for (j = 0; j < nvals; j++) {
      if (list[j] > maxval)
        croak("keyval list value overflow");
      vptr[j] = list[j] * mult;
    }
    L->keylist[h].size = nvals;
    L->keylist[h].key = key;
    L->size++;
    if (_keyval_too_full(L->size, L->maxsize))
      _setlist_expand(L);
  }
}

static void setlist_addval(set_list_t *L, UV key, UV val) {
  setlist_addlist(L, key, 1, &val, 1);
}

static UV* setlist_getlist(Size_t *nvals, set_list_t L, UV key) {
  SSize_t i = setlist_search(L, key);
  if (i == -1) {
    *nvals = 0;
    return 0;
  }
  *nvals = L.keylist[i].size;
  return L.keylist[i].vals;
}

static void setlist_merge(set_list_t *L, set_list_t T) {
  Size_t j;
  for (j = 0; j < T.maxsize; j++) {
    if (T.keylist[j].key > 0) {
      UV key   = T.keylist[j].key;
      Size_t nvals = T.keylist[j].size;
      const UV *vals = T.keylist[j].vals;
      setlist_addlist(L, key, nvals, vals, 1);
    }
  }
}

#if 0
static void setlist_zerolist(set_list_t *L, UV key) {
  SSize_t i = setlist_search(*L, key);
  if (i != -1) {
    Safefree(L->keylist[i].vals);
    L->keylist[i].vals = 0;
    L->keylist[i].size = L->keylist[i].maxsize = 0;
  }
}
#endif

#endif
