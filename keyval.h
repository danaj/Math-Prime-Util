#ifndef MPU_KEYVAL_H
#define MPU_KEYVAL_H

#include "ptypes.h"

typedef struct {
  UV key;
  UV val;
} keyval_t;

typedef struct {
  keyval_t *keyval;
  UV mask;
  long maxsize;
  long size;
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


/******************************************************************************/

static void init_set(set_t *S, UV isize) {
  int bits = 0;
  while (isize > 0) {
    bits++;
    isize >>= 1;
  }
  S->size = 0;
  S->maxsize = UVCONST(1) << ((bits < 3) ? 3 : bits);
  S->mask = S->maxsize - 1;
  Newz(0,S->keyval,S->maxsize,keyval_t);
}

static void free_set(set_t *S) {
  S->size = S->maxsize = 0;
  Safefree(S->keyval);
}

static void _set_expand(set_t *S) {
  long i, max = S->maxsize, newmax = max*2, newsize = 0, newmask = newmax-1;
  keyval_t *nkv;
  Newz(0, nkv, newmax, keyval_t);
  for (i = 0; i < max; i++) {
    UV key = S->keyval[i].key;
    if (key != 0) {
      UV h = _hash(key) & newmask;
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

static long set_search(set_t S, UV key) {
  long h = _hash(key) & S.mask;
  while (S.keyval[h].key > 0 && S.keyval[h].key != key)
    h = (h+1) & S.mask;   /* Linear probe */
  return (S.keyval[h].key == key) ? h : -1;
}

static UV set_getval(set_t S, UV key) {
  long i = set_search(S, key);
  return (i == -1) ? 0 : S.keyval[i].val;
}

static void set_addsum(set_t *S, keyval_t kv) {
  UV h = _hash(kv.key) & S->mask;
  while (S->keyval[h].key > 0 && S->keyval[h].key != kv.key)
    h = (h+1) & S->mask;
  if (S->keyval[h].key == kv.key) {
    /* if (kv.val > UV_MAX - S->keyval[h].val) croak("add overflow\n"); */
    S->keyval[h].val += kv.val;
  } else {
    S->keyval[h] = kv;
    if (S->size++ > 0.65 * S->maxsize)
      _set_expand(S);
  }
}

static void set_merge(set_t *S, set_t T) {
  long j;
  for (j = 0; j < T.maxsize; j++)
    if (T.keyval[j].key > 0)
      set_addsum(S, T.keyval[j]);
}

/******************************************************************************/

typedef struct {
  UV key;
  UV *vals;
  long size;
  long maxsize;
} keylist_t;

typedef struct {
  keylist_t *keylist;
  UV mask;
  long maxsize;
  long size;
} set_list_t;

static void init_setlist(set_list_t *L, UV isize) {
  int bits = 0;
  while (isize > 0) {
    bits++;
    isize >>= 1;
  }
  L->size = 0;
  L->maxsize = UVCONST(1) << ((bits < 3) ? 3 : bits);
  L->mask = L->maxsize - 1;
  Newz(0, L->keylist, L->maxsize, keylist_t);
}

static void free_setlist(set_list_t *L) {
  long i;
  for (i = 0; i < L->maxsize; i++)
    if (L->keylist[i].size > 0)
      Safefree(L->keylist[i].vals);
  Safefree(L->keylist);
  L->size = L->maxsize = 0;
}

static void _setlist_expand(set_list_t *L) {
  long i, max = L->maxsize, newmax = max*2, newsize = 0, newmask = newmax-1;
  keylist_t *nlist;
  Newz(0, nlist, newmax, keylist_t);
  for (i = 0; i < max; i++) {
    UV key = L->keylist[i].key;
    if (key != 0) {
      UV h = _hash(key) & newmask;
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

static long setlist_search(set_list_t L, UV key) {
  long h = _hash(key) & L.mask;
  while (L.keylist[h].key > 0 && L.keylist[h].key != key)
    h = (h+1) & L.mask;   /* Linear probe */
  return (L.keylist[h].key == key) ? h : -1;
}

static void setlist_addlist(set_list_t *L, UV key, long nvals, UV* list, UV mult) {
  UV *vptr;
  long j, h = _hash(key) & L->mask;
  while (L->keylist[h].key > 0 && L->keylist[h].key != key)
    h = (h+1) & L->mask;
  if (L->keylist[h].key == key) {
    long size = L->keylist[h].size;
    long maxsize = L->keylist[h].maxsize;
    if (size + nvals > maxsize) {
      maxsize = 2 * (size+nvals);
      Renew(L->keylist[h].vals, maxsize, UV);
      L->keylist[h].maxsize = maxsize;
    }
    vptr = L->keylist[h].vals + size;
    for (j = 0; j < nvals; j++)
      vptr[j] = list[j] * mult;
    L->keylist[h].size = size + nvals;
  } else {
    long maxsize = (nvals < 5) ? 12 : (nvals+1) * 2;
    New(0, L->keylist[h].vals, maxsize, UV);
    L->keylist[h].maxsize = maxsize;
    vptr = L->keylist[h].vals;
    for (j = 0; j < nvals; j++)
      vptr[j] = list[j] * mult;
    L->keylist[h].size = nvals;
    L->keylist[h].key = key;
    if (L->size++ > 0.65 * L->maxsize)
      _setlist_expand(L);
  }
}

static void setlist_addval(set_list_t *L, UV key, UV val) {
  setlist_addlist(L, key, 1, &val, 1);
}

static UV* setlist_getlist(UV *nvals, set_list_t L, UV key) {
  long i = setlist_search(L, key);
  if (i == -1) {
    *nvals = 0;
    return 0;
  }
  *nvals = L.keylist[i].size;
  return L.keylist[i].vals;
}

static void setlist_merge(set_list_t *L, set_list_t T) {
  long j;
  for (j = 0; j < T.maxsize; j++) {
    if (T.keylist[j].key > 0) {
      UV key   = T.keylist[j].key;
      UV nvals = T.keylist[j].size;
      UV *vals = T.keylist[j].vals;
      setlist_addlist(L, key, nvals, vals, 1);
    }
  }
}

#endif
