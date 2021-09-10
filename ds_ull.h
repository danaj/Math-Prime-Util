#ifndef MPU_DS_ULL_H
#define MPU_DS_ULL_H

#include "ptypes.h"

/******************************************************************************/
/*                       UNROLLED LINKED LIST DATA STRUCTURE                  */
/******************************************************************************/

/*
 * I made this data type for lucky numbers, and called it 'pagelist'.  After
 * some searching I found it is an 'unrolled linked list' published in 1994.
 * There are other names for similar structures, * e.g. Space Efficient List.
 * This is a singly linked list.
 *
 * ull32   for uint32_t
 * ull     for UV (or your type)
 */

/* By default we're storing UVs, but you can use something else */
#ifndef ULLTYPE
#define ULLTYPE UV
#endif
#ifndef ULLPAGESIZE
#define ULLPAGESIZE (2303*sizeof(ULLTYPE))
#endif
#ifndef ULL32PAGESIZE
#define ULL32PAGESIZE (1303*sizeof(uint32_t))
#endif

/******************************************************************************/

typedef struct _page32_t {
  struct _page32_t *next;
  uint32_t size;
  uint32_t data[ULL32PAGESIZE];
} _page32_t;

static _page32_t* _page32_create(void) {
  _page32_t *p;
  Newz(0, p, 1, _page32_t);
  return p;
}

typedef struct ull32_t {
  /* skip is the last-used page, fast is a settable index */
  uint32_t nelems, skip_amount, fast_amount;
  _page32_t *head, *tail, *skip_page, *fast_page;
} ull32_t;

static _page32_t* _page32_for_index(ull32_t *pl, uint32_t *idx) {
  uint32_t oidx = *idx;
  _page32_t *p;

  if (*idx >= pl->nelems) croak("ull32 index out of range");

  if (oidx >= pl->fast_amount) {
    if (oidx >= pl->skip_amount && pl->skip_amount >= pl->fast_amount)
                                 { p = pl->skip_page; *idx -= pl->skip_amount; }
    else                         { p = pl->fast_page; *idx -= pl->fast_amount; }
  } else {
    if (oidx >= pl->skip_amount) { p = pl->skip_page; *idx -= pl->skip_amount; }
    else                         { p = pl->head; }
  }

  for ( ;  p != 0;  p = p->next) {  /* Walk to correct page */
    if (*idx < p->size)  break;
    *idx -= p->size;
  }
  if (p == 0) croak("walked off list");
  pl->skip_amount = oidx - *idx;
  pl->skip_page = p;
  return p;
}

static void ull32_set_fastpage(ull32_t *pl, uint32_t idx) {
  uint32_t oidx = idx;
  _page32_t *p;
  if (oidx >= pl->nelems) croak("ull32 index out of range");
  if (oidx >= pl->fast_amount) {
    if (oidx >= pl->skip_amount && pl->skip_amount >= pl->fast_amount)
                                 { p = pl->skip_page; idx -= pl->skip_amount; }
    else                         { p = pl->fast_page; idx -= pl->fast_amount; }
  } else {
    if (oidx >= pl->skip_amount) { p = pl->skip_page; idx -= pl->skip_amount; }
    else                         { p = pl->head; }
  }
  for ( ;  p != 0;  p = p->next) {  /* Walk to correct page */
    if (idx < p->size)  break;
    idx -= p->size;
  }
  pl->fast_amount = oidx - idx;
  pl->fast_page = p;
}

static ull32_t* ull32_create(void) {
  ull32_t *pl;
  New(0, pl, 1, ull32_t);
  pl->head = pl->tail = _page32_create();
  pl->nelems = 0;
  pl->skip_amount = pl->fast_amount = 0;
  pl->skip_page = pl->fast_page = pl->head;
  return pl;
}

static void ull32_destroy(ull32_t *pl) {
  _page32_t *head = pl->head;
  while (head != 0) {
    _page32_t *next = head->next;
    Safefree(head);
    head = next;
  }
  pl->head = pl->tail = 0;
  pl->skip_amount = pl->fast_amount = 0;
  pl->skip_page = pl->fast_page = 0;
}

/* For insert function, do it by splitting this page into two if needed */

static void ull32_append(ull32_t *pl, uint32_t n) {
  _page32_t *p = pl->tail;
  pl->nelems++;
  if (p->size < ULL32PAGESIZE) {
    p->data[p->size++] = n;
  } else {
    _page32_t *q = _page32_create();
    q->size = 1;
    q->data[0] = n;
    p->next = q;
    pl->tail = q;
  }
}

static void ull32_delete(ull32_t *pl, uint32_t idx) {   /* Index 0,1,... */
  _page32_t *q, *p;
  /* Invalidate cached pages if needed. */
  if (pl->skip_amount > idx) { pl->skip_amount = 0; pl->skip_page = pl->head; }
  if (pl->fast_amount > idx) { pl->fast_amount = 0; pl->fast_page = pl->head; }
  /* Delete entry */
  p = _page32_for_index(pl, &idx);
  if (idx < p->size - 1)
    memmove(p->data+idx, p->data+idx+1, (p->size-idx-1) * sizeof(uint32_t));
  p->size--;
  pl->nelems--;
  /* Try to merge with next page */
  q = p->next;
  if (q != 0 && (p->size + q->size <= ULL32PAGESIZE)) {
    memcpy( p->data + p->size, q->data, q->size * sizeof(uint32_t) );
    p->size += q->size;
    p->next = q->next;
    Safefree(q);
  }
}

static uint32_t ull32_val(ull32_t *pl, uint32_t idx) {
  _page32_t *p = _page32_for_index(pl, &idx);
  return p->data[idx];
}

static uint32_t ull32_npages(ull32_t *pl) {
  _page32_t *p;
  uint32_t npages;
  for (npages = 0, p = pl->head;  p != 0;  p = p->next)
    npages++;
  return npages;
}

static uint32_t* ull32_toarray(UV *size, ull32_t *pl) {
  uint32_t i, *arr;
  _page32_t *p;
  New(0, arr, pl->nelems, uint32_t);
  for (i = 0, p = pl->head;  p != 0;  p = p->next) {
    if (pl->nelems - i < p->size) croak("ull32: bad number of elements in page");
    memcpy(arr + i,  p->data,  p->size * sizeof(uint32_t));
    i += p->size;
  }
  if (i != pl->nelems) croak("ull31: bad number of elements in list");
  *size = i;
  return arr;
}

typedef struct ull32_iter_t {
  _page32_t *p;
  uint32_t idx;
} ull32_iter_t;

static ull32_iter_t ull32_iterator_create(ull32_t *pl, uint32_t idx) {
  ull32_iter_t iter;
  iter.p = _page32_for_index(pl, &idx);
  iter.idx = idx;
  return iter;
}

static uint32_t ull32_iterator_next(ull32_iter_t *iter) {
  uint32_t v;
  _page32_t *P = iter->p;
  if (P == 0 || iter->idx >= P->size) return 0;
  v = P->data[iter->idx++];
  if (iter->idx >= P->size) {
    iter->idx = 0;
    iter->p = P->next;
  }
  return v;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

typedef struct _page_t {
  struct _page_t *next;
  UV size;
  ULLTYPE data[ULLPAGESIZE];
} _page_t;

static _page_t* _page_create(void) {
  _page_t *p;
  Newz(0, p, 1, _page_t);
  return p;
}

typedef struct ull_t {
  /* skip is the last-used page, fast is a settable index */
  UV nelems, skip_amount, fast_amount;
  _page_t *head, *tail, *skip_page, *fast_page;
} ull_t;

static _page_t* _page_for_index(ull_t *pl, UV *idx) {
  UV oidx = *idx;
  _page_t *p;

  if (*idx >= pl->nelems) croak("ull index out of range");

  if (oidx >= pl->fast_amount) {
    if (oidx >= pl->skip_amount && pl->skip_amount >= pl->fast_amount)
                                 { p = pl->skip_page; *idx -= pl->skip_amount; }
    else                         { p = pl->fast_page; *idx -= pl->fast_amount; }
  } else {
    if (oidx >= pl->skip_amount) { p = pl->skip_page; *idx -= pl->skip_amount; }
    else                         { p = pl->head; }
  }
  for ( ;  p != 0;  p = p->next) {  /* Walk to correct page */
    if (*idx < p->size)  break;
    *idx -= p->size;
  }
  if (p == 0) croak("walked off list");
  pl->skip_amount = oidx - *idx;
  pl->skip_page = p;
  return p;
}

static void ull_set_fastpage(ull_t *pl, UV idx) {
  UV oidx = idx;
  _page_t *p;
  if (oidx >= pl->nelems) croak("ull32 index out of range");
  if (oidx >= pl->fast_amount) {
    if (oidx >= pl->skip_amount && pl->skip_amount >= pl->fast_amount)
                                 { p = pl->skip_page; idx -= pl->skip_amount; }
    else                         { p = pl->fast_page; idx -= pl->fast_amount; }
  } else {
    if (oidx >= pl->skip_amount) { p = pl->skip_page; idx -= pl->skip_amount; }
    else                         { p = pl->head; }
  }
  for ( ;  p != 0;  p = p->next) {  /* Walk to correct page */
    if (idx < p->size)  break;
    idx -= p->size;
  }
  pl->fast_amount = oidx - idx;
  pl->fast_page = p;
}

static ull_t* ull_create(void) {
  ull_t *pl;
  New(0, pl, 1, ull_t);
  pl->head = pl->tail = _page_create();
  pl->nelems = 0;
  pl->skip_amount = pl->fast_amount = 0;
  pl->skip_page = pl->fast_page = pl->head;
  return pl;
}

static void ull_destroy(ull_t *pl) {
  _page_t *head = pl->head;
  while (head != 0) {
    _page_t *next = head->next;
    Safefree(head);
    head = next;
  }
  pl->head = pl->tail = 0;
  pl->skip_amount = pl->fast_amount = 0;
  pl->skip_page = pl->fast_page = 0;
}

/* For insert function, do it by splitting this page into two if needed */

static void ull_append(ull_t *pl, ULLTYPE n) {
  _page_t *p = pl->tail;
  pl->nelems++;
  if (p->size < ULLPAGESIZE) {
    p->data[p->size++] = n;
  } else {
    _page_t *q = _page_create();
    q->size = 1;
    q->data[0] = n;
    p->next = q;
    pl->tail = q;
  }
}

static void ull_delete(ull_t *pl, UV idx) {   /* Index 0,1,... */
  _page_t *q, *p;
  /* Invalidate cached pages if needed. */
  if (pl->skip_amount > idx) { pl->skip_amount = 0; pl->skip_page = pl->head; }
  if (pl->fast_amount > idx) { pl->fast_amount = 0; pl->fast_page = pl->head; }
  /* Delete entry */
  p = _page_for_index(pl, &idx);
  if (idx < p->size - 1)
    memmove(p->data+idx, p->data+idx+1, (p->size-idx-1) * sizeof(ULLTYPE));
  p->size--;
  pl->nelems--;
  /* Try to merge with next page */
  q = p->next;
  if (q != 0 && (p->size + q->size <= ULLPAGESIZE)) {
    memcpy( p->data + p->size, q->data, q->size * sizeof(ULLTYPE) );
    p->size += q->size;
    p->next = q->next;
    Safefree(q);
  }
}

static ULLTYPE ull_val(ull_t *pl, UV idx) {
  _page_t *p = _page_for_index(pl, &idx);
  return p->data[idx];
}

static UV ull_npages(ull_t *pl) {
  _page_t *p;
  UV npages;
  for (npages = 0, p = pl->head;  p != 0;  p = p->next)
    npages++;
  return npages;
}

static ULLTYPE* ull_toarray(UV *size, ull_t *pl) {
  ULLTYPE *arr;
  UV i;
  _page_t *p;
  New(0, arr, pl->nelems, ULLTYPE);
  for (i = 0, p = pl->head;  p != 0;  p = p->next) {
    if (pl->nelems - i < p->size) croak("ull: bad number of elements in page");
    memcpy(arr + i,  p->data,  p->size * sizeof(ULLTYPE));
    i += p->size;
  }
  if (i != pl->nelems) croak("ull: bad number of elements in list");
  *size = i;
  return arr;
}

typedef struct ull_iter_t {
  _page_t *p;
  UV idx;
} ull_iter_t;

static ull_iter_t ull_iterator_create(ull_t *pl, UV idx) {
  ull_iter_t iter;
  iter.p = _page_for_index(pl, &idx);
  iter.idx = idx;
  return iter;
}

static ULLTYPE ull_iterator_next(ull_iter_t *iter) {
  ULLTYPE v;
  _page_t *P = iter->p;
  if (P == 0 || iter->idx >= P->size) return 0;
  v = P->data[iter->idx++];
  if (iter->idx >= P->size) {
    iter->idx = 0;
    iter->p = P->next;
  }
  return v;
}

/******************************************************************************/
/******************************************************************************/

#endif
