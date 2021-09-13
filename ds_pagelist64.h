#ifndef MPU_DS_PAGELIST64_H
#define MPU_DS_PAGELIST64_H

#include "ptypes.h"

/******************************************************************************/
/*                           PAGELIST DATA STRUCTURE                          */
/******************************************************************************/

/*
 * This is a form of unrolled skip list.  However, for the sake of vastly
 * improved cache hits to find the right page, we use two arrays rather
 * than a linked list.
 */

#ifndef PLTYPE
#define PLTYPE UV
#endif

#define PLDEBUG    0
#define PAGESIZE   256
#define ELEV       7
#define ESHIFT     4
#define DELFACTOR  0.66  /* 1.0 means always merge */
#define NEPG(n,i)  (((n) + (1UL<<((i)*ESHIFT)) - 1) >> ((i)*ESHIFT))
#define ADDSIZE(pl, p, n) \
  { int _i; \
    for (_i = 0; _i < ELEV; _i++) \
      { pl->pgsize[_i][(p) >> (_i*ESHIFT)] += n; } }


typedef struct pagelist64_t {
  uint32_t pagesize;
  PLTYPE   nelems;
  PLTYPE   npages_allocated[ELEV];
  PLTYPE   npages[ELEV];
  PLTYPE   *pgsize[ELEV];
  PLTYPE   **pgdata;
} pagelist64_t;

static pagelist64_t* pagelist64_create(void) {
  pagelist64_t *pl;
  int i;
  New(0, pl, 1, pagelist64_t);
  pl->pagesize = PAGESIZE;
  pl->nelems = 0;

  for (i = 0; i < ELEV; i++) {
    pl->npages[i] = 1;
    pl->npages_allocated[i] = (i == 0)  ?  (1UL << ESHIFT)  :  1;
    Newz(0, pl->pgsize[i], pl->npages_allocated[i], PLTYPE);
  }
  New(0, pl->pgdata, pl->npages_allocated[0], PLTYPE *);
  New(0, pl->pgdata[0], pl->pagesize, PLTYPE);
  return pl;
}

static void pagelist64_destroy(pagelist64_t *pl) {
  PLTYPE p;
  int i;

  for (p = 0; p < pl->npages[0]; p++)
    Safefree(pl->pgdata[p]);
  Safefree(pl->pgdata);

  for (i = 0; i < ELEV; i++) {
    Safefree(pl->pgsize[i]);
    pl->npages[i] = pl->npages_allocated[i] = 0;
  }
  pl->nelems = 0;
}

#if PLDEBUG
#define CHECKPL(pl, msg)  _pagelist64_check_epg(pl, msg)
static void _pagelist64_check_epg(pagelist64_t *pl, const char* msg) {
  PLTYPE p, npages0 = pl->npages[0],  *pgsz0 = pl->pgsize[0],  sum[ELEV] = {0};
  int i;

  for (i = 0; i < ELEV; i++)
    if (pl->npages[i] > pl->npages_allocated[i])
      croak("level %u more pages in use than allocated\n", i);

  for (i = 1; i < ELEV; i++)
    if ( NEPG(npages0,i) > pl->npages[i] )
      croak("%s: level %u not enough pages nepg\n",msg,i);

  for (i = 1; i < ELEV; i++)
    if ( ((npages0-1) >> (i*ESHIFT)) >= pl->npages[i] )
      croak("%s: level %u not enough pages\n",msg,i);

  for (p = 0; p < npages0; p++) {
    for (i = 1; i < ELEV; i++) {
      PLTYPE pagesper = 1UL << (i*ESHIFT);
      sum[i] += pgsz0[p];
      if (p == npages0-1 || (p % pagesper) == (pagesper-1)) {
        PLTYPE ep = p >> (i*ESHIFT);
        if (sum[i] != pl->pgsize[i][ep])
          croak("%s: bad epg:  sum %u  pgsize[%u][%u] %u\n", msg, sum[i], i, ep, pl->pgsize[i][ep]);
        sum[i] = 0;
      }
    }
  }
}
#else
#define CHECKPL(pl, msg)
#endif

#if 0
static void _pagelist64_remake_epg(pagelist64_t *pl) {
  PLTYPE i, p,  npages0 = pl->npages[0],  *pgsz0 = pl->pgsize[0];

  for (i = 1; i < ELEV; i++) {
    PLTYPE* pgszi = pl->pgsize[i];
    memset( pgszi, 0, pl->npages[i] * sizeof(PLTYPE) );
    for (p = 0; p < npages0; p++)
      pgszi[ p >> (i*ESHIFT) ] += pgsz0[p];
  }
  CHECKPL(pl, "remake");
}
#endif

static PLTYPE _pagelist64_find_page(pagelist64_t *pl, PLTYPE *n) {
  PLTYPE p;
  int i;

  if (*n >= pl->nelems) croak("pagelist64 index out of range");

  CHECKPL(pl, "find page");
  for (i = ELEV-1, p = 0; i >= 0; i--) {
    PLTYPE npagesi = pl->npages[i],  *pgszi  = pl->pgsize[i];
    for (p = p << ESHIFT;  p < npagesi && *n >= pgszi[p];  p++)
      *n -= pgszi[p];
  }
  return p;
}

static void _pagelist64_add_page(pagelist64_t *pl) {
  PLTYPE npages0;
  int i;
  if (pl->npages[0] == pl->npages_allocated[0]) {
    pl->npages_allocated[0] += (1UL << ESHIFT);
    Renew(pl->pgsize[0], pl->npages_allocated[0], PLTYPE);
    Renew(pl->pgdata, pl->npages_allocated[0], PLTYPE*);
  }
  /* Mark one more page in use */
  npages0 = pl->npages[0];
  pl->pgsize[0][npages0] = 0;
  New(0, pl->pgdata[npages0], pl->pagesize, PLTYPE);
  pl->npages[0]++;

  /* Enlarge the higher levels if needed */
  for (i = 1; i < ELEV; i++) {
    if (NEPG(pl->npages_allocated[0], i) > pl->npages_allocated[i]) {
      pl->npages_allocated[i] += 4 + NEPG(pl->npages_allocated[0], i);
      Renew(pl->pgsize[i], pl->npages_allocated[i], PLTYPE);
    }
    if (NEPG(pl->npages[0], i) > pl->npages[i]) {
      pl->pgsize[i][pl->npages[i]] = 0;
      pl->npages[i]++;
    }
  }
  CHECKPL(pl, "add page");
}
static void pagelist64_append(pagelist64_t *pl, PLTYPE v) {
  PLTYPE *pgsz0  = pl->pgsize[0],  ptail = pl->npages[0] - 1;

  if (pgsz0[ptail] >= pl->pagesize) {
    _pagelist64_add_page(pl);
    pgsz0   = pl->pgsize[0];
    ptail   = pl->npages[0] - 1;
  }
  pl->pgdata[ptail][pgsz0[ptail]] = v;
  ADDSIZE(pl, ptail, 1);
  pl->nelems++;
  CHECKPL(pl, "append");
}
static PLTYPE pagelist64_val(pagelist64_t *pl, PLTYPE idx) {
  PLTYPE p = _pagelist64_find_page(pl, &idx);
#if PLDEBUG
  if (p >= pl->npages[0]) croak("pagelist64: bad page in val");
  if (idx >= pl->pgsize[0][p]) croak("pagelist64: bad index in val");
#endif
  return pl->pgdata[p][idx];
}
static PLTYPE* pagelist64_to_array(UV *size, pagelist64_t *pl) {
  PLTYPE *pgsz0  = pl->pgsize[0];
  PLTYPE d, p, *arr;
  New(0, arr, pl->nelems, PLTYPE);
  for (d = 0, p = 0; p < pl->npages[0]; p++) {
    memcpy(arr + d, pl->pgdata[p], pgsz0[p] * sizeof(PLTYPE));
    d += pgsz0[p];
  }
  if (d != pl->nelems) croak("pagelist64: bad number of elements in list");
  *size = d;
  return arr;
}

static void pagelist64_delete(pagelist64_t *pl, PLTYPE idx) { /* idx 0,1,... */
  PLTYPE npages0 = pl->npages[0],  *pgsz0 = pl->pgsize[0];
  PLTYPE p = _pagelist64_find_page(pl, &idx);
  if (p >= npages0) croak("pagelist64: bad page in delete");
  if (idx < pgsz0[p]-1)
    memmove(pl->pgdata[p]+idx, pl->pgdata[p]+idx+1, (pgsz0[p]-1-idx) * sizeof(PLTYPE));
  ADDSIZE(pl, p, -1);
  pl->nelems--;
#if 1  /* Merge with leveled add/subtract. */
  if ((p+1) < npages0 && (pgsz0[p] + pgsz0[p+1] <= DELFACTOR * pl->pagesize)) {
    int i;
    /* 1 copy data to end of this page, and remove next page data */
    memcpy(pl->pgdata[p] + pgsz0[p], pl->pgdata[p+1], pgsz0[p+1] * sizeof(PLTYPE));
    Safefree(pl->pgdata[p+1]);
    if ( (p+1) < (npages0-1) )
      memmove(pl->pgdata + p + 1, pl->pgdata + p + 2, (npages0-1-p) * sizeof(PLTYPE*));
    /* 2 adjust upper levels, moving sizes on boundaries */
    for (i = 1; i < ELEV; i++) {
      PLTYPE ep, npagesi = pl->npages[i], *pgszi = pl->pgsize[i];
      for (ep = 1 + ((p+0) >> (i*ESHIFT)); ep < npagesi; ep++) {
        PLTYPE amt = pgsz0[ep << (i*ESHIFT)];
        pgszi[ep] -= amt;
        pgszi[ep-1] += amt;
      }
      pl->npages[i] = NEPG(npages0-1, i); /* Possibly lower now */
    }
    /* 3 Move sizes at base level over efficiently */
    pgsz0[p] += pgsz0[p+1];
    if ( (p+1) < (npages0-1) )
      memmove(pgsz0 + p + 1,  pgsz0 + p + 2, (npages0-1-p) * sizeof(PLTYPE));
    pl->npages[0]--;
  }
#endif
  CHECKPL(pl, "delete");
}


typedef struct pagelist64_iter_t {
  pagelist64_t *pl;
  PLTYPE p;
  PLTYPE idx;
} pagelist64_iter_t;

static pagelist64_iter_t pagelist64_iterator_create(pagelist64_t *pl, PLTYPE idx) {
  pagelist64_iter_t iter;
  iter.pl = pl;
  iter.p = _pagelist64_find_page(pl, &idx);
  iter.idx = idx;
  return iter;
}

static PLTYPE pagelist64_iterator_next(pagelist64_iter_t *iter) {
  PLTYPE v, p = iter->p;

  if (p >= iter->pl->npages[0]) return 0;

  v = iter->pl->pgdata[p][iter->idx];

  if (++iter->idx >= iter->pl->pgsize[0][p]) {
    iter->p++;
    iter->idx = 0;
  }
  return v;
}

#undef PLTYPE
#undef PLDEBUG
#undef PAGESIZE
#undef ELEV
#undef ESHIFT
#undef DELFACTOR
#undef NEPG
#undef ADDSIZE
#undef CHECKPL

#endif
