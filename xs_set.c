#define PERL_NO_GET_CONTEXT 1

#include <string.h>
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#define NEED_sv_2pv_flags
#include "ppport.h"

#include "ptypes.h"
#include "xs_internal.h"
#include "xs_set.h"

#if PERL_VERSION_LT(5,17,7)
#  define SvREFCNT_dec_NN(sv)    SvREFCNT_dec(sv)
#endif

#define NEWSVINT(sign,v) (((sign) > 0) ? newSVuv(v) : newSViv(v))

int type_of_sumset(int typea, int typeb, UV amin, UV amax, UV bmin, UV bmax) {
  if (typea == IARR_TYPE_BAD || typeb == IARR_TYPE_BAD)
    return IARR_TYPE_BAD;
  if (typea != IARR_TYPE_NEG && typeb != IARR_TYPE_NEG) {
    /* ANY+ANY  ANY+POS  POS+ANY  POS+POS */
    if (UV_MAX-amax < bmax)     return IARR_TYPE_BAD;
    if (amax+bmax > (UV)IV_MAX) return IARR_TYPE_POS;
    return IARR_TYPE_ANY;
  }
  /* For simplicity, throw out NEG+POS to avoid UV+IV */
  if (typea == IARR_TYPE_POS || typeb == IARR_TYPE_POS)
    return IARR_TYPE_BAD;
  /* NEG+NEG  NEG+ANY  ANY+NEG */
  if ((IV)amax > 0 && (IV)bmax > 0 && amax + bmax > (UV)IV_MAX)
    return IARR_TYPE_BAD;  /* overflow */
  if ((IV)amin < 0 && (IV)bmin < 0 && (UV)(-(IV)amin) + (UV)(-(IV)bmin) > (UV)IV_MAX)
    return IARR_TYPE_BAD;  /* underflow */
  if (((IV)amin > 0 || (IV)bmin > 0) && (IV)(amin+bmin) >= 0)
    return IARR_TYPE_ANY;  /* Result is all positive */
  return IARR_TYPE_NEG;
}

void _sc_clear_cache(set_data_t *cache) {
  memset(cache->status, 0, sizeof(signed char) * (2+MPU_SC_SIZE));
}

#define _SC_GET_VALUE(statvar, var, arr, i) \
  statvar = _validate_and_set(&var, aTHX_ arr[i], IFLAG_ANY); \
  if (statvar == 0) return -1;

#define SC_SET_MID_VALUE(statvar, var, arr, i, cache) \
  do { \
    if (cache == 0) { \
      _SC_GET_VALUE(statvar, var, arr, i) \
    } else { \
      unsigned int imod_ = 2 + ((i) % MPU_SC_SIZE); \
      if (cache->status[imod_] != 0 && cache->index[imod_] == i) { \
        statvar = cache->status[imod_]; \
        var     = cache->value[imod_]; \
      } else { \
        _SC_GET_VALUE(statvar, var, arr, i) \
        cache->status[imod_] = statvar; \
        cache->value[imod_]  = var; \
        cache->index[imod_]  = i; \
      } \
    } \
  } while (0)

int _sc_set_lohi(pTHX_ SV** avarr, set_data_t *cache, int loindex, int hiindex, int *lostatus, int *histatus, UV *loval, UV *hival)
{
  if (cache && cache->status[0] != 0) {
    *lostatus = cache->status[0];  *loval = cache->value[0];
  } else {
    _SC_GET_VALUE(*lostatus, *loval, avarr, loindex);
    if (cache) {
      cache->status[0] = *lostatus;
      cache->value[0]  = *loval;
    }
  }
  if (cache && cache->status[1] != 0) {
    *histatus = cache->status[1];  *hival = cache->value[1];
  } else {
    _SC_GET_VALUE(*histatus, *hival, avarr, hiindex);
    if (cache) {
      cache->status[1] = *histatus;
      cache->value[1] = *hival;
    }
  }
  return 1;
}



/* index of val in a set (array ref of sorted unique integers)
 *    -1 bigint
 *     n nth-position (0 .. count-1)
 *  eq will be set to 1 if the element in that position is the input value.
 */
static int index_for_set(pTHX_ AV* av, set_data_t *cache, int sign, UV val, int *eq)
{
  Size_t len;
  int lo, hi, lostatus, histatus, midstatus, cmp;
  UV  rlo, rhi, rmid;
  SV** arr;

  if (sign != 1 && sign != -1)
    return -1;
  len = av_count(av);
  if (len == 0) {
    *eq = 0;
    return 0;
  }
  arr = AvARRAY(av);

  lo = 0;
  hi = len-1;
  if (_sc_set_lohi(aTHX_ arr, cache, lo, hi, &lostatus, &histatus, &rlo, &rhi) < 0)
    return -1;

  cmp = _sign_cmp(sign, val, lostatus, rlo);
  if (cmp <= 0) { *eq = cmp==0; return lo; }
  /* val > rlo */
  cmp = _sign_cmp(sign, val, histatus, rhi);
  if (cmp >= 0) { *eq = cmp==0; return hi + (cmp>0); }
  /* val < rhi */

  while (hi-lo > 1) {
    int mid = lo + ((hi-lo) >> 1);
    SC_SET_MID_VALUE(midstatus, rmid, arr, (size_t)mid, cache);
    cmp = _sign_cmp(midstatus, rmid, sign, val);
    if (cmp == 0) { *eq = 1; return mid; }
    if (cmp < 0) { lo = mid; rlo = rmid; lostatus = midstatus; }
    else         { hi = mid; rhi = rmid; histatus = midstatus; }
  }
  if (sign == histatus && rhi == val)
    *eq = 1;
  else if (_sign_cmp(sign,val, histatus,rhi) > 0)
    croak("internal index error");
  return hi;
}

/* Find index to insert in a set (array ref of sorted unique integers)
 *    -1 bigint
 *    0  already in set
 *    n  should be in n-th position (1 means should be first element)
 */
int insert_index_in_set(pTHX_ AV* av, set_data_t *cache, int sign, UV val) {
  int eq = 0;
  int index = index_for_set(aTHX_ av, cache, sign, val, &eq);
  return (index < 0) ? index : eq ? 0 : index+1;
}

/* Find index of element in a set (array ref of sorted unique integers)
 *    -1 bigint
 *    0  not in set
 *    n  in n-th position (1 means first element)
 */
int index_in_set(pTHX_ AV* av, set_data_t *cache, int sign, UV val) {
  int eq = 0;
  int index = index_for_set(aTHX_ av, cache, sign, val, &eq);
  return (index < 0) ? index : eq ? index+1 : 0;
}

/* See if an element is in a set (array ref of sorted unique integers) */
/* -1 = bigint, 0 = not found, 1 = found */
int is_in_set(pTHX_ AV* av, set_data_t *cache, int sign, UV val)
{
  int eq = 0;
  int index = index_for_set(aTHX_ av, cache, sign, val, &eq);
  return (index < 0) ? index : eq ? 1 : 0;
}

/* 1 if deleted, 0 if not deleted, -1 if need to punt to PP */
int del_from_set(pTHX_ AV* ava, int bstatus, UV b) {
  int index = index_in_set(aTHX_ ava, 0, bstatus, b);
  if (index <= 0)
    return index;
  {
    SV **arr = AvARRAY(ava);
    SV *savep = arr[index-1];
    Size_t pos = index, alen = av_count(ava);
    if (pos > alen/2) {
      if (pos < alen) {
        memmove(arr+pos-1, arr+pos, sizeof(SV*) * (alen-pos));
        arr[alen-1] = savep;
      }
      SvREFCNT_dec_NN(av_pop(ava));
    } else {
      if (pos > 1) {
        memmove(arr+1, arr+0, sizeof(SV*) * (pos-1));
        arr[0] = savep;
      }
      SvREFCNT_dec_NN(av_shift(ava));
    }
  }
  return 1;
}

/* 1 if inserted, 0 if not inserted, -1 if need to punt to PP */
int ins_into_set(pTHX_ AV* ava, int bstatus, UV b) {
  int index = insert_index_in_set(aTHX_ ava, 0, bstatus, b);
  if (index <= 0)
    return index;
  {
    SV *newb, **arr;
    SV* newsvb = NEWSVINT(bstatus, b);
    Size_t alen = av_count(ava);
    if ((Size_t)index > alen/2) {
      av_push(ava, newsvb);
      if ((Size_t)index <= alen) {
        arr = AvARRAY(ava);
        newb = arr[alen];
        memmove(arr+index, arr+index-1, sizeof(SV*) * (alen-(index-1)));
        arr[index-1] = newb;
      }
    } else {
      av_unshift(ava, 1);
      av_store(ava, 0, newsvb);
      if (index > 1) {
        arr = AvARRAY(ava);
        newb = arr[0];
        memmove(arr+0, arr+1, sizeof(SV*) * index);
        arr[index-1] = newb;
      }
    }
  }
  return 1;
}

#undef NEWSVINT
