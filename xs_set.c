#define PERL_NO_GET_CONTEXT 1

#include <string.h>
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#define NEED_sv_2pv_flags
#include "ppport.h"

#include "ptypes.h"
#include "util.h"
#include "ds_iset.h"
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
  if ((IV)amax > 0 && (IV)bmax > 0) {
    if (UV_MAX-amax < bmax || amax + bmax > (UV)IV_MAX)
      return IARR_TYPE_BAD;  /* overflow */
  }
  if ((IV)amin < 0 && (IV)bmin < 0) {
    UV aneg = neg_iv(amin);
    UV bneg = neg_iv(bmin);
    if (UV_MAX-aneg < bneg || aneg + bneg > (UV)IV_MAX)
      return IARR_TYPE_BAD;  /* underflow */
  }
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

bool xs_set_relation(pTHX_ SV* sva, SV* svb, set_relation_op_t op, int *ret, const char *name)
{
  int atype, btype;
  UV *ra, *rb;
  size_t alen, blen, inalen, inblen;

  if (op < SET_REL_DISJOINT || op > SET_REL_PROPER_INTERSECTION)
    croak("%s: unknown set relation", name);

  /* If one set is much smaller than the other, it would be faster using
   * is_in_set().  We'll keep things simple and slurp in both sets. */

  /* THIS ASSUMES THE INPUT LISTS HAVE NO DUPLICATES */
  inalen = inblen = 0;
  if (SvROK(sva) && SvTYPE(SvRV(sva)) == SVt_PVAV && SvROK(svb) && SvTYPE(SvRV(svb)) == SVt_PVAV) {
    /* Shortcut on length if we can to skip intersection. */
    inalen = av_count((AV*) SvRV(sva));
    inblen = av_count((AV*) SvRV(svb));
    if ( (op == SET_REL_EQUAL           && inalen != inblen) ||
         (op == SET_REL_SUBSET          && inalen <  inblen) ||
         (op == SET_REL_PROPER_SUBSET   && inalen <= inblen) ||
         (op == SET_REL_SUPERSET        && inalen >  inblen) ||
         (op == SET_REL_PROPER_SUPERSET && inalen >= inblen) ) {
      *ret = 0;
      return 1;
    }
  }

  /* Get the integers as sorted arrays of IV or UV */
  atype = arrayref_to_int_array(aTHX_ &alen, &ra, 1, sva, name);
  btype = arrayref_to_int_array(aTHX_ &blen, &rb, 1, svb, name);

  if (CAN_COMBINE_IARR_TYPES(atype,btype)) {
    size_t rlen = 0, ia = 0, ib = 0;
    int pcmp = (atype == IARR_TYPE_NEG || btype == IARR_TYPE_NEG) ? 0 : 1;

    while (ia < alen && ib < blen) {
      if (ra[ia] == rb[ib]) {
        rlen++;
        ia++; ib++;
      } else {
        if (SIGNED_CMP_LT(pcmp, ra[ia], rb[ib])) ia++;
        else                                     ib++;
      }
    }
    Safefree(ra);
    Safefree(rb);

    *ret = 0;
    switch (op) {
      case SET_REL_DISJOINT:
        *ret = (rlen == 0);
        break;
      case SET_REL_EQUAL:
        *ret = (alen == blen && rlen == blen);
        break;
      case SET_REL_SUBSET:
        *ret = (alen >= blen && rlen == blen);
        break;
      case SET_REL_PROPER_SUBSET:
        *ret = (alen >  blen && rlen == blen);
        break;
      case SET_REL_SUPERSET:
        *ret = (alen <= blen && rlen == alen);
        break;
      case SET_REL_PROPER_SUPERSET:
        *ret = (alen <  blen && rlen == alen);
        break;
      case SET_REL_PROPER_INTERSECTION:
        *ret = (rlen > 0 && rlen < alen && rlen < blen);
        break;
    }
    return 1;
  }
  Safefree(ra);
  Safefree(rb);
  return 0;
}

static SV* set_arrayref_from_uv_array(pTHX_ size_t len, UV *data, int sign)
{
  size_t i;
  AV *av = newAV();
  if (len > 0) {
    SV **ar;
    av_extend(av, (SSize_t)len - 1);
    ar = AvARRAY(av);
    for (i = 0; i < len; i++)
      ar[i] = NEWSVINT(sign, data[i]);
    AvFILLp(av) = (SSize_t)len - 1;
  }
  Safefree(data);
  return newRV_noinc((SV*)av);
}

static SV* set_arrayref_from_sv_merge(pTHX_ SV **aa, size_t alen, SV **bb, size_t blen, set_op_t op)
{
  int inc_eq = (op == SET_OP_UNION || op == SET_OP_INTERSECT);
  int inc_lt = (op != SET_OP_INTERSECT);
  int inc_gt = (op == SET_OP_UNION || op == SET_OP_DELTA);
  size_t maxlen = (op == SET_OP_INTERSECT) ? (alen < blen ? alen : blen) : alen + blen;
  AV *av = newAV();
  size_t rlen = 0, ia = 0, ib = 0;
  SV **ar;

  if (maxlen > 0)
    av_extend(av, (SSize_t)maxlen - 1);
  ar = AvARRAY(av);

  while (ia < alen && ib < blen) {
    UV va = SvUVX(aa[ia]), vb = SvUVX(bb[ib]);
    if (va == vb) {
      if (inc_eq) ar[rlen++] = SvREFCNT_inc(aa[ia]);
      ia++;
      ib++;
    } else if (va < vb) {
      if (inc_lt) ar[rlen++] = SvREFCNT_inc(aa[ia]);
      ia++;
    } else {
      if (inc_gt) ar[rlen++] = SvREFCNT_inc(bb[ib]);
      ib++;
    }
  }
  if (inc_lt) while (ia < alen) ar[rlen++] = SvREFCNT_inc(aa[ia++]);
  if (inc_gt) while (ib < blen) ar[rlen++] = SvREFCNT_inc(bb[ib++]);
  AvFILLp(av) = (SSize_t)rlen - 1;

  return newRV_noinc((SV*)av);
}

bool xs_set_op(pTHX_ SV* sva, SV* svb, set_op_t op, SV **ret, const char *name)
{
  int atype, btype;
  UV *ra, *rb;
  size_t alen, blen;

  if (op < SET_OP_UNION || op > SET_OP_DELTA)
    croak("%s: unknown set operation", name);

  /* Fast path: both inputs are non-magical arrayrefs of native non-negative
   * sorted unique integers.  Merge SV* directly, preserving existing values
   * and skipping intermediate UV arrays and per-element newSVuv calls. */
  {
    size_t fa, fb;
    SV **aa = _check_sorted_nonneg_arrayref(aTHX_ sva, &fa);
    SV **bb = aa ? _check_sorted_nonneg_arrayref(aTHX_ svb, &fb) : NULL;
    if (aa && bb) {
      *ret = set_arrayref_from_sv_merge(aTHX_ aa, fa, bb, fb, op);
      return 1;
    }
  }

  /* Get the integers and ensure they are sorted unique integers first. */
  atype = arrayref_to_int_array(aTHX_ &alen, &ra, 1, sva, name);
  btype = arrayref_to_int_array(aTHX_ &blen, &rb, 1, svb, name);

  if (CAN_COMBINE_IARR_TYPES(atype,btype)) {
    UV *r = 0;
    size_t rlen = 0, ia = 0, ib = 0;
    int pcmp = (atype == IARR_TYPE_NEG || btype == IARR_TYPE_NEG) ? 0 : 1;

    if (op == SET_OP_UNION) {
      New(0, r, alen + blen, UV);
      while (ia < alen && ib < blen) {
        if (ra[ia] == rb[ib]) {
          r[rlen++] = ra[ia];
          ia++;
          ib++;
        } else {
          if (SIGNED_CMP_LT(pcmp, ra[ia], rb[ib])) r[rlen++] = ra[ia++];
          else                                     r[rlen++] = rb[ib++];
        }
      }
      if (ia < alen) { Copy(ra+ia, r+rlen, alen-ia, UV); rlen += alen-ia; }
      if (ib < blen) { Copy(rb+ib, r+rlen, blen-ib, UV); rlen += blen-ib; }
    } else if (op == SET_OP_INTERSECT) {
      New(0, r, (alen < blen) ? alen : blen, UV);
      while (ia < alen && ib < blen) {
        if (ra[ia] == rb[ib]) {
          r[rlen++] = ra[ia];
          ia++;
          ib++;
        } else {
          if (SIGNED_CMP_LT(pcmp, ra[ia], rb[ib])) ia++;
          else                                     ib++;
        }
      }
    } else if (op == SET_OP_MINUS) {
      New(0, r, alen, UV);
      while (ia < alen && ib < blen) {
        if (ra[ia] == rb[ib]) {
          ia++;
          ib++;
        } else {
          if (SIGNED_CMP_LT(pcmp, ra[ia], rb[ib])) r[rlen++] = ra[ia++];
          else                                     ib++;
        }
      }
      if (ia < alen) { Copy(ra+ia, r+rlen, alen-ia, UV); rlen += alen-ia; }
    } else if (op == SET_OP_DELTA) {
      New(0, r, alen + blen, UV);
      while (ia < alen && ib < blen) {
        if (ra[ia] == rb[ib]) {
          ia++;
          ib++;
        } else {
          if (SIGNED_CMP_LT(pcmp, ra[ia], rb[ib])) r[rlen++] = ra[ia++];
          else                                     r[rlen++] = rb[ib++];
        }
      }
      if (ia < alen) { Copy(ra+ia, r+rlen, alen-ia, UV); rlen += alen-ia; }
      if (ib < blen) { Copy(rb+ib, r+rlen, blen-ib, UV); rlen += blen-ib; }
    }
    Safefree(ra);
    Safefree(rb);
    *ret = set_arrayref_from_uv_array(aTHX_ rlen, r, pcmp);
    return 1;
  }
  Safefree(ra);
  Safefree(rb);
  return 0;
}

bool xs_is_sidon_set(pTHX_ SV* sva, int *ret)
{
  int itype, is_sidon;
  size_t len, i, j;
  UV *data;
  iset_t s;

  itype = arrayref_to_int_array(aTHX_ &len, &data, 1, sva, "is_sidon_set");
  if (itype == IARR_TYPE_NEG) {  /* All elements must be non-negative. */
    Safefree(data);
    *ret = 0;
    return 1;
  }
  /* If any bigints or we cannot add the values in 64-bits, call PP. */
  if (itype == IARR_TYPE_BAD || itype == IARR_TYPE_POS) {
    Safefree(data);
    return 0;
  }

  is_sidon = 1;
  s = iset_create(20UL * len);
  for (i = 0; i < len && is_sidon; i++)
    for (j = i; j < len; j++)
      if (!iset_add(&s, data[i] + data[j], 1))
        { is_sidon = 0; break; }
  Safefree(data);
  iset_destroy(&s);

  *ret = is_sidon;
  return 1;
}

bool xs_is_sumfree_set(pTHX_ SV* sva, int *ret)
{
  UV *data;
  size_t len, i, j;
  int itype;
  bool is_sumfree;

  itype = arrayref_to_int_array(aTHX_ &len, &data, 1, sva, "is_sumfree_set");
  if (itype != IARR_TYPE_BAD && len <= 1) { /* Degenerate cases: len 0 or 1 */
    *ret = (len == 0 || data[0] != 0);
    Safefree(data);
    return 1;
  }
  /* Check for IV overflow on sum */
  if (itype == IARR_TYPE_NEG) {
    IV min = data[0], max = data[len-1];  /* Array is sorted */
    if (min < IV_MIN/2 || max > IV_MAX/2)
      itype = IARR_TYPE_BAD;
  }
  is_sumfree = 1;
  if (itype == IARR_TYPE_ANY) {
    for (i = 0; i < len && is_sumfree; i++)
      for (j = i; j < len; j++)
        if (is_in_sorted_uv_array(data[i]+data[j], data, len))
          { is_sumfree = 0; break; }
  } else if (itype == IARR_TYPE_NEG) {
    for (i = 0; i < len && is_sumfree; i++)
      for (j = i; j < len; j++)
        if (is_in_sorted_iv_array((IV)data[i]+(IV)data[j], (IV*)data, len))
          { is_sumfree = 0; break; }
  }
  Safefree(data);

  if (itype == IARR_TYPE_ANY || itype == IARR_TYPE_NEG) {
    *ret = is_sumfree;
    return 1;
  }

  /* We're here because one of:
   *   1) itype is TYPE_BAD because there were bigints.
   *   2) itype is TYPE_BAD because summed IVs would overflow.
   *   3) itype is TYPE_POS.
   *      At least one element is >= 2^63, so we would overflow on sum.
   */
  return 0;
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
