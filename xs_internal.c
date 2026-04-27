#define PERL_NO_GET_CONTEXT 1

#include <stdlib.h>
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#define NEED_sv_2pv_flags
#define NEED_HvNAME_get
#include "ppport.h"

#include "ptypes.h"
#include "sort.h"
#include "xs_internal.h"

#if 0
UV neg_iv(UV n) {
  if ((IV)n == IV_MIN)  return (UV_MAX >> 1) + 1;
  else                  return (UV) (-(IV)n);
}
#endif

/* Given 'a' and astatus (-1 means 'a' is an IV), properly mod with n */
void _mod_with(UV *a, int astatus, UV n) {
  if (n == 0) return;
  if (astatus != -1) {
    *a %= n;
  } else {
    UV r = neg_iv(*a) % n;
    *a = (r == 0) ? 0 : n-r;
  }
}

/* Given values and a sign indicating IV or UV, returns -1 (<), 0 (eq), 1 (>) */
int _sign_cmp(int xsign, UV x, int ysign, UV y) {
  /* Convert sign to -1 (neg), 0 (small pos), 1 (big pos) */
  if (x <= (UV)IV_MAX) xsign = 0;
  if (y <= (UV)IV_MAX) ysign = 0;
  if (xsign == ysign && x == y) return 0;
  /* neg < small pos < big pos */
  if (xsign != ysign) return (xsign < ysign) ? -1 : 1;
  /* Numerical comparison as IV or UV */
  return ((xsign == -1 && (IV)x < (IV)y) || (xsign != -1 && x < y)) ? -1 : 1;
}

int _sv_is_bigint(pTHX_ SV* n) {
  if (SvROK(n) && SvOBJECT(SvRV(n))) {
    const char *hvname = HvNAME_get(SvSTASH(SvRV(n)));
    if (hvname != 0 &&
        (strEQ(hvname, "Math::BigInt") || /* BigFloat not here, force to PP */
         strEQ(hvname, "Math::GMPz")   || strEQ(hvname, "Math::GMP") ||
         strEQ(hvname, "Math::GMPq")   || strEQ(hvname, "Math::AnyNum") ||
         strEQ(hvname, "Math::Pari")   || strEQ(hvname, "Math::BigInt::Lite")))
      return 1;
  }
  return 0;
}

int _sv_is_bigint_fast(pTHX_ SV* n) {
  if (SvROK(n) && SvOBJECT(SvRV(n))) {
    const char *hvname = HvNAME_get(SvSTASH(SvRV(n)));
    if (hvname != 0 &&
        (strEQ(hvname, "Math::GMPz")   || strEQ(hvname, "Math::GMP")))
      return 1;
  }
  return 0;
}

#if BITS_PER_WORD == 32
  static const unsigned int uvmax_maxlen = 10;
  static const unsigned int ivmax_maxlen = 10;
  static const char uvmax_str[] = "4294967295";
  /* static const char ivmax_str[] = "2147483648"; */
  static const char ivmin_str[] = "2147483648";
#else
  static const unsigned int uvmax_maxlen = 20;
  static const unsigned int ivmax_maxlen = 19;
  static const char uvmax_str[] = "18446744073709551615";
  /* static const char ivmax_str[] =  "9223372036854775808"; */
  static const char ivmin_str[] =  "9223372036854775808";
#endif

/* Parse any numeric string.  Returns:
 *   SNUMFLAG_UV               non-negative integer, fits in UV
 *   SNUMFLAG_NEG              negative integer, fits in IV
 *   SNUMFLAG_BIGINT [|NEG]    integer too large for native word
 *   SNUMFLAG_FP    [|NEG]     valid floating-point (non-integer) number
 *   SNUMFLAG_INVALID          not a valid number
 *
 * There is a good argument that we should be using Perl's numeric.c functions.
 */
uint32_t _parse_strnum(const char* s, STRLEN len)
{
  STRLEN i = 0, sig_start;
  STRLEN maxlen;
  const char* maxstr;
  uint32_t flag = 0;
  int had_zeros = 0;

  if (s == 0 || len == 0) return SNUMFLAG_UV;   /* null/empty -> 0 */

  /* Sign */
  if      (s[i] == '-') { flag |= SNUMFLAG_NEG; i++; }
  else if (s[i] == '+') { i++; }
  if (i > 0 && i == len) return SNUMFLAG_INVALID;  /* lone sign character */

  /* Strip leading zeros, noting if any existed */
  while (i < len && s[i] == '0') { i++; had_zeros = 1; }
  if (i == len) return SNUMFLAG_UV;      /* zero */
  sig_start = i;

  /* Scan integer digits */
  while (i < len && isDIGIT(s[i])) i++;

  if (i == len) {
    /* Pure integer: range check */
    STRLEN sig_len = i - sig_start;
    if (flag & SNUMFLAG_NEG) { maxlen = ivmax_maxlen; maxstr = ivmin_str; }
    else                     { maxlen = uvmax_maxlen; maxstr = uvmax_str; }
    if (sig_len > maxlen) return flag | SNUMFLAG_BIGINT;
    if (sig_len == maxlen) {
      STRLEN j;
      for (j = 0; j < maxlen; j++)
        if (s[sig_start + j] != maxstr[j]) break;
      if (j < maxlen && s[sig_start + j] > maxstr[j])
        return flag | SNUMFLAG_BIGINT;
    }
    return flag;   /* SNUMFLAG_UV or SNUMFLAG_NEG */
  }

  /* Not a pure integer - try to parse as float */
  /*   [+-]? digit* (. digit*)? ([eE] [+-]? digit+)? */
  {
    int has_frac = 0;
    int has_int  = (i > sig_start || had_zeros);  /* had any integer digits? */

    if (i < len && s[i] == '.') {
      i++;
      while (i < len && isDIGIT(s[i])) { i++; has_frac = 1; }
    }

    /* Reject lone ".", ".e5", bare "e5", etc. */
    if (!has_int && !has_frac) return SNUMFLAG_INVALID;

    if (i < len && (s[i] == 'e' || s[i] == 'E')) {
      i++;
      if (i < len && (s[i] == '+' || s[i] == '-')) i++;
      if (i >= len || !isDIGIT(s[i])) return SNUMFLAG_INVALID;
      while (i < len && isDIGIT(s[i])) i++;
    }

    return (i == len) ? flag | SNUMFLAG_FP : SNUMFLAG_INVALID;
  }
}

/* Is this a pedantically valid integer?
 * Croaks if undefined or invalid.
 * Returns 0 if it is an object or a string too large for a UV.
 * Returns 1/-1 if it is good to process by XS.
 * TODO: it would be useful to know the sign even if returning 0 for bigint.
 */
int _validate_int(pTHX_ SV* n, int negok)
{
  const char* mustbe = (negok) ? "must be an integer" : "must be a non-negative integer";
  const char* sptr;
  STRLEN len;
  uint32_t stype, isbignum = 0;

  /* TODO: magic, grok_number, etc. */
  if (SVNUMTEST(n)) { /* If defined as number, use it */
    if (SvIsUV(n) || SvIVX(n) >= 0)  return 1; /* The normal case */
    if (negok)  return -1;
    else croak("Parameter '%" SVf "' %s", n, mustbe);
  }
  if (sv_isobject(n)) {
    isbignum = _sv_is_bigint(aTHX_ n);
    if (!isbignum) return 0;
  }
  if (!SvOK(n))  croak("Parameter must be defined");
  if (SvGAMAGIC(n) && !isbignum)   sptr = SvPV(n, len);
  else                             sptr = SvPV_nomg(n, len);
  if (len == 0 || sptr == 0)  croak("Parameter %s", mustbe);
  stype = _parse_strnum(sptr, len);
  if (stype == SNUMFLAG_UV)      return 1;
  if (negok || !(stype & SNUMFLAG_NEG)) {
    if (stype == SNUMFLAG_NEG)   return -1;
    if (stype & SNUMFLAG_BIGINT) return 0;
  }
  croak("Parameter '%" SVf "' %s", n, mustbe);
}

int _validate_and_set(UV* val, pTHX_ SV* svn, uint32_t mask) {
  int status;

  if (svn == 0) croak("Parameter must be defined");
  /* Streamline the typical path of input being a native integer. */
  if (SVNUMTEST(svn)) {
    IV n = SvIVX(svn);
    if (n >= 0) {
      if (n == 0 && (mask & IFLAG_POS))
        croak("Parameter '%" SVf "' must be a positive integer", svn);
      *val = (UV)n;
      return 1;
    }
    if (SvIsUV(svn)) {
      if (mask & IFLAG_IV)
        return 0;
      *val = (UV)n;
      return 1;
    }
    if (mask & IFLAG_ABS)    { *val = (UV)(-n); return 1; }
    if (mask & IFLAG_POS)    croak("Parameter '%" SVf "' must be a positive integer", svn);
    if (mask & IFLAG_NONNEG) croak("Parameter '%" SVf "' must be a non-negative integer", svn);
    *val = n;
    return -1;
  }

  status = _validate_int( aTHX_ svn, !(mask & (IFLAG_NONNEG|IFLAG_POS)) );
  if (status == 1) {
    UV n = my_svuv(svn);
    if (n == 0 && (mask & IFLAG_POS))
      croak("Parameter '%" SVf "' must be a positive integer", svn);
    if (n > (UV)IV_MAX && (mask & IFLAG_IV))
      return 0;
    *val = n;
  } else if (status == -1) {
    IV n = my_sviv(svn);
    if (mask & IFLAG_ABS) { *val = (UV)(-n); status = 1; }
    else                  { *val = (UV)n; }
  }
  return status;
}

SV* _fetch_arref(pTHX_ AV* av, SV** svarr, size_t i) {
  if (svarr == 0) {
    SV **svp = av_fetch(av, i, 0);
    return svp ? *svp : &PL_sv_undef;
  }
  return svarr[i];
}

#define READ_UV_IARR(dst, src, itype) \
  { \
    UV n; \
    int istatus = _validate_and_set(&n, aTHX_ src, IFLAG_ANY); \
    if      (istatus == -1)                  itype |= IARR_TYPE_NEG; \
    else if (istatus == 1 && n > (UV)IV_MAX) itype |= IARR_TYPE_POS; \
    if (istatus == 0 || itype == IARR_TYPE_BAD) break; \
    dst = n; \
  }

int arrayref_to_int_array(pTHX_ size_t *retlen, UV** ret, bool want_sort, SV* sva, const char* fstr)
{
  Size_t len, i;
  int itype = IARR_TYPE_ANY;
  UV  *r;
  DECL_ARREF(avp);

  USE_ARREF(avp, sva, fstr, AR_READ);
  len = len_avp;
  *retlen = len;
  if (len == 0) {
    *ret = 0;
    return itype;
  }
  New(0, r, len, UV);
  for (i = 0; i < len; i++) {
    SV *iv = FETCH_ARREF(avp,i);
    if (iv == 0) continue;
    if (SVNUMTEST(iv)) {
      IV n = SvIVX(iv);
      if (n < 0) {
        if (SvIsUV(iv))  itype |= IARR_TYPE_POS;
        else             itype |= IARR_TYPE_NEG;
        if (itype == IARR_TYPE_BAD) break;
      }
      r[i] = (UV)n;
    } else {
      READ_UV_IARR(r[i], iv, itype);
    }
  }
  if (i < len) {
    Safefree(r);
    *ret = 0;
    return IARR_TYPE_BAD;
  }
  *ret = r;
  if (want_sort) {
    if (itype == IARR_TYPE_NEG) {
      for (i = 1; i < len; i++)
        if ( (IV)r[i] <= (IV)r[i-1] )
          break;
    } else {
      for (i = 1; i < len; i++)
        if (r[i] <= r[i-1])
          break;
    }
    if (i < len)
      sort_dedup_uv_array(r, itype == IARR_TYPE_NEG, retlen);
  }
  return itype;
}

/* Check whether an SV is a non-magical arrayref whose elements are all native
 * non-negative integers in strictly increasing order (i.e. sorted and unique).
 * On success returns the AvARRAY pointer and sets *lenp; otherwise NULL.
 * Used by the set-op fast path to skip intermediate UV array allocation. */
SV** _check_sorted_nonneg_arrayref(pTHX_ SV *sv, size_t *lenp)
{
  AV *av;
  SV **arr;
  size_t len, i;
  if (!SvROK(sv) || SvTYPE(SvRV(sv)) != SVt_PVAV) return NULL;
  av = (AV*)SvRV(sv);
  if (SvMAGICAL(av)) return NULL;
  arr = AvARRAY(av);
  len = av_count(av);
  for (i = 0; i < len; i++) {
    SV *elem = arr[i];
    if (!SVNUMTEST(elem) || (!SvIsUV(elem) && SvIVX(elem) < 0)) return NULL;
    if (i > 0 && SvUVX(elem) <= SvUVX(arr[i-1]))                return NULL;
  }
  *lenp = len;
  return arr;
}

int array_to_int_array(pTHX_ size_t *retlen, UV** ret, bool want_sort, SV** svbase, size_t len)
{
  size_t i;
  int itype = IARR_TYPE_ANY;
  UV  *r;
  *retlen = len;
  if (len == 0) {
    *ret = 0;
    return itype;
  }
  New(0, r, len, UV);
  for (i = 0; i < len; i++) {
    SV *iv = svbase[i];
    if (SVNUMTEST(iv)) {
      IV n = SvIVX(iv);
      if (n < 0) {
        if (SvIsUV(iv))  itype |= IARR_TYPE_POS;
        else             itype |= IARR_TYPE_NEG;
        if (itype == IARR_TYPE_BAD) break;
      }
      r[i] = (UV)n;
    } else {
      READ_UV_IARR(r[i], iv, itype);
    }
  }
  if (i < len) {
    Safefree(r);
    *ret = 0;
    return IARR_TYPE_BAD;
  }
  *ret = r;
  if (want_sort) {
    if (itype == IARR_TYPE_NEG) {
      for (i = 1; i < len; i++)
        if ( (IV)r[i] <= (IV)r[i-1] )
          break;
    } else {
      for (i = 1; i < len; i++)
        if (r[i] <= r[i-1])
          break;
    }
    if (i < len)
      sort_dedup_uv_array(r, itype == IARR_TYPE_NEG, retlen);
  }
  return itype;
}

int arrayref_to_digit_array(pTHX_ UV** ret, AV* av, int base)
{
  SSize_t len, i;
  UV *r, carry = 0;
  if (SvTYPE((SV*)av) != SVt_PVAV)
    croak("fromdigits first argument must be a string or array reference");
  len = av_count(av);
  New(0, r, len, UV);
  for (i = len-1; i >= 0; i--) {
    SV** psvd = av_fetch(av, i, 0);
    if (_validate_and_set(r+i, aTHX_ *psvd, IFLAG_ANY) != 1) break;
    r[i] += carry;
    if (r[i] >= (UV)base && i > 0) {
      carry = r[i] / base;
      r[i] -= carry * base;
    } else {
      carry = 0;
    }
  }
  if (i >= 0) {
    Safefree(r);
    return -1;
  }
  /* printf("array is ["); for(i=0;i<len;i++)printf("%lu,",r[i]); printf("]\n"); */
  *ret = r;
  return len;
}

#undef READ_UV_IARR

int _compare_array_refs(pTHX_ SV* a, SV* b)
{
  AV *ava, *avb;
  SSize_t i, alen, blen;
  if ( ((!SvROK(a)) || (SvTYPE(SvRV(a)) != SVt_PVAV)) ||
       ((!SvROK(b)) || (SvTYPE(SvRV(b)) != SVt_PVAV)) )
    return -1;
  ava = (AV*) SvRV(a);
  avb = (AV*) SvRV(b);
  alen = av_len(ava);
  blen = av_len(avb);
  if (alen != blen)
    return 0;
  for (i = 0; i <= alen; i++) {
    SV** iva = av_fetch(ava, i, 0);
    SV** ivb = av_fetch(avb, i, 0);
    SV *sva, *svb;
    int res;

    if (!iva || !ivb)  return -1;
    sva = *iva;
    svb = *ivb;

    if (!SvOK(sva) && !SvOK(svb))  /* Two undefs are fine. */
      continue;
    if (!SvOK(sva) || !SvOK(svb))  /* One undef isn't ok. */
      return 0;
    /* Hashes, I/O, etc. are not ok. */
    if (SvTYPE(sva) >= SVt_PVAV || SvTYPE(svb) >= SVt_PVAV)
      return -1;

    /* One of them is a non-object reference */
    if ( (SvROK(sva) && !sv_isobject(sva)) ||
         (SvROK(svb) && !sv_isobject(svb)) ) {
      /* Always error if either one is not an array reference. */
      if ( (SvROK(sva) && SvTYPE(SvRV(sva)) != SVt_PVAV) ||
           (SvROK(svb) && SvTYPE(SvRV(svb)) != SVt_PVAV) )
        return -1;
      /* One reference, one non-reference = not equal */
      if (SvROK(sva) != SvROK(svb))
        return 0;
      /* Now we know both are array references.  Compare. */
      res = _compare_array_refs(aTHX_ sva, svb);
      if (res == 1) continue;
      return res;
    }

    /* Common case: two simple integers */
    if (    SVNUMTEST(sva) && SVNUMTEST(svb)
         && (SvTYPE(sva) == SVt_IV || SvTYPE(sva) == SVt_PVIV)
         && (SvTYPE(svb) == SVt_IV || SvTYPE(svb) == SVt_PVIV) ) {
      UV va = my_svuv(sva), vb = my_svuv(svb);
      if (va != vb) return 0;
      continue;
    }

    /* This function is more useful if we allow more than strictly integers */
    {  /* Compare the string representation */
      STRLEN alen, blen;
      const char* stra = SvPV(sva, alen);
      const char* strb = SvPV(svb, blen);
      if (alen != blen || strcmp(stra,strb) != 0)
        return 0;
    }
  }
  return 1;
}
