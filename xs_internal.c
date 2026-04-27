#define PERL_NO_GET_CONTEXT 1

#include <stdlib.h>
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#define NEED_sv_2pv_flags
#define NEED_HvNAME_get
#include "ppport.h"

#include "ptypes.h"
#include "xs_internal.h"

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
