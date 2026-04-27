#ifndef MPU_XS_INTERNAL_H
#define MPU_XS_INTERNAL_H

#if BITS_PER_WORD == 64
  #if defined(_MSC_VER)
    #include <stdlib.h>
    #define strtoull _strtoui64
    #define strtoll  _strtoi64
  #endif
  #define PSTRTOULL(str, end, base) strtoull (str, end, base)
  #define PSTRTOLL(str, end, base)  strtoll (str, end, base)
#else
  #define PSTRTOULL(str, end, base) strtoul (str, end, base)
  #define PSTRTOLL(str, end, base)  strtol (str, end, base)
#endif

#if defined(_MSC_VER) && !defined(strtold)
  #define strtold strtod
#endif

#ifdef USE_QUADMATH
  #define STRTONV(t)  strtoflt128(t,NULL)
#elif defined(USE_LONG_DOUBLE) && defined(HAS_LONG_DOUBLE)
  #define STRTONV(t)  strtold(t,NULL)
#else
  #define STRTONV(t)  strtod(t,NULL)
#endif

#if PERL_VERSION_LT(5,7,0) && BITS_PER_WORD == 64
 /* Workaround perl 5.6 UVs and bigints */
 #define my_svuv(sv)  PSTRTOULL(SvPV_nolen(sv), NULL, 10)
 #define my_sviv(sv)  PSTRTOLL(SvPV_nolen(sv), NULL, 10)
#elif PERL_VERSION_LT(5,14,0) && BITS_PER_WORD == 64
 /* Workaround RT 49569 in Math::BigInt::FastCalc (pre 5.14.0) */
 /* TODO: Math::BigInt::Pari has the same problem with negs pre-5.18.0 */
 #define my_svuv(sv) ( (!SvROK(sv)) ? SvUV(sv) : PSTRTOULL(SvPV_nolen(sv),NULL,10) )
 #define my_sviv(sv) ( (!SvROK(sv)) ? SvIV(sv) : PSTRTOLL(SvPV_nolen(sv),NULL,10) )
#else
 #define my_svuv(sv) SvUV(sv)
 #define my_sviv(sv) SvIV(sv)
#endif

#if PERL_VERSION_GE(5,9,4) || PERL_VERSION_EQ(5,8,9)
  #define SVf_MAGTEST  SVf_ROK
#else
  #define SVf_MAGTEST  SVf_AMAGIC
  #define GV_NOTQUAL   0
#endif

#define SVNUMTEST(n) \
  ((SvFLAGS(n) & (SVf_IOK | SVf_MAGTEST | SVs_GMG )) == SVf_IOK)

#define SNUMFLAG_UV       0x00000000U
#define SNUMFLAG_NEG      0x00000001U
#define SNUMFLAG_BIGINT   0x00000002U
#define SNUMFLAG_FP       0x00000004U
#define SNUMFLAG_INVALID  0x00000008U

#define IFLAG_ANY      0x00000000U
#define IFLAG_NONNEG   0x00000001U  /* Must be non-negative */
#define IFLAG_POS      0x00000002U  /* Must be positive (n > 0) */
#define IFLAG_ABS      0x00000004U  /* Absolute value returned */
#define IFLAG_IV       0x00000008U  /* Value returned as IV */

#if 1
  /* This is NEGATE_2UV(iv) from handy.h */
  #define neg_iv(n) ((UV)-((n)+1) + 1U)
#else
UV neg_iv(UV n);
#endif

void _mod_with(UV *a, int astatus, UV n);
int _sign_cmp(int xsign, UV x, int ysign, UV y);

int _sv_is_bigint(pTHX_ SV* n);
int _sv_is_bigint_fast(pTHX_ SV* n);
uint32_t _parse_strnum(const char* s, STRLEN len);
int _validate_int(pTHX_ SV* n, int negok);
int _validate_and_set(UV* val, pTHX_ SV* svn, uint32_t mask);

#endif
