#define PERL_NO_GET_CONTEXT 1 /* Define at top for more efficiency. */

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "multicall.h"  /* only works in 5.6 and newer */
#include <stdio.h>      /* For fileno and stdout */

#define NEED_newCONSTSUB
#define NEED_newRV_noinc
#define NEED_sv_2pv_flags
#define NEED_HvNAME_get
#include "ppport.h"

#define FUNC_gcd_ui 1
#define FUNC_isqrt 1
#define FUNC_ipow 1
#define FUNC_popcnt 1
#define FUNC_is_strong_pseudoprime 1
#include "ptypes.h"
#include "cache.h"
#include "sieve.h"
#include "sieve_cluster.h"
#include "util.h"
#include "primality.h"
#include "lucas_seq.h"
#include "factor.h"
#include "totients.h"
#include "lehmer.h"
#include "lmo.h"
#include "legendre_phi.h"
#include "aks.h"
#include "constants.h"
#include "mulmod.h"
#include "entropy.h"
#include "csprng.h"
#include "random_prime.h"
#include "perfect_powers.h"
#include "prime_powers.h"
#include "ramanujan_primes.h"
#include "semi_primes.h"
#include "twin_primes.h"
#include "almost_primes.h"
#include "omega_primes.h"
#include "prime_counts.h"
#include "prime_sums.h"
#include "powerfree.h"
#include "powerful.h"
#include "lucky_numbers.h"
#include "rootmod.h"
#include "real.h"
#include "ds_iset.h"

#ifdef FACTORING_HARNESSES
#include <sys/time.h>
static double my_difftime (struct timeval * start, struct timeval * end) {
  double secs, usecs;
  if (start->tv_sec == end->tv_sec) {
    secs = 0;
    usecs = end->tv_usec - start->tv_usec;
  } else {
    usecs = 1000000 - start->tv_usec;
    secs = end->tv_sec - (start->tv_sec + 1);
    usecs += end->tv_usec;
    if (usecs >= 1000000) {
      usecs -= 1000000;
      secs += 1;
    }
  }
  return secs + usecs / 1000000.;
}
#endif

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

#if PERL_REVISION <= 5 && PERL_VERSION <= 6 && BITS_PER_WORD == 64
 /* Workaround perl 5.6 UVs and bigints */
 #define my_svuv(sv)  PSTRTOULL(SvPV_nolen(sv), NULL, 10)
 #define my_sviv(sv)  PSTRTOLL(SvPV_nolen(sv), NULL, 10)
#elif PERL_REVISION <= 5 && PERL_VERSION < 14 && BITS_PER_WORD == 64
 /* Workaround RT 49569 in Math::BigInt::FastCalc (pre 5.14.0) */
 /* TODO: Math::BigInt::Pari has the same problem with negs pre-5.18.0 */
 #define my_svuv(sv) ( (!SvROK(sv)) ? SvUV(sv) : PSTRTOULL(SvPV_nolen(sv),NULL,10) )
 #define my_sviv(sv) ( (!SvROK(sv)) ? SvIV(sv) : PSTRTOLL(SvPV_nolen(sv),NULL,10) )
#else
 #define my_svuv(sv) SvUV(sv)
 #define my_sviv(sv) SvIV(sv)
#endif

#if PERL_REVISION >=5 && PERL_VERSION >= 9 && PERL_SUBVERSION >= 4
  #define SVf_MAGTEST  SVf_ROK
#else
  #define SVf_MAGTEST  SVf_AMAGIC
#endif

#define SVNUMTEST(n) \
  ((SvFLAGS(n) & (SVf_IOK | SVf_MAGTEST | SVs_GMG )) == SVf_IOK)

/* multicall compatibility stuff */
#if (PERL_REVISION <= 5 && PERL_VERSION < 7) || !defined(dMULTICALL)
# define USE_MULTICALL 0   /* Too much trouble to work around it */
#else
# define USE_MULTICALL 1
#endif
#if PERL_VERSION < 13 || (PERL_VERSION == 13 && PERL_SUBVERSION < 9)
#  define FIX_MULTICALL_REFCOUNT \
      if (CvDEPTH(multicall_cv) > 1) SvREFCNT_inc(multicall_cv);
#else
#  define FIX_MULTICALL_REFCOUNT
#endif

#ifndef CvISXSUB
#  define CvISXSUB(cv) CvXSUB(cv)
#endif

/* Not right, but close */
#if !defined cxinc && ( (PERL_VERSION == 8 && PERL_SUBVERSION >= 1) || (PERL_VERSION == 10 && PERL_SUBVERSION <= 1) )
# define cxinc() Perl_cxinc(aTHX)
#endif

#if PERL_VERSION < 17 || (PERL_VERSION == 17 && PERL_SUBVERSION < 7)
#  define SvREFCNT_dec_NN(sv)    SvREFCNT_dec(sv)
#endif

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

#define MY_CXT_KEY "Math::Prime::Util::API_guts"
#define CINTS 100
typedef struct {
  HV* MPUroot;
  HV* MPUGMP;
  HV* MPUPP;
  SV* const_int[CINTS+1];   /* -1, 0, 1, ..., 99 */
  void* randcxt;            /* per-thread csprng context */
  uint16_t forcount;
  char     forexit;
} my_cxt_t;

START_MY_CXT

static int _is_sv_bigint(pTHX_ SV* n)
{
  if (sv_isobject(n)) {
    const char *hvname = HvNAME_get(SvSTASH(SvRV(n)));
    if (hvname != 0) {
      if (strEQ(hvname, "Math::BigInt") || strEQ(hvname, "Math::BigFloat") ||
          strEQ(hvname, "Math::GMPz")   || strEQ(hvname, "Math::GMP") ||
          strEQ(hvname, "Math::GMPq")   || strEQ(hvname, "Math::AnyNum") ||
          strEQ(hvname, "Math::Pari")   || strEQ(hvname, "Math::BigInt::Lite"))
        return 1;
    }
  }
  return 0;
}

/* Is this a pedantically valid integer?
 * Croaks if undefined or invalid.
 * Returns 0 if it is an object or a string too large for a UV.
 * Returns 1 if it is good to process by XS.
 */
static int _validate_int(pTHX_ SV* n, int negok)
{
  const char* mustbe = (negok) ? "must be an integer" : "must be a non-negative integer";
  const char* maxstr;
  char* ptr;
  STRLEN i, len, maxlen;
  int ret, isbignum = 0, isneg = 0;

  /* TODO: magic, grok_number, etc. */
  if (SVNUMTEST(n)) { /* If defined as number, use it */
    if (SvIsUV(n) || SvIVX(n) >= 0)  return 1; /* The normal case */
    if (negok)  return -1;
    else croak("Parameter '%" SVf "' %s", n, mustbe);
  }
  if (sv_isobject(n)) {
    isbignum = _is_sv_bigint(aTHX_ n);
    if (!isbignum) return 0;
  }
  /* Without being very careful, don't process magic variables here */
  if (SvGAMAGIC(n) && !isbignum) return 0;
  if (!SvOK(n))  croak("Parameter must be defined");
  ptr = SvPV_nomg(n, len);             /* Includes stringifying bigints */
  if (len == 0 || ptr == 0)  croak("Parameter %s", mustbe);
  if (ptr[0] == '-' && negok) {
    isneg = 1; ptr++; len--;           /* Read negative sign */
  } else if (ptr[0] == '+') {
    ptr++; len--;                      /* Allow a single plus sign */
  }
  /* Empty string or non-numeric */
  if (len == 0 || !isDIGIT(ptr[0])) croak("Parameter '%" SVf "' %s", n, mustbe);
  /* Leading zeros and if left with only zero */
  while (len > 0 && *ptr == '0')       /* Strip all leading zeros */
    { ptr++; len--; }
  if (len == 0)                        /* 0 or -0 */
    return 1;
  /* We're going to look more carefully at the string to ensure it's a number */
  if (isneg) { ret = -1;  maxlen = ivmax_maxlen;  maxstr = ivmin_str; }
  else       { ret =  1;  maxlen = uvmax_maxlen;  maxstr = uvmax_str; }
  if (len > maxlen)
    return 0;                          /* Huge number, don't even look at it */
  for (i = 0; i < len; i++)            /* Ensure all characters are digits */
    if (!isDIGIT(ptr[i]))
      croak("Parameter '%" SVf "' %s", n, mustbe);
  if (len < maxlen)                    /* Valid small integer */
    return ret;
  for (i = 0; i < maxlen; i++) {       /* Check if in range */
    if (ptr[i] < maxstr[i]) return ret;
    if (ptr[i] > maxstr[i]) return 0;
  }
  return ret;                          /* value = UV_MAX/UV_MIN.  That's ok */
}

#define IFLAG_ANY      0x00000000U
#define IFLAG_POS      0x00000001U  /* Must be non-negative */
#define IFLAG_NONZERO  0x00000002U  /* Must not be zero */
#define IFLAG_ABS      0x00000004U  /* Absolute value returned */
#define IFLAG_IV       0x00000008U  /* Value returned as IV */

static int _validate_and_set(UV* val, pTHX_ SV* svn, uint32_t mask) {
  int status;

  if (svn == 0) croak("Parameter must be defined");
  /* Streamline the typical path of input being a native integer. */
  if (SVNUMTEST(svn)) {
    IV n = SvIVX(svn);
    if (n >= 0) {
      if (n == 0 && (mask & IFLAG_NONZERO))
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
    if (mask & IFLAG_ABS) { *val = (UV)(-n); return 1; }
    if (mask & IFLAG_POS) croak("Parameter '%" SVf "' must be a non-negative integer", svn);
    *val = n;
    return -1;
  }

  status = _validate_int( aTHX_ svn, !(mask & IFLAG_POS) );
  if (status == 1) {
    UV n = my_svuv(svn);
    if (n == 0 && (mask & IFLAG_NONZERO))
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

static UV neg_iv(UV n) {
  if ((IV)n == IV_MIN)  return (UV_MAX >> 1) + 1;
  else                  return (UV) (-(IV)n);
}

/* Given 'a' and astatus (-1 means 'a' is an IV), properly mod with n */
static void _mod_with(UV *a, int astatus, UV n) {
  if (n == 0) return;
  if (astatus != -1) {
    *a %= n;
  } else {
    UV r = neg_iv(*a) % n;
    *a = (r == 0) ? 0 : n-r;
  }
}

#define VCALL_ROOT 0x0
#define VCALL_PP   0x1
#define VCALL_GMP  0x2
/* Call a Perl sub to handle work for us. */
static int _vcallsubn(pTHX_ I32 flags, I32 stashflags, const char* name, int nargs, int minversion)
{
    GV* gv = NULL;
    dMY_CXT;
    Size_t namelen = strlen(name);
    /* If given a GMP function, and GMP enabled, and function exists, use it. */
    int use_gmp = stashflags & VCALL_GMP && _XS_get_callgmp() && _XS_get_callgmp() >= minversion;
    assert(!(stashflags & ~(VCALL_PP|VCALL_GMP)));
    if (use_gmp && hv_exists(MY_CXT.MPUGMP,name,namelen)) {
      GV ** gvp = (GV**)hv_fetch(MY_CXT.MPUGMP,name,namelen,0);
      if (gvp) gv = *gvp;
    }
    if (!gv && (stashflags & VCALL_PP))
      perl_require_pv("Math/Prime/Util/PP.pm");
    if (!gv) {
      GV ** gvp = (GV**)hv_fetch(stashflags & VCALL_PP? MY_CXT.MPUPP : MY_CXT.MPUroot, name,namelen,0);
      if (gvp) gv = *gvp;
    }
    /* use PL_stack_sp in PUSHMARK macro directly it will be read after
      the possible mark stack extend */
    PUSHMARK(PL_stack_sp-nargs);
    /* no PUTBACK bc we didn't move global SP */
    return call_sv((SV*)gv, flags);
}
#define _vcallsub(func) (void)_vcallsubn(aTHX_ G_SCALAR, VCALL_ROOT, func, items,0)
#define _vcallsub_with_gmp(ver,func) (void)_vcallsubn(aTHX_ G_SCALAR, VCALL_GMP|VCALL_PP, func, items,(int)(100*(ver)))
#define _vcallsub_with_pp(func) (void)_vcallsubn(aTHX_ G_SCALAR, VCALL_PP, func, items,0)
#define CALLPPSUB(func) (void)_vcallsubn(aTHX_ GIMME_V,VCALL_PP,func,items,0)
/* #define _vcallsub_with_gmpobj(ver,func) (void)_vcallsubn(aTHX_ G_SCALAR, (PERL_REVISION >= 5 && PERL_VERSION > 8) ? VCALL_GMP|VCALL_PP : VCALL_PP, func, items,(int)(100*(ver))) */
#define _vcallsub_with_gmpobj(ver,func) _vcallsub_with_gmp(ver,func)

#if 0
static int _vcallgmpsubn(pTHX_ I32 flags, const char* name, int nargs, int minversion)
{
  Size_t namelen = strlen(name);
  int gmpver = _XS_get_callgmp();
  dMY_CXT;
  if (gmpver && gmpver >= minversion && hv_exists(MY_CXT.MPUGMP,name,namelen)) {
    GV ** gvp = (GV**)hv_fetch(MY_CXT.MPUGMP,name,namelen,0);
    if (gvp) {
      GV* gv = *gvp;
      PUSHMARK(PL_stack_sp-nargs);
      return call_sv((SV*)gv, flags);
    }
  }
  return 0;
}
#endif

/* In my testing, this constant return works fine with threads, but to be
 * correct (see perlxs) one has to make a context, store separate copies in
 * each one, then retrieve them from a struct using a hash index.  This
 * defeats the purpose if only done once. */
#define RETURN_NPARITY(ret) \
  do { int r_ = ret; \
       dMY_CXT; \
       if (r_ >= -1 && r_<CINTS) { ST(0) = MY_CXT.const_int[r_+1]; XSRETURN(1); } \
       else                      { XSRETURN_IV(r_);                      } \
  } while (0)
#define PUSH_NPARITY(ret) \
  do { int r_ = ret; \
       if (r_ >= -1 && r_<CINTS) { PUSHs( MY_CXT.const_int[r_+1] );       } \
       else                      { PUSHs(sv_2mortal(newSViv(r_))); } \
  } while (0)

static void objectify_result(pTHX_ SV* input, SV* output) {
  /* Leave unchanged: undef, objects, small integers */
  if (!SvOK(output) || sv_isobject(output) || SVNUMTEST(output))
    return;
  /* If they didn't give us a bigint, then try to be smart */
  if (!input || !sv_isobject(input)) {
    (void)_vcallsubn(aTHX_ G_SCALAR, VCALL_ROOT, "_to_bigint_if_needed", 1, 0);
  } else {
    const char *iname = HvNAME_get(SvSTASH(SvRV(input)));
    if (strEQ(iname, "Math::BigInt")) {
      (void)_vcallsubn(aTHX_ G_SCALAR, VCALL_ROOT, "_to_bigint", 1, 0);
    } else if (strEQ(iname, "Math::GMPz")) {
      (void)_vcallsubn(aTHX_ G_SCALAR, VCALL_ROOT, "_to_gmpz", 1, 0);
    } else if (strEQ(iname, "Math::GMP")) {
      (void)_vcallsubn(aTHX_ G_SCALAR, VCALL_ROOT, "_to_gmp", 1, 0);
    } else { /* Return it as: ref(input)->new(result) */
      dSP;  (void)POPs;  ENTER;  PUSHMARK(SP);
      XPUSHs(sv_2mortal(newSVpv(iname, 0)));  XPUSHs(output);
      PUTBACK;  call_method("new", G_SCALAR);  LEAVE;
    }
  }
}

static SV* call_sv_to_func(pTHX_ SV* r, const char* name) {
  dSP;  ENTER;  PUSHMARK(SP);
  XPUSHs(r);
  PUTBACK;
  call_pv(name, G_SCALAR);
  SPAGAIN;
  r = POPs;
  PUTBACK; LEAVE;
  return r;
}
static SV* sv_to_bigint(pTHX_ SV* r) {
  return call_sv_to_func(aTHX_ r, "Math::Prime::Util::_to_bigint");
}
static SV* sv_to_bigint_abs(pTHX_ SV* r) {
  return call_sv_to_func(aTHX_ r, "Math::Prime::Util::_to_bigint_abs");
}
static SV* sv_to_bigint_nonneg(pTHX_ SV* r) {
  return call_sv_to_func(aTHX_ r, "Math::Prime::Util::_to_bigint_nonneg");
}

#define NEWSVINT(sign,v) (((sign) > 0) ? newSVuv(v) : newSViv(v))

#define RETURN_128(hi,lo) \
  do { char str_[40]; \
       int slen_ = to_string_128(str_, hi, lo); \
       ST(0) = sv_to_bigint( aTHX_ sv_2mortal(newSVpv(str_,slen_)) ); \
       XSRETURN(1); } while(0)

#define CREATE_AV(av) \
  do { \
       av = newAV(); \
       { \
         SV * retsv = sv_2mortal(newRV_noinc( (SV*) av )); \
         PUSHs(retsv); \
         PUTBACK; \
         SP = NULL; /* never use SP again, poison */ \
       } \
  } while(0)

#define RETURN_LIST_VALS(in_alen,arr,sign)   /* Return array values */ \
  { \
    unsigned long k_, alen_ = in_alen; \
    EXTEND(SP,(long)alen_); \
    for (k_ = 0; k_ < alen_; k_++) \
      ST(k_) = sv_2mortal(NEWSVINT(sign,arr[k_])); \
    Safefree(arr); \
    XSRETURN(alen_); \
  }

#define RETURN_LIST_REF(in_alen,arr,sign)   /* Return array values as ref */ \
  { \
    AV* av_; \
    unsigned long k_, alen_ = in_alen; \
    CREATE_AV(av_); \
    av_extend(av_, (SSize_t)alen_); \
    for (k_ = 0; k_ < alen_; k_++) \
      av_push(av_, NEWSVINT(sign,arr[k_])); \
    Safefree(arr); \
    XSRETURN(1); \
  }

static int arrayref_to_digit_array(pTHX_ UV** ret, AV* av, int base)
{
  int len, i;
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


#define IARR_TYPE_ANY 0x00
#define IARR_TYPE_NEG 0x01
#define IARR_TYPE_POS 0x02
#define IARR_TYPE_BAD 0x03

/* BAD always bad, ANY with ANY/POS/NEG, POS and NEG only with ANY and self. */
#define CAN_COMBINE_IARR_TYPES(t1,t2)  ( ((t1) | (t2)) != IARR_TYPE_BAD )
/* Convert to 0/1/-1 status */
#define IARR_TYPE_TO_STATUS(t) \
  (((t) == IARR_TYPE_BAD) ? 0 : ((t) == IARR_TYPE_NEG) ? -1 : 1)
#define STATUS_TO_IARR_TYPE(s,n) \
  (((s) == 0) ? IARR_TYPE_BAD : ((s) == -1) ? IARR_TYPE_NEG : ((n) > (UV)IV_MAX) ? IARR_TYPE_POS : IARR_TYPE_ANY)

/* Compare using first argument non-zero to indicate UV, otherwise IV */
#define SIGNED_CMP_LE(pos,x,y) ((pos)  ?  (x <= y)  :  ((IV)x <= (IV)y))
#define SIGNED_CMP_LT(pos,x,y) ((pos)  ?  (x <  y)  :  ((IV)x <  (IV)y))
#define SIGNED_CMP_GT(pos,x,y) ((pos)  ?  (x >  y)  :  ((IV)x >  (IV)y))

/* Given values and a sign indicating IV or UV, returns -1 (<), 0 (eq), 1 (>) */
static int _sign_cmp(int xsign, UV x, int ysign, UV y) {
  if (xsign == ysign && x == y) return 0;
  /* Convert sign to -1 (neg), 0 (small pos), 1 (big pos) */
  if (x <= (UV)IV_MAX) xsign = 0;
  if (y <= (UV)IV_MAX) ysign = 0;
  if (xsign != ysign) return (xsign < ysign) ? -1 : 1;
  if (xsign == -1) return ((IV)x < (IV)y) ? -1 : 1;
  return (x < y) ? -1 : 1;
}

#define DEBUG_PRINT_ARRAY(name,av) \
  { Size_t j_; SV** arr_ = AvARRAY(av);  printf("%s: [",name);  for(j_=0; j_<av_count(av); j_++) printf("%lu ",SvUV(arr_[j_])); printf("]\n"); }

static int arrayref_to_int_array(pTHX_ unsigned long *retlen, UV** ret, int want_sort, SV* sva, const char* fstr)
{
  AV *av;
  int len, i, itype = IARR_TYPE_ANY;
  UV  *r;

  if ((!SvROK(sva)) || (SvTYPE(SvRV(sva)) != SVt_PVAV))
    croak("%s argument must be an array reference", fstr);
  av = (AV*) SvRV(sva);
  len = av_len(av);
  if (len < 0) {
    *retlen = 0;  *ret = 0;  return itype;
  }
  New(0, r, len+1, UV);
  for (i = 0; i <= len; i++) {
    SV **iv = av_fetch(av, i, 0);
    if (iv == 0) croak("Array entries must be defined");
    if (SVNUMTEST(*iv)) {
      IV n = SvIVX(*iv);
      if (n < 0) {
        if (SvIsUV(*iv))  itype |= IARR_TYPE_POS;
        else              itype |= IARR_TYPE_NEG;
        if (itype == IARR_TYPE_BAD) break;
      }
      r[i] = (UV)n;
    } else {
      UV n;
      int istatus = _validate_and_set(&n, aTHX_ *iv, IFLAG_ANY);
      if (istatus == -1) {
        itype |= IARR_TYPE_NEG;
      } else if (istatus == 1 && n > (UV)IV_MAX) {
        itype |= IARR_TYPE_POS;
      }
      if (istatus == 0 || itype == IARR_TYPE_BAD) break;
      r[i] = n;
    }
  }
  if (i <= len) {
    Safefree(r);
    *retlen = 0;  *ret = 0;  return IARR_TYPE_BAD;
  }
  *retlen = len+1;
  *ret = r;
  if (want_sort) {
    if (itype == IARR_TYPE_NEG) {
      for (i = 1; i <= len; i++)
        if ( (IV)r[i] <= (IV)r[i-1] )
          break;
    } else {
      for (i = 1; i <= len; i++)
        if (r[i] <= r[i-1])
          break;
    }
    if (i <= len)
      sort_dedup_uv_array(*ret, itype == IARR_TYPE_NEG, retlen);
  }
  return itype;
}

static int type_of_sumset(int typea, int typeb, UV amin, UV amax, UV bmin, UV bmax) {
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

#define SC_SIZE  257   /* Choose 131, 257, 521, 1031, 2053 */
typedef struct {
  char lostatus, histatus, *midstatus;
  UV loval, hival, *midval;
  int *midindex;
} set_data_t;
static set_data_t init_set_lookup_cache(pTHX_ AV *av) {
  set_data_t d;
  int len = av_count(av);
  if (len > SC_SIZE) len = SC_SIZE;
  d.lostatus = d.histatus = 0;
  Newz(0, d.midstatus, len, char);
  New(0, d.midval, len, UV);
  New(0, d.midindex, len, int);
  return d;
}
static void free_set_lookup_cache(set_data_t *d) {
  Safefree(d->midstatus);
  Safefree(d->midval);
  Safefree(d->midindex);
}
#define _SC_GET_VALUE(statvar, var, av, i) \
  SV** iv = av_fetch(av, i, 0); \
  statvar = _validate_and_set(&var, aTHX_ iv ? *iv : 0, IFLAG_ANY); \
  if (statvar == 0) return -1;

#define SC_SET_MID_VALUE(statvar, var, av, i, cache) \
  do { \
    unsigned int imod_ = (i) % SC_SIZE; \
    if (cache != 0 && cache->midstatus[imod_] != 0 && cache->midindex[imod_] == i) { \
      statvar = cache->midstatus[imod_];  var = cache->midval[imod_]; \
    } else { \
      _SC_GET_VALUE(statvar, var, av, i) \
      if (cache) { cache->midstatus[imod_] = statvar;  cache->midval[imod_] = var; cache->midindex[imod_] = i; } \
    } \
  } while (0)

static int _sc_set_lohi(pTHX_ AV* av, set_data_t *cache, int loindex, int hiindex, int *lostatus, int *histatus, UV *loval, UV *hival)
{
  if (cache != 0 && cache->lostatus != 0) {
    *lostatus = cache->lostatus;  *loval = cache->loval;
  } else {
    _SC_GET_VALUE(*lostatus, *loval, av, loindex);
    if (cache) { cache->lostatus = *lostatus;  cache->loval = *loval; }
  }
  if (cache != 0 && cache->histatus != 0) {
    *histatus = cache->histatus;  *hival = cache->hival;
  } else {
    _SC_GET_VALUE(*histatus, *hival, av, hiindex);
    if (cache) { cache->histatus = *histatus;  cache->hival = *hival; }
  }
  return 1;
}

/* Find index in a set (array ref of sorted unique integers)
 *    -1 bigint
 *    0  already in set
 *    n  should be in n-th position (1 means should be first element)
 */
static int index_in_set(pTHX_ AV* av, set_data_t *cache, int sign, UV val)
{
  int last_index, lo, hi, lostatus, histatus, midstatus, cmp;
  UV  rlo, rhi, rmid;

  if (sign != 1 && sign != -1)
    return -1;
  last_index = av_len(av);
  if (last_index < 0)
    return 1;

  lo = 0;
  hi = last_index;
  if (_sc_set_lohi(aTHX_ av, cache, lo, hi, &lostatus, &histatus, &rlo, &rhi) < 0)
    return -1;

  cmp = _sign_cmp(sign, val, lostatus, rlo);
  if (cmp <= 0) return (cmp == 0) ? 0 : lo+1;
  /* val > rlo */
  cmp = _sign_cmp(sign, val, histatus, rhi);
  if (cmp >= 0) return (cmp == 0) ? 0 : hi+2;
  /* val < rhi */

  while (hi-lo > 1) {
    int mid = lo + ((hi-lo) >> 1);
    SC_SET_MID_VALUE(midstatus, rmid, av, mid, cache);
    cmp = _sign_cmp(midstatus, rmid, sign, val);
    if (cmp == 0) return 0;
    if (cmp < 0) { lo = mid; rlo = rmid; lostatus = midstatus; }
    else         { hi = mid; rhi = rmid; histatus = midstatus; }
  }
  if (sign == histatus && rhi == val) return 0;
  if (_sign_cmp(sign,val, histatus, rhi) > 0) croak("internal index error");
  return hi+1;
}

/* See if an element is in a set (array ref of sorted unique integers) */
/* -1 = bigint, 0 = not found, 1 = found */
static int is_in_set(pTHX_ AV* av, set_data_t *cache, int sign, UV val)
{
  int index = index_in_set(aTHX_ av, cache, sign, val);
  return (index < 0) ? index : !index;
}

#if 0
/* Take two SV pointers to array references, slurp them all in.
 * Returns -1 if cannot answer, 0 if not a subset, 1 if a subset. */
static int is_subset_read_all(pTHX_ SV* svset, SV* svsub)
{
  int atype, btype, bstatus, res;
  UV *ra, *rb;
  unsigned long alen, blen, ia, ib;

  atype = arrayref_to_int_array(aTHX_ &alen, &ra, 1, svset, "setcontains arg 1");
  if (atype != IARR_TYPE_BAD) {
    if (SvROK(svsub) && SvTYPE(SvRV(svsub)) == SVt_PVAV) {
      btype = arrayref_to_int_array(aTHX_ &blen, &rb, 1, svsub, "setcontains arg 1");
    } else {   /* We got a single integer */
      blen = 1;
      New(0, rb, blen, UV);
      bstatus = _validate_and_set(rb+0, aTHX_ svsub, IFLAG_ANY);
      btype = STATUS_TO_IARR_TYPE(bstatus, rb[0]);
    }
  }
  if (atype == IARR_TYPE_BAD || btype == IARR_TYPE_BAD) {
    Safefree(ra);
    Safefree(rb);
    return -1;
  }
  res = 1;
  /* Check for empty sets */
  if (alen == 0 || blen == 0)
    { Safefree(ra);  Safefree(rb);  return (blen == 0); }
  /* Check if impossible via types */
  if ( (btype == IARR_TYPE_NEG && atype != IARR_TYPE_NEG) ||
       (btype == IARR_TYPE_POS && atype != IARR_TYPE_POS) ||
       (blen > alen) )
    res = 0;
  /* Check to see if the subset is outside the range */
  if ( (atype != IARR_TYPE_NEG &&
        (rb[0] < ra[0] || rb[blen-1] > ra[alen-1])) ||
       (atype == IARR_TYPE_NEG &&
        ((IV)rb[0] < (IV)ra[0] || (IV)rb[blen-1] > (IV)ra[alen-1])) )
    res = 0;
  for (ia = 0, ib = 0;  ia < alen && ib < blen && res; ia++) {
    if (ra[ia] == rb[ib])
      ib++;
    else if (SIGNED_CMP_GT(btype != IARR_TYPE_NEG, ra[ia], rb[ib]))
      res = 0;
  }
  Safefree(ra);
  Safefree(rb);
  return res;
}
#endif

static int _compare_array_refs(pTHX_ SV* a, SV* b)
{
  AV *ava, *avb;
  int i, alen, blen;
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

static void csprng_init_seed(void* ctx) {
  unsigned char* data;
  New(0, data, 64, unsigned char);
  get_entropy_bytes(64, data);
  csprng_seed(ctx, 64, data);
  Safefree(data);
}

static void _comb_init(UV* cm, UV k, int derangement) {
  UV i;
  cm[0] = UV_MAX;
  for (i = 0; i < k; i++)
    cm[i] = k-i;
  if (derangement && k >= 2) {   /* Make derangements start deranged */
    for (i = 0; i < k; i++)
      cm[k-i-1] = (i&1) ? i : i+2;
    if (k & 1) {
      cm[0] = k-2;
      cm[1] = k;
    }
  }
}

static int _comb_iterate(UV* cm, UV k, UV n, int ix) {
  UV i, j, m;
  if (ix == 0) {
    if (cm[0]++ < n)  return 0;                /* Increment last value */
    for (i = 1; i < k && cm[i] >= n-i; i++) ;  /* Find next index to incr */
    if (i >= k)  return 1;                     /* Done! */
    cm[i]++;                                   /* Increment this one */
    while (i-- > 0)  cm[i] = cm[i+1] + 1;      /* Set the rest */
  } else if (ix == 1) {
    for (j = 1; j < k && cm[j] > cm[j-1]; j++) ;  /* Find last decrease */
    if (j >= k) return 1;                         /* Done! */
    for (m = 0; cm[j] > cm[m]; m++) ;             /* Find next greater */
    { UV t = cm[j];  cm[j] = cm[m];  cm[m] = t; } /* Swap */
    for (i = j-1, m = 0;  m < i;  i--, m++)       /* Reverse the end */
      { UV t = cm[i];  cm[i] = cm[m];  cm[m] = t; }
  } else {
    REDERANGE:
    for (j = 1; j < k && cm[j] > cm[j-1]; j++) ;  /* Find last decrease */
    if (j >= k) return 1;                         /* Done! */
    for (m = 0; cm[j] > cm[m]; m++) ;             /* Find next greater */
    { UV t = cm[j];  cm[j] = cm[m];  cm[m] = t; } /* Swap */
    if (cm[j] == k-j) goto REDERANGE;             /* Skip? */
    for (i = j-1, m = 0;  m < i;  i--, m++)       /* Reverse the end */
      { UV t = cm[i];  cm[i] = cm[m];  cm[m] = t; }
    for (i = 0; i < k; i++)                       /* Check deranged */
      if (cm[k-i-1]-1 == i)
        break;
    if (i != k) goto REDERANGE;
  }
  return 0;
}



MODULE = Math::Prime::Util	PACKAGE = Math::Prime::Util

PROTOTYPES: ENABLE

BOOT:
{
    int i;
    SV * sv = newSViv(BITS_PER_WORD);
    HV * stash = gv_stashpv("Math::Prime::Util", TRUE);
    newCONSTSUB(stash, "_XS_prime_maxbits", sv);

    {
      MY_CXT_INIT;
      MY_CXT.MPUroot = stash;
      MY_CXT.MPUGMP = gv_stashpv("Math::Prime::Util::GMP", TRUE);
      MY_CXT.MPUPP = gv_stashpv("Math::Prime::Util::PP", TRUE);
      for (i = 0; i <= CINTS; i++) {
        MY_CXT.const_int[i] = newSViv(i-1);
        SvREADONLY_on(MY_CXT.const_int[i]);
      }
      New(0, MY_CXT.randcxt, csprng_context_size(), char);
      csprng_init_seed(MY_CXT.randcxt);
      MY_CXT.forcount = 0;
      MY_CXT.forexit = 0;
   }
}

#if defined(USE_ITHREADS) && defined(MY_CXT_KEY)

void
CLONE(...)
PREINIT:
  int i;
PPCODE:
  {
    MY_CXT_CLONE; /* possible declaration */
    MY_CXT.MPUroot = gv_stashpv("Math::Prime::Util", TRUE);
    MY_CXT.MPUGMP = gv_stashpv("Math::Prime::Util::GMP", TRUE);
    MY_CXT.MPUPP = gv_stashpv("Math::Prime::Util::PP", TRUE);
    /* These should be shared between threads, but that's dodgy. */
    for (i = 0; i <= CINTS; i++) {
      MY_CXT.const_int[i] = newSViv(i-1);
      SvREADONLY_on(MY_CXT.const_int[i]);
    }
    /* Make a new CSPRNG context for this thread */
    New(0, MY_CXT.randcxt, csprng_context_size(), char);
    csprng_init_seed(MY_CXT.randcxt);
    /* NOTE:  There is no thread destroy, so these never get freed... */
    MY_CXT.forcount = 0;
    MY_CXT.forexit = 0;
  }
  return; /* skip implicit PUTBACK, returning @_ to caller, more efficient*/

#endif

void
END(...)
PREINIT:
  dMY_CXT;
  int i;
PPCODE:
  _prime_memfreeall();
  MY_CXT.MPUroot = NULL;
  MY_CXT.MPUGMP = NULL;
  MY_CXT.MPUPP = NULL;
  for (i = 0; i <= CINTS; i++) {
    SV * const sv = MY_CXT.const_int[i];
    MY_CXT.const_int[i] = NULL;
    SvREFCNT_dec_NN(sv);
  } /* stashes are owned by stash tree, no refcount on them in MY_CXT */
  Safefree(MY_CXT.randcxt); MY_CXT.randcxt = 0;
  return; /* skip implicit PUTBACK, returning @_ to caller, more efficient*/


void csrand(IN SV* seed = 0)
  PREINIT:
    unsigned char* data;
    STRLEN size;
    dMY_CXT;
  PPCODE:
    if (items == 0) {
      csprng_init_seed(MY_CXT.randcxt);
    } else if (_XS_get_secure()) {
      croak("secure option set, manual seeding disabled");
    } else {
      data = (unsigned char*) SvPV(seed, size);
      csprng_seed(MY_CXT.randcxt, size, data);
    }
    if (_XS_get_callgmp() >= 42) _vcallsub("_csrand_p");
    return;

UV srand(IN UV seedval = 0)
  PREINIT:
    dMY_CXT;
  CODE:
    if (_XS_get_secure())
      croak("secure option set, manual seeding disabled");
    if (items == 0)
      get_entropy_bytes(sizeof(UV), (unsigned char*) &seedval);
    csprng_srand(MY_CXT.randcxt, seedval);
    if (_XS_get_callgmp() >= 42) _vcallsub("_srand_p");
    RETVAL = seedval;
  OUTPUT:
    RETVAL

UV irand()
  ALIAS:
    irand64 = 1
  PREINIT:
    dMY_CXT;
  CODE:
    if (ix == 0)
      RETVAL = irand32(MY_CXT.randcxt);
    else
#if BITS_PER_WORD >= 64
      RETVAL = irand64(MY_CXT.randcxt);
#else /* TODO: should irand64 on 32-bit perl (1) croak, (2) return 32-bits */
      RETVAL = irand32(MY_CXT.randcxt);
#endif
  OUTPUT:
    RETVAL

NV drand(NV m = 0.0)
  ALIAS:
    rand = 1
  PREINIT:
    dMY_CXT;
  CODE:
    PERL_UNUSED_VAR(ix);
    RETVAL = drand64(MY_CXT.randcxt);
    if (m != 0) RETVAL *= m;
  OUTPUT:
    RETVAL

SV* random_bytes(IN UV n)
  PREINIT:
    char* sptr;
    dMY_CXT;
  CODE:
    RETVAL = newSV(n == 0 ? 1 : n);
    SvPOK_only(RETVAL);
    SvCUR_set(RETVAL, n);
    sptr = SvPVX(RETVAL);
    csprng_rand_bytes(MY_CXT.randcxt, n, (unsigned char*)sptr);
    sptr[n] = '\0';
  OUTPUT:
    RETVAL

SV* entropy_bytes(IN UV n)
  PREINIT:
    char* sptr;
  CODE:
    RETVAL = newSV(n == 0 ? 1 : n);
    SvPOK_only(RETVAL);
    SvCUR_set(RETVAL, n);
    sptr = SvPVX(RETVAL);
    get_entropy_bytes(n, (unsigned char*)sptr);
    sptr[n] = '\0';
  OUTPUT:
    RETVAL

UV _is_csprng_well_seeded()
  ALIAS:
    _XS_get_verbose = 1
    _XS_get_callgmp = 2
    _XS_get_secure = 3
    _XS_set_secure = 4
    _get_forexit = 5
    _start_for_loop = 6
    _get_prime_cache_size = 7
  CODE:
    switch (ix) {
      case 0:  { dMY_CXT; RETVAL = is_csprng_well_seeded(MY_CXT.randcxt); } break;
      case 1:  RETVAL = _XS_get_verbose(); break;
      case 2:  RETVAL = _XS_get_callgmp(); break;
      case 3:  RETVAL = _XS_get_secure(); break;
      case 4:  _XS_set_secure(); RETVAL = 1; break;
      case 5:  { dMY_CXT; RETVAL = MY_CXT.forexit; } break;
      case 6:  { dMY_CXT; MY_CXT.forcount++; RETVAL = MY_CXT.forexit; MY_CXT.forexit = 0; } break;
      case 7:
      default: RETVAL = get_prime_cache(0,0); break;
    }
  OUTPUT:
    RETVAL

bool _validate_integer(SV* svn)
  ALIAS:
    _validate_integer_nonneg = 1
    _validate_integer_positive = 2
    _validate_integer_abs = 3
  PREINIT:
    uint32_t mask;
    int status;
    UV n;
  CODE:
    /* Flag:  0 neg ok,  1 neg err,  2 zero or neg err,  3 abs */
    switch (ix) {
      case 0: mask = IFLAG_ANY; break;
      case 1: mask = IFLAG_POS; break;
      case 2: mask = IFLAG_POS | IFLAG_NONZERO; break;
      case 3: mask = IFLAG_ABS; break;
      default: croak("_validate_integer unknown flag value");
    }
    status = _validate_and_set(&n, aTHX_ svn, mask);
    if (status != 0) {
#if PERL_REVISION <= 5 && PERL_VERSION < 8 && BITS_PER_WORD == 64
      sv_setpviv(svn, n);
#else
      if (status == 1)  sv_setuv(svn, n);
      else              sv_setiv(svn, n);
#endif
    } else {  /* Status 0 = bigint */
      if (mask & IFLAG_ABS) {
        /* TODO: if given a positive bigint, no need for this */
        sv_setsv(svn, sv_to_bigint_abs(aTHX_ svn));
      } else if (mask & IFLAG_POS) {
        if (!_is_sv_bigint(aTHX_ svn))
          sv_setsv(svn, sv_to_bigint_nonneg(aTHX_ svn));
      } else {
        if (!_is_sv_bigint(aTHX_ svn))
          sv_setsv(svn, sv_to_bigint(aTHX_ svn));
      }
    }
    RETVAL = TRUE;
  OUTPUT:
    RETVAL

void prime_memfree()
  PREINIT:
    dMY_CXT;
  PPCODE:
    prime_memfree();
    /* (void) _vcallgmpsubn(aTHX_ G_VOID|G_DISCARD, "_GMP_memfree", 0, 49); */
    if (MY_CXT.MPUPP != NULL) _vcallsub_with_pp("prime_memfree");
    return;

void
prime_precalc(IN UV n)
  ALIAS:
    _XS_set_verbose = 1
    _XS_set_callgmp = 2
    _end_for_loop = 3
  PPCODE:
    PUTBACK; /* SP is never used again, the 3 next func calls are tailcall
    friendly since this XSUB has nothing to do after the 3 calls return */
    switch (ix) {
      case 0:  prime_precalc(n);    break;
      case 1:  _XS_set_verbose(n);  break;
      case 2:  _XS_set_callgmp(n);  break;
      case 3:
      default: { dMY_CXT; MY_CXT.forcount--; MY_CXT.forexit = n; } break;
    }
    return; /* skip implicit PUTBACK */


void prime_count(IN SV* svlo, IN SV* svhi = 0)
  ALIAS:
    semiprime_count = 1
    twin_prime_count = 2
    ramanujan_prime_count = 3
    perfect_power_count = 4
    prime_power_count = 5
    lucky_count = 6
  PREINIT:
    UV lo = 0, hi, count = 0;
  PPCODE:
    if ((items == 1 && _validate_and_set(&hi, aTHX_ svlo, IFLAG_POS)) ||
        (items == 2 && _validate_and_set(&lo, aTHX_ svlo, IFLAG_POS) && _validate_and_set(&hi, aTHX_ svhi, IFLAG_POS))) {
      if (lo <= hi) {
        switch (ix) {
          case 0:  count = prime_count_range(lo, hi);           break;
          case 1:  count = semiprime_count_range(lo, hi);       break;
          case 2:  count = twin_prime_count_range(lo, hi);      break;
          case 3:  count = ramanujan_prime_count_range(lo, hi); break;
          case 4:  count = perfect_power_count_range(lo, hi);   break;
          case 5:  count = prime_power_count_range(lo, hi);     break;
          case 6:  count = lucky_count_range(lo, hi);     break;
        }
      }
      XSRETURN_UV(count);
    }
    switch (ix) {
      case 0: _vcallsubn(aTHX_ GIMME_V, VCALL_ROOT, "_generic_prime_count", items, 0); break;
      case 1: _vcallsub_with_pp("semiprime_count");  break;
      case 2: _vcallsub_with_pp("twin_prime_count");  break;
      case 3: _vcallsub_with_pp("ramanujan_prime_count");  break;
      case 4: _vcallsub_with_gmpobj(0.53,"perfect_power_count"); break;
      case 5: _vcallsub_with_gmpobj(0.53,"prime_power_count"); break;
      case 6: _vcallsub_with_pp("lucky_count"); break;
    }
    return;


void prime_count_upper(IN SV* svn)
  ALIAS:
    prime_count_lower = 1
    prime_count_approx = 2
    prime_power_count_upper = 3
    prime_power_count_lower = 4
    prime_power_count_approx = 5
    perfect_power_count_upper = 6
    perfect_power_count_lower = 7
    perfect_power_count_approx = 8
    ramanujan_prime_count_upper = 9
    ramanujan_prime_count_lower = 10
    ramanujan_prime_count_approx = 11
    twin_prime_count_approx = 12
    semiprime_count_approx = 13
    lucky_count_upper = 14
    lucky_count_lower = 15
    lucky_count_approx = 16
  PREINIT:
    UV n, ret;
  PPCODE:
    if (_validate_and_set(&n, aTHX_ svn, IFLAG_POS)) {
      switch (ix) {
        case  0: ret = prime_count_upper(n); break;
        case  1: ret = prime_count_lower(n); break;
        case  2: ret = prime_count_approx(n); break;
        case  3: ret = prime_power_count_upper(n); break;
        case  4: ret = prime_power_count_lower(n); break;
        case  5: ret = prime_power_count_approx(n); break;
        case  6: ret = perfect_power_count_upper(n); break;
        case  7: ret = perfect_power_count_lower(n); break;
        case  8: ret = perfect_power_count_approx(n); break;
        case  9: ret = ramanujan_prime_count_upper(n); break;
        case 10: ret = ramanujan_prime_count_lower(n); break;
        case 11: ret = ramanujan_prime_count_approx(n); break;
        case 12: ret = twin_prime_count_approx(n); break;
        case 13: ret = semiprime_count_approx(n); break;
        case 14: ret = lucky_count_upper(n); break;
        case 15: ret = lucky_count_lower(n); break;
        case 16:
        default: ret = lucky_count_approx(n); break;
      }
      XSRETURN_UV(ret);
    }
    switch (ix) {
      case  0: _vcallsub_with_pp("prime_count_upper");   break;
      case  1: _vcallsub_with_pp("prime_count_lower");   break;
      case  2: _vcallsub_with_pp("prime_count_approx");   break;
      case  3: _vcallsub_with_pp("prime_power_count_upper");   break;
      case  4: _vcallsub_with_pp("prime_power_count_lower");   break;
      case  5: _vcallsub_with_pp("prime_power_count_approx");   break;
      case  6: _vcallsub_with_pp("perfect_power_count_upper");   break;
      case  7: _vcallsub_with_pp("perfect_power_count_lower");   break;
      case  8: _vcallsub_with_pp("perfect_power_count_approx");   break;
      case  9: _vcallsub_with_pp("ramanujan_prime_count_upper");  break;
      case 10: _vcallsub_with_pp("ramanujan_prime_count_lower");  break;
      case 11: _vcallsub_with_pp("ramanujan_prime_count_approx");  break;
      case 12: _vcallsub_with_pp("twin_prime_count_approx"); break;
      case 13: _vcallsub_with_pp("semiprime_count_approx"); break;
      case 14: _vcallsub_with_pp("lucky_count_upper"); break;
      case 15: _vcallsub_with_pp("lucky_count_lower"); break;
      case 16:
      default: _vcallsub_with_pp("lucky_count_approx"); break;
    }
    objectify_result(aTHX_ svn, ST(0));
    return; /* skip implicit PUTBACK */


void sum_primes(IN SV* svlo, IN SV* svhi = 0)
  PREINIT:
    UV lo = 2, hi;
  PPCODE:
    if ((items == 1 && _validate_and_set(&hi, aTHX_ svlo, IFLAG_POS)) ||
        (items == 2 && _validate_and_set(&lo, aTHX_ svlo, IFLAG_POS) && _validate_and_set(&hi, aTHX_ svhi, IFLAG_POS))) {
      UV count = 0;
      int retok = 1;
      /* 32/64-bit, Legendre or table-accelerated sieving. */
      retok = sum_primes(lo, hi, &count);
      /* If that didn't work, try the 128-bit version if supported. */
      if (retok == 0 && HAVE_SUM_PRIMES128) {
        UV hicount, lo_hic, lo_loc;
        retok = sum_primes128(hi, &hicount, &count);
        if (retok == 1 && lo > 2) {
          retok = sum_primes128(lo-1, &lo_hic, &lo_loc);
          hicount -= lo_hic;
          if (count < lo_loc) hicount--;
          count -= lo_loc;
        }
        if (retok == 1 && hicount > 0)
          RETURN_128(hicount, count);
      }
      if (retok == 1)
        XSRETURN_UV(count);
    }
    _vcallsub_with_pp("sum_primes");
    return;

void random_prime(IN SV* svlo, IN SV* svhi = 0)
  PREINIT:
    UV lo = 2, hi, ret;
    dMY_CXT;
  PPCODE:
    if ((items == 1 && _validate_and_set(&hi, aTHX_ svlo, IFLAG_POS)) ||
        (items == 2 && _validate_and_set(&lo, aTHX_ svlo, IFLAG_POS) && _validate_and_set(&hi, aTHX_ svhi, IFLAG_POS))) {
      ret = random_prime(MY_CXT.randcxt,lo,hi);
      if (ret) XSRETURN_UV(ret);
      else     XSRETURN_UNDEF;
    }
    _vcallsub_with_gmpobj(0.44,"random_prime");
    objectify_result(aTHX_ ST(0), ST(0));
    XSRETURN(1);

void print_primes(IN SV* svlo, IN SV* svhi = 0, IN int infd = -1)
  PREINIT:
    UV lo = 2, hi;
  PPCODE:
    if ((items == 1 && _validate_and_set(&hi, aTHX_ svlo, IFLAG_POS)) ||
        (items >= 2 && _validate_and_set(&lo, aTHX_ svlo, IFLAG_POS) && _validate_and_set(&hi, aTHX_ svhi, IFLAG_POS))) {
      if (lo <= hi) {
        int fd = (infd == -1) ? fileno(stdout) : infd;
        print_primes(lo, hi, fd);
      }
    } else {
      _vcallsub_with_pp("print_primes");
    }
    return;

UV
_LMO_pi(IN UV n)
  ALIAS:
    _legendre_pi = 1
    _meissel_pi = 2
    _lehmer_pi = 3
    _LMOS_pi = 4
    _segment_pi = 5
  PREINIT:
    UV ret;
  CODE:
    switch (ix) {
      case 0: ret = LMO_prime_count(n); break;
      case 1: ret = legendre_prime_count(n); break;
      case 2: ret = meissel_prime_count(n); break;
      case 3: ret = lehmer_prime_count(n); break;
      case 4: ret = LMOS_prime_count(n); break;
      default:ret = segment_prime_count(2,n); break;
    }
    RETVAL = ret;
  OUTPUT:
    RETVAL




void
sieve_primes(IN UV low, IN UV high)
  ALIAS:
    trial_primes = 1
    erat_primes = 2
    segment_primes = 3
  PREINIT:
    AV* av;
  PPCODE:
    CREATE_AV(av);
    if ((low <= 2) && (high >= 2)) av_push(av, newSVuv( 2 ));
    if ((low <= 3) && (high >= 3)) av_push(av, newSVuv( 3 ));
    if ((low <= 5) && (high >= 5)) av_push(av, newSVuv( 5 ));
    if (low < 7)  low = 7;
    if (low <= high) {
      if (ix == 0) {                          /* Sieve with primary cache */
        START_DO_FOR_EACH_PRIME(low, high) {
          av_push(av,newSVuv(p));
        } END_DO_FOR_EACH_PRIME
      } else if (ix == 1) {                   /* Trial */
        for (low = next_prime(low-1);
             low <= high && low != 0;
             low = next_prime(low) ) {
          av_push(av,newSVuv(low));
        }
      } else if (ix == 2) {                   /* Erat with private memory */
        unsigned char* sieve = sieve_erat30(high);
        START_DO_FOR_EACH_SIEVE_PRIME( sieve, 0, low, high ) {
           av_push(av,newSVuv(p));
        } END_DO_FOR_EACH_SIEVE_PRIME
        Safefree(sieve);
      } else if (ix == 3) {        /* Segment */
        unsigned char* segment;
        UV seg_base, seg_low, seg_high;
        void* ctx = start_segment_primes(low, high, &segment);
        while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
          START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
            av_push(av,newSVuv( p ));
          END_DO_FOR_EACH_SIEVE_PRIME
        }
        end_segment_primes(ctx);
      }
    }
    return; /* skip implicit PUTBACK */


void primes(IN SV* svlo, IN SV* svhi = 0)
  PREINIT:
    AV* av;
    UV lo = 0, hi, i;
  PPCODE:
    if ((items == 1 && _validate_and_set(&hi, aTHX_ svlo, IFLAG_POS)) ||
        (items == 2 && _validate_and_set(&lo, aTHX_ svlo, IFLAG_POS) && _validate_and_set(&hi, aTHX_ svhi, IFLAG_POS))) {
      CREATE_AV(av);
      if ((lo <= 2) && (hi >= 2)) av_push(av, newSVuv( 2 ));
      if ((lo <= 3) && (hi >= 3)) av_push(av, newSVuv( 3 ));
      if ((lo <= 5) && (hi >= 5)) av_push(av, newSVuv( 5 ));
      if (lo < 7)  lo = 7;
      if (lo <= hi) {
        if ( hi-lo <= 10
             || (hi >  100000000UL && hi-lo <=  330)
             || (hi > 4000000000UL && hi-lo <= 1500)
           ) {
          for (i = !(lo&1); i <= hi-lo; i += 2)
            if (is_prime(lo+i))
              av_push(av,newSVuv(lo+i));
        } else if (hi < (65536*30) ||  hi <= get_prime_cache(0,0)) {
          START_DO_FOR_EACH_PRIME(lo, hi) {
            av_push(av,newSVuv(p));
          } END_DO_FOR_EACH_PRIME
        } else {
          unsigned char* segment;
          UV seg_base, seg_low, seg_high;
          void* ctx = start_segment_primes(lo, hi, &segment);
          while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
            START_DO_FOR_EACH_SIEVE_PRIME(segment, seg_base, seg_low, seg_high)
              av_push(av,newSVuv( p ));
            END_DO_FOR_EACH_SIEVE_PRIME
          }
          end_segment_primes(ctx);
        }
      }
    } else {
      _vcallsub_with_pp("primes");
    }
    return;

void almost_primes(IN UV k, IN SV* svlo, IN SV* svhi = 0)
  ALIAS:
    omega_primes = 1
  PREINIT:
    AV* av;
    UV lo = 1, hi, i, n, *S;
  PPCODE:
    if ((items == 1 && _validate_and_set(&hi, aTHX_ svlo, IFLAG_POS)) ||
        (items >= 2 && _validate_and_set(&lo, aTHX_ svlo, IFLAG_POS) && _validate_and_set(&hi, aTHX_ svhi, IFLAG_POS))) {
      CREATE_AV(av);
      S = 0;
      if (ix == 0) n = generate_almost_primes(&S, k, lo, hi);
      else         n = range_omega_prime_sieve(&S, k, lo, hi);
      for (i = 0; i < n; i++)
        av_push(av, newSVuv(S[i]));
      if (S != 0) Safefree(S);
    } else {
      if (ix == 0) _vcallsub_with_pp("almost_primes");
      else         _vcallsub_with_pp("omega_primes");
    }
    return;


void prime_powers(IN SV* svlo, IN SV* svhi = 0)
  ALIAS:
    twin_primes = 1
    semi_primes = 2
    ramanujan_primes = 3
  PREINIT:
    AV* av;
    UV lo = 0, hi, i, num, *L;
  PPCODE:
    if ((items == 1 && _validate_and_set(&hi, aTHX_ svlo, IFLAG_POS)) ||
        (items == 2 && _validate_and_set(&lo, aTHX_ svlo, IFLAG_POS) && _validate_and_set(&hi, aTHX_ svhi, IFLAG_POS))) {
      CREATE_AV(av);
      if (ix == 0) {         /* Prime power */
        if ((lo <= 2) && (hi >= 2)) av_push(av, newSVuv( 2 ));
        if ((lo <= 3) && (hi >= 3)) av_push(av, newSVuv( 3 ));
        if ((lo <= 4) && (hi >= 4)) av_push(av, newSVuv( 4 ));
        if ((lo <= 5) && (hi >= 5)) av_push(av, newSVuv( 5 ));
      } else if (ix == 1) {  /* Twin */
        if ((lo <= 3) && (hi >= 3)) av_push(av, newSVuv( 3 ));
        if ((lo <= 5) && (hi >= 5)) av_push(av, newSVuv( 5 ));
      } else if (ix == 2) {  /* Semi */
        if ((lo <= 4) && (hi >= 4)) av_push(av, newSVuv( 4 ));
        if ((lo <= 6) && (hi >= 6)) av_push(av, newSVuv( 6 ));
      } else if (ix == 3) {  /* Ramanujan */
        if ((lo <= 2) && (hi >= 2)) av_push(av, newSVuv( 2 ));
      }
      if (lo < 7)  lo = 7;
      if (lo <= hi) {
        switch (ix) {
          case  0: num = prime_power_sieve(&L,lo,hi);           break;
          case  1: num = range_twin_prime_sieve(&L,lo,hi);      break;
          case  2: num = range_semiprime_sieve(&L,lo,hi);       break;
          case  3: num = range_ramanujan_prime_sieve(&L,lo,hi); break;
          default: num = 0; L = 0; break;
        }
        for (i = 0; i < num; i++)
          av_push(av,newSVuv(L[i]));
        Safefree(L);
      }
    } else {
      switch (ix) {
        case  0: _vcallsub_with_pp("prime_powers"); break;
        case  1: _vcallsub_with_pp("twin_primes"); break;
        case  2: _vcallsub_with_pp("semi_primes"); break;
        case  3: _vcallsub_with_pp("ramanujan_primes"); break;
        default: break;
      }
    }
    return;

void
lucky_numbers(IN SV* svlo, IN SV* svhi = 0)
  PREINIT:
    AV* av;
    UV lo = 0, hi, i, nlucky = 0;
  PPCODE:
    if ((items == 1 && _validate_and_set(&hi, aTHX_ svlo, IFLAG_POS)) ||
        (items == 2 && _validate_and_set(&lo, aTHX_ svlo, IFLAG_POS) && _validate_and_set(&hi, aTHX_ svhi, IFLAG_POS))) {
      CREATE_AV(av);
      if (lo == 0 && hi <= UVCONST(4000000000)) {
        uint32_t* lucky = lucky_sieve32(&nlucky, hi);
        for (i = 0; i < nlucky; i++)
          av_push(av,newSVuv(lucky[i]));
        Safefree(lucky);
      } else {
        UV* lucky = lucky_sieve_range(&nlucky, lo, hi);
        for (i = 0; i < nlucky; i++)
          av_push(av,newSVuv(lucky[i]));
        Safefree(lucky);
      }
    } else {
      _vcallsub_with_pp("lucky_numbers");
    }
    return;

void powerful_numbers(IN SV* svlo, IN SV* svhi = 0, IN UV k = 2)
  PREINIT:
    AV* av;
    UV lo = 1, hi, i, npowerful, *powerful;
  PPCODE:
    if ((items == 1 && _validate_and_set(&hi, aTHX_ svlo, IFLAG_POS)) ||
        (items >= 2 && _validate_and_set(&lo, aTHX_ svlo, IFLAG_POS) && _validate_and_set(&hi, aTHX_ svhi, IFLAG_POS))) {
      CREATE_AV(av);
      powerful = powerful_numbers_range(&npowerful, lo, hi, k);
      for (i = 0; i < npowerful; i++)
        av_push(av,newSVuv(powerful[i]));
      Safefree(powerful);
    } else {
      _vcallsub_with_pp("powerful_numbers");
    }
    return;

void
sieve_range(IN SV* svn, IN UV width, IN UV depth)
  PREINIT:
    int status;
    UV i, n;
  PPCODE:
    /* Return index of every n unless it is a composite with factor > depth */
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_POS);
    if (status == 1) {
      if ((n+width) < n) {
        status = 0;   /* range will overflow */
      } else { /* TODO: actually sieve */
        for (i = (n<2)?2-n:0; i < width; i++)
          if (is_rough(n+i, (depth+1) >= (n+i) ? n+i : depth+1))
            XPUSHs(sv_2mortal(newSVuv( i )));
      }
    }
    if (status != 1) {
      _vcallsubn(aTHX_ GIMME_V, VCALL_GMP|VCALL_PP, "sieve_range", items, 36);
      return;
    }

void
sieve_prime_cluster(IN SV* svlo, IN SV* svhi, ...)
  PREINIT:
    uint32_t nc, cl[100];
    UV i, lo, hi, cval, nprimes, *list;
    int done;
  PPCODE:
    nc = items-1;
    if (items > 100) croak("sieve_prime_cluster: too many entries");
    cl[0] = 0;
    for (i = 1; i < nc; i++) {
      if (!_validate_and_set(&cval, aTHX_ ST(1+i), IFLAG_POS))
        croak("sieve_prime_cluster: cluster values must be standard integers");
      if (cval & 1) croak("sieve_prime_cluster: values must be even");
      if (cval > 2147483647UL) croak("sieve_prime_cluster: values must be 31-bit");
      if (cval <= cl[i-1]) croak("sieve_prime_cluster: values must be increasing");
      cl[i] = cval;
    }
    done = 0;
    if (_validate_and_set(&lo, aTHX_ svlo, IFLAG_POS) &&
        _validate_and_set(&hi, aTHX_ svhi, IFLAG_POS)) {
      list = sieve_cluster(lo, hi, nc, cl, &nprimes);
      if (list != 0) {
        done = 1;
        EXTEND(SP, (IV)nprimes);
        for (i = 0; i < nprimes; i++)
          PUSHs(sv_2mortal(newSVuv( list[i] )));
        Safefree(list);
      }
    }
    if (!done) {
      _vcallsubn(aTHX_ GIMME_V, VCALL_GMP|VCALL_PP, "sieve_prime_cluster", items, 34);
      return;
    }

void is_pseudoprime(IN SV* svn, ...)
  ALIAS:
    is_euler_pseudoprime = 1
    is_strong_pseudoprime = 2
  PREINIT:
    int i, status, ret = 0;
    UV n, base;
  PPCODE:
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
    if (status == 1) {
      if (n < 3) {
        ret = (n == 2);
      } else if (ix == 2 && !(n&1)) {
        ret = 0;
      } else if (items == 1) {
        ret = (ix == 0) ? is_pseudoprime(n, 2) :
              (ix == 1) ? is_euler_pseudoprime(n, 2) :
                          is_strong_pseudoprime(n, 2);
      } else {
        for (i = 1, ret = 1;  i < items && ret == 1; i++) {
          status = _validate_and_set(&base, aTHX_ ST(i), IFLAG_POS);
          if (status != 1) break;
          ret = (ix == 0) ? is_pseudoprime(n, base) :
                (ix == 1) ? is_euler_pseudoprime(n, base) :
                            is_strong_pseudoprime(n, base);
        }
      }
    }
    if (status != 0)  RETURN_NPARITY(ret);
    _vcallsub_with_gmp( (items == 1) ? 0.74 : 0.41,
                        (ix == 0) ? "is_pseudoprime" :
                        (ix == 1) ? "is_euler_pseudoprime" :
                                    "is_strong_pseudoprime");
    return; /* skip implicit PUTBACK */


void is_prime(IN SV* svn)
  ALIAS:
    is_prob_prime = 1
    is_provable_prime = 2
    is_bpsw_prime = 3
    is_aks_prime = 4
    is_lucas_pseudoprime = 5
    is_strong_lucas_pseudoprime = 6
    is_extra_strong_lucas_pseudoprime = 7
    is_frobenius_underwood_pseudoprime = 8
    is_frobenius_khashin_pseudoprime = 9
    is_catalan_pseudoprime = 10
    is_euler_plumb_pseudoprime = 11
    is_ramanujan_prime = 12
    is_semiprime = 13
    is_chen_prime = 14
    is_mersenne_prime = 15
  PREINIT:
    int status, ret;
    UV n;
  PPCODE:
    ret = 0;
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
    if (status == 1) {
      switch (ix) {
        case 0:
        case 1:
        case 2:  ret = is_prime(n); break;
        case 3:  ret = BPSW(n); break;
        case 4:  ret = is_aks_prime(n); break;
        case 5:  ret = is_lucas_pseudoprime(n, 0); break;
        case 6:  ret = is_lucas_pseudoprime(n, 1); break;
        case 7:  ret = is_lucas_pseudoprime(n, 3); break;
        case 8:  ret = is_frobenius_underwood_pseudoprime(n); break;
        case 9:  ret = is_frobenius_khashin_pseudoprime(n); break;
        case 10: ret = is_catalan_pseudoprime(n); break;
        case 11: ret = is_euler_plumb_pseudoprime(n); break;
        case 12: ret = is_ramanujan_prime(n); break;
        case 13: ret = is_semiprime(n); break;
        case 14: ret = is_chen_prime(n); break;
        case 15: ret = is_mersenne_prime(n);  if (ret == -1) status = 0; break;
        default: break;
      }
    }
    if (status != 0)  RETURN_NPARITY(ret);
    switch (ix) {
      case 0: _vcallsub_with_gmp(0.01,"is_prime");       break;
      case 1: _vcallsub_with_gmp(0.01,"is_prob_prime");  break;
      case 2: _vcallsub_with_gmp(0.04,"is_provable_prime");  break;
      case 3: _vcallsub_with_gmp(0.17,"is_bpsw_prime");  break;
      case 4: _vcallsub_with_gmp(0.16,"is_aks_prime"); break;
      case 5: _vcallsub_with_gmp(0.01,"is_lucas_pseudoprime"); break;
      case 6: _vcallsub_with_gmp(0.01,"is_strong_lucas_pseudoprime"); break;
      case 7: _vcallsub_with_gmp(0.01,"is_extra_strong_lucas_pseudoprime"); break;
      case 8: _vcallsub_with_gmp(0.13,"is_frobenius_underwood_pseudoprime"); break;
      case 9: _vcallsub_with_gmp(0.30,"is_frobenius_khashin_pseudoprime"); break;
      case 10:_vcallsub_with_gmp(0.00,"is_catalan_pseudoprime"); break;
      case 11:_vcallsub_with_gmp(0.39,"is_euler_plumb_pseudoprime"); break;
      case 12:_vcallsub_with_gmp(0.00,"is_ramanujan_prime"); break;
      case 13:_vcallsub_with_gmp(0.42,"is_semiprime"); break;
      case 14:_vcallsub_with_gmp(0.00,"is_chen_prime"); break;
      case 15:_vcallsub_with_gmp(0.28,"is_mersenne_prime"); break;
      default: break;
    }
    return; /* skip implicit PUTBACK */

void
is_perrin_pseudoprime(IN SV* svn, IN UV k = 0)
  ALIAS:
    is_almost_extra_strong_lucas_pseudoprime = 1
    is_delicate_prime = 2
  PREINIT:
    int status, ret;
    UV n;
  PPCODE:
    if (k < 0) croak("second argument must not be negative");
    /*  ix = 0    k = 0 - 3       n below 2 returns 0 for all k
     *  ix = 1    k = 0 - 256     n below 2 returns 0 for all k
     *  ix = 2    k = 0 - 2^32    n below 2 returns 0 for all k
     */
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
    ret = 0;
    if (status == 1) {
      switch (ix) {
        case 0:  ret = is_perrin_pseudoprime(n, k); break;
        case 1:  ret = is_almost_extra_strong_lucas_pseudoprime(n, (k < 1) ? 1 : k); break;
        case 2:  ret = is_delicate_prime(n, (k<1) ? 10 : k);
                 if (ret < 0) status = 0; break;
        default: break;
      }
    }
    if (status != 0) RETURN_NPARITY(ret);
    switch (ix) {
      case 0: _vcallsub_with_gmp( (k == 0) ? 0.20 : 0.40, "is_perrin_pseudoprime"); break;
      case 1: _vcallsub_with_gmp(0.13,"is_almost_extra_strong_lucas_pseudoprime"); break;
      case 2: _vcallsub_with_gmp(0.00,"is_delicate_prime"); break;
      default: break;
    }
    return;

void
is_frobenius_pseudoprime(IN SV* svn, IN IV P = 0, IN IV Q = 0)
  PREINIT:
    int status;
    UV n;
  PPCODE:
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
    if (status != 0)
      RETURN_NPARITY((status == 1) ?  is_frobenius_pseudoprime(n, P, Q)  :  0);
    _vcallsub_with_gmp(0.24,"is_frobenius_pseudoprime");
    return;

void
miller_rabin_random(IN SV* svn, IN IV bases = 1, IN char* seed = 0)
  PREINIT:
    int status;
    UV n;
    dMY_CXT;
  PPCODE:
    if (bases < 0) croak("miller_rabin_random: number of bases must be positive");
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
    if (status == -1) RETURN_NPARITY(0);
    if (seed == 0 && status == 1)
      RETURN_NPARITY( is_mr_random(MY_CXT.randcxt, n, bases) );
    _vcallsub_with_gmp(0.46,"miller_rabin_random");
    return;

void is_gaussian_prime(IN SV* sva, IN SV* svb)
  PREINIT:
    UV a, b;
  PPCODE:
    if (_validate_and_set(&a, aTHX_ sva, IFLAG_ABS) &&
        _validate_and_set(&b, aTHX_ svb, IFLAG_ABS)) {
      if (a == 0) RETURN_NPARITY( ((b % 4) == 3) ? is_prime(b) : 0 );
      if (b == 0) RETURN_NPARITY( ((a % 4) == 3) ? is_prime(a) : 0 );
      if (a < HALF_WORD && b < HALF_WORD) {
        UV aa = a*a, bb = b*b;
        if (UV_MAX-aa >= bb)
          RETURN_NPARITY( is_prime(aa+bb) );
      }
    }
    _vcallsub_with_gmp(0.52,"is_gaussian_prime");
    return;


void
gcd(...)
  PROTOTYPE: @
  ALIAS:
    lcm = 1
    vecmin = 2
    vecmax = 3
    vecsum = 4
    vecprod = 5
  PREINIT:
    int i, status = 1;
    UV ret, nullv, n;
  PPCODE:
    if (ix == 2 || ix == 3) {
      UV retindex = 0;
      int sign, minmax = (ix == 2);
      if (items == 0) XSRETURN_UNDEF;
      if (items == 1) XSRETURN(1);
      if (items > 1 && (status = _validate_and_set(&ret, aTHX_ ST(0), IFLAG_ANY))) {
        sign = status;
        for (i = 1; i < items; i++) {
          status = _validate_and_set(&n, aTHX_ ST(i), IFLAG_ANY);
          if (status == 0) break;
          if (( (sign == -1 && status == 1) ||
                (n >= ret && sign == status)
              ) ? !minmax : minmax ) {
            sign = status;
            ret = n;
            retindex = i;
          }
        }
      }
      if (status != 0) {
        ST(0) = ST(retindex);
        XSRETURN(1);
      }
    } else if (ix == 4) {
      UV lo = 0;
      IV hi = 0;
      for (ret = i = 0; i < items; i++) {
        status = _validate_and_set(&n, aTHX_ ST(i), IFLAG_ANY);
        if (status == 0) break;
        if (status == 1) {
          hi += (n > (UV_MAX - lo));
        } else {
          /* TODO: This shouldn't be needed with validate_and_set */
          if (UV_MAX-n == (UV)IV_MAX) { status = 0; break; }  /* IV Overflow */
          hi -= ((UV_MAX-n) >= lo);
        }
        lo += n;
      }
      if (status != 0 && hi != 0) {
        if (hi == -1 && lo > IV_MAX) XSRETURN_IV((IV)lo);
        else                         RETURN_128(hi, lo);
      }
      ret = lo;
    } else if (ix == 5) {
      int sign = 1;
      ret = 1;
      for (i = 0; i < items; i++) {
        status = _validate_and_set(&n, aTHX_ ST(i), IFLAG_ANY);
        if (status == 0) break;
        if (ret > 0 && n > UV_MAX/ret) { status = 0; break; }
        sign *= status;
        ret *= n;
      }
      if (sign == -1 && status != 0) {
        if (ret <= (UV)IV_MAX)  XSRETURN_IV(neg_iv(ret));
        else                    status = 0;
      }
    } else {
      /* For each arg, while valid input, validate+gcd/lcm.  Shortcut stop. */
      if (ix == 0) { ret = 0; nullv = 1; }
      else         { ret = 1; nullv = 0; }
      for (i = 0; i < items && ret != nullv && status != 0; i++) {
        status = _validate_and_set(&n, aTHX_ ST(i), IFLAG_ABS);
        if (status == 0) break;
        if (i == 0) {
          ret = n;
        } else {
          UV gcd = gcd_ui(ret, n);
          if (ix == 0) {
            ret = gcd;
          } else {
            n /= gcd;
            if (n <= (UV_MAX / ret) )    ret *= n;
            else                         status = 0;   /* Overflow */
          }
        }
      }
    }
    if (status != 0)
      XSRETURN_UV(ret);
    /* For min/max, use string compare if not an object */
    if ((ix == 2 || ix == 3) && !sv_isobject(ST(0))) {
      int retindex = 0;
      int minmax = (ix == 2);
      STRLEN alen, blen;
      char *aptr, *bptr;
      aptr = SvPV(ST(0), alen);
      (void) strnum_minmax(minmax, 0, 0, aptr, alen);
      for (i = 1; i < items; i++) {
        bptr = SvPV(ST(i), blen);
        if (strnum_minmax(minmax, aptr, alen, bptr, blen)) {
          aptr = bptr;
          alen = blen;
          retindex = i;
        }
      }
      ST(0) = ST(retindex);
      XSRETURN(1);
    }
    switch (ix) {
      case 0: _vcallsub_with_gmp(0.17,"gcd");   break;
      case 1: _vcallsub_with_gmp(0.17,"lcm");   break;
      case 2: _vcallsub_with_gmp(0.00,"vecmin"); break;
      case 3: _vcallsub_with_gmp(0.00,"vecmax"); break;
      case 4: _vcallsub_with_pp("vecsum");  break;
      case 5:
      default:_vcallsub_with_pp("vecprod");  break;
    }
    return; /* skip implicit PUTBACK */

void
vecextract(IN SV* x, IN SV* svm)
  PREINIT:
    AV* av;
    UV mask, i = 0;
  PPCODE:
    if ((!SvROK(x)) || (SvTYPE(SvRV(x)) != SVt_PVAV))
      croak("vecextract first argument must be an array reference");
    av = (AV*) SvRV(x);
    if (SvROK(svm) && SvTYPE(SvRV(svm)) == SVt_PVAV) {
      AV* avm = (AV*) SvRV(svm);
      int j, mlen = av_len(avm);
      for (j = 0; j <= mlen; j++) {
        SV** iv = av_fetch(avm, j, 0);
        if (iv && SvTYPE(*iv) == SVt_IV) {
          SV **v = av_fetch(av, SvIV(*iv), 0);
          if (v) XPUSHs(*v);
        }
      }
    } else if (_validate_and_set(&mask, aTHX_ svm, IFLAG_POS)) {
      while (mask) {
        if (mask & 1) {
          SV** v = av_fetch(av, i, 0);
          if (v) XPUSHs(*v);
        }
        i++;
        mask >>= 1;
      }
    } else {
      _vcallsubn(aTHX_ GIMME_V, VCALL_PP, "vecextract", items, 0);
      return;
    }

void
vecequal(IN SV* a, IN SV* b)
  PREINIT:
    int res;
  PPCODE:
    res = _compare_array_refs(aTHX_ a, b);
    if (res == -1)
      croak("vecequal element not scalar or array reference");
    RETURN_NPARITY(res);
    return;

void
vecmex(...)
  ALIAS:
    vecpmex = 1
  PROTOTYPE: @
  PREINIT:
    char *setv;
    int i, status = 1;
    UV min, n;
    uint32_t mask;
  PPCODE:
    if (ix == 0) {
      min = 0;
      mask = IFLAG_POS;
    } else {
      min = 1;
      mask = IFLAG_POS | IFLAG_NONZERO;
    }
    if (items == 0)
      XSRETURN_UV(min);
    Newz(0, setv, items, char);
    for (i = 0; i < items; i++) {
      status = _validate_and_set(&n, aTHX_ ST(i), mask);
      /* Ignore any bigint */
      if (status == 1 && n-min < (UV)items)
        setv[n-min] = 1;
    }
    for (i = 0; i < items; i++)
      if (setv[i] == 0)
        break;
    Safefree(setv);
    XSRETURN_UV(i+min);

void
frobenius_number(...)
  PROTOTYPE: @
  PREINIT:
    int i, found1 = 0;
    UV fn, n, *A;
  PPCODE:
    if (items == 0) XSRETURN_UNDEF;
    Newz(0, A, items, UV);
    for (i = 0; i < items; i++) {
      if (!_validate_and_set(&n, aTHX_ ST(i), IFLAG_POS | IFLAG_NONZERO)) break;
      if (n == 1) { found1 = 1; break; }
      A[i] = n;
    }
    if (i == items) {
      fn = frobenius_number(A, i);
      Safefree(A);
      if (fn == 0) XSRETURN_UNDEF;
      if (fn != UV_MAX) XSRETURN_UV(fn);
    } else {
      Safefree(A);
      if (found1) XSRETURN_IV(-1);
    }
    (void)_vcallsubn(aTHX_ GIMME_V, VCALL_PP, "frobenius_number", items, 0);
    return;

void
chinese(...)
  ALIAS:
    chinese2 = 1
  PROTOTYPE: @
  PREINIT:
    int i, status, astatus, nstatus;
    UV ret, lcm, *an;
    SV **psva, **psvn;
  PPCODE:
    status = 1;
    New(0, an, 2*items, UV);
    ret = 0;
    for (i = 0; i < items; i++) {
      AV* av;
      if (!SvROK(ST(i)) || SvTYPE(SvRV(ST(i))) != SVt_PVAV || av_len((AV*)SvRV(ST(i))) != 1)
        croak("chinese arguments are two-element array references");
      av = (AV*) SvRV(ST(i));
      psva = av_fetch(av, 0, 0);
      psvn = av_fetch(av, 1, 0);
      if (psva == 0 || psvn == 0) { status = 0; break; }
      astatus = _validate_and_set(an+i, aTHX_ *psva, IFLAG_ANY);
      nstatus = _validate_and_set(an+i+items, aTHX_ *psvn, IFLAG_ABS);
      if (astatus == 0 || nstatus == 0) { status = 0; break; }
      if (an[i+items] == 0) {
        XPUSHs(&PL_sv_undef);
        if (ix == 1) XPUSHs(&PL_sv_undef);
        XSRETURN(1 + ix);
      }
      _mod_with(an+i, astatus, an[i+items]);
    }
    if (status)
      status = chinese(&ret, &lcm, an, an+items, items);
    Safefree(an);
    if (status) {
      if (ix == 0) {
        if (status < 0)  XSRETURN_UNDEF;
        else             XSRETURN_UV(ret);
      } else {
        if (status < 0) {
          XPUSHs(&PL_sv_undef);
          XPUSHs(&PL_sv_undef);
        } else {
          XPUSHs(sv_2mortal(newSVuv( ret )));
          XPUSHs(sv_2mortal(newSVuv( lcm )));
        }
        XSRETURN(2);
      }
    }
    if (ix == 0) {
      psvn = av_fetch((AV*) SvRV(ST(0)), 1, 0);
      _vcallsub_with_gmpobj(0.32,"chinese");
      objectify_result(aTHX_  (psvn ? *psvn : 0), ST(0));
    } else {
      (void)_vcallsubn(aTHX_ GIMME_V, VCALL_PP, "chinese2", 2, 0);
    }
    return; /* skip implicit PUTBACK */

void cornacchia(IN SV* svd, IN SV* svn)
  PREINIT:
    UV d, n, x, y;
  PPCODE:
    if (_validate_and_set(&d, aTHX_ svd, IFLAG_POS) &&
        _validate_and_set(&n, aTHX_ svn, IFLAG_POS) ) {
      if (!cornacchia(&x, &y, d, n))  XSRETURN_UNDEF;
      PUSHs(sv_2mortal(newSVuv( x )));
      PUSHs(sv_2mortal(newSVuv( y )));
    } else {
      _vcallsubn(aTHX_ GIMME_V, VCALL_PP, "cornacchia", items, 0);
      return;
    }

void lucas_sequence(...)
  PREINIT:
    UV U, V, Qk,  n, P, Q, k;
  PPCODE:
    if (items != 4) croak("lucas_sequence: n, P, Q, k");
    if (_validate_and_set(&n, aTHX_ ST(0), IFLAG_POS | IFLAG_NONZERO) &&
        _validate_and_set(&P, aTHX_ ST(1), IFLAG_ANY | IFLAG_IV) &&
        _validate_and_set(&Q, aTHX_ ST(2), IFLAG_ANY | IFLAG_IV) &&
        _validate_and_set(&k, aTHX_ ST(3), IFLAG_POS)) {
      lucas_seq(&U, &V, &Qk, n, (IV)P, (IV)Q, k);
      PUSHs(sv_2mortal(newSVuv( U )));  /* 4 args in, 3 out, no EXTEND needed */
      PUSHs(sv_2mortal(newSVuv( V )));
      PUSHs(sv_2mortal(newSVuv( Qk )));
    } else {
      _vcallsubn(aTHX_ GIMME_V, VCALL_PP, "lucas_sequence", items, 0);
      return;
    }

void lucasuvmod(IN SV* svp, IN SV* svq, IN SV* svk, IN SV* svn)
  ALIAS:
    lucasumod = 1
    lucasvmod = 2
  PREINIT:
    int pstatus, qstatus;
    UV P, Q, k, n, U, V;
  PPCODE:
    pstatus = _validate_and_set(&P, aTHX_ svp, IFLAG_ANY);
    qstatus = _validate_and_set(&Q, aTHX_ svq, IFLAG_ANY);
    if ((pstatus != 0) && (qstatus != 0) &&
        _validate_and_set(&k, aTHX_ svk, IFLAG_POS) &&
        _validate_and_set(&n, aTHX_ svn, IFLAG_ABS)
        ) {
      if (n == 0) XSRETURN_UNDEF;
      P = (pstatus == 1)  ?  P % n  :  ivmod((IV)P,n);
      Q = (qstatus == 1)  ?  Q % n  :  ivmod((IV)Q,n);
      switch (ix) {
        case 0:  lucasuvmod(&U, &V, P, Q, k, n);
                 PUSHs(sv_2mortal(newSVuv( U )));
                 PUSHs(sv_2mortal(newSVuv( V )));
                 break;
        case 1:  XSRETURN_UV(lucasumod(P, Q, k, n)); break;
        case 2:
        default: XSRETURN_UV(lucasvmod(P, Q, k, n)); break;
      }
    } else {
      if (ix == 0) {
        _vcallsubn(aTHX_ GIMME_V, VCALL_PP, "lucasuvmod", items, 0);
      } else {
        _vcallsub_with_gmp(0.53, (ix == 1) ? "lucasumod" : "lucasvmod");
        objectify_result(aTHX_ svn, ST(0));
      }
      return;
    }

void lucasuv(IN SV* svp, IN SV* svq, IN SV* svk)
  ALIAS:
    lucasu = 1
    lucasv = 2
  PREINIT:
    UV k;
    IV P, Q, U, V;
  PPCODE:
    if (_validate_and_set((UV*)&P, aTHX_ svp, IFLAG_IV) &&
        _validate_and_set((UV*)&Q, aTHX_ svq, IFLAG_IV) &&
        _validate_and_set(&k, aTHX_ svk, IFLAG_POS) &&
        lucasuv(&U, &V, P, Q, k)) {
      if (ix == 1)  XSRETURN_IV(U);     /* U = lucasu(P,Q,k) */
      if (ix == 2)  XSRETURN_IV(V);     /* V = lucasv(P,Q,k) */
      PUSHs(sv_2mortal(newSViv( U )));  /* (U,V) = lucasuv(P,Q,k) */
      PUSHs(sv_2mortal(newSViv( V )));
    } else {
      if (ix == 0) {
        _vcallsubn(aTHX_ GIMME_V, VCALL_PP, "lucasuv", items, 0);
      } else {
        _vcallsub_with_gmpobj(0.29, (ix==1) ? "lucasu" : "lucasv");
        objectify_result(aTHX_ ST(0), ST(0));
      }
      return;
    }


void is_sum_of_squares(IN SV* svn, IN UV k = 2)
  PREINIT:
    int status, ret;
    UV n;
  PPCODE:
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_ABS);
    if (status != 0) {
      switch (k) {
        case 0:  ret = (n==0);                     break;
        case 1:  ret = is_power(n,2);              break;
        case 2:  ret = is_sum_of_two_squares(n);   break;
        case 3:  ret = is_sum_of_three_squares(n); break;
        default: ret = 1;                          break;
      }
      RETURN_NPARITY(ret);
    }
    _vcallsub_with_pp("is_sum_of_squares");
    return;


void is_square(IN SV* svn)
  ALIAS:
    is_carmichael = 1
    is_quasi_carmichael = 2
    is_perfect_power = 3
    is_fundamental = 4
    is_lucky = 5
    is_practical = 6
    is_perfect_number = 7
    is_totient = 8
  PREINIT:
    int status, ret;
    UV n;
  PPCODE:
    ret = 0;
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
    if (ix == 3 && status == -1) {  /* is_perfect_power special behavior */
      n = neg_iv(n);
      ret = is_power(n,0);
      /* A power with exponent other than 0,1,2,4,8,... is ok */
      ret = (n == 1 || (ret > 2 && (ret & (ret-1)) != 0));
    }
    if (ix == 4 && status == -1) {  /* is_fundamental special behavior */
      n = neg_iv(n);
      ret = is_fundamental(n,1);
    }
    if (status == 1) {
      switch (ix) {
        case 0: ret = is_power(n,2); break;
        case 1: ret = is_carmichael(n); break;
        case 2: ret = is_quasi_carmichael(n); break;
        case 3: ret = is_perfect_power(n); break;
        case 4: ret = is_fundamental(n,0); break;
        case 5: ret = is_lucky(n); break;
        case 6: ret = is_practical(n); break;
        case 7: ret = is_perfect_number(n); break;
        case 8:
        default:ret = is_totient(n); break;
      }
    }
    if (status != 0)  RETURN_NPARITY(ret);
    switch (ix) {
      case  0: _vcallsub_with_gmp(0.47,"is_square"); break;
      case  1: _vcallsub_with_gmp(0.47,"is_carmichael"); break;
      case  2: _vcallsub_with_gmp(0.00,"is_quasi_carmichael"); break;
      case  3: _vcallsub_with_gmp(0.47,"is_perfect_power"); break;
      case  4: _vcallsub_with_gmp(0.47,"is_fundamental"); break;
      case  5: _vcallsub_with_gmp(0.48,"is_lucky"); break;
      case  6: _vcallsub_with_gmp(0.53,"is_practical"); break;
      case  7: _vcallsub_with_gmp(0.53,"is_perfect_number"); break;
      case  8:
      default: _vcallsub_with_gmp(0.47,"is_totient"); break;
    }
    return; /* skip implicit PUTBACK */

void squarefree_kernel(IN SV* svn)
  PREINIT:
    int status;
    UV n;
  PPCODE:
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
    if (status == -1)
      XSRETURN_IV( neg_iv(squarefree_kernel(neg_iv(n))) );
    if (status == 1)
      XSRETURN_UV( squarefree_kernel(n) );
    _vcallsub_with_pp("squarefree_kernel");
    return;

void is_powerfree(IN SV* svn, IN int k = 2)
  ALIAS:
    powerfree_sum = 1
    powerfree_part = 2
    powerfree_part_sum = 3
  PREINIT:
    int status;
    UV n, res;
  PPCODE:
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
    if (status == -1) {
      if (ix == 2)
        XSRETURN_IV( neg_iv(powerfree_part(neg_iv(n),k)) );
      RETURN_NPARITY(0);
    }
    if (status == 1) {
      switch (ix) {
        case 0:  res = is_powerfree(n,k);    break;
        case 1:  res = powerfree_sum(n,k);   break;
        case 2:  res = powerfree_part(n,k);  break;
        case 3:
        default: res = powerfree_part_sum(n,k);  break;
      }
      if (ix == 0)
        RETURN_NPARITY(res);
      if (res != 0 || n == 0)
        XSRETURN_UV(res);
      /* res is 0 and n > 0, so we overflowed.  Fall through to PP. */
    }
    switch (ix) {
      case  0: _vcallsub_with_gmp(0.53,"is_powerfree"); break;
      case  1: _vcallsub_with_pp("powerfree_sum"); break;
      case  2: _vcallsub_with_pp("powerfree_part"); break;
      case  3:
      default: _vcallsub_with_pp("powerfree_part_sum"); break;
    }
    return;

void powerfree_count(IN SV* svn, IN int k = 2)
  ALIAS:
    nth_powerfree = 1
  PREINIT:
    int status;
    UV n, res;
  PPCODE:
    status = _validate_and_set(&n, aTHX_ svn, (ix==0) ? IFLAG_ANY : IFLAG_POS);
    if (status != 0) {
      if (status == -1)
        XSRETURN_UV(0);
      if (ix == 0) {
        res = powerfree_count(n,k);
        XSRETURN_UV(res);
      } else {
        if (n == 0 || k < 2)
          XSRETURN_UNDEF;
        res = nth_powerfree(n,k);
        if (res != 0)
          XSRETURN_UV(res);
        /* if res = 0, overflow */
      }
    }
    switch (ix) {
      case  0: _vcallsub_with_gmp(0.53,"powerfree_count"); break;
      case  1: _vcallsub_with_gmp(0.53,"nth_powerfree"); break;
      default: break;
    }
    objectify_result(aTHX_ svn, ST(0));
    return;

void
is_power(IN SV* svn, IN UV k = 0, IN SV* svroot = 0)
  PREINIT:
    int status, ret;
    UV n;
  PPCODE:
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
    if (status != 0) {
      if (status == -1) {
        /* Negative n with even positive k return 0. */
        if (k > 0 && !(k & 1))  RETURN_NPARITY(0);
        n = neg_iv(n);
      }
      ret = is_power(n, k);
      if (status == -1 && k == 0) {
        ret >>= valuation(ret,2);
        if (ret == 1) ret = 0;
      }
      if (ret && svroot != 0) {
        UV root = rootint(n, k ? k : (UV)ret);
        if (!SvROK(svroot)) croak("is_power: third argument not a scalar reference");
        if (status == 1) sv_setuv(SvRV(svroot),  root);
        else             sv_setiv(SvRV(svroot), -root);
      }
      RETURN_NPARITY(ret);
    }
    if (svroot == 0) { _vcallsub_with_gmp(0.28, "is_power"); }
    else             { _vcallsub_with_pp("is_power"); }
    return;

void
is_prime_power(IN SV* svn, IN SV* svroot = 0)
  PREINIT:
    int status, ret;
    UV n, root;
  PPCODE:
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
    if (status != 0) {
      ret = (status == 1)  ?  prime_power(n, &root)  :  0;
      if (ret && svroot != 0) {
        if (!SvROK(svroot))croak("is_prime_power: second argument not a scalar reference");
        sv_setuv(SvRV(svroot), root);
      }
      RETURN_NPARITY(ret);
    }
    (void)_vcallsubn(aTHX_ G_SCALAR, (svroot == 0) ? (VCALL_GMP|VCALL_PP) : (VCALL_PP), "is_prime_power", items, 40);
    return;

void
is_polygonal(IN SV* svn, IN UV k, IN SV* svroot = 0)
  PREINIT:
    int status, result, overflow;
    UV n, root;
  PPCODE:
    if (svroot != 0 && !SvROK(svroot))
      croak("is_polygonal: third argument not a scalar reference");
    if (k < 3)
      croak("is_polygonal: k must be >= 3");

    status = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
    if (status == -1)
      RETURN_NPARITY(0);
    if (status == 1) {
      overflow = 0;
      root = polygonal_root(n, k, &overflow);
      result = (n == 0) || root;
      if (!overflow) {
        if (result && svroot != 0)
          sv_setuv(SvRV(svroot), root);
        RETURN_NPARITY(result);
      }
    }
    if (items != 3) { _vcallsub_with_gmp(0.47, "is_polygonal"); }
    else            { _vcallsub_with_pp("is_polygonal"); }
    return;


void inverse_li(IN SV* svn)
  PREINIT:
    UV n;
  PPCODE:
    if (_validate_and_set(&n, aTHX_ svn, IFLAG_POS)) {
      if (n < MPU_MAX_PRIME_IDX) /* Fall through to Perl if out of range. */
        XSRETURN_UV(inverse_li(n));
    }
    _vcallsub_with_pp("inverse_li");
    XSRETURN(1);

void nth_prime(IN SV* svn)
  ALIAS:
    nth_prime_upper = 1
    nth_prime_lower = 2
    nth_prime_approx = 3
  PREINIT:
    UV n, ret;
  PPCODE:
    if ( _validate_and_set(&n, aTHX_ svn, IFLAG_POS) &&
         n <= MPU_MAX_PRIME_IDX ) {
      if (n == 0) XSRETURN_UNDEF;
      switch (ix) {
        case 0:  ret = nth_prime(n); break;
        case 1:  ret = nth_prime_upper(n); break;
        case 2:  ret = nth_prime_lower(n); break;
        case 3:
        default: ret = nth_prime_approx(n); break;
      }
      XSRETURN_UV(ret);
    }
    switch (ix) {
      case 0:  _vcallsub_with_pp("nth_prime");          break;
      case 1:  _vcallsub_with_pp("nth_prime_upper");    break;
      case 2:  _vcallsub_with_pp("nth_prime_lower");    break;
      case 3:
      default: _vcallsub_with_pp("nth_prime_approx");   break;
    }
    return;

void nth_prime_power(IN SV* svn)
  ALIAS:
    nth_prime_power_upper = 1
    nth_prime_power_lower = 2
    nth_prime_power_approx = 3
  PREINIT:
    UV n, ret;
  PPCODE:
    if ( _validate_and_set(&n, aTHX_ svn, IFLAG_POS) &&
         n <= MPU_MAX_PRIME_IDX ) {
      if (n == 0) XSRETURN_UNDEF;
      switch (ix) {
        case 0:  ret = nth_prime_power(n); break;
        case 1:  ret = nth_prime_power_upper(n); break;
        case 2:  ret = nth_prime_power_lower(n); break;
        case 3:
        default: ret = nth_prime_power_approx(n); break;
      }
      XSRETURN_UV(ret);
    }
    switch (ix) {
      case 0:  _vcallsub_with_pp("nth_prime_power");          break;
      case 1:  _vcallsub_with_pp("nth_prime_power_upper");    break;
      case 2:  _vcallsub_with_pp("nth_prime_power_lower");    break;
      case 3:
      default: _vcallsub_with_pp("nth_prime_power_approx");   break;
    }
    return;

void nth_perfect_power(IN SV* svn)
  ALIAS:
    nth_perfect_power_upper = 1
    nth_perfect_power_lower = 2
    nth_perfect_power_approx = 3
  PREINIT:
    UV n, ret;
  PPCODE:
    if ( _validate_and_set(&n, aTHX_ svn, IFLAG_POS) &&
         n <= MPU_MAX_PERFECT_POW_IDX ) {
      if (n == 0) XSRETURN_UNDEF;
      switch (ix) {
        case 0:  ret = nth_perfect_power(n); break;
        case 1:  ret = nth_perfect_power_upper(n); break;
        case 2:  ret = nth_perfect_power_lower(n); break;
        case 3:
        default: ret = nth_perfect_power_approx(n); break;
      }
      XSRETURN_UV(ret);
    }
    switch (ix) {
      case 0:  _vcallsub_with_gmp(0.53,"nth_perfect_power");        break;
      case 1:  _vcallsub_with_pp("nth_perfect_power_upper");        break;
      case 2:  _vcallsub_with_pp("nth_perfect_power_lower");        break;
      case 3:
      default: _vcallsub_with_gmp(0.53,"nth_perfect_power_approx"); break;
    }
    objectify_result(aTHX_ svn, ST(0));
    return;

void nth_ramanujan_prime(IN SV* svn)
  ALIAS:
    nth_ramanujan_prime_upper = 1
    nth_ramanujan_prime_lower = 2
    nth_ramanujan_prime_approx = 3
  PREINIT:
    UV n, ret;
  PPCODE:
    if ( _validate_and_set(&n, aTHX_ svn, IFLAG_POS) &&
         n <= MPU_MAX_RMJN_PRIME_IDX ) {
      if (n == 0) XSRETURN_UNDEF;
      switch (ix) {
        case 0:  ret = nth_ramanujan_prime(n); break;
        case 1:  ret = nth_ramanujan_prime_upper(n); break;
        case 2:  ret = nth_ramanujan_prime_lower(n); break;
        case 3:
        default: ret = nth_ramanujan_prime_approx(n); break;
      }
      XSRETURN_UV(ret);
    }
    switch (ix) {
      case 0:  _vcallsub_with_pp("nth_ramanujan_prime");          break;
      case 1:  _vcallsub_with_pp("nth_ramanujan_prime_upper");    break;
      case 2:  _vcallsub_with_pp("nth_ramanujan_prime_lower");    break;
      case 3:
      default: _vcallsub_with_pp("nth_ramanujan_prime_approx");   break;
    }
    return;

void nth_twin_prime(IN SV* svn)
  ALIAS:
    nth_twin_prime_approx = 1
  PREINIT:
    UV n, ret;
  PPCODE:
    if ( _validate_and_set(&n, aTHX_ svn, IFLAG_POS) &&
         n <= MPU_MAX_TWIN_PRIME_IDX ) {
      if (n == 0) XSRETURN_UNDEF;
      switch (ix) {
        case 0:  ret = nth_twin_prime(n); break;
        case 1:
        default: ret = nth_twin_prime_approx(n); break;
      }
      XSRETURN_UV(ret);
    }
    switch (ix) {
      case 0:  _vcallsub_with_pp("nth_twin_prime");          break;
      case 1:
      default: _vcallsub_with_pp("nth_twin_prime_approx");   break;
    }
    return;

void nth_semiprime(IN SV* svn)
  ALIAS:
    nth_semiprime_approx = 1
  PREINIT:
    UV n, ret;
  PPCODE:
    if ( _validate_and_set(&n, aTHX_ svn, IFLAG_POS) &&
         n <= MPU_MAX_SEMI_PRIME_IDX ) {
      if (n == 0) XSRETURN_UNDEF;
      switch (ix) {
        case 0:  ret = nth_semiprime(n); break;
        case 1:
        default: ret = nth_semiprime_approx(n); break;
      }
      XSRETURN_UV(ret);
    }
    switch (ix) {
      case 0:  _vcallsub_with_pp("nth_semiprime");          break;
      case 1:
      default: _vcallsub_with_pp("nth_semiprime_approx");   break;
    }
    return;

void nth_lucky(IN SV* svn)
  ALIAS:
    nth_lucky_upper = 1
    nth_lucky_lower = 2
    nth_lucky_approx = 3
  PREINIT:
    UV n, ret;
  PPCODE:
    if ( _validate_and_set(&n, aTHX_ svn, IFLAG_POS) &&
         n <= MPU_MAX_LUCKY_IDX ) {
      if (n == 0) XSRETURN_UNDEF;
      switch (ix) {
        case 0:  ret = nth_lucky(n); break;
        case 1:  ret = nth_lucky_upper(n); break;
        case 2:  ret = nth_lucky_lower(n); break;
        case 3:
        default: ret = nth_lucky_approx(n); break;
      }
      XSRETURN_UV(ret);
    }
    switch (ix) {
      case 0:  _vcallsub_with_pp("nth_lucky");          break;
      case 1:  _vcallsub_with_pp("nth_lucky_upper");    break;
      case 2:  _vcallsub_with_pp("nth_lucky_lower");    break;
      case 3:
      default: _vcallsub_with_pp("nth_lucky_approx");   break;
    }
    return;


void next_prime(IN SV* svn)
  ALIAS:
    prev_prime = 1
  PREINIT:
    UV n, ret;
  PPCODE:
    if (_validate_and_set(&n, aTHX_ svn, IFLAG_POS)
        && !(ix == 0 && n >= MPU_MAX_PRIME)) {
      ret = 0;
      switch (ix) {
        case 0:  ret = next_prime(n); break;
        case 1:  ret = prev_prime(n); break;
        default: break;
      }
      if (ret == 0) XSRETURN_UNDEF;
      XSRETURN_UV(ret);
    }
    switch (ix) {
      case 0:  _vcallsub_with_gmpobj(0.01,"next_prime");   break;
      case 1:
      default: _vcallsub_with_gmpobj(0.01,"prev_prime");   break;
    }
    objectify_result(aTHX_ svn, ST(0));
    return;

void next_prime_power(IN SV* svn)
  ALIAS:
    prev_prime_power = 1
  PREINIT:
    UV n, ret;
  PPCODE:
    if (_validate_and_set(&n, aTHX_ svn, IFLAG_ABS)
        && !(ix == 0 && n >= MPU_MAX_PRIME)) {
      ret = 0;
      switch (ix) {
        case 0:  ret = next_prime_power(n); break;
        case 1:  ret = prev_prime_power(n); break;
        default: break;
      }
      if (ret == 0) XSRETURN_UNDEF;
      XSRETURN_UV(ret);
    }
    switch (ix) {
      case 0:  _vcallsub_with_pp("next_prime_power");  break;
      case 1:
      default: _vcallsub_with_pp("prev_prime_power");  break;
    }
    /* objectify_result(aTHX_ svn, ST(0)); */
    return;

void next_perfect_power(IN SV* svn)
  ALIAS:
    prev_perfect_power = 1
  PREINIT:
    UV n, ret;
    int status, power;
  PPCODE:
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
    if (status == -1) n = neg_iv(n);
    if        (status == 1 && ix == 0) {
      ret = next_perfect_power(n);
      if (ret != 0) XSRETURN_UV(ret);
    } else if (status == 1 && ix == 1) {
      if (n <= 1) XSRETURN_IV(-1);
      ret = prev_perfect_power(n);
      if (ret != 0) XSRETURN_UV(ret);
    } else if (status == -1 && ix == 0) { /* next perfect power: negative n */
      if (n == 1) XSRETURN_UV(1);
      do {
        n = prev_perfect_power(n);
        power = is_power(n,0);
      } while (n > 1 && (power <= 2 || (power & (power-1)) == 0));
      XSRETURN_IV(neg_iv(n));
    } else if (status == -1 && ix == 1) { /* prev perfect power: negative n */
      do {
        n = next_perfect_power(n);
        power = is_power(n,0);
      } while (n > 1 && (power <= 2 || (power & (power-1)) == 0));
      if (n <= (UV)IV_MAX)
        XSRETURN_IV(neg_iv(n));
    }
    switch (ix) {
      case 0:  _vcallsub_with_gmp(0.53,"next_perfect_power"); break;
      case 1:
      default: _vcallsub_with_gmp(0.53,"prev_perfect_power"); break;
    }
    objectify_result(aTHX_ svn, ST(0));
    return;

void next_chen_prime(IN SV* svn)
  PREINIT:
    UV n, ret;
  PPCODE:
    if (_validate_and_set(&n, aTHX_ svn, IFLAG_POS)) {
      ret = next_chen_prime(n);
      if (ret != 0) XSRETURN_UV(ret);
    }
    _vcallsub_with_pp("next_chen_prime");
    /* No objectify because we're calling PP, not GMP */
    return;

void urandomb(IN UV bits)
  ALIAS:
    random_ndigit_prime = 1
    random_semiprime = 2
    random_unrestricted_semiprime = 3
    random_safe_prime = 4
    random_nbit_prime = 5
    random_shawe_taylor_prime = 6
    random_maurer_prime = 7
    random_proven_prime = 8
    random_strong_prime = 9
  PREINIT:
    UV res, minarg;
    dMY_CXT;
    void* cs;
  PPCODE:
    switch (ix) {
      case 1:  minarg =   1; break;
      case 2:  minarg =   4; break;
      case 3:  minarg =   3; break;
      case 4:  minarg =   3; break;
      case 5:
      case 6:
      case 7:
      case 8:  minarg =   2; break;
      case 9:  minarg = 128; break;
      default: minarg =   0; break;
    }
    if (minarg > 0 && bits < minarg)
      croak("Parameter '%d' must be >= %d", (int)bits, (int)minarg);
    cs = MY_CXT.randcxt;
    if (bits <= BITS_PER_WORD) {
      switch (ix) {
        case 0:  res = urandomb(cs,bits); break;
        case 1:  res = random_ndigit_prime(cs,bits); break;
        case 2:  res = random_semiprime(cs,bits); break;
        case 3:  res = random_unrestricted_semiprime(cs,bits); break;
        case 4:  res = random_safe_prime(cs,bits); break;
        case 5:
        case 6:
        case 7:
        case 8:
        case 9:
        default: res = random_nbit_prime(cs,bits); break;
      }
      if (res || ix == 0) XSRETURN_UV(res);
    }
    switch (ix) {
      case 0:  _vcallsub_with_gmpobj(0.43,"urandomb"); break;
      case 1:  _vcallsub_with_gmpobj(0.42,"random_ndigit_prime"); break;
      case 2:  _vcallsub_with_gmpobj(0.00,"random_semiprime"); break;
      case 3:  _vcallsub_with_gmpobj(0.00,"random_unrestricted_semiprime"); break;
      case 4:  _vcallsub_with_gmpobj(0.52,"random_safe_prime"); break;
      case 5:  _vcallsub_with_gmpobj(0.42,"random_nbit_prime"); break;
      case 6:  _vcallsub_with_gmpobj(0.43,"random_shawe_taylor_prime"); break;
      case 7:
      case 8:  _vcallsub_with_gmpobj(0.43,"random_maurer_prime"); break;
      case 9:
      default: _vcallsub_with_gmpobj(0.43,"random_strong_prime"); break;
    }
    objectify_result(aTHX_ ST(0), ST(0));
    XSRETURN(1);

void urandomm(IN SV* svn)
  PREINIT:
    UV n, ret;
  PPCODE:
    if (_validate_and_set(&n, aTHX_ svn, IFLAG_POS)) {
      dMY_CXT;
      ret = urandomm64(MY_CXT.randcxt, n);
      XSRETURN_UV(ret);
    }
    _vcallsub_with_gmpobj(0.44,"urandomm");
    objectify_result(aTHX_ svn, ST(0));
    XSRETURN(1);

void pisano_period(IN SV* svn)
  ALIAS:
    partitions = 1
    consecutive_integer_lcm = 2
  PREINIT:
    UV n, r = 0;
  PPCODE:
    if (_validate_and_set(&n, aTHX_ svn, IFLAG_POS)) {
      switch (ix) {
        case  0: r = pisano_period(n); break;
        case  1: r = npartitions(n); break;
        case  2: r = consecutive_integer_lcm(n); break;
        default: break;
      }
      /* Returns 0 if n=0 or result overflows */
      if (r != 0 || n == 0)
        XSRETURN_UV(r);
    }
    switch (ix) {
      case 0: _vcallsub_with_gmp(0.53,"pisano_period"); break;
      case 1: _vcallsub_with_gmp(0.16,"partitions"); break;
      case 2: _vcallsub_with_gmp(0.04,"consecutive_integer_lcm"); break;
      default: break;
    }
    objectify_result(aTHX_ svn, ST(0));
    XSRETURN(1);

void random_factored_integer(IN SV* svn)
  PREINIT:
    UV n;
  PPCODE:
    if (_validate_and_set(&n, aTHX_ svn, IFLAG_POS | IFLAG_NONZERO)) {
      dMY_CXT;
      int f, nf, flip;
      UV r, F[MPU_MAX_FACTORS+1];
      AV* av = newAV();
      r = random_factored_integer(MY_CXT.randcxt, n, &nf, F);
      flip = (F[0] >= F[nf-1]);  /* Handle results in either sort order */
      for (f = 0; f < nf; f++)
        av_push(av, newSVuv(F[flip ? nf-1-f : f]));
      XPUSHs(sv_2mortal(newSVuv( r )));
      XPUSHs(sv_2mortal(newRV_noinc( (SV*) av )));
    } else {
      (void)_vcallsubn(aTHX_ G_ARRAY, VCALL_PP, "random_factored_integer",items,0);
      return;
    }

void Pi(IN UV digits = 0)
  PREINIT:
#ifdef USE_QUADMATH
#define STRTONV(t)  strtoflt128(t,NULL)
    const UV mantsize = FLT128_DIG;
    const NV pival = 3.141592653589793238462643383279502884197169Q;
#elif defined(USE_LONG_DOUBLE) && defined(HAS_LONG_DOUBLE)
#define STRTONV(t)  strtold(t,NULL)
    const UV mantsize = LDBL_DIG;
    const NV pival = 3.141592653589793238462643383279502884197169L;
#else
#define STRTONV(t)  strtod(t,NULL)
    const UV mantsize = DBL_DIG;
    const NV pival = 3.141592653589793238462643383279502884197169;
#endif
  PPCODE:
    if (digits == 0) {
      XSRETURN_NV( pival );
    } else if (digits <= mantsize) {
      char* out = pidigits(digits);
      NV pi = STRTONV(out);
      Safefree(out);
      XSRETURN_NV( pi );
    } else {
      _vcallsub_with_pp("Pi");
      return;
    }

void
_pidigits(IN int digits)
  PREINIT:
    char* out;
  PPCODE:
    if (digits <= 0) XSRETURN_EMPTY;
    out = pidigits(digits);
    XPUSHs(sv_2mortal(newSVpvn(out, digits+1)));
    Safefree(out);

void inverse_totient(IN SV* svn)
  PREINIT:
    U32 gimme_v;
    int status, it_overflow;
    UV i, n, ntotients;
  PPCODE:
    gimme_v = GIMME_V;
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_POS);
    it_overflow = (status == 1 && gimme_v == G_ARRAY && n > UV_MAX/7.5);
    if (status == 1 && !it_overflow) {
      if (gimme_v == G_SCALAR) {
        XSRETURN_UV( inverse_totient_count(n) );
      } else if (gimme_v == G_ARRAY) {
        UV* tots = inverse_totient_list(&ntotients, n);
        EXTEND(SP, (IV)ntotients);
        for (i = 0; i < ntotients; i++)
          PUSHs(sv_2mortal(newSVuv( tots[i] )));
        Safefree(tots);
      }
    } else {
      _vcallsubn(aTHX_ gimme_v, VCALL_GMP|VCALL_PP, "inverse_totient", 1, 0);
      return; /* skip implicit PUTBACK */
    }

void
factor(IN SV* svn)
  ALIAS:
    factor_exp = 1
  PREINIT:
    U32 gimme_v;
    int status, i, nfactors;
    UV n;
  PPCODE:
    gimme_v = GIMME_V;
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_POS);
    if (status == 1) {
      UV factors[MPU_MAX_FACTORS+1];
      UV exponents[MPU_MAX_FACTORS+1];
      if (gimme_v == G_SCALAR) {
        UV res = (ix == 0)  ?  factor(n, factors)  :  factor_exp(n, factors, 0);
        XSRETURN_UV(res);
      } else if (gimme_v == G_ARRAY) {
        if (ix == 0) {
          nfactors = factor(n, factors);
          EXTEND(SP, nfactors);
          for (i = 0; i < nfactors; i++)
            PUSHs(sv_2mortal(newSVuv( factors[i] )));
        } else {
          nfactors = factor_exp(n, factors, exponents);
          /* if (n == 1)  XSRETURN_EMPTY; */
          EXTEND(SP, nfactors);
          for (i = 0; i < nfactors; i++) {
            AV* av = newAV();
            av_push(av, newSVuv(factors[i]));
            av_push(av, newSVuv(exponents[i]));
            PUSHs( sv_2mortal(newRV_noinc( (SV*) av )) );
          }
        }
      }
    } else {
      if (ix == 0)
        _vcallsubn(aTHX_ gimme_v, VCALL_ROOT, "_generic_factor", 1, 0);
      else
        _vcallsubn(aTHX_ gimme_v, VCALL_ROOT, "_generic_factor_exp", 1, 0);
      return; /* skip implicit PUTBACK */
    }

void divisors(IN SV* svn, IN SV* svk = 0)
  PREINIT:
    U32 gimme_v;
    int status;
    UV n, k, i, ndivisors, *divs;
  PPCODE:
    gimme_v = GIMME_V;
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_POS);
    k = n;
    if (status == 1 && svk != 0) {
      status = _validate_and_set(&k, aTHX_ svk, IFLAG_POS);
      if (k > n)  k = n;
    }
    if (status != 1) {
      _vcallsubn(aTHX_ gimme_v, VCALL_GMP|VCALL_PP, "divisors", items, 53);
      return; /* skip implicit PUTBACK */
    }
    if (GIMME_V == G_VOID) {
      /* Nothing */
    } else if (GIMME_V == G_SCALAR && k >= n) {
      ndivisors = divisor_sum(n, 0);
      PUSHs(sv_2mortal(newSVuv( ndivisors )));
    } else {
      divs = divisor_list(n, &ndivisors, k);
      if (GIMME_V == G_SCALAR) {
        PUSHs(sv_2mortal(newSVuv( ndivisors )));
      } else {
        EXTEND(SP, (IV)ndivisors);
        for (i = 0; i < ndivisors; i++)
          PUSHs(sv_2mortal(newSVuv( divs[i] )));
      }
      Safefree(divs);
    }

void
trial_factor(IN UV n, ...)
  ALIAS:
    fermat_factor = 1
    holf_factor = 2
    squfof_factor = 3
    lehman_factor = 4
    prho_factor = 5
    cheb_factor = 6
    pplus1_factor = 7
    pbrent_factor = 8
    pminus1_factor = 9
    ecm_factor = 10
  PREINIT:
    UV arg1, arg2;
    static const UV default_arg1[] =
       {0,     64000000, 8000000, 4000000, 1,   4000000, 0,    200, 4000000, 1000000};
     /* Trial, Fermat,   Holf,    SQUFOF,  Lmn, PRHO,    Cheb, P+1, Brent,    P-1 */
  PPCODE:
    if (n == 0)  XSRETURN_UV(0);
    if (ix == 10) {  /* We don't have an ecm_factor, call PP. */
      _vcallsubn(aTHX_ GIMME_V, VCALL_PP, "ecm_factor", 1, 0);
      return;
    }
    /* Must read arguments before pushing anything */
    arg1 = (items >= 2) ? my_svuv(ST(1)) : default_arg1[ix];
    arg2 = (items >= 3) ? my_svuv(ST(2)) : 0;
    /* Small factors */
    while ( (n% 2) == 0 ) {  n /=  2;  XPUSHs(sv_2mortal(newSVuv( 2 ))); }
    while ( (n% 3) == 0 ) {  n /=  3;  XPUSHs(sv_2mortal(newSVuv( 3 ))); }
    while ( (n% 5) == 0 ) {  n /=  5;  XPUSHs(sv_2mortal(newSVuv( 5 ))); }
    if (n == 1) {  /* done */ }
    else if (is_prime(n)) { XPUSHs(sv_2mortal(newSVuv( n ))); }
    else {
      UV factors[MPU_MAX_FACTORS+1];
      int i, nfactors = 0;
      switch (ix) {
        case 0:  nfactors = trial_factor  (n, factors, 2, arg1);  break;
        case 1:  nfactors = fermat_factor (n, factors, arg1);  break;
        case 2:  nfactors = holf_factor   (n, factors, arg1);  break;
        case 3:  nfactors = squfof_factor (n, factors, arg1);  break;
        case 4:  nfactors = lehman_factor (n, factors, arg1);  break;
        case 5:  nfactors = prho_factor   (n, factors, arg1);  break;
        case 6:  nfactors = cheb_factor   (n, factors, arg1, arg2);  break;
        case 7:  nfactors = pplus1_factor (n, factors, arg1);  break;
        case 8:  if (items < 3) arg2 = 1;
                 nfactors = pbrent_factor (n, factors, arg1, arg2);  break;
        case 9:
        default: if (items < 3) arg2 = 10*arg1;
                 nfactors = pminus1_factor(n, factors, arg1, arg2);  break;
      }
      EXTEND(SP, nfactors);
      for (i = 0; i < nfactors; i++)
        PUSHs(sv_2mortal(newSVuv( factors[i] )));
    }


void
divisor_sum(IN SV* svn, ...)
  PREINIT:
    UV n, k, sigma;
  PPCODE:
    if (items == 1 && _validate_and_set(&n, aTHX_ svn, IFLAG_POS)) {
      sigma = divisor_sum(n, 1);
      if (n <= 1 || sigma != 0)
        XSRETURN_UV(sigma);
    } else {
      SV* svk = ST(1);
      if ( (!SvROK(svk) || (SvROK(svk) && SvTYPE(SvRV(svk)) != SVt_PVCV)) &&
           _validate_and_set(&n, aTHX_ svn, IFLAG_POS) &&
           _validate_and_set(&k, aTHX_ svk, IFLAG_POS) ) {
        sigma = divisor_sum(n, k);
        if (n <= 1 || sigma != 0)
          XSRETURN_UV(sigma);
      }
    }
    _vcallsub_with_pp("divisor_sum");
    return; /* skip implicit PUTBACK */

void
jordan_totient(IN SV* sva, IN SV* svn)
  ALIAS:
    powersum = 1
    ramanujan_sum = 2
    legendre_phi = 3
    smooth_count = 4
    rough_count = 5
  PREINIT:
    int astatus, nstatus;
    UV a, n, ret;
  PPCODE:
    astatus = _validate_and_set(&a, aTHX_ sva, IFLAG_POS);
    nstatus = _validate_and_set(&n, aTHX_ svn, IFLAG_POS);
    if (astatus != 0 && nstatus != 0) {
      switch (ix) {
        case 0:  ret = jordan_totient(a, n);
                 if (ret == 0 && n > 1)
                   goto overflow;
                 break;
        case 1:  ret = powersum(a, n);
                 if (ret == 0 && a > 0)
                   goto overflow;
                 break;
        case 2:  if (a < 1 || n < 1) XSRETURN_IV(0);
                 {
                   UV g = a / gcd_ui(a,n);
                   int m = moebius(g);
                   if (m == 0 || a == g) RETURN_NPARITY(m);
                   XSRETURN_IV( m * (totient(a) / totient(g)) );
                 }
                 break;
        case 3:  ret = legendre_phi(a, n); break;
        case 4:  ret = debruijn_psi(a, n); break;
        case 5:
        default: ret = buchstab_phi(a, n); break;
      }
      XSRETURN_UV(ret);
    }
    overflow:
    switch (ix) {
      case 0:  _vcallsub_with_gmp(0.22,"jordan_totient");  break;
      case 1:  _vcallsub_with_gmp(0.53,"powersum");  break;
      case 2:  _vcallsub_with_pp("ramanujan_sum");  break;
      case 3:  _vcallsub_with_pp("legendre_phi");  break;
      case 4:  _vcallsub_with_pp("smooth_count"); break;
      case 5:
      default: _vcallsub_with_pp("rough_count"); break;
    }
    objectify_result(aTHX_ sva, ST(0));
    return; /* skip implicit PUTBACK */

void almost_prime_count(IN SV* svk, IN SV* svn)
  ALIAS:
    almost_prime_count_approx = 1
    almost_prime_count_lower = 2
    almost_prime_count_upper = 3
    omega_prime_count = 4
  PREINIT:
    UV k, n, ret;
  PPCODE:
    if (_validate_and_set(&k, aTHX_ svk, IFLAG_ABS) &&
        _validate_and_set(&n, aTHX_ svn, IFLAG_ABS) &&
        k < BITS_PER_WORD) {
      ret = 0;
      switch (ix) {
        case 0:  ret = almost_prime_count(k, n); break;
        case 1:  ret = almost_prime_count_approx(k, n); break;
        case 2:  ret = almost_prime_count_lower(k, n); break;
        case 3:  ret = almost_prime_count_upper(k, n); break;
        case 4:  ret = omega_prime_count(k, n); break;
        default: break;
      }
      XSRETURN_UV(ret);
    }
    switch (ix) {
      case 0:  _vcallsub_with_pp("almost_prime_count");  break;
      case 1:  _vcallsub_with_pp("almost_prime_count_approx");  break;
      case 2:  _vcallsub_with_pp("almost_prime_count_lower");  break;
      case 3:  _vcallsub_with_pp("almost_prime_count_upper");  break;
      case 4:
      default: _vcallsub_with_pp("omega_prime_count");  break;
    }
    objectify_result(aTHX_ svn, ST(0));
    return;

void nth_almost_prime(IN SV* svk, IN SV* svn)
  ALIAS:
    nth_almost_prime_approx = 1
    nth_almost_prime_lower = 2
    nth_almost_prime_upper = 3
  PREINIT:
    UV k, n, max;
  PPCODE:
    if (_validate_and_set(&k, aTHX_ svk, IFLAG_ABS) &&
        _validate_and_set(&n, aTHX_ svn, IFLAG_ABS) &&
        k < BITS_PER_WORD) {
      UV ret = 0;
      if (n == 0 || (k == 0 && n > 1)) XSRETURN_UNDEF;
      max = max_almost_prime_count(k);
      if (max > 0  &&  n <= max) {
        switch (ix) {
          case 0: ret = nth_almost_prime(k, n); break;
          case 1: ret = nth_almost_prime_approx(k, n); break;
          case 2: ret = nth_almost_prime_lower(k, n); break;
          case 3: ret = nth_almost_prime_upper(k, n); break;
        }
        XSRETURN_UV(ret);
      }
    }
    switch (ix) {
      case 0:  _vcallsub_with_pp("nth_almost_prime");  break;
      case 1:  _vcallsub_with_pp("nth_almost_prime_approx");  break;
      case 2:  _vcallsub_with_pp("nth_almost_prime_lower");  break;
      case 3:
      default: _vcallsub_with_pp("nth_almost_prime_upper");  break;
    }
    objectify_result(aTHX_ svn, ST(0));
    return;

void nth_omega_prime(IN SV* svk, IN SV* svn)
  PREINIT:
    UV k, n, max, ret;
  PPCODE:
    if (_validate_and_set(&k, aTHX_ svk, IFLAG_ABS) &&
        _validate_and_set(&n, aTHX_ svn, IFLAG_ABS) &&
        k < 16) {
      if (n == 0 || (k == 0 && n > 1)) XSRETURN_UNDEF;
      max = max_omega_prime_count(k);
      if (max > 0  &&  n <= max) {
        ret = nth_omega_prime(k, n);
        XSRETURN_UV(ret);
      }
    }
    _vcallsub_with_pp("nth_omega_prime");
    objectify_result(aTHX_ svn, ST(0));
    return;


void powmod(IN SV* sva, IN SV* svg, IN SV* svn)
  ALIAS:
    rootmod = 1
  PREINIT:
    int astatus, gstatus, nstatus, retundef;
    UV a, g, n, ret;
  PPCODE:
    astatus = _validate_and_set(&a, aTHX_ sva, IFLAG_ANY);
    gstatus = _validate_and_set(&g, aTHX_ svg, IFLAG_ANY);
    nstatus = _validate_and_set(&n, aTHX_ svn, IFLAG_ABS);
    if (astatus != 0 && gstatus != 0 && nstatus != 0) {
      if (n == 0) XSRETURN_UNDEF;
      if (n == 1) XSRETURN_UV(0);
      _mod_with(&a, astatus, n);
      retundef = ret = 0;
      if (ix == 0) {
        retundef = !prep_pow_inv(&a,&g,gstatus,n);
        if (!retundef) ret = powmod(a, g, n);
      } else {
        retundef = !(prep_pow_inv(&a,&g,gstatus,n) && rootmod(&ret,a,g,n));
      }
      if (retundef) XSRETURN_UNDEF;
      XSRETURN_UV(ret);
    }
    if (ix == 0) _vcallsub_with_gmpobj(0.36,"powmod");
    else         _vcallsub_with_gmpobj(0.00,"rootmod");
    objectify_result(aTHX_ svn, ST(0));
    return;

void addmod(IN SV* sva, IN SV* svb, IN SV* svn)
  ALIAS:
    submod = 1
    mulmod = 2
    divmod = 3
    znlog = 4
  PREINIT:
    int astatus, bstatus, nstatus, retundef;
    UV a, b, n, ret;
  PPCODE:
    astatus = _validate_and_set(&a, aTHX_ sva, IFLAG_ANY);
    bstatus = _validate_and_set(&b, aTHX_ svb, IFLAG_ANY);
    nstatus = _validate_and_set(&n, aTHX_ svn, IFLAG_ABS);
    if (astatus != 0 && bstatus != 0 && nstatus != 0) {
      if (n == 0) XSRETURN_UNDEF;
      if (n == 1) XSRETURN_UV(0);
      _mod_with(&a, astatus, n);
      _mod_with(&b, bstatus, n);
      retundef = ret = 0;
      switch (ix) {
        case 0:  ret = addmod(a, b, n); break;
        case 1:  ret = submod(a, b, n); break;
        case 2:  ret = mulmod(a, b, n); break;
        case 3:  b = modinverse(b, n);
                 if (b == 0) retundef = 1;
                 else        ret = mulmod(a, b, n);
                 break;
        case 4:  ret = znlog(a, b, n);
                 if (ret == 0 && (b == 0 || a != 1))  retundef = 1;
        default: break;
      }
      if (retundef) XSRETURN_UNDEF;
      XSRETURN_UV(ret);
    }
    switch (ix) {
      case 0: _vcallsub_with_gmpobj(0.36,"addmod"); break;
      case 1: _vcallsub_with_gmpobj(0.53,"submod"); break;
      case 2: _vcallsub_with_gmpobj(0.36,"mulmod"); break;
      case 3: _vcallsub_with_gmpobj(0.36,"divmod"); break;
      case 4: _vcallsub_with_gmpobj(0.00,"znlog"); break;
      default:break;
    }
    objectify_result(aTHX_ svn, ST(0));
    return;

void muladdmod(IN SV* sva, IN SV* svb, IN SV* svc, IN SV* svn)
  ALIAS:
    mulsubmod = 1
  PREINIT:
    int astatus, bstatus, cstatus, nstatus;
    UV a, b, c, n, ret;
  PPCODE:
    astatus = _validate_and_set(&a, aTHX_ sva, IFLAG_ANY);
    bstatus = _validate_and_set(&b, aTHX_ svb, IFLAG_ANY);
    cstatus = _validate_and_set(&c, aTHX_ svc, IFLAG_ANY);
    nstatus = _validate_and_set(&n, aTHX_ svn, IFLAG_ABS);
    if (astatus != 0 && bstatus != 0 && cstatus != 0 && nstatus != 0) {
      if (n == 0) XSRETURN_UNDEF;
      if (n == 1) XSRETURN_UV(0);
      _mod_with(&a, astatus, n);
      _mod_with(&b, bstatus, n);
      _mod_with(&c, cstatus, n);
      ret = (ix==0)  ?  muladdmod(a,b,c,n)  :  mulsubmod(a,b,c,n);
      XSRETURN_UV(ret);
    }
    switch (ix) {
      case 0: _vcallsub_with_gmpobj(0.53,"muladdmod"); break;
      case 1: _vcallsub_with_gmpobj(0.53,"mulsubmod"); break;
      default:break;
    }
    objectify_result(aTHX_ svn, ST(0));
    return;

void binomialmod(IN SV* svn, IN SV* svk, IN SV* svm)
  PREINIT:
    int nstatus, kstatus, mstatus;
    UV ret, n, k, m;
  PPCODE:
    nstatus = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
    kstatus = _validate_and_set(&k, aTHX_ svk, IFLAG_ANY);
    mstatus = _validate_and_set(&m, aTHX_ svm, IFLAG_ABS);
    if (nstatus != 0 && kstatus != 0 && mstatus != 0) {
      if (m == 0) XSRETURN_UNDEF;
      if (m == 1) XSRETURN_UV(0);
      if ( (nstatus == 1 && (kstatus == -1 || k > n)) ||
           (nstatus ==-1 && (kstatus == -1 && k > n)) )
         XSRETURN_UV(0);
      if (kstatus == -1) k = n - k;
      if (nstatus == -1) n = neg_iv(n) + k - 1;
      if (binomialmod(&ret, n, k, m)) {
        if ((nstatus == -1) && (k & 1)) ret = (ret==0) ? 0 : m-ret;
        XSRETURN_UV(ret);
      }
    }
    _vcallsub_with_gmpobj(0.00,"binomialmod");
    objectify_result(aTHX_ svm, ST(0));
    return;

void factorialmod(IN SV* sva, IN SV* svn)
  PREINIT:
    int astatus, nstatus;
    UV a, n;
  PPCODE:
    astatus = _validate_and_set(&a, aTHX_ sva, IFLAG_POS);
    nstatus = _validate_and_set(&n, aTHX_ svn, IFLAG_ABS);
    if (astatus != 0 && nstatus != 0) {
      if (n == 0) XSRETURN_UNDEF;
      if (n == 1) XSRETURN_UV(0);
      XSRETURN_UV( factorialmod(a, n) );
    }
    _vcallsub_with_gmp(0.47,"factorialmod");
    objectify_result(aTHX_ svn, ST(0));
    return;

void invmod(IN SV* sva, IN SV* svn)
  ALIAS:
    znorder = 1
    sqrtmod = 2
    negmod = 3
  PREINIT:
    int astatus, nstatus;
    UV a, n, r, retok;
  PPCODE:
    astatus = _validate_and_set(&a, aTHX_ sva, IFLAG_ANY);
    nstatus = _validate_and_set(&n, aTHX_ svn, IFLAG_ABS);
    if (astatus != 0 && nstatus != 0) {
      if (n == 0) XSRETURN_UNDEF;
      if (n == 1) XSRETURN_UV((ix==1) ? 1 : 0); /* znorder different */
      _mod_with(&a, astatus, n);
      retok = 1;
      switch (ix) {
        case 0:  retok = r = modinverse(a, n); break;
        case 1:  retok = r = znorder(a, n);    break;
        case 2:  retok = sqrtmod(&r, a, n);    break;
        case 3:
        default: r = negmod(a, n); break;
      }
      if (retok == 0) XSRETURN_UNDEF;
      XSRETURN_UV(r);
    }
    /* GMP sqrtmod didn't work with composites until 0.53. */
    switch (ix) {
      case 0:  _vcallsub_with_gmp(0.20,"invmod"); break;
      case 1:  _vcallsub_with_gmp(0.22,"znorder"); break;
      case 2:  _vcallsub_with_gmp(0.53,"sqrtmod"); break;
      case 3:
      default: _vcallsub_with_gmp(0.00,"negmod");  break;
    }
    objectify_result(aTHX_ ST(0), ST(0));
    return;

void allsqrtmod(IN SV* sva, IN SV* svn)
  PREINIT:
    int astatus, nstatus;
    UV a, n, i, numr, *roots;
  PPCODE:
    astatus = _validate_and_set(&a, aTHX_ sva, IFLAG_ANY);
    nstatus = _validate_and_set(&n, aTHX_ svn, IFLAG_ABS);
    if (astatus != 0 && nstatus != 0) {
      if (n == 0) XSRETURN_EMPTY;
      _mod_with(&a, astatus, n);
      roots = allsqrtmod(&numr, a, n);
      if (roots != 0) {
        EXTEND(SP, (IV)numr);
        for (i = 0; i < numr; i++)
          PUSHs(sv_2mortal(newSVuv(roots[i])));
        Safefree(roots);
      }
    } else {
      (void)_vcallsubn(aTHX_ GIMME_V, VCALL_PP, "allsqrtmod", items, 0);
      return;
    }

void allrootmod(IN SV* sva, IN SV* svg, IN SV* svn)
  PREINIT:
    int astatus, gstatus, nstatus;
    UV a, g, n, i, numr, *roots;
  PPCODE:
    astatus = _validate_and_set(&a, aTHX_ sva, IFLAG_ANY);
    gstatus = _validate_and_set(&g, aTHX_ svg, IFLAG_ANY);
    nstatus = _validate_and_set(&n, aTHX_ svn, IFLAG_ABS);
    if (astatus != 0 && gstatus != 0 && nstatus != 0) {
      if (n == 0) XSRETURN_EMPTY;
      _mod_with(&a, astatus, n);
      if (!prep_pow_inv(&a,&g,gstatus,n)) XSRETURN_EMPTY;
      roots = allrootmod(&numr, a, g, n);
      if (roots != 0) {
        EXTEND(SP, (IV)numr);
        for (i = 0; i < numr; i++)
          PUSHs(sv_2mortal(newSVuv(roots[i])));
        Safefree(roots);
      }
    } else {
      (void)_vcallsubn(aTHX_ GIMME_V, VCALL_PP, "allrootmod", items, 0);
      return;
    }

void is_primitive_root(IN SV* sva, IN SV* svn)
  PREINIT:
    int astatus, nstatus;
    UV a, n;
  PPCODE:
    astatus = _validate_and_set(&a, aTHX_ sva, IFLAG_ANY);
    nstatus = _validate_and_set(&n, aTHX_ svn, IFLAG_ABS);
    if (astatus != 0 && nstatus != 0) {
      if (n == 0) XSRETURN_UNDEF;
      _mod_with(&a, astatus, n);
      RETURN_NPARITY( is_primitive_root(a,n,0) );
    }
    _vcallsub_with_gmp(0.36,"is_primitive_root");
    return;

void qnr(IN SV* svn)
  ALIAS:
    znprimroot = 1
  PREINIT:
    UV n, r;
  PPCODE:
    if (_validate_and_set(&n, aTHX_ svn, IFLAG_ABS)) {
      if (n == 0) XSRETURN_UNDEF;
      if (ix == 0) {
        r = qnr(n);
      } else {
        r = znprimroot(n);
        if (r == 0 && n != 1)  XSRETURN_UNDEF;
      }
      if (r < 100)  RETURN_NPARITY(r);
      else          XSRETURN_UV(r);
    }
    if (ix == 0)  _vcallsub_with_gmp(0.00,"qnr");
    else          _vcallsub_with_gmp(0.22,"znprimroot");
    objectify_result(aTHX_ svn, ST(0));
    return;

void
is_smooth(IN SV* svn, IN SV* svk)
  ALIAS:
    is_rough = 1
    is_omega_prime = 2
    is_almost_prime = 3
  PREINIT:
    UV n, k;
  PPCODE:
    if (_validate_and_set(&n, aTHX_ svn, IFLAG_POS) &&
        _validate_and_set(&k, aTHX_ svk, IFLAG_POS)) {
      int res;
      switch (ix) {
        case 0:  res = is_smooth(n,k); break;
        case 1:  res = is_rough(n,k); break;
        case 2:  res = is_omega_prime(n,k); break; /* Note order */
        case 3:
        default: res = is_almost_prime(n,k); break; /* Note order */
      }
      RETURN_NPARITY(res);
    }
    switch (ix) {
      case 0:  _vcallsub_with_gmp(0.53,"is_smooth");  break;
      case 1:  _vcallsub_with_gmp(0.53,"is_rough");  break;
      case 2:  _vcallsub_with_gmp(0.00,"is_omega_prime");  break;
      case 3:
      default: _vcallsub_with_gmp(0.53,"is_almost_prime"); break;
    }
    return;

void is_divisible(IN SV* svn, IN SV* svd)
  PREINIT:
    UV n, d;
  PPCODE:
    if (_validate_and_set(&n, aTHX_ svn, IFLAG_ABS) &&
        _validate_and_set(&d, aTHX_ svd, IFLAG_ABS)) {
      if (d == 0) RETURN_NPARITY(n == 0);
      RETURN_NPARITY( (n % d) == 0 );
    }
    _vcallsub_with_gmp(0.53,"is_divisible");
    return;

void is_congruent(IN SV* svn, IN SV* svc, IN SV* svd)
  PREINIT:
    UV n, c, d;
    int nstatus, cstatus, dstatus;
  PPCODE:
    nstatus = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
    cstatus = _validate_and_set(&c, aTHX_ svc, IFLAG_ANY);
    dstatus = _validate_and_set(&d, aTHX_ svd, IFLAG_ABS);
    if (nstatus != 0 && cstatus != 0 && dstatus != 0) {
      if (d != 0) {
        _mod_with(&n, nstatus, d);
        _mod_with(&c, cstatus, d);
      }
      RETURN_NPARITY( n == c );
    }
    _vcallsub_with_gmp(0.53,"is_congruent");
    return;

void valuation(IN SV* svn, IN SV* svk)
  PREINIT:
    UV n, k;
  PPCODE:
    if (_validate_and_set(&n, aTHX_ svn, IFLAG_ABS) &&
        _validate_and_set(&k, aTHX_ svk, IFLAG_POS)) {
      if (k <= 1)  croak("valuation: k must be > 1");
      if (n == 0) XSRETURN_UNDEF;
      RETURN_NPARITY(valuation(n, k));
    }
    _vcallsub_with_gmp(0.20,"valuation");
    return;

void is_powerful(IN SV* svn, IN SV* svk = 0);
  ALIAS:
    powerful_count = 1
    sumpowerful = 2
    nth_powerful = 3
  PREINIT:
    int nstatus;
    UV n, ret, k = 2;
  PPCODE:
    nstatus = _validate_and_set(&n, aTHX_ svn, (ix < 3) ? IFLAG_ANY: IFLAG_POS);
    if (nstatus != 0 && (!svk || _validate_and_set(&k, aTHX_ svk, IFLAG_POS))) {
      if (nstatus == -1) RETURN_NPARITY(0);
      if (ix == 0) RETURN_NPARITY( is_powerful(n, k) );
      if (ix == 1) XSRETURN_UV( powerful_count(n, k) );
      if (ix == 2) {
        if (n == 0) XSRETURN_UV(0);
        ret = sumpowerful(n, k);
      } else {
        if (n == 0) XSRETURN_UNDEF;
        ret = nth_powerful(n, k);
      }
      /* ret=0: nth_powerful / sumpowerful result > UV_MAX, so go to PP/GMP */
      if (ret > 0) XSRETURN_UV(ret);
    }
    switch (ix) {
      case 0: _vcallsub_with_gmp(0.53, "is_powerful"); break;
      case 1: _vcallsub_with_gmp(0.53, "powerful_count"); break;
      case 2: _vcallsub_with_gmp(0.00, "sumpowerful"); break;
      case 3: _vcallsub_with_gmp(0.00, "nth_powerful"); break;
      default: break;
    }
    return;


void kronecker(IN SV* sva, IN SV* svb)
  PREINIT:
    int astatus, bstatus;
    UV a, b;
  PPCODE:
    astatus = _validate_and_set(&a, aTHX_ sva, IFLAG_ANY);
    bstatus = _validate_and_set(&b, aTHX_ svb, IFLAG_ANY);
    if (astatus != 0 && bstatus != 0) {
      int k;
      if (bstatus == 1)
        k = (astatus==1) ? kronecker_uu(a,b)         :  kronecker_su((IV)a,b);
      else
        k = (astatus==1) ? kronecker_uu(a,neg_iv(b)) : -kronecker_su((IV)a,neg_iv(b));
      RETURN_NPARITY( k );
    }
    _vcallsub_with_gmp(0.17,"kronecker");
    return;

void is_qr(IN SV* sva, IN SV* svn)
  PREINIT:
    int astatus, nstatus;
    UV a, n;
  PPCODE:
    astatus = _validate_and_set(&a, aTHX_ sva, IFLAG_ANY);
    nstatus = _validate_and_set(&n, aTHX_ svn, IFLAG_ABS);
    if (astatus != 0 && nstatus != 0) {
      if (n == 0) XSRETURN_UNDEF;
      if (n == 1) RETURN_NPARITY(1);
      _mod_with(&a, astatus, n);
      RETURN_NPARITY( is_qr(a,n) );
    }
    _vcallsub_with_gmp(0.53,"is_qr");
    return;

void addint(IN SV* sva, IN SV* svb)
  ALIAS:
    subint = 1
    mulint = 2
    divint = 3
    modint = 4
    cdivint = 5
    powint = 6
  PREINIT:
    int astatus, bstatus, overflow, postneg, nix, smask;
    UV a, b, t, ret;
  PPCODE:
    astatus = _validate_and_set(&a, aTHX_ sva, IFLAG_ANY);
    bstatus = _validate_and_set(&b, aTHX_ svb, (ix == 6) ? IFLAG_POS : IFLAG_ANY);

    if (astatus != 0 && bstatus != 0) {
      /* We will try to do everything with non-negative integers, with overflow
       * detection.  This means some pre-processing and post-processing for
       * negative inputs. */
      nix = ix;  /* So we can modify */
      ret = overflow = postneg = 0;
      smask = ((astatus == -1) << 1) + (bstatus == -1);
      /* smask=0: +a +b  smask=1: +a -b  smask=2: -a +b  smask=3: -a -b */

      if (ix == 0 && smask != 0) {
        switch (smask) {
          case 1: nix=1; b = neg_iv(b); break;                /* a - |b| */
          case 2: nix=1; t = neg_iv(a); a = b; b = t; break;  /* b - |a| */
          case 3: a=neg_iv(a); b=neg_iv(b); postneg=1; break; /* -(|a| + |b|) */
          default: break;
        }
        astatus = bstatus = 1;
      }
      if (ix == 1 && smask != 0) {
        switch (smask) {
          case 1: nix=0; b = neg_iv(b); break;                 /* a + |b| */
          case 2: nix=0; a = neg_iv(a); postneg=1; break;      /* -(|a| + b) */
          case 3: t=neg_iv(a); a=neg_iv(b); b = t; break;      /* |b| - |a| */
          default: break;
        }
        astatus = bstatus = 1;
      }
      if (ix == 2 && smask != 0) {
        switch (smask) {
          case 1: b = neg_iv(b); postneg = 1; break;
          case 2: a = neg_iv(a); postneg = 1; break;
          case 3: a = neg_iv(a); b = neg_iv(b); break;
          default: break;
        }
        astatus = bstatus = 1;
      }
      if (ix == 3 || ix == 4) {
        if (ix == 3 && b == 0) croak("divint: divide by zero");
        if (ix == 4 && b == 0) croak("modint: divide by zero");
        if (smask != 0 && (astatus == -1 || a <= (UV)IV_MAX) && (bstatus == -1 || b <= (UV)IV_MAX)) {
          IV q, r;
          (void) fdivrem(&q, &r, (IV)a, (IV)b);
          XSRETURN_IV( (ix == 3) ? q : r );
        }
      }
      if (ix == 5) {
        if (b == 0) croak("cdivint: divide by zero");
        if (smask != 0 && (astatus == -1 || a <= (UV)IV_MAX) && (bstatus == -1 || b <= (UV)IV_MAX)) {
          IV q, r;
          (void) cdivrem(&q, &r, (IV)a, (IV)b);
          XSRETURN_IV( q );
        }
      }
      if (ix == 6 && astatus != 1) {  /* bstatus is never -1 for powint */
        a = neg_iv(a);
        postneg = (b & 1);
        astatus = 1;
      }
      if (astatus == 1 && bstatus == 1) {
        switch (nix) {
          case 0:  ret = a + b;                  /* addint */
                   overflow = UV_MAX-a < b;
                   break;
          case 1:  ret = a - b;                  /* subint */
                   if (b > a && (IV)ret < 0) XSRETURN_IV((IV)ret);
                   overflow = (b > a);
                   break;
          case 2:  ret = a * b;                  /* mulint */
                   overflow = a > 0 && UV_MAX/a < b;
                   break;
          case 3:  ret = a / b; break;           /* divint */
          case 4:  ret = a % b; break;           /* modint */
          case 5:  ret = a / b + (a % b != 0);   /* cdivint */
                   break;
          case 6:
          default: ret = ipowsafe(a, b);
                   overflow = (a > 1 && ret == UV_MAX);
                   break;
        }
        if (!overflow) {
          if (!postneg)
            XSRETURN_UV(ret);
          if (ret <= (UV)IV_MAX)
            XSRETURN_IV(neg_iv(ret));
        }
      }
    }
    switch (ix) {
      case 0:  _vcallsub_with_gmp(0.52,"addint"); break;
      case 1:  _vcallsub_with_gmp(0.52,"subint"); break;
      case 2:  _vcallsub_with_gmp(0.52,"mulint"); break;
      case 3:  _vcallsub_with_gmp(0.52,"divint"); break;
      case 4:  _vcallsub_with_gmp(0.52,"modint"); break;
      case 5:  _vcallsub_with_gmp(0.53,"cdivint"); break;
      case 6:
      default: _vcallsub_with_gmp(0.52,"powint"); break;
    }
    objectify_result(aTHX_ ST(0), ST(0));
    return; /* skip implicit PUTBACK */

void add1int(IN SV* svn)
  ALIAS:
    sub1int = 1
  PREINIT:
    int status;
    UV n;
  PPCODE:
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
    if (status == 1) {
      if (ix == 1 && n == 0)  XSRETURN_IV(-1);
      if (ix == 1 || (ix == 0 && n < UV_MAX))
        XSRETURN_UV( (ix==0) ? n+1 : n-1 );
    } else if (status == -1) {
      if (ix == 0 || (ix == 1 && (IV)n > IV_MIN))
        XSRETURN_IV( (ix==0) ? (IV)n+1 : (IV)n-1 );
    }
    _vcallsub_with_gmp(0.53, (ix == 0) ? "add1int" : "sub1int");
    objectify_result(aTHX_ svn, ST(0));
    return;

void absint(IN SV* svn)
  ALIAS:
    negint = 1
  PREINIT:
    UV n;
  PPCODE:
    if (ix == 0) {
      if (_validate_and_set(&n, aTHX_ svn, IFLAG_ABS))
        XSRETURN_UV(n);
      _vcallsub_with_gmp(0.52,"absint");
    } else {
      int status = _validate_and_set(&n, aTHX_ svn, IFLAG_IV);
      if      (status == -1) XSRETURN_UV(neg_iv(n));
      else if (status ==  1) XSRETURN_IV(neg_iv(n));
      _vcallsub_with_gmp(0.52,"negint");
    }
    objectify_result(aTHX_ svn, ST(0));
    return;

void signint(IN SV* svn)
  ALIAS:
    is_odd = 1
    is_even = 2
  PREINIT:
    int status, sign, isodd;
    UV n;
    const char* s;
    STRLEN len;
  PPCODE:
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
    if (status == 0) {  /* Look at the string input */
      s = SvPV(svn, len);
      if (len == 0 || s == 0) croak("%s: invalid parameter", (ix==0) ? "signint" : (ix==1) ? "is_odd" : "is_even");
      sign = (s[0] == '-')  ?  -1  : (s[0] == '0')  ?  0  :  1;
      isodd = (s[len-1] == '1' || s[len-1] == '3' || s[len-1] == '5' || s[len-1] == '7' || s[len-1] == '9');
    } else {
      sign = (status == -1)  ?  -1  :  (n == 0)  ?  0  :  1;
      isodd = n & 1;
    }
    RETURN_NPARITY( (ix==0) ? sign : (ix==1) ? isodd : !isodd );

void cmpint(IN SV* sva, IN SV* svb)
  PREINIT:
    int astatus, bstatus, ret = 0;
    UV a, b;
  PPCODE:
    astatus = _validate_and_set(&a, aTHX_ sva, IFLAG_ANY);
    bstatus = _validate_and_set(&b, aTHX_ svb, IFLAG_ANY);
    if (astatus != 0 && bstatus != 0) {
      if      (astatus > bstatus) ret = 1;
      else if (astatus < bstatus) ret = -1;
      else if (a == b)            ret = 0;
      else                        ret = ((astatus == 1 && a > b) || (astatus == -1 && (IV)a > (IV)b)) ? 1 : -1;
    } else {
      STRLEN alen, blen;
      char *aptr, *bptr;
      aptr = SvPV(sva, alen);
      bptr = SvPV(svb, blen);
      ret = strnum_cmp(aptr, alen, bptr, blen);
    }
    RETURN_NPARITY(ret);

void logint(IN SV* svn, IN UV k, IN SV* svret = 0)
  ALIAS:
    rootint = 1
  PREINIT:
    UV n, root;
  PPCODE:
    if (ix == 0 && k <= 1)  croak("logint: base must be > 1");
    if (ix == 1 && k <= 0)  croak("rootint: k must be > 0");
    if (svret != 0 && !SvROK(svret))
      croak("%s: third argument not a scalar reference",(ix==0)?"logint":"rootint");
    if (_validate_and_set(&n, aTHX_ svn, ix == 0 ? IFLAG_POS | IFLAG_NONZERO : IFLAG_POS)) {
      root = (ix == 0) ? logint(n, k) : rootint(n, k);
      if (svret) sv_setuv(SvRV(svret), ix == 0 ? ipow(k,root) : ipow(root,k));
      XSRETURN_UV(root);
    }
    switch (ix) {
      case 0: (void)_vcallsubn(aTHX_ G_SCALAR, (svret == 0) ? (VCALL_GMP|VCALL_PP) : (VCALL_PP), "logint", items, 47); break;
      case 1: (void)_vcallsubn(aTHX_ G_SCALAR, (svret == 0) ? (VCALL_GMP|VCALL_PP) : (VCALL_PP), "rootint", items, 40); break;
      default: break;
    }
    return;

void divrem(IN SV* sva, IN SV* svb)
  ALIAS:
    fdivrem = 1
    cdivrem = 2
    tdivrem = 3
  PREINIT:
    int astatus, bstatus;
    UV D, d;
    IV iD, id;
  PPCODE:
    astatus = _validate_and_set(&D, aTHX_ sva, IFLAG_ANY);
    bstatus = _validate_and_set(&d, aTHX_ svb, IFLAG_ANY);
    if (astatus == 1 && bstatus == 1 && (ix != 2 || (D % d) == 0)) {
      if (d == 0) croak("divrem: divide by zero");
      XPUSHs(sv_2mortal(newSVuv( D / d )));
      XPUSHs(sv_2mortal(newSVuv( D % d )));
      XSRETURN(2);
    } else if (ix == 2 && astatus == 1 && bstatus == 1 && d <= (UV)IV_MAX) {
      if (d == 0) croak("cdivrem: divide by zero");
      if ( (D % d) == 0 ) {
        XPUSHs(sv_2mortal(newSVuv( D / d )));
        XPUSHs(sv_2mortal(newSVuv( D % d )));
      } else {
        XPUSHs(sv_2mortal(newSVuv( D/d + 1 )));
        XPUSHs(sv_2mortal(newSViv( ((IV)D%d) - d )));
      }
      XSRETURN(2);
    } else if (astatus != 0 && bstatus != 0 &&
               _validate_and_set((UV*)&iD, aTHX_ sva, IFLAG_IV) != 0 &&
               _validate_and_set((UV*)&id, aTHX_ svb, IFLAG_IV) != 0) {
      /* Both values fit in an IV */
      IV q, r;
      if (id == 0) croak("divrem: divide by zero");
      switch (ix) {
        case 0:  edivrem(&q, &r, iD, id); break;
        case 1:  fdivrem(&q, &r, iD, id); break;
        case 2:  cdivrem(&q, &r, iD, id); break;
        case 3:
        default: tdivrem(&q, &r, D, d); break;
      }
      XPUSHs(sv_2mortal(newSViv( q )));
      XPUSHs(sv_2mortal(newSViv( r )));
      XSRETURN(2);
    }
    switch (ix) {
      case 0:  _vcallsubn(aTHX_ GIMME_V, VCALL_PP|VCALL_GMP, "divrem", items, 52); break;
      case 1:  _vcallsubn(aTHX_ GIMME_V, VCALL_PP|VCALL_GMP, "fdivrem", items, 53); break;
      case 2:  _vcallsubn(aTHX_ GIMME_V, VCALL_PP|VCALL_GMP, "cdivrem", items, 53); break;
      case 3:
      default: _vcallsubn(aTHX_ GIMME_V, VCALL_PP|VCALL_GMP, "tdivrem", items, 52); break;
    }
    return; /* skip implicit PUTBACK */

void lshiftint(IN SV* svn, IN unsigned long k = 1)
  ALIAS:
    rshiftint = 1
    rashiftint = 2
  PREINIT:
    int status;
    UV n, nk;
  PPCODE:
    status = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
    if (k == 0) {
      ST(0) = svn;
      XSRETURN(1);
    }
    if (status != 0 && ix > 0 && k >= BITS_PER_WORD)  /* Big right shift */
      XSRETURN_UV(0);
    if (status == 1 && k < BITS_PER_WORD) {
      if (ix > 0)                       XSRETURN_UV(n >> k);  /* Right shift */
      if ( ((n << k) >> k) == n)        XSRETURN_UV(n << k);  /* Left shift */
      /* Fall through -- left shift needs more bits */
    } else if (status == -1 && k < BITS_PER_WORD) {
      n = neg_iv(n);
      nk = n >> k;
      if (ix == 1)                      XSRETURN_IV(-nk);
      if (ix == 2)                      XSRETURN_IV(((nk<<k)==n) ? -nk : -nk-1);
      if (((n << (k+1)) >> (k+1)) == n) XSRETURN_IV(-(n << k));
      /* Fall through -- left shift needs more bits */
    }
    switch (ix) {
      case 0:  _vcallsub_with_gmpobj(0.53,"lshiftint");  break;
      case 1:  _vcallsub_with_gmpobj(0.53,"rshiftint");  break;
      case 2:
      default: _vcallsub_with_gmpobj(0.53,"rashiftint"); break;
    }
    objectify_result(aTHX_ ST(0), ST(0));
    return;

void
gcdext(IN SV* sva, IN SV* svb)
  PREINIT:
    IV u, v, d, a, b;
  PPCODE:
    if (_validate_and_set((UV*)&a, aTHX_ sva, IFLAG_IV) &&
        _validate_and_set((UV*)&b, aTHX_ svb, IFLAG_IV)) {
      d = gcdext(a, b, &u, &v, 0, 0);
      XPUSHs(sv_2mortal(newSViv( u )));
      XPUSHs(sv_2mortal(newSViv( v )));
      XPUSHs(sv_2mortal(newSViv( d )));
    } else {
      _vcallsubn(aTHX_ GIMME_V, VCALL_PP, "gcdext", items, 0);
      return; /* skip implicit PUTBACK */
    }

void
stirling(IN UV n, IN UV m, IN UV type = 1)
  PPCODE:
    if (type != 1 && type != 2 && type != 3)
      croak("stirling type must be 1, 2, or 3");
    if (n == m)
      XSRETURN_UV(1);
    else if (n == 0 || m == 0 || m > n)
      XSRETURN_UV(0);
    else if (type == 3) {
      UV s = stirling3(n, m);
      if (s != 0) XSRETURN_UV(s);
    } else if (type == 2) {
      IV s = stirling2(n, m);
      if (s != 0) XSRETURN_IV(s);
    } else if (type == 1) {
      IV s = stirling1(n, m);
      if (s != 0) XSRETURN_IV(s);
    }
    _vcallsub_with_gmpobj(0.26,"stirling");
    objectify_result(aTHX_ ST(0), ST(0));
    return;

NV
_XS_ExponentialIntegral(IN SV* x)
  ALIAS:
    _XS_LogarithmicIntegral = 1
    _XS_RiemannZeta = 2
    _XS_RiemannR = 3
    _XS_LambertW = 4
  PREINIT:
    NV nv, ret;
  CODE:
    nv = SvNV(x);
    switch (ix) {
      case 0: ret = Ei(nv); break;
      case 1: ret = Li(nv); break;
      case 2: ret = (NV) ld_riemann_zeta(nv); break;
      case 3: ret = (NV) RiemannR(nv,0); break;
      case 4:
      default:ret = lambertw(nv); break;
    }
    RETVAL = ret;
  OUTPUT:
    RETVAL


void euler_phi(IN SV* svlo, IN SV* svhi = 0)
  ALIAS:
    moebius = 1
  PREINIT:
    UV lo, hi;
    int lostatus, histatus;
    uint32_t mask;
  PPCODE:
    mask = (ix == 1 && items == 1)  ?  IFLAG_ABS  :  IFLAG_ANY;
    lostatus = _validate_and_set(&lo, aTHX_ svlo, mask);
    histatus = (svhi == 0) || _validate_and_set(&hi, aTHX_ svhi, IFLAG_ANY);
    if (svhi == 0 && lostatus != 0) {
      if (ix == 0) XSRETURN_UV( (lostatus == -1) ? 0 : totient(lo) );
      else         RETURN_NPARITY( moebius(lo) );
    } else if (items == 2 && lostatus == 1 && histatus == 1) {
      if (lo <= hi) {
        UV i, count = hi - lo + 1;
        EXTEND(SP, (IV)count);
        if (ix == 0) {
          UV arrlo = (lo < 100) ?  0 : lo;
          UV *totients = range_totient(arrlo, hi);
          for (i = 0; i < count; i++)
            PUSHs(sv_2mortal(newSVuv(totients[i+lo-arrlo])));
          Safefree(totients);
        } else {
          signed char* mu = range_moebius(lo, hi);
          dMY_CXT;
          for (i = 0; i < count; i++)
            PUSH_NPARITY(mu[i]);
          Safefree(mu);
        }
      }
    } else {
      I32 flags = VCALL_PP;
      if (ix == 0) {
        _vcallsubn(aTHX_ GIMME_V, flags, "euler_phi", items, 22);
      } else {
        if (lostatus == 1 && histatus == 1)  flags |= VCALL_GMP;
        _vcallsubn(aTHX_ GIMME_V, flags, "moebius", items, 22);
      }
      return;
    }

void sqrtint(IN SV* svn)
  ALIAS:
    carmichael_lambda = 1
    exp_mangoldt = 2
  PREINIT:
    UV n, r;
  PPCODE:
    if (_validate_and_set(&n, aTHX_ svn, IFLAG_POS)) {
      r = 0;
      switch (ix) {
        case 0:  r = isqrt(n);  break;
        case 1:  r = carmichael_lambda(n);  break;
        case 2:  r = exp_mangoldt(n);  break;
        default: break;
      }
      XSRETURN_UV(r);
    }
    switch (ix) {
      case 0:  _vcallsub_with_gmp(0.40,"sqrtint"); break;
      case 1:  _vcallsub_with_gmp(0.22,"carmichael_lambda"); break;
      case 2:  _vcallsub_with_gmp(0.19,"exp_mangoldt"); break;
      default: break;
    }
    objectify_result(aTHX_ svn, ST(0));
    return; /* skip implicit PUTBACK */

void prime_omega(IN SV* svn)
  ALIAS:
    prime_bigomega = 1
    hammingweight = 2
    is_square_free = 3
  PREINIT:
    UV n, ret;
  PPCODE:
    if (_validate_and_set(&n, aTHX_ svn, IFLAG_ABS)) {
      ret = 0;
      switch (ix) {
        case 0:  ret = prime_omega(n);    break;
        case 1:  ret = prime_bigomega(n); break;
        case 2:  ret = popcnt(n);         break;
        case 3:  ret = is_square_free(n); break;
        default: break;
      }
      RETURN_NPARITY(ret);
    }
    switch (ix) {
      case 0:  _vcallsub_with_gmp(0.53,"prime_omega"); break;
      case 1:  _vcallsub_with_gmp(0.53,"prime_bigomega"); break;
      case 2:  if (_XS_get_callgmp() >= 47) { /* Very fast */
                 _vcallsub_with_gmp(0.47,"hammingweight");
               } else {                       /* Better than PP */
                 char* ptr;  STRLEN len;  ptr = SvPV(svn, len);
                 XSRETURN_UV(mpu_popcount_string(ptr, len));
               }
               break;
      case 3: _vcallsub_with_gmp(0.53,"is_square_free"); break;
      default: break;
    }
    /* No objectify */
    return;

void factorial(IN SV* svn)
  ALIAS:
    subfactorial = 1
    fubini = 2
    primorial = 3
    pn_primorial = 4
    sumtotient = 5
  PREINIT:
    UV n, r;
  PPCODE:
    if (_validate_and_set(&n, aTHX_ svn, IFLAG_POS)) {
      r = 0;
      switch(ix) {
        case 0:  r = factorial(n);    break;
        case 1:  r = subfactorial(n); break;
        case 2:  r = fubini(n);       break;
        case 3:  r = primorial(n);    break;
        case 4:  r = pn_primorial(n); break;
        case 5:  r = sumtotient(n);   break;
        default: break;
      }
      if (n == 0 || r > 0) XSRETURN_UV(r);
      if (ix == 5) {  /* Probably an overflow, try 128-bit. */
        UV hicount, count;
        int retok = sumtotient128(n, &hicount, &count);
        if (retok == 1 && hicount > 0)
          RETURN_128(hicount, count);
        if (retok == 1)
          XSRETURN_UV(count);
      }
    }
    switch (ix) {
      case 0:  _vcallsub_with_pp("factorial"); break;  /* use PP */
      case 1:  _vcallsub_with_gmp(0.51,"subfactorial"); break;
      case 2:  _vcallsub_with_gmp(0.00,"fubini"); break;
      case 3:  _vcallsub_with_gmp(0.37,"primorial"); break;
      case 4:  _vcallsub_with_gmp(0.37,"pn_primorial"); break;
      case 5:  _vcallsub_with_gmp(0.00,"sumtotient"); break;
      default: break;
    }
    objectify_result(aTHX_ svn, ST(0));
    return;

void binomial(IN SV* svn, IN SV* svk)
  PREINIT:
    int nstatus, kstatus;
    UV n, k, ret;
  PPCODE:
    nstatus = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY);
    kstatus = _validate_and_set(&k, aTHX_ svk, IFLAG_ANY);
    if (nstatus != 0 && kstatus != 0) {
      if ( (nstatus == 1 && (kstatus == -1 || k > n)) ||
           (nstatus ==-1 && (kstatus == -1 && k > n)) )
         XSRETURN_UV(0);
      if (kstatus == -1)
        k = n - k; /* n<0,k<=n:  (-1)^(n-k) * binomial(-k-1,n-k) */
      if (nstatus == -1) {
        ret = binomial( neg_iv(n)+k-1, k );
        if (ret > 0 && ret <= (UV)IV_MAX)
          XSRETURN_IV( (IV)ret * ((k&1) ? -1 : 1) );
      } else if (nstatus == 1) {
        ret = binomial(n, k);
        if (ret != 0) XSRETURN_UV(ret);
      }
    }
    _vcallsub_with_gmp(0.22,"binomial");
    objectify_result(aTHX_ svn, ST(0));
    return;

void falling_factorial(IN SV* svn, IN SV* svk)
  ALIAS:
    rising_factorial = 1
  PREINIT:
    int nstatus, kstatus;
    UV n, k;
  PPCODE:
    nstatus = _validate_and_set(&n, aTHX_ svn, IFLAG_ANY | IFLAG_IV);
    kstatus = _validate_and_set(&k, aTHX_ svk, IFLAG_POS);
    if (nstatus == 1 && kstatus == 1) {
      UV ret = (ix==0) ? falling_factorial(n,k) : rising_factorial(n,k);
      if (ret != UV_MAX) XSRETURN_UV(ret);
    } else {
      IV in = (IV)n;
      IV ret = (ix==0) ? falling_factorial_s(in,k) : rising_factorial_s(in,k);
      if (ret != IV_MAX) XSRETURN_IV(ret);
    }
    if (ix == 0) _vcallsub_with_gmp(0.51,"falling_factorial");
    else         _vcallsub_with_gmp(0.51,"rising_factorial");
    objectify_result(aTHX_ svn, ST(0));
    return;

void mertens(IN SV* svn)
  ALIAS:
    liouville = 1
    sumliouville = 2
    is_pillai = 3
    is_congruent_number = 4
    hclassno = 5
    ramanujan_tau = 6
  PREINIT:
    UV n;
    int status;
  PPCODE:
    status = _validate_and_set(&n, aTHX_ svn, (ix < 5) ? IFLAG_POS : IFLAG_ANY);
    if (status == -1)
      XSRETURN_IV(0);
    if (status == 1) {
      IV r = 0;
      switch(ix) {
        case 0:  r = mertens(n); break;
        case 1:  r = liouville(n); break;
        case 2:  r = sumliouville(n); break;
        case 3:  r = pillai_v(n); break;
        case 4:  r = is_congruent_number(n); break;
        case 5:  r = hclassno(n); break;
        case 6:  r = ramanujan_tau(n);
                 if (r == 0 && n != 0)
                   status = 0;
                 break;
        default: break;
      }
      if (status != 0) RETURN_NPARITY(r);
    }
    switch (ix) {
      case 0:  _vcallsub_with_pp("mertens"); break;
      case 1:  _vcallsub_with_gmp(0.22,"liouville"); break;
      case 2:  _vcallsub_with_gmp(0.00,"sumliouville"); break;
      case 3:  _vcallsub_with_gmp(0.00,"is_pillai"); break;
      case 4:  _vcallsub_with_gmp(0.00,"is_congruent_number"); break;
      case 5:  _vcallsub_with_pp("hclassno"); break;
      case 6:  _vcallsub_with_gmp(0.32,"ramanujan_tau"); break;
      default: break;
    }
    objectify_result(aTHX_ svn, ST(0));
    return;

void chebyshev_theta(IN SV* svn)
  ALIAS:
    chebyshev_psi = 1
  PREINIT:
    UV n;
  PPCODE:
    if (_validate_and_set(&n, aTHX_ svn, IFLAG_POS)) {
      NV r = (ix==0)  ?  chebyshev_theta(n)  :  chebyshev_psi(n);
      XSRETURN_NV(r);
    }
    _vcallsub_with_pp( (ix==0) ? "chebyshev_theta" : "chebyshev_psi" );
    /* Result is FP, don't objectify */
    return;


#define RETURN_SET_VALS(s)   /* Return sorted set values */ \
  { \
    UV *sdata; \
    unsigned long slen = iset_size(s); \
    int sign = iset_sign(s); \
    if (GIMME_V == G_SCALAR) { \
      iset_destroy(&s); \
      XSRETURN_UV(slen); \
    } \
    New(0, sdata, slen, UV); \
    iset_allvals(s, sdata); \
    iset_destroy(&s); \
    RETURN_LIST_VALS( slen, sdata, sign ); \
  }

void sumset(IN SV* sva, IN SV* svb = 0)
  PROTOTYPE: $;$
  PREINIT:
    int atype, btype, stype, sign;
    UV *ra, *rb;
    unsigned long alen, blen,  i, j;
    iset_t s;
  PPCODE:
    atype = arrayref_to_int_array(aTHX_ &alen, &ra, 1, sva, "sumset arg 1");
    if (svb == 0 || atype == IARR_TYPE_BAD) {
      rb = ra;
      blen = alen;
      btype = atype;
    } else {
      btype = arrayref_to_int_array(aTHX_ &blen, &rb, 1, svb, "sumset arg 2");
    }
    if (atype == IARR_TYPE_BAD || btype == IARR_TYPE_BAD)
      stype = IARR_TYPE_BAD;
    else
      stype = type_of_sumset(atype, btype, ra[0],ra[alen-1], rb[0],rb[blen-1]);
    if (stype == IARR_TYPE_BAD) {
      if (rb != ra) Safefree(rb);
      Safefree(ra);
      _vcallsubn(aTHX_ GIMME_V, VCALL_PP, "sumset", items, 0);
      return;
    }
    if (alen == 0 || blen == 0) {
      if (rb != ra) Safefree(rb);
      Safefree(ra);
      XSRETURN_EMPTY;
    }
    sign = IARR_TYPE_TO_STATUS(stype);
    /* Sumset */
    s = iset_create( 10UL * (alen+blen) );
    for (i = 0; i < alen; i++)
      for (j = 0; j < blen; j++)
        iset_add(&s, ra[i]+rb[j], sign);
    if (rb != ra) Safefree(rb);
    Safefree(ra);
    RETURN_SET_VALS(s);

void setbinop(SV* block, IN SV* sva, IN SV* svb = 0)
  PROTOTYPE: &$;$
  CODE:
  {
    AV *ava, *avb;
    GV *agv, *bgv, *gv;
    SV **asv, **bsv;
    iset_t s;
    int alen, blen, i, j, status;
    UV ret;
    HV *stash;
    CV *cv = sv_2cv(block, &stash, &gv, 0);

    if (cv == Nullcv) croak("Not a subroutine reference");

    /* ====== Get and store the input array references ====== */
    if ((!SvROK(sva)) || (SvTYPE(SvRV(sva)) != SVt_PVAV))
      croak("setbinop: first array arg must be an array reference");
    ava = (AV*) SvRV(sva);
    alen = av_len(ava);
    if (alen < 0) XSRETURN_EMPTY;

    if (svb == 0) {
      avb = ava;
    } else {
      if ((!SvROK(svb)) || (SvTYPE(SvRV(svb)) != SVt_PVAV))
        croak("setbinop: second array arg must be an array reference");
      avb = (AV*) SvRV(svb);
    }
    blen = av_len(avb);
    if (blen < 0) XSRETURN_EMPTY;

    /* ====== Walk the array references storing the SV pointers ====== */
    New(0, asv, alen+1, SV*);
    for (i = 0; i <= alen; i++) {
      SV** iv = av_fetch(ava, i, 0);
      if (iv == 0) break;
      asv[i] = *iv;
    }
    if (i > alen) {
      if (avb == ava) {
        bsv = asv;
        j = blen+1;
      } else {
        New(0, bsv, blen+1, SV*);
        for (j = 0; j <= blen; j++) {
          SV** iv = av_fetch(avb, j, 0);
          if (iv == 0) break;
          bsv[j] = *iv;
        }
      }
    }
    if (i <= alen || j <= blen) {
      _vcallsubn(aTHX_ GIMME_V, VCALL_PP, "setbinop", items, 0);
      return;
    }

    /* ====== Call block on cross product, insert result into set ====== */
    agv = gv_fetchpv("a", GV_ADD, SVt_PV);
    bgv = gv_fetchpv("b", GV_ADD, SVt_PV);
    SAVESPTR(GvSV(agv));
    SAVESPTR(GvSV(bgv));
    s = iset_create( 4UL *  ((unsigned long)alen + (unsigned long)blen + 2) );
#ifdef dMULTICALL
    if (!CvISXSUB(cv)) {
      dMULTICALL;
      I32 gimme = G_SCALAR;
      PUSH_MULTICALL(cv);
      for (i = 0; i <= alen; i++) {
        GvSV(agv) = asv[i];
        for (j = 0; j <= blen; j++) {
          GvSV(bgv) = bsv[j];
          { ENTER; MULTICALL; LEAVE; }
          status = _validate_and_set(&ret, aTHX_ *PL_stack_sp, IFLAG_ANY);
          if (status != 0)  iset_add(&s, ret, status);
          if (status == 0 || iset_is_invalid(s)) break;
        }
        if (j <= blen) break;
      }
      FIX_MULTICALL_REFCOUNT;
      POP_MULTICALL;
    }
    else
#endif
    {
      for (i = 0; i <= alen; i++) {
        GvSV(agv) = asv[i];
        for (j = 0; j <= blen; j++) {
          dSP;
          GvSV(bgv) = bsv[j];
          PUSHMARK(SP);
          call_sv((SV*)cv, G_SCALAR);
          status = _validate_and_set(&ret, aTHX_ *PL_stack_sp, IFLAG_ANY);
          if (status != 0)  iset_add(&s, ret, status);
          if (status == 0 || iset_is_invalid(s)) break;
        }
        if (j <= blen) break;
      }
    }

    /* ====== Free cached SV pointers.  Call PP if we not complete. ====== */
    if (bsv != asv) Safefree(bsv);
    Safefree(asv);
    /* ====== Call PP if not finished. Return scalar if all needed. ====== */
    if (i <= alen || j <= blen) {
      iset_destroy(&s);
      _vcallsubn(aTHX_ GIMME_V, VCALL_PP, "setbinop", items, 0);
      return;
    }
    /* ====== Get sorted set values.  Put on return stack. ====== */
    RETURN_SET_VALS(s);
  }

void setunion(IN SV* sva, IN SV* svb)
  PROTOTYPE: $$
  ALIAS:
    setintersect = 1
    setminus = 2
    setdelta = 3
  PREINIT:
    int atype, btype;
    UV *ra, *rb;
    unsigned long alen, blen;
  PPCODE:
    /* Get the integers and check if they are sorted unique integers first. */
    atype = arrayref_to_int_array(aTHX_ &alen, &ra, 1, sva, "setunion arg 1");
    btype = arrayref_to_int_array(aTHX_ &blen, &rb, 1, svb, "setunion arg 2");

    if (CAN_COMBINE_IARR_TYPES(atype,btype)) {
      UV *r;
      unsigned long rlen = 0, ia = 0, ib = 0;
      int pcmp = (atype == IARR_TYPE_NEG || btype == IARR_TYPE_NEG) ? 0 : 1;

      if (ix == 0) {        /* union */
        New(0, r, alen + blen, UV);
        while (ia < alen && ib < blen) {
          if (ra[ia] == rb[ib]) {
            r[rlen++] = ra[ia];
            ia++; ib++;
          } else {
            if (SIGNED_CMP_LT(pcmp, ra[ia], rb[ib])) r[rlen++] = ra[ia++];
            else                                     r[rlen++] = rb[ib++];
          }
        }
        while (ia < alen)   r[rlen++] = ra[ia++];
        while (ib < blen)   r[rlen++] = rb[ib++];
      } else if (ix == 1) { /* intersect */
        New(0, r, (alen > blen) ? alen : blen, UV);
        while (ia < alen && ib < blen) {
          if (ra[ia] == rb[ib]) {
            r[rlen++] = ra[ia];
            ia++; ib++;
          } else {
            if (SIGNED_CMP_LT(pcmp, ra[ia], rb[ib])) ia++;
            else                                     ib++;
          }
        }
      } else if (ix == 2) { /* minus (difference) */
        New(0, r, alen, UV);
        while (ia < alen && ib < blen) {
          if (ra[ia] == rb[ib]) {
            ia++; ib++;
          } else {
            if (SIGNED_CMP_LT(pcmp, ra[ia], rb[ib])) r[rlen++] = ra[ia++];
            else                                     ib++;
          }
        }
        while (ia < alen)   r[rlen++] = ra[ia++];
      } else if (ix == 3) { /* delta (symmetric difference) */
        New(0, r, alen + blen, UV);
        while (ia < alen && ib < blen) {
          if (ra[ia] == rb[ib]) {
            ia++; ib++;
          } else {
            if (SIGNED_CMP_LT(pcmp, ra[ia], rb[ib])) r[rlen++] = ra[ia++];
            else                                     r[rlen++] = rb[ib++];
          }
        }
        while (ia < alen)   r[rlen++] = ra[ia++];
        while (ib < blen)   r[rlen++] = rb[ib++];
      }
      Safefree(ra);
      Safefree(rb);
      if (GIMME_V == G_SCALAR) {
        Safefree(r);
        XSRETURN_UV(rlen);
      }
      RETURN_LIST_VALS(rlen, r, pcmp);
    }
    /* if (atype != IARR_TYPE_BAD && btype != IARR_TYPE_BAD) { .. isets .. } */
    Safefree(ra);
    Safefree(rb);
    switch (ix) {
      case 0: CALLPPSUB("setunion");     break;
      case 1: CALLPPSUB("setintersect"); break;
      case 2: CALLPPSUB("setminus");     break;
      case 3:
      default:CALLPPSUB("setdelta");     break;
    }
    return;

void set_is_disjoint(IN SV* sva, IN SV* svb)
  PROTOTYPE: $$
  ALIAS:
    set_is_equal = 1
    set_is_subset = 2
    set_is_proper_subset = 3
    set_is_superset = 4
    set_is_proper_superset = 5
    set_is_proper_intersection = 6
  PREINIT:
    int atype, btype, ret;
    UV *ra, *rb;
    unsigned long alen, blen, inalen, inblen;
  PPCODE:
    /* If one set is much smaller than the other, it would be faster using
     * is_in_set().  We'll keep things simple and slurp in both sets. */

    /* THIS ASSUMES THE INPUT LISTS HAVE NO DUPLICATES */
    inalen = inblen = 0;
    if (SvROK(sva) && SvTYPE(SvRV(sva)) == SVt_PVAV && SvROK(svb) && SvTYPE(SvRV(svb)) == SVt_PVAV) {
      /* Shortcut on length if we can to skip intersection. */
      inalen = av_count((AV*) SvRV(sva));
      inblen = av_count((AV*) SvRV(svb));
      if ( (ix == 1 && inalen != inblen) ||
           (ix == 2 && inalen <  inblen) || (ix == 3 && inalen <= inblen) ||
           (ix == 4 && inalen >  inblen) || (ix == 5 && inalen >= inblen) )
        RETURN_NPARITY(0);
    }

    /* Get the integers as sorted arrays of IV or UV */
    atype = arrayref_to_int_array(aTHX_ &alen, &ra, 1, sva, "setunion arg 1");
    btype = arrayref_to_int_array(aTHX_ &blen, &rb, 1, svb, "setunion arg 2");

    if (CAN_COMBINE_IARR_TYPES(atype,btype)) {
      unsigned long rlen = 0, ia = 0, ib = 0;
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
      ret = 0;
      switch (ix) {
        case 0: if (rlen == 0) ret = 1;   break;
        case 1: if (alen == blen && rlen == blen) ret = 1;  break;
        case 2: if (alen >= blen && rlen == blen) ret = 1;  break;
        case 3: if (alen >  blen && rlen == blen) ret = 1;  break;
        case 4: if (alen <= blen && rlen == alen) ret = 1;  break;
        case 5: if (alen <  blen && rlen == alen) ret = 1;  break;
        case 6:
        default:if (rlen > 1 && rlen < alen && rlen < blen) ret = 1; break;
      }
      RETURN_NPARITY(ret);
    }
    Safefree(ra);
    Safefree(rb);
    switch (ix) {
      case 0: CALLPPSUB("set_is_disjoint");       break;
      case 1: CALLPPSUB("set_is_equal");          break;
      case 2: CALLPPSUB("set_is_subset");         break;
      case 3: CALLPPSUB("set_is_proper_subset");  break;
      case 4: CALLPPSUB("set_is_superset");       break;
      case 5: CALLPPSUB("set_is_proper_superset");break;
      case 6:
      default:CALLPPSUB("set_is_proper_intersection"); break;
    }
    return;

void setcontains(IN SV* sva, IN SV* svb)
  PROTOTYPE: $$
  PREINIT:
    UV b;
    AV *ava, *avb;
    int bstatus, subset, alast, blast, i;
  PPCODE:
    /* First arg must be an array reference (set) */
    if ((!SvROK(sva)) || (SvTYPE(SvRV(sva)) != SVt_PVAV))
      croak("setcontains first argument must be an array reference");
    ava = (AV*) SvRV(sva);
    alast = av_len(ava);
    /* Case of the second argument being a single integer. */
    if (!SvROK(svb) || SvTYPE(SvRV(svb)) != SVt_PVAV) {
      bstatus = _validate_and_set(&b, aTHX_ svb, IFLAG_ANY);
      subset = is_in_set(aTHX_ ava, 0, bstatus, b);
    } else { /* The second argument is an array reference (set) */
      avb = (AV*) SvRV(svb);
      blast = av_len(avb);
      if (blast > alast) {
        subset = 0;  /* cannot fit */
      } else {
        set_data_t cache = init_set_lookup_cache(aTHX_ ava);
        for (i = 0, subset = 1; i <= blast && subset == 1; i++) {
          SV **iv = av_fetch(avb, i, 0);
          if (iv != 0) {
            bstatus = _validate_and_set(&b, aTHX_ *iv, IFLAG_ANY);
            subset = is_in_set(aTHX_ ava, &cache, bstatus, b);
          }
        }
        free_set_lookup_cache(&cache);
      }
      /* Much slower in most cases:
       * subset = is_subset_read_all(aTHX_ sva, svb);
       */
    }
    if (subset != -1)
      RETURN_NPARITY(subset);
    CALLPPSUB("setcontains");
    return;

void setinsert(IN SV* sva, IN SV* svb)
  PROTOTYPE: $$
  PREINIT:
    UV b;
    AV *ava;
    int bstatus;
  PPCODE:
    /* First arg must be an array reference (set) */
    if ((!SvROK(sva)) || (SvTYPE(SvRV(sva)) != SVt_PVAV))
      croak("setcontains first argument must be an array reference");
    ava = (AV*) SvRV(sva);
    /* Case of the second argument being a single integer. */
    if (!SvROK(svb) || SvTYPE(SvRV(svb)) != SVt_PVAV) {
      bstatus = _validate_and_set(&b, aTHX_ svb, IFLAG_ANY);
      if (bstatus != 0) {
        Size_t count = av_count(ava);
        int index = index_in_set(aTHX_ ava, 0, bstatus, b);
        if (index == 0) RETURN_NPARITY(0);
        if (index > 0) { /* Insert into array reference */
          SV *newb, **arr;
          SV* newsvb = NEWSVINT(bstatus, b);
          if ((Size_t)index > count/2) {
            av_push(ava, newsvb);
            if ((Size_t)index <= count) {
              arr = AvARRAY(ava);
              newb = arr[count];
              memmove(arr+index, arr+index-1, sizeof(SV*) * (count-(index-1)));
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
          RETURN_NPARITY(1);
        }
      }
    } else {
      UV *rb;
      unsigned long i, blen, nbeg, nmid, nend, nmidcheck;
      int btype = arrayref_to_int_array(aTHX_ &blen, &rb, 1, svb, "setinsert arg 2");
      bstatus = IARR_TYPE_TO_STATUS(btype);
      if (bstatus != 0) {
        Size_t alen = av_count(ava);
        int alostatus, ahistatus;
        UV  alo, ahi;
        set_data_t cache;

        /* 1. ava is empty.  push everything and we're done. */
        if (alen == 0) {
          for (i = 0; i < blen; i++)
            av_push(ava, NEWSVINT(bstatus, rb[i]));
          Safefree(rb);
          RETURN_NPARITY(blen);
        }
        /* Get ready to insert */
        cache = init_set_lookup_cache(aTHX_ ava);
        /* Get hi and lo values of set */
        if (_sc_set_lohi(aTHX_ ava, &cache, 0, alen-1, &alostatus, &ahistatus, &alo, &ahi) >= 0) {
          nbeg = nend = nmid = 0;
          /* 1. Find out how many elements go in front. */
          while (nbeg < blen && _sign_cmp(bstatus,rb[nbeg],alostatus,alo) < 0)
            nbeg++;
            ;
          /* 2. Find out how many elements go at the end. */
          while (nend < blen-nbeg && _sign_cmp(bstatus,rb[blen-1-nend],ahistatus,ahi) > 0)
            nend++;
          /* 3. In-place insert everything in the middle. */
          nmidcheck = blen - nbeg - nend;
          if (nmidcheck > 0) {
            unsigned long *insert_idx;
            SV **insert_sv;
            New(0, insert_idx, nmidcheck, unsigned long);
            New(0, insert_sv,  nmidcheck, SV*);
            for (i = nbeg; bstatus != 0 && i < blen-nend; i++) {
              int index = index_in_set(aTHX_ ava, &cache, bstatus, rb[i]);
              if (index < 0) croak("bigint value found in interior of non-bigint sorted set");
              if (index > 0) {
                insert_sv[nmid]  = NEWSVINT(bstatus,rb[i]);/* Value to insert */
                insert_idx[nmid] = index-1;                /* Where to insert */
                nmid++;
              }
            }
            if (nmid > 0) {
              SV** arr;
              Size_t newalen = alen + nmid;
              unsigned long index_lastorig = alen-1;
              unsigned long index_moveto   = index_lastorig + nmid;

              av_extend(ava, newalen);
              /* Push new values on end so Perl calculates array right. */
              for (i = 0; i < nmid; i++)
                av_push(ava, insert_sv[i]);
              arr = AvARRAY(ava);
              /* SV* pointer manipulation to insert new values in place. */
              for (i = 0; i < nmid; i++) {
                unsigned long j   = nmid-1-i;
                unsigned long idx = insert_idx[j];
                if (i == 0 || idx != insert_idx[j+1]) {
                  unsigned long nmove  = index_lastorig - idx + 1;
                  unsigned long moveto = index_moveto - nmove + 1;
                  memmove(arr+moveto, arr+idx, sizeof(SV*) * nmove);
                  index_lastorig -= nmove;
                  index_moveto -= nmove;
                }
                arr[index_moveto--] = insert_sv[j];
              }
            }
            Safefree(insert_sv);
            Safefree(insert_idx);
          }
          /* 4. Insert at front */
          if (nbeg > 0) {
            av_unshift(ava, nbeg);
            for (i = 0; i < nbeg; i++)
              av_store(ava, i, NEWSVINT(bstatus, rb[i]));
          }
          /* 5. Push onto back */
          if (nend > 0) {
            for (i = 0; i < nend; i++)
              av_push(ava, NEWSVINT(bstatus, rb[blen-nend+i]));
          }
          Safefree(rb);
          RETURN_NPARITY(nbeg+nmid+nend);
        }
        free_set_lookup_cache(&cache);
      }
      Safefree(rb);
    }
    CALLPPSUB("setinsert");
    return;


void is_sidon_set(IN SV* sva)
  PROTOTYPE: $
  PREINIT:
    int itype, is_sidon;
    unsigned long len, i, j;
    UV *data;
    iset_t s;
  PPCODE:
    itype = arrayref_to_int_array(aTHX_ &len, &data, 1, sva,"is_sidon_set arg");
    if (itype == IARR_TYPE_NEG) {  /* All elements must be non-negative. */
      Safefree(data);
      RETURN_NPARITY(0);
    }
    /* If any bigints or we cannot add the values in 64-bits, call PP. */
    if (itype == IARR_TYPE_BAD || itype == IARR_TYPE_POS) {
      Safefree(data);
      CALLPPSUB("is_sidon_set");
      return;
    }
    /* Check if the set is a Sidon set. */
    is_sidon = 1;
    s = iset_create( 20UL * len );
    for (i = 0; i < len && is_sidon; i++)
      for (j = i; j < len; j++)
        if (!iset_add(&s, data[i] + data[j], 1))
          { is_sidon = 0; break; }
    Safefree(data);
    iset_destroy(&s);
    RETURN_NPARITY(is_sidon);

void is_sumfree_set(IN SV* sva)
  PROTOTYPE: $
  PREINIT:
    int itype, is_sumfree;
    unsigned long len, i, j;
    UV *data;
    iset_t sa;
  PPCODE:
    itype = arrayref_to_int_array(aTHX_ &len, &data,1,sva,"is_sumfree_set arg");
    /* Check for UV overflow on sum */
    if (itype == IARR_TYPE_POS)  itype = IARR_TYPE_BAD;
    /* Check for IV overflow on sum */
    if (itype == IARR_TYPE_NEG) {
      for (i = 0; i < len; i++)
        if ((IV)data[i] > (IV_MAX/2) || (IV)data[i] < (IV_MIN/2))
          break;
      if (i < len)  itype = IARR_TYPE_BAD;
    }
    /* Call PP code if bigints or summation overflow */
    if (itype == IARR_TYPE_BAD) {
      Safefree(data);
      CALLPPSUB("is_sumfree_set");
      return;
    }
    sa = iset_create_from_array(data, len, IARR_TYPE_TO_STATUS(itype));
    is_sumfree = 1;
    for (i = 0; i < len && is_sumfree; i++)
      for (j = i; j < len; j++)
        if (iset_contains(sa, data[i] + data[j]))
          { is_sumfree = 0; break; }
    iset_destroy(&sa);
    Safefree(data);
    RETURN_NPARITY(is_sumfree);

void toset(IN SV* sva)
  PREINIT:
    int atype;
    UV *ra;
    unsigned long alen;
  PPCODE:
    atype = arrayref_to_int_array(aTHX_ &alen, &ra, 1, sva, "toset");
    if (atype == IARR_TYPE_BAD) {
      CALLPPSUB("toset");
      return;
    }
    RETURN_LIST_VALS(alen, ra, atype != IARR_TYPE_NEG);


void vecsortr(IN SV* sva)
  PROTOTYPE: $
  ALIAS:
    vecsortrr = 1
  PREINIT:
    int type;
    unsigned long len;
    UV *L;
  PPCODE:
    type = arrayref_to_int_array(aTHX_ &len, &L, 0, sva, "vecsort arg");
    if (type == IARR_TYPE_ANY || type == IARR_TYPE_POS) {
      sort_uv_array(L, len);
    } else if (type == IARR_TYPE_NEG) {
      sort_iv_array((IV*)L, len);
    } else {
      Safefree(L);
      CALLPPSUB((ix == 0) ? "vecsortr" : "vecsortrr");
      return;
    }
    if (ix == 0) { RETURN_LIST_VALS( len, L, (type != IARR_TYPE_NEG) ); }
    else         { RETURN_LIST_REF(  len, L, (type != IARR_TYPE_NEG) ); }

void vecsort(...)
  PROTOTYPE: @
  PREINIT:
    UV *L;
    int type;
    I32 i;
  PPCODE:
    if (items == 0) XSRETURN_EMPTY;
    New(0, L, items, UV);
    type = IARR_TYPE_ANY;
    for (i = 0; i < items; i++) {
      UV n;
      int istatus = _validate_and_set(&n, aTHX_ ST(i), IFLAG_ANY);
      if (istatus == -1) {
        type |= IARR_TYPE_NEG;
      } else if (istatus == 1 && n > (UV)IV_MAX) {
        type |= IARR_TYPE_POS;
      }
      if (istatus == 0 || type == IARR_TYPE_BAD) break;
      L[i] = n;
    }
    if (i < items) {
      Safefree(L);
      CALLPPSUB("vecsort");
      return;
    }
    if (type == IARR_TYPE_ANY || type == IARR_TYPE_POS)
      sort_uv_array(L, items);
    else
      sort_iv_array((IV*)L, items);
    RETURN_LIST_VALS( items, L, (type != IARR_TYPE_NEG) );


void
numtoperm(IN UV n, IN SV* svk)
  PREINIT:
    UV k;
    int i, S[32];
  PPCODE:
    if (n == 0)
      XSRETURN_EMPTY;
    if (n < 32 && _validate_and_set(&k, aTHX_ svk, IFLAG_ABS) == 1) {
      if (num_to_perm(k, n, S)) {
        dMY_CXT;
        EXTEND(SP, (IV)n);
        for (i = 0; i < (int)n; i++)
          PUSH_NPARITY( S[i] );
        XSRETURN(n);
      }
    }
    _vcallsubn(aTHX_ GIMME_V, VCALL_PP|VCALL_GMP, "numtoperm", items, 47);
    return;

void
permtonum(IN SV* svp)
  PREINIT:
    AV *av;
    UV val, num;
    int plen, i;
  PPCODE:
    if ((!SvROK(svp)) || (SvTYPE(SvRV(svp)) != SVt_PVAV))
      croak("permtonum argument must be an array reference");
    av = (AV*) SvRV(svp);
    plen = av_len(av);
    if (plen < 32) {
      int V[32], A[32] = {0};
      for (i = 0; i <= plen; i++) {
        SV **iv = av_fetch(av, i, 0);
        if (_validate_and_set(&val, aTHX_ iv ? *iv : 0, IFLAG_POS) != 1)
          break;
        if (val > (UV)plen || A[val] != 0) break;
        A[val] = i+1;
        V[i] = val;
      }
      if (i > plen && perm_to_num(plen+1, V, &num))
        XSRETURN_UV(num);
    }
    _vcallsub_with_gmpobj(0.47,"permtonum");
    objectify_result(aTHX_ ST(0), ST(0));
    XSRETURN(1);

void
randperm(IN UV n, IN UV k = 0)
  PREINIT:
    UV i, *S;
    dMY_CXT;
  PPCODE:
    if (items == 1) k = n;
    if (k > n) k = n;
    if (k == 0) XSRETURN_EMPTY;
    New(0, S, k, UV);
    randperm(MY_CXT.randcxt, n, k, S);
    EXTEND(SP, (IV)k);
    for (i = 0; i < k; i++) {
      if (n < 2*CINTS)  PUSH_NPARITY(S[i]);
      else              PUSHs(sv_2mortal(newSVuv(S[i])));
    }
    Safefree(S);

void shuffle(...)
  PROTOTYPE: @
  PREINIT:
    int i, j;
    void* randcxt;
    dMY_CXT;
  PPCODE:
    if (items == 0)
      XSRETURN_EMPTY;
    for (i = 0, randcxt = MY_CXT.randcxt; i < items-1; i++) {
      j = urandomm32(randcxt, items-i);
      { SV* t = ST(i); ST(i) = ST(i+j); ST(i+j) = t; }
    }
    XSRETURN(items);

void
sumdigits(SV* svn, UV ibase = 255)
  PREINIT:
    UV base, sum;
    STRLEN i, len;
    const char* s;
  PPCODE:
    base = (ibase == 255) ? 10 : ibase;
    if (base < 2 || base > 36) croak("sumdigits: invalid base %"UVuf, base);
    sum = 0;
    /* faster for integer input in base 10 */
    if (base == 10 && SVNUMTEST(svn) && (SvIsUV(svn) || SvIVX(svn) >= 0)) {
      UV n, t = my_svuv(svn);
      while ((n=t)) {
        t = n / base;
        sum += n - base*t;
      }
      XSRETURN_UV(sum);
    }
    s = SvPV(svn, len);
    /* If no base given and input is 0x... or 0b..., select base. */
    if (ibase == 255 && len > 2 && s[0] == '0' && (s[1] == 'x' || s[1] == 'b')){
      base = (s[1] == 'x') ? 16 : 2;
      s += 2;
      len -= 2;
    }
    for (i = 0; i < len; i++) {
      UV d = 0;
      const char c = s[i];
      if      (c >= '0' && c <= '9') { d = c - '0';      }
      else if (c >= 'a' && c <= 'z') { d = c - 'a' + 10; }
      else if (c >= 'A' && c <= 'Z') { d = c - 'A' + 10; }
      if (d < base)
        sum += d;
    }
    XSRETURN_UV(sum);

void todigits(SV* svn, int base=10, int length=-1)
  ALIAS:
    todigitstring = 1
    fromdigits = 2
  PREINIT:
    int i, status;
    UV n;
    char *str;
  PPCODE:
    if (base < 2) croak("invalid base: %d", base);
    status = 0;
    if (ix == 0 || ix == 1) {
      status = _validate_and_set(&n, aTHX_ svn, IFLAG_ABS);
    }
    /* todigits with native input */
    if (ix == 0 && status != 0 && length < 128) {
      int digits[128];
      IV len = to_digit_array(digits, n, base, length);
      if (len >= 0) {
        dMY_CXT;
        EXTEND(SP, len);
        for (i = 0; i < len; i++)
          PUSH_NPARITY( digits[len-i-1] );
        XSRETURN(len);
      }
    }
    /* todigitstring with native input */
    if (ix == 1 && status != 0 && length < 128) {
      char s[128+1];
      IV len = to_digit_string(s, n, base, length);
      if (len >= 0) {
        XPUSHs(sv_2mortal(newSVpv(s, len)));
        XSRETURN(1);
      }
    }
    /* todigits or todigitstring base 10 (large size) */
    if ((ix == 0 || ix == 1) && base == 10 && length < 0) {
      STRLEN len;
      str = SvPV(svn, len);
      if (ix == 1) {
        XPUSHs(sv_2mortal(newSVpv(str, len)));
        XSRETURN(1);
      }
      if (len == 1 && str[0] == '0') XSRETURN(0);
      {
        dMY_CXT;
        EXTEND(SP, (IV)len);
        for (i = 0; i < (int)len; i++)
          PUSH_NPARITY(str[i]-'0');
      }
      XSRETURN(len);
    }
    if (ix == 2) { /* fromdigits */
      if (!SvROK(svn)) {  /* string */
        if (from_digit_string(&n, SvPV_nolen(svn), base)) {
          XSRETURN_UV(n);
        }
      } else if (!_is_sv_bigint(aTHX_ svn)) {     /* array ref of digits */
        UV* r = 0;
        int len = arrayref_to_digit_array(aTHX_ &r, (AV*) SvRV(svn), base);
        if (from_digit_to_UV(&n, r, len, base)) {
          Safefree(r);
          XSRETURN_UV(n);
        } else if (from_digit_to_str(&str, r, len, base)){
          Safefree(r);
          XPUSHs( sv_to_bigint(aTHX_ sv_2mortal(newSVpv(str,0))) );
          Safefree(str);
          XSRETURN(1);
        }
        Safefree(r);
      }
    }
    switch (ix) {
      case 0:  _vcallsubn(aTHX_ GIMME_V, VCALL_GMP|VCALL_PP, "todigits", items, 41); break;
      case 1:  _vcallsub_with_gmp(0.00,"todigitstring"); break;
      case 2:
      default: _vcallsub_with_gmp(0.53,"fromdigits");
               objectify_result(aTHX_ 0, ST(0));
               break;
    }
    return;

void tozeckendorf(SV* svn)
  PREINIT:
    UV n;
  PPCODE:
    if (_validate_and_set(&n, aTHX_ svn, IFLAG_POS)) {
      char *str = to_zeckendorf(n);
      XPUSHs(sv_2mortal(newSVpv(str, 0)));
      Safefree(str);
      XSRETURN(1);
    } else {
      _vcallsub_with_pp("tozeckendorf");
    }
    return;

void fromzeckendorf(IN char* str)
  PREINIT:
    int status;
  PPCODE:
    status = validate_zeckendorf(str);
    if (status == 0) {
      croak("fromzeckendorf takes a binary string as input");
    } else if (status == -1) {
      croak("fromzeckendorf binary input not in canonical Zeckendorf form");
    } else if (status == 1) {
      XSRETURN_UV(from_zeckendorf(str));
    } else {
      _vcallsub_with_pp("fromzeckendorf");
      objectify_result(aTHX_ 0, ST(0));
    }
    return;

void
lastfor()
  PREINIT:
    dMY_CXT;
  PPCODE:
    /* printf("last for with count = %u\n", MY_CXT.forcount); */
    if (MY_CXT.forcount == 0) croak("lastfor called outside a loop");
    MY_CXT.forexit = 1;
    /* In some ideal world this would also act like a last */
    return;

#define START_FORCOUNT \
    do { \
      oldforloop = ++MY_CXT.forcount; \
      oldforexit = MY_CXT.forexit; \
      forexit = &MY_CXT.forexit; \
      *forexit = 0; \
    } while(0)

#define CHECK_FORCOUNT \
    if (*forexit) break;

#define END_FORCOUNT \
    do { \
      /* Put back outer loop's exit request, if any. */ \
      *forexit = oldforexit; \
      /* Ensure loops are nested and not woven. */ \
      if (MY_CXT.forcount-- != oldforloop) croak("for loop mismatch"); \
    } while (0)

#define DECL_FORCOUNT \
    uint16_t oldforloop; \
    char     oldforexit; \
    char    *forexit

void
forprimes (SV* block, IN SV* svbeg, IN SV* svend = 0)
  PROTOTYPE: &$;$
  PREINIT:
    GV *gv;
    HV *stash;
    SV* svarg;
    CV *cv;
    unsigned char* segment;
    UV beg, end, seg_base, seg_low, seg_high;
    DECL_FORCOUNT;
    dMY_CXT;
  PPCODE:
    cv = sv_2cv(block, &stash, &gv, 0);
    if (cv == Nullcv)
      croak("Not a subroutine reference");

    if (!_validate_and_set(&beg, aTHX_ svbeg, IFLAG_POS) ||
        (svend && !_validate_and_set(&end, aTHX_ svend, IFLAG_POS))) {
      _vcallsubn(aTHX_ G_VOID|G_DISCARD, VCALL_ROOT, "_generic_forprimes", items, 0);
      return;
    }
    if (!svend) { end = beg; beg = 2; }

    START_FORCOUNT;
    SAVESPTR(GvSV(PL_defgv));
    svarg = newSVuv(beg);
    GvSV(PL_defgv) = svarg;
    /* Handle early part */
#if USE_MULTICALL
    if (!CvISXSUB(cv) && beg <= end) {
      dMULTICALL;
      I32 gimme = G_VOID;
      PUSH_MULTICALL(cv);
      if (beg < 6) {
        beg = (beg <= 2) ? 2 : (beg <= 3) ? 3 : 5;
        for ( ; beg < 6 && beg <= end; beg += 1+(beg>2) ) {
          CHECK_FORCOUNT;
          sv_setuv(svarg, beg);
          { ENTER; MULTICALL; LEAVE; }
        }
      }
      if (beg <= end) {
       if (
#if BITS_PER_WORD == 64
          (beg >= UVCONST(     100000000000000) && end-beg <    100000) ||
          (beg >= UVCONST(      10000000000000) && end-beg <     40000) ||
          (beg >= UVCONST(       1000000000000) && end-beg <     17000) ||
#endif
          ((end-beg) < 500) ) {     /* MULTICALL next prime */
        for (beg = next_prime(beg-1); beg <= end && beg != 0; beg = next_prime(beg)) {
          CHECK_FORCOUNT;
          sv_setuv(svarg, beg);
          { ENTER; MULTICALL; LEAVE; }
        }
       } else {                      /* MULTICALL segment sieve */
        void* ctx = start_segment_primes(beg, end, &segment);
        while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
          int crossuv = (seg_high > IV_MAX) && !SvIsUV(svarg);
          START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
            CHECK_FORCOUNT;
            /* sv_setuv(svarg, p); */
            if      (SvTYPE(svarg) != SVt_IV) { sv_setuv(svarg, p);            }
            else if (crossuv && p > IV_MAX)   { sv_setuv(svarg, p); crossuv=0; }
            else                              { SvUV_set(svarg, p);            }
            { ENTER; MULTICALL; LEAVE; }
          END_DO_FOR_EACH_SIEVE_PRIME
          CHECK_FORCOUNT;
        }
        end_segment_primes(ctx);
       }
      }
      FIX_MULTICALL_REFCOUNT;
      POP_MULTICALL;
    }
    else
#endif
    {
      if (beg < 6) {
        beg = (beg <= 2) ? 2 : (beg <= 3) ? 3 : 5;
        for ( ; beg < 6 && beg <= end; beg += 1+(beg>2) ) {
          sv_setuv(svarg, beg);
          PUSHMARK(SP); PUTBACK; call_sv((SV*)cv, G_VOID|G_DISCARD); SPAGAIN;
          CHECK_FORCOUNT;
        }
      }
      if (beg <= end) {               /* NO-MULTICALL segment sieve */
        void* ctx = start_segment_primes(beg, end, &segment);
        while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
          START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
            CHECK_FORCOUNT;
            sv_setuv(svarg, p);
            PUSHMARK(SP); PUTBACK; call_sv((SV*)cv, G_VOID|G_DISCARD); SPAGAIN;
          END_DO_FOR_EACH_SIEVE_PRIME
          CHECK_FORCOUNT;
        }
        end_segment_primes(ctx);
      }
    }
    SvREFCNT_dec(svarg);
    END_FORCOUNT;

#define FORCOMPTEST(ix,n) \
  ( (ix==1) || (ix==0 && n&1) )

void
foroddcomposites (SV* block, IN SV* svbeg, IN SV* svend = 0)
  ALIAS:
    forcomposites = 1
  PROTOTYPE: &$;$
  PREINIT:
    UV beg, end;
    GV *gv;
    HV *stash;
    SV* svarg;  /* We use svarg to prevent clobbering $_ outside the block */
    CV *cv;
    DECL_FORCOUNT;
    dMY_CXT;
  PPCODE:
    cv = sv_2cv(block, &stash, &gv, 0);
    if (cv == Nullcv)
      croak("Not a subroutine reference");

    if (!_validate_and_set(&beg, aTHX_ svbeg, IFLAG_POS) ||
        (svend && !_validate_and_set(&end, aTHX_ svend, IFLAG_POS))) {
      _vcallsubn(aTHX_ G_VOID|G_DISCARD, VCALL_ROOT,
         (ix == 0) ? "_generic_foroddcomposites"
       :             "_generic_forcomposites", items, 0);
      return;
    }
    if (!svend) { end = beg; beg = ix ? 4 : 9; }

    START_FORCOUNT;
    SAVESPTR(GvSV(PL_defgv));
    svarg = newSVuv(0);
    GvSV(PL_defgv) = svarg;
#if USE_MULTICALL
    if (!CvISXSUB(cv) && end >= beg) {
      unsigned char* segment;
      UV seg_base, seg_low, seg_high, c, cbeg, cend, cinc, prevprime, nextprime;
      void* ctx;
      dMULTICALL;
      I32 gimme = G_VOID;
      PUSH_MULTICALL(cv);
      if (beg >= MPU_MAX_PRIME ||
#if BITS_PER_WORD == 64
          (beg >= UVCONST(     100000000000000) && end-beg <    120000) ||
          (beg >= UVCONST(      10000000000000) && end-beg <     50000) ||
          (beg >= UVCONST(       1000000000000) && end-beg <     20000) ||
#endif
          end-beg < 1000 ) {
        beg = (beg <= 4) ? 3 : beg-1;
        nextprime = next_prime(beg);
        while (beg++ < end) {
          if      (beg == nextprime)
            nextprime = next_prime(beg);
          else if (FORCOMPTEST(ix,beg)) {
            sv_setuv(svarg, beg);
            { ENTER; MULTICALL; LEAVE; }
          }
          CHECK_FORCOUNT;
        }
      } else {
        if (!ix) {
          if (beg < 8)  beg = 8;
        } else if (beg <= 4) { /* sieve starts at 7, so handle this here */
          sv_setuv(svarg, 4);
          { ENTER; MULTICALL; LEAVE; }
          beg = 6;
        }
        /* Find the two primes that bound their interval. */
        /* beg must be < max_prime, and end >= max_prime is special. */
        prevprime = prev_prime(beg);
        nextprime = (end >= MPU_MAX_PRIME) ? MPU_MAX_PRIME : next_prime(end);
        ctx = start_segment_primes(beg, nextprime, &segment);
        while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
          int crossuv = (seg_high > IV_MAX) && !SvIsUV(svarg);
          START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
            cbeg = prevprime+1;
            if (cbeg < beg)
              cbeg = beg - (ix == 0 && (beg % 2));
            prevprime = p;
            cend = prevprime-1;  if (cend > end) cend = end;
            /* If ix=0, skip evens by starting 1 farther and skipping by 2 */
            cinc = 1 + (ix==0);
            for (c = cbeg + (ix==0); c <= cend; c += cinc) {
              CHECK_FORCOUNT;
              if      (SvTYPE(svarg) != SVt_IV) { sv_setuv(svarg,c); }
              else if (crossuv && c > IV_MAX)   { sv_setuv(svarg,c); crossuv=0;}
              else                              { SvUV_set(svarg,c); }
              { ENTER; MULTICALL; LEAVE; }
            }
          END_DO_FOR_EACH_SIEVE_PRIME
        }
        end_segment_primes(ctx);
        if (end > nextprime)   /* Complete the case where end > max_prime */
          while (nextprime++ < end)
            if (FORCOMPTEST(ix,nextprime)) {
              CHECK_FORCOUNT;
              sv_setuv(svarg, nextprime);
              { ENTER; MULTICALL; LEAVE; }
            }
      }
      FIX_MULTICALL_REFCOUNT;
      POP_MULTICALL;
    }
    else
#endif
    if (beg <= end) {
      beg = (beg <= 4) ? 3 : beg-1;
      while (beg++ < end) {
        if (FORCOMPTEST(ix,beg) && !is_prob_prime(beg)) {
          sv_setuv(svarg, beg);
          PUSHMARK(SP); PUTBACK; call_sv((SV*)cv, G_VOID|G_DISCARD); SPAGAIN;
          CHECK_FORCOUNT;
        }
      }
    }
    SvREFCNT_dec(svarg);
    END_FORCOUNT;

void
forsemiprimes (SV* block, IN SV* svbeg, IN SV* svend = 0)
  PROTOTYPE: &$;$
  PREINIT:
    UV beg, end;
    GV *gv;
    HV *stash;
    SV* svarg;  /* We use svarg to prevent clobbering $_ outside the block */
    CV *cv;
    DECL_FORCOUNT;
    dMY_CXT;
  PPCODE:
    cv = sv_2cv(block, &stash, &gv, 0);
    if (cv == Nullcv)
      croak("Not a subroutine reference");

    if (!_validate_and_set(&beg, aTHX_ svbeg, IFLAG_POS) ||
        (svend && !_validate_and_set(&end, aTHX_ svend, IFLAG_POS))) {
      _vcallsubn(aTHX_ G_VOID|G_DISCARD, VCALL_ROOT, "_generic_forsemiprimes", items, 0);
      return;
    }
    if (!svend) { end = beg; beg = 4; }

    if (beg < 4) beg = 4;
    if (end > MPU_MAX_SEMI_PRIME) end = MPU_MAX_SEMI_PRIME;

    START_FORCOUNT;
    SAVESPTR(GvSV(PL_defgv));
    svarg = newSVuv(0);
    GvSV(PL_defgv) = svarg;
#if USE_MULTICALL
    if (!CvISXSUB(cv) && end >= beg) {
      UV c, seg_beg, seg_end, *S, count;
      dMULTICALL;
      I32 gimme = G_VOID;
      PUSH_MULTICALL(cv);
      if (beg >= MPU_MAX_SEMI_PRIME ||
#if BITS_PER_WORD == 64
          (beg >= UVCONST(10000000000000000000) && end-beg <  1400000) ||
          (beg >= UVCONST( 1000000000000000000) && end-beg <   950000) ||
          (beg >= UVCONST(  100000000000000000) && end-beg <   440000) ||
          (beg >= UVCONST(   10000000000000000) && end-beg <   240000) ||
          (beg >= UVCONST(    1000000000000000) && end-beg <    65000) ||
          (beg >= UVCONST(     100000000000000) && end-beg <    29000) ||
          (beg >= UVCONST(      10000000000000) && end-beg <    11000) ||
          (beg >= UVCONST(       1000000000000) && end-beg <     5000) ||
#endif
          end-beg < 200 ) {
        for (c = beg; c <= end && c >= beg; c++) {
          if (is_semiprime(c)) {
            sv_setuv(svarg, c);
            { ENTER; MULTICALL; LEAVE; }
          }
          CHECK_FORCOUNT;
        }
      } else {
        while (beg < end) {
          seg_beg = beg;
          seg_end = end;
          if ((seg_end - seg_beg) > 50000000) seg_end = seg_beg + 50000000 - 1;
          count = range_semiprime_sieve(&S, seg_beg, seg_end);
          for (c = 0; c < count; c++) {
            sv_setuv(svarg, S[c]);
            { ENTER; MULTICALL; LEAVE; }
            CHECK_FORCOUNT;
          }
          Safefree(S);
          beg = seg_end+1;
          CHECK_FORCOUNT;
        }
      }
      FIX_MULTICALL_REFCOUNT;
      POP_MULTICALL;
    }
    else
#endif
    if (beg <= end) {
      beg = (beg <= 4) ? 3 : beg-1;
      while (beg++ < end) {
        if (is_semiprime(beg)) {
          sv_setuv(svarg, beg);
          PUSHMARK(SP); PUTBACK; call_sv((SV*)cv, G_VOID|G_DISCARD); SPAGAIN;
          CHECK_FORCOUNT;
        }
      }
    }
    SvREFCNT_dec(svarg);
    END_FORCOUNT;

void
foralmostprimes (SV* block, IN UV k, IN SV* svbeg, IN SV* svend = 0)
  PROTOTYPE: &$$;$
  PREINIT:
    UV c, beg, end, shiftres;
    GV *gv;
    HV *stash;
    SV* svarg;  /* We use svarg to prevent clobbering $_ outside the block */
    CV *cv;
    DECL_FORCOUNT;
    dMY_CXT;
  PPCODE:
    cv = sv_2cv(block, &stash, &gv, 0);
    if (cv == Nullcv)
      croak("Not a subroutine reference");

    if (!_validate_and_set(&beg, aTHX_ svbeg, IFLAG_POS) ||
        (svend && !_validate_and_set(&end, aTHX_ svend, IFLAG_POS))) {
      _vcallsub_with_pp("foralmostprimes");
      return;
    }
    if (!svend) { end = beg; beg = 1; }

    /* If k is over 63 but the beg/end points are UVs, then we're empty. */
    if (k == 0 || k >= BITS_PER_WORD) return;

    if (beg < (UVCONST(1) << k)) beg = UVCONST(1) << k;
    if (end > max_nth_almost_prime(k)) end = max_nth_almost_prime(k);
    if (beg > end) return;

    /* We might be able to reduce the k value. */
    shiftres = 0;
    if (k > MPU_MAX_POW3)
      shiftres = k - MPU_MAX_POW3;
    while ((k-shiftres) > 1 && (end >> shiftres) < ipow(3, k - shiftres))
      shiftres++;
    beg = (beg >> shiftres) + (((beg >> shiftres) << shiftres) < beg);
    end = end >> shiftres;
    k -= shiftres;
    /* k <= 40 (64-bit) or 20 (32-bit). */

    START_FORCOUNT;
    SAVESPTR(GvSV(PL_defgv));
    svarg = newSVuv(0);
    GvSV(PL_defgv) = svarg;
#if USE_MULTICALL
    if (!CvISXSUB(cv) && end >= beg) {
      UV seg_beg, seg_end, *S, count, k3 = ipow(3,k);
      dMULTICALL;
      I32 gimme = G_VOID;
      PUSH_MULTICALL(cv);
      while (beg <= end) {
        /* TODO: Tuning this better would be nice */
        UV ssize = 65536 * 256;
        seg_beg = beg;
        seg_end = end;
        if (k > 12)                    ssize *= 16;
        if (k > 18 || seg_beg >  9*k3) ssize *= 4;
        if (k > 24 || seg_beg > 81*k3) ssize *= 3;
        if ((seg_end - seg_beg) > ssize) seg_end = seg_beg + ssize - 1;
        count = generate_almost_primes(&S, k, seg_beg, seg_end);
        for (c = 0; c < count; c++) {
          sv_setuv(svarg, S[c] << shiftres);
          { ENTER; MULTICALL; LEAVE; }
          CHECK_FORCOUNT;
        }
        Safefree(S);
        if (seg_end == UV_MAX) break;
        beg = seg_end+1;
        CHECK_FORCOUNT;
      }
      FIX_MULTICALL_REFCOUNT;
      POP_MULTICALL;
    }
    else
#endif
    if (beg <= end) {
      for (c = beg; c <= end && c >= beg; c++) {
        if (is_almost_prime(k,c)) {
          sv_setuv(svarg, c << shiftres);
          PUSHMARK(SP); PUTBACK; call_sv((SV*)cv, G_VOID|G_DISCARD); SPAGAIN;
          CHECK_FORCOUNT;
        }
      }
    }
    SvREFCNT_dec(svarg);
    END_FORCOUNT;

void
fordivisors (SV* block, IN SV* svn)
  PROTOTYPE: &$
  PREINIT:
    UV i, n, ndivisors;
    UV *divs;
    GV *gv;
    HV *stash;
    SV* svarg;  /* We use svarg to prevent clobbering $_ outside the block */
    CV *cv;
    DECL_FORCOUNT;
    dMY_CXT;
  PPCODE:
    cv = sv_2cv(block, &stash, &gv, 0);
    if (cv == Nullcv)
      croak("Not a subroutine reference");

    if (!_validate_and_set(&n, aTHX_ svn, IFLAG_POS)) {
      _vcallsubn(aTHX_ G_VOID|G_DISCARD, VCALL_ROOT, "_generic_fordivisors", 2, 0);
      return;
    }

    divs = divisor_list(n, &ndivisors, UV_MAX);

    START_FORCOUNT;
    SAVESPTR(GvSV(PL_defgv));
    svarg = newSVuv(0);
    GvSV(PL_defgv) = svarg;
#if USE_MULTICALL
    if (!CvISXSUB(cv)) {
      dMULTICALL;
      I32 gimme = G_VOID;
      PUSH_MULTICALL(cv);
      for (i = 0; i < ndivisors; i++) {
        sv_setuv(svarg, divs[i]);
        { ENTER; MULTICALL; LEAVE; }
        CHECK_FORCOUNT;
      }
      FIX_MULTICALL_REFCOUNT;
      POP_MULTICALL;
    }
    else
#endif
    {
      for (i = 0; i < ndivisors; i++) {
        sv_setuv(svarg, divs[i]);
        PUSHMARK(SP); PUTBACK; call_sv((SV*)cv, G_VOID|G_DISCARD); SPAGAIN;
        CHECK_FORCOUNT;
      }
    }
    SvREFCNT_dec(svarg);
    Safefree(divs);
    END_FORCOUNT;

void
forpart (SV* block, IN SV* svn, IN SV* svh = 0)
  ALIAS:
    forcomp = 1
  PROTOTYPE: &$;$
  PREINIT:
    UV i, n, amin, amax, nmin, nmax;
    int primeq;
    GV *gv;
    HV *stash;
    CV *cv;
    SV** svals;
    DECL_FORCOUNT;
    dMY_CXT;
  PPCODE:
    cv = sv_2cv(block, &stash, &gv, 0);
    if (cv == Nullcv)
      croak("Not a subroutine reference");
    if (!_validate_and_set(&n, aTHX_ svn, IFLAG_POS)) {
      _vcallsub_with_pp("forpart");
      return;
    }
    if (n > (UV_MAX-2)) croak("forpart argument overflow");

    New(0, svals, n+1, SV*);
    for (i = 0; i <= n; i++) {
      svals[i] = newSVuv(i);
      SvREADONLY_on(svals[i]);
    }

    amin = 1;  amax = n;  nmin = 1;  nmax = n;  primeq = -1;
    if (svh != 0) {
      HV* rhash;
      SV** svp;
      if (!SvROK(svh) || SvTYPE(SvRV(svh)) != SVt_PVHV)
        croak("forpart second argument must be a hash reference");
      rhash = (HV*) SvRV(svh);
      if ((svp = hv_fetchs(rhash, "n", 0)) != NULL)
        { nmin = my_svuv(*svp);  nmax = nmin; }
      if ((svp = hv_fetchs(rhash, "amin", 0)) != NULL) amin = my_svuv(*svp);
      if ((svp = hv_fetchs(rhash, "amax", 0)) != NULL) amax = my_svuv(*svp);
      if ((svp = hv_fetchs(rhash, "nmin", 0)) != NULL) nmin = my_svuv(*svp);
      if ((svp = hv_fetchs(rhash, "nmax", 0)) != NULL) nmax = my_svuv(*svp);
      if ((svp = hv_fetchs(rhash, "prime",0)) != NULL) primeq=my_svuv(*svp);

      if (amin < 1) amin = 1;
      if (amax > n) amax = n;
      if (nmin < 1) nmin = 1;
      if (nmax > n) nmax = n;
      if (primeq != 0 && primeq != -1) primeq = 2;  /* -1, 0, or 2 */
    }

    if (primeq == 2) {
      UV prev =                 prev_prime(amax+1);
      UV next = amin <= 2 ? 2 : next_prime(amin-1);
      if (amin < next)  amin = next;
      if (amax > prev)  amax = prev;
    }

    if (n==0 && nmin <= 1) {
      { PUSHMARK(SP);
        /* Nothing */
        PUTBACK; call_sv((SV*)cv, G_VOID|G_DISCARD); SPAGAIN;
      }
    }
    if (n >= nmin && nmin <= nmax && amin <= amax && nmax > 0 && amax > 0)
    { /* RuleAsc algorithm from Kelleher and O'Sullivan 2009/2014) */
      UV *a, k, x, y, r;
      New(0, a, n+1, UV);
      k = 1;
      a[0] = amin-1;
      a[1] = n-amin+1;
      START_FORCOUNT;
      while (k != 0) {
        x = a[k-1]+1;
        y = a[k]-1;
        k--;
        r = (ix == 0) ? x : 1;
        while (r <= y) {
          a[k++] = x;
          x = r;
          y -= x;
        }
        a[k] = x + y;

        /* ------ length restrictions ------ */
        while (k+1 > nmax) {   /* Skip range if over max size */
          a[k-1] += a[k];
          k--;
        }
        /* Look into: quick skip over nmin range */
        if (k+1 < nmin) {      /* Skip if not over min size */
          if (a[0] >= n-nmin+1 && a[k] > 1) break; /* early exit check */
          continue;
        }

        /* ------ value restrictions ------ */
        if (amin > 1 || amax < n) {
          /* Lexical order allows us to start at amin, and exit early */
          if (a[0] > amax) break;

          if (ix == 0) {  /* value restrictions for partitions */
            if (a[k] > amax) continue;
          } else {  /* restrictions for compositions */
            /* TODO: maybe skip forward? */
            for (i = 0; i <= k; i++)
              if (a[i] < amin || a[i] > amax)
                break;
            if (i <= k) continue;
          }
        }
        if (primeq != -1) {
          for (i = 0; i <= k; i++) if (is_prime(a[i]) != primeq) break;
          if (i <= k) continue;
        }

        PUSHMARK(SP); EXTEND(SP, (IV)k);
        for (i = 0; i <= k; i++) { PUSHs(svals[a[i]]); }
        PUTBACK; call_sv((SV*)cv, G_VOID|G_DISCARD); SPAGAIN;
        CHECK_FORCOUNT;
      }
      Safefree(a);
      END_FORCOUNT;
    }
    for (i = 0; i <= n; i++)
      SvREFCNT_dec(svals[i]);
    Safefree(svals);

void
forcomb (SV* block, IN SV* svn, IN SV* svk = 0)
  ALIAS:
    forperm = 1
    forderange = 2
  PROTOTYPE: &$;$
  PREINIT:
    UV i, n, k, begk, endk;
    GV *gv;
    HV *stash;
    CV *cv;
    SV** svals;
    UV*  cm;
    DECL_FORCOUNT;
    dMY_CXT;
  PPCODE:
    cv = sv_2cv(block, &stash, &gv, 0);
    if (cv == Nullcv)
      croak("Not a subroutine reference");
    if (ix > 0 && svk != 0)
      croak("Too many arguments for forperm");

    if (!_validate_and_set(&n, aTHX_ svn, IFLAG_POS) ||
        (svk && !_validate_and_set(&k, aTHX_ svk, IFLAG_POS))) {
      _vcallsub_with_pp(   (ix == 0) ? "forcomb"
                         : (ix == 1) ? "forperm"
                                     : "forderange" );
      return;
    }

    if (svk == 0) {
      begk = (ix == 0) ? 0 : n;
      endk = n;
    } else {
      begk = endk = k;
      if (begk > n)
        return;
    }

    New(0, svals, n, SV*);
    for (i = 0; i < n; i++) {
      svals[i] = newSVuv(i);
      SvREADONLY_on(svals[i]);
    }
    New(0, cm, endk+1, UV);

    START_FORCOUNT;
#if USE_MULTICALL
    if (!CvISXSUB(cv)) {
      dMULTICALL;
      I32 gimme = G_VOID;
      AV *av = save_ary(PL_defgv);
      AvREAL_off(av);
      PUSH_MULTICALL(cv);
      for (k = begk; k <= endk; k++) {
        _comb_init(cm, k, ix == 2);
        while (1) {
          if (ix < 2 || k != 1) {
            IV j;
            av_extend(av, k-1);
            av_fill(av, k-1);
            for (j = k-1; j >= 0; j--)
              AvARRAY(av)[j] = svals[ cm[k-j-1]-1 ];
            { ENTER; MULTICALL; LEAVE; }
          }
          CHECK_FORCOUNT;
          if (_comb_iterate(cm, k, n, ix)) break;
        }
        CHECK_FORCOUNT;
      }
      FIX_MULTICALL_REFCOUNT;
      POP_MULTICALL;
    } else
#endif
    {
      for (k = begk; k <= endk; k++) {
        _comb_init(cm, k, ix == 2);
        while (1) {
          if (ix < 2 || k != 1) {
            PUSHMARK(SP); EXTEND(SP, ((IV)k));
            for (i = 0; i < k; i++) { PUSHs(svals[ cm[k-i-1]-1 ]); }
            PUTBACK; call_sv((SV*)cv, G_VOID|G_DISCARD); SPAGAIN;
          }
          CHECK_FORCOUNT;
          if (_comb_iterate(cm, k, n, ix)) break;
        }
        CHECK_FORCOUNT;
      }
    }

    Safefree(cm);
    for (i = 0; i < n; i++)
      SvREFCNT_dec(svals[i]);
    Safefree(svals);
    END_FORCOUNT;

void
forsetproduct (SV* block, ...)
  PROTOTYPE: &@
  PREINIT:
    IV narrays, i, *arlen, *arcnt;
    AV **arptr;
    SV **arout;
    GV *gv;
    HV *stash;
    CV *cv;
    DECL_FORCOUNT;
    dMY_CXT;
  PPCODE:
    cv = sv_2cv(block, &stash, &gv, 0);
    if (cv == Nullcv) croak("Not a subroutine reference");

    narrays = items-1;
    if (narrays < 1) return;

    for (i = 1; i <= narrays; i++) {
      SvGETMAGIC(ST(i));
      if ((!SvROK(ST(i))) || (SvTYPE(SvRV(ST(i))) != SVt_PVAV))
        croak("forsetproduct arguments must be array references");
      if (av_len((AV *)SvRV(ST(i))) < 0)
        return;
    }

    Newz(0, arcnt, narrays, IV);
    New(0, arlen, narrays, IV);
    New(0, arptr, narrays, AV*);
    New(0, arout, narrays, SV*);
    for (i = 0; i < narrays; i++) {
      arptr[i] = (AV*) SvRV(ST(i+1));
      arlen[i] = 1 + av_len(arptr[i]);
      arout[i] = AvARRAY(arptr[i])[0];
    }

    START_FORCOUNT;
#if USE_MULTICALL
    if (!CvISXSUB(cv)) {
      dMULTICALL;
      I32 gimme = G_VOID;
      AV *av = save_ary(PL_defgv);
      AvREAL_off(av);
      PUSH_MULTICALL(cv);
      do {
        av_extend(av, narrays-1);
        av_fill(av, narrays-1);
        for (i = narrays-1; i >= 0; i--)  /* Faster to fill backwards */
          AvARRAY(av)[i] = arout[i];
        { ENTER; MULTICALL; LEAVE; }
        CHECK_FORCOUNT;
        for (i = narrays-1; i >= 0; i--) {
          if (++arcnt[i] >= arlen[i])  arcnt[i] = 0;
          arout[i] = AvARRAY(arptr[i])[arcnt[i]];
          if (arcnt[i] > 0)  break;
        }
      } while (i >= 0);
      FIX_MULTICALL_REFCOUNT;
      POP_MULTICALL;
    }
    else
#endif
    do {
      PUSHMARK(SP); EXTEND(SP, narrays);
      for (i = 0; i < narrays; i++) { PUSHs(arout[i]); }
      PUTBACK; call_sv((SV*)cv, G_VOID|G_DISCARD); SPAGAIN;
      CHECK_FORCOUNT;
      for (i = narrays-1; i >= 0; i--) {
        if (++arcnt[i] >= arlen[i])  arcnt[i] = 0;
        arout[i] = AvARRAY(arptr[i])[arcnt[i]];
        if (arcnt[i] > 0)  break;
      }
    } while (i >= 0);
    Safefree(arout);
    Safefree(arptr);
    Safefree(arlen);
    Safefree(arcnt);
    END_FORCOUNT;

void
forfactored (SV* block, IN SV* svbeg, IN SV* svend = 0)
  ALIAS:
    forsquarefree = 1
  PROTOTYPE: &$;$
  PREINIT:
    UV beg, end, n, *factors;
    int i, nfactors, maxfactors;
    factor_range_context_t fctx;
    GV *gv;
    HV *stash;
    SV* svarg;  /* We use svarg to prevent clobbering $_ outside the block */
    CV *cv;
    SV* svals[64];
    DECL_FORCOUNT;
    dMY_CXT;
  PPCODE:
    cv = sv_2cv(block, &stash, &gv, 0);
    if (cv == Nullcv)
      croak("Not a subroutine reference");

    if (!_validate_and_set(&beg, aTHX_ svbeg, IFLAG_POS) ||
        (svend && !_validate_and_set(&end, aTHX_ svend, IFLAG_POS))) {
      _vcallsubn(aTHX_ G_VOID|G_DISCARD, VCALL_ROOT, (ix == 0) ? "_generic_forfactored" : "_generic_forsquarefree", items, 0);
      return;
    }
    if (!svend) { end = beg; beg = 1; }
    if (beg > end) return;

    for (maxfactors = 0, n = end >> 1;  n;  n >>= 1)
      maxfactors++;
    for (i = 0; i < maxfactors; i++) {
      svals[i] = newSVuv(UV_MAX);
      SvREADONLY_on(svals[i]);
    }

    SAVESPTR(GvSV(PL_defgv));
    svarg = newSVuv(0);
    GvSV(PL_defgv) = svarg;
    START_FORCOUNT;
    if (beg <= 1) {
      PUSHMARK(SP);
      sv_setuv(svarg, 1);
      PUTBACK; call_sv((SV*)cv, G_VOID|G_DISCARD); SPAGAIN;
      beg = 2;
    }
    fctx = factor_range_init(beg, end, ix);
#if USE_MULTICALL
    if (!CvISXSUB(cv)) {
      dMULTICALL;
      I32 gimme = G_VOID;
      AV *av = save_ary(PL_defgv);
      AvREAL_off(av);
      PUSH_MULTICALL(cv);
      for (n = 0; n < end-beg+1; n++) {
        CHECK_FORCOUNT;
        nfactors = factor_range_next(&fctx);
        if (nfactors > 0) {
          sv_setuv(svarg, fctx.n);
          factors = fctx.factors;
          av_extend(av, nfactors-1);
          av_fill(av, nfactors-1);
          for (i = nfactors-1; i >= 0; i--) {
            SV* sv = svals[i];
            SvREADONLY_off(sv);
            sv_setuv(sv, factors[i]);
            SvREADONLY_on(sv);
            AvARRAY(av)[i] = sv;
          }
          { ENTER; MULTICALL; LEAVE; }
        }
      }
      FIX_MULTICALL_REFCOUNT;
      POP_MULTICALL;
    }
    else
#endif
    for (n = 0; n < end-beg+1; n++) {
      CHECK_FORCOUNT;
      nfactors = factor_range_next(&fctx);
      if (nfactors > 0) {
        PUSHMARK(SP); EXTEND(SP, nfactors);
        sv_setuv(svarg, fctx.n);
        factors = fctx.factors;
        for (i = 0; i < nfactors; i++) {
          SV* sv = svals[i];
          SvREADONLY_off(sv);
          sv_setuv(sv, factors[i]);
          SvREADONLY_on(sv);
          PUSHs(sv);
        }
        PUTBACK; call_sv((SV*)cv, G_VOID|G_DISCARD); SPAGAIN;
      }
    }
    factor_range_destroy(&fctx);
    SvREFCNT_dec(svarg);
    for (i = 0; i < maxfactors; i++)
      SvREFCNT_dec(svals[i]);
    END_FORCOUNT;

void
vecreduce(SV* block, ...)
PROTOTYPE: &@
CODE:
{   /* This is basically reduce from List::Util.  Try to maintain compat. */
    SV *ret = sv_newmortal();
    int i;
    GV *agv,*bgv,*gv;
    HV *stash;
    SV **args = &PL_stack_base[ax];
    CV *cv = sv_2cv(block, &stash, &gv, 0);

    if (cv == Nullcv) croak("Not a subroutine reference");
    if (items <= 1) XSRETURN_UNDEF;

    agv = gv_fetchpv("a", GV_ADD, SVt_PV);
    bgv = gv_fetchpv("b", GV_ADD, SVt_PV);
    SAVESPTR(GvSV(agv));
    SAVESPTR(GvSV(bgv));
    GvSV(agv) = ret;
    SvSetMagicSV(ret, args[1]);
#ifdef dMULTICALL
    if (!CvISXSUB(cv)) {
      dMULTICALL;
      I32 gimme = G_SCALAR;
      PUSH_MULTICALL(cv);
      for (i = 2; i < items; i++) {
        GvSV(bgv) = args[i];
        { ENTER; MULTICALL; LEAVE; }
        SvSetMagicSV(ret, *PL_stack_sp);
      }
      FIX_MULTICALL_REFCOUNT;
      POP_MULTICALL;
    }
    else
#endif
    {
      for (i = 2; i < items; i++) {
        dSP;
        GvSV(bgv) = args[i];
        PUSHMARK(SP);
        call_sv((SV*)cv, G_SCALAR);
        SvSetMagicSV(ret, *PL_stack_sp);
      }
    }
    ST(0) = ret;
    XSRETURN(1);
}

void
vecnone(SV* block, ...)
ALIAS:
    vecall    = 1
    vecany    = 2
    vecnotall = 3
    vecfirst  = 4
    vecfirstidx = 6
PROTOTYPE: &@
PPCODE:
{   /* This is very similar to List::Util.  Try to maintain compat. */
    int ret_true = !(ix & 2); /* return true at end of loop for none/all; false for any/notall */
    int invert   =  (ix & 1); /* invert block test for all/notall */
    int index;
    GV *gv;
    HV *stash;
    SV **args = &PL_stack_base[ax];
    CV *cv    = sv_2cv(block, &stash, &gv, 0);

    if (cv == Nullcv) croak("Not a subroutine reference");

    SAVESPTR(GvSV(PL_defgv));
#ifdef dMULTICALL
    if (!CvISXSUB(cv)) {
      dMULTICALL;
      I32 gimme = G_SCALAR;

      PUSH_MULTICALL(cv);
      for (index = 1; index < items; index++) {
        GvSV(PL_defgv) = args[index];
        { ENTER; MULTICALL; LEAVE; }
        if (SvTRUEx(*PL_stack_sp) ^ invert)
          break;
      }
      FIX_MULTICALL_REFCOUNT;
      POP_MULTICALL;
    }
    else
#endif
    {
      for (index = 1; index < items; index++) {
        dSP;
        GvSV(PL_defgv) = args[index];
        PUSHMARK(SP);
        call_sv((SV*)cv, G_SCALAR);
        if (SvTRUEx(*PL_stack_sp) ^ invert)
          break;
      }
    }

    if (ix == 4) {
      if (index == items)
        XSRETURN_UNDEF;
      ST(0) = ST(index);
      XSRETURN(1);
    }
    if (ix == 6) {
      if (index == items)
        XSRETURN_IV(-1);
      XSRETURN_UV(index-1);
    }

    if (index != items)           /* We exited the loop early */
      ret_true = !ret_true;

    if (ret_true)  XSRETURN_YES;
    else           XSRETURN_NO;
}

void vecuniq(...)
  PROTOTYPE: @
  PREINIT:
    iset_t s;
    int i, status, retvals;
    UV n;
    unsigned long count;
  PPCODE:
    retvals = (GIMME_V != G_SCALAR);
    s = iset_create(items);
    for (status = 1, i = 0; i < items; i++) {
      status = _validate_and_set(&n, aTHX_ ST(i), IFLAG_ANY);
      if (status == 0) break;
      if (iset_add(&s, n, status) == 0)
        continue;
      if (iset_sign(s) == 0) { status = 0; break; }
      if (retvals)
        PUSHs(sv_2mortal(NEWSVINT(status,n)));
    }
    count = iset_size(s);
    iset_destroy(&s);
    if (status != 0 && retvals) {
      XSRETURN(count);
    } else if (status != 0) {
      ST(0) = sv_2mortal(newSVuv(count));
      XSRETURN(1);
    } else {
      /* This is 100% from List::MoreUtils::XS by Parseval and Rehsack */
      I32 i;
      IV count = 0, seen_undef = 0;
      HV *hv = newHV();
      SV **args = &PL_stack_base[ax];
      SV *tmp = sv_newmortal();
      sv_2mortal(newRV_noinc((SV*)hv));

      if (GIMME_V == G_SCALAR) { /* don't build return list if not needed */
        for (i = 0; i < items; i++) {
          SvGETMAGIC(args[i]);
          if (SvOK(args[i])) {
            sv_setsv_nomg(tmp, args[i]);
            if (!hv_exists_ent(hv, tmp, 0)) {
              ++count;
              hv_store_ent(hv, tmp, &PL_sv_yes, 0);
            }
          } else if (0 == seen_undef++)
            ++count;
        }
        ST(0) = sv_2mortal(newSVuv(count));
        XSRETURN(1);
      }
      /* list context: populate SP with mortal copies */
      for (i = 0; i < items; i++) {
        SvGETMAGIC(args[i]);
        if (SvOK(args[i])) {
          SvSetSV_nosteal(tmp, args[i]);
          if (!hv_exists_ent(hv, tmp, 0)) {
            args[count++] = args[i];
            hv_store_ent(hv, tmp, &PL_sv_yes, 0);
          }
        } else if (0 == seen_undef++)
          args[count++] = args[i];
      }
      XSRETURN(count);
    }

#ifdef FACTORING_HARNESSES
void
factor_test_harness1(...)
  PROTOTYPE: @
  PPCODE:
    /* Pass in a big array of numbers, we factor them in a timed loop */
    {
      UV res, factors[MPU_MAX_FACTORS+1], *comp;
      struct timeval gstart, gstop;
      double t_time;
      int i, j, k, correct, nf, num = items;

      //num = (items > 100000) ? 100000 : items;
      New(0, comp, num, UV);
      for (i = 0; i < num; i++)
        comp[i] = my_svuv(ST(i));
      gettimeofday(&gstart, NULL);
      for (j = 0; j < 1; j++) {
        correct = 0;
        for (i = 0; i < num; i++) {
          nf = factor(comp[i], factors);
          //nf = squfof_factor(comp[i], factors, 140000);
          //nf = pbrent_factor(comp[i], factors, 500000, 1);
          //nf = holf_factor(comp[i], factors, 1000000);
          //nf = lehman_factor(comp[i], factors, 1);
          //nf = lehman_factor(comp[i], factors, 0);  if (nf < 2) nf=pbrent_factor(comp[i], factors, 500000, 3);
          //nf = factor63(comp[i], factors);
          //nf = pminus1_factor(comp[i], factors, 1000,10000);
          //nf = prho_factor(comp[i], factors, 10000);
          if (nf >= 2) {
            for (res = factors[0], k = 1; k < nf; k++)
              res *= factors[k];
            if (res == comp[i])
              correct++;
          }
          //printf("%lu:",comp[i]);for(k=0;k<nf;k++)printf(" %lu",factors[k]);printf("\n");
        }
      }
      gettimeofday(&gstop, NULL);
      t_time = my_difftime(&gstart, &gstop);
      Safefree(comp);
      printf("factoring got %d of %d correct in %2.2f sec\n", correct, num, t_time);
      printf("percent correct = %.3f\n", 100.0*(double)correct / (double)num);
      printf("average time per input = %1.4f ms\n", 1000 * t_time / (double)num);
      XSRETURN_UV(correct);
    }

void
factor_test_harness2(IN int count, IN int bits = 63)
  PPCODE:
    /* We'll factor <count> <bits>-bit numbers */
    {
      UV factors[MPU_MAX_FACTORS+1], exponents[MPU_MAX_FACTORS+1];
      FILE *fid = 0;  // fopen("results.txt", "w");
      uint64_t n, state = 28953;
      int i, nfactors, totfactors = 0;
      /* Use Knuth MMIX -- simple and no worse than Chacha20 for this */
      for (i = 0; i < count; i++) {
        state = 6364136223846793005ULL * state + 1442695040888963407ULL;
        n = state >> (64-bits);
        nfactors = factor_exp(n, factors, exponents);
        if (fid) fprintf(fid, "%llu has %d factors\n", n, nfactors);
        totfactors += nfactors;
      }
      if (fid) fclose(fid);
      XSRETURN_IV(totfactors);
    }

void
factor_test_harness3(IN UV start, IN UV end)
  PPCODE:
    /* We'll factor <count> <bits>-bit numbers */
    {
      UV totf = 0, i, factors[MPU_MAX_FACTORS];
      for (i = start; i < end; i++) {
        totf += factor(i, factors);
      }
      XSRETURN_UV(totf);
    }

#endif
