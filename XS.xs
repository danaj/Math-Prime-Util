
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

#include "ptypes.h"
#include "cache.h"
#include "sieve.h"
#include "sieve_cluster.h"
#define FUNC_gcd_ui 1
#define FUNC_isqrt 1
#define FUNC_ipow 1
#include "util.h"
#include "primality.h"
#include "factor.h"
#include "lehmer.h"
#include "lmo.h"
#include "aks.h"
#include "constants.h"
#include "mulmod.h"

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
#if !defined cxinc && ( (PERL_VERSION == 8 && PERL_SUBVERSION >= 2) || (PERL_VERSION == 10 && PERL_SUBVERSION <= 1) )
# define cxinc() Perl_cxinc(aTHX)
#endif

#if PERL_VERSION < 17 || (PERL_VERSION == 17 && PERL_SUBVERSION < 7)
#  define SvREFCNT_dec_NN(sv)    SvREFCNT_dec(sv)
#endif

#if BITS_PER_WORD == 32
  static const unsigned int uvmax_maxlen = 10;
  static const unsigned int ivmax_maxlen = 10;
  static const char uvmax_str[] = "4294967295";
  static const char ivmax_str[] = "2147483648";
#else
  static const unsigned int uvmax_maxlen = 20;
  static const unsigned int ivmax_maxlen = 19;
  static const char uvmax_str[] = "18446744073709551615";
  static const char ivmax_str[] =  "9223372036854775808";
#endif

#define MY_CXT_KEY "Math::Prime::Util::API_guts"
typedef struct {
  SV* const_int[11];   /* -1, 0, 1, ..., 9 */
  HV* MPUroot;
  HV* MPUGMP;
  HV* MPUPP;
} my_cxt_t;

START_MY_CXT

/* Is this a pedantically valid integer?
 * Croaks if undefined or invalid.
 * Returns 0 if it is an object or a string too large for a UV.
 * Returns 1 if it is good to process by XS.
 */
static int _validate_int(pTHX_ SV* n, int negok)
{
  const char* maxstr;
  char* ptr;
  STRLEN i, len, maxlen;
  int ret, isbignum = 0, isneg = 0;

  /* TODO: magic, grok_number, etc. */
  if (SVNUMTEST(n)) { /* If defined as number, use it */
    if (SvIsUV(n) || SvIVX(n) >= 0)  return 1; /* The normal case */
    if (negok)  return -1;
    else croak("Parameter '%" SVf "' must be a positive integer", n);
  }
  if (sv_isobject(n)) {
    const char *hvname = HvNAME_get(SvSTASH(SvRV(n)));
    if (hvname == 0)
      return 0;
    if (strEQ(hvname, "Math::BigInt") || strEQ(hvname, "Math::BigFloat") ||
        strEQ(hvname, "Math::GMPz")   || strEQ(hvname, "Math::GMP") ||
        strEQ(hvname, "Math::Pari") )
      isbignum = 1;
    else
      return 0;
  }
  /* Without being very careful, don't process magic variables here */
  if (SvGAMAGIC(n) && !isbignum) return 0;
  if (!SvOK(n))  croak("Parameter must be defined");
  ptr = SvPV_nomg(n, len);             /* Includes stringifying bigints */
  if (len == 0 || ptr == 0)  croak("Parameter must be a positive integer");
  if (ptr[0] == '-' && negok) {
    isneg = 1; ptr++; len--;           /* Read negative sign */
  } else if (ptr[0] == '+') {
    ptr++; len--;                      /* Allow a single plus sign */
  }
  if (len == 0 || !isDIGIT(ptr[0]))
    croak("Parameter '%" SVf "' must be a positive integer", n);
  while (len > 0 && *ptr == '0')       /* Strip all leading zeros */
    { ptr++; len--; }
  if (len > uvmax_maxlen)              /* Huge number, don't even look at it */
    return 0;
  for (i = 0; i < len; i++)            /* Ensure all characters are digits */
    if (!isDIGIT(ptr[i]))
      croak("Parameter '%" SVf "' must be a positive integer", n);
  if (isneg == 1)                      /* Negative number (ignore overflow) */
    return -1;
  ret    = isneg ? -1           : 1;
  maxlen = isneg ? ivmax_maxlen : uvmax_maxlen;
  maxstr = isneg ? ivmax_str    : uvmax_str;
  if (len < maxlen)                    /* Valid small integer */
    return ret;
  for (i = 0; i < maxlen; i++) {       /* Check if in range */
    if (ptr[i] < maxstr[i]) return ret;
    if (ptr[i] > maxstr[i]) return 0;
  }
  return ret;                          /* value = UV_MAX/UV_MIN.  That's ok */
}

#define VCALL_ROOT 0x0
#define VCALL_PP 0x1
#define VCALL_GMP 0x2
/* Call a Perl sub to handle work for us. */
static int _vcallsubn(pTHX_ I32 flags, I32 stashflags, const char* name, int nargs, int version)
{
    GV* gv = NULL;
    dMY_CXT;
    Size_t namelen = strlen(name);
    /* If given a GMP function, and GMP enabled, and function exists, use it. */
    int use_gmp = stashflags & VCALL_GMP && _XS_get_callgmp() && _XS_get_callgmp() >= version;
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

/* In my testing, this constant return works fine with threads, but to be
 * correct (see perlxs) one has to make a context, store separate copies in
 * each one, then retrieve them from a struct using a hash index.  This
 * defeats the purpose if only done once. */
#define RETURN_NPARITY(ret) \
  do { int r_ = ret; \
       dMY_CXT; \
       if (r_ >= -1 && r_ <= 9) { ST(0) = MY_CXT.const_int[r_+1]; XSRETURN(1); } \
       else                     { XSRETURN_IV(r_);                      } \
  } while (0)
#define PUSH_NPARITY(ret) \
  do { int r_ = ret; \
       if (r_ >= -1 && r_ <= 9) { PUSHs( MY_CXT.const_int[r_+1] );       } \
       else                     { PUSHs(sv_2mortal(newSViv(r_))); } \
  } while (0)

#define OBJECTIFY_RESULT(input, output) \
  if (!sv_isobject(output)) { \
    SV* resptr = output; \
    const char *iname = (input && sv_isobject(input)) \
                      ? HvNAME_get(SvSTASH(SvRV(input))) : 0; \
    if (iname == 0 || strEQ(iname, "Math::BigInt")) { \
      (void)_vcallsubn(aTHX_ G_SCALAR, VCALL_ROOT, "_to_bigint", 1, 0); \
    } else if (iname == 0 || strEQ(iname, "Math::GMPz")) { \
      (void)_vcallsubn(aTHX_ G_SCALAR, VCALL_ROOT, "_to_gmpz", 1, 0); \
    } else if (iname == 0 || strEQ(iname, "Math::GMP")) { \
      (void)_vcallsubn(aTHX_ G_SCALAR, VCALL_ROOT, "_to_gmp", 1, 0); \
    } else { /* Return it as: ref(input)->new(result) */ \
      dSP;  ENTER;  PUSHMARK(SP); \
      XPUSHs(sv_2mortal(newSVpv(iname, 0)));  XPUSHs(resptr); \
      PUTBACK;  call_method("new", G_SCALAR);  LEAVE; \
    } \
  }

static SV* sv_to_bigint(pTHX_ SV* r) {
  dSP;  ENTER;  PUSHMARK(SP);
  XPUSHs(r);
  PUTBACK;
  call_pv("Math::Prime::Util::_to_bigint", G_SCALAR);
  SPAGAIN;
  r = POPs;
  PUTBACK; LEAVE;
  return r;
}

static int arrayref_to_int_array(pTHX_ UV** ret, AV* av, int base)
{
  int len, i;
  UV *r, carry = 0;
  if (SvTYPE((SV*)av) != SVt_PVAV)
    croak("fromdigits first argument must be a string or array reference");
  len = 1 + av_len(av);
  New(0, r, len, UV);
  for (i = len-1; i >= 0; i--) {
    SV** psvd = av_fetch(av, i, 0);
    if (_validate_int(aTHX_ *psvd, 1) != 1) break;
    r[i] = my_svuv(*psvd) + carry;
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

static UV negmod(IV a, UV n) {
  UV negamod = ((UV)(-a)) % n;
  return (negamod == 0) ? 0 : n-negamod;
}


MODULE = Math::Prime::Util	PACKAGE = Math::Prime::Util

PROTOTYPES: ENABLE

BOOT:
{
    SV * sv = newSViv(BITS_PER_WORD);
    HV * stash = gv_stashpv("Math::Prime::Util", TRUE);
    newCONSTSUB(stash, "_XS_prime_maxbits", sv);
    { int i;
      MY_CXT_INIT;
      MY_CXT.MPUroot = stash;
      for (i = 0; i <= 10; i++) {
        MY_CXT.const_int[i] = newSViv(i-1);
        SvREADONLY_on(MY_CXT.const_int[i]);
      }
      MY_CXT.MPUGMP = gv_stashpv("Math::Prime::Util::GMP", TRUE);
      MY_CXT.MPUPP = gv_stashpv("Math::Prime::Util::PP", TRUE);
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
    for (i = 0; i <= 10; i++) {
      MY_CXT.const_int[i] = newSViv(i-1);
      SvREADONLY_on(MY_CXT.const_int[i]);
    }
    MY_CXT.MPUroot = gv_stashpv("Math::Prime::Util", TRUE);
    MY_CXT.MPUGMP = gv_stashpv("Math::Prime::Util::GMP", TRUE);
    MY_CXT.MPUPP = gv_stashpv("Math::Prime::Util::PP", TRUE);
  }
  return; /* skip implicit PUTBACK, returning @_ to caller, more efficient*/

#endif

void
END(...)
PREINIT:
  dMY_CXT;
  int i;
PPCODE:
  for (i = 0; i <= 10; i++) {
    SV * const sv = MY_CXT.const_int[i];
    MY_CXT.const_int[i] = NULL;
    SvREFCNT_dec_NN(sv);
  } /* stashes are owned by stash tree, no refcount on them in MY_CXT */
  MY_CXT.MPUroot = NULL;
  MY_CXT.MPUGMP = NULL;
  MY_CXT.MPUPP = NULL;
  _prime_memfreeall();
  return; /* skip implicit PUTBACK, returning @_ to caller, more efficient*/

void
prime_memfree()
  ALIAS:
    _XS_get_verbose = 1
    _XS_get_callgmp = 2
    _get_prime_cache_size = 3
  PREINIT:
    UV ret;
  PPCODE:
    switch (ix) {
      case 0:  prime_memfree(); goto return_nothing;
      case 1:  ret = _XS_get_verbose(); break;
      case 2:  ret = _XS_get_callgmp(); break;
      case 3:
      default: ret = get_prime_cache(0,0); break;
    }
    XSRETURN_UV(ret);
    return_nothing:

void
prime_precalc(IN UV n)
  ALIAS:
    _XS_set_verbose = 1
    _XS_set_callgmp = 2
  PPCODE:
    PUTBACK; /* SP is never used again, the 3 next func calls are tailcall
    friendly since this XSUB has nothing to do after the 3 calls return */
    switch (ix) {
      case 0:  prime_precalc(n);    break;
      case 1:  _XS_set_verbose(n);  break;
      default: _XS_set_callgmp(n);  break;
    }
    return; /* skip implicit PUTBACK */

void
prime_count(IN SV* svlo, ...)
  ALIAS:
    _XS_segment_pi = 1
    twin_prime_count = 2
    ramanujan_prime_count = 3
    sum_primes = 4
    print_primes = 5
  PREINIT:
    int lostatus, histatus;
    UV lo, hi;
  PPCODE:
    lostatus = _validate_int(aTHX_ svlo, 0);
    histatus = (items == 1 || _validate_int(aTHX_ ST(1), 0));
    if (lostatus == 1 && histatus == 1) {
      UV count = 0;
      if (items == 1) {
        lo = 2;
        hi = my_svuv(svlo);
      } else {
        lo = my_svuv(svlo);
        hi = my_svuv(ST(1));
      }
      if (lo <= hi) {
        if (ix == 2) {
          count = twin_prime_count(lo, hi);
        } else if (ix == 3) {
          count = ramanujan_prime_count(lo, hi);
        } else if (ix == 4) {
          lostatus = sum_primes(lo, hi, &count);
        } else if (ix == 5) {
          int fd = (items < 3) ? fileno(stdout) : my_sviv(ST(2));
          print_primes(lo, hi, fd);
          XSRETURN_EMPTY;
        } else if (ix == 1 || (hi / (hi-lo+1)) > 100) {
          count = _XS_prime_count(lo, hi);
        } else {
          count = _XS_LMO_pi(hi);
          if (lo > 2)
            count -= _XS_LMO_pi(lo-1);
        }
      }
      if (lostatus == 1) XSRETURN_UV(count);
    }
    switch (ix) {
      case 0:
      case 1: _vcallsubn(aTHX_ GIMME_V, VCALL_ROOT, "_generic_prime_count", items, 0); break;
      case 2:_vcallsub_with_pp("twin_prime_count");  break;
      case 3:_vcallsub_with_pp("ramanujan_prime_count");  break;
      case 4:_vcallsub_with_pp("sum_primes");  break;
      case 5:
      default:_vcallsub_with_pp("print_primes");  break;
    }
    return; /* skip implicit PUTBACK */

UV
_XS_LMO_pi(IN UV n)
  ALIAS:
    _XS_legendre_pi = 1
    _XS_meissel_pi = 2
    _XS_lehmer_pi = 3
    _XS_LMOS_pi = 4
  PREINIT:
    UV ret;
  CODE:
    switch (ix) {
      case 0: ret = _XS_LMO_pi(n); break;
      case 1: ret = _XS_legendre_pi(n); break;
      case 2: ret = _XS_meissel_pi(n); break;
      case 3: ret = _XS_lehmer_pi(n); break;
      default:ret = _XS_LMOS_pi(n); break;
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
    segment_twin_primes = 4
    _ramanujan_primes = 5
  PREINIT:
    AV* av;
  PPCODE:
    av = newAV();
    {
      SV * retsv = sv_2mortal(newRV_noinc( (SV*) av ));
      PUSHs(retsv);
      PUTBACK;
      SP = NULL; /* never use SP again, poison */
    }
    if ((low <= 2) && (high >= 2) && ix != 4) { av_push(av, newSVuv( 2 )); }
    if ((low <= 3) && (high >= 3) && ix != 5) { av_push(av, newSVuv( 3 )); }
    if ((low <= 5) && (high >= 5) && ix != 5) { av_push(av, newSVuv( 5 )); }
    if (low < 7)  low = 7;
    if (low <= high) {
      if (ix == 4) high += 2;
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
      } else if (ix == 3 || ix == 4) {        /* Segment */
        unsigned char* segment;
        UV seg_base, seg_low, seg_high, lastp = 0;
        void* ctx = start_segment_primes(low, high, &segment);
        while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
          START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
            if (ix == 3)            av_push(av,newSVuv( p ));
            else if (lastp+2 == p)  av_push(av,newSVuv( lastp ));
            lastp = p;
          END_DO_FOR_EACH_SIEVE_PRIME
        }
        end_segment_primes(ctx);
      } else if (ix == 5) {                   /* Ramanujan primes */
        UV i, beg, end, *L;
        L = ramanujan_primes(&beg, &end, low, high);
        if (L && end >= beg)
          for (i = beg; i <= end; i++)
            av_push(av,newSVuv(L[i]));
        Safefree(L);
      }
    }
    return; /* skip implicit PUTBACK */

void
sieve_range(IN SV* svn, IN UV width, IN UV depth)
  PREINIT:
    int status;
  PPCODE:
    status = _validate_int(aTHX_ svn, 0);
    if (status == 1) {
      /* TODO: actually sieve */
      UV factors[MPU_MAX_FACTORS+1], i, n = my_svuv(svn);
      if (depth == 0) depth = 1; /* Trial factor takes 0 to means sqrt(n) */
      if ( (n + width) < n) {  /* Overflow */
        status = 0;
      } else if (depth <= 100) { /* trial division for each value */
        for (i = (n<2)?2-n:0; i < width; i++)
          if (trial_factor(n+i, factors, depth) < 2)
            XPUSHs(sv_2mortal(newSVuv( i )));
      } else {                   /* small trial + factor for each value */
        for (i = (n<2)?2-n:0; i < width; i++)
          if (trial_factor(n+i, factors, 100) < 2)
            if (factor(n+i,factors) < 2 || factors[0] > depth)
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
    UV i, cval, nprimes, *list;
    int lostatus, histatus, done;
  PPCODE:
    nc = items-1;
    if (items > 100) croak("sieve_prime_cluster: too many entries");
    cl[0] = 0;
    for (i = 1; i < nc; i++) {
      if (!_validate_int(aTHX_ ST(1+i), 0)) croak("sieve_prime_cluster: cluster values must be standard integers");
      cval = my_svuv(ST(1+i));
      if (cval & 1) croak("sieve_prime_cluster: values must be even");
      if (cval > 2147483647UL) croak("sieve_prime_cluster: values must be 31-bit");
      if (cval <= cl[i-1]) croak("sieve_prime_cluster: values must be increasing");
      cl[i] = cval;
    }
    lostatus = _validate_int(aTHX_ svlo, 1);
    histatus = _validate_int(aTHX_ svhi, 1);
    done = 0;
    if (lostatus == 1 && histatus == 1) {
      UV low = my_svuv(svlo);
      UV high = my_svuv(svhi);
      list = sieve_cluster(low, high, nc, cl, &nprimes);
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

void
trial_factor(IN UV n, ...)
  ALIAS:
    fermat_factor = 1
    holf_factor = 2
    squfof_factor = 3
    prho_factor = 4
    pplus1_factor = 5
    pbrent_factor = 6
    pminus1_factor = 7
    ecm_factor = 8
  PREINIT:
    UV arg1, arg2;
    static const UV default_arg1[] =
       {0,     64000000, 8000000, 4000000, 4000000, 200, 4000000, 1000000};
     /* Trial, Fermat,   Holf,    SQUFOF,  PRHO,    P+1, Brent,    P-1 */
  PPCODE:
    if (n == 0)  XSRETURN_UV(0);
    if (ix == 8) {  /* We don't have an ecm_factor, call PP. */
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
    else if (_XS_is_prime(n)) { XPUSHs(sv_2mortal(newSVuv( n ))); }
    else {
      UV factors[MPU_MAX_FACTORS+1];
      int i, nfactors = 0;
      switch (ix) {
        case 0:  nfactors = trial_factor  (n, factors, arg1);  break;
        case 1:  nfactors = fermat_factor (n, factors, arg1);  break;
        case 2:  nfactors = holf_factor   (n, factors, arg1);  break;
        case 3:  nfactors = squfof_factor (n, factors, arg1);  break;
        case 4:  nfactors = prho_factor   (n, factors, arg1);  break;
        case 5:  nfactors = pplus1_factor (n, factors, arg1);  break;
        case 6:  if (items < 3) arg2 = 1;
                 nfactors = pbrent_factor (n, factors, arg1, arg2);  break;
        case 7:
        default: if (items < 3) arg2 = 10*arg1;
                 nfactors = pminus1_factor(n, factors, arg1, arg2);  break;
      }
      EXTEND(SP, nfactors);
      for (i = 0; i < nfactors; i++)
        PUSHs(sv_2mortal(newSVuv( factors[i] )));
    }

void
is_strong_pseudoprime(IN SV* svn, ...)
  ALIAS:
    is_pseudoprime = 1
    is_euler_pseudoprime = 2
  PREINIT:
    int c, status = 1;
  PPCODE:
    if (items < 2)
      croak("No bases given to is_strong_pseudoprime");
    /* Check all arguments */
    for (c = 0; c < items && status == 1; c++)
      if (_validate_int(aTHX_ ST(c), 0) != 1)
        status = 0;
    if (status == 1) {
      UV n = my_svuv(svn);
      int b, ret = 1;
      if        (n < 4) {                        /* 0,1 composite; 2,3 prime */
        ret = (n >= 2);
      } else if (ix == 1) {                      /* Fermat test */
        for (c = 1; c < items && ret == 1; c++)
          ret = is_pseudoprime(n, my_svuv(ST(c)));
      } else if (ix == 2) {                      /* Euler test */
        for (c = 1; c < items && ret == 1; c++)
          ret = is_euler_pseudoprime(n, my_svuv(ST(c)));
      } else if ((n % 2) == 0) {                 /* evens composite */
         ret = 0;
      } else {
        UV bases[32];
        for (c = 1; c < items && ret == 1; ) {
          for (b = 0; b < 32 && c < items; c++)
            bases[b++] = my_svuv(ST(c));
          ret = miller_rabin(n, bases, b);
        }
      }
      RETURN_NPARITY(ret);
    }
    switch (ix) {
      case 0: _vcallsub_with_gmp(0.00,"is_strong_pseudoprime"); break;
      case 1: _vcallsub_with_gmp(0.20,"is_pseudoprime"); break;
      case 2:
      default:_vcallsub_with_gmp(0.00,"is_euler_pseudoprime");  break;
    }
    return; /* skip implicit PUTBACK */

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
      status = _validate_int(aTHX_ ST(0), 2);
      if (status != 0 && items > 1) {
        sign = status;
        ret = my_svuv(ST(0));
        for (i = 1; i < items; i++) {
          status = _validate_int(aTHX_ ST(i), 2);
          if (status == 0) break;
          n = my_svuv(ST(i));
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
        status = _validate_int(aTHX_ ST(i), 2);
        if (status == 0) break;
        n = my_svuv(ST(i));
        if (status == 1) {
          hi += (n > (UV_MAX - lo));
        } else {
          if (UV_MAX-n == (UV)IV_MAX) { status = 0; break; }  /* IV Overflow */
          hi -= ((UV_MAX-n) >= lo);
        }
        lo += n;
      }
      if (status != 0 && hi == -1 && lo > IV_MAX)  XSRETURN_IV((IV)lo);
      if (status != 0 && hi != 0) {
        /* If status != 0 then the 128-bit result is:
         *   result = (hi << 64) + lo    (hi is signed, lo is unsigned)
         * Call a Perl routine that will take these and make a bigint for us. */
        XPUSHs(sv_2mortal(newSViv( hi )));
        XPUSHs(sv_2mortal(newSVuv( lo )));
        _vcallsub_with_pp("_from_128");
        XSRETURN(1);
      }
      if (hi != 0) status = 0;  /* Overflow */
      ret = lo;
    } else if (ix == 5) {
      int sign = 1;
      ret = 1;
      for (i = 0; i < items; i++) {
        status = _validate_int(aTHX_ ST(i), 2);
        if (status == 0) break;
        n = (status == 1) ? my_svuv(ST(i)) : (UV)-my_sviv(ST(i));
        if (ret > 0 && n > UV_MAX/ret) { status = 0; break; }
        sign *= status;
        ret *= n;
      }
      if (sign == -1 && status != 0) {
        if (ret <= (UV)IV_MAX)  XSRETURN_IV(-(IV)ret);
        else                    status = 0;
      }
    } else {
      /* For each arg, while valid input, validate+gcd/lcm.  Shortcut stop. */
      if (ix == 0) { ret = 0; nullv = 1; }
      else         { ret = (items == 0) ? 0 : 1; nullv = 0; }
      for (i = 0; i < items && ret != nullv && status != 0; i++) {
        status = _validate_int(aTHX_ ST(i), 2);
        if (status == 0)
          break;
        n = status * my_svuv(ST(i));  /* n = abs(arg) */
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
    UV i = 0;
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
    } else if (_validate_int(aTHX_ svm, 0)) {
      UV mask = my_svuv(svm);
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
chinese(...)
  PROTOTYPE: @
  PREINIT:
    int i, status;
    UV* an;
    UV ret;
  PPCODE:
    status = 1;
    New(0, an, 2*items, UV);
    ret = 0;
    for (i = 0; i < items; i++) {
      AV* av;
      SV** psva;
      SV** psvn;
      if (!SvROK(ST(i)) || SvTYPE(SvRV(ST(i))) != SVt_PVAV || av_len((AV*)SvRV(ST(i))) != 1)
        croak("chinese arguments are two-element array references");
      av = (AV*) SvRV(ST(i));
      psva = av_fetch(av, 0, 0);
      psvn = av_fetch(av, 1, 0);
      if (psva == 0 || psvn == 0 || _validate_int(aTHX_ *psva, 1) != 1 || !_validate_int(aTHX_ *psvn, 0)) {
        status = 0;
        break;
      }
      an[i+0]     = my_svuv(*psva);
      an[i+items] = my_svuv(*psvn);
    }
    if (status)
      ret = chinese(an, an+items, items, &status);
    Safefree(an);
    if (status == -1) XSRETURN_UNDEF;
    if (status)       XSRETURN_UV(ret);
    _vcallsub_with_pp("chinese");
    return; /* skip implicit PUTBACK */

void
lucas_sequence(...)
  ALIAS:
    lucasu = 1
    lucasv = 2
  PREINIT:
    UV U, V, Qk;
  PPCODE:
    if (ix == 1 || ix == 2) {
      if (items != 3) croak("lucasu: P, Q, k");
      if (_validate_int(aTHX_ ST(0), 1) && _validate_int(aTHX_ ST(1), 1) &&
          _validate_int(aTHX_ ST(2), 0)) {
        IV P = my_sviv(ST(0));
        IV Q = my_sviv(ST(1));
        UV k = my_svuv(ST(2));
        IV ret;
        int ok = (ix == 1) ? lucasu(&ret, P, Q, k) : lucasv(&ret, P, Q, k);
        if (ok) XSRETURN_IV(ret);
      }
      _vcallsub_with_gmp(0.29,(ix==1) ? "lucasu" : "lucasv");
      OBJECTIFY_RESULT(ST(2), ST(0));
      return;
    }
    if (items != 4) croak("lucas_sequence: n, P, Q, k");
    if (_validate_int(aTHX_ ST(0), 0) && _validate_int(aTHX_ ST(1), 1) &&
        _validate_int(aTHX_ ST(2), 1) && _validate_int(aTHX_ ST(3), 0)) {
      lucas_seq(&U, &V, &Qk,
                my_svuv(ST(0)), my_sviv(ST(1)), my_sviv(ST(2)), my_svuv(ST(3)));
      PUSHs(sv_2mortal(newSVuv( U )));  /* 4 args in, 3 out, no EXTEND needed */
      PUSHs(sv_2mortal(newSVuv( V )));
      PUSHs(sv_2mortal(newSVuv( Qk )));
    } else {
      _vcallsubn(aTHX_ GIMME_V, VCALL_PP, "lucas_sequence", items, 0);
      return;
    }

void
is_prime(IN SV* svn, ...)
  ALIAS:
    is_prob_prime = 1
    is_provable_prime = 2
    is_bpsw_prime = 3
    is_aks_prime = 4
    is_lucas_pseudoprime = 5
    is_strong_lucas_pseudoprime = 6
    is_extra_strong_lucas_pseudoprime = 7
    is_frobenius_pseudoprime = 8
    is_frobenius_underwood_pseudoprime = 9
    is_frobenius_khashin_pseudoprime = 10
    is_perrin_pseudoprime = 11
    is_catalan_pseudoprime = 12
    is_almost_extra_strong_lucas_pseudoprime = 13
    is_euler_plumb_pseudoprime = 14
    is_ramanujan_prime = 15
    is_square_free = 16
    is_carmichael = 17
    is_quasi_carmichael = 18
    is_primitive_root = 19
    is_mersenne_prime = 20
    is_power = 21
    is_prime_power = 22
    logint = 23
    rootint = 24
  PREINIT:
    int status, astatus;
  PPCODE:
    status = _validate_int(aTHX_ svn, 1);
    astatus = (items == 1 || ix==22) ? 1 : _validate_int(aTHX_ ST(1), 1);
    if (status != 0 && astatus != 0) {
      int ret = 0;
      UV n = my_svuv(svn);
      if (status == 1 && astatus == 1 && ix < 21) {
        UV a = (items == 1) ? 0 : my_svuv(ST(1));
        switch (ix) {
          case 0:
          case 1:
          case 2:  ret = _XS_is_prime(n); break;
          case 3:  ret = BPSW(n); break;
          case 4:  ret = is_aks_prime(n); break;
          case 5:  ret = is_lucas_pseudoprime(n, 0); break;
          case 6:  ret = is_lucas_pseudoprime(n, 1); break;
          case 7:  ret = is_lucas_pseudoprime(n, 2); break;
          case 8:  {
                     /* IV P = 1, Q = -1; */ /* Fibonacci polynomial */
                     IV P = 0, Q = 0;        /* Q=2,P=least odd s.t. (D|n)=-1 */
                     if (items == 3) { P = my_sviv(ST(1)); Q = my_sviv(ST(2)); }
                     else if (items != 1) croak("is_frobenius_pseudoprime takes P,Q");
                     ret = is_frobenius_pseudoprime(n, P, Q);
                   } break;
          case 9:  ret = is_frobenius_underwood_pseudoprime(n); break;
          case 10: ret = is_frobenius_khashin_pseudoprime(n); break;
          case 11: ret = is_perrin_pseudoprime(n,a); break;
          case 12: ret = is_catalan_pseudoprime(n); break;
          case 13: ret = is_almost_extra_strong_lucas_pseudoprime
                         (n, (items == 1) ? 1 : a); break;
          case 14: ret = is_euler_plumb_pseudoprime(n); break;
          case 15: ret = is_ramanujan_prime(n); break;
          case 16: ret = is_square_free(n); break;
          case 17: ret = is_carmichael(n); break;
          case 18: ret = is_quasi_carmichael(n); break;
          case 19: ret = is_primitive_root(n,a,0); break;
          case 20:
          default: ret = is_mersenne_prime(n);
                   if (ret == -1) status = 0;
                   break;
        }
        if (status != 0 && astatus == 1) RETURN_NPARITY(ret);
      } else if (ix == 21) {
        UV a = (items == 1) ? 0 : my_svuv(ST(1));
        if (status == -1) {
          if (n == UV_MAX || (UV_MAX-n) == (UV)IV_MAX)  status = 0;
          else                                          n = (UV) -my_sviv(svn);
        }
        if (status == 1 || (status == -1 && (a == 0 || a & 1))) {
          ret = is_power(n, a);
          if (status == -1 && a == 0) {
            ret >>= valuation(ret,2);
            if (ret == 1) ret = 0;
          }
          if (ret && items == 3) {
            UV root = rootof(n, a ? a : (UV)ret);
            if (!SvROK(ST(2))) croak("is_power: third argument not a scalar reference");
            if (status == 1) sv_setuv(SvRV(ST(2)),  root);
            else             sv_setiv(SvRV(ST(2)), -root);
          }
        }
        if (status != 0 && astatus != 0) RETURN_NPARITY(ret);
      } else if (ix == 22) {
        if (status == 1) {
          UV root;
          ret = primepower(n, &root);
          if (ret && items > 1) {
            if (!SvROK(ST(1))) croak("is_prime_power: second argument not a scalar reference");
            sv_setuv(SvRV(ST(1)), root);
          }
        }
        if (status != 0) RETURN_NPARITY(ret);
      } else if (ix == 23) {
        UV e, a = (items == 1) ? 0 : my_svuv(ST(1));
        if (status != 1 || n <= 0)   croak("logint: n must be > 0");
        if (items == 1)              croak("logint: missing base");
        if (astatus != 1 || a <= 1)  croak("logint: base must be > 1");
        e = logint(n, a);
        if (items == 3) {
          if (!SvROK(ST(2))) croak("logint: third argument not a scalar reference");
          sv_setuv(SvRV(ST(2)), ipow(a,e));
        }
        XSRETURN_UV(e);
      } else if (ix == 24) {
        UV r, a = (items == 1) ? 0 : my_svuv(ST(1));
        if (items == 1)              croak("rootint: missing exponent");
        if (astatus != 1 || a == 0)  croak("rootint: k must be > 0");
        if (status == -1)            croak("rootint: n must be >= 0");
        if (status == 1) {
          r = rootof(n, a);
          if (items == 3) {
            if (!SvROK(ST(2))) croak("rootint: third argument not a scalar reference");
            sv_setuv(SvRV(ST(2)), ipow(r,a));
          }
          XSRETURN_UV(r);
        }
      }
    }
    switch (ix) {
      case 0: _vcallsub_with_gmp(0.01,"is_prime");       break;
      case 1: _vcallsub_with_gmp(0.01,"is_prob_prime");  break;
      case 2: _vcallsub_with_gmp(0.04,"is_provable_prime");  break;
      case 3: _vcallsub_with_gmp(0.17,"is_bpsw_prime");  break;
      case 4: _vcallsub_with_gmp(0.16,"is_aks_prime"); break;
      case 5: _vcallsub_with_gmp(0.01,"is_lucas_pseudoprime"); break;
      case 6: _vcallsub_with_gmp(0.01,"is_strong_lucas_pseudoprime"); break;
      case 7: _vcallsub_with_gmp(0.01,"is_extra_strong_lucas_pseudoprime"); break;
      case 8: _vcallsub_with_gmp(0.24,"is_frobenius_pseudoprime"); break;
      case 9: _vcallsub_with_gmp(0.13,"is_frobenius_underwood_pseudoprime"); break;
      case 10:_vcallsub_with_gmp(0.30,"is_frobenius_khashin_pseudoprime"); break;
      case 11:_vcallsub_with_gmp(
                (items==1 || (astatus==0 && my_svuv(ST(1)) == 0)) ? 0.20 : 0.40,
                "is_perrin_pseudoprime"); break;
      case 12:_vcallsub_with_gmp(0.00,"is_catalan_pseudoprime"); break;
      case 13:_vcallsub_with_gmp(0.13,"is_almost_extra_strong_lucas_pseudoprime"); break;
      case 14:_vcallsub_with_gmp(0.39,"is_euler_plumb_pseudoprime"); break;
      case 15:_vcallsub_with_gmp(0.00,"is_ramanujan_prime"); break;
      case 16:_vcallsub_with_gmp(0.00,"is_square_free"); break;
      case 17:_vcallsub_with_gmp(0.00,"is_carmichael"); break;
      case 18:_vcallsub_with_gmp(0.00,"is_quasi_carmichael"); break;
      case 19:_vcallsub_with_gmp(0.36,"is_primitive_root"); break;
      case 20:_vcallsub_with_gmp(0.28,"is_mersenne_prime"); break;
      case 21:if (items != 3) {
                _vcallsub_with_gmp((status == -1) ? 0.28 : 0.19, "is_power");
              } else {
                _vcallsub_with_pp("is_power");
              }
              break;
      case 22:(void)_vcallsubn(aTHX_ G_SCALAR, (items == 1) ? (VCALL_GMP|VCALL_PP) : (VCALL_PP), "is_prime_power", items, 40); break;
      case 23:_vcallsub_with_gmp(0.00,"logint"); break;
      case 24:
      default:(void)_vcallsubn(aTHX_ G_SCALAR, (items == 2) ? (VCALL_GMP|VCALL_PP) : (VCALL_PP), "rootint", items, 40); break;
    }
    return; /* skip implicit PUTBACK */

void
next_prime(IN SV* svn)
  ALIAS:
    prev_prime = 1
    nth_prime = 2
    nth_prime_upper = 3
    nth_prime_lower = 4
    nth_prime_approx = 5
    nth_twin_prime = 6
    nth_twin_prime_approx = 7
    nth_ramanujan_prime = 8
    prime_count_upper = 9
    prime_count_lower = 10
    prime_count_approx = 11
    twin_prime_count_approx = 12
  PPCODE:
    if (_validate_int(aTHX_ svn, 0)) {
      UV n = my_svuv(svn);
      if ( (n >= MPU_MAX_PRIME     && ix == 0) ||
           (n >= MPU_MAX_PRIME_IDX && (ix==2 || ix==3 || ix==4 || ix==5)) ||
           (n >= MPU_MAX_TWIN_PRIME_IDX && (ix==6 || ix==7)) ||
           (n >= MPU_MAX_RMJN_PRIME_IDX && (ix==8)) ) {
        /* Out of range.  Fall through to Perl. */
      } else {
        UV ret;
        /* Prev prime of 2 or less should return undef */
        if (ix == 1 && n < 3) XSRETURN_UNDEF;
        /* nth_prime(0) and similar should return undef */
        if (n == 0 && (ix >= 2 && ix <= 8)) XSRETURN_UNDEF;
        switch (ix) {
          case 0: ret = next_prime(n);  break;
          case 1: ret = (n < 3) ? 0 : prev_prime(n);  break;
          case 2: ret = nth_prime(n); break;
          case 3: ret = nth_prime_upper(n); break;
          case 4: ret = nth_prime_lower(n); break;
          case 5: ret = nth_prime_approx(n); break;
          case 6: ret = nth_twin_prime(n); break;
          case 7: ret = nth_twin_prime_approx(n); break;
          case 8: ret = nth_ramanujan_prime(n); break;
          case 9: ret = prime_count_upper(n); break;
          case 10:ret = prime_count_lower(n); break;
          case 11:ret = prime_count_approx(n); break;
          case 12:
          default:ret = twin_prime_count_approx(n); break;
        }
        XSRETURN_UV(ret);
      }
    }
    if ((ix == 0 || ix == 1) && _XS_get_callgmp() && PERL_REVISION >= 5 && PERL_VERSION > 8) {
      _vcallsub_with_gmp(0.01, ix ? "prev_prime" : "next_prime");
      OBJECTIFY_RESULT(svn, ST(0));
      return;
    }
    switch (ix) {
      case 0:  _vcallsub_with_pp("next_prime");         break;
      case 1:  _vcallsub_with_pp("prev_prime");         break;
      case 2:  _vcallsub_with_pp("nth_prime");          break;
      case 3:  _vcallsub_with_pp("nth_prime_upper");    break;
      case 4:  _vcallsub_with_pp("nth_prime_lower");    break;
      case 5:  _vcallsub_with_pp("nth_prime_approx");   break;
      case 6:  _vcallsub_with_pp("nth_twin_prime");     break;
      case 7:  _vcallsub_with_pp("nth_twin_prime_approx"); break;
      case 8:  _vcallsub_with_pp("nth_ramanujan_prime"); break;
      case 9:  _vcallsub_with_pp("prime_count_upper");  break;
      case 10: _vcallsub_with_pp("prime_count_lower");  break;
      case 11: _vcallsub_with_pp("prime_count_approx"); break;
      case 12:
      default: _vcallsub_with_pp("twin_prime_count_approx"); break;
    }
    return; /* skip implicit PUTBACK */

void Pi(IN UV digits = 0)
  PREINIT:
#ifdef USE_QUADMATH
    NV pival = 3.141592653589793238462643383279502884197169Q;
    UV mantsize = FLT128_MANT_DIG / 3.322 - 1;
#else
    NV pival = 3.141592653589793238462643383279502884197169L;
    UV mantsize = DBL_MANT_DIG / 3.322;   /* Let long doubles go to BF */
#endif
  PPCODE:
    if (digits == 0) {
      XSRETURN_NV( pival );
    } else if (digits <= mantsize && digits <= 40) {
      char t[40+2];
      NV pi;
      (void)sprintf(t, "%.*"NVff, (int)(digits-1), pival);
#ifdef USE_QUADMATH
      pi = strtoflt128(t, NULL);
#elif defined(USE_LONG_DOUBLE) && defined(HAS_LONG_DOUBLE)
      pi = strtold(t, NULL);
#else
      pi = strtod(t, NULL);
#endif
      XSRETURN_NV( pi );
    } else {
      _vcallsub_with_pp("Pi");
      return;
    }

void
_pidigits(IN int digits)
  PPCODE:
    if (digits == 1) {
      XSRETURN_UV(3);
    } else {
      char *out;
      IV  *a;
      IV b, c, d, e, f, g, i,  d4, d3, d2, d1;

      digits++;   /* For rounding */
      b = d = e = g = i = 0;  f = 10000;
      c = 14*(digits/4 + 2);
      New(0, a, c, IV);
      New(0, out, digits+5+1, char);
      *out++ = '3';  /* We'll turn "31415..." into "3.1415..." */
      for (b = 0; b < c; b++)  a[b] = 20000000;

      while ((b = c -= 14) > 0 && i < digits) {
        d = e = d % f;
        while (--b > 0) {
          d = d * b + a[b];
          g = (b << 1) - 1;
          a[b] = (d % g) * f;
          d /= g;
        }
        /* sprintf(out+i, "%04d", e+d/f);   i += 4; */
        d4 = e+d/f;
        if (d4 > 9999) {
          d4 -= 10000;
          out[i-1]++;
          for (b=i-1; out[b] == '0'+1; b--) { out[b]='0'; out[b-1]++; }
        }
        d3 = d4/10;  d2 = d3/10;  d1 = d2/10;
        out[i++] = '0' + d1;
        out[i++] = '0' + d2-d1*10;
        out[i++] = '0' + d3-d2*10;
        out[i++] = '0' + d4-d3*10;
      }
      Safefree(a);
      if (out[digits-1] >= '5') out[digits-2]++;  /* Round */
      for (i = digits-2; out[i] == '9'+1; i--)    /* Keep rounding */
        { out[i] = '0';  out[i-1]++; }
      digits--;  /* Undo the extra digit we used for rounding */
      out[digits] = '\0';
      *out-- = '.';
      XPUSHs(sv_2mortal(newSVpvn(out, digits+1)));
      Safefree(out);
    }


void
factor(IN SV* svn)
  ALIAS:
    factor_exp = 1
    divisors = 2
  PREINIT:
    U32 gimme_v;
    int status, i, nfactors;
  PPCODE:
    gimme_v = GIMME_V;
    status = _validate_int(aTHX_ svn, 0);
    if (status == 1) {
      UV factors[MPU_MAX_FACTORS+1];
      UV exponents[MPU_MAX_FACTORS+1];
      UV n = my_svuv(svn);
      if (gimme_v == G_SCALAR) {
        switch (ix) {
          case 0:  nfactors = factor(n, factors);        break;
          case 1:  nfactors = factor_exp(n, factors, 0); break;
          default: nfactors = divisor_sum(n, 0);         break;
        }
        PUSHs(sv_2mortal(newSVuv( nfactors )));
      } else if (gimme_v == G_ARRAY) {
        switch (ix) {
          case 0:  nfactors = factor(n, factors);
                   EXTEND(SP, nfactors);
                   for (i = 0; i < nfactors; i++)
                     PUSHs(sv_2mortal(newSVuv( factors[i] )));
                   break;
          case 1:  nfactors = factor_exp(n, factors, exponents);
                   /* if (n == 1)  XSRETURN_EMPTY; */
                   EXTEND(SP, nfactors);
                   for (i = 0; i < nfactors; i++) {
                     AV* av = newAV();
                     av_push(av, newSVuv(factors[i]));
                     av_push(av, newSVuv(exponents[i]));
                     PUSHs( sv_2mortal(newRV_noinc( (SV*) av )) );
                   }
                   break;
          default: {
                     UV ndivisors;
                     UV* divs = _divisor_list(n, &ndivisors);
                     EXTEND(SP, (IV)ndivisors);
                     for (i = 0; (UV)i < ndivisors; i++)
                       PUSHs(sv_2mortal(newSVuv( divs[i] )));
                     Safefree(divs);
                   }
                   break;
        }
      }
    } else {
      switch (ix) {
        case 0:  _vcallsubn(aTHX_ gimme_v, VCALL_ROOT, "_generic_factor", 1, 0);     break;
        case 1:  _vcallsubn(aTHX_ gimme_v, VCALL_ROOT, "_generic_factor_exp", 1, 0); break;
        default: _vcallsubn(aTHX_ gimme_v, VCALL_GMP|VCALL_PP, "divisors", 1, 0);   break;
      }
      return; /* skip implicit PUTBACK */
    }

void
divisor_sum(IN SV* svn, ...)
  PREINIT:
    SV* svk;
    int nstatus, kstatus;
  PPCODE:
    svk = (items > 1) ? ST(1) : 0;
    nstatus = _validate_int(aTHX_ svn, 0);
    kstatus = (items == 1 || (SvIOK(svk) && SvIV(svk) >= 0))  ?  1  :  0;
    /* The above doesn't understand small bigints */
    if (nstatus == 1 && kstatus == 0 && SvROK(svk) && (sv_isa(svk, "Math::BigInt") || sv_isa(svk, "Math::GMP") || sv_isa(svk, "Math::GMPz")))
      kstatus = _validate_int(aTHX_ svk, 0);
    if (nstatus == 1 && kstatus == 1) {
      UV n = my_svuv(svn);
      UV k = (items > 1) ? my_svuv(svk) : 1;
      UV sigma = divisor_sum(n, k);
      if (sigma != 0)  XSRETURN_UV(sigma);   /* sigma 0 means overflow */
    }
    _vcallsub_with_gmp(0.00,"divisor_sum");
    return; /* skip implicit PUTBACK */

void
znorder(IN SV* sva, IN SV* svn)
  ALIAS:
    binomial = 1
    jordan_totient = 2
    ramanujan_sum = 3
    legendre_phi = 4
  PREINIT:
    int astatus, nstatus;
  PPCODE:
    astatus = _validate_int(aTHX_ sva, (ix==1) ? 2 : 0);
    nstatus = _validate_int(aTHX_ svn, (ix==1) ? 2 : 0);
    if (astatus != 0 && nstatus != 0) {
      UV a = my_svuv(sva);
      UV n = my_svuv(svn);
      UV ret;
      switch (ix) {
        case 0:  ret = znorder(a, n);
                 break;
        case 1:  if ( (astatus == 1 && (nstatus == -1 || n > a)) ||
                      (astatus ==-1 && (nstatus == -1 && n > a)) )
                   { ret = 0; break; }
                 if (nstatus == -1)
                   n = a - n; /* n<0,k<=n:  (-1)^(n-k) * binomial(-k-1,n-k) */
                 if (astatus == -1) {
                   ret = binomial( -my_sviv(sva)+n-1, n );
                   if (ret > 0 && ret <= (UV)IV_MAX)
                     XSRETURN_IV( (IV)ret * ((n&1) ? -1 : 1) );
                   goto overflow;
                 } else {
                   ret = binomial(a, n);
                   if (ret == 0)
                     goto overflow;
                 }
                 break;
        case 2:  ret = jordan_totient(a, n);
                 if (ret == 0 && n > 1)
                   goto overflow;
                 break;
        case 3:  if (a < 1 || n < 1) XSRETURN_IV(0);
                 {
                   UV g = a / gcd_ui(a,n);
                   int m = moebius(g);
                   if (m == 0 || a == g) RETURN_NPARITY(m);
                   XSRETURN_IV( m * (totient(a) / totient(g)) );
                 }
                 break;
        case 4:
        default: ret = legendre_phi(a, n);
                 break;
      }
      if (ret == 0 && ix == 0)  XSRETURN_UNDEF;  /* not defined */
      XSRETURN_UV(ret);
    }
    overflow:
    switch (ix) {
      case 0:  _vcallsub_with_pp("znorder");  break;
      case 1:  _vcallsub_with_pp("binomial");  break;
      case 2:  _vcallsub_with_pp("jordan_totient");  break;
      case 3:  _vcallsub_with_pp("ramanujan_sum");  break;
      case 4:
      default: _vcallsub_with_pp("legendre_phi"); break;
    }
    return; /* skip implicit PUTBACK */

void
znlog(IN SV* sva, IN SV* svg, IN SV* svp)
  ALIAS:
    addmod = 1
    mulmod = 2
    divmod = 3
    powmod = 4
  PREINIT:
    int astatus, gstatus, pstatus, retundef;
    UV ret;
  PPCODE:
    astatus = _validate_int(aTHX_ sva, (ix == 0) ? 0 : 1);
    gstatus = _validate_int(aTHX_ svg, (ix == 0) ? 0 : 1);
    pstatus = _validate_int(aTHX_ svp, 0);
    if (astatus != 0 && gstatus != 0 && pstatus == 1) {
      UV a, g, p = my_svuv(svp);
      if (p <= 1) XSRETURN_UV(0);
      ret = 0;
      retundef = 0;
      a = (astatus == 1) ? my_svuv(sva) : negmod(my_sviv(sva), p);
      g = (gstatus == 1) ? my_svuv(svg) : negmod(my_sviv(svg), p);
      if (a >= p) a %= p;
      if (g >= p && ix != 4) g %= p;
      switch (ix) {
        case 0: ret = znlog(a, g, p);
                if (ret == 0 && a > 1) retundef = 1;
                if (ret == 0 && (a == 0 || g == 0)) retundef = 1;
                break;
        case 1: ret = addmod(a, g, p); break;
        case 2: ret = mulmod(a, g, p); break;
        case 3: g = modinverse(g, p);
                if (g == 0) retundef = 1;
                else        ret = mulmod(a, g, p);
                break;
        case 4:
        default:if (a == 0) {
                  ret = (g == 0);
                  retundef = (gstatus == -1);
                } else {
                  if (gstatus == -1) {
                    a = modinverse(a, p);
                    if (a == 0) retundef = 1;
                    else        g = -my_sviv(svg);
                  }
                  ret = powmod(a, g, p);
                }
                break;
      }
      if (retundef) XSRETURN_UNDEF;
      XSRETURN_UV(ret);
    }
    switch (ix) {
      case 0: _vcallsub_with_gmp(0.00,"znlog"); break;
      case 1: _vcallsub_with_gmp(0.36,"addmod"); break;
      case 2: _vcallsub_with_gmp(0.36,"mulmod"); break;
      case 3: _vcallsub_with_gmp(0.36,"divmod"); break;
      case 4:
      default:_vcallsub_with_gmp(0.36,"powmod"); break;
    }
    OBJECTIFY_RESULT(svp, ST(items-1));
    return; /* skip implicit PUTBACK */

void
kronecker(IN SV* sva, IN SV* svb)
  ALIAS:
    valuation = 1
    invmod = 2
    sqrtmod = 3
  PREINIT:
    int astatus, bstatus, abpositive, abnegative;
  PPCODE:
    astatus = _validate_int(aTHX_ sva, 2);
    bstatus = _validate_int(aTHX_ svb, 2);
    if (astatus != 0 && bstatus != 0) {
      if (ix == 0) {
        /* Are both a and b positive? */
        abpositive = astatus == 1 && bstatus == 1;
        /* Will both fit in IVs?  We should use a bitmask return. */
        abnegative = !abpositive
                     && (SvIOK(sva) && !SvIsUV(sva))
                     && (SvIOK(svb) && !SvIsUV(svb));
        if (abpositive || abnegative) {
          UV a = my_svuv(sva);
          UV b = my_svuv(svb);
          int k = (abpositive) ? kronecker_uu(a,b) : kronecker_ss(a,b);
          RETURN_NPARITY(k);
        }
      } else if (ix == 1) {
        UV n = (astatus == -1) ? (UV)(-(my_sviv(sva))) : my_svuv(sva);
        UV k = (bstatus == -1) ? (UV)(-(my_sviv(svb))) : my_svuv(svb);
        /* valuation of 0-2 is very common, so return a constant if possible */
        RETURN_NPARITY( valuation(n, k) );
      } else if (ix == 2) {
        UV a, n, ret = 0;
        n = (bstatus != -1) ? my_svuv(svb) : (UV)(-(my_sviv(svb)));
        if (n > 0) {
          a = (astatus == 1) ? my_svuv(sva) : negmod(my_sviv(sva), n);
          if (a > 0) {
            if (n == 1) XSRETURN_UV(0);
            ret = modinverse(a, n);
          }
        }
        if (ret == 0) XSRETURN_UNDEF;
        XSRETURN_UV(ret);
      } else {
        UV a, n, s;
        n = (bstatus != -1) ? my_svuv(svb) : (UV)(-(my_sviv(svb)));
        a = (n == 0) ? 0 : (astatus != -1) ? my_svuv(sva) % n : negmod(my_sviv(sva), n);
        if (is_prob_prime(n)) {
          if (!sqrtmod(&s, a, n)) XSRETURN_UNDEF;
        } else {
          if (!sqrtmod_composite(&s, a, n)) XSRETURN_UNDEF;
        }
        XSRETURN_UV(s);
      }
    }
    switch (ix) {
      case 0:  _vcallsub_with_gmp(0.17,"kronecker");  break;
      case 1:  _vcallsub_with_gmp(0.20,"valuation"); break;
      case 2:  _vcallsub_with_gmp(0.20,"invmod"); break;
      case 3:
      default: _vcallsub_with_gmp(0.36,"sqrtmod"); break;
    }
    return; /* skip implicit PUTBACK */

void
gcdext(IN SV* sva, IN SV* svb)
  PREINIT:
    int astatus, bstatus;
  PPCODE:
    astatus = _validate_int(aTHX_ sva, 2);
    bstatus = _validate_int(aTHX_ svb, 2);
    /* TODO: These should be built into validate_int */
    if ( (astatus == 1 && SvIsUV(sva)) || (astatus == -1 && !SvIOK(sva)) )
      astatus = 0;  /* too large */
    if ( (bstatus == 1 && SvIsUV(svb)) || (bstatus == -1 && !SvIOK(svb)) )
      bstatus = 0;  /* too large */
    if (astatus != 0 && bstatus != 0) {
      IV u, v, d;
      IV a = my_sviv(sva);
      IV b = my_sviv(svb);
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
    if (type != 1 && type != 2)
      croak("stirling type must be 1 or 2");
    if (n == m)
      XSRETURN_UV(1);
    else if (n == 0 || m == 0 || m > n)
      XSRETURN_UV(0);
    else if (m == 1 && type == 2)
      XSRETURN_UV(1);
    else if (m == 1 && type == 1) {
      UV f = factorial(n-1);
      if (f != 0 && f <= (UV)IV_MAX) XSRETURN_IV( f * ((n&1) ? 1 : -1) );
    } else if (type == 2) {
      IV s = stirling2(n, m);
      if (s != 0) XSRETURN_IV(s);
    } else if (type == 1) {
      IV s = stirling1(n, m);
      if (s != 0) XSRETURN_IV(s);
    }
    _vcallsub_with_pp("stirling");
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
      case 0: ret = (NV) _XS_ExponentialIntegral(nv); break;
      case 1: ret = (NV) _XS_LogarithmicIntegral(nv); break;
      case 2: ret = (NV) ld_riemann_zeta(nv); break;
      case 3: ret = (NV) _XS_RiemannR(nv); break;
      case 4:
      default:ret = (NV) lambertw(nv); break;
    }
    RETVAL = ret;
  OUTPUT:
    RETVAL

void
euler_phi(IN SV* svlo, ...)
  ALIAS:
    moebius = 1
  PREINIT:
    int lostatus, histatus;
  PPCODE:
    lostatus = _validate_int(aTHX_ svlo, 2);
    histatus = (items == 1 || _validate_int(aTHX_ ST(1), 0));
    if (items == 1 && lostatus != 0) {
      /* input is a single value and in UV/IV range */
      if (ix == 0) {
        UV n = (lostatus == -1) ? 0 : my_svuv(svlo);
        XSRETURN_UV(totient(n));
      } else {
        UV n = (lostatus == -1) ? (UV)(-(my_sviv(svlo))) : my_svuv(svlo);
        RETURN_NPARITY(moebius(n));
      }
    } else if (items == 2 && lostatus == 1 && histatus == 1) {
      /* input is a range and both lo and hi are non-negative */
      UV lo = my_svuv(svlo);
      UV hi = my_svuv(ST(1));
      if (lo <= hi) {
        UV i;
        EXTEND(SP, (IV)(hi-lo+1));
        if (ix == 0) {
          UV  arraylo = (lo < 100)  ?  0  :  lo;
          UV* totients = _totient_range(arraylo, hi);
          for (i = lo; i <= hi; i++)
            PUSHs(sv_2mortal(newSVuv(totients[i-arraylo])));
          Safefree(totients);
        } else {
          signed char* mu = _moebius_range(lo, hi);
          dMY_CXT;
          for (i = lo; i <= hi; i++)
            PUSH_NPARITY(mu[i-lo]);
          Safefree(mu);
        }
      }
    } else {
      /* Whatever we didn't handle above */
      U32 gimme_v = GIMME_V;
      switch (ix) {
        case 0:  _vcallsubn(aTHX_ gimme_v, VCALL_PP, "euler_phi", items, 0.22);break;
        case 1:
        default: _vcallsubn(aTHX_ gimme_v, VCALL_GMP|VCALL_PP, "moebius", items, 22);  break;
      }
      return;
    }

void
carmichael_lambda(IN SV* svn)
  ALIAS:
    mertens = 1
    liouville = 2
    chebyshev_theta = 3
    chebyshev_psi = 4
    factorial = 5
    sqrtint = 6
    exp_mangoldt = 7
    znprimroot = 8
    hammingweight = 9
    hclassno = 10
    ramanujan_tau = 11
  PREINIT:
    int status;
  PPCODE:
    status = _validate_int(aTHX_ svn, (ix >= 7) ? 1 : 0);
    if (status != 0) {
      UV r, n = my_svuv(svn);
      switch (ix) {
        case 0:  XSRETURN_UV(carmichael_lambda(n)); break;
        case 1:  XSRETURN_IV(mertens(n)); break;
        case 2:  { UV factors[MPU_MAX_FACTORS+1];
                   int nfactors = factor(my_svuv(svn), factors);
                   RETURN_NPARITY( (nfactors & 1) ? -1 : 1 ); }
                 break;
        case 3:  XSRETURN_NV(chebyshev_function(n, 0)); break;
        case 4:  XSRETURN_NV(chebyshev_function(n, 1)); break;
        case 5:  r = factorial(n);
                 if (r != 0) XSRETURN_UV(r);
                 status = 0; break;
        case 6:  XSRETURN_UV(isqrt(n)); break;
        case 7:  XSRETURN_UV( (status == -1) ? 1 : exp_mangoldt(n) ); break;
        case 8:  if (status == -1) n = -(IV)n;
                 r = znprimroot(n);
                 if (r == 0 && n != 1)  XSRETURN_UNDEF;  /* No root */
                 XSRETURN_UV(r);  break;
        case 9:  if (status == -1) n = -(IV)n;
                 XSRETURN_UV(mpu_popcount(n));  break;
        case 10: XSRETURN_IV( (status == -1) ? 0 : hclassno(n) ); break;
        case 11:
        default: { IV tau = (status == 1) ? ramanujan_tau(n) : 0;
                   if (tau != 0 || status == -1 || n == 0)
                     XSRETURN_IV(tau);
                 } /* Fall through if n > 0 and we got 0 back */
                 break;
      }
    }
    switch (ix) {
      case 0:  _vcallsub_with_gmp(0.22,"carmichael_lambda"); break;
      case 1:  _vcallsub_with_pp("mertens"); break;
      case 2:  _vcallsub_with_gmp(0.22,"liouville"); break;
      case 3:  _vcallsub_with_pp("chebyshev_theta"); break;
      case 4:  _vcallsub_with_pp("chebyshev_psi"); break;
      case 5:  _vcallsub_with_pp("factorial"); break;
      case 6:  _vcallsub_with_pp("sqrtint"); break;
      case 7:  _vcallsub_with_gmp(0.19,"exp_mangoldt"); break;
      case 8:  _vcallsub_with_gmp(0.22,"znprimroot"); break;
      case 9:  { char* ptr;  STRLEN len;  ptr = SvPV(svn, len);
                 XSRETURN_UV(mpu_popcount_string(ptr, len)); } break;
      case 10:  _vcallsub_with_pp("hclassno"); break;
      case 11:
      default:_vcallsub_with_pp("ramanujan_tau"); break;
    }
    return; /* skip implicit PUTBACK */

void
sumdigits(SV* svn, UV ibase = 255)
  PREINIT:
    UV base, sum;
    STRLEN i, len;
    const char* s;
  PPCODE:
    base = (ibase == 255) ? 10 : ibase;
    if (base < 2 || base > 36) croak("sumdigits: invalid base %lu", base);
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
      status = _validate_int(aTHX_ svn, 1);
      n = (status == 0) ? 0 : status * my_svuv(svn);
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
      } else {            /* array ref of digits */
        UV* r = 0;
        int len = arrayref_to_int_array(aTHX_ &r, (AV*) SvRV(svn), base);
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
      case 0:  _vcallsubn(aTHX_ GIMME_V, VCALL_GMP|VCALL_PP, "todigits", items, 0); break;
      case 1:  _vcallsub_with_gmp(0.00,"todigitstring"); break;
      case 2:
      default: _vcallsub_with_gmp(0.00,"fromdigits"); break;
    }
    return;

bool
_validate_num(SV* svn, ...)
  PREINIT:
    SV* sv1;
    SV* sv2;
  CODE:
    /* Internal function.  Emulate the PP version of this:
     *   $is_valid = _validate_num( $n [, $min [, $max] ] )
     * Return 0 if we're befuddled by the input.
     * Otherwise croak if n isn't >= 0 and integer, n < min, or n > max.
     * Small bigints will be converted to scalars.
     */
    RETVAL = FALSE;
    if (_validate_int(aTHX_ svn, 0)) {
      if (SvROK(svn)) {  /* Convert small Math::BigInt object into scalar */
        UV n = my_svuv(svn);
#if PERL_REVISION <= 5 && PERL_VERSION < 8 && BITS_PER_WORD == 64
        sv_setpviv(svn, n);
#else
        sv_setuv(svn, n);
#endif
      }
      if (items > 1 && ((sv1 = ST(1)), SvOK(sv1))) {
        UV n = my_svuv(svn);
        UV min = my_svuv(sv1);
        if (n < min)
          croak("Parameter '%"UVuf"' must be >= %"UVuf, n, min);
        if (items > 2 && ((sv2 = ST(2)), SvOK(sv2))) {
          UV max = my_svuv(sv2);
          if (n > max)
            croak("Parameter '%"UVuf"' must be <= %"UVuf, n, max);
          MPUassert( items <= 3, "_validate_num takes at most 3 parameters");
        }
      }
      RETVAL = TRUE;
    }
  OUTPUT:
    RETVAL

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
  PPCODE:
    cv = sv_2cv(block, &stash, &gv, 0);
    if (cv == Nullcv)
      croak("Not a subroutine reference");

    if (!_validate_int(aTHX_ svbeg, 0) || (items >= 3 && !_validate_int(aTHX_ svend,0))) {
      _vcallsubn(aTHX_ G_VOID|G_DISCARD, VCALL_ROOT, "_generic_forprimes", items, 0);
      return;
    }

    if (items < 3) {
      beg = 2;
      end = my_svuv(svbeg);
    } else {
      beg = my_svuv(svbeg);
      end = my_svuv(svend);
    }

    SAVESPTR(GvSV(PL_defgv));
    svarg = newSVuv(beg);
    GvSV(PL_defgv) = svarg;
    /* Handle early part */
    while (beg < 6) {
      beg = (beg <= 2) ? 2 : (beg <= 3) ? 3 : 5;
      if (beg <= end) {
        sv_setuv(svarg, beg);
        PUSHMARK(SP);
        call_sv((SV*)cv, G_VOID|G_DISCARD);
      }
      beg += 1 + (beg > 2);
    }
#if USE_MULTICALL
    if (!CvISXSUB(cv) && beg <= end) {
      dMULTICALL;
      I32 gimme = G_VOID;
      PUSH_MULTICALL(cv);
      if (
#if BITS_PER_WORD == 64
          (beg >= UVCONST(     100000000000000) && end-beg <    100000) ||
          (beg >= UVCONST(      10000000000000) && end-beg <     40000) ||
          (beg >= UVCONST(       1000000000000) && end-beg <     17000) ||
#endif
          ((end-beg) < 500) ) {     /* MULTICALL next prime */
        for (beg = next_prime(beg-1); beg <= end && beg != 0; beg = next_prime(beg)) {
          sv_setuv(svarg, beg);
          MULTICALL;
        }
      } else {                      /* MULTICALL segment sieve */
        void* ctx = start_segment_primes(beg, end, &segment);
        while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
          int crossuv = (seg_high > IV_MAX) && !SvIsUV(svarg);
          START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
            /* sv_setuv(svarg, p); */
            if      (SvTYPE(svarg) != SVt_IV) { sv_setuv(svarg, p);            }
            else if (crossuv && p > IV_MAX)   { sv_setuv(svarg, p); crossuv=0; }
            else                              { SvUV_set(svarg, p);            }
            MULTICALL;
          END_DO_FOR_EACH_SIEVE_PRIME
        }
        end_segment_primes(ctx);
      }
      FIX_MULTICALL_REFCOUNT;
      POP_MULTICALL;
    }
    else
#endif
    if (beg <= end) {               /* NO-MULTICALL segment sieve */
      void* ctx = start_segment_primes(beg, end, &segment);
      while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
        START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
          sv_setuv(svarg, p);
          PUSHMARK(SP);
          call_sv((SV*)cv, G_VOID|G_DISCARD);
        END_DO_FOR_EACH_SIEVE_PRIME
      }
      end_segment_primes(ctx);
    }
    SvREFCNT_dec(svarg);

void
forcomposites (SV* block, IN SV* svbeg, IN SV* svend = 0)
  ALIAS:
    foroddcomposites = 1
  PROTOTYPE: &$;$
  PREINIT:
    UV beg, end;
    GV *gv;
    HV *stash;
    SV* svarg;  /* We use svarg to prevent clobbering $_ outside the block */
    CV *cv;
  PPCODE:
    cv = sv_2cv(block, &stash, &gv, 0);
    if (cv == Nullcv)
      croak("Not a subroutine reference");

    if (!_validate_int(aTHX_ svbeg, 0) || (items >= 3 && !_validate_int(aTHX_ svend,0))) {
      _vcallsubn(aTHX_ G_VOID|G_DISCARD, VCALL_ROOT, (ix == 0) ? "_generic_forcomposites" : "_generic_foroddcomposites", items, 0);
      return;
    }

    if (items < 3) {
      beg = ix ? 8 : 4;
      end = my_svuv(svbeg);
    } else {
      beg = my_svuv(svbeg);
      end = my_svuv(svend);
    }

    SAVESPTR(GvSV(PL_defgv));
    svarg = newSVuv(0);
    GvSV(PL_defgv) = svarg;
#if USE_MULTICALL
    if (!CvISXSUB(cv) && end >= beg) {
      unsigned char* segment;
      UV seg_base, seg_low, seg_high, c, cbeg, cend, prevprime, nextprime;
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
          if (beg == nextprime)     nextprime = next_prime(beg);
          else if (!ix || beg & 1)  { sv_setuv(svarg, beg); MULTICALL; }
        }
      } else {
        if (ix) {
          if (beg < 8)  beg = 8;
        } else if (beg <= 4) { /* sieve starts at 7, so handle this here */
          sv_setuv(svarg, 4);  MULTICALL;
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
              cbeg = beg - (ix == 1 && (beg % 2));
            prevprime = p;
            cend = prevprime-1;  if (cend > end) cend = end;
            /* If ix=1, skip evens by starting 1 farther and skipping by 2 */
            for (c = cbeg + ix; c <= cend; c=c+1+ix) {
              if      (SvTYPE(svarg) != SVt_IV) { sv_setuv(svarg,c); }
              else if (crossuv && c > IV_MAX)   { sv_setuv(svarg,c); crossuv=0;}
              else                              { SvUV_set(svarg,c); }
              MULTICALL;
            }
          END_DO_FOR_EACH_SIEVE_PRIME
        }
        end_segment_primes(ctx);
        if (end > nextprime)   /* Complete the case where end > max_prime */
          while (nextprime++ < end)
            if (!ix || nextprime & 1)
              { sv_setuv(svarg, nextprime);  MULTICALL; }
      }
      FIX_MULTICALL_REFCOUNT;
      POP_MULTICALL;
    }
    else
#endif
    if (beg <= end) {
      beg = (beg <= 4) ? 3 : beg-1;
      while (beg++ < end) {
        if ((!ix || beg&1) && !is_prob_prime(beg)) {
          sv_setuv(svarg, beg);
          PUSHMARK(SP);
          call_sv((SV*)cv, G_VOID|G_DISCARD);
        }
      }
    }
    SvREFCNT_dec(svarg);

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
  PPCODE:
    cv = sv_2cv(block, &stash, &gv, 0);
    if (cv == Nullcv)
      croak("Not a subroutine reference");

    if (!_validate_int(aTHX_ svn, 0)) {
      _vcallsubn(aTHX_ G_VOID|G_DISCARD, VCALL_ROOT, "_generic_fordivisors", 2, 0);
      return;
    }

    n = my_svuv(svn);
    divs = _divisor_list(n, &ndivisors);

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
        MULTICALL;
      }
      FIX_MULTICALL_REFCOUNT;
      POP_MULTICALL;
    }
    else
#endif
    {
      for (i = 0; i < ndivisors; i++) {
        sv_setuv(svarg, divs[i]);
        PUSHMARK(SP);
        call_sv((SV*)cv, G_VOID|G_DISCARD);
      }
    }
    SvREFCNT_dec(svarg);
    Safefree(divs);

void
forpart (SV* block, IN SV* svn, IN SV* svh = 0)
  ALIAS:
    forcomp = 1
  PROTOTYPE: &$;$
  PREINIT:
    UV i, n, amin, amax, nmin, nmax;
    GV *gv;
    HV *stash;
    CV *cv;
    SV** svals;
  PPCODE:
    cv = sv_2cv(block, &stash, &gv, 0);
    if (cv == Nullcv)
      croak("Not a subroutine reference");
    if (!_validate_int(aTHX_ svn, 0)) {
      _vcallsub_with_pp("forpart");
      return;
    }
    n = my_svuv(svn);
    if (n > (UV_MAX-2)) croak("forpart argument overflow");

    New(0, svals, n+1, SV*);
    for (i = 0; i <= n; i++) {
      svals[i] = newSVuv(i);
      SvREADONLY_on(svals[i]);
    }

    amin = 1;  amax = n;  nmin = 1;  nmax = n;
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

      if (amin < 1) amin = 1;
      if (amax > n) amax = n;
      if (nmin < 1) nmin = 1;
      if (nmax > n) nmax = n;
    }

    if (n==0 && nmin <= 1) {
      { dSP; ENTER; PUSHMARK(SP);
        /* Nothing */
        PUTBACK; call_sv((SV*)cv, G_VOID|G_DISCARD); LEAVE;
      }
    }
    if (n >= nmin && nmin <= nmax && amin <= amax && nmax > 0 && amax > 0)
    { /* RuleAsc algorithm from Kelleher and O'Sullivan 2009/2014) */
      UV *a, k, x, y, r;
      New(0, a, n+1, UV);
      k = 1;
      a[0] = amin-1;
      a[1] = n-amin+1;
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

        { dSP; ENTER; PUSHMARK(SP);
          EXTEND(SP, (IV)k); for (i = 0; i <= k; i++) { PUSHs(svals[a[i]]); }
          PUTBACK; call_sv((SV*)cv, G_VOID|G_DISCARD); LEAVE;
        }
      }
      Safefree(a);
    }
    for (i = 0; i <= n; i++)
      SvREFCNT_dec(svals[i]);
    Safefree(svals);

void
forcomb (SV* block, IN SV* svn, IN SV* svk = 0)
  ALIAS:
    forperm = 1
  PROTOTYPE: &$;$
  PREINIT:
    UV i, n, k, j, m;
    GV *gv;
    HV *stash;
    CV *cv;
    SV** svals;
    UV*  cm;
  PPCODE:
    cv = sv_2cv(block, &stash, &gv, 0);
    if (cv == Nullcv)
      croak("Not a subroutine reference");
    if (ix == 1 && svk != 0)
      croak("Too many arguments for forperm");

    if (!_validate_int(aTHX_ svn, 0) || (svk != 0 && !_validate_int(aTHX_ svk, 0))) {
      _vcallsub_with_pp( (ix == 0) ? "forcomb" : "forperm" );
      return;
    }

    n = my_svuv(svn);
    k = (svk == 0) ? n : my_svuv(svk);
    if (k > n)
      return;

    New(0, cm, k+1, UV);
    cm[0] = UV_MAX;
    for (i = 0; i < k; i++)
      cm[i] = k-i;

    New(0, svals, n, SV*);
    for (i = 0; i < n; i++) {
      svals[i] = newSVuv(i);
      SvREADONLY_on(svals[i]);
    }

    while (1) {
      { dSP; ENTER; PUSHMARK(SP);                /* Send the values */
        EXTEND(SP, ((IV)k));
        for (i = 0; i < k; i++) { PUSHs(svals[ cm[k-i-1]-1 ]); }
        PUTBACK; call_sv((SV*)cv, G_VOID|G_DISCARD); LEAVE;
      }
      if (ix == 0) {
        if (cm[0]++ < n)  continue;                /* Increment last value */
        for (i = 1; i < k && cm[i] >= n-i; i++) ;  /* Find next index to incr */
        if (i >= k)  break;                        /* Done! */
        cm[i]++;                                   /* Increment this one */
        while (i-- > 0)  cm[i] = cm[i+1] + 1;      /* Set the rest */
      } else {
        for (j = 1; j < k && cm[j] > cm[j-1]; j++) ;    /* Find last decrease */
        if (j >= k) break;                              /* Done! */
        for (m = 0; cm[j] > cm[m]; m++)                 /* Find next greater */
          ;
        { UV t = cm[j];  cm[j] = cm[m];  cm[m] = t; }   /* Swap */
        for (i = j-1, m = 0;  m < i;  i--, m++)         /* Reverse the end */
          { UV t = cm[i];  cm[i] = cm[m];  cm[m] = t; }
      }
    }
    Safefree(cm);
    for (i = 0; i < n; i++)
      SvREFCNT_dec(svals[i]);
    Safefree(svals);

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
        MULTICALL;
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
        MULTICALL;
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

    if (index != items)           /* We exited the loop early */
      ret_true = !ret_true;

    if (ret_true)  XSRETURN_YES;
    else           XSRETURN_NO;
}
