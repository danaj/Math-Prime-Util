
#define PERL_NO_GET_CONTEXT 1 /* Define at top for more efficiency. */

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "multicall.h"  /* only works in 5.6 and newer */

#define NEED_sv_2pv_flags
#include "ppport.h"

#include "ptypes.h"
#include "cache.h"
#include "sieve.h"
#define FUNC_gcd_ui 1
#include "util.h"
#include "primality.h"
#include "factor.h"
#include "lehmer.h"
#include "lmo.h"
#include "aks.h"
#include "constants.h"

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

#if PERL_REVISION <= 5 && PERL_VERSION <= 6 && BITS_PER_WORD == 64
 /* Workaround perl 5.6 UVs and bigints */
 #define my_svuv(sv)  PSTRTOULL(SvPV_nolen(sv), NULL, 10)
 #define my_sviv(sv)  PSTRTOLL(SvPV_nolen(sv), NULL, 10)
#elif PERL_REVISION <= 5 && PERL_VERSION < 14 && BITS_PER_WORD == 64
 /* Workaround RT 49569 in Math::BigInt::FastCalc (pre 5.14.0) */
 #define my_svuv(sv) ( (!SvROK(sv)) ? SvUV(sv) : PSTRTOULL(SvPV_nolen(sv),NULL,10) )
 #define my_sviv(sv) ( (!SvROK(sv)) ? SvIV(sv) : PSTRTOLL(SvPV_nolen(sv),NULL,10) )
#else
 #define my_svuv(sv) SvUV(sv)
 #define my_sviv(sv) SvIV(sv)
#endif

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
  SV* const_int[4];   /* -1, 0, 1, 2 */
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
  if ((SvFLAGS(n) & (SVf_IOK |
#if PERL_REVISION >=5 && PERL_VERSION >= 9 && PERL_SUBVERSION >= 4
                     SVf_ROK |
#else
                     SVf_AMAGIC |
#endif
                     SVs_GMG )) == SVf_IOK) { /* If defined as number, use it */
    if (SvIsUV(n) || SvIVX(n) >= 0)  return 1; /* The normal case */
    if (negok)  return -1;
    else croak("Parameter '%" SVf "' must be a positive integer", n);
  }
  if (SvROK(n)) {
    if (sv_isa(n, "Math::BigInt") || sv_isa(n, "Math::BigFloat") ||
        sv_isa(n, "Math::Pari") || sv_isa(n, "Math::GMP") ||
        sv_isa(n, "Math::GMPz") )
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
static int _vcallsubn(pTHX_ I32 flags, I32 stashflags, const char* name, int nargs)
{
    GV* gv = NULL;
    dMY_CXT;
    Size_t namelen = strlen(name);
    /* If given a GMP function, and GMP enabled, and function exists, use it. */
    int use_gmp = stashflags & VCALL_GMP && _XS_get_callgmp();
    assert(!(stashflags & ~(VCALL_PP|VCALL_GMP)));
    if (use_gmp) {
      GV ** gvp = (GV**)hv_fetch(MY_CXT.MPUGMP,name,namelen,0);
      if (gvp) gv = *gvp;
    }
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
#define _vcallsub(func) (void)_vcallsubn(aTHX_ G_SCALAR, VCALL_ROOT, func, items)
#define _vcallsub_with_gmp(func) (void)_vcallsubn(aTHX_ G_SCALAR, VCALL_GMP|VCALL_PP, func, items)
#define _vcallsub_with_pp(func) (void)_vcallsubn(aTHX_ G_SCALAR, VCALL_PP, func, items)

/* In my testing, this constant return works fine with threads, but to be
 * correct (see perlxs) one has to make a context, store separate copies in
 * each one, then retrieve them from a struct using a hash index.  This
 * defeats the purpose if only done once. */
#define RETURN_NPARITY(ret) \
  do { int r_ = ret; \
       dMY_CXT; \
       if (r_ >= -1 && r_ <= 2) { ST(0) = MY_CXT.const_int[r_+1]; XSRETURN(1); } \
       else                     { XSRETURN_IV(r_);                      } \
  } while (0)
#define PUSH_NPARITY(ret) \
  do { int r_ = ret; \
       if (r_ >= -1 && r_ <= 2) { PUSHs( MY_CXT.const_int[r_+1] );       } \
       else                     { PUSHs(sv_2mortal(newSViv(r_))); } \
  } while (0)

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
      for (i = 0; i <= 3; i++) {
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
    for (i = 0; i <= 3; i++) {
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
  for (i = 0; i <= 3; i++) {
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
        if (ix == 1 || (hi / (hi-lo+1)) > 100) {
          count = _XS_prime_count(lo, hi);
        } else {
          count = _XS_LMO_pi(hi);
          if (lo > 2)
            count -= _XS_LMO_pi(lo-1);
        }
      }
      XSRETURN_UV(count);
    }
    _vcallsubn(aTHX_ GIMME_V, VCALL_ROOT, "_generic_prime_count", items);
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
    if ((low <= 2) && (high >= 2)) { av_push(av, newSVuv( 2 )); }
    if ((low <= 3) && (high >= 3)) { av_push(av, newSVuv( 3 )); }
    if ((low <= 5) && (high >= 5)) { av_push(av, newSVuv( 5 )); }
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
        START_DO_FOR_EACH_SIEVE_PRIME( sieve, low, high ) {
           av_push(av,newSVuv(p));
        } END_DO_FOR_EACH_SIEVE_PRIME
        Safefree(sieve);
      } else if (ix == 3) {                   /* Segment */
        unsigned char* segment;
        UV seg_base, seg_low, seg_high;
        void* ctx = start_segment_primes(low, high, &segment);
        while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
          START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_low - seg_base, seg_high - seg_base )
            av_push(av,newSVuv( seg_base + p ));
          END_DO_FOR_EACH_SIEVE_PRIME
        }
        end_segment_primes(ctx);
      }
    }
    return; /* skip implicit PUTBACK */

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
  PREINIT:
    UV arg1, arg2;
    static const UV default_arg1[] =
       {0,     64000000, 8000000, 4000000, 4000000, 200, 4000000, 1000000};
     /* Trial, Fermat,   Holf,    SQUFOF,  PRHO,    P+1, Brent,    P-1 */
  PPCODE:
    if (n == 0)  XSRETURN_UV(0);
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
  PREINIT:
    int c, status = 1;
  PPCODE:
    if (items < 2)
      croak("No bases given to miller_rabin");
    /* Check all arguments */
    for (c = 0; c < items && status == 1; c++)
      if (_validate_int(aTHX_ ST(c), 0) != 1)
        status = 0;
    if (status == 1) {
      UV n = my_svuv(svn);
      int b, ret = 1;
      if      (n < 4)        { ret = (n >= 2); } /* 0,1 composite; 2,3 prime */
      else if ((n % 2) == 0) { ret = 0; }        /* evens composite */
      else {
        UV bases[32];
        for (c = 1; c < items && ret == 1; ) {
          for (b = 0; b < 32 && c < items; c++)
            bases[b++] = my_svuv(ST(c));
          ret = _XS_miller_rabin(n, bases, b);
        }
      }
      RETURN_NPARITY(ret);
    }
    _vcallsub_with_gmp("is_strong_pseudoprime");
    return; /* skip implicit PUTBACK */

void
gcd(...)
  PROTOTYPE: @
  ALIAS:
    lcm = 1
  PREINIT:
    int i, status = 1;
    UV ret, nullv, n;
  PPCODE:
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
    if (status != 0)
      XSRETURN_UV(ret);
    switch (ix) {
      case 0: _vcallsub_with_gmp("gcd");  break;
      case 1:
      default:_vcallsub_with_gmp("lcm");  break;
    }
    return; /* skip implicit PUTBACK */

void
_XS_lucas_sequence(IN UV n, IN IV P, IN IV Q, IN UV k)
  PREINIT:
    UV U, V, Qk;
  PPCODE:
    lucas_seq(&U, &V, &Qk,  n, P, Q, k);
    PUSHs(sv_2mortal(newSVuv( U )));    /* 4 args in, 3 out, no EXTEND needed */
    PUSHs(sv_2mortal(newSVuv( V )));
    PUSHs(sv_2mortal(newSVuv( Qk )));

void
is_prime(IN SV* svn, ...)
  ALIAS:
    is_prob_prime = 1
    is_bpsw_prime = 2
    is_lucas_pseudoprime = 3
    is_strong_lucas_pseudoprime = 4
    is_extra_strong_lucas_pseudoprime = 5
    is_frobenius_underwood_pseudoprime = 6
    is_aks_prime = 7
    is_pseudoprime = 8
    is_almost_extra_strong_lucas_pseudoprime = 9
  PREINIT:
    int status;
  PPCODE:
    status = _validate_int(aTHX_ svn, 1);
    if (status != 0) {
      int ret = 0;
      if (status == 1) {
        UV n = my_svuv(svn);
        UV a = (items == 1) ? 0 : my_svuv(ST(1));
        switch (ix) {
          case 0:
          case 1:  ret = _XS_is_prime(n);  break;
          case 2:  ret = _XS_BPSW(n);      break;
          case 3:  ret = _XS_is_lucas_pseudoprime(n, 0); break;
          case 4:  ret = _XS_is_lucas_pseudoprime(n, 1); break;
          case 5:  ret = _XS_is_lucas_pseudoprime(n, 2); break;
          case 6:  ret = _XS_is_frobenius_underwood_pseudoprime(n); break;
          case 7:  ret = _XS_is_aks_prime(n); break;
          case 8:  ret = _XS_is_pseudoprime(n, (items == 1) ? 2 : a); break;
          case 9:
          default: ret = _XS_is_almost_extra_strong_lucas_pseudoprime
                         (n, (items == 1) ? 1 : a); break;
        }
      }
      RETURN_NPARITY(ret);
    }
    switch (ix) {
      case 0: _vcallsub_with_gmp("is_prime");       break;
      case 1: _vcallsub_with_gmp("is_prob_prime");  break;
      case 2: _vcallsub_with_gmp("is_bpsw_prime");  break;
      case 3: _vcallsub_with_gmp("is_lucas_pseudoprime"); break;
      case 4: _vcallsub_with_gmp("is_strong_lucas_pseudoprime"); break;
      case 5: _vcallsub_with_gmp("is_extra_strong_lucas_pseudoprime"); break;
      case 6: _vcallsub_with_gmp("is_frobenius_underwood_pseudoprime"); break;
      case 7: _vcallsub_with_gmp("is_aks_prime"); break;
      case 8: _vcallsub_with_gmp("is_pseudoprime"); break;
      case 9:
      default:_vcallsub_with_gmp("is_almost_extra_strong_lucas_pseudoprime"); break;
    }
    return; /* skip implicit PUTBACK */

void
next_prime(IN SV* svn)
  ALIAS:
    prev_prime = 1
    nth_prime = 2
  PPCODE:
    if (_validate_int(aTHX_ svn, 0)) {
      UV n = my_svuv(svn);
      if ( ((ix == 0) && (n >= MPU_MAX_PRIME))    ||
           ((ix == 2) && (n >= MPU_MAX_PRIME_IDX)) ) {
        /* Out of range.  Fall through to Perl. */
      } else {
        UV ret;
        switch (ix) {
          case 0: ret = next_prime(n);  break;
          case 1: ret = (n < 3) ? 0 : prev_prime(n);  break;
          case 2:
          default:ret = _XS_nth_prime(n);  break;
        }
        XSRETURN_UV(ret);
      }
    }
    switch (ix) {
      case 0:  _vcallsub("_generic_next_prime");     break;
      case 1:  _vcallsub("_generic_prev_prime");     break;
      default: _vcallsub_with_pp("nth_prime");           break;
    }
    return; /* skip implicit PUTBACK */

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
                     EXTEND(SP, ndivisors);
                     for (i = 0; (UV)i < ndivisors; i++)
                       PUSHs(sv_2mortal(newSVuv( divs[i] )));
                     Safefree(divs);
                   }
                   break;
        }
      }
    } else {
      switch (ix) {
        case 0:  _vcallsubn(aTHX_ gimme_v, VCALL_ROOT, "_generic_factor", 1);     break;
        case 1:  _vcallsubn(aTHX_ gimme_v, VCALL_ROOT, "_generic_factor_exp", 1); break;
        default: _vcallsubn(aTHX_ gimme_v, VCALL_ROOT, "_generic_divisors", 1);   break;
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
    kstatus = (items == 1 || (SvIOK(svk) && SvIV(svk)))  ?  1  :  0;
    if (nstatus == 1 && kstatus == 1) {
      UV n = my_svuv(svn);
      UV k = (items > 1) ? my_svuv(svk) : 1;
      UV sigma = divisor_sum(n, k);
      if (sigma != 0)  XSRETURN_UV(sigma);   /* sigma 0 means overflow */
    }
    _vcallsub("_generic_divisor_sum");
    return; /* skip implicit PUTBACK */

void
znorder(IN SV* sva, IN SV* svn)
  ALIAS:
    legendre_phi = 1
  PREINIT:
    int astatus, nstatus;
  PPCODE:
    astatus = _validate_int(aTHX_ sva, 0);
    nstatus = _validate_int(aTHX_ svn, 0);
    if (astatus == 1 && nstatus == 1) {
      UV a = my_svuv(sva);
      UV n = my_svuv(svn);
      UV ret;
      switch (ix) {
        case 0:  ret = znorder(a, n);
                 if (ret == 0) XSRETURN_UNDEF;  /* not defined */
                 break;
        case 1:
        default: ret = legendre_phi(a, n); break;
      }
      XSRETURN_UV(ret);
    }
    switch (ix) {
      case 0:  _vcallsub_with_pp("znorder");  break;
      case 1:
      default: _vcallsub_with_pp("legendre_phi"); break;
    }
    return; /* skip implicit PUTBACK */

void
znlog(IN SV* sva, IN SV* svg, IN SV* svp)
  PREINIT:
    int astatus, gstatus, pstatus;
  PPCODE:
    astatus = _validate_int(aTHX_ sva, 0);
    gstatus = _validate_int(aTHX_ svg, 0);
    pstatus = _validate_int(aTHX_ svp, 0);
    if (astatus == 1 && gstatus == 1 && pstatus == 1) {
      UV ret = znlog(my_svuv(sva), my_svuv(svg), my_svuv(svp));
      if (ret == 0) XSRETURN_UNDEF;
      XSRETURN_UV(ret);
    }
    _vcallsub_with_pp("znlog");
    return; /* skip implicit PUTBACK */

void
kronecker(IN SV* sva, IN SV* svb)
  PREINIT:
    int astatus, bstatus, abpositive, abnegative;
  PPCODE:
    astatus = _validate_int(aTHX_ sva, 2);
    bstatus = _validate_int(aTHX_ svb, 2);
    /* Are both a and b positive? */
    abpositive = astatus == 1 && bstatus == 1;
    /* Will both fit in IVs?  We should use a bitmask return. */
    abnegative = !abpositive
                 && (astatus != 0 && SvIOK(sva) && !SvIsUV(sva))
                 && (bstatus != 0 && SvIOK(svb) && !SvIsUV(svb));
    if (abpositive || abnegative) {
      UV a = my_svuv(sva);
      UV b = my_svuv(svb);
      int k = (abpositive) ? kronecker_uu(a,b) : kronecker_ss(a,b);
      RETURN_NPARITY(k);
    }
    _vcallsub("_generic_kronecker");
    return; /* skip implicit PUTBACK */

double
_XS_ExponentialIntegral(IN SV* x)
  ALIAS:
    _XS_LogarithmicIntegral = 1
    _XS_RiemannZeta = 2
    _XS_RiemannR = 3
    _XS_chebyshev_theta = 4
    _XS_chebyshev_psi = 5
  PREINIT:
    double ret;
  CODE:
    if (ix < 4) {
      NV nv = SvNV(x);
      switch (ix) {
        case 0: ret = (NV) _XS_ExponentialIntegral(nv); break;
        case 1: ret = (NV) _XS_LogarithmicIntegral(nv); break;
        case 2: ret = (NV) ld_riemann_zeta(nv); break;
        case 3:
        default:ret = (NV) _XS_RiemannR(nv); break;
      }
    } else {
      UV uv = SvUV(x);
      ret = (NV) chebyshev_function(uv, ix-4);
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
        EXTEND(SP, hi-lo+1);
        if (ix == 0) {
          UV* totients = _totient_range(lo, hi);
          for (i = lo; i <= hi; i++)
            PUSHs(sv_2mortal(newSVuv(totients[i-lo])));
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
        case 0:  _vcallsubn(aTHX_ gimme_v, VCALL_ROOT,"_generic_euler_phi", items);break;
        case 1:
        default: _vcallsubn(aTHX_ gimme_v, VCALL_ROOT,"_generic_moebius", items);  break;
      }
      return;
    }

void
carmichael_lambda(IN SV* svn)
  ALIAS:
    mertens = 1
    exp_mangoldt = 2
    znprimroot = 3
  PREINIT:
    int status;
  PPCODE:
    status = _validate_int(aTHX_ svn, (ix > 1) ? 1 : 0);
    switch (ix) {
      case 0: if (status == 1) XSRETURN_UV(carmichael_lambda(my_svuv(svn)));
              _vcallsub("_generic_carmichael_lambda");
              break;
      case 1: if (status == 1) XSRETURN_IV(mertens(my_svuv(svn)));
              _vcallsub("_generic_mertens");
              break;
      case 2: if (status ==-1) XSRETURN_UV(1);
              if (status == 1) XSRETURN_UV(exp_mangoldt(my_svuv(svn)));
              _vcallsub("_generic_exp_mangoldt");
              break;
      case 3:
      default:if (status != 0) {
                UV r, n = my_svuv(svn);
                if (status == -1) n = -(IV)n;
                r = znprimroot(n);
                if (r == 0 && n != 1)
                  XSRETURN_UNDEF;   /* No root, return undef */
                else
                  XSRETURN_UV(r);
              }
              _vcallsub("_generic_znprimroot");
              break;
    }
    return; /* skip implicit PUTBACK */

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
      _vcallsubn(aTHX_ G_VOID|G_DISCARD, VCALL_ROOT, "_generic_forprimes", items);
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
    svarg = newSVuv(0);
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
          START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_low - seg_base, seg_high - seg_base ) {
            sv_setuv(svarg, seg_base + p);
            MULTICALL;
          } END_DO_FOR_EACH_SIEVE_PRIME
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
        START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_low - seg_base, seg_high - seg_base ) {
          sv_setuv(svarg, seg_base + p);
          PUSHMARK(SP);
          call_sv((SV*)cv, G_VOID|G_DISCARD);
        } END_DO_FOR_EACH_SIEVE_PRIME
      }
      end_segment_primes(ctx);
    }
    SvREFCNT_dec(svarg);

void
forcomposites (SV* block, IN SV* svbeg, IN SV* svend = 0)
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
      _vcallsubn(aTHX_ G_VOID|G_DISCARD, VCALL_ROOT, "_generic_forcomposites", items);
      return;
    }

    if (items < 3) {
      beg = 4;
      end = my_svuv(svbeg);
    } else {
      beg = my_svuv(svbeg);
      end = my_svuv(svend);
    }

    SAVESPTR(GvSV(PL_defgv));
    svarg = newSVuv(0);
    GvSV(PL_defgv) = svarg;
#if USE_MULTICALL
    if (!CvISXSUB(cv) && (end-beg) > 200) {
      unsigned char* segment;
      UV seg_base, seg_low, seg_high, c, cbeg, cend, prevprime, nextprime;
      void* ctx;
      dMULTICALL;
      I32 gimme = G_VOID;
      PUSH_MULTICALL(cv);
      if (beg <= 4) { /* sieve starts at 7, so handle this here */
        sv_setuv(svarg, 4);  MULTICALL;
        beg = 6;
      }
      /* Find the two primes that bound their interval. */
      /* If beg or end are >= max_prime, then this will die. */
      prevprime = prev_prime(beg);
      nextprime = next_prime(end);
      ctx = start_segment_primes(beg, nextprime, &segment);
      while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
        START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_low - seg_base, seg_high - seg_base ) {
          cbeg = prevprime+1;  if (cbeg < beg) cbeg = beg;
          prevprime = seg_base + p;
          cend = prevprime-1;  if (cend > end) cend = end;
          for (c = cbeg; c <= cend; c++) {
            sv_setuv(svarg, c);  MULTICALL;
          }
        } END_DO_FOR_EACH_SIEVE_PRIME
      }
      end_segment_primes(ctx);
      MPUassert( nextprime >= end, "composite sieve skipped end numbers" );
      FIX_MULTICALL_REFCOUNT;
      POP_MULTICALL;
    }
    else
#endif
    if (beg <= end) {
      beg = (beg <= 4) ? 3 : beg-1;
      while (beg++ < end) {
        if (!is_prob_prime(beg)) {
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
      _vcallsubn(aTHX_ G_VOID|G_DISCARD, VCALL_ROOT, "_generic_fordivisors", 2);
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
