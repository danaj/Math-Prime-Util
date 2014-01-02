
#define PERL_NO_GET_CONTEXT 1 /* Define at top for more efficiency. */

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "multicall.h"  /* only works in 5.6 and newer */
/* Perhaps we should use ppport.h */
#ifndef XSRETURN_UV   /* Fix 21086 from Sep 2003 */
  #define XST_mUV(i,v)  (ST(i) = sv_2mortal(newSVuv(v))  )
  #define XSRETURN_UV(v) STMT_START { XST_mUV(0,v);  XSRETURN(1); } STMT_END
#endif
#if PERL_REVISION <= 5 && (PERL_VERSION < 7 || (PERL_VERSION == 7 && PERL_SUBVERSION <= 2))
  #define SvPV_nomg SvPV
#endif
#include "ptypes.h"
#include "cache.h"
#include "sieve.h"
#include "util.h"
#include "primality.h"
#include "factor.h"
#include "lehmer.h"
#include "lmo.h"
#include "aks.h"

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

/* Workaround perl 5.6 UVs and bigints in later */
#if PERL_REVISION <= 5 && PERL_VERSION <= 6 && BITS_PER_WORD == 64
 #define my_svuv(sv)  PSTRTOULL(SvPV_nolen(sv), NULL, 10)
 #define my_sviv(sv)  PSTRTOLL(SvPV_nolen(sv), NULL, 10)
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
        sv_isa(n, "Math::GMP") || sv_isa(n, "Math::GMPz") )
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

/* Call a Perl sub to handle work for us. */
static int _vcallsubn(pTHX_ I32 flags, const char* name, int nargs)
{
    dSP;
    char fullname[80] = "Math::Prime::Util::";
    strncat(fullname, name, 60);
    PUSHMARK(SP-nargs);
    PUTBACK;
    return call_pv(fullname, flags);
}
#define _vcallsub(func) (void)_vcallsubn(aTHX_ G_SCALAR, func, 1)

#if BITS_PER_WORD == 64
static const UV _max_prime = UVCONST(18446744073709551557);
#else
static const UV _max_prime = UVCONST(4294967291);
#endif


MODULE = Math::Prime::Util	PACKAGE = Math::Prime::Util

PROTOTYPES: ENABLE

BOOT:
{
    SV * sv = newSViv(BITS_PER_WORD);
    HV * stash = gv_stashpv("Math::Prime::Util", TRUE);
    newCONSTSUB(stash, "_XS_prime_maxbits", sv);
}

void
prime_memfree()
  ALIAS:
    _prime_memfreeall = 1
    _XS_get_verbose = 2
    _XS_get_callgmp = 3
    _get_prime_cache_size = 4
  PPCODE:
    UV ret;
    switch (ix) {
      case 0:  prime_memfree(); goto return_nothing;
      case 1:  _prime_memfreeall(); goto return_nothing;
      case 2:  ret = _XS_get_verbose(); break;
      case 3:  ret = _XS_get_callgmp(); break;
      case 4:
      default:  ret = get_prime_cache(0,0); break;
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
_XS_prime_count(IN UV low, IN UV high = 0)
  PPCODE:
    if (high == 0) {   /* Without a Perl layer in front of this, we'll have */
      high = low;      /* the pathological case of a-0 turning into 0-a.    */
      low = 0;
    }
    if (GIMME_V == G_VOID) {
      prime_precalc(high);
    } else {
      PUSHs(sv_2mortal(newSVuv( _XS_prime_count(low, high) )));
    }

UV
_XS_nth_prime(IN UV n)
  ALIAS:
    _XS_next_prime = 1
    _XS_prev_prime = 2
    _XS_legendre_pi = 3
    _XS_meissel_pi = 4
    _XS_lehmer_pi = 5
    _XS_LMOS_pi = 6
    _XS_LMO_pi = 7
  PREINIT:
    UV ret;
  CODE:
    switch (ix) {
      case 0: ret = _XS_nth_prime(n); break;
      case 1: ret = _XS_next_prime(n); break;
      case 2: ret = _XS_prev_prime(n); break;
      case 3: ret = _XS_legendre_pi(n); break;
      case 4: ret = _XS_meissel_pi(n); break;
      case 5: ret = _XS_lehmer_pi(n); break;
      case 6: ret = _XS_LMOS_pi(n); break;
      default:ret = _XS_LMO_pi(n); break;
    }
    RETVAL = ret;
  OUTPUT:
    RETVAL

SV*
sieve_primes(IN UV low, IN UV high)
  ALIAS:
    trial_primes = 1
    erat_primes = 2
    segment_primes = 3
  PREINIT:
    AV* av = newAV();
  CODE:
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
        for (low = _XS_next_prime(low-1);
             low <= high && low != 0;
             low = _XS_next_prime(low) ) {
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
    RETVAL = newRV_noinc( (SV*) av );
  OUTPUT:
    RETVAL

void
trial_factor(IN UV n, ...)
  ALIAS:
    fermat_factor = 1
    holf_factor = 2
    squfof_factor = 3
    rsqufof_factor = 4
    pbrent_factor = 5
    prho_factor = 6
    pplus1_factor = 7
    pminus1_factor = 8
  PPCODE:
    if (n == 0)  XSRETURN_UV(0);
    while ( (n% 2) == 0 ) {  n /=  2;  XPUSHs(sv_2mortal(newSVuv( 2 ))); }
    while ( (n% 3) == 0 ) {  n /=  3;  XPUSHs(sv_2mortal(newSVuv( 3 ))); }
    while ( (n% 5) == 0 ) {  n /=  5;  XPUSHs(sv_2mortal(newSVuv( 5 ))); }
    if (n == 1) {  /* done */ }
    else if (_XS_is_prime(n)) { XPUSHs(sv_2mortal(newSVuv( n ))); }
    else {
      UV arg1, arg2, factors[MPU_MAX_FACTORS+1];
      int i, nfactors = 0;
      switch (ix) {
        case 0:  arg1 = (items < 2)  ?  0             : SvUV(ST(1));
                 nfactors = trial_factor  (n, factors, arg1);  break;
        case 1:  arg1 = (items < 2)  ?  64*1024*1024  : SvUV(ST(1));
                 nfactors = fermat_factor (n, factors, arg1);  break;
        case 2:  arg1 = (items < 2)  ?   8*1024*1024  : SvUV(ST(1));
                 nfactors = holf_factor   (n, factors, arg1);  break;
        case 3:  arg1 = (items < 2)  ?   4*1024*1024  : SvUV(ST(1));
                 nfactors = squfof_factor (n, factors, arg1);  break;
        case 4:  arg1 = (items < 2)  ?   4*1024*1024  : SvUV(ST(1));
                 nfactors = racing_squfof_factor(n, factors, arg1);  break;
        case 5:  arg1 = (items < 2)  ?   4*1024*1024  : SvUV(ST(1));
                 arg2 = (items < 3)  ?             1  : SvUV(ST(2));
                 nfactors = pbrent_factor (n, factors, arg1, arg2);  break;
        case 6:  arg1 = (items < 2)  ?   4*1024*1024  : SvUV(ST(1));
                 nfactors = prho_factor   (n, factors, arg1);  break;
        case 7:  arg1 = (items < 2)  ?           200  : SvUV(ST(1));
                 nfactors = pplus1_factor (n, factors, arg1);  break;
        case 8:  arg1 = (items < 2)  ?   1*1024*1024  : SvUV(ST(1));
                 arg2 = (items < 3)  ?             0  : SvUV(ST(2));
                 if (arg2 == 0) arg2 = 10*arg1;  /* default B2 */
                 nfactors = pminus1_factor(n, factors, arg1, arg2);  break;
        default: break;
      }
      EXTEND(SP, nfactors);
      for (i = 0; i < nfactors; i++)
        PUSHs(sv_2mortal(newSVuv( factors[i] )));
    }

void
is_strong_pseudoprime(IN SV* svn, ...)
  PREINIT:
    int status;
  PPCODE:
    if (items < 2)
      croak("No bases given to miller_rabin");
    status = _validate_int(aTHX_ svn, 0);
    if (status == 1) {
      UV n = my_svuv(svn);
      int b, c, ret = 1;
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
      XSRETURN_UV(ret);
    } else {
      const char* sub = _XS_get_callgmp() ? "GMP::is_strong_pseudoprime"
                                          : "_generic_is_strong_pseudoprime";
      _vcallsubn(aTHX_ G_SCALAR, sub, items);
      return; /* skip implicit PUTBACK */
    }

void
_XS_lucas_sequence(IN UV n, IN IV P, IN IV Q, IN UV k)
  PREINIT:
    UV U, V, Qk;
  PPCODE:
    lucas_seq(&U, &V, &Qk,  n, P, Q, k);
    PUSHs(sv_2mortal(newSVuv( U )));    /* 4 args in, 3 out, no EXTEND needed */
    PUSHs(sv_2mortal(newSVuv( V )));
    PUSHs(sv_2mortal(newSVuv( Qk )));

int
_XS_is_prime(IN UV n, ...)
  ALIAS:
    _XS_is_prob_prime = 1
    _XS_is_lucas_pseudoprime = 2
    _XS_is_strong_lucas_pseudoprime = 3
    _XS_is_extra_strong_lucas_pseudoprime = 4
    _XS_is_frobenius_underwood_pseudoprime = 5
    _XS_is_almost_extra_strong_lucas_pseudoprime = 6
    _XS_is_pseudoprime = 7
    _XS_is_aks_prime = 8
    _XS_BPSW = 9
  PREINIT:
    int ret;
  CODE:
    switch (ix) {
      case 0: ret = _XS_is_prime(n); break;
      case 1: ret = _XS_is_prob_prime(n); break;
      case 2: ret = _XS_is_lucas_pseudoprime(n, 0); break;
      case 3: ret = _XS_is_lucas_pseudoprime(n, 1); break;
      case 4: ret = _XS_is_lucas_pseudoprime(n, 2); break;
      case 5: ret = _XS_is_frobenius_underwood_pseudoprime(n); break;
      case 6: ret = _XS_is_almost_extra_strong_lucas_pseudoprime
                    ( n, (items == 1) ? 1 : my_svuv(ST(1)) );  break;
      case 7: ret = _XS_is_pseudoprime(n, my_svuv(ST(1)));     break;
      case 8: ret = _XS_is_aks_prime(n); break;
      default:ret = _XS_BPSW(n); break;
    }
    RETVAL = ret;
  OUTPUT:
    RETVAL


void
is_prime(IN SV* svn)
  ALIAS:
    is_prob_prime = 1
  PREINIT:
    int status;
  PPCODE:
    status = _validate_int(aTHX_ svn, 1);
    if (status == -1) {
      XSRETURN_UV(0);
    } else if (status == 1) {
      UV n = my_svuv(svn);
      XSRETURN_UV(_XS_is_prime(n));
    } else {
      const char* sub = 0;
      if (_XS_get_callgmp())
        sub = (ix == 0) ? "GMP::is_prime"
                        : "GMP::is_prob_prime";
      else
        sub = (ix == 0) ? "_generic_is_prime"
                        : "_generic_is_prob_prime";
      _vcallsub(sub);
      return; /* skip implicit PUTBACK */
    }

void
next_prime(IN SV* svn)
  ALIAS:
    prev_prime = 1
  PPCODE:
    if (_validate_int(aTHX_ svn, 0)) {
      UV n = my_svuv(svn);
      if (ix) {
        XSRETURN_UV( (n < 3) ? 0 : _XS_prev_prime(n));
      } else {
        if (n < _max_prime)
          XSRETURN_UV(_XS_next_prime(n));
      }
    }
    _vcallsub((ix == 0) ?  "_generic_next_prime" :
                           "_generic_prev_prime" );
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
        case 0:  _vcallsubn(aTHX_ gimme_v, "_generic_factor", 1);     break;
        case 1:  _vcallsubn(aTHX_ gimme_v, "_generic_factor_exp", 1); break;
        default: _vcallsubn(aTHX_ gimme_v, "_generic_divisors", 1);   break;
      }
      return; /* skip implicit PUTBACK */
    }


void
divisor_sum(IN SV* svn, ...)
  PPCODE:
    SV* svk = (items > 1) ? ST(1) : 0;
    int nstatus = _validate_int(aTHX_ svn, 0);
    int kstatus = (items == 1 || (SvIOK(svk) && SvIV(svk)))  ?  1  :  0;
    if (nstatus == 1 && kstatus == 1) {
      UV n = my_svuv(svn);
      UV k = (items > 1) ? my_svuv(svk) : 1;
      UV sigma = divisor_sum(n, k);
      if (sigma != 0)  XSRETURN_UV(sigma);   /* sigma 0 means overflow */
    }
    _vcallsubn(aTHX_ G_SCALAR, "_generic_divisor_sum",items);
    return; /* skip implicit PUTBACK */


void
znorder(IN SV* sva, IN SV* svn)
  ALIAS:
    legendre_phi = 1
  PPCODE:
    int astatus = _validate_int(aTHX_ sva, 0);
    int nstatus = _validate_int(aTHX_ svn, 0);
    if (astatus == 1 && nstatus == 1) {
      UV a = my_svuv(sva);
      UV n = my_svuv(svn);
      UV ret;
      switch (ix) {
        case 0:  ret = znorder(a,n);
                 if (ret == 0) XSRETURN_UNDEF;  /* not defined */
                 break;
        default: ret = legendre_phi(a, n);
                 break;
      }
      XSRETURN_UV(ret);
    }
    switch (ix) {
      case 0:  _vcallsubn(aTHX_ G_SCALAR, "_generic_znorder", 2);  break;
      /* TODO: Fixup public PP legendre_phi */
      default: _vcallsubn(aTHX_ G_SCALAR, "PP::_legendre_phi", 2); break;
    }
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
      IV k = (abpositive) ? kronecker_uu(a,b) : kronecker_ss(a,b);
      XSRETURN_IV( k );
    }
    _vcallsubn(aTHX_ G_SCALAR, "_generic_kronecker", 2);
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
        case 0: ret = _XS_ExponentialIntegral(nv); break;
        case 1: ret = _XS_LogarithmicIntegral(nv); break;
        case 2: ret = (double) ld_riemann_zeta(nv); break;
        case 3:
        default:ret = _XS_RiemannR(nv); break;
      }
    } else {
      UV uv = SvUV(x);
      switch (ix) {
        case 4: ret = _XS_chebyshev_theta(uv); break;
        case 5:
        default:ret = _XS_chebyshev_psi(uv); break;
      }
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
        XSRETURN_IV(moebius(n));
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
          for (i = lo; i <= hi; i++)
            PUSHs(sv_2mortal(newSViv(mu[i-lo])));
          Safefree(mu);
        }
      }
    } else {
      /* Whatever we didn't handle above */
      U32 gimme_v = GIMME_V;
      switch (ix) {
        case 0:  _vcallsubn(aTHX_ gimme_v, "_generic_euler_phi", items); break;
        case 1:
        default: _vcallsubn(aTHX_ gimme_v, "_generic_moebius", items);   break;
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
        sv_setuv(svn, n);
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
  PPCODE:
  {
#if !USE_MULTICALL
    _vcallsubn(aTHX_ G_VOID|G_DISCARD, "_generic_forprimes",items);
#else
    UV beg, end;
    GV *gv;
    HV *stash;
    CV *cv;
    SV* svarg;  /* We use svarg to prevent clobbering $_ outside the block */
    void* ctx;
    unsigned char* segment;
    UV seg_base, seg_low, seg_high;

    if (!_validate_int(aTHX_ svbeg, 0) || (items >= 3 && !_validate_int(aTHX_ svend,0))) {
      _vcallsubn(aTHX_ G_VOID|G_DISCARD, "_generic_forprimes",items);
      return;
    }
    if (items < 3) {
      beg = 2;
      end = my_svuv(svbeg);
    } else {
      beg = my_svuv(svbeg);
      end = my_svuv(svend);
    }
    if (beg > end)
      XSRETURN(0);

    cv = sv_2cv(block, &stash, &gv, 0);
    if (cv == Nullcv)
      croak("Not a subroutine reference");
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
    /* For small ranges with large bases or tiny ranges, it will be faster
     * and less memory to just iterate through the primes in range.  The
     * exact limits will change based on the sieve vs. next_prime speed. */
    if (beg <= end && !CvISXSUB(cv) && (
#if BITS_PER_WORD == 64
          (beg >= UVCONST(     100000000000000) && end-beg <    100000) ||
          (beg >= UVCONST(      10000000000000) && end-beg <     40000) ||
          (beg >= UVCONST(       1000000000000) && end-beg <     17000) ||
#endif
          (                                        end-beg <       500) ) ) {
      dMULTICALL;
      I32 gimme = G_VOID;
      PUSH_MULTICALL(cv);
      for (beg = _XS_next_prime(beg-1); beg <= end && beg != 0; beg = _XS_next_prime(beg)) {
        sv_setuv(svarg, beg);
        MULTICALL;
      }
      FIX_MULTICALL_REFCOUNT;
      POP_MULTICALL;
    } else if (beg <= end) {
      ctx = start_segment_primes(beg, end, &segment);
      while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
        if (!CvISXSUB(cv)) {
          dMULTICALL;
          I32 gimme = G_VOID;
          PUSH_MULTICALL(cv);
          START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_low - seg_base, seg_high - seg_base ) {
            sv_setuv(svarg, seg_base + p);
            MULTICALL;
          } END_DO_FOR_EACH_SIEVE_PRIME
          FIX_MULTICALL_REFCOUNT;
          POP_MULTICALL;
        } else {
          START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_low - seg_base, seg_high - seg_base ) {
            sv_setuv(svarg, seg_base + p);
            PUSHMARK(SP);
            call_sv((SV*)cv, G_VOID|G_DISCARD);
          } END_DO_FOR_EACH_SIEVE_PRIME
        }
      }
      end_segment_primes(ctx);
    }
    SvREFCNT_dec(svarg);
#endif
  }

void
forcomposites (SV* block, IN SV* svbeg, IN SV* svend = 0)
  PROTOTYPE: &$;$
  PPCODE:
  {
    UV beg, end;
    GV *gv;
    HV *stash;
    SV* svarg;  /* We use svarg to prevent clobbering $_ outside the block */
    CV *cv = sv_2cv(block, &stash, &gv, 0);

    if (cv == Nullcv)
      croak("Not a subroutine reference");

    if (!_validate_int(aTHX_ svbeg, 0) || (items >= 3 && !_validate_int(aTHX_ svend,0))) {
      _vcallsubn(aTHX_ G_VOID|G_DISCARD, "_generic_forcomposites",items);
      return;
    }

    if (items < 3) {
      beg = 4;
      end = my_svuv(svbeg);
    } else {
      beg = my_svuv(svbeg);
      end = my_svuv(svend);
    }
    if (beg > end)
      return;

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
      /* If beg or end are >= _max_prime, then this will die. */
      prevprime = _XS_prev_prime(beg);
      nextprime = _XS_next_prime(end);
      ctx = start_segment_primes(beg, nextprime, &segment);
      while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
        START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_low - seg_base, seg_high - seg_base ) {
          cbeg = prevprime+1;  if (cbeg < beg) cbeg = beg;
          prevprime = seg_base + p;
          cend = prevprime-1;  if (cend > end) cbeg = end;
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
    {
      beg = (beg <= 4) ? 3 : beg-1;
      while (beg++ < end) {
        if (!_XS_is_prob_prime(beg)) {
          sv_setuv(svarg, beg);
          PUSHMARK(SP);
          call_sv((SV*)cv, G_VOID|G_DISCARD);
        }
      }
    }
    SvREFCNT_dec(svarg);
  }

void
fordivisors (SV* block, IN SV* svn)
  PROTOTYPE: &$
  PPCODE:
  {
    UV i, n, ndivisors;
    UV *divs;
    GV *gv;
    HV *stash;
    SV* svarg;  /* We use svarg to prevent clobbering $_ outside the block */
    CV *cv = sv_2cv(block, &stash, &gv, 0);

    if (cv == Nullcv)
      croak("Not a subroutine reference");

    if (!_validate_int(aTHX_ svn, 0)) {
      _vcallsubn(aTHX_ G_VOID|G_DISCARD, "_generic_fordivisors",2);
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
  }
