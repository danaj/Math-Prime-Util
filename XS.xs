
#define PERL_NO_GET_CONTEXT 1 /* Define at top for more efficiency. */

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "multicall.h"  /* only works in 5.6 and newer */
/* We're not using anything for which we need ppport.h */
#ifndef XSRETURN_UV   /* Er, almost.  Fix 21086 from Sep 2003 */
  #define XST_mUV(i,v)  (ST(i) = sv_2mortal(newSVuv(v))  )
  #define XSRETURN_UV(v) STMT_START { XST_mUV(0,v);  XSRETURN(1); } STMT_END
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
 #define set_val_from_sv(val, sv) \
   do { val = PSTRTOULL(SvPV_nolen(sv), NULL, 10); } while(0)
 #define set_sval_from_sv(val, sv) \
   do { val = PSTRTOLL(SvPV_nolen(sv), NULL, 10); } while(0)
#else
 #define set_val_from_sv(val, sv) \
   do { val = (!SvROK(sv)) ? SvUV(sv) : PSTRTOULL(SvPV_nolen(sv), NULL, 10); } while(0)
 #define set_sval_from_sv(val, sv) \
   do { val = (!SvROK(sv)) ? SvIV(sv) : PSTRTOLL(SvPV_nolen(sv), NULL, 10); } while(0)
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


static int pbrent_factor_a1(UV n, UV *factors, UV maxrounds) {
  return pbrent_factor(n, factors, maxrounds, 1);
}

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
static int _validate_int(SV* n, int negok)
{
  dTHX;
  const char* maxstr;
  char* ptr;
  STRLEN i, len, maxlen;
  int ret, isneg = 0;

  if (!SvOK(n))  croak("Parameter must be defined");
  if (SvIOK(n)) {                      /* If defined as number, use that */
    if (SvIsUV(n) || SvIV(n) >= 0)  return 1;
    if (negok)  return -1;
    else      croak("Parameter '%s' must be a positive integer", SvPV_nolen(n));
  }
  if (SvROK(n) && !sv_isa(n, "Math::BigInt"))  return 0;
  ptr = SvPV(n, len);                  /* Includes stringifying bigints */
  if (len == 0 || ptr == 0)  croak("Parameter '' must be a positive integer");
  if (ptr[0] == '-') {                 /* Read negative sign */
    if (negok) { isneg = 1; ptr++; len--; }
    else       croak("Parameter '%s' must be a positive integer", ptr);
  }
  if (ptr[0] == '+') { ptr++; len--; } /* Allow a single plus sign */
  while (len > 0 && *ptr == '0')       /* Strip all leading zeros */
    { ptr++; len--; }
  if (len > uvmax_maxlen)              /* Huge number, don't even look at it */
    return 0;
  for (i = 0; i < len; i++)            /* Ensure all characters are digits */
    if (!isDIGIT(ptr[i]))
      croak("Parameter '%s' must be a positive integer", ptr);
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

/* Call a Perl sub to handle work for us.
 *   The input is a single SV on the top of the stack.
 *   The output is a single mortal SV that is on the stack.
 */
static void _vcallsub(const char* name)
{
  dTHX;
  dSP;                               /* Local copy of stack pointer         */
  int count;
  SV* arg;

  ENTER;                             /* Start wrapper                       */
  SAVETMPS;                          /* Start (2)                           */

  arg = POPs;                        /* Get argument value from stack       */
  PUSHMARK(SP);                      /* Start args: note our SP             */
  XPUSHs(arg);
  PUTBACK;                           /* End args:   set global SP to ours   */

  count = call_pv(name, G_SCALAR);   /* Call the sub                        */
  SPAGAIN;                           /* refresh local stack pointer         */

  if (count != 1)
    croak("callback sub should return one value");

  TOPs = SvREFCNT_inc(TOPs);         /* Make sure FREETMPS doesn't kill it  */

  PUTBACK;
  FREETMPS;                          /* End wrapper                         */
  LEAVE;                             /* End (2)                             */
  TOPs = sv_2mortal(TOPs);           /* mortalize it.  refcnt will be 1.    */
}

#if BITS_PER_WORD == 64
static const UV _max_prime = UVCONST(18446744073709551557);
#else
static const UV _max_prime = UVCONST(4294967291);
#endif


MODULE = Math::Prime::Util	PACKAGE = Math::Prime::Util

PROTOTYPES: ENABLE


void
_XS_set_verbose(IN int verbose)

int
_XS_get_verbose()

void
_XS_set_callgmp(IN int call_gmp)

int
_XS_get_callgmp()


void
prime_precalc(IN UV n)

void
prime_memfree()

void
_prime_memfreeall()

UV
_XS_prime_count(IN UV low, IN UV high = 0)
  CODE:
    if (high == 0) {   /* Without a Perl layer in front of this, we'll have */
      high = low;      /* the pathological case of a-0 turning into 0-a.    */
      low = 0;
    }
    if (GIMME_V == G_VOID) {
      prime_precalc(high);
      RETVAL = 0;
    } else {
      RETVAL = _XS_prime_count(low, high);
    }
  OUTPUT:
    RETVAL

UV
_XS_nth_prime(IN UV n)
  ALIAS:
    _XS_next_prime = 1
    _XS_prev_prime = 2
    _XS_legendre_pi = 3
    _XS_meissel_pi = 4
    _XS_lehmer_pi = 5
    _XS_LMO_pi = 6
    _XS_LMOS_pi = 7
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
      case 6: ret = _XS_LMO_pi(n); break;
      case 7: ret = _XS_LMOS_pi(n); break;
      default: croak("_XS_nth_prime: Unknown function alias"); break;
    }
    RETVAL = ret;
  OUTPUT:
    RETVAL

UV
_XS_divisor_sum(IN UV n, IN UV k)

UV
_XS_legendre_phi(IN UV x, IN UV a)


UV
_get_prime_cache_size()
  CODE:
    RETVAL = get_prime_cache(0, 0);
  OUTPUT:
    RETVAL

int
_XS_prime_maxbits()
  CODE:
    RETVAL = BITS_PER_WORD;
  OUTPUT:
    RETVAL


SV*
sieve_primes(IN UV low, IN UV high)
  PREINIT:
    AV* av = newAV();
  CODE:
    if (low <= high) {
      START_DO_FOR_EACH_PRIME(low, high) {
        av_push(av,newSVuv(p));
      } END_DO_FOR_EACH_PRIME
    }
    RETVAL = newRV_noinc( (SV*) av );
  OUTPUT:
    RETVAL


SV*
trial_primes(IN UV low, IN UV high)
  PREINIT:
    UV  p;
    AV* av = newAV();
  CODE:
    if (low <= high) {
      if (low >= 2) low--;   /* Make sure low gets included */
      for (p = _XS_next_prime(low); p <= high && p != 0; p = _XS_next_prime(p)) {
        av_push(av,newSVuv(p));
      }
    }
    RETVAL = newRV_noinc( (SV*) av );
  OUTPUT:
    RETVAL

SV*
segment_primes(IN UV low, IN UV high);
  PREINIT:
    AV* av = newAV();
  CODE:
    if ((low <= 2) && (high >= 2)) { av_push(av, newSVuv( 2 )); }
    if ((low <= 3) && (high >= 3)) { av_push(av, newSVuv( 3 )); }
    if ((low <= 5) && (high >= 5)) { av_push(av, newSVuv( 5 )); }
    if (low < 7)  low = 7;
    if (low <= high) {
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
    RETVAL = newRV_noinc( (SV*) av );
  OUTPUT:
    RETVAL

SV*
erat_primes(IN UV low, IN UV high)
  PREINIT:
    unsigned char* sieve;
    AV* av = newAV();
  CODE:
    if (low <= high) {
      sieve = sieve_erat30(high);
      if (sieve == 0) {
        croak("Could not generate sieve for %"UVuf, high);
      } else {
        if ((low <= 2) && (high >= 2)) { av_push(av, newSVuv( 2 )); }
        if ((low <= 3) && (high >= 3)) { av_push(av, newSVuv( 3 )); }
        if ((low <= 5) && (high >= 5)) { av_push(av, newSVuv( 5 )); }
        if (low < 7) { low = 7; }
        START_DO_FOR_EACH_SIEVE_PRIME( sieve, low, high ) {
           av_push(av,newSVuv(p));
        } END_DO_FOR_EACH_SIEVE_PRIME
        Safefree(sieve);
      }
    }
    RETVAL = newRV_noinc( (SV*) av );
  OUTPUT:
    RETVAL


#define SIMPLE_FACTOR(func, n, arg1) \
    if (n <= 1) { \
      if (n == 0) \
        XPUSHs(sv_2mortal(newSVuv( n ))); \
    } else { \
      while ( (n% 2) == 0 ) {  n /=  2;  XPUSHs(sv_2mortal(newSVuv( 2 ))); } \
      while ( (n% 3) == 0 ) {  n /=  3;  XPUSHs(sv_2mortal(newSVuv( 3 ))); } \
      while ( (n% 5) == 0 ) {  n /=  5;  XPUSHs(sv_2mortal(newSVuv( 5 ))); } \
      if (n == 1) {  /* done */ } \
      else if (_XS_is_prime(n)) { XPUSHs(sv_2mortal(newSVuv( n ))); } \
      else { \
        UV factors[MPU_MAX_FACTORS+1]; \
        int i, nfactors; \
        nfactors = func(n, factors, arg1); \
        for (i = 0; i < nfactors; i++) { \
          XPUSHs(sv_2mortal(newSVuv( factors[i] ))); \
        } \
      } \
    }
#define SIMPLE_FACTOR_2ARG(func, n, arg1, arg2) \
    /* Stupid MSVC won't bring its C compiler out of the 1980s. */ \
    if (n <= 1) { \
      if (n == 0) \
        XPUSHs(sv_2mortal(newSVuv( n ))); \
    } else { \
      while ( (n% 2) == 0 ) {  n /=  2;  XPUSHs(sv_2mortal(newSVuv( 2 ))); } \
      while ( (n% 3) == 0 ) {  n /=  3;  XPUSHs(sv_2mortal(newSVuv( 3 ))); } \
      while ( (n% 5) == 0 ) {  n /=  5;  XPUSHs(sv_2mortal(newSVuv( 5 ))); } \
      if (n == 1) {  /* done */ } \
      else if (_XS_is_prime(n)) { XPUSHs(sv_2mortal(newSVuv( n ))); } \
      else { \
        UV factors[MPU_MAX_FACTORS+1]; \
        int i, nfactors; \
        nfactors = func(n, factors, arg1, arg2); \
        for (i = 0; i < nfactors; i++) { \
          XPUSHs(sv_2mortal(newSVuv( factors[i] ))); \
        } \
      } \
    }

void
_XS_factor(IN UV n)
  PREINIT:
    UV factors[MPU_MAX_FACTORS+1];
    int i, nfactors;
  PPCODE:
    nfactors = factor(n, factors);
    if (GIMME_V == G_SCALAR) {
      PUSHs(sv_2mortal(newSVuv(nfactors)));
    } else if (GIMME_V == G_ARRAY) {
      EXTEND(SP, nfactors);
      for (i = 0; i < nfactors; i++) {
        PUSHs(sv_2mortal(newSVuv( factors[i] )));
      }
    }

void
_XS_factor_exp(IN UV n)
  PREINIT:
    UV factors[MPU_MAX_FACTORS+1];
    UV exponents[MPU_MAX_FACTORS+1];
    int i, nfactors;
  PPCODE:
    nfactors = factor_exp(n, factors, exponents);
    if (GIMME_V == G_SCALAR) {
      PUSHs(sv_2mortal(newSVuv(nfactors)));
    } else {
      /* Return ( [p1,e1], [p2,e2], [p3,e3], ... ) */
      if (n == 1)  XSRETURN_EMPTY;
      for (i = 0; i < nfactors; i++) {
        AV* av = newAV();
        av_push(av, newSVuv(factors[i]));
        av_push(av, newSVuv(exponents[i]));
        XPUSHs( sv_2mortal(newRV_noinc( (SV*) av )) );
      }
    }

void
_XS_divisors(IN UV n)
  PPCODE:
    if (GIMME_V == G_SCALAR) {
      PUSHs(sv_2mortal(newSVuv( _XS_divisor_sum(n, 0) )));
    } else {
      UV i, ndivisors;
      UV* divs = _divisor_list(n, &ndivisors);
      for (i = 0; i < ndivisors; i++)
        XPUSHs(sv_2mortal(newSVuv(divs[i])));
      Safefree(divs);
    }

void
trial_factor(IN UV n, IN UV maxfactor = 0)
  PPCODE:
    SIMPLE_FACTOR(trial_factor, n, maxfactor);

void
fermat_factor(IN UV n, IN UV maxrounds = 64*1024*1024)
  PPCODE:
    SIMPLE_FACTOR(fermat_factor, n, maxrounds);

void
holf_factor(IN UV n, IN UV maxrounds = 8*1024*1024)
  PPCODE:
    SIMPLE_FACTOR(holf_factor, n, maxrounds);

void
squfof_factor(IN UV n, IN UV maxrounds = 4*1024*1024)
  PPCODE:
    SIMPLE_FACTOR(squfof_factor, n, maxrounds);

void
rsqufof_factor(IN UV n, IN UV maxrounds = 4*1024*1024)
  PPCODE:
    SIMPLE_FACTOR(racing_squfof_factor, n, maxrounds);

void
pbrent_factor(IN UV n, IN UV maxrounds = 4*1024*1024)
  PPCODE:
    SIMPLE_FACTOR(pbrent_factor_a1, n, maxrounds);

void
prho_factor(IN UV n, IN UV maxrounds = 4*1024*1024)
  PPCODE:
    SIMPLE_FACTOR(prho_factor, n, maxrounds);

void
pminus1_factor(IN UV n, IN UV B1 = 1*1024*1024, IN UV B2 = 0)
  PPCODE:
    if (B2 == 0)
      B2 = 10*B1;
    SIMPLE_FACTOR_2ARG(pminus1_factor, n, B1, B2);

void
pplus1_factor(IN UV n, IN UV B = 200)
  PPCODE:
    SIMPLE_FACTOR(pplus1_factor, n, B);

int
_XS_is_pseudoprime(IN UV n, IN UV a)

int
_XS_miller_rabin(IN UV n, ...)
  PREINIT:
    UV bases[64];
    int prob_prime = 1;
    int c = 1;
  CODE:
    if (items < 2)
      croak("No bases given to miller_rabin");
    if ( (n == 0) || (n == 1) ) XSRETURN_IV(0);   /* 0 and 1 are composite */
    if ( (n == 2) || (n == 3) ) XSRETURN_IV(1);   /* 2 and 3 are prime */
    if (( n % 2 ) == 0)  XSRETURN_IV(0);          /* MR works with odd n */
    while (c < items) {
      int b = 0;
      while (c < items) {
        bases[b++] = SvUV(ST(c));
        c++;
        if (b == 64) break;
      }
      prob_prime = _XS_miller_rabin(n, bases, b);
      if (prob_prime != 1)
        break;
    }
    RETVAL = prob_prime;
  OUTPUT:
    RETVAL

void
_XS_lucas_sequence(IN UV n, IN IV P, IN IV Q, IN UV k)
  PREINIT:
    UV U, V, Qk;
  PPCODE:
    lucas_seq(&U, &V, &Qk,  n, P, Q, k);
    XPUSHs(sv_2mortal(newSVuv( U )));
    XPUSHs(sv_2mortal(newSVuv( V )));
    XPUSHs(sv_2mortal(newSVuv( Qk )));

int
_XS_is_prime(IN UV n)
  ALIAS:
    _XS_is_prob_prime = 1
    _XS_is_lucas_pseudoprime = 2
    _XS_is_strong_lucas_pseudoprime = 3
    _XS_is_extra_strong_lucas_pseudoprime = 4
    _XS_is_frobenius_underwood_pseudoprime = 5
    _XS_is_aks_prime = 6
    _XS_BPSW = 7
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
      case 6: ret = _XS_is_aks_prime(n); break;
      case 7: ret = _XS_BPSW(n); break;
      default: croak("_XS_is_prime: Unknown function alias"); break;
    }
    RETVAL = ret;
  OUTPUT:
    RETVAL

int
_XS_is_almost_extra_strong_lucas_pseudoprime(IN UV n, IN UV increment = 1)
  CODE:
    RETVAL = _XS_is_almost_extra_strong_lucas_pseudoprime(n, increment);
  OUTPUT:
    RETVAL

int
is_prime(IN SV* svn)
  ALIAS:
    is_prob_prime = 1
  PREINIT:
    int status;
  PPCODE:
    status = _validate_int(svn, 1);
    if (status == -1) {
      XSRETURN_UV(0);
    } else if (status == 1) {
      UV n;
      set_val_from_sv(n, svn);
      XSRETURN_UV(_XS_is_prime(n));
    } else {
      const char* sub = 0;
      if (_XS_get_callgmp())
        sub = (ix == 0) ? "Math::Prime::Util::GMP::is_prime"
                        : "Math::Prime::Util::GMP::is_prob_prime";
      else
        sub = (ix == 0) ? "Math::Prime::Util::_generic_is_prime"
                        : "Math::Prime::Util::_generic_is_prob_prime";
      _vcallsub(sub);
      XSRETURN(1);
    }

void
next_prime(IN SV* n)
  ALIAS:
    prev_prime = 1
  PREINIT:
    UV val;
  PPCODE:
    if (ix) {
      if (_validate_int(n, 0)) {
        set_val_from_sv(val, n);
        XSRETURN_UV( (val < 3) ? 0 : _XS_prev_prime(val));
      }
    } else {
      if (_validate_int(n, 0)) {
        set_val_from_sv(val, n);
        if (val < _max_prime)
          XSRETURN_UV(_XS_next_prime(val));
      }
    }
    _vcallsub((ix == 0) ?  "Math::Prime::Util::_generic_next_prime" :
                           "Math::Prime::Util::_generic_prev_prime" );
    XSRETURN(1);

void
znprimroot(IN SV* svn)
  PREINIT:
    int status;
  PPCODE:
    status = _validate_int(svn, 1);
    if (status != 0) {
      UV n, r;
      set_val_from_sv(n, svn);
      if (status == -1)
        n = -(IV)n;
      r = znprimroot(n);
      if (r == 0 && n != 1) XSRETURN_EMPTY;
      XSRETURN_UV(r);
    }
    _vcallsub("Math::Prime::Util::_generic_znprimroot");
    XSRETURN(1);

int
kronecker(IN SV* sva, IN SV* svb)
  PREINIT:
    int astatus, bstatus;
  CODE:
    astatus = _validate_int(sva, 2);
    bstatus = _validate_int(svb, 2);
    if (astatus == 1 && bstatus == 1) {
      UV a, b;
      set_val_from_sv(a, sva);
      set_val_from_sv(b, svb);
      RETVAL = kronecker_uu(a, b);
    } else if (astatus != 0 && SvIOK(sva) && !SvIsUV(sva) &&
               bstatus != 0 && SvIOK(svb) && !SvIsUV(svb) ) {
      IV a, b;
      set_sval_from_sv(a, sva);
      set_sval_from_sv(b, svb);
      RETVAL = kronecker_ss(a, b);
    } else {
      dTHX;
      dSP;
      int count;
      ENTER;
      SAVETMPS;
      PUSHMARK(SP);
      XPUSHs(sva);
      XPUSHs(svb);
      PUTBACK;
      count = call_pv("Math::Prime::Util::_generic_kronecker", G_SCALAR);
      SPAGAIN;
      if (count != 1) croak("callback sub should return one value");
      RETVAL = POPi;
      PUTBACK;
      FREETMPS;
      LEAVE;
    }
  OUTPUT:
    RETVAL

double
_XS_ExponentialIntegral(IN double x)
  ALIAS:
    _XS_LogarithmicIntegral = 1
    _XS_RiemannZeta = 2
    _XS_RiemannR = 3
  PREINIT:
    double ret;
  CODE:
    switch (ix) {
      case 0: ret = _XS_ExponentialIntegral(x); break;
      case 1: ret = _XS_LogarithmicIntegral(x); break;
      case 2: ret = (double) ld_riemann_zeta(x); break;
      case 3: ret = _XS_RiemannR(x); break;
      default: croak("_XS_ExponentialIntegral: Unknown function alias"); break;
    }
    RETVAL = ret;
  OUTPUT:
    RETVAL

double
_XS_chebyshev_theta(IN UV n)
  ALIAS:
    _XS_chebyshev_psi = 1
  PREINIT:
    double ret;
  CODE:
    switch (ix) {
      case 0: ret = _XS_chebyshev_theta(n); break;
      case 1: ret = _XS_chebyshev_psi(n); break;
      default: croak("_XS_chebyshev_theta: Unknown function alias"); break;
    }
    RETVAL = ret;
  OUTPUT:
    RETVAL

void
_XS_totient(IN UV lo, IN UV hi = 0)
  PREINIT:
    UV i;
  PPCODE:
    if (hi != lo && hi != 0) { /* Totients in a range, returns array */
      UV* totients;
      if (hi < lo) XSRETURN_EMPTY;
      if (lo < 2) {
        if (lo <= 0           ) XPUSHs(sv_2mortal(newSVuv(0)));
        if (lo <= 1 && hi >= 1) XPUSHs(sv_2mortal(newSVuv(1)));
        lo = 2;
      }
      if (hi >= lo) {
        totients = _totient_range(lo, hi);
        /* Extend the stack to handle how many items we'll return */
        EXTEND(SP, hi-lo+1);
        for (i = lo; i <= hi; i++)
          PUSHs(sv_2mortal(newSVuv(totients[i-lo])));
        Safefree(totients);
      }
    } else {
      XSRETURN_UV(totient(lo));
    }

void
carmichael_lambda(IN SV* svn)
  PREINIT:
    int status;
  PPCODE:
    status = _validate_int(svn, 0);
    if (status == 1) {
      UV n;
      set_val_from_sv(n, svn);
      XSRETURN_UV(carmichael_lambda(n));
    }
    _vcallsub("Math::Prime::Util::_generic_carmichael_lambda");
    XSRETURN(1);

void
_XS_moebius(IN UV lo, IN UV hi = 0)
  PREINIT:
    UV i;
  PPCODE:
    if (hi != lo && hi != 0) {   /* mobius in a range */
      signed char* mu = _moebius_range(lo, hi);
      MPUassert( mu != 0, "_moebius_range returned 0" );
      EXTEND(SP, hi-lo+1);
      for (i = lo; i <= hi; i++)
        PUSHs(sv_2mortal(newSViv(mu[i-lo])));
      Safefree(mu);
    } else {
      UV factors[MPU_MAX_FACTORS+1];
      UV nfactors;
      UV n = lo;

      if (n <= 1)
        XSRETURN_IV(n);

      if ( (!(n% 4) && n >=  4) || (!(n% 9) && n >=  9) ||
           (!(n%25) && n >= 25) || (!(n%49) && n >= 49) )
        XSRETURN_IV(0);

      nfactors = factor(n, factors);
      for (i = 1; i < nfactors; i++)
        if (factors[i] == factors[i-1])
          XSRETURN_IV(0);

      XSRETURN_IV( (nfactors % 2) ? -1 : 1 );
    }

IV
_XS_mertens(IN UV n)

void
exp_mangoldt(IN SV* svn)
  PREINIT:
    int status;
    UV n;
  PPCODE:
    status = _validate_int(svn, 1);
    if (status == -1) {
      XSRETURN_UV(1);
    } else if (status == 1) {
      set_val_from_sv(n, svn);
      if (n <= 1)
        XSRETURN_UV(1);
      else if ((n & (n-1)) == 0)  /* Power of 2 */
        XSRETURN_UV(2);
      else if ((n & 1) == 0)      /* Even number */
        XSRETURN_UV(1);
      else {
        UV factors[MPU_MAX_FACTORS+1];
        UV nfactors, i;
        /* We could try a partial factor, e.g. looking for two small factors */
        /* We could also check powers of primes searching for n */
        nfactors = factor(n, factors);
        for (i = 1; i < nfactors; i++) {
          if (factors[i] != factors[0])
            XSRETURN_UV(1);
        }
        XSRETURN_UV(factors[0]);
      }
    } else {
      _vcallsub("Math::Prime::Util::_generic_exp_mangoldt");
      XSRETURN(1);
    }

int
_validate_num(SV* n, ...)
  CODE:
    RETVAL = 0;
    if (_validate_int(n, 0)) {
      if (SvROK(n)) {  /* Convert small Math::BigInt object into scalar */
        UV val;
        set_val_from_sv(val, n);
        sv_setuv(n, val);
      }
      if (items > 1 && SvOK(ST(1))) {
        UV val, min, max;
        set_val_from_sv(val, n);
        set_val_from_sv(min, ST(1));
        if (val < min)
          croak("Parameter '%"UVuf"' must be >= %"UVuf, val, min);
        if (items > 2 && SvOK(ST(2))) {
          set_val_from_sv(max, ST(2));
          if (val > max)
            croak("Parameter '%"UVuf"' must be <= %"UVuf, val, max);
          MPUassert( items <= 3, "_validate_num takes at most 3 parameters");
        }
      }
      RETVAL = 1;
    }
  OUTPUT:
    RETVAL

void
forprimes (SV* block, IN SV* svbeg, IN SV* svend = 0)
  PROTOTYPE: &$;$
  CODE:
  {
#if !USE_MULTICALL
    dSP;
    PUSHMARK(SP);
    XPUSHs(block); XPUSHs(svbeg); XPUSHs(svend);
    PUTBACK;
    (void) call_pv("Math::Prime::Util::_generic_forprimes", G_VOID|G_DISCARD);
    SPAGAIN;
#else
    UV beg, end;
    GV *gv;
    HV *stash;
    CV *cv;
    SV* svarg;  /* We use svarg to prevent clobbering $_ outside the block */
    void* ctx;
    unsigned char* segment;
    UV seg_base, seg_low, seg_high;

    if (!_validate_int(svbeg, 0) || (items >= 3 && !_validate_int(svend,0))) {
      dSP;
      PUSHMARK(SP);
      XPUSHs(block); XPUSHs(svbeg); XPUSHs(svend);
      PUTBACK;
      (void) call_pv("Math::Prime::Util::_generic_forprimes", G_VOID|G_DISCARD);
      SPAGAIN;
      XSRETURN_UNDEF;
    }
    if (items < 3) {
      beg = 2;
      set_val_from_sv(end, svbeg);
    } else {
      set_val_from_sv(beg, svbeg);
      set_val_from_sv(end, svend);
    }
    if (beg > end)
      XSRETURN_UNDEF;

    cv = sv_2cv(block, &stash, &gv, 0);
    if (cv == Nullcv)
      croak("Not a subroutine reference");
    SAVESPTR(GvSV(PL_defgv));
    svarg = newSVuv(0);
    /* Handle early part */
    while (beg < 6) {
      dSP;
      beg = (beg <= 2) ? 2 : (beg <= 3) ? 3 : 5;
      if (beg <= end) {
        sv_setuv(svarg, beg);
        GvSV(PL_defgv) = svarg;
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
        GvSV(PL_defgv) = svarg;
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
            GvSV(PL_defgv) = svarg;
            MULTICALL;
          } END_DO_FOR_EACH_SIEVE_PRIME
          FIX_MULTICALL_REFCOUNT;
          POP_MULTICALL;
        } else {
          START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_low - seg_base, seg_high - seg_base ) {
            dSP;
            sv_setuv(svarg, seg_base + p);
            GvSV(PL_defgv) = svarg;
            PUSHMARK(SP);
            call_sv((SV*)cv, G_VOID|G_DISCARD);
          } END_DO_FOR_EACH_SIEVE_PRIME
        }
      }
      end_segment_primes(ctx);
    }
    SvREFCNT_dec(svarg);
#endif
    XSRETURN_UNDEF;
  }

void
forcomposites (SV* block, IN SV* svbeg, IN SV* svend = 0)
  PROTOTYPE: &$;$
  CODE:
  {
    UV beg, end;
    GV *gv;
    HV *stash;
    SV* svarg;  /* We use svarg to prevent clobbering $_ outside the block */
    CV *cv = sv_2cv(block, &stash, &gv, 0);

    if (cv == Nullcv)
      croak("Not a subroutine reference");
    if (items <= 1) XSRETURN_UNDEF;

    if (!_validate_int(svbeg, 0) || (items >= 3 && !_validate_int(svend,0))) {
      dSP;
      PUSHMARK(SP);
      XPUSHs(block); XPUSHs(svbeg); XPUSHs(svend);
      PUTBACK;
      (void) call_pv("Math::Prime::Util::_generic_forcomposites", G_VOID|G_DISCARD);
      SPAGAIN;
      XSRETURN_UNDEF;
    }

    if (items < 3) {
      beg = 4;
      set_val_from_sv(end, svbeg);
    } else {
      set_val_from_sv(beg, svbeg);
      set_val_from_sv(end, svend);
    }
    if (beg > end)
      XSRETURN_UNDEF;

    SAVESPTR(GvSV(PL_defgv));
    svarg = newSVuv(0);
#if USE_MULTICALL
    if (!CvISXSUB(cv) && (end-beg) > 200) {
      unsigned char* segment;
      UV seg_base, seg_low, seg_high, c, cbeg, cend, prevprime, nextprime;
      void* ctx;
      dMULTICALL;
      I32 gimme = G_VOID;
      PUSH_MULTICALL(cv);
      if (beg <= 4) { /* sieve starts at 7, so handle this here */
        sv_setuv(svarg, 4);  GvSV(PL_defgv) = svarg;  MULTICALL;
        beg = 6;
      }
      /* Find the two primes that bound their interval. */
      /* If beg or end are >= _max_prime, then this will die. */
      prevprime = _XS_prev_prime(beg);
      nextprime = _XS_next_prime(end);
      ctx = start_segment_primes(beg, nextprime, &segment);
      while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
        START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_low - seg_base, seg_high - seg_base ) {
          cbeg = (prevprime+1 < beg)   ? beg : prevprime+1;
          prevprime = seg_base + p;
          cend = (prevprime-1 > end)   ? end : prevprime-1;
          for (c = cbeg; c <= cend; c++) {
            sv_setuv(svarg, c);  GvSV(PL_defgv) = svarg;  MULTICALL;
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
          dSP;
          sv_setuv(svarg, beg);
          GvSV(PL_defgv) = svarg;
          PUSHMARK(SP);
          call_sv((SV*)cv, G_VOID|G_DISCARD);
        }
      }
    }
    SvREFCNT_dec(svarg);
    XSRETURN_UNDEF;
  }

void
fordivisors (SV* block, IN SV* svn)
  PROTOTYPE: &$
  CODE:
  {
    UV i, n, ndivisors;
    UV *divs;
    GV *gv;
    HV *stash;
    SV* svarg;  /* We use svarg to prevent clobbering $_ outside the block */
    CV *cv = sv_2cv(block, &stash, &gv, 0);

    if (cv == Nullcv)
      croak("Not a subroutine reference");
    if (items <= 1) XSRETURN_UNDEF;

    if (!_validate_int(svn, 0)) {
      dSP;
      PUSHMARK(SP);
      XPUSHs(block); XPUSHs(svn);
      PUTBACK;
      (void) call_pv("Math::Prime::Util::_generic_fordivisors", G_VOID|G_DISCARD);
      SPAGAIN;
      XSRETURN_UNDEF;
    }

    set_val_from_sv(n, svn);
    divs = _divisor_list(n, &ndivisors);

    SAVESPTR(GvSV(PL_defgv));
    svarg = newSVuv(0);
#if USE_MULTICALL
    if (!CvISXSUB(cv)) {
      dMULTICALL;
      I32 gimme = G_VOID;
      PUSH_MULTICALL(cv);
      for (i = 0; i < ndivisors; i++) {
        sv_setuv(svarg, divs[i]);
        GvSV(PL_defgv) = svarg;
        MULTICALL;
      }
      FIX_MULTICALL_REFCOUNT;
      POP_MULTICALL;
    }
    else
#endif
    {
      for (i = 0; i < ndivisors; i++) {
        dSP;
        sv_setuv(svarg, divs[i]);
        GvSV(PL_defgv) = svarg;
        PUSHMARK(SP);
        call_sv((SV*)cv, G_VOID|G_DISCARD);
      }
    }
    SvREFCNT_dec(svarg);
    Safefree(divs);
    XSRETURN_UNDEF;
  }
