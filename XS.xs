
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
#include "aks.h"

/* Workaround perl 5.6 UVs */
#if PERL_REVISION <= 5 && PERL_VERSION <= 6 && BITS_PER_WORD == 64
 /* This could be blown up with a wacky string, but it's just for 5.6 */
 #define set_val_from_sv(val, sv) \
   { char*ptr = SvPV_nolen(sv); val = Strtoul(ptr, NULL, 10); }
#else
 #define set_val_from_sv(val, sv) \
   val = SvUV(sv)
#endif

/* multicall compatibility stuff */
#if PERL_REVISION <= 5 && PERL_VERSION < 7
# define USE_MULTICALL 0   /* Too much trouble to work around it */
#else
# define USE_MULTICALL 1
#endif

#if PERL_VERSION < 13 || (PERL_VERSION == 13 && PERL_SUBVERSION < 9)
#  define PERL_HAS_BAD_MULTICALL_REFCOUNT
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

/* Is this a pedantically valid integer?
 * Croaks if undefined or invalid.
 * Returns 0 if it is an object or a string too large for a UV.
 * Returns 1 if it is good to process by XS.
 */
static int _validate_int(SV* n, int negok)
{
  dTHX;
  char* ptr;
  STRLEN i, len;
  UV val;
  int isneg = 0;

  if (!SvOK(n))  croak("Parameter must be defined");
  /* aside: to detect bigint: if ( SvROK(n) && sv_isa(n, "Math::BigInt") ) */
  if (SvROK(n))  return 0;
  /* Perhaps SvPVbyte, or other UTF8 stuff? */
  ptr = SvPV(n, len);
  if (len == 0)  croak("Parameter '' must be a positive integer");
  for (i = 0; i < len; i++) {
    if (!isDIGIT(ptr[i])) {
      if (i == 0 && ptr[i] == '-' && negok)
        isneg = 1;
      else if (i == 0 && ptr[i] == '+')
        /* Allowed */ ;
      else
        croak("Parameter '%s' must be a positive integer", ptr); /* TODO NULL */
    }
  }
  if (isneg) return -1;  /* It's a valid negative number */
  set_val_from_sv(val, n);
  if (val == UV_MAX) { /* Could be bigger than UV_MAX.  Need to find out. */
    char vstr[40];
    sprintf(vstr, "%"UVuf, val);
    /* Skip possible leading zeros */
    while (len > 0 && (*ptr == '0' || *ptr == '+'))
      { ptr++; len--; }
    for (i = 0; i < len; i++)
      if (vstr[i] != ptr[i])
        return 0;
  }
  return 1;
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
      default: croak("_XS_nth_prime: Unknown function alias"); break;
    }
    RETVAL = ret;
  OUTPUT:
    RETVAL

UV
_XS_divisor_sum(IN UV n, IN UV k)


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
    UV  curprime;
    AV* av = newAV();
  CODE:
    if (low <= high) {
      if (low >= 2) low--;   /* Make sure low gets included */
      curprime = _XS_next_prime(low);
      while (curprime <= high) {
        av_push(av,newSVuv(curprime));
        curprime = _XS_next_prime(curprime);
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
    /* Could rewrite using {start/next/end}_segment_primes functions */
    if ((low <= 2) && (high >= 2)) { av_push(av, newSVuv( 2 )); }
    if ((low <= 3) && (high >= 3)) { av_push(av, newSVuv( 3 )); }
    if ((low <= 5) && (high >= 5)) { av_push(av, newSVuv( 5 )); }
    if (low < 7)  low = 7;
    if (low <= high) {
      /* Call the segment siever one or more times */
      UV low_d, high_d, segment_size;
      unsigned char* sieve = get_prime_segment(&segment_size);
      if (sieve == 0)
        croak("Could not get segment cache");

      /* To protect vs. overflow, work entirely with d. */
      low_d  = low  / 30;
      high_d = high / 30;

      {  /* Avoid recalculations of this */
        UV endp = (high_d >= (UV_MAX/30))  ?  UV_MAX-2  :  30*high_d+29;
        prime_precalc(isqrt(endp) + 1 );
      }

      while ( low_d <= high_d ) {
        UV seghigh_d = ((high_d - low_d) < segment_size)
                       ? high_d
                       : (low_d + segment_size-1);
        UV range_d = seghigh_d - low_d + 1;
        UV seghigh = (seghigh_d == high_d) ? high : (seghigh_d*30+29);
        UV segbase = low_d * 30;
        /* printf("  startd = %"UVuf"  endd = %"UVuf"\n", startd, endd); */

        MPUassert( seghigh_d >= low_d, "segment_primes highd < lowd");
        MPUassert( range_d <= segment_size, "segment_primes range > segment size");

        /* Sieve from startd*30+1 to endd*30+29.  */
        if (sieve_segment(sieve, low_d, seghigh_d) == 0) {
          release_prime_segment(sieve);
          croak("Could not segment sieve from %"UVuf" to %"UVuf, segbase+1, seghigh);
        }

        START_DO_FOR_EACH_SIEVE_PRIME( sieve, low - segbase, seghigh - segbase )
          av_push(av,newSVuv( segbase + p ));
        END_DO_FOR_EACH_SIEVE_PRIME

        low_d += range_d;
        low = seghigh+2;
      }
      release_prime_segment(sieve);
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
    EXTEND(SP, nfactors);
    for (i = 0; i < nfactors; i++) {
      PUSHs(sv_2mortal(newSVuv( factors[i] )));
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
is_prime(IN SV* n)
  ALIAS:
    is_prob_prime = 1
  PREINIT:
    int status;
  PPCODE:
    status = _validate_int(n, 1);
    if (status == -1) {
      XSRETURN_UV(0);
    } else if (status == 1) {
      UV val;
      set_val_from_sv(val, n);
      XSRETURN_UV(_XS_is_prime(val));
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
  PPCODE:
    if (_validate_int(n, 0)) {
      UV val;
      set_val_from_sv(val, n);
      if ( (ix && val < 3) || (!ix && val >= _max_prime) )  XSRETURN_UV(0);
      if (ix) XSRETURN_UV(_XS_prev_prime(val));
      else    XSRETURN_UV(_XS_next_prime(val));
    } else {
      _vcallsub((ix == 0) ?  "Math::Prime::Util::_generic_next_prime" :
                             "Math::Prime::Util::_generic_prev_prime" );
      XSRETURN(1);
    }

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
      UV facs[MPU_MAX_FACTORS+1];  /* maximum number of factors is log2n */
      UV nfacs, totient, lastf;
      UV n = lo;
      if (n <= 1) XSRETURN_UV(n);
      nfacs = factor(n, facs);
      totient = 1;
      lastf = 0;
      for (i = 0; i < nfacs; i++) {
        UV f = facs[i];
        if (f == lastf) { totient *= f;               }
        else            { totient *= f-1;  lastf = f; }
      }
      PUSHs(sv_2mortal(newSVuv(totient)));
    }

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
      UV nfactors, lastf;
      UV n = lo;

      if (n <= 1)
        XSRETURN_IV(n);
      if ( (n >= 25) && ( !(n%4) || !(n%9) || !(n%25) ) )
        XSRETURN_IV(0);

      nfactors = factor(n, factors);
      lastf = 0;
      for (i = 0; i < nfactors; i++) {
        if (factors[i] == lastf)
          XSRETURN_IV(0);
        lastf = factors[i];
      }
      XSRETURN_IV( (nfactors % 2) ? -1 : 1 );
    }

IV
_XS_mertens(IN UV n)

UV
_XS_exp_mangoldt(IN UV n)
  CODE:
    if (n <= 1)
      RETVAL = 1;
    else if ((n & (n-1)) == 0)  /* Power of 2 */
      RETVAL = 2;
    else if ((n & 1) == 0)      /* Even number */
      RETVAL = 1;
    /* else if (_XS_is_prime(n))  RETVAL = n; */
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
      RETVAL = factors[0];
    }
  OUTPUT:
    RETVAL

int
_validate_num(SV* n, ...)
  CODE:
    RETVAL = 0;
    if (_validate_int(n, 0)) {
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
    SV* svarg;
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
    if (beg <= end) {
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
#ifdef PERL_HAS_BAD_MULTICALL_REFCOUNT
          if (CvDEPTH(multicall_cv) > 1)
            SvREFCNT_inc(multicall_cv);
#endif
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
