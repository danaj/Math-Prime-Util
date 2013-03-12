
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
/* We're not using anything for which we need ppport.h */
#ifndef XSRETURN_UV   /* Er, almost.  Fix 21086 from Sep 2003 */
  #define XST_mUV(i,v)  (ST(i) = sv_2mortal(newSVuv(v))  )
  #define XSRETURN_UV(v) STMT_START { XST_mUV(0,v);  XSRETURN(1); } STMT_END
#endif
#include "ptypes.h"
#include "cache.h"
#include "sieve.h"
#include "util.h"
#include "factor.h"
#include "lehmer.h"
#include "aks.h"

MODULE = Math::Prime::Util	PACKAGE = Math::Prime::Util

PROTOTYPES: ENABLE


void
_XS_set_verbose(IN int verbose)

int
_XS_get_verbose()


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
_XS_legendre_pi(IN UV n)

UV
_XS_meissel_pi(IN UV n)

UV
_XS_lehmer_pi(IN UV n)

UV
_XS_LMO_pi(IN UV n)

UV
_XS_nth_prime(IN UV n)

int
_XS_is_prime(IN UV n)

int
_XS_is_aks_prime(IN UV n)

UV
_XS_next_prime(IN UV n)

UV
_XS_prev_prime(IN UV n)


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
    const unsigned char* sieve;
    AV* av = newAV();
  CODE:
    if (low <= high) {
      if (get_prime_cache(high, &sieve) < high) {
        release_prime_cache(sieve);
        croak("Could not generate sieve for %"UVuf, high);
      } else {
        if ((low <= 2) && (high >= 2)) { av_push(av, newSVuv( 2 )); }
        if ((low <= 3) && (high >= 3)) { av_push(av, newSVuv( 3 )); }
        if ((low <= 5) && (high >= 5)) { av_push(av, newSVuv( 5 )); }
        if (low < 7) { low = 7; }
        START_DO_FOR_EACH_SIEVE_PRIME( sieve, low, high ) {
           av_push(av,newSVuv(p));
        } END_DO_FOR_EACH_SIEVE_PRIME
        release_prime_cache(sieve);
      }
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


#define SIMPLE_FACTOR(func, n, rounds) \
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
        nfactors = func(n, factors, rounds); \
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
    if (n <= 1) {
      XPUSHs(sv_2mortal(newSVuv( n )));
    } else {
      while ( (n% 2) == 0 ) {  n /=  2;  XPUSHs(sv_2mortal(newSVuv( 2 ))); }
      while ( (n% 3) == 0 ) {  n /=  3;  XPUSHs(sv_2mortal(newSVuv( 3 ))); }
      while ( (n% 5) == 0 ) {  n /=  5;  XPUSHs(sv_2mortal(newSVuv( 5 ))); }
      if (n == 1) {  /* done */ }
      else if (_XS_is_prime(n)) { XPUSHs(sv_2mortal(newSVuv( n ))); }
      else {
        UV factors[MPU_MAX_FACTORS+1];
        int i, nfactors;
        nfactors = pminus1_factor(n, factors, B1, B2);
        for (i = 0; i < nfactors; i++) {
          XPUSHs(sv_2mortal(newSVuv( factors[i] )));
        }
      }
    }

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

int
_XS_is_prob_prime(IN UV n)

double
_XS_ExponentialIntegral(double x)

double
_XS_LogarithmicIntegral(double x)

double
_XS_RiemannZeta(double x)
  CODE:
    RETVAL = (double) ld_riemann_zeta(x);
  OUTPUT:
    RETVAL

double
_XS_RiemannR(double x)


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
      char* mu = _moebius_range(lo, hi);
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

double
_XS_chebyshev_theta(IN UV n)

double
_XS_chebyshev_psi(IN UV n)

UV
_XS_divisor_sum(IN UV n)
