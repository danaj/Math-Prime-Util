
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
/* We're not using anything for which we need ppport.h */
#include "sieve.h"
#include "util.h"
#include "bitarray.h"
#include "factor.h"

MODULE = Math::Prime::Util	PACKAGE = Math::Prime::Util

PROTOTYPES: ENABLE


void
prime_precalc(IN UV n)

void
prime_free()

UV
prime_count(IN UV n)
  CODE:
    if (GIMME_V == G_VOID) {
      prime_precalc(n);
      RETVAL = 0;
    } else {
      RETVAL = prime_count(n);
    }
  OUTPUT:
    RETVAL

UV
prime_count_lower(IN UV n)

UV
prime_count_upper(IN UV n)

UV
prime_count_approx(IN UV n)

UV
nth_prime(IN UV n)

UV
nth_prime_lower(IN UV n)

UV
nth_prime_upper(IN UV n)

UV
nth_prime_approx(IN UV n)

int
is_prime(IN UV n)

UV
next_prime(IN UV n)

UV
prev_prime(IN UV n)


UV
_get_prime_cache_size()
  CODE:
    RETVAL = get_prime_cache_size();
  OUTPUT:
    RETVAL

int
_maxbits()
  CODE:
    RETVAL = BITS_PER_WORD;
  OUTPUT:
    RETVAL


SV*
sieve_primes(IN UV low, IN UV high)
  PREINIT:
    UV  s;
    const unsigned char* sieve;
    AV* av = newAV();
  CODE:
    if (low <= high) {
      if (get_prime_cache(high, &sieve) < high) {
        croak("Could not generate sieve for %"UVuf, high);
      } else {
        if ((low <= 2) && (high >= 2)) { av_push(av, newSVuv( 2 )); }
        if ((low <= 3) && (high >= 3)) { av_push(av, newSVuv( 3 )); }
        if ((low <= 5) && (high >= 5)) { av_push(av, newSVuv( 5 )); }
        if (low < 7) { low = 7; }
        START_DO_FOR_EACH_SIEVE_PRIME( sieve, low, high ) {
           av_push(av,newSVuv(p));
        } END_DO_FOR_EACH_SIEVE_PRIME
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
      curprime = next_trial_prime(low);
      while (curprime <= high) {
        av_push(av,newSVuv(curprime));
        curprime = next_trial_prime(curprime);
      }
    }
    RETVAL = newRV_noinc( (SV*) av );
  OUTPUT:
    RETVAL

SV*
segment_primes(IN UV low, IN UV high, IN UV segment_size = 65536UL)
  PREINIT:
    AV* av = newAV();
    unsigned char* sieve;
    UV low_d, high_d;
  CODE:
    if ((low <= 2) && (high >= 2)) { av_push(av, newSVuv( 2 )); }
    if ((low <= 3) && (high >= 3)) { av_push(av, newSVuv( 3 )); }
    if ((low <= 5) && (high >= 5)) { av_push(av, newSVuv( 5 )); }
    if (low < 7)  low = 7;
    if (low <= high) {
      /* Call the segment siever one or more times */
      sieve = (unsigned char*) malloc( segment_size );
      if (sieve == 0)
        croak("Could not allocate %"UVuf" bytes for segment sieve", segment_size);

      /* To protect vs. overflow, work entirely with d. */
      low_d  = low  / 30;
      high_d = high / 30;
      while ( low_d <= high_d ) {
        UV seghigh_d = ((high_d - low_d) < segment_size)
                       ? high_d
                       : (low_d + segment_size-1);
        UV range_d = seghigh_d - low_d + 1;
        UV seghigh = (seghigh_d == high_d) ? high : (seghigh_d*30+29);
        UV segbase = low_d * 30;
        /* printf("  startd = %"UVuf"  endd = %"UVuf"\n", startd, endd); */
  
        assert( seghigh_d >= low_d );
        assert( range_d <= segment_size );

        /* Sieve from startd*30+1 to endd*30+29.  */
        if (sieve_segment(sieve, low_d, seghigh_d) == 0) {
          croak("Could not segment sieve from %"UVuf" to %"UVuf, segbase+1, seghigh);
          break;
        }

        START_DO_FOR_EACH_SIEVE_PRIME( sieve, low - segbase, seghigh - segbase )
          av_push(av,newSVuv( segbase + p ));
        END_DO_FOR_EACH_SIEVE_PRIME

        low_d += range_d;
        low = seghigh+2;
      }
      free(sieve);
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
        free(sieve);
      }
    }
    RETVAL = newRV_noinc( (SV*) av );
  OUTPUT:
    RETVAL


SV*
erat_simple_primes(IN UV low, IN UV high)
  PREINIT:
    UV* sieve;
    UV s;
    AV* av = newAV();
  CODE:
    if (low <= high) {
      sieve = sieve_erat(high);
      if (sieve == 0) {
        croak("Could not generate sieve for %"UVuf, high);
      } else {
        if (low <= 2) { av_push(av, newSVuv( 2 )); low = 3; }
        low  = low/2;
        high = (high-1)/2;
        for (s = low; s <= high; s++) {
          if ( ! IS_SET_ARRAY_BIT(sieve, s) ) {
            av_push(av,newSVuv( 2*s+1 ));
          }
        }
        free(sieve);
      }
    }
    RETVAL = newRV_noinc( (SV*) av );
  OUTPUT:
    RETVAL

void
factor(IN UV n)
  PREINIT:
    UV const maxtrials = UV_MAX;
    UV factors[MPU_MAX_FACTORS+1];
    int nfactors;
    int i;
  PPCODE:
#if BITS_PER_WORD == 32
    nfactors = trial_factor(n, factors, maxtrials);
    if (nfactors < 1)
      croak("No factors from trial_factor");
    for (i = 0; i < nfactors; i++) {
      XPUSHs(sv_2mortal(newSVuv( factors[i] )));
    }
#else
    /* TODO: worry about squfof overflow */
    if ( n < UVCONST(0xFFFFFFFF) ) {
      /* small number */
      nfactors = trial_factor(n, factors, maxtrials);
    } else {
      UV sqfactors, f1, f2;
      /* First factor out small numbers */
      nfactors = trial_factor(n, factors, 29);
      /* Use SQUFOF to separate the remainder */
      n = factors[--nfactors];
      sqfactors = squfof_factor(n, factors+nfactors, 800000);
      assert(sqfactors <= 2);
      if (sqfactors == 1) {
        n = factors[nfactors];
        nfactors += trial_factor(n, factors+nfactors, maxtrials);
      } else {
        UV n1 = factors[nfactors];
        n = factors[nfactors+1];
        nfactors += trial_factor(n1, factors+nfactors, maxtrials);
        nfactors += trial_factor(n, factors+nfactors, maxtrials);
      }
    }
    if (nfactors < 1) croak("No factors");
    for (i = 0; i < nfactors; i++) {
      XPUSHs(sv_2mortal(newSVuv( factors[i] )));
    }
#endif

void
fermat_factor(IN UV n)
  PREINIT:
    UV factors[MPU_MAX_FACTORS+1];
    int nfactors;
    int i;
  PPCODE:
    nfactors = fermat_factor(n, factors);
    if (nfactors < 1)
      croak("No factors from fermat_factor");
    for (i = 0; i < nfactors; i++) {
      XPUSHs(sv_2mortal(newSVuv( factors[i] )));
    }

void
squfof_factor(IN UV n)
  PREINIT:
    UV factors[MPU_MAX_FACTORS+1];
    int nfactors;
    int i;
  PPCODE:
    nfactors = squfof_factor(n, factors, 800000);
    if (nfactors < 1)
      croak("No factors from squfof_factor");
    for (i = 0; i < nfactors; i++) {
      XPUSHs(sv_2mortal(newSVuv( factors[i] )));
    }


void
pbrent_factor(IN UV n)
  PREINIT:
    UV factors[MPU_MAX_FACTORS+1];
    int nfactors;
    int i;
  PPCODE:
    nfactors = pbrent_factor(n, factors, 50000000);
    if (nfactors < 1)
      croak("No factors from pbrent_factor");
    for (i = 0; i < nfactors; i++) {
      XPUSHs(sv_2mortal(newSVuv( factors[i] )));
    }

void
prho_factor(IN UV n)
  PREINIT:
    UV factors[MPU_MAX_FACTORS+1];
    int nfactors;
    int i;
  PPCODE:
    nfactors = prho_factor(n, factors, 50000000);
    if (nfactors < 1)
      croak("No factors from prho_factor");
    for (i = 0; i < nfactors; i++) {
      XPUSHs(sv_2mortal(newSVuv( factors[i] )));
    }

void
pminus1_factor(IN UV n)
  PREINIT:
    UV factors[MPU_MAX_FACTORS+1];
    int nfactors;
    int i;
  PPCODE:
    nfactors = pminus1_factor(n, factors, 50000000);
    if (nfactors < 1)
      croak("No factors from pminus1_factor");
    for (i = 0; i < nfactors; i++) {
      XPUSHs(sv_2mortal(newSVuv( factors[i] )));
    }
