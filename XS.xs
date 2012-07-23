
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
/* We're not using anything for which we need ppport.h */
#include "ptypes.h"
#include "cache.h"
#include "sieve.h"
#include "util.h"
#include "factor.h"

MODULE = Math::Prime::Util	PACKAGE = Math::Prime::Util

PROTOTYPES: ENABLE


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

int
_XS_is_prime(IN UV n)

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
        prime_precalc( sqrt(endp) + 0.1 + 1 );
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


void
_XS_factor(IN UV n)
  PPCODE:
    if (n < 4) {
      XPUSHs(sv_2mortal(newSVuv( n ))); /* If n is 0-3, we're done. */
    } else {
      int const verbose = 0;
      UV tlim = 101;  /* Below this we've checked with trial division */
      UV tofac_stack[MPU_MAX_FACTORS+1];
      UV factored_stack[MPU_MAX_FACTORS+1];
      int ntofac = 0;
      int nfactored = 0;
      /* Quick trial divisions.  Crude use of GCD to hopefully go faster. */
      while ( (n% 2) == 0 ) {  n /=  2;  XPUSHs(sv_2mortal(newSVuv(  2 ))); }
      if ( (n >= UVCONST(3*3)) && (gcd_ui(n, UVCONST(3234846615) != 1)) ) {
        while ( (n% 3) == 0 ) {  n /=  3;  XPUSHs(sv_2mortal(newSVuv(  3 ))); }
        while ( (n% 5) == 0 ) {  n /=  5;  XPUSHs(sv_2mortal(newSVuv(  5 ))); }
        while ( (n% 7) == 0 ) {  n /=  7;  XPUSHs(sv_2mortal(newSVuv(  7 ))); }
        while ( (n%11) == 0 ) {  n /= 11;  XPUSHs(sv_2mortal(newSVuv( 11 ))); }
        while ( (n%13) == 0 ) {  n /= 13;  XPUSHs(sv_2mortal(newSVuv( 13 ))); }
        while ( (n%17) == 0 ) {  n /= 17;  XPUSHs(sv_2mortal(newSVuv( 17 ))); }
        while ( (n%19) == 0 ) {  n /= 19;  XPUSHs(sv_2mortal(newSVuv( 19 ))); }
        while ( (n%23) == 0 ) {  n /= 23;  XPUSHs(sv_2mortal(newSVuv( 23 ))); }
        while ( (n%29) == 0 ) {  n /= 29;  XPUSHs(sv_2mortal(newSVuv( 29 ))); }
      }
      if ( (n >= UVCONST(31*31)) && (gcd_ui(n, UVCONST(95041567) != 1)) ) {
        while ( (n%31) == 0 ) {  n /= 31;  XPUSHs(sv_2mortal(newSVuv( 31 ))); }
        while ( (n%37) == 0 ) {  n /= 37;  XPUSHs(sv_2mortal(newSVuv( 37 ))); }
        while ( (n%41) == 0 ) {  n /= 41;  XPUSHs(sv_2mortal(newSVuv( 41 ))); }
        while ( (n%43) == 0 ) {  n /= 43;  XPUSHs(sv_2mortal(newSVuv( 43 ))); }
        while ( (n%47) == 0 ) {  n /= 47;  XPUSHs(sv_2mortal(newSVuv( 47 ))); }
      }
      if ( (n >= UVCONST(53*53)) && (gcd_ui(n, UVCONST(907383479) != 1)) ) {
        while ( (n%53) == 0 ) {  n /= 53;  XPUSHs(sv_2mortal(newSVuv( 53 ))); }
        while ( (n%59) == 0 ) {  n /= 59;  XPUSHs(sv_2mortal(newSVuv( 59 ))); }
        while ( (n%61) == 0 ) {  n /= 61;  XPUSHs(sv_2mortal(newSVuv( 61 ))); }
        while ( (n%67) == 0 ) {  n /= 67;  XPUSHs(sv_2mortal(newSVuv( 67 ))); }
        while ( (n%71) == 0 ) {  n /= 71;  XPUSHs(sv_2mortal(newSVuv( 71 ))); }
      }
      if ( (n >= UVCONST(73*73)) && (gcd_ui(n, UVCONST(4132280413) != 1)) ) {
        while ( (n%73) == 0 ) {  n /= 73;  XPUSHs(sv_2mortal(newSVuv( 73 ))); }
        while ( (n%79) == 0 ) {  n /= 79;  XPUSHs(sv_2mortal(newSVuv( 79 ))); }
        while ( (n%83) == 0 ) {  n /= 83;  XPUSHs(sv_2mortal(newSVuv( 83 ))); }
        while ( (n%89) == 0 ) {  n /= 89;  XPUSHs(sv_2mortal(newSVuv( 89 ))); }
        while ( (n%97) == 0 ) {  n /= 97;  XPUSHs(sv_2mortal(newSVuv( 97 ))); }
      }

      do { /* loop over each remaining factor */
        /* In theory we can try to minimize work using is_definitely_prime(n)
         * but in practice it seems slower. */
        while ( (n >= (tlim*tlim)) && (!_XS_is_prime(n)) ) {
          int split_success = 0;
          /* How many rounds of SQUFOF to try depends on the number size */
          UV sq_rounds = ((n>>29) ==     0) ? 100000 :
                         ((n>>29) < 100000) ? 250000 :
                                              600000;

          /* Small factors will be found quite rapidly with this */
          if (!split_success) {
            split_success = pbrent_factor(n, tofac_stack+ntofac, 1500)-1;
            if (verbose) { if (split_success) printf("pbrent 1:  %"UVuf" %"UVuf"\n", tofac_stack[ntofac], tofac_stack[ntofac+1]); else printf("pbrent 0\n"); }
          }

          if (!split_success) {
            /* SQUFOF does very well with what's left after TD and Rho. */
            /* Racing SQUFOF is about 40% faster and has better success, but
             * has input size restrictions and I'm seeing cases where it gets
             * stuck in stage 2.  For now just stick with the old one.  */
            split_success = squfof_factor(n, tofac_stack+ntofac, sq_rounds)-1;
            if (verbose) printf("squfof %d\n", split_success);
          }

          /* Perhaps prho using different parameters will find it */
          if (!split_success) {
            split_success = prho_factor(n, tofac_stack+ntofac, 800)-1;
            if (verbose) printf("prho %d\n", split_success);
          }

          /* Some rounds of HOLF, good for close to perfect squares */
          if (!split_success) {
            split_success = holf_factor(n, tofac_stack+ntofac, 2000)-1;
            if (verbose) printf("holf %d\n", split_success);
          }

          /* A miniscule fraction of numbers make it here */
          if (!split_success) {
            split_success = prho_factor(n, tofac_stack+ntofac, 256*1024)-1;
            if (verbose) printf("long prho %d\n", split_success);
          }

          if (split_success) {
            MPUassert( split_success == 1, "split factor returned more than 2 factors");
            ntofac++; /* Leave one on the to-be-factored stack */
            if ((tofac_stack[ntofac] == n) || (tofac_stack[ntofac] == 1))
              croak("bad factor\n");
            n = tofac_stack[ntofac];  /* Set n to the other one */
          } else {
            /* trial divisions */
            UV f = tlim;
            UV m = tlim % 30;
            UV limit = (UV) (sqrt(n)+0.1);
            if (verbose) printf("doing trial on %"UVuf"\n", n);
            while (f <= limit) {
              if ( (n%f) == 0 ) {
                do {
                  n /= f;
                  factored_stack[nfactored++] = f;
                } while ( (n%f) == 0 );
                limit = (UV) (sqrt(n)+0.1);
              }
              f += wheeladvance30[m];
              m =  nextwheel30[m];
            }
            break;  /* We just factored n via trial division.  Exit loop. */
          }
        }
        /* n is now prime (or 1), so add to already-factored stack */
        if (n != 1)  factored_stack[nfactored++] = n;
        /* Pop the next number off the to-factor stack */
        if (ntofac > 0)  n = tofac_stack[ntofac-1];
      } while (ntofac-- > 0);
      /* Now push all the factored results in sorted order */
      {
        int i, j;
        for (i = 0; i < nfactored; i++) {
          int mini = i;
          for (j = i+1; j < nfactored; j++)
            if (factored_stack[j] < factored_stack[mini])
              mini = j;
          if (mini != i) {
            UV t = factored_stack[mini];
            factored_stack[mini] = factored_stack[i];
            factored_stack[i] = t;
          }
          XPUSHs(sv_2mortal(newSVuv( factored_stack[i] )));
        }
      }
    }

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
rsqufof_factor(IN UV n)
  PPCODE:
    SIMPLE_FACTOR(racing_squfof_factor, n, 0);

void
pbrent_factor(IN UV n, IN UV maxrounds = 4*1024*1024)
  PPCODE:
    SIMPLE_FACTOR(pbrent_factor, n, maxrounds);

void
prho_factor(IN UV n, IN UV maxrounds = 4*1024*1024)
  PPCODE:
    SIMPLE_FACTOR(prho_factor, n, maxrounds);

void
pminus1_factor(IN UV n, IN UV maxrounds = 1*1024*1024)
  PPCODE:
    SIMPLE_FACTOR(pminus1_factor, n, maxrounds);

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
_XS_RiemannR(double x)
