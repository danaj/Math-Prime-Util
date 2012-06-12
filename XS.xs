
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
/* We're not using anything for which we need ppport.h */
#include "ptypes.h"
#include "sieve.h"
#include "util.h"
#include "bitarray.h"
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
prime_count(IN UV low, IN UV high = 0)
  CODE:
    if (high == 0) {   /* Without a Perl layer in front of this, we'll have */
      high = low;      /* the pathological case of a-0 turning into 0-a.    */
      low = 0;
    }
    if (GIMME_V == G_VOID) {
      prime_precalc(high);
      RETVAL = 0;
    } else {
      RETVAL = prime_count(low, high);
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
  PPCODE:
    if (n < 4) {
      XPUSHs(sv_2mortal(newSVuv( n ))); /* If n is 0-3, we're done. */
    } else {
      UV tlim = 19;  /* Below this we've checked */
      UV factor_stack[MPU_MAX_FACTORS+1];
      int nstack = 0;
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
        tlim = 31;
      }
      if ( (n >= UVCONST(31*31)) && (gcd_ui(n, UVCONST(95041567) != 1)) ) {
        while ( (n%31) == 0 ) {  n /= 31;  XPUSHs(sv_2mortal(newSVuv( 31 ))); }
        while ( (n%37) == 0 ) {  n /= 37;  XPUSHs(sv_2mortal(newSVuv( 37 ))); }
        while ( (n%41) == 0 ) {  n /= 41;  XPUSHs(sv_2mortal(newSVuv( 41 ))); }
        while ( (n%43) == 0 ) {  n /= 43;  XPUSHs(sv_2mortal(newSVuv( 43 ))); }
        while ( (n%47) == 0 ) {  n /= 47;  XPUSHs(sv_2mortal(newSVuv( 47 ))); }
        tlim = 53;
      }
      do { /* loop over each remaining factor */
        while ( (n >= (tlim*tlim)) && (!is_definitely_prime(n)) ) {
          int split_success = 0;
          if (n > UVCONST(10000000) ) {  /* tune this */
            /* For sufficiently large n, try more complex methods. */
            /* SQUFOF (succeeds 98-99.9%) */
            split_success = squfof_factor(n, factor_stack+nstack, 256*1024)-1;
            /* A few rounds of Pollard rho (succeeds most of the rest) */
            if (!split_success) {
              split_success = prho_factor(n, factor_stack+nstack, 400)-1;
            }
            /* Some rounds of HOLF, good for close to perfect squares */
            if (!split_success) {
              split_success = holf_factor(n, factor_stack+nstack, 2000)-1;
            }
          }
          if (split_success) {
            MPUassert( split_success == 1, "split factor returned more than 2 factors");
            nstack++;
            n = factor_stack[nstack];
          } else {
            /* trial divisions */
            UV f = tlim;
            UV m = tlim % 30;
            UV limit = (UV) (sqrt(n)+0.1);
            while (f <= limit) {
              if ( (n%f) == 0 ) {
                do {
                  n /= f;  XPUSHs(sv_2mortal(newSVuv( f )));
                } while ( (n%f) == 0 );
                limit = (UV) (sqrt(n)+0.1);
              }
              f += wheeladvance30[m];
              m =  nextwheel30[m];
            }
            break;  /* We just factored n via trial division.  Exit loop. */
          }
        }
        if (n != 1)  XPUSHs(sv_2mortal(newSVuv( n )));
        if (nstack > 0)  n = factor_stack[nstack-1];
      } while (nstack-- > 0);
    }

#define SIMPLE_FACTOR(func, n, rounds) \
    if (n <= 1) { \
      XPUSHs(sv_2mortal(newSVuv( n ))); \
    } else { \
      while ( (n% 2) == 0 ) {  n /=  2;  XPUSHs(sv_2mortal(newSVuv( 2 ))); } \
      while ( (n% 3) == 0 ) {  n /=  3;  XPUSHs(sv_2mortal(newSVuv( 3 ))); } \
      while ( (n% 5) == 0 ) {  n /=  5;  XPUSHs(sv_2mortal(newSVuv( 5 ))); } \
      if (n == 1) {  /* done */ } \
      else if (is_prime(n)) { XPUSHs(sv_2mortal(newSVuv( n ))); } \
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
trial_factor(IN UV n)
  PPCODE:
    SIMPLE_FACTOR(trial_factor, n, UV_MAX);

void
fermat_factor(IN UV n, IN UV maxrounds = 64*1024*1024)
  PPCODE:
    SIMPLE_FACTOR(fermat_factor, n, maxrounds);

void
holf_factor(IN UV n, IN UV maxrounds = 64*1024*1024)
  PPCODE:
    SIMPLE_FACTOR(holf_factor, n, maxrounds);

void
squfof_factor(IN UV n, IN UV maxrounds = 4*1024*1024)
  PPCODE:
    SIMPLE_FACTOR(squfof_factor, n, maxrounds);

void
pbrent_factor(IN UV n, IN UV maxrounds = 4*1024*1024)
  PPCODE:
    SIMPLE_FACTOR(pbrent_factor, n, maxrounds);

void
prho_factor(IN UV n, IN UV maxrounds = 4*1024*1024)
  PPCODE:
    SIMPLE_FACTOR(prho_factor, n, maxrounds);

void
pminus1_factor(IN UV n, IN UV maxrounds = 4*1024*1024)
  PPCODE:
    SIMPLE_FACTOR(pminus1_factor, n, maxrounds);

int
miller_rabin(IN UV n, ...)
  PREINIT:
    UV bases[64];
    int prob_prime = 1;
    int c = 1;
  CODE:
    if (items < 2)
      croak("No bases given to miller_rabin");
    if ( (n == 0) || (n == 1) ) XSRETURN_IV(0);   /* 0 and 1 are composite */
    if ( (n == 2) || (n == 3) ) XSRETURN_IV(2);   /* 2 and 3 are prime */
    while (c < items) {
      int b = 0;
      while (c < items) {
        bases[b++] = SvUV(ST(c));
        c++;
        if (b == 64) break;
      }
      prob_prime = miller_rabin(n, bases, b);
      if (prob_prime != 1)
        break;
    }
    RETVAL = prob_prime;
  OUTPUT:
    RETVAL

int
is_prob_prime(IN UV n)

double
ExponentialIntegral(double x)

double
LogarithmicIntegral(double x)

double
RiemannR(double x)
