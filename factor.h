#ifndef MPU_FACTOR_H
#define MPU_FACTOR_H

#include "ptypes.h"

#if BITS_PER_WORD == 64
#define MPU_MAX_FACTORS 64
#define MPU_MAX_DFACTORS 15
#else
#define MPU_MAX_FACTORS 32
#define MPU_MAX_DFACTORS 9
#endif

/* These all return the number of factors set in factors[].
 * Nothing found:     returns 1 and factors[0] = n
 * One factor found:  returns 2 and factors[0] = f, factors[1] = n/f
 * ...
 */

extern int factor(UV n, UV *factors);
extern int factor_exp(UV n, UV *factors, UV* exponents);
extern int factor_one(UV n, UV *factors, bool primality, bool trial);
extern UV  divisor_sum(UV n, UV k);

extern int trial_factor(UV n, UV *factors, UV first, UV last);
extern int fermat_factor(UV n, UV *factors, UV rounds);
extern int holf_factor(UV n, UV *factors, UV rounds);
extern int pbrent_factor(UV n, UV *factors, UV maxrounds, UV a);
extern int prho_factor(UV n, UV *factors, UV maxrounds);
extern int pminus1_factor(UV n, UV *factors, UV B1, UV B2);
extern int pplus1_factor(UV n, UV *factors, UV B);
extern int squfof_factor(UV n, UV *factors, UV rounds);
extern int lehman_factor(UV n, UV *factors, bool dotrial);
extern int cheb_factor(UV n, UV *factors, UV B, UV initx);

extern UV* divisor_list(UV n, UV *num_divisors, UV maxd);

extern int prime_omega(UV n);     /* number of distinct prime factors */
extern int prime_bigomega(UV n);  /* number of prime factors w/ multiplicity */
/* bigomega => with_multiplicity=1        omega => with_multiplicity=0 */
extern unsigned char* range_nfactor_sieve(UV lo, UV hi, bool with_multiplicity);

/* Factoring all numbers in a range. */
typedef struct {
  UV   lo;
  UV   hi;
  UV   n;
  bool is_square_free;
  UV  *factors;
  UV   _coffset;
  UV   _noffset;
  UV  *_farray;
  UV  *_nfactors;
} factor_range_context_t;
extern factor_range_context_t factor_range_init(UV lo, UV hi, bool square_free);
extern int factor_range_next(factor_range_context_t *ctx);
extern void factor_range_destroy(factor_range_context_t *ctx);

/*
extern UV dlp_trial(UV a, UV g, UV p, UV maxrounds);
extern UV dlp_prho(UV a, UV g, UV p, UV n, UV maxrounds);
extern UV dlp_bsgs(UV a, UV g, UV p, UV n, UV maxent);
*/
/* Generic znlog returns k that solves a = g^k mod p */
extern UV znlog(UV a, UV g, UV p);
/* znlog given prime gorder = znorder(g,p) */
extern UV znlog_solve(UV a, UV g, UV p, UV gorder);

#endif
