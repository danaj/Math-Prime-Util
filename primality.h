#ifndef MPU_PRIMALITY_H
#define MPU_PRIMALITY_H

#include "ptypes.h"

#if BITS_PER_WORD == 64 && HAVE_STD_U64 && defined(__GNUC__) && defined(__x86_64__)
#define USE_MONT_PRIMALITY 1
#else
#define USE_MONT_PRIMALITY 0
#endif

extern int is_pseudoprime(UV const n, UV a);
extern int is_euler_pseudoprime(UV const n, UV a);
extern int miller_rabin(UV const n, const UV *bases, int nbases);
extern void lucas_seq(UV* U, UV* V, UV* Qk,  UV n, IV P, IV Q, UV k);
extern int is_lucas_pseudoprime(UV n, int strength);
extern int is_almost_extra_strong_lucas_pseudoprime(UV n, UV increment);
extern int is_frobenius_pseudoprime(UV n, IV P, IV Q);
extern int is_frobenius_underwood_pseudoprime(UV n);
extern int is_frobenius_khashin_pseudoprime(UV n);
extern int is_perrin_pseudoprime(UV n, int restricted);
extern int is_mersenne_prime(UV p);
extern int lucas_lehmer(UV p);
extern int lucasu(IV* U, IV P, IV Q, UV k);
extern int lucasv(IV* V, IV P, IV Q, UV k);

extern int BPSW(UV const n);
extern int is_prob_prime(UV n);

#endif
