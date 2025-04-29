#ifndef MPU_PRIMALITY_H
#define MPU_PRIMALITY_H

#include "ptypes.h"

extern bool is_pseudoprime(UV const n, UV a);
extern bool is_euler_pseudoprime(UV const n, UV a);
extern bool is_euler_plumb_pseudoprime(UV const n);
extern bool miller_rabin(UV const n, const UV *bases, int nbases);
extern bool is_lucas_pseudoprime(UV n, int strength);
extern bool is_almost_extra_strong_lucas_pseudoprime(UV n, UV increment);
extern bool is_frobenius_pseudoprime(UV n, IV P, IV Q);
extern bool is_frobenius_underwood_pseudoprime(UV n);
extern bool is_frobenius_khashin_pseudoprime(UV n);
extern bool is_perrin_pseudoprime(UV n, uint32_t restricted);

/* 2^p-1:  0 composite, 1 prime, -1 don't know */
extern int is_mersenne_prime(UV p);
extern bool lucas_lehmer(UV p);

#if defined(FUNC_is_strong_pseudoprime)
static bool is_strong_pseudoprime(UV n, UV base) {
  return miller_rabin(n, &base, 1);
}
#endif

extern bool BPSW(UV const n);
extern bool MR32(uint32_t n);

/* General purpose primality test.  Does small-prime divisibility. */
extern bool is_prob_prime(UV n);

/* General purpose primality test without small divisibility tests. */
#if BITS_PER_WORD == 32
#define is_def_prime(n)  MR32(n)
#else
#define is_def_prime(n)  ((n <= 4294967295U) ? MR32(n) : BPSW(n))
#endif

#endif
