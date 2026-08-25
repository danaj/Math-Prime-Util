#ifndef MPU_UTIL_H
#define MPU_UTIL_H

#include "ptypes.h"

extern int  _XS_get_verbose(void);
extern void _XS_set_verbose(int v);
extern int  _XS_get_callgmp(void);
extern void _XS_set_callgmp(int v);
extern int  _XS_get_nobigint(void);
extern void _XS_set_nobigint(int v);
/* Disable all manual seeding */
extern bool  _XS_get_secure(void);
extern void _XS_set_secure(void);

#if HAVE_UINT128
extern void* mpu_aligned_alloc(UV count, Size_t size, Size_t alignment);
extern void  mpu_aligned_free(void* ptr);
#endif

extern bool is_prime(UV x);
extern UV   next_prime(UV x);
extern UV   prev_prime(UV x);

/* Simple estimate for upper limit:  max_nprimes(n) >= prime_count(n) */
extern UV   max_nprimes(UV n) ISCONSTFUNC;

/* If failure is returned, errno should be set. */
extern bool print_primes(UV low, UV high, int fd);

/* Returns maximal k for c^k = n for k > 1, n > 1.  0 otherwise. */
extern uint32_t powerof_ret(UV n, uint32_t *root);
#define powerof(n) powerof_ret(n,0)

/* Return true if n = r^k for the given k, sets root if given */
extern bool is_power_ret(UV n, UV k, uint32_t *root);
#define is_power(n,k) is_power_ret(n,k,0)

extern uint32_t icbrt(UV n) ISCONSTFUNC;
extern UV rootint(UV n, UV k) ISCONSTFUNC;
extern UV crootint(UV n, UV k) ISCONSTFUNC;
extern UV ipowsafe(UV n, UV k) ISCONSTFUNC;  /* returns UV_MAX if overflows */
extern UV lcmsafe(UV x, UV u) ISCONSTFUNC;   /* returns 0 if overflows */
extern UV valuation(UV n, UV k) ISCONSTFUNC;
extern UV valuation_remainder(UV n, UV k, UV *r);
extern UV logint(UV n, UV b) ISCONSTFUNC;

extern UV powersum(UV n, UV k) ISCONSTFUNC;

extern int kronecker_uu(UV a, UV b) ISCONSTFUNC;
extern int kronecker_su(IV a, UV b) ISCONSTFUNC;
extern int kronecker_ss(IV a, IV b) ISCONSTFUNC;

extern UV pn_primorial(UV n) ISCONSTFUNC;
extern UV primorial(UV n) ISCONSTFUNC;
extern UV factorial(UV n) ISCONSTFUNC;
extern UV subfactorial(UV n) ISCONSTFUNC;
extern UV bell_number(UV n) ISCONSTFUNC;
extern UV fubini(UV n) ISCONSTFUNC;
extern UV binomial(UV n, UV k) ISCONSTFUNC;
extern UV catalan_number(UV n) ISCONSTFUNC;
extern UV multifactorial(UV n, UV k) ISCONSTFUNC;
extern UV falling_factorial(UV n, UV m) ISCONSTFUNC;
extern UV rising_factorial(UV n, UV m) ISCONSTFUNC;
extern IV falling_factorial_s(IV n, UV m) ISCONSTFUNC;
extern IV rising_factorial_s(IV n, UV m) ISCONSTFUNC;
extern IV gcdext(IV a, IV b, IV* u, IV* v, IV* s, IV* t); /* Ext Euclidean */
extern UV modinverse(UV a, UV p) ISCONSTFUNC;  /* Returns 1/a mod p */
extern UV divmod(UV a, UV b, UV n) ISCONSTFUNC;/* Returns a/b mod n */
extern UV gcddivmod(UV a, UV b, UV n) ISCONSTFUNC; /* divmod(a/gcd,b/gcd,n) */

/* UV_MAX overflow */
extern UV pisano_period(UV n);
extern UV floor_sum(UV n, UV m, UV a, UV b);

/* 0 overflow, -1 no inverse, 1 ok */
/* The a/n arrays will be sorted by descending n. */
extern int chinese(UV *r, UV *lcm, UV* a, UV* n, UV num);/* Chinese Remainder */

/* Do the inverse for a negative modular power / root. a^-k => (1/a)^k mod n */
extern bool prep_pow_inv(UV *a, UV *k, int kstatus, UV n);

/* Signed division and remainder.  Returns remainder.*/
extern IV tdivrem(IV *q, IV *r, IV D, IV d);   /* divrem trunc */
extern IV fdivrem(IV *q, IV *r, IV D, IV d);   /* divrem floor */
extern IV cdivrem(IV *q, IV *r, IV D, IV d);   /* divrem ceiling */
extern IV edivrem(IV *q, IV *r, IV D, IV d);   /* divrem Euclidian */
extern UV ivmod(IV a, UV n) ISCONSTFUNC;       /* Returns a mod n (trunc) */

extern bool muladd_uv_signmag(int *sign, UV *hi, UV *lo,
                              UV a, UV b, UV c,
                              int asign, int bsign, int csign);

extern UV   exp_mangoldt(UV n);
extern UV   znprimroot(UV n);
extern UV   znorder(UV a, UV n);
/* nprime says to assume n = p or n = 2p.  Skips power and primality tests. */
extern bool is_primitive_root(UV a, UV n, bool nprime);

extern bool is_square_free(UV n);
extern bool is_perfect_number(UV n) ISCONSTFUNC;
extern bool is_fundamental(UV n, bool neg);
extern bool is_semiprime(UV n);
extern bool is_almost_prime(UV k, UV n);
extern bool is_cyclic(UV n);
extern bool is_carmichael(UV n);
extern UV   is_quasi_carmichael(UV n);  /* Returns number of bases */
extern UV   pillai_v(UV n) ISCONSTFUNC; /* v: v! % n == n-1 && n % v != 1 */
extern UV   qnr(UV n);
extern bool is_qr(UV a, UV n);         /* kronecker that works for composites */
extern bool is_practical(UV n);
extern int  is_delicate_prime(UV n, uint32_t b);
/* Returns -1 if an intermediate digit sum overflows. */
extern int  happy_height(UV n, uint32_t base, uint32_t exponent) ISCONSTFUNC;

extern bool is_smooth(UV n, UV k);
extern bool is_rough(UV n, UV k);

extern bool is_sum_of_two_squares(UV n);
extern bool is_sum_of_three_squares(UV n) ISCONSTFUNC;
extern bool cornacchia(UV *x, UV *y, UV d, UV p);

extern UV debruijn_psi(UV x, UV y);
extern UV buchstab_phi(UV x, UV y);

extern UV stirling3(UV n, UV m) ISCONSTFUNC;
extern IV stirling2(UV n, UV m) ISCONSTFUNC;
extern IV stirling1(UV n, UV m) ISCONSTFUNC;

extern bool bernfrac(IV *num, UV *den, UV n);
extern bool harmfrac(UV *num, UV *den, UV n);

extern UV hclassno(UV n);
extern IV ramanujan_tau(UV n) ISCONSTFUNC;

extern char* pidigits(uint32_t digits);

/* Returns 1 if good, 0 if bad, -1 if non canon, 2 ok but out of range */
extern int validate_zeckendorf(const char* str, size_t len);
extern UV  from_zeckendorf(const char* str, size_t len);
extern char* to_zeckendorf(UV n);

extern bool is_catalan_pseudoprime(UV n);
extern bool is_higgs_prime(UV n);

extern UV  polygonal_root(UV n, UV k, bool* overflow);

extern UV  npartitions(UV n);
extern UV  npartitionsq(UV n);
extern UV  consecutive_integer_lcm(UV n);

extern UV  frobenius_number(UV* A, uint32_t alen);

extern bool num_to_perm(UV rank, int n, int *vec);
extern bool perm_to_num(int n, int *vec, UV *rank);
extern void randperm(void* ctx, UV n, UV k, UV *S);

extern UV random_factored_integer(void* ctx, UV n, int *nf, UV *factors);

extern UV gcdz(UV x, UV y) ISCONSTFUNC;

#endif
