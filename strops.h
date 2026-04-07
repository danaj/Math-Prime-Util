#ifndef MPU_STROPS_H
#define MPU_STROPS_H

#include "ptypes.h"

/* Compare two signed decimal integer strings.
 * Returns -1, 0, or 1 like strcmp. */
extern int  strint_cmp(const char* a, STRLEN alen, const char* b, STRLEN blen);

/* Compare a/alen against b/blen for min or max.
 * b/blen is validated; a/alen must already be canonical.
 * Returns 0 to select a, 1 to select b. */
extern bool strint_minmax(bool min, const char* a, STRLEN alen, const char* b, STRLEN blen);


/* Add or subtract two signed decimal integer strings (a+b or a-b).
 * negate_b=false for addition, negate_b=true for subtraction.
 * Result is written to out; caller must ensure out has at least
 * max(alen,blen)+2 bytes.  out may alias a; out must not alias b.
 * Returns the length of the result. */
extern STRLEN strint_add_s(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen, bool negate_b);
#define strint_add(out, a, alen, b, blen) strint_add_s(out, a, alen, b, blen, 0)
#define strint_sub(out, a, alen, b, blen) strint_add_s(out, a, alen, b, blen, 1)

/* Multiply two signed decimal integer strings.
 * Result is written to out; caller must ensure out has at least
 * alen+blen bytes.  Returns the length of the result (no NUL written). */
extern STRLEN strint_mul(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen);

/* Return a*b+c or a*b-c to out.
 * out must have at least max(alen+blen,clen)+2 bytes. */
extern STRLEN strint_muladd_s(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen, const char* c, STRLEN clen, bool negate_c);

/* Return (a*b+c) mod m or (a*b-c) mod m.
 * negate_c=false for addition, negate_c=true for subtraction.
 * Result is in [0, m-1].  out must have at least mlen bytes. */
extern STRLEN strint_muladdmod_s(char* out,
                                 const char* a, STRLEN alen,
                                 const char* b, STRLEN blen,
                                 const char* c, STRLEN clen, bool negate_c,
                                 const char* m, STRLEN mlen);
#define strint_muladdmod(out,a,al,b,bl,c,cl,m,ml) strint_muladdmod_s(out,a,al,b,bl,c,cl,0,m,ml)
#define strint_mulsubmod(out,a,al,b,bl,c,cl,m,ml) strint_muladdmod_s(out,a,al,b,bl,c,cl,1,m,ml)

/* Raise a to non-negative integer power exp.
 * out must have at least limit bytes.
 * Returns 0 if result would not fit in limit characters.
 * Returns the result length. */
extern STRLEN strint_pow(char* out, const char* a, STRLEN alen, UV exp, STRLEN limit);

/* floor(log_base(a)) for a >= 1, base >= 2.
 * Returns UV_MAX on error (a < 1 or base < 2). */
extern UV strint_logint(const char* a, STRLEN alen, UV base);

/* floor(a^(1/k)) for a >= 0, k >= 1.
 * out must have at least alen+2 bytes.
 * Returns 0 on error (k==0, a<0). */
extern STRLEN strint_rootint(char* out, const char* a, STRLEN alen, UV k);

/* Signed floor division and remainder.
 * Either qout or rout (and its companion length pointer) may be NULL.
 * qout needs alen+1 bytes; rout needs blen bytes.
 * Remainder carries the sign of b (floor convention).
 * Returns false if b = 0, true otherwise. */
extern bool strint_divmod(char* qout, STRLEN* qlen, char* rout, STRLEN* rlen,
                          const char* a, STRLEN alen,
                          const char* b, STRLEN blen);

/* Three division helpers, all using strint_divmod but more convenient. */

/* out needs alen+1 bytes.  Returns 0 on failure. */
extern STRLEN strint_divint(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen);
/* out needs blen bytes.  Returns 0 on failure. */
extern STRLEN strint_modint(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen);
/* out needs alen+2 bytes.  Returns 0 on failure. */
extern STRLEN strint_cdivint(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen);

/* Return a mod b as a UV.  b > 0.  For negative a, returns |a| mod b. */
extern UV strint_moduv(const char* a, STRLEN alen, UV b);

/* Extract all small prime factors (up to 47 on 64-bit, 23 on 32-bit) from a.
 * Factors found, with repetition, are appended to out_f[0..*nf); *nf is updated.
 * Returns 0 if the remaining cofactor fits in a UV (*uv_out is set).
 * Returns the cofactor string length otherwise (written to str_out, may alias a;
 * str_out needs alen+1 bytes). */
extern STRLEN strint_remove_small_factors(char* str_out, UV* uv_out,
                                          UV* out_f, int* nf,
                                          const char* a, STRLEN alen);

/* Trial-factor a decimal integer string by primes in primes[0..nprimes).
 * Factors found (with repetition) are appended to out_f[0..*nf); *nf is updated.
 * Caller must ensure out_f has room for at least alen*4 entries.
 * Returns 0 if the remaining cofactor fits in a UV (*uv_out is set).
 * Returns the cofactor string length otherwise (written to str_out, may alias a;
 * str_out needs alen+1 bytes). */
extern STRLEN strint_trial_factor(char* str_out, UV* uv_out,
                                  UV* out_f, int* nf,
                                  const char* a, STRLEN alen,
                                  const uint32_t* primes, uint32_t nprimes);

/* Absolute value.  out needs alen bytes. */
extern STRLEN strint_abs(char* out, const char* a, STRLEN alen);
/* Negation.  out needs alen+1 bytes. */
extern STRLEN strint_neg(char* out, const char* a, STRLEN alen);

/* Left shift: a * 2^k.
 * out needs alen + ceil(k * log10(2)) + 2 bytes. */
extern STRLEN strint_lshiftint(char* out, const char* a, STRLEN alen, UV k);
/* Right shift truncating toward zero: sign(a) * floor(|a| / 2^k).
 * out needs alen bytes. */
extern STRLEN strint_rshiftint(char* out, const char* a, STRLEN alen, UV k);
/* Arithmetic right shift (floor division): floor(a / 2^k).
 * out needs alen bytes. */
extern STRLEN strint_rashiftint(char* out, const char* a, STRLEN alen, UV k);

#endif
