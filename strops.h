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

/* Add 1 to or subtract 1 from a signed decimal integer string s(len).
 * Result is written to out; caller must ensure out has at
 * least len+1 bytes.  Returns the length of the result. */
extern STRLEN strint_incr(char* out, const char* s, STRLEN len);
extern STRLEN strint_decr(char* out, const char* s, STRLEN len);

/* Add or subtract two signed decimal integer strings (a+b or a-b).
 * negate_b=false for addition, negate_b=true for subtraction.
 * Result is written to out; caller must ensure out has at least
 * max(alen,blen)+2 bytes.  Returns the length of the result. */
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

#endif
