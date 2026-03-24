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

/* Add 1 to (strint_incr) or subtract 1 from (strint_decr) a signed decimal integer
 * string s/len.  Result is written to out; caller must ensure out has at
 * least len+2 bytes.  Returns the length of the result (no NUL written). */
extern STRLEN strint_incr(char* out, const char* s, STRLEN len);
extern STRLEN strint_decr(char* out, const char* s, STRLEN len);

/* Add (strint_add) or subtract (strint_sub) two signed decimal integer strings.
 * Result is written to out; caller must ensure out has at least
 * max(alen,blen)+1 bytes.  Returns the length of the result (no NUL written). */
extern STRLEN strint_add(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen);
extern STRLEN strint_sub(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen);

#endif
