#ifndef MPU_STROPS_H
#define MPU_STROPS_H

#include "ptypes.h"

/* Compare two signed decimal integer strings.
 * Returns -1, 0, or 1 like strcmp. */
extern int  strnum_cmp(const char* a, STRLEN alen, const char* b, STRLEN blen);

/* Compare a/alen against b/blen for min or max.
 * b/blen is validated; a/alen must already be canonical.
 * Returns 0 to select a, 1 to select b. */
extern bool strnum_minmax(bool min, const char* a, STRLEN alen, const char* b, STRLEN blen);

/* Add 1 to (strincr) or subtract 1 from (strdecr) a signed decimal integer
 * string s/len.  Result is written to out; caller must ensure out has at
 * least len+2 bytes.  Returns the length of the result (no NUL written). */
extern STRLEN strincr(char* out, const char* s, STRLEN len);
extern STRLEN strdecr(char* out, const char* s, STRLEN len);

#endif
