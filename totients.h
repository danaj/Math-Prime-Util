#ifndef MPU_TOTIENTS_H
#define MPU_TOTIENTS_H

#include "ptypes.h"

extern UV  totient(UV n);

extern UV* range_totient(UV low, UV high);

extern UV  sumtotient(UV n);
#define HAVE_SUMTOTIENT128 (HAVE_UINT64 && HAVE_UINT128)
#if HAVE_SUMTOTIENT128
extern bool sumtotient128(uint64_t n, uint128_t *sum);
#endif

extern UV  jordan_totient(UV k, UV n);

extern UV  carmichael_lambda(UV n);

extern UV  dedekind_psi(UV n);

extern bool is_totient(UV n);

extern UV  inverse_totient_count(UV n);
extern UV* inverse_totient_list(UV *ntotients, UV n);

#endif
