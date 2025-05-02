#ifndef MPU_TOTIENTS_H
#define MPU_TOTIENTS_H

#include "ptypes.h"

extern UV  totient(UV n);

extern UV* range_totient(UV low, UV high);

extern UV  sumtotient(UV n);
extern int sumtotient128(UV n, UV *hisum, UV *losum);

extern UV  jordan_totient(UV k, UV n);

extern bool is_totient(UV n);

extern UV  inverse_totient_count(UV n);
extern UV* inverse_totient_list(UV *ntotients, UV n);

#endif
