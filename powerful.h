#ifndef MPU_POWERFUL_H
#define MPU_POWERFUL_H

#include "ptypes.h"

extern int is_powerful(UV n, UV k);

extern UV powerful_count(UV n, UV k);

extern UV nth_powerful(UV n, UV k);

extern UV sumpowerful(UV n, UV k);

extern UV* powerful_numbers_range(UV* npowerful, UV lo, UV hi, UV k);

#endif
