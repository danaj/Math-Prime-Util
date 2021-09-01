#ifndef MPU_LUCKY_NUMBERS_H
#define MPU_LUCKY_NUMBERS_H

#include "ptypes.h"

extern int is_lucky(UV n);

extern UV* lucky_sieve(UV *size, UV n);
extern uint32_t* lucky_sieve32(UV *size, UV n);

extern UV lucky_count_range(UV lo, UV hi);

extern UV lucky_count(UV n);
extern UV lucky_count_upper(UV n);
extern UV lucky_count_lower(UV n);
extern UV lucky_count_approx(UV n);

extern UV nth_lucky(UV n);
extern UV nth_lucky_lower(UV n);
extern UV nth_lucky_upper(UV n);
extern UV nth_lucky_approx(UV n);

#endif
