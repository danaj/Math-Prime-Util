#ifndef MPU_GOLDBACH_H
#define MPU_GOLDBACH_H

#include "ptypes.h"

extern UV minimal_goldbach_pair(UV n);

extern UV* goldbach_pairs(size_t *size, UV n);

extern UV goldbach_pair_count(UV n);

#endif
