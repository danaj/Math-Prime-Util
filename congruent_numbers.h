#ifndef MPU_CONGRUENT_NUMBERS_H
#define MPU_CONGRUENT_NUMBERS_H

#include "ptypes.h"


extern bool is_congruent_number(UV n);


/* We supply two functions to allow testing. */

/* Answers 0: non-congruent, 1: congruent, -1: don't know */
extern int  is_congruent_number_filter(UV n);

extern bool is_congruent_number_tunnell(UV n);

#endif
