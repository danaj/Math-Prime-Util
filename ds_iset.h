#ifndef MPU_DS_ISET_H
#define MPU_DS_ISET_H

#include "ptypes.h"

/******************************************************************************/
/*                        INTEGER SET DATA STRUCTURE                          */
/******************************************************************************/

typedef struct {
  UV *arr;
  UV mask;
  unsigned long maxsize;
  unsigned long size;
  char contains_zero;
} iset_t;

iset_t iset_create(unsigned long init_size);
void  iset_destroy(iset_t *set);

static unsigned long iset_size(iset_t set) { return set.size; }

int   iset_contains(const iset_t set, UV val);   /* returns 0 or 1 */

int iset_add(iset_t *set, UV val);  /* Returns 1 if added, 0 otherwise */
/* void  iset_remove(iset_t *set, UV val); */

/* We could make an iterator */

void  iset_allvals(const iset_t set, UV* array); /* array must have room */

void iset_test(void);

#endif
