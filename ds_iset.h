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
  char seen_uv;
  char sign;
} iset_t;

iset_t iset_create(unsigned long init_size);
void  iset_destroy(iset_t *set);

static unsigned long iset_size(const iset_t set) { return set.size; }

/* Returns 1 if unsigned, -1 if signed, 0 if messed up. */
static int iset_sign(const iset_t set) { return set.sign; }

int iset_contains(const iset_t set, UV val);   /* returns 0 or 1 */

/* sign indicates: val is a UV (1) or IV (-1) */
int iset_add(iset_t *set, UV val, int sign); /* Returns 1 if added, 0 if not */
/* void  iset_remove(iset_t *set, UV val); */

/* We could make an iterator */

/* caller supplied array must have room */
void iset_allvals(const iset_t set, UV* array);

void iset_union_with(iset_t *set, const iset_t L);
void iset_intersect_with(iset_t *set, const iset_t L);

void iset_test(void);

#endif
