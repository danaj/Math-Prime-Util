#ifndef MPU_DS_ISET_H
#define MPU_DS_ISET_H

#include "ptypes.h"

/******************************************************************************/
/*                        INTEGER SET DATA STRUCTURE                          */
/******************************************************************************/

#define ISET_TYPE_ANY      0
#define ISET_TYPE_UV       1
#define ISET_TYPE_IV       2
#define ISET_TYPE_INVALID  3

typedef struct {
  UV *arr;
  UV mask;
  unsigned long maxsize;
  unsigned long size;
  char contains_zero;
  unsigned char type;
} iset_t;

iset_t iset_create(unsigned long init_size);
void  iset_destroy(iset_t *set);

iset_t iset_create_from_array(UV* d, unsigned long dlen, int dsign);

/* Returns 1 if unsigned, -1 if signed, 0 if messed up. */
static int iset_sign(const iset_t set) {
  static const char _iset_typeret[4] = {1,1,-1,0};
  return _iset_typeret[set.type];
}
static int iset_is_invalid(const iset_t set)
  { return set.type == ISET_TYPE_INVALID; }

static unsigned long iset_size(const iset_t set) { return set.size; }

int iset_contains(const iset_t set, UV val);   /* returns 0 or 1 */

/* void iset_minmax(const iset_t set, UV *min, UV *max); */

/* sign indicates: val is a UV (1) or IV (-1) */
int iset_add(iset_t *set, UV val, int sign); /* Returns 1 if added, 0 if not */
/* void  iset_remove(iset_t *set, UV val); */

/* We could make an iterator */

/* caller supplied array must have room */
void iset_allvals(const iset_t set, UV* array);

void iset_union_with(iset_t *set, const iset_t L);
void iset_intersect_with(iset_t *set, const iset_t L);
void iset_difference_with(iset_t *set, const iset_t L);
void iset_symdiff_with(iset_t *set, const iset_t L);

iset_t iset_union_of(const iset_t A, const iset_t B);
iset_t iset_intersection_of(const iset_t A, const iset_t B);
iset_t iset_difference_of(const iset_t A, const iset_t B);
iset_t iset_symdiff_of(const iset_t A, const iset_t B);
int iset_is_subset_of(const iset_t A, const iset_t B);  /* A subset of B? */

void iset_test(void);

#endif
