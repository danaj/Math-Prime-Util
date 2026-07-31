#ifndef MPU_XS_SET_H
#define MPU_XS_SET_H

typedef enum {
  SET_REL_DISJOINT = 0,
  SET_REL_EQUAL,
  SET_REL_SUBSET,
  SET_REL_PROPER_SUBSET,
  SET_REL_SUPERSET,
  SET_REL_PROPER_SUPERSET,
  SET_REL_PROPER_INTERSECTION
} set_relation_op_t;
bool xs_set_relation(pTHX_ SV* sva, SV* svb, set_relation_op_t op, int *ret, const char *name);

typedef enum {
  SET_OP_UNION = 0,
  SET_OP_INTERSECT,
  SET_OP_MINUS,
  SET_OP_DELTA
} set_op_t;
bool xs_set_op(pTHX_ SV* sva, SV* svb, set_op_t op, SV **ret, const char *name);
bool xs_is_sidon_set(pTHX_ SV* sva, int *ret);
bool xs_is_sumfree_set(pTHX_ SV* sva, int *ret);

#define MPU_SC_TREE_SIZE  511  /* Top nine levels of the binary-search tree. */
#define MPU_SC_TAIL_SIZE   64  /* Power-of-two cache for deeper tree nodes. */
#if MPU_SC_TAIL_SIZE == 0 || \
    (MPU_SC_TAIL_SIZE & (MPU_SC_TAIL_SIZE-1)) != 0
#error "MPU_SC_TAIL_SIZE must be a power of two"
#endif
typedef struct {
  AV *av;
  Size_t len;
  UV lo_value, hi_value;
  signed char lo_status, hi_status;
  UV tree_value[MPU_SC_TREE_SIZE];
  signed char tree_status[MPU_SC_TREE_SIZE];
  UV tail_value[MPU_SC_TAIL_SIZE];
  Size_t tail_index[MPU_SC_TAIL_SIZE];
  signed char tail_status[MPU_SC_TAIL_SIZE];
} set_cache_t;

void _sc_init_cache(set_cache_t *cache);
int _sc_set_bounds(pTHX_ AV* av, set_cache_t *cache, Size_t len,
                   int *lostatus, int *histatus, UV *loval, UV *hival);

SSize_t insert_index_in_set(pTHX_ AV* av, set_cache_t *cache, int sign, UV val);
SSize_t index_in_set(pTHX_ AV* av, set_cache_t *cache, int sign, UV val);
int is_in_set(pTHX_ AV* av, set_cache_t *cache, int sign, UV val);
int del_from_set(pTHX_ AV* ava, int bstatus, UV b);
int ins_into_set(pTHX_ AV* ava, int bstatus, UV b);
int type_of_sumset(int typea, int typeb, UV amin, UV amax, UV bmin, UV bmax);

#endif
