#ifndef MPU_XS_SET_H
#define MPU_XS_SET_H

#define MPU_SC_SIZE  257   /* Choose 131, 257, 521, 1031, 2053 */

typedef struct { /* lo in 0, hi in 1, cached values in rest */
  UV     value[2+MPU_SC_SIZE];
  size_t index[2+MPU_SC_SIZE];
  signed char status[2+MPU_SC_SIZE];
} set_data_t;

int type_of_sumset(int typea, int typeb, UV amin, UV amax, UV bmin, UV bmax);
void _sc_clear_cache(set_data_t *cache);
int _sc_set_lohi(pTHX_ SV** avarr, set_data_t *cache, int loindex, int hiindex, int *lostatus, int *histatus, UV *loval, UV *hival);
int insert_index_in_set(pTHX_ AV* av, set_data_t *cache, int sign, UV val);
int index_in_set(pTHX_ AV* av, set_data_t *cache, int sign, UV val);
int is_in_set(pTHX_ AV* av, set_data_t *cache, int sign, UV val);
int del_from_set(pTHX_ AV* ava, int bstatus, UV b);
int ins_into_set(pTHX_ AV* ava, int bstatus, UV b);

#endif
