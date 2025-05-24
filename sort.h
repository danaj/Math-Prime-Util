#ifndef MPU_SORT_H
#define MPU_SORT_H

#include "ptypes.h"

extern void sort_uv_array(UV* L, size_t len);
extern void sort_iv_array(IV* L, size_t len);

extern void sort_dedup_uv_array(UV* L, bool data_is_signed, size_t *len);

#endif
