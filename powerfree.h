#ifndef MPU_POWERFREE_H
#define MPU_POWERFREE_H

extern bool is_powerfree(UV n, uint32_t k);

extern UV powerfree_count(UV n, uint32_t k);
extern UV powerfree_sum(UV n, uint32_t k);

extern UV powerfree_part(UV n, uint32_t k);

extern UV powerfree_part_sum(UV n, uint32_t k);

extern UV nth_powerfree(UV n, uint32_t k);

extern UV  squarefree_kernel(UV n);

#endif
