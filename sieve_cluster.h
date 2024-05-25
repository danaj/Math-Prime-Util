#ifndef MPU_SIEVE_CLUSTER_H
#define MPU_SIEVE_CLUSTER_H

#include "ptypes.h"

extern UV* sieve_cluster_simple(UV beg, UV end, uint32_t nc, const uint32_t* cl, UV* numret);
extern UV* sieve_cluster(UV beg, UV end, uint32_t nc, const uint32_t* cl, UV* numret);

#endif
