#ifndef MPU_ISAAC_H
#define MPU_ISAAC_H

#include "ptypes.h"

extern void isaac_seed(uint32_t bytes, const unsigned char* data);
extern void isaac_rand_bytes(uint32_t bytes, unsigned char* data);

extern uint32_t isaac_irand32(void);
extern UV       isaac_irand64(void);

extern int isaac_selftest(void);
#endif
