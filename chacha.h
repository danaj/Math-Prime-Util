#ifndef MPU_CHACHA_H
#define MPU_CHACHA_H

#include "ptypes.h"

extern void chacha_seed(uint32_t bytes, const unsigned char* data);
extern void chacha_rand_bytes(uint32_t bytes, unsigned char* data);

extern uint32_t chacha_irand32(void);
extern UV       chacha_irand64(void);

#endif
