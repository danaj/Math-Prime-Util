#ifndef MPU_ISAAC_H
#define MPU_ISAAC_H

#include "ptypes.h"

extern void isaac_seed(uint32_t bytes, const unsigned char* data);
extern int isaac_well_seeded(void);

/* get random bytes */
extern void isaac_rand_bytes(uint32_t bytes, unsigned char* data);

/* fast 32-bit, 0- to 64-bit, and float */
extern uint32_t irand32(void);
extern UV irand64(void);
extern NV drand64(void);

/* integer less than n */
extern uint32_t isaac_rand(uint32_t n);
extern UV isaac_rand64(UV n);
/* integer with n bits */
extern UV irandb(int nbits);

#endif
