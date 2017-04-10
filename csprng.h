#ifndef MPU_CSPRNG_H
#define MPU_CSPRNG_H

#include "ptypes.h"

/*****************************************************************************/

/* Seed and init if needed */
extern void csprng_seed(uint32_t bytes, const unsigned char* data);

/* Simple seed */
extern void csprng_srand(UV insecure_seed);

/* Fill buffer with this many bytes of random data */
extern void csprng_rand_bytes(uint32_t bytes, unsigned char* data);

extern uint32_t irand32(void);
extern UV       irand64(void);

/*****************************************************************************/

extern int is_csprng_well_seeded(void);

extern NV drand64(void);

extern uint32_t urandomm32(uint32_t n);   /* integer less than n */
extern UV       urandomm64(UV n);
extern UV       urandomb(int nbits);      /* integer with n bits */

/*****************************************************************************/

/* Very simple PRNG for other use. */
extern char* prng_new(uint32_t a, uint32_t b, uint32_t c, uint32_t d);
extern uint32_t prng_next(char* rng);

#endif
