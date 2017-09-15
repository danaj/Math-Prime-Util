#ifndef MPU_CSPRNG_H
#define MPU_CSPRNG_H

#include "ptypes.h"

/*****************************************************************************/

extern uint32_t csprng_context_size(void);

/* Seed and init if needed */
extern void csprng_seed(void *ctx, uint32_t bytes, const unsigned char* data);

/* Simple seed */
extern void csprng_srand(void *ctx, UV insecure_seed);

/* Fill buffer with this many bytes of random data */
extern void csprng_rand_bytes(void *ctx, uint32_t bytes, unsigned char* data);

extern uint32_t irand32(void *ctx);
extern UV       irand64(void *ctx);

/*****************************************************************************/

extern int is_csprng_well_seeded(void *ctx);

extern NV drand64(void *ctx);

extern uint32_t urandomm32(void* ctx, uint32_t n);   /* integer less than n */
extern UV       urandomm64(void* ctx, UV n);
extern UV       urandomb(void* ctx, int nbits);      /* integer with n bits */

/*****************************************************************************/

/* Very simple PRNG for other use. */
extern char* prng_new(uint32_t a, uint32_t b, uint32_t c, uint32_t d);
extern uint32_t prng_next(char* rng);

#endif
