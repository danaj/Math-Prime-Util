#ifndef MPU_CSPRNG_H
#define MPU_CSPRNG_H

#include "ptypes.h"

/*****************************************************************************/

extern uint32_t csprng_context_size(void);

typedef UV (*csprng_entropy_f)(UV bytes, unsigned char *buf);

/* Initialize with an entropy source and immediately seed the context. */
extern void csprng_init(void *ctx, csprng_entropy_f get_entropy);

/* Immediately reseed the context from its configured entropy source. */
extern void csprng_reseed(void *ctx);

/* Clear all sensitive state.  The caller still owns the context memory. */
extern void csprng_clear(void *ctx);

/* Discard buffered output and require reseeding before further use. */
extern void csprng_require_reseed(void *ctx);

/* Seed an initialized context, expanding short input if needed. */
extern void csprng_seed(void *ctx, uint32_t bytes, const unsigned char* data);

/* Simple seed.  Returns the number of bytes written to seed[8]. */
extern uint32_t csprng_srand(void *ctx, UV insecure_seed,
                            unsigned char seed[8]);

/* Fill buffer with this many bytes of random data */
extern void csprng_rand_bytes(void *ctx, uint32_t bytes, unsigned char* data);

extern uint32_t irand32(void *ctx);
extern UV       irand64(void *ctx);

/*****************************************************************************/

extern bool is_csprng_well_seeded(void *ctx);

extern NV drand64(void *ctx);

extern uint32_t urandomm32(void* ctx, uint32_t n);   /* integer less than n */
extern UV       urandomm64(void* ctx, UV n);
extern UV       urandomb(void* ctx, uint32_t nbits); /* integer with n bits */

/*****************************************************************************/

/* Very simple PRNG for other use. */
extern void* prng_new(uint32_t a, uint32_t b, uint32_t c, uint32_t d);
extern uint32_t prng_next(void* rng);

#endif
