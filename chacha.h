#ifndef MPU_CHACHA_H
#define MPU_CHACHA_H

#include "ptypes.h"

/*  State info  */

#define STATESZ  16        /* words: 4 constant, 8 key, 2 counter, 2 nonce */
#define KEYSZ    40        /* bytes of user supplied key+nonce */
#define CORESZ   64        /* bytes output by core */
#define BUFSZ    16*CORESZ /* bytes we get at a time (1024) */
typedef struct {
  uint32_t      state[STATESZ];
  unsigned char buf[BUFSZ];
  uint16_t      have;
  char          goodseed;
} chacha_context_t;

/*  API  */

extern void chacha_seed(chacha_context_t *cs, uint32_t bytes, const unsigned char* data, char good);
extern void chacha_rand_bytes(chacha_context_t *cs, uint32_t bytes, unsigned char* data);

extern uint32_t chacha_irand32(chacha_context_t *cs);
extern UV       chacha_irand64(chacha_context_t *cs);

extern int chacha_selftest(void);

#endif
