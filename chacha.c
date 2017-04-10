/*
 * The ChaCha(20) CSPRNG interface.
 * Using OpenBSD version of merged core, and ideas from their arc4random.
 */

#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include "ptypes.h"
#include "chacha.h"

#define CHACHA_ROUNDS 20
#define KEYSTREAM_ONLY
#include "chacha_private.h"

#define KEYSZ    32
#define IVSZ      8
#define BLOCKSZ  64
#define RSBUFSZ  (16*BLOCKSZ)
#define minimum(a,b)  ((a) < (b) ? (a) : (b))

typedef unsigned char u8;
typedef unsigned int u32;

static struct _rs {
    size_t        rs_have;            /* valid bytes at end of rs_buf */
    size_t        rs_count;           /* bytes till reseed */
} *rs;
static struct _rsx {
    chacha_ctx    rs_chacha;          /* chacha context for random keystream */
    u8            rs_buf[RSBUFSZ];    /* keystream blocks */
} *rsx;
static struct {
    struct _rs  rs;
    struct _rsx rsx;
} static_ctx;

/*****************************************************************************/

static inline void _rs_init(const u8 *buf, size_t n)
{
  if (n < KEYSZ + IVSZ)
    croak("ChaCha bad key size: %lu", (unsigned long)n);

  if (rs == NULL) {
    rs  = &(static_ctx.rs);
    rsx = &(static_ctx.rsx);
    rs->rs_have  = 0;
    rs->rs_count = 0;
  }

  chacha_keysetup(&rsx->rs_chacha, buf, KEYSZ*8, 0);
  chacha_ivsetup(&rsx->rs_chacha, buf + KEYSZ);
}

/* We don't do periodic reseeds with more entropy */
static inline void _rs_stir(void) { }
static inline void _rs_stir_if_needed(size_t len) { }

static inline void _rs_rekey(const u8 *dat, size_t datlen)
{
#ifndef KEYSTREAM_ONLY
  memset(rsx->rs_buf, 0, RSBUFSZ);
#endif
  /* fill rs_buf with the keystream */
  chacha_encrypt_bytes(&rsx->rs_chacha, rsx->rs_buf, rsx->rs_buf, RSBUFSZ);
  /* mix in optional user provided data */
  if (dat) {
    size_t i, m = minimum(datlen, KEYSZ + IVSZ);
    for (i = 0; i < m; i++)
        rsx->rs_buf[i] ^= dat[i];
  }
  /* immediately reinit for backtracking resistance */
  _rs_init(rsx->rs_buf, KEYSZ + IVSZ);
  memset(rsx->rs_buf, 0, KEYSZ + IVSZ);
  rs->rs_have = RSBUFSZ - KEYSZ - IVSZ;
}

void chacha_seed(uint32_t bytes, const unsigned char* data)
{
  _rs_init(data, bytes);
}

void chacha_rand_bytes(uint32_t bytes, unsigned char* data)
{
  u8 *keystream;
  uint32_t m;

  _rs_stir_if_needed(bytes);
  while (bytes > 0) {
    if (rs->rs_have > 0) {
      m = minimum(bytes, rs->rs_have);
      keystream = rsx->rs_buf + RSBUFSZ - rs->rs_have;
      memcpy(data, keystream, m);
      memset(keystream, 0, m);
      data += m;
      bytes -= m;
      rs->rs_have -= m;
    }
    if (rs->rs_have == 0)
      _rs_rekey(NULL, 0);
  }
}

uint32_t chacha_irand32(void)
{
  u8 *keystream;
  uint32_t a;
  _rs_stir_if_needed(4);
  if (rs->rs_have < 4)
    _rs_rekey(NULL, 0);
  keystream = rsx->rs_buf + RSBUFSZ - rs->rs_have;
  memcpy(&a, keystream, 4);
  memset(keystream, 0, 4);
  rs->rs_have -= 4;
  return a;
}
UV chacha_irand64(void)
{
  u8 *keystream;
  UV a;
  _rs_stir_if_needed(8);
  if (rs->rs_have < 8)
    _rs_rekey(NULL, 0);
  keystream = rsx->rs_buf + RSBUFSZ - rs->rs_have;
  memcpy(&a, keystream, 8);
  memset(keystream, 0, 8);
  rs->rs_have -= 8;
  return a;
}
