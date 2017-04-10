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
#define RUN_INTERNAL_TESTS 0
#define RESEED_ON_REFILL 1

#if 0  /* My code */

#define STATESZ  16        /* words: 4 constant, 8 key, 2 counter, 2 nonce */
#define KEYSZ    40        /* bytes of user supplied key+nonce */
#define CORESZ   64        /* bytes output by core */
#define BUFSZ    16*CORESZ /* bytes we get at a time (1024) */
typedef struct {
  uint32_t      state[STATESZ];
  unsigned char buf[BUFSZ];
  uint32_t      have;
} chacha_context_t;

/*****************************************************************************/
/* Chacha routines: init, quarter round, core, keystream                     */
/*****************************************************************************/

#ifndef U8TO32_LE
#define U8TO32_LE(p) \
  (((uint32_t)((p)[0])      ) | \
   ((uint32_t)((p)[1]) <<  8) | \
   ((uint32_t)((p)[2]) << 16) | \
   ((uint32_t)((p)[3]) << 24))
#endif

static void init_context(chacha_context_t *ctx, const unsigned char *seed, int init_buffer)
{
  uint32_t *x = ctx->state;

  x[ 0] = 0x61707865;
  x[ 1] = 0x3320646e;
  x[ 2] = 0x79622d32;
  x[ 3] = 0x6b206574;
  x[ 4] = U8TO32_LE(seed+ 0);
  x[ 5] = U8TO32_LE(seed+ 4);
  x[ 6] = U8TO32_LE(seed+ 8);
  x[ 7] = U8TO32_LE(seed+12);
  x[ 8] = U8TO32_LE(seed+16);
  x[ 9] = U8TO32_LE(seed+20);
  x[10] = U8TO32_LE(seed+24);
  x[11] = U8TO32_LE(seed+28);
  x[12] = 0;
  x[13] = 0;
  x[14] = U8TO32_LE(seed+32);
  x[15] = U8TO32_LE(seed+36);

  if (init_buffer) {
    memset(ctx->buf, 0, BUFSZ);
    ctx->have = 0;
  }
}

static INLINE uint32_t rotl32(uint32_t x, const unsigned int n) {
  return (x << n) | (x >> (32u - n));
}
#define QUARTERROUND(a,b,c,d) \
  a += b;  d = rotl32(d ^ a, 16); \
  c += d;  b = rotl32(b ^ c, 12); \
  a += b;  d = rotl32(d ^ a,  8); \
  c += d;  b = rotl32(b ^ c,  7); \

/* Produces buffer from state, does not change state */
static void chacha_core(unsigned char* buf, const chacha_context_t *ctx) {
  uint32_t i, x[16];
  const uint32_t *s = ctx->state;

  memcpy(x, s, 16*sizeof(uint32_t));

  for (i = 0; i < CHACHA_ROUNDS; i += 2) {
    QUARTERROUND( x[ 0], x[ 4], x[ 8], x[12] );
    QUARTERROUND( x[ 1], x[ 5], x[ 9], x[13] );
    QUARTERROUND( x[ 2], x[ 6], x[10], x[14] );
    QUARTERROUND( x[ 3], x[ 7], x[11], x[15] );
    QUARTERROUND( x[ 0], x[ 5], x[10], x[15] );
    QUARTERROUND( x[ 1], x[ 6], x[11], x[12] );
    QUARTERROUND( x[ 2], x[ 7], x[ 8], x[13] );
    QUARTERROUND( x[ 3], x[ 4], x[ 9], x[14] );
  }

  for (i = 0; i < 16; i++)
    x[i] += s[i];

  /* TODO: endian */
  for (i = 0; i < 16; i++)
    memcpy(buf+(4*i), x+i, sizeof(uint32_t));
}

static INLINE void increment_chacha_counter(chacha_context_t *ctx) {
  if (++ctx->state[12] == 0)
    ctx->state[13]++;
}

static uint32_t chacha_keystream(unsigned char* buf, uint32_t n, chacha_context_t *ctx) {
  uint32_t r = n;
  while (r >= CORESZ) {
    chacha_core(buf, ctx);
    increment_chacha_counter(ctx);
    buf += CORESZ;
    r -= CORESZ;
  }
  if (r > 0) {
    unsigned char sbuf[CORESZ];
    chacha_core(sbuf, ctx);
    increment_chacha_counter(ctx);
    memcpy(buf, sbuf, r);
  }
  return n;
}

/* The method for refilling our buffer.
 * This includes reseeding policy.
 */
static uint32_t  _refill_buffer(chacha_context_t *ctx) {
#if RESEED_ON_REFILL
  ctx->have = chacha_keystream(ctx->buf, BUFSZ, ctx);
  init_context(ctx, ctx->buf, 0);
  memset(ctx->buf, 0, KEYSZ);
  ctx->have = BUFSZ - KEYSZ;
#else
  ctx->have = chacha_keystream(ctx->buf, BUFSZ, ctx);
#endif
  return ctx->have;
}


/*****************************************************************************/
/*   Test vectors                                                            */
/*****************************************************************************/
#if RUN_INTERNAL_TESTS
static int _test_qr(void) {
  uint32_t i;
  uint32_t tv1i[4] = {0x11111111, 0x01020304, 0x9b8d6f43, 0x01234567};
  uint32_t tv1o[4] = {0xea2a92f4, 0xcb1cf8ce, 0x4581472e, 0x5881c4bb};
  uint32_t tv2i[4] = {0x516461b1, 0x2a5f714c, 0x53372767, 0x3d631689};
  uint32_t tv2o[4] = {0xbdb886dc, 0xcfacafd2, 0xe46bea80, 0xccc07c79};
  if (CHACHA_ROUNDS != 20) return 0;
  QUARTERROUND(tv1i[0],tv1i[1],tv1i[2],tv1i[3]);
  QUARTERROUND(tv2i[0],tv2i[1],tv2i[2],tv2i[3]);
  for (i = 0; i < 4; i++) {
    if (tv1i[i] != tv1o[i]) croak("QR test 2.1.1 fail %u\n",i);
    if (tv2i[i] != tv2o[i]) croak("QR test 2.2.1 fail %u\n",i);
  }
  return 1;
}
static int _test_core(void) {
  uint32_t test, i;
  unsigned char keys[6][40] = { {0},{0},{0},{0},{0} };
  char ebuf[6][129] = {
   "76b8e0ada0f13d90405d6ae55386bd28bdd219b8a08ded1aa836efcc8b770dc7da41597c5157488d7724e03fb8d84a376a43b8f41518a11cc387b669b2ee6586",
   "4540f05a9f1fb296d7736e7b208e3c96eb4fe1834688d2604f450952ed432d41bbe2a0b6ea7566d2a5d1e7e20d42af2c53d792b1c43fea817e9ad275ae546963",
   "de9cba7bf3d69ef5e786dc63973f653a0b49e015adbff7134fcb7df137821031e85a050278a7084527214f73efc7fa5b5277062eb7a0433e445f41e31afab757",
   "ef3fdfd6c61578fbf5cf35bd3dd33b8009631634d21e42ac33960bd138e50d32111e4caf237ee53ca8ad6426194a88545ddc497a0b466e7d6bbdb0041b2f586b",
   "f798a189f195e66982105ffb640bb7757f579da31602fc93ec01ac56f85ac3c134a4547b733b46413042c9440049176905d3be59ea1c53f15916155c2be8241a",
   "10f1e7e4d13b5915500fdd1fa32071c4c7d1f4c733c068030422aa9ac3d46c4ed2826446079faa0914c2d705d98b02a2b5129cd1de164eb9cbd083e8a2503c4e",
  };
  keys[1][31] = 1;
  keys[2][39] = 1;
  keys[3][32] = 1;
  for (i = 0; i < 32; i++) keys[4][ 0+i] = i;
  for (i = 0; i <  8; i++) keys[4][32+i] = i;
  for (i = 0; i < 32; i++) keys[5][ 0+i] = i;
  keys[5][35] = 0x4a;

  if (CHACHA_ROUNDS != 20) return 0;
  for (test = 0; test < 6; test++) {
    unsigned char* key = keys[test];
    char* expout = ebuf[test];
    char got[129];
    chacha_context_t ctx;
    init_context(&ctx, key, 1);
    if (test == 5) { ctx.state[12]=1; ctx.state[13]=0x09000000; }
    chacha_core(ctx.buf, &ctx);
    if (test == 0) {
      for (i = 5; i < 16; i++)
        if (ctx.state[i] != 0)
          croak("core modified state");
    }
    for (i = 0; i < 64; i++)
      sprintf(got+2*i,"%02x", ctx.buf[i]);
    got[128] = '\0';
    if (memcmp(got, expout, 128))
      croak("fail core test vector %u:\n  exp %s\n  got %s\n",test,expout,got);
  }
  return 1;
}
static int _test_keystream(void) {
  uint32_t test, i;
  unsigned char keys[2][40] = { {0},{0} };
  char ebuf[2][1024+1] = {
   "f798a189f195e66982105ffb640bb7757f579da31602fc93ec01ac56f85ac3c134a4547b733b46413042c9440049176905d3be59ea1c53f15916155c2be8241a38008b9a26bc35941e2444177c8ade6689de95264986d95889fb60e84629c9bd9a5acb1cc118be563eb9b3a4a472f82e09a7e778492b562ef7130e88dfe031c79db9d4f7c7a899151b9a475032b63fc385245fe054e3dd5a97a5f576fe064025d3ce042c566ab2c507b138db853e3d6959660996546cc9c4a6eafdc777c040d70eaf46f76dad3979e5c5360c3317166a1c894c94a371876a94df7628fe4eaaf2ccb27d5aaae0ad7ad0f9d4b6ad3b54098746d4524d38407a6deb3ab78fab78c9",
   "af051e40bba0354981329a806a140eafd258a22a6dcb4bb9f6569cb3efe2deaf837bd87ca20b5ba12081a306af0eb35c41a239d20dfc74c81771560d9c9c1e4b224f51f3401bd9e12fde276fb8631ded8c131f823d2c06e27e4fcaec9ef3cf788a3b0aa372600a92b57974cded2b9334794cba40c63e34cdea212c4cf07d41b769a6749f3f630f4122cafe28ec4dc47e26d4346d70b98c73f3e9c53ac40c5945398b6eda1a832c89c167eacd901d7e2bf363",
  };
  for (i = 0; i < 32; i++) keys[0][ 0+i] = i;
  for (i = 0; i <  8; i++) keys[0][32+i] = i;
  for (i = 0; i < 32; i++) keys[1][ 0+i] = i;
  keys[1][35] = 0x4a;

  if (CHACHA_ROUNDS != 20) return 0;

  for (test = 0; test < 2; test++) {
    unsigned char* key = keys[test];
    char* expout = ebuf[test];
    unsigned char kbuf[512];
    char got[1024+1];
    uint32_t gen, len = strlen(expout) / 2;
    chacha_context_t ctx;

    if (len > 512) croak("Test vector too large");
    init_context(&ctx, key, 1);
    gen = chacha_keystream(kbuf, len, &ctx);
    /* Check state block counter */
    for (i = 0; i < len; i++)
      sprintf(got+2*i,"%02x", kbuf[i]);
    got[2*len] = '\0';
    if (memcmp(got, expout, 2*len))
      croak("fail keystream test vector %u:\n  exp %s\n  got %s\n",test,expout,got);
  }
  return 1;
}

int chacha_selftest(void) {
  if (_test_qr() && _test_core() && _test_keystream())
    return 1;
  return 0;
}
#else
int chacha_selftest(void) { return 1; }
#endif

/*****************************************************************************/
/*   API                                                                     */
/*****************************************************************************/

static  chacha_context_t  _cs;

void chacha_seed(uint32_t bytes, const unsigned char* data)
{
  if (bytes < 40) croak("Not enough seed bytes given to ChaCha\n");
  init_context(&_cs, data, 1);
}
void chacha_rand_bytes(uint32_t bytes, unsigned char* data)
{
  while (bytes > 0) {
    uint32_t copybytes;
    if (_cs.have == 0)  _refill_buffer(&_cs);
    copybytes = (bytes > _cs.have)  ?  _cs.have  :  bytes;
    memcpy(data, _cs.buf + BUFSZ - _cs.have, copybytes);
    data += copybytes;
    _cs.have -= copybytes;
    bytes -= copybytes;
  }
}
uint32_t chacha_irand32(void)
{
  uint32_t a;
  unsigned char* ptr;
  if (_cs.have < 4)
    _refill_buffer(&_cs);
  ptr = _cs.buf + BUFSZ - _cs.have;
  _cs.have -= 4;
  a = U8TO32_LE(ptr);
  return a;
}
UV chacha_irand64(void)
{
  uint32_t a,b;
  unsigned char* ptr;
  if (_cs.have < 8)
    _refill_buffer(&_cs);
  ptr = _cs.buf + BUFSZ - _cs.have;
  _cs.have -= 8;
  a = U8TO32_LE(ptr);
  b = U8TO32_LE(ptr);
  return (((UV)b) << 32) | a;
}

/*****************************************************************************/
/*****************************************************************************/
#else  /* OpenBSD style */
/*****************************************************************************/
/*****************************************************************************/

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
#if RUN_INTERNAL_TESTS
  rs->rs_have = RSBUFSZ;
#else
  /* immediately reinit for backtracking resistance */
  _rs_init(rsx->rs_buf, KEYSZ + IVSZ);
  memset(rsx->rs_buf, 0, KEYSZ + IVSZ);
  rs->rs_have = RSBUFSZ - KEYSZ - IVSZ;
#endif
}

int chacha_selftest(void) { return 1; }

void chacha_seed(uint32_t bytes, const unsigned char* data)
{
  _rs_init(data, bytes);
}

void chacha_rand_bytes(uint32_t bytes, unsigned char* data)
{
  u8 *keystream;
  uint32_t m;

#if RUN_INTERNAL_TESTS
  {
    uint32_t i;
    unsigned char key[40] = {0};
    //char expout[] = "76b8e0ada0f13d90405d6ae55386bd28bdd219b8a08ded1aa836efcc8b770dc7da41597c5157488d7724e03fb8d84a376a43b8f41518a11cc387b669b2ee6586";
    char expout[] = "f798a189f195e66982105ffb640bb7757f579da31602fc93ec01ac56f85ac3c134a4547b733b46413042c9440049176905d3be59ea1c53f15916155c2be8241a38008b9a26bc35941e2444177c8ade6689de95264986d95889fb60e84629c9bd9a5acb1cc118be563eb9b3a4a472f82e09a7e778492b562ef7130e88dfe031c79db9d4f7c7a899151b9a475032b63fc385245fe054e3dd5a97a5f576fe064025d3ce042c566ab2c507b138db853e3d6959660996546cc9c4a6eafdc777c040d70eaf46f76dad3979e5c5360c3317166a1c894c94a371876a94df7628fe4eaaf2ccb27d5aaae0ad7ad0f9d4b6ad3b54098746d4524d38407a6deb3ab78fab78c9";
    for (i = 0; i < 32; i++) key[ 0+i] = i;
    for (i = 0; i <  8; i++) key[32+i] = i;

    unsigned char kbuf[512];
    unsigned char* kptr = kbuf;
    char got[1024+1];
    uint32_t gen, len = strlen(expout) / 2;

    if (len > 512) croak("Test vector too large");
    _rs_init(key, KEYSZ + IVSZ);

    //  gen = chacha_keystream(kbuf, len, &ctx);
    gen = len;
    while (gen > 0) {
      if (rs->rs_have > 0) {
        m = minimum(gen, rs->rs_have);
        keystream = rsx->rs_buf + RSBUFSZ - rs->rs_have;
        memcpy(kptr, keystream, m);
        memset(keystream, 0, m);
        kptr += m;
        gen -= m;
        rs->rs_have -= m;
      }
      if (rs->rs_have == 0)
        _rs_rekey(NULL, 0);
    }

    for (i = 0; i < len; i++)
      sprintf(got+2*i,"%02x", kbuf[i]);
    got[2*len] = '\0';
    if (memcmp(got, expout, 2*len))
      croak("fail keystream test vector %u:\n  exp %s\n  got %s\n",0,expout,got);
    printf("Test vector length %d passed\n",(int)len);
  }
#endif

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
#endif
