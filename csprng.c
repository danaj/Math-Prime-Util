
/* Our API for random numbers.
 *
 * We can use ISAAC, ChaCha20, or something else.
 *
 * 4798    ns/word  ChaCha20 in Perl
 * 2106    ns/word  ChaCha8 in Perl
 *  850    ns/word  ISAAC in Perl
 *
 *   11.20 ns/word  ChaCha20 (openbsd)
 *   10.31 ns/word  ChaCha20 (dj)
 *    6.85 ns/word  ChaCha12 (dj)
 *    5.99 ns/word  Tyche
 *    5.11 ns/word  ChaCha8 (dj)
 *    4.37 ns/word  MT19937 (Cokus)
 *    4.14 ns/word  Tyche-i
 *    3.26 ns/word  ISAAC
 *    3.18 ns/word  PCG32 basic
 *    1.95 ns/word  PCG64 opt
 *    1.84 ns/word  ChaCha20 (AVX2 chacha-opt)
 *    1.48 ns/word  Xoroshiro128+
 *    1.16 ns/word  SplitMix64
 *
 * These functions do locking, the underlying library does not.
 */

#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include "ptypes.h"
#include "threadlock.h"
#include "csprng.h"

#define USE_ISAAC      0
#define USE_CHACHA20   1

static int mutex_init = 0;
static int good_seed = 0;
MUTEX_DECL(state);

#if (USE_ISAAC + USE_CHACHA20) != 1
#error Exactly one CSPRNG should been selected
#endif

#if USE_ISAAC

#include "isaac.h"
#define SEED_BYTES 1024
#define CSEED     isaac_seed
#define CRBYTES   isaac_rand_bytes
#define CIRAND32  isaac_irand32
#define CIRAND64  isaac_irand64
#define CSELFTEST isaac_selftest

#elif USE_CHACHA20

#include "chacha.h"
#define SEED_BYTES (32+8)
#define CSEED     chacha_seed
#define CRBYTES   chacha_rand_bytes
#define CIRAND32  chacha_irand32
#define CIRAND64  chacha_irand64
#define CSELFTEST chacha_selftest

#endif

/* Helper macros, similar to ChaCha, so we're consistent. */
#ifndef U8TO32_LE
#define U8TO32_LE(p) \
  (((uint32_t)((p)[0])      ) | \
   ((uint32_t)((p)[1]) <<  8) | \
   ((uint32_t)((p)[2]) << 16) | \
   ((uint32_t)((p)[3]) << 24))
#endif
#define U32TO8_LE(p, v) \
  do { \
    uint32_t _v = v; \
    (p)[0] = (((_v)      ) & 0xFFU); \
    (p)[1] = (((_v) >>  8) & 0xFFU); \
    (p)[2] = (((_v) >> 16) & 0xFFU); \
    (p)[3] = (((_v) >> 24) & 0xFFU); \
  } while (0)

/*****************************************************************************/

/* We put a simple non-CS PRNG here to help fill small seeds.
 * I'd like to use PCG64 or Xoroshiro128+ here, but to be maximally
 * portable and reproducible, we'll use PCG (RXS M XS 32).
 */
#if 0
uint32_t prng_next(char* rng) {
  uint32_t xorshifted, rot;
  uint64_t *state = (uint64_t*) rng;
  uint64_t oldstate = state[0];
  pcg_state = oldstate * 6364136223846793005ULL + state[1];
  xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
  rot = oldstate >> 59u;
  return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}
char* prng_new(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
  uint64_t *state, initstate, initseq;
  New(0, state, 2, uint64_t);
  initstate = (((uint64_t)a) << 32) | b;
  initseq = (((uint64_t)c) << 32) | d;
  state[0] = 0U;
  state[1] = (initseq << 1u) | 1u;
  (void) prng_next((char*)state);
  state[0] += initstate;
  (void) prng_next((char*)state);
  return (char*) state;
}
#else
uint32_t prng_next(char* ctx) {
  uint32_t *rng = (uint32_t*) ctx;
  uint32_t word, oldstate = rng[0];
  rng[0] = rng[0] * 747796405U + rng[1];
  word = ((oldstate >> ((oldstate >> 28u) + 4u)) ^ oldstate) * 277803737u;
  return (word >> 22u) ^ word;
}
char* prng_new(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
  uint32_t *state;
  New(0, state, 2, uint32_t);
  state[0] = 0U;
  state[1] = (b << 1u) | 1u;
  (void) prng_next((char*)state);
  state[0] += a;
  (void) prng_next((char*)state);
  /* TODO: We really should do *something* with c and d */
  return (char*) state;
}
#endif

/*****************************************************************************/

void csprng_seed(uint32_t bytes, const unsigned char* data)
{
  unsigned char seed[SEED_BYTES + 4];

  /* If given a short seed, minimize zeros in state */
  if (bytes >= SEED_BYTES) {
    memcpy(seed, data, SEED_BYTES);
  } else {
    char* rng;
    uint32_t a, b, c, d, i;
    memcpy(seed, data, bytes);
    memset(seed+bytes, 0, sizeof(seed)-bytes);
    a = U8TO32_LE(seed +  0);
    b = U8TO32_LE(seed +  4);
    c = U8TO32_LE(seed +  8);
    d = U8TO32_LE(seed + 12);
    rng = prng_new(a,b,c,d);
    for (i = bytes; i < SEED_BYTES; i += 4)
      U32TO8_LE(seed + i, prng_next(rng));
    Safefree(rng);
#if 0
    printf("got %u bytes in expanded to %u\n", bytes, SEED_BYTES);
    printf("from: ");for(i=0;i<bytes;i++)printf("%02x",data[i]);printf("\n");
    printf("to:   ");for(i=0;i<SEED_BYTES;i++)printf("%02x",seed[i]);printf("\n");
#endif
  }

  if (!mutex_init) {
    MUTEX_INIT(&state_mutex);
    mutex_init = 1;
    MUTEX_LOCK(&state_mutex);
    CSELFTEST();
    MUTEX_UNLOCK(&state_mutex);
  }
  MUTEX_LOCK(&state_mutex);
  CSEED(SEED_BYTES, seed);
  good_seed = (bytes >= 16);
  MUTEX_UNLOCK(&state_mutex);
}

extern void csprng_srand(UV insecure_seed)
{
  unsigned char seed[8];
  U32TO8_LE(seed, insecure_seed);
#if BITS_PER_WORD == 64
  U32TO8_LE(seed + 4, (insecure_seed >> 32));
#endif
  csprng_seed((BITS_PER_WORD == 64) ? 8 : 4, seed);
}

void csprng_rand_bytes(uint32_t bytes, unsigned char* data)
{
  if (!mutex_init) croak("CSPRNG used before init");
  MUTEX_LOCK(&state_mutex);
  CRBYTES(bytes, data);
  MUTEX_UNLOCK(&state_mutex);
}

uint32_t irand32(void)
{
  uint32_t a;
  MUTEX_LOCK(&state_mutex);
  a = CIRAND32();
  MUTEX_UNLOCK(&state_mutex);
  return a;
}
UV irand64(void)
{
#if BITS_PER_WORD < 64
  croak("irand64 too many bits for UV");
#else
  UV a;
  MUTEX_LOCK(&state_mutex);
  a = CIRAND64();
  MUTEX_UNLOCK(&state_mutex);
  return a;
#endif
}


/*****************************************************************************/

int is_csprng_well_seeded(void) { return good_seed; }

/* Many ways to do this.  2 or 3 good-but-different ways, many bad. */
static NV _tonv_32 = 2.3283064365386962890625000000000000000E-10L;
static NV _tonv_64 = 5.4210108624275221700372640043497085571E-20L;
NV drand64(void)
{
  NV r;
  MUTEX_LOCK(&state_mutex);
  if (sizeof(NV) == 4)
    r = CIRAND32() * _tonv_32;
  else if (BITS_PER_WORD == 64)
    r = CIRAND64() * _tonv_64;
  else
    r = CIRAND32() * _tonv_64 + CIRAND32() * _tonv_32;
  MUTEX_UNLOCK(&state_mutex);
  return r;
}

/* Return rand 32-bit integer between 0 to n-1 inclusive */
uint32_t urandomm32(uint32_t n)
{
  uint32_t r, rmin;

  if (n <= 1)
    return 0;

  rmin = -n % n;
  MUTEX_LOCK(&state_mutex);
  while (1) {
    r = CIRAND32();
    if (r >= rmin)
      break;
  }
  MUTEX_UNLOCK(&state_mutex);
  return r % n;
}

UV urandomm64(UV n)
{
  UV r, rmin;

  if (n <= 1)
    return 0;

  rmin = -n % n;
  MUTEX_LOCK(&state_mutex);
  while (1) {
    r = CIRAND64();
    if (r >= rmin)
      break;
  }
  MUTEX_UNLOCK(&state_mutex);
  return r % n;
}

UV urandomb(int nbits)
{
  if (nbits == 0) {
    return 0;
  } else if (nbits <= 32) {
    return irand32() >> (32-nbits);
#if BITS_PER_WORD == 64
  } else if (nbits <= 64) {
    return irand64() >> (64-nbits);
#endif
  }
  croak("irand64 too many bits for UV");
}
