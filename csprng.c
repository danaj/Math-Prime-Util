
/* Our API for random numbers.
 *
 * We can use ISAAC, ChaCha20, or something else.
 *
 * 4770    ns/word  ChaCha20 in Perl
 * 2050    ns/word  ChaCha8 in Perl
 *  760    ns/word  ISAAC in Perl
 *
 *   11.20 ns/word  ChaCha20 (openbsd)
 *   10.31 ns/word  ChaCha20 (dj)
 *    6.85 ns/word  ChaCha12 (dj)
 *    5.99 ns/word  Tyche
 *    5.11 ns/word  ChaCha8 (dj)
 *    4.37 ns/word  MT19937 (Cokus)
 *    4.14 ns/word  Tyche-i
 *    3.26 ns/word  ISAAC
 *    3.18 ns/word  PCG64 (64-bit state, 64-bit types)
 *    1.95 ns/word  PCG64 (64-bit state, 128-bit types)
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
/* PCG XSH RR 64/32.  32-bit output, 64-bit state and types. */
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
/* PCG RXS M XS 32.  32-bit output, 32-bit state and types. */
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
  state[0] ^= c;
  (void) prng_next((char*)state);
  state[0] ^= d;
  (void) prng_next((char*)state);
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
  unsigned char seed[8] = {0};
  U32TO8_LE(seed, insecure_seed);
#if BITS_PER_WORD == 64
  U32TO8_LE(seed + 4, (insecure_seed >> 32));
#endif
  csprng_seed(8, seed);
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

/* There are many ways to get floats from integers.  A few good, many bad.
 *
 * Vigna recommends (x64 >> 11) * (1.0 / (1ULL<<53)).
 * http://xoroshiro.di.unimi.it
 * Also see alternatives discussed:
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/speed-up-real.html
 *
 * Melissa O'Neill notes the problem is harder than it looks, doesn't address.
 * http://www.pcg-random.org/pdf/toms-oneill-pcg-family-v1.02.pdf
 *
 * randomstate for numpy uses separate code for each generator.
 * With the exception of dSFMT, they each one one of:
 *     (x64 >> 11) * (1 / 9007199254740992.0)
 *     ((x32 >> 5) * 67108864.0 + (y32 >> 6)) / 9007199254740992.0
 * where the first one is identical to Vigna.
 *
 * David Jones recommends the minor 32-bit variant:
 *     ((x32 >> 6) * 134217728.0 + (y32 >> 5)) / 9007199254740992.0
 * http://www0.cs.ucl.ac.uk/staff/d.jones/GoodPracticeRNG.pdf
 *
 * Taylor Campbell discusses this in:
 * http://mumble.net/~campbell/tmp/random_real.c
 * He points out that there are two common non-broken choices,
 * div by 2^-53  or   div by 2^-64, and each are slightly flawed in
 * different ways.  He shows a theoretically better method.
 */

/*
 * We choose the x64 / 2^-64 method here.  This way is simple and yields
 * 64-bit results if NV's are long doubles.  It is similar to what
 * Geoff Kuenning does in MTwist, though he has the host calculate the
 * constants to ensure a dodgy compiler won't munge them.
 */
static const NV _tonv_32  = 2.3283064365386962890625000000000000000E-10L;
static const NV _tonv_64  = 5.4210108624275221700372640043497085571E-20L;
static const NV _tonv_96  = 1.2621774483536188886587657044524579675E-29L;
static const NV _tonv_128 = 2.9387358770557187699218413430556141945E-39L;
#ifndef NVMANTBITS
#define NVMANTBITS ((sizeof(NV) <= 4) ? 24 : (sizeof(NV) <= 8) ? 52 : 112)
#endif
NV drand64(void)
{
  NV r;
  MUTEX_LOCK(&state_mutex);
#if NVMANTBITS <= 32
  r = CIRAND32() * _tonv_32;
#elif NVMANTBITS <= 64
  #if BITS_PER_WORD > 32
    r = CIRAND64() * _tonv_64;
  #else
    r = ((CIRAND32() >> 5) * 67108864.0 + (CIRAND32() >> 6)) / 9007199254740992.0;
  #endif
#else
  #if BITS_PER_WORD > 32
    r = CIRAND64() * _tonv_64 + CIRAND64() * _tonv_128;
  #else
    r = CIRAND32() * _tonv_32 + CIRAND32() * _tonv_64 + CIRAND32() * _tonv_96 + CIRAND32() * _tonv_128;
  #endif
#endif
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
