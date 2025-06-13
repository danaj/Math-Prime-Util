
/* Our API for random numbers.
 *
 * We can use ISAAC, ChaCha20, or something else.
 *
 * 3700    ns/word  ChaCha20 in Perl
 * 3100    ns/word  Salsa20 in Perl
 * 1600    ns/word  ChaCha8 in Perl
 *  760    ns/word  ISAAC in Perl
 *
 *   11.20 ns/word  ChaCha20 (openbsd)
 *   10.31 ns/word  ChaCha20 (dj)
 *    8.66 ns/word  ChaCha20 (sse2 Peters)
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
#include "csprng.h"

#include "chacha.h"
#define SEED_BYTES (32+8)
#define CSEED(ctx,bytes,data,good)  chacha_seed(ctx,bytes,data,good)
#define CRBYTES(ctx,bytes,data)     chacha_rand_bytes(ctx,bytes,data)
#define CIRAND32(ctx)               chacha_irand32(ctx)
#define CIRAND64(ctx)               chacha_irand64(ctx)
#define CSELFTEST()                 chacha_selftest()

/* Helper macros, similar to ChaCha, so we're consistent. */
#if !defined(__x86_64__)
#undef U8TO32_LE
#undef U32TO8_LE
#endif
#ifndef U8TO32_LE
#define U8TO32_LE(p) \
  ((uint32_t)(p)[0]       | \
   (uint32_t)(p)[1] <<  8 | \
   (uint32_t)(p)[2] << 16 | \
   (uint32_t)(p)[3] << 24)
#endif
#ifndef U32TO8_LE
#define U32TO8_LE(p, v) \
  do { uint32_t _v = v; \
       (p)[0] = _v       & 0xFF; \
       (p)[1] = _v >>  8 & 0xFF; \
       (p)[2] = _v >> 16 & 0xFF; \
       (p)[3] = _v >> 24 & 0xFF; \
  } while (0)
#endif

/*****************************************************************************/

/* We put a simple 32-bit non-CS PRNG here to help fill small seeds. */
#if 0
/* XOSHIRO128**  32-bit output, 32-bit types, 128-bit state */
static INLINE uint32_t rotl(const uint32_t x, int k) {
  return (x << k) | (x >> (32 - k));
}
uint32_t prng_next(void* ctx) {
  uint32_t *s = (uint32_t*) ctx;
  const uint32_t result_starstar = rotl(s[0] * 5, 7) * 9;
  const uint32_t t = s[1] << 9;
  s[2] ^= s[0];  s[3] ^= s[1];  s[1] ^= s[2];  s[0] ^= s[3];
  s[2] ^= t;
  s[3] = rotl(s[3], 11);
  return result_starstar;
}
void* prng_new(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
  uint32_t *state;
  New(0, state, 4, uint32_t);
  state[0] = 1;  state[1] = b;  state[2] = c;  state[3] = d;
  (void) prng_next((void*)state);
  state[0] += a;
  (void) prng_next((void*)state);
  return (void*) state;
}
#else
/* PCG RXS M XS 32.  32-bit output, 32-bit state and types. */
uint32_t prng_next(void* ctx) {
  uint32_t *rng = (uint32_t*) ctx;
  uint32_t word, oldstate = rng[0];
  rng[0] = rng[0] * 747796405U + rng[1];
  word = ((oldstate >> ((oldstate >> 28u) + 4u)) ^ oldstate) * 277803737u;
  return (word >> 22u) ^ word;
}
void* prng_new(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
  uint32_t *state;
  New(0, state, 2, uint32_t);
  state[0] = 0U;
  state[1] = (b << 1u) | 1u;
  (void) prng_next((void*)state);
  state[0] += a;
  (void) prng_next((void*)state);
  state[0] ^= c;
  (void) prng_next((void*)state);
  state[0] ^= d;
  (void) prng_next((void*)state);
  return (void*) state;
}
#endif

/*****************************************************************************/

uint32_t csprng_context_size(void)
{
  return sizeof(chacha_context_t);
}
static char _has_selftest_run = 0;

void csprng_seed(void *ctx, uint32_t bytes, const unsigned char* data)
{
  unsigned char seed[SEED_BYTES + 4];

  /* If given a short seed, minimize zeros in state */
  if (bytes >= SEED_BYTES) {
    memcpy(seed, data, SEED_BYTES);
  } else {
    void* rng;
    uint32_t a, b, c, d, i;
    memcpy(seed, data, bytes);
    memset(seed+bytes, 0, sizeof(seed)-bytes);
    a = U8TO32_LE(seed +  0);
    b = U8TO32_LE(seed +  4);
    c = U8TO32_LE(seed +  8);
    d = U8TO32_LE(seed + 12);
    rng = prng_new(a,b,c,d);
    for (i = 4*((bytes+3)/4); i < SEED_BYTES; i += 4)
      U32TO8_LE(seed + i, prng_next(rng));
    Safefree(rng);
#if 0
    printf("got %u bytes in expanded to %u\n", bytes, SEED_BYTES);
    printf("from: ");for(i=0;i<bytes;i++)printf("%02x",data[i]);printf("\n");
    printf("to:   ");for(i=0;i<SEED_BYTES;i++)printf("%02x",seed[i]);printf("\n");
#endif
  }

  if (!_has_selftest_run) {
    _has_selftest_run = 1;
    CSELFTEST();
  }
  CSEED(ctx, SEED_BYTES, seed, (bytes >= 16));
}

extern void csprng_srand(void* ctx, UV insecure_seed)
{
#if BITS_PER_WORD == 32
  unsigned char seed[4] = {0};
  U32TO8_LE(seed, insecure_seed);
  csprng_seed(ctx, 4, seed);
#else
  unsigned char seed[8] = {0};
  if (insecure_seed <= UVCONST(4294967295)) {
    U32TO8_LE(seed, insecure_seed);
    csprng_seed(ctx, 4, seed);
  } else {
    U32TO8_LE(seed, insecure_seed);
    U32TO8_LE(seed + 4, (insecure_seed >> 32));
    csprng_seed(ctx, 8, seed);
  }
#endif
}

void csprng_rand_bytes(void* ctx, uint32_t bytes, unsigned char* data)
{
  CRBYTES(ctx, bytes, data);
}

uint32_t irand32(void* ctx)
{
  return CIRAND32(ctx);
}
UV irand64(void* ctx)
{
#if BITS_PER_WORD < 64
  croak("irand64 too many bits for UV");
#else
  return CIRAND64(ctx);
#endif
}


/*****************************************************************************/

bool is_csprng_well_seeded(void *ctx)
{
  chacha_context_t *cs = ctx;
  return cs->goodseed;
}

/* There are many ways to get floats from integers.  A few good, many bad.
 *
 * Vigna in https://prng.di.unimi.it recommends this C99:
 *    #include <stdint.h>
 *    (x64 >> 11) * 0x1.0p-53
 * Or the older:
 *    (x64 >> 11) * (1.0 / (1ULL<<53)).
 *
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
 * We prefer the x64 / 2^-64 method.  It seems to produce the best results
 * and is easiest for ensuring we fill up all the bits.
 * It is similar to what Geoff Kuenning does in MTwist, though he computes
 * the constants at runtime to ensure a dodgy compiler won't munge them.
 */
#define TO_NV_32    2.3283064365386962890625000000000000000E-10L
#define TO_NV_64    5.4210108624275221700372640043497085571E-20L
#define TO_NV_96    1.2621774483536188886587657044524579675E-29L
#define TO_NV_128   2.9387358770557187699218413430556141945E-39L

#define DRAND_32_32  (CIRAND32(ctx) * TO_NV_32)
#define DRAND_64_32  (((CIRAND32(ctx)>>5) * 67108864.0 + (CIRAND32(ctx)>>6)) / 9007199254740992.0)
#define DRAND_64_64  (CIRAND64(ctx) * TO_NV_64)
#define DRAND_128_32 (CIRAND32(ctx) * TO_NV_32 + CIRAND32(ctx) * TO_NV_64 + CIRAND32(ctx) * TO_NV_96 + CIRAND32(ctx) * TO_NV_128)
#define DRAND_128_64 (CIRAND64(ctx) * TO_NV_64 + CIRAND64(ctx) * TO_NV_128)

NV drand64(void* ctx)
{
  NV r;
#if NVMANTBITS <= 32
  r = DRAND_32_32;
#elif NVMANTBITS <= 64
  r = (BITS_PER_WORD <= 32)  ?  DRAND_64_32  :  DRAND_64_64;
#else
  r = (BITS_PER_WORD <= 32)  ?  DRAND_128_32 :  DRAND_128_64;
#endif
  return r;
}

/* Return rand 32-bit integer between 0 to n-1 inclusive */
uint32_t urandomm32(void *ctx, uint32_t n)
{
  uint32_t r, rmin;

  if (n <= 1)
    return 0;

  rmin = -n % n;
  while (1) {
    r = CIRAND32(ctx);
    if (r >= rmin)
      break;
  }
  return r % n;
}

UV urandomm64(void* ctx, UV n)
{
  UV r, rmin;

  if (n <= 4294967295UL)
    return urandomm32(ctx,n);
  if (n-1 == 4294967295UL)
    return irand32(ctx);

  rmin = -n % n;
  while (1) {
    r = CIRAND64(ctx);
    if (r >= rmin)
      break;
  }
  return r % n;
}

UV urandomb(void* ctx, int nbits)
{
  if (nbits == 0) {
    return 0;
  } else if (nbits <= 32) {
    return irand32(ctx) >> (32-nbits);
#if BITS_PER_WORD == 64
  } else if (nbits <= 64) {
    return irand64(ctx) >> (64-nbits);
#endif
  }
  croak("irand64 too many bits for UV");
}
