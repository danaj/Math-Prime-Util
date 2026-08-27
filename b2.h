#ifndef MPU_B2_H
#define MPU_B2_H

#include "ptypes.h"

/* Little-endian binary limbs with a double-width accumulator. */
#if HAVE_UINT64
  typedef uint32_t b2limb_t;
  typedef uint64_t b2acc_t;
  #define B2_LIMB_BITS       32
  #define B2_DEC_BASE        UINT32_C(1000000000)
  #define B2_DEC_DIGITS      9
  #define B2_DEC_BOUND_DIGITS 10
#else
  typedef uint16_t b2limb_t;
  typedef uint32_t b2acc_t;
  #define B2_LIMB_BITS       16
  #define B2_DEC_BASE        UINT32_C(10000)
  #define B2_DEC_DIGITS      4
  #define B2_DEC_BOUND_DIGITS 5
#endif

#define B2_LIMB_MAX     ((b2limb_t)~(b2limb_t)0)
#define B2_INLINE_LIMBS ((2 * BITS_PER_WORD) / B2_LIMB_BITS)

typedef struct {
  b2limb_t *d;          /* d_small when inline, otherwise heap storage */
  b2limb_t  d_small[B2_INLINE_LIMBS];
  uint32_t  alloc;
  uint32_t  n;          /* significant limbs; zero means value is zero */
  bool      neg;
} b2_t;

MAYBE_UNUSED static INLINE void b2_init(b2_t *x)
{
  x->d = x->d_small;
  x->alloc = B2_INLINE_LIMBS;
  x->n = 0;
  x->neg = 0;
}

MAYBE_UNUSED static INLINE bool b2_is_zero(const b2_t *x)
{
  return x->n == 0;
}

MAYBE_UNUSED static INLINE bool b2_is_one(const b2_t *x)
{
  return !x->neg && x->n == 1 && x->d[0] == 1;
}

MAYBE_UNUSED static INLINE bool b2_is_negative(const b2_t *x)
{
  return x->neg;
}

MAYBE_UNUSED static INLINE bool b2_is_even(const b2_t *x)
{
  return x->n == 0 || !(x->d[0] & 1);
}

MAYBE_UNUSED static INLINE void b2_abs(b2_t *x)
{
  x->neg = 0;
}

MAYBE_UNUSED static INLINE void b2_set_negative(b2_t *x, bool neg)
{
  x->neg = x->n != 0 && neg;
}

void b2_free(b2_t *x);
void b2_neg(b2_t *x);
void b2_ensure(b2_t *x, uint32_t need);
void b2_set(b2_t *dst, const b2_t *src);
void b2_move(b2_t *dst, b2_t *src);
void b2_swap(b2_t *a, b2_t *b);

void b2_set_uv(b2_t *x, UV v);
/* The caller must establish that x fits in a UV. */
UV b2_get_uv(const b2_t *x);
void b2_set_str(b2_t *x, const char *s, STRLEN len);
/* Set from an unsigned, validated base-2, base-8, or base-16 string. */
void b2_set_radix(b2_t *x, const char *s, STRLEN len, uint32_t base);
/* Upper bound excluding the trailing NUL; b2_get_str returns actual length. */
STRLEN b2_string_size(const b2_t *x);
STRLEN b2_get_str(char *buf, const b2_t *x);
/* Return 0 without writing when the canonical result exceeds limit. */
STRLEN b2_get_str_limit(char *buf, const b2_t *x, STRLEN limit);

void b2_init_set(b2_t *x, const b2_t *y);
void b2_init_set_uv(b2_t *x, UV v);
void b2_init_set_str(b2_t *x, const char *s, STRLEN len);
void b2_init_set_radix(b2_t *x, const char *s, STRLEN len, uint32_t base);
void b2_init_set_pow2(b2_t *x, uint32_t k);

int b2_cmp_abs(const b2_t *a, const b2_t *b);
int b2_cmp(const b2_t *a, const b2_t *b);
size_t b2_bit_length(const b2_t *x);

void b2_add(b2_t *out, const b2_t *a, const b2_t *b);
void b2_mul(b2_t *out, const b2_t *a, const b2_t *b);
void b2_sqr(b2_t *out, const b2_t *a);
void b2_pow(b2_t *out, const b2_t *a, UV exp);
void b2_product(b2_t A[], size_t a, size_t b);
void b2_product_u32(b2_t *out, const uint32_t A[], size_t len);

uint32_t b2_mod_u32(const b2_t *a, uint32_t p);
void b2_divexact_u32(b2_t *a, uint32_t p);

void b2_mul_2exp(b2_t *out, const b2_t *a, UV bits);
void b2_tdiv_2exp(b2_t *out, const b2_t *a, UV bits);
void b2_fdiv_2exp(b2_t *out, const b2_t *a, UV bits);

UV b2_popcount(const b2_t *x);
bool b2_testbit(const b2_t *x, uint32_t k);
void b2_setbit(b2_t *x, uint32_t k);
void b2_clearbit(b2_t *x, uint32_t k);
void b2_flipbit(b2_t *x, uint32_t k);
bool b2_bitscan0(uint32_t *result, const b2_t *x, uint32_t start);
bool b2_bitscan1(uint32_t *result, const b2_t *x, uint32_t start);

void b2_bitand(b2_t *out, const b2_t *a, const b2_t *b);
void b2_bitor(b2_t *out, const b2_t *a, const b2_t *b);
void b2_bitxor(b2_t *out, const b2_t *a, const b2_t *b);
void b2_bitandnot(b2_t *out, const b2_t *a, const b2_t *b);
void b2_bitnot(b2_t *out, const b2_t *a, bool fixed_width,
               uint32_t width);

#endif
