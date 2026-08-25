#ifndef MPU_B9_H
#define MPU_B9_H

#include <stdlib.h>

#include "ptypes.h"

/* Little-endian base-10^N limbs, selected for the available accumulator. */
#if 0 && HAVE_UINT128
  /* This works, but benchmarking shows it is rarely faster when we have
   * good uint128_t support, and is much worse when it is _BitInt(128). */
  #define B9_DIGS  9
  #define B9_BASE  UINT32_C(1000000000)
  typedef uint32_t  b9limb_t;
  typedef uint128_t b9acc_t;
  #define B9_ACC_DEC_DIGS  39
#elif HAVE_UINT64
  #define B9_DIGS  6
  #define B9_BASE  UINT32_C(1000000)
  typedef uint32_t b9limb_t;
  typedef uint64_t b9acc_t;
  #define B9_ACC_DEC_DIGS  20
#else
  #define B9_DIGS  4
  #define B9_BASE  UINT32_C(10000)
  typedef uint32_t b9limb_t;
  typedef uint32_t b9acc_t;
  #define B9_ACC_DEC_DIGS  10
#endif

#define B9_NLIMBS(len)  (((len) + B9_DIGS - 1) / B9_DIGS)
#define B9_UV_DEC_DIGS  ((BITS_PER_WORD == 64) ? 20 : 10)
#define B9ACC_MAX        ((b9acc_t)~(b9acc_t)0)
#define B9_INLINE_LIMBS  (2 * sizeof(UV) / sizeof(b9limb_t))

typedef struct {
  b9limb_t *d;          /* d_small when inline, otherwise heap storage */
  UV        d_small[2]; /* two UVs of inline limb storage */
  uint32_t  alloc;
  uint32_t  n;          /* significant limbs; zero means value is zero */
  int       neg;
} b9_t;

MAYBE_UNUSED static INLINE void b9_init(b9_t *x)
{
  x->d = (b9limb_t*)x->d_small;
  x->alloc = B9_INLINE_LIMBS;
  x->n = 0;
  x->neg = 0;
}

MAYBE_UNUSED static INLINE bool b9_is_zero(const b9_t *x)
{
  return x->n == 0;
}

MAYBE_UNUSED static INLINE bool b9_is_one(const b9_t *x)
{
  return !x->neg && x->n == 1 && x->d[0] == 1;
}

MAYBE_UNUSED static INLINE bool b9_is_negative(const b9_t *x)
{
  return x->neg != 0;
}

MAYBE_UNUSED static INLINE bool b9_is_even(const b9_t *x)
{
  return x->n == 0 || !(x->d[0] & 1);
}

/* Replace x with |x| in place. */
MAYBE_UNUSED static INLINE void b9_abs(b9_t *x)
{
  x->neg = 0;
}

/* Set the sign while preserving canonical zero. */
MAYBE_UNUSED static INLINE void b9_set_negative(b9_t *x, bool neg)
{
  x->neg = (x->n != 0 && neg);
}

/* Return the canonical residue modulo 2, 4, 8, or 16. */
MAYBE_UNUSED static INLINE uint32_t b9_mod_power2(const b9_t *x,
                                                  uint32_t modulus)
{
  uint32_t r;
  MPUassert(modulus == 2 || modulus == 4 ||
            modulus == 8 || modulus == 16,
            "b9_mod_power2 bad modulus");
  if (x->n == 0) return 0;
  r = x->d[0] & (modulus - 1);
  return (x->neg && r != 0) ? modulus - r : r;
}

void b9_free(b9_t *x);
void b9_neg(b9_t *x);

void b9_ensure(b9_t *x, uint32_t need);
void b9_swap(b9_t *a, b9_t *b);
void b9_move(b9_t *dst, b9_t *src);

STRLEN b9_length(const b9_t *x);
void b9_set_str(b9_t *x, const char *s, STRLEN len);
STRLEN b9_get_str(char *buf, const b9_t *x);
void b9_set_uv(b9_t *x, UV v);
UV b9_get_uv(const b9_t *x);

void b9_set(b9_t *dst, const b9_t *src);
void b9_init_set(b9_t *x, const b9_t *y);
void b9_init_set_uv(b9_t *x, UV v);
void b9_init_set_str(b9_t *x, const char *s, STRLEN len);
void b9_init_set_pow2(b9_t *x, UV k);

int b9_cmp(const b9_t *a, const b9_t *b);
void b9_add(b9_t *out, const b9_t *a, const b9_t *b);
void b9_add_u32(b9_t *out, const b9_t *a, uint32_t v);
void b9_add_uv(b9_t *out, const b9_t *a, UV v);
void b9_mul(b9_t *out, const b9_t *a, const b9_t *b);
void b9_mul_u32(b9_t *out, const b9_t *a, uint32_t v);
void b9_mul_2exp(b9_t *out, const b9_t *a, UV bits);
void b9_pow(b9_t *out, const b9_t *a, UV exp);
void b9_fdivrem(b9_t *q, b9_t *r, const b9_t *a, const b9_t *b);
uint32_t b9_mod_u32(const b9_t *a, uint32_t p);
void b9_divexact_u32(b9_t *a, uint32_t p);

/* B9_BASE is divisible by 2 and 16 in every representation tier. */
MAYBE_UNUSED static INLINE void b9_tdiv2(b9_t *a)
{
  uint32_t carry = 0;
  uint32_t i;
  for (i = a->n; i-- > 0; ) {
    uint32_t limb = a->d[i];
    a->d[i] = (b9limb_t)((limb >> 1) + carry * (B9_BASE >> 1));
    carry = limb & 1;
  }
  while (a->n > 0 && a->d[a->n-1] == 0) a->n--;
  if (a->n == 0) a->neg = 0;
}

MAYBE_UNUSED static INLINE void b9_tdiv16(b9_t *a)
{
  uint32_t carry = 0;
  uint32_t i;
  for (i = a->n; i-- > 0; ) {
    uint32_t limb = a->d[i];
    a->d[i] = (b9limb_t)((limb >> 4) + carry * (B9_BASE >> 4));
    carry = limb & 15;
  }
  while (a->n > 0 && a->d[a->n-1] == 0) a->n--;
  if (a->n == 0) a->neg = 0;
}

uint32_t b9_divrem_u32_inplace(b9_t *a, uint32_t p);
int b9_cmp_abs(const b9_t *a, const b9_t *b);
void b9_addmod(b9_t *out, b9_t *tmp, const b9_t *a, const b9_t *b,
               const b9_t *m);
void b9_submod(b9_t *out, b9_t *tmp, const b9_t *a, const b9_t *b,
               const b9_t *m);
void b9_mulmod(b9_t *out, b9_t *tmp, const b9_t *a, const b9_t *b,
               const b9_t *m);
void b9_powmod(b9_t *out, b9_t *tmp, const b9_t *a, const b9_t *e,
               const b9_t *m);
void b9_product(b9_t A[], size_t a, size_t b);
void b9_product_u32(b9_t *out, const uint32_t A[], size_t len);

#endif
