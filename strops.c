#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "constants.h"
#include "mulmod.h"
#include "strops.h"

/******************************************************************************/
/* Helper functions. */
/******************************************************************************/

/* Write UV v as a decimal string; return length.  buf must have >= 21 bytes. */
static STRLEN uv_to_str(char* buf, UV v)
{
  STRLEN slen = 0;
  if (v == 0) { buf[0] = '0'; return 1; }
  while (v > 0) { buf[slen++] = '0' + (char)(v % 10); v /= 10; }
  if (slen > 1) {
    char *L = buf, *R = buf+slen;
    while (--R > L) { char t = *R; *R = *L; *L++ = t; }
  }
  return slen;
}

static UV gcduv(UV x, UV y) {
  while (y) { UV t = y; y = x % y; x = t; }
  return x;
}

/* Parse a canonical decimal string (no sign, no leading zeros) to UV.
 * Returns 1 if it fits in UV, 0 on overflow. */
static int str_to_uv_s(const char* s, STRLEN slen, UV* out) {
  UV v = 0;
  STRLEN i;
  for (i = 0; i < slen; i++) {
    UV d = (UV)(s[i] - '0');
    if (v > (UV_MAX - d) / 10) return 0;
    v = v * 10 + d;
  }
  *out = v;
  return 1;
}


/******************************************************************************/
/*               BASE-B9 INTERNAL REPRESENTATION FOR FAST ARITHMETIC          */
/*                                                                            */
/* Numbers are stored internally as little-endian arrays of b9limb_t, each    */
/* holding B9_DIGS decimal digits (one "digit" in base B9_BASE).  The         */
/* accumulator type b9acc_t is wide enough that the schoolbook inner-product  */
/* loop never overflows for any practical input length.                       */
/*                                                                            */
/* Tier selection is based on compile-time capability flags from ptypes.h:    */
/*   HAVE_UINT128: base 10^9, 128-bit accumulator (no flush needed)           */
/*   HAVE_UINT64:  base 10^6,  64-bit accumulator (no flush needed in practice) */
/*   else:         base 10^4,  32-bit accumulator (per-row carry flush needed)*/
/******************************************************************************/

#if 0 && HAVE_UINT128
  /* This works, but benchmarking shows it is rarely faster when we have
   * good uint128_t support, and is much worse when it is _BitInt(128).
   * Disable this path for now. */
  #define B9_DIGS  9
  #define B9_BASE  UINT32_C(1000000000)
  typedef uint32_t  b9limb_t;
  typedef uint128_t b9acc_t;
  /* Decimal digits needed for any b9acc_t value. */
  #define B9_ACC_DEC_DIGS  39
#elif HAVE_UINT64
  #define B9_DIGS  6
  #define B9_BASE  UINT32_C(1000000)
  typedef uint32_t  b9limb_t;
  typedef uint64_t  b9acc_t;
  #define B9_ACC_DEC_DIGS  20
#else
  #define B9_DIGS  4
  #define B9_BASE  UINT32_C(10000)
  typedef uint32_t  b9limb_t;
  typedef uint32_t  b9acc_t;
  #define B9_ACC_DEC_DIGS  10
#endif

/* Limbs needed for a decimal string of len digits. */
#define B9_NLIMBS(len)  (((len) + B9_DIGS - 1) / B9_DIGS)

/* Decimal digits needed for any UV value. */
#define B9_UV_DEC_DIGS  ((BITS_PER_WORD == 64) ? 20 : 10)

#define B9ACC_MAX        ((b9acc_t)~(b9acc_t)0)

/* Inline limb slots: 2 UVs worth, so 4 limbs on 64-bit (36 digits)
 * or 2 limbs on 32-bit with uint64_t (12 digits).
 * Avoids malloc for many real-world inputs. */
#define B9_INLINE_LIMBS  (2 * sizeof(UV) / sizeof(b9limb_t))

/******************************************************************************/
/*                       B9 ARITHMETIC  (GMP-STYLE INTERFACE)                 */
/*                                                                            */
/* b9_t holds a signed big integer as a managed little-endian limb array.     */
/* d points to d_small when the value fits in B9_INLINE_LIMBS limbs,          */
/* otherwise to a heap allocation.  All b9_* functions are static.            */
/******************************************************************************/

typedef struct {
    b9limb_t *d;         /* (b9limb_t*)d_small when inline, else heap */
    UV        d_small[2];/* inline storage: 2 UVs = B9_INLINE_LIMBS limbs */
    uint32_t  alloc;     /* allocated limb slots */
    uint32_t  n;         /* significant limbs (0 means value is zero) */
    int       neg;       /* 1 if negative, 0 otherwise */
} b9_t;

static void b9_init(b9_t *x)
  { x->d = (b9limb_t*)x->d_small;  x->alloc = B9_INLINE_LIMBS;  x->n = 0;  x->neg = 0; }
/* b9_init_set, b9_init_set_uv, b9_init_set_str are defined later */

static void b9_free(b9_t *x)
{
  if (x->d != (b9limb_t*)x->d_small) free(x->d);
  x->d = (b9limb_t*)x->d_small;  x->alloc = B9_INLINE_LIMBS;  x->n = 0;  x->neg = 0;
}

static void b9_ensure(b9_t *x, uint32_t need)
{
  if (x->alloc >= need) return;
  if (x->d == (b9limb_t*)x->d_small) {
    x->d = (b9limb_t *) malloc((size_t)need * sizeof(b9limb_t));
    if (x->n) memcpy(x->d, x->d_small, x->n * sizeof(b9limb_t));
  } else {
    x->d = (b9limb_t *) realloc(x->d, (size_t)need * sizeof(b9limb_t));
  }
  x->alloc = need;
}

/* Exchange two b9_t values, correctly handling inline storage. */
static void b9_swap(b9_t *a, b9_t *b)
{
  int a_inline = (a->d == (b9limb_t*)a->d_small);
  int b_inline = (b->d == (b9limb_t*)b->d_small);
  UV  a0 = a->d_small[0],  a1 = a->d_small[1];
  UV  b0 = b->d_small[0],  b1 = b->d_small[1];
  b9_t tmp = *a;  *a = *b;  *b = tmp;
  if (b_inline) { a->d = (b9limb_t*)a->d_small;  a->d_small[0] = b0;  a->d_small[1] = b1; }
  if (a_inline) { b->d = (b9limb_t*)b->d_small;  b->d_small[0] = a0;  b->d_small[1] = a1; }
}

/* Free dst, move src into dst.  src is blanked. */
static void b9_move(b9_t *dst, b9_t *src)
{
  b9_free(dst);
  *dst = *src;
  if (src->d == (b9limb_t*)src->d_small)
      dst->d = (b9limb_t*)dst->d_small;
  b9_init(src);
}

static void b9_neg(b9_t *x)
  { x->neg = x->n == 0 ? 0 : !x->neg; }

/* Exact length in decimal digits including possible sign character */
static STRLEN b9_length(const b9_t *x)
{
  STRLEN digits;
  uint32_t v;
  if (x->n == 0)
    return 1;
  digits = (STRLEN)(x->n - 1) * B9_DIGS;
  v      = x->d[x->n - 1];
  do { digits++; v /= 10; } while (v > 0);
  return digits + (x->neg ? 1 : 0);
}

/* Set x from a signed decimal string (may have leading sign and/or zeros). */
static void b9_set_str(b9_t *x, const char *s, STRLEN len)
{
  int neg;
  uint32_t nlimbs;
  STRLEN pos;

  neg = (len > 0 && s[0] == '-');
  if (len > 0 && (s[0] == '-' || s[0] == '+')) { s++;  len--; }
  while (len > 1 && s[0] == '0') { s++;  len--; }

  if (len == 0 || (len == 1 && s[0] == '0')) { x->n = 0;  x->neg = 0;  return; }

  nlimbs = (uint32_t)B9_NLIMBS(len);
  b9_ensure(x, nlimbs);
  x->n = 0;
  pos = len;
  while (pos > 0) {
    STRLEN start = (pos > (STRLEN)B9_DIGS) ? pos - B9_DIGS : 0;
    STRLEN i;
    uint32_t v = 0;
    for (i = start; i < pos; i++)
      v = v * 10 + (uint32_t)(s[i] - '0');
    x->d[x->n++] = (b9limb_t)v;
    pos = start;
  }
  x->neg = neg;
}

/* Write x to buf as a signed decimal string (no NUL).
 * buf must have at least x->n * B9_DIGS + 2 bytes.
 * Returns string length. */
static STRLEN b9_get_str(char *buf, const b9_t *x)
{
  STRLEN pos = 0;
  uint32_t i;

  if (x->n == 0) { buf[0] = '0';  return 1; }
  if (x->neg)    buf[pos++] = '-';

  { /* Most-significant limb: no leading zeros */
    char tmp[12];
    uint32_t v = x->d[x->n - 1];
    STRLEN j = (STRLEN)B9_DIGS;
    do { tmp[--j] = '0' + (char)(v % 10); v /= 10; } while (v > 0);
    while (j < (STRLEN)B9_DIGS) buf[pos++] = tmp[j++];
  }
  /* Remaining limbs: exactly B9_DIGS digits, zero-padded */
  for (i = x->n - 1; i-- > 0; ) {
    char tmp[12];
    uint32_t v = x->d[i];
    STRLEN j = (STRLEN)B9_DIGS;
    while (j > 0) { tmp[--j] = '0' + (char)(v % 10); v /= 10; }
    memcpy(buf + pos, tmp, (size_t)B9_DIGS);
    pos += (STRLEN)B9_DIGS;
  }
  return pos;
}

/* Set x from an unsigned machine word. */
static void b9_set_uv(b9_t *x, UV v)
{
  uint32_t n = 0;
  b9_ensure(x, (uint32_t)B9_NLIMBS(B9_UV_DEC_DIGS));
  x->neg = 0;
  if (v == 0) { x->n = 0;  return; }
  while (v > 0) { x->d[n++] = (b9limb_t)((UV)(v % B9_BASE)); v /= (UV)B9_BASE; }
  x->n = n;
}

/* Set x from an unsigned accumulator value. */
static void b9_set_acc(b9_t *x, b9acc_t v)
{
  uint32_t n = 0;
  b9_ensure(x, (uint32_t)B9_NLIMBS(B9_ACC_DEC_DIGS));
  x->neg = 0;
  if (v == 0) { x->n = 0;  return; }
  while (v > 0) { x->d[n++] = (b9limb_t)(v % B9_BASE); v /= B9_BASE; }
  x->n = n;
}

/* Convert b9 to u32, no size check. */
MAYBE_UNUSED static uint32_t b9_get_u32(const b9_t *x)
{
  uint32_t i, v;
  for (v = 0, i = x->n; i-- > 0; )
    v = v * B9_BASE + (uint32_t)x->d[i];
  return v;
}

/* Convert b9 to UV, no size check. [n * B9_DIGS < sizeof(UV)*3] */
static UV b9_get_uv(const b9_t* x) {
  UV v = 0;
  uint32_t i;
  for (i = x->n; i-- > 0; )
    v = v * (UV)B9_BASE + (UV)x->d[i];
  return v;
}

/* Simpler helpers */
static void b9_init_set_uv(b9_t *x, UV v)
  { b9_init(x);  b9_set_uv(x,v); }
static void b9_init_set_acc(b9_t *x, b9acc_t v)
  { b9_init(x);  b9_set_acc(x,v); }
static void b9_init_set_str(b9_t *x, const char *s, STRLEN len)
  { b9_init(x);  b9_set_str(x,s,len); }
static void b9_init_set(b9_t *x, const b9_t *y)
{
  b9_init(x);
  b9_ensure(x, y->n);
  if (y->n) memcpy(x->d, y->d, y->n * sizeof(b9limb_t));
  x->n = y->n;
  x->neg = y->neg;
}


/* Signed comparison.  Returns -1, 0, or 1. */
static int b9_cmp(const b9_t *a, const b9_t *b)
{
  uint32_t i;
  int gt;
  if (a->n == 0 && b->n == 0) return 0;
  if (a->neg != b->neg) return a->neg ? -1 : 1;
  if (a->n != b->n) {
    gt = (a->n > b->n) ? 1 : -1;
    return a->neg ? -gt : gt;
  }
  for (i = a->n; i-- > 0; ) {
    if (a->d[i] != b->d[i]) {
      gt = (a->d[i] > b->d[i]) ? 1 : -1;
      return a->neg ? -gt : gt;
    }
  }
  return 0;
}

/*--- Unsigned-magnitude helpers (raw limb arrays, no alloc) ---*/

static int b9mag_cmp(const b9limb_t *a, uint32_t na,
                     const b9limb_t *b, uint32_t nb)
{
  uint32_t i;
  if (na != nb) return (na > nb) ? 1 : -1;
  for (i = na; i-- > 0; )
    if (a[i] != b[i]) return (a[i] > b[i]) ? 1 : -1;
  return 0;
}

/* out = |a| + |b|.  out needs max(na,nb)+1 limbs.  Returns limb count. */
static uint32_t b9mag_add(b9limb_t *out,
                          const b9limb_t *a, uint32_t na,
                          const b9limb_t *b, uint32_t nb)
{
  uint32_t i, rn;
  b9acc_t carry = 0;
  if (na < nb) {
    const b9limb_t *t = a;  uint32_t tn = na;
    a = b;  na = nb;  b = t;  nb = tn;
  }
  rn = na;
  for (i = 0; i < nb; i++) {
    b9acc_t s = (b9acc_t)a[i] + b[i] + carry;
    out[i] = (b9limb_t)(s % B9_BASE);
    carry  = s / B9_BASE;
  }
  for (; i < na; i++) {
    b9acc_t s = (b9acc_t)a[i] + carry;
    out[i] = (b9limb_t)(s % B9_BASE);
    carry  = s / B9_BASE;
  }
  if (carry) out[rn++] = (b9limb_t)carry;
  return rn;
}

/* out = |a| - |b|  (|a| >= |b| required).
 * out needs na limbs.  Returns significant limb count. */
static uint32_t b9mag_sub(b9limb_t *out,
                          const b9limb_t *a, uint32_t na,
                          const b9limb_t *b, uint32_t nb)
{
  uint32_t i, rn;
  int borrow = 0;
  for (i = 0; i < nb; i++) {
    b9acc_t bi_b = (b9acc_t)b[i] + borrow;
    if ((b9acc_t)a[i] >= bi_b) {
      out[i] = (b9limb_t)((b9acc_t)a[i] - bi_b);
      borrow = 0;
    } else {
      out[i] = (b9limb_t)(B9_BASE + (b9acc_t)a[i] - bi_b);
      borrow = 1;
    }
  }
  for (; i < na; i++) {
    if (!borrow)   { out[i] = a[i]; }
    else if (a[i]) { out[i] = a[i] - 1;  borrow = 0; }
    else           { out[i] = (b9limb_t)(B9_BASE - 1); }
  }
  rn = na;
  while (rn > 1 && out[rn-1] == 0) rn--;
  return rn;
}

/* out = a + b (signed).  out may alias a or b. */
static void b9_add(b9_t *out, const b9_t *a, const b9_t *b)
{
  int cmp, use_tmp;
  b9limb_t *dst;
  uint32_t need, rn;

  if (a->n == 0) {
    if (out != b) {
      b9_ensure(out, b->n);
      if (b->n) memcpy(out->d, b->d, b->n * sizeof(b9limb_t));
      out->n = b->n;  out->neg = b->neg;
    }
    return;
  }
  if (b->n == 0) {
    if (out != a) {
      b9_ensure(out, a->n);
      if (a->n) memcpy(out->d, a->d, a->n * sizeof(b9limb_t));
      out->n = a->n;  out->neg = a->neg;
    }
    return;
  }

  use_tmp = (out->d == a->d || out->d == b->d);
  if (a->neg == b->neg) {
    need = (a->n > b->n ? a->n : b->n) + 1;
    if (use_tmp) {
      dst = (b9limb_t*) malloc(need * sizeof(b9limb_t));
    } else {
      b9_ensure(out, need);
      dst = out->d;
    }
    rn       = b9mag_add(dst, a->d, a->n, b->d, b->n);
    out->neg = a->neg;
  } else {
    cmp = b9mag_cmp(a->d, a->n, b->d, b->n);
    if (cmp == 0) { out->n = 0;  out->neg = 0;  return; }
    need = (cmp > 0) ? a->n : b->n;
    if (use_tmp) {
      dst = (b9limb_t*) malloc(need * sizeof(b9limb_t));
    } else {
      b9_ensure(out, need);
      dst = out->d;
    }
    if (cmp > 0) {
      rn = b9mag_sub(dst, a->d, a->n, b->d, b->n);
      out->neg = a->neg;
    } else {
      rn = b9mag_sub(dst, b->d, b->n, a->d, a->n);
      out->neg = b->neg;
    }
  }
  if (use_tmp) {
    b9_ensure(out, rn);
    memcpy(out->d, dst, rn * sizeof(b9limb_t));
    free(dst);
  }
  out->n = rn;
  if (out->n == 0) out->neg = 0;
}

#if 0
/* out = a - b (signed).  out may alias a or b. */
static void b9_sub(b9_t *out, const b9_t *a, const b9_t *b)
{
  b9_t nb = *b;   /* shallow copy: nb.d aliases b->d, only neg changes */
  b9_neg(&nb);
  b9_add(out, a, &nb);
}
#endif

static void b9_add_u32(b9_t *out, const b9_t *a, uint32_t v)
{
  b9acc_t carry;
  uint32_t i, need;

  if (v == 0) {
    if (out != a) {
      b9_ensure(out, a->n);
      if (a->n) memcpy(out->d, a->d, a->n * sizeof(b9limb_t));
      out->n = a->n;  out->neg = a->neg;
    }
    return;
  }

  if (a->neg) {
    b9_t b;
    b9_init_set_uv(&b, v);
    b9_add(out, a, &b);
    b9_free(&b);
    return;
  }

  need = a->n + (uint32_t)B9_NLIMBS(B9_UV_DEC_DIGS);
  b9_ensure(out, need);
  if (out != a && a->n)
    memcpy(out->d, a->d, a->n * sizeof(b9limb_t));

  carry = v;
  i = 0;
  while (carry) {
    b9acc_t s = carry;
    if (i < a->n) s += out->d[i];
    out->d[i++] = (b9limb_t)(s % B9_BASE);
    carry = s / B9_BASE;
  }
  out->n = (a->n > i) ? a->n : i;
  out->neg = 0;
}

static void b9_add_uv(b9_t *out, const b9_t *a, UV v)
{
  if (v <= (UV)UINT32_MAX) {
    b9_add_u32(out, a, (uint32_t)v);
  } else {
    b9_t b;
    b9_init_set_uv(&b, v);
    b9_add(out, a, &b);
    b9_free(&b);
  }
}

/* out = a * b (signed).
 * Reads all of a->d and b->d before writing to out->d, so out may alias
 * a or b (b9_ensure is called only after the multiply loop). */

/* Small multiplications use a stack buffer to avoid malloc overhead.
 * 64 limbs covers 576 digits (base 10^9), 384 (10^6), 256 (10^4). */
#define B9_MUL_STACK_LIMBS 64

static void b9_mul(b9_t *out, const b9_t *a, const b9_t *b)
{
  uint32_t i, j, rn;
  b9acc_t  stack_acc[B9_MUL_STACK_LIMBS + 1];
  b9acc_t *acc;
  int neg;

  if (a->n == 0 || b->n == 0) { out->n = 0;  out->neg = 0;  return; }

  neg = (a->neg != b->neg) ? 1 : 0;
  rn  = a->n + b->n;

  if (rn <= B9_MUL_STACK_LIMBS) {
    acc = stack_acc;
    memset(acc, 0, (rn + 1) * sizeof(b9acc_t));
  } else {
    acc = (b9acc_t*) calloc((size_t)(rn + 1), sizeof(b9acc_t));
  }

  for (i = 0; i < a->n; i++) {
    if (a->d[i] == 0) continue;
    for (j = 0; j < b->n; j++)
      acc[i + j] += (b9acc_t)a->d[i] * b->d[j];
    /* Per-row carry flush: only needed for base 10^4 (32-bit uint32_t acc).
     * With uint64_t (base 10^6) overflow needs 18M products/position. */
#if !HAVE_UINT128 && !HAVE_UINT64
    for (j = i; j <= i + b->n && j < rn; j++) {
      acc[j + 1] += acc[j] / B9_BASE;
      acc[j]     %= B9_BASE;
    }
#endif
  }

  /* a->d and b->d fully consumed; safe to resize out even if aliased */
  b9_ensure(out, rn);
  for (i = 0; i < rn; i++) {
    acc[i + 1] += acc[i] / B9_BASE;
    acc[i]     %= B9_BASE;
    out->d[i]   = (b9limb_t)acc[i];
  }
  if (acc != stack_acc) free(acc);

  while (rn > 1 && out->d[rn-1] == 0) rn--;
  out->n   = rn;
  out->neg = (out->n == 1 && out->d[0] == 0) ? 0 : neg;
}

MAYBE_UNUSED static void b9_mul_uv(b9_t *out, const b9_t *a, UV v)
{
  b9_t b;
  b9_init_set_uv(&b, v);
  b9_mul(out, a, &b);
  b9_free(&b);
}

static void b9_mul_u32(b9_t *out, const b9_t *a, uint32_t v)
{
#if HAVE_UINT64
  uint64_t carry = 0;
  uint32_t i, need;
  int neg = a->neg;

  if (v == 0 || a->n == 0) {
    out->n = 0;
    out->neg = 0;
    return;
  }

  need = a->n + (uint32_t)B9_NLIMBS(10) + 1;
  b9_ensure(out, need);
  for (i = 0; i < a->n; i++) {
    uint64_t t = (uint64_t)a->d[i] * v + carry;
    out->d[i] = (b9limb_t)(t % B9_BASE);
    carry = t / B9_BASE;
  }
  while (carry > 0) {
    out->d[i++] = (b9limb_t)(carry % B9_BASE);
    carry /= B9_BASE;
  }
  out->n = i;
  out->neg = neg;
#else
  b9_mul_uv(out, a, (UV)v);
#endif
}

/* out = a^exp.  out can alias a. */
static void b9_pow(b9_t *out, const b9_t *a, UV exp)
{
  b9_t sq, tmp;
  b9_init_set(&sq, a);  /* sq is a copy of a */
  b9_set_uv(out, 1);
  b9_init(&tmp);
  while (exp > 0) {
    if (exp & 1) {
      b9_mul(&tmp, out, &sq);
      b9_move(out, &tmp);
    }
    exp >>= 1;
    if (exp > 0) {
      b9_mul(&tmp, &sq, &sq);
      b9_move(&sq, &tmp);
    }
  }
  b9_free(&sq);  b9_free(&tmp);
}

/* There are a few places that want a power of two.  This wraps it. */
static void b9_init_set_pow2(b9_t *x, UV k)
{
  if (k < BITS_PER_WORD) {
    b9_init_set_uv(x, UVCONST(1) << k);
  } else {
    b9_init_set_uv(x, 2);
    b9_pow(x, x, k);
  }
}

/* Signed floor division and remainder (floor convention: rem has sign of b).
 * q or r may be NULL (that output is skipped).
 * b must be non-zero.  q and r must not alias a or b. */
static void b9_fdivrem(b9_t *q, b9_t *r, const b9_t *a, const b9_t *b)
{
  int na = (int)a->n, nb = (int)b->n;
  int m, j, k, qn;
  int neg_q = (a->neg != b->neg) ? 1 : 0;
  int neg_r = b->neg;
  b9_t q_tmp, r_tmp;
  b9_t *qp, *rp;

  if (!q) { b9_init(&q_tmp);  qp = &q_tmp; } else qp = q;
  if (!r) { b9_init(&r_tmp);  rp = &r_tmp; } else rp = r;

  /* a = 0 */
  if (na == 0) { qp->n = 0;  qp->neg = 0;  rp->n = 0;  rp->neg = 0;  goto done; }

  m = na - nb;

  /* |a| < |b|: q = 0, r = |a| */
  if (m < 0) {
    qp->n = 0;  qp->neg = 0;
    b9_ensure(rp, (uint32_t)na);
    memcpy(rp->d, a->d, (size_t)na * sizeof(b9limb_t));
    rp->n = (uint32_t)na;  rp->neg = 0;
    goto adjust;
  }

  b9_ensure(qp, (uint32_t)(m + 1));
  b9_ensure(rp, (uint32_t)nb);

  /* Single-limb divisor */
  if (nb == 1) {
    b9acc_t rem = 0;
    for (k = na - 1; k >= 0; k--) {
      b9acc_t cur = rem * B9_BASE + a->d[k];
      qp->d[k] = (b9limb_t)(cur / b->d[0]);
      rem       = cur % b->d[0];
    }
    rp->d[0] = (b9limb_t)rem;
    rp->n = 1;  rp->neg = 0;
    qn = na;
    while (qn > 1 && qp->d[qn-1] == 0) qn--;
    qp->n = (uint32_t)qn;  qp->neg = 0;
    goto adjust;
  }

  /* Multi-limb: Knuth Algorithm D (TAOCP Vol.2 §4.3.1) */
  {
    b9limb_t *u, *v;
    uint32_t d;

    u = (b9limb_t*) malloc((size_t)(na + 1) * sizeof(b9limb_t));
    v = (b9limb_t*) malloc((size_t)nb        * sizeof(b9limb_t));

    d = (uint32_t)(B9_BASE / ((b9acc_t)b->d[nb-1] + 1));

    { b9acc_t carry = 0;
      for (k = 0; k < na; k++) {
        b9acc_t t = (b9acc_t)a->d[k] * d + carry;
        u[k] = (b9limb_t)(t % B9_BASE);  carry = t / B9_BASE;
      }
      u[na] = (b9limb_t)carry;
    }
    { b9acc_t carry = 0;
      for (k = 0; k < nb; k++) {
        b9acc_t t = (b9acc_t)b->d[k] * d + carry;
        v[k] = (b9limb_t)(t % B9_BASE);  carry = t / B9_BASE;
      }
    }

    for (j = m; j >= 0; j--) {
      b9acc_t q_hat, r_hat, borrow;
      { b9acc_t top2 = (b9acc_t)u[j+nb] * B9_BASE + u[j+nb-1];
        q_hat = top2 / v[nb-1];
        r_hat = top2 % v[nb-1];
      }
      while (q_hat >= B9_BASE ||
             q_hat * v[nb-2] > B9_BASE * r_hat + u[j+nb-2]) {
        q_hat--;  r_hat += v[nb-1];
        if (r_hat >= B9_BASE) break;
      }
      borrow = 0;
      for (k = 0; k < nb; k++) {
        b9acc_t prod = (b9acc_t)q_hat * v[k] + borrow;
        b9limb_t sub = (b9limb_t)(prod % B9_BASE);
        borrow = prod / B9_BASE;
        if (u[j+k] >= sub) { u[j+k] -= sub; }
        else { u[j+k] = (b9limb_t)(B9_BASE + u[j+k] - sub);  borrow++; }
      }
      qp->d[j] = (b9limb_t)q_hat;
      if (u[j+nb] < (b9limb_t)borrow) {
        b9acc_t carry = 0;
        qp->d[j]--;
        for (k = 0; k < nb; k++) {
          b9acc_t t = (b9acc_t)u[j+k] + v[k] + carry;
          u[j+k] = (b9limb_t)(t % B9_BASE);  carry = t / B9_BASE;
        }
        u[j+nb] = (b9limb_t)((b9acc_t)u[j+nb] + B9_BASE - borrow + carry);
        if (u[j+nb] >= (b9limb_t)B9_BASE) u[j+nb] -= (b9limb_t)B9_BASE;
      } else { u[j+nb] -= (b9limb_t)borrow; }
    }

    if (d == 1) {
      for (k = 0; k < nb; k++) rp->d[k] = u[k];
    } else {
      b9acc_t rem = 0;
      for (k = nb - 1; k >= 0; k--) {
        b9acc_t t = rem * B9_BASE + u[k];
        rp->d[k] = (b9limb_t)(t / d);  rem = t % d;
      }
    }
    rp->n = (uint32_t)nb;
    while (rp->n > 1 && rp->d[rp->n-1] == 0) rp->n--;
    rp->neg = 0;

    qn = m + 1;
    while (qn > 1 && qp->d[qn-1] == 0) qn--;
    qp->n = (uint32_t)qn;  qp->neg = 0;

    free(u);  free(v);
  }

adjust:
  /* Normalize zero representation */
  if (qp->n == 1 && qp->d[0] == 0) qp->n = 0;
  if (rp->n == 1 && rp->d[0] == 0) rp->n = 0;

  /* Floor adjustment: when signs differ and remainder != 0 */
  if (neg_q && rp->n != 0) {
    /* q_floor = q_trunc + 1 */
    { uint32_t ii;
      b9acc_t carry = 1;
      b9_ensure(qp, qp->n + 1);
      for (ii = 0; ii < qp->n && carry; ii++) {
        b9acc_t s = (b9acc_t)qp->d[ii] + carry;
        qp->d[ii] = (b9limb_t)(s % B9_BASE);
        carry = s / B9_BASE;
      }
      if (carry) { qp->d[qp->n] = (b9limb_t)carry;  qp->n++; }
    }
    /* r_floor = |b| - r_trunc */
    { b9limb_t *tmp = (b9limb_t*) malloc((size_t)nb * sizeof(b9limb_t));
      uint32_t rn_new = b9mag_sub(tmp, b->d, (uint32_t)nb, rp->d, rp->n);
      b9_ensure(rp, rn_new);
      memcpy(rp->d, tmp, rn_new * sizeof(b9limb_t));
      rp->n = rn_new;
      free(tmp);
      if (rp->n == 1 && rp->d[0] == 0) rp->n = 0;
    }
  }

  qp->neg = (qp->n != 0) ? neg_q : 0;
  rp->neg = (rp->n != 0) ? neg_r : 0;

done:
  if (!q) b9_free(&q_tmp);
  if (!r) b9_free(&r_tmp);
}

/* Compute b9 value a mod small UV p (read-only). */
/* Compute a mod p for a non-zero uint32_t p.
 *
 * With p < 2^32, rem < p < 2^32 always, so rem*B9_BASE + d[i] is at most
 * (2^32-1)*B9_BASE + B9_BASE-1 < 2^32 * B9_BASE.  Since B9_BASE <= 10^9 < 2^30,
 * this is under 2^62, fitting in uint64_t for all three B9 tiers.
 * On the rare platform with neither uint64_t nor uint128_t (B9_BASE=10^4),
 * fall back to b9_fdivrem. */
static uint32_t b9_mod_u32(const b9_t* a, uint32_t p) {
#if HAVE_UINT64
  uint32_t r = 0;
  uint32_t i;
  for (i = a->n; i-- > 0; )
    r = (uint32_t)(((uint64_t)r * B9_BASE + a->d[i]) % p);
  return r;
#else
  b9_t bp, bq, br;
  uint32_t result;
  b9_init_set_uv(&bp, (UV)p);
  b9_init(&bq);  b9_init(&br);
  b9_fdivrem(&bq, &br, a, &bp);
  result = b9_get_u32(&br);
  b9_free(&bp);  b9_free(&bq);  b9_free(&br);
  return result;
#endif
}

/* Divide b9 value a in-place by uint32_t p.  Caller must ensure p | a exactly.
 * Same overflow analysis as b9_mod_u32: uint64_t always suffices when available. */
static void b9_divexact_u32(b9_t* a, uint32_t p) {
#if HAVE_UINT64
  uint32_t rem = 0;
  uint32_t i;
  for (i = a->n; i-- > 0; ) {
    uint64_t cur = (uint64_t)rem * B9_BASE + a->d[i];
    a->d[i] = (b9limb_t)(cur / p);
    rem = (uint32_t)(cur % p);
  }
  while (a->n > 0 && a->d[a->n-1] == 0) a->n--;
#else
  b9_t bp, bq, br;
  b9_init_set_uv(&bp, (UV)p);
  b9_init(&bq);  b9_init(&br);
  b9_fdivrem(&bq, &br, a, &bp);
  b9_move(a, &bq);
  b9_free(&bp);  b9_free(&br);
#endif
}

/* Divide positive b9 value a in-place by 2, truncating toward zero. */
static void b9_tdiv2(b9_t* a) {
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

static void b9_tdiv16(b9_t* a) {
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

static uint32_t b9_divrem_u32_inplace(b9_t* a, uint32_t p)
{
#if HAVE_UINT64
  uint32_t rem = 0;
  uint32_t i;

  if (p == 2) {
    rem = (a->n == 0) ? 0 : (a->d[0] & 1);
    b9_tdiv2(a);
    return rem;
  }
  if (p == 16) {
    rem = (a->n == 0) ? 0 : (a->d[0] & 15);
    b9_tdiv16(a);
    return rem;
  }

  for (i = a->n; i-- > 0; ) {
    uint64_t cur = (uint64_t)rem * B9_BASE + a->d[i];
    a->d[i] = (b9limb_t)(cur / p);
    rem = (uint32_t)(cur % p);
  }
  while (a->n > 0 && a->d[a->n-1] == 0) a->n--;
  if (a->n == 0) a->neg = 0;
  return rem;
#else
  b9_t bp, bq, br;
  uint32_t rem;
  b9_init_set_uv(&bp, (UV)p);
  b9_init(&bq);  b9_init(&br);
  b9_fdivrem(&bq, &br, a, &bp);
  rem = b9_get_u32(&br);
  b9_move(a, &bq);
  b9_free(&bp);  b9_free(&br);
  return rem;
#endif
}


static int b9_cmp_abs(const b9_t *a, const b9_t *b)
{
  return b9mag_cmp(a->d, a->n, b->d, b->n);
}

static void b9_reduce_mod(b9_t *out, b9_t *tmp, const b9_t *a, const b9_t *m)
{
  b9_fdivrem(NULL, tmp, a, m);
  b9_move(out, tmp);
}

static void b9_addmod(b9_t *out, b9_t *tmp, const b9_t *a, const b9_t *b, const b9_t *m)
{
  b9_t ar, br;
  const b9_t *ap = a, *bp = b;

  b9_init(&ar);  b9_init(&br);
  if (a->neg || b9_cmp_abs(a, m) >= 0) {
    b9_reduce_mod(&ar, tmp, a, m);
    ap = &ar;
  }
  if (b->neg || b9_cmp_abs(b, m) >= 0) {
    b9_reduce_mod(&br, tmp, b, m);
    bp = &br;
  }

  b9_add(out, ap, bp);
  b9_reduce_mod(out, tmp, out, m);
  b9_free(&ar);  b9_free(&br);
}

static void b9_submod(b9_t *out, b9_t *tmp, const b9_t *a, const b9_t *b, const b9_t *m)
{
  b9_t nb = *b;   /* shallow copy: nb.d aliases b->d, only neg changes */
  b9_neg(&nb);
  b9_addmod(out, tmp, a, &nb, m);
}

static void b9_mulmod_reduced(b9_t *out, b9_t *tmp, const b9_t *a, const b9_t *b, const b9_t *m)
{
  b9_mul(out, a, b);
  b9_reduce_mod(out, tmp, out, m);
}

static void b9_mulmod(b9_t *out, b9_t *tmp, const b9_t *a, const b9_t *b, const b9_t *m)
{
  b9_t ar, br;
  const b9_t *ap = a, *bp = b;

  b9_init(&ar);  b9_init(&br);
  if (a->neg || b9_cmp_abs(a, m) >= 0) {
    b9_reduce_mod(&ar, tmp, a, m);
    ap = &ar;
  }
  if (b->neg || b9_cmp_abs(b, m) >= 0) {
    b9_reduce_mod(&br, tmp, b, m);
    bp = &br;
  }

  b9_mulmod_reduced(out, tmp, ap, bp, m);
  b9_free(&ar);  b9_free(&br);
}

static void b9_powmod(b9_t *out, b9_t *tmp, const b9_t *a, const b9_t *e, const b9_t *m)
{
  b9_t base, exp;

  b9_init_set(&base, a);
  b9_init_set(&exp, e);
  b9_set_uv(out, 1);

  if (base.neg || b9_cmp_abs(&base, m) >= 0)
    b9_reduce_mod(&base, tmp, &base, m);

  while (exp.n != 0) {
    if (exp.d[0] & 1)
      b9_mulmod_reduced(out, tmp, out, &base, m);
    b9_tdiv2(&exp);
    if (exp.n != 0)
      b9_mulmod_reduced(&base, tmp, &base, &base, m);
  }

  b9_free(&base);  b9_free(&exp);
}

static void b9_product(b9_t A[], size_t a, size_t b) {
  if (b <= a) {
    /* A[a] already correct */
  } else if (b == a+1) {
    b9_mul(&A[a], &A[a], &A[b]);
  } else if (b == a+2) {
    b9_mul(&A[a+1], &A[a+1], &A[a+2]);
    b9_mul(&A[a], &A[a], &A[a+1]);
  } else {
    size_t c = a + (b-a+1)/2;
    b9_product(A, a, c-1);
    b9_product(A, c, b);
    b9_mul(&A[a], &A[a], &A[c]);
  }
}

static void b9_product_u32(b9_t *out, const uint32_t A[], size_t len)
{
  b9_t *B = 0;
  b9acc_t prod = 1;
  size_t i, max_chunks, nprod = 0;
  const size_t n_u32_per_acc = sizeof(b9acc_t) / sizeof(uint32_t);

  if (len == 0) {
    b9_set_uv(out, 1);
    return;
  }

  max_chunks = (len + n_u32_per_acc - 1) / n_u32_per_acc;
  if (max_chunks > (size_t)MAX_SIZET / sizeof(b9_t))
    croak("internal: b9_product_u32 cannot allocate");
  B = (b9_t*) malloc(max_chunks * sizeof(b9_t));
  if (B == 0)
    croak("internal: b9_product_u32 allocation failed");

  for (i = 0; i < len; i++) {
    if (A[i] == 0) { prod = 0; break; }
    if (A[i] == 1) continue;
    if (prod > B9ACC_MAX / A[i]) {
      b9_init_set_acc(&B[nprod++], prod);
      prod = A[i];
    } else {
      prod *= A[i];
    }
  }

  if (prod == 0) {
    b9_set_uv(out, 0);
  } else {
    if (prod != 1)
      b9_init_set_acc(&B[nprod++], prod);
    if (nprod == 0) {
      b9_set_uv(out, 1);
    } else {
      b9_product(B, 0, nprod-1);
      b9_move(out, &B[0]);
    }
  }

  for (i = 0; i < nprod; i++)
    b9_free(&B[i]);
  free(B);
}


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


/* Strip sign and leading zeros from an already-validated decimal string.
 * No character-by-character validation — caller guarantees valid input.
 * Returns 1 if negative, 0 otherwise. */
static int strint_strip(const char **sp, STRLEN *slen)
{
  const char *s = *sp;
  STRLEN len = *slen;
  int neg = (len > 0 && s[0] == '-');
  if (len > 0 && (s[0] == '-' || s[0] == '+')) { s++;  len--; }
  while (len > 1 && s[0] == '0') { s++;  len--; }
  if (neg && len == 1 && s[0] == '0') neg = 0;   /* -0 → 0 */
  *sp = s;
  *slen = len;
  return neg;
}

/* Parse a signed decimal integer string: strip optional sign and leading
 * zeros, validate that all remaining characters are digits, and update *sp
 * and *slen to point to the canonical digit sequence.  Returns 1 if negative,
 * 0 otherwise.  Croaks on invalid input. */
static int strint_strip_and_validate(const char **sp, STRLEN *slen)
{
  STRLEN i = 0;
  int neg = 0;
  if (*sp != NULL) {
    neg = strint_strip(sp, slen);
    for (i = 0; i < *slen; i++)
      if (!isDIGIT((*sp)[i]))
        break;
  }
  if (*sp == NULL || *slen == 0 || i < *slen)
    croak("Parameter must be an integer");
  return neg;
}

char* strint_reverse_digits(const char* n, STRLEN nlen, UV base, STRLEN* rlen)
{
  b9_t N, R;
  char *out;
  STRLEN outlen;

  if (rlen) *rlen = 0;
  if (base < 2 || base > (UV)UINT32_MAX)
    return 0;

  strint_strip(&n, &nlen);  /* Ignore sign, remove leading zeros. */

  if (nlen == 1 && n[0] == '0') {
    out = (char*) malloc(2);
    if (out == 0) return 0;
    out[0] = '0';
    out[1] = '\0';
    if (rlen) *rlen = 1;
    return out;
  }

  if (base == 10) {
    STRLEN i;
    while (nlen > 1 && n[nlen-1] == '0')
      nlen--;
    out = (char*) malloc((size_t)nlen + 1);
    if (out == 0) return 0;
    for (i = 0; i < nlen; i++)
      out[i] = n[nlen-i-1];
    out[nlen] = '\0';
    if (rlen) *rlen = nlen;
    return out;
  }

  b9_init_set_str(&N, n, nlen);
  b9_init_set_uv(&R, 0);
  while (N.n != 0) {
    uint32_t d = b9_divrem_u32_inplace(&N, (uint32_t)base);
    b9_mul_u32(&R, &R, (uint32_t)base);
    b9_add_u32(&R, &R, d);
  }

  out = (char*) malloc((size_t)b9_length(&R) + 1);
  outlen = 0;
  if (out != 0) {
    outlen = b9_get_str(out, &R);
    out[outlen] = '\0';
  }
  if (rlen) *rlen = outlen;
  b9_free(&N);
  b9_free(&R);
  return out;
}

int strint_cmp(const char* a, STRLEN alen, const char* b, STRLEN blen) {
  STRLEN i;
  int aneg = strint_strip(&a, &alen);
  int bneg = strint_strip(&b, &blen);
  if (aneg != bneg)  return (bneg) ? 1 : -1;
  if (aneg) { /* swap a and b if both negative */
    const char* t = a;  STRLEN tlen = alen;
    a = b; b = t;  alen = blen;  blen = tlen;
  }
  if (alen != blen)  return (alen > blen) ? 1 : -1;
  for (i = 0; i < blen; i++)
    if (a[i] != b[i])
      return  (a[i] > b[i]) ? 1 : -1;
  return 0;
}

/* 1. Perform signed integer validation on b/blen.
 * 2. Compare to a/alen using min or max based on first arg.
 * 3. Return 0 to select a, 1 to select b.
 */
bool strint_minmax(bool min, const char* a, STRLEN alen, const char* b, STRLEN blen)
{
  int aneg, bneg;
  STRLEN i;

  /* a is checked, process b */
  bneg = strint_strip_and_validate(&b, &blen);

  if (a == 0) return 1;

  aneg = (a[0] == '-');
  if (a[0] == '-' || a[0] == '+') { a++; alen--; }
  while (alen > 0 && *a == '0') { a++; alen--; }

  if (aneg != bneg)  return  min  ?  (bneg == 1)  :  (aneg == 1);
  if (aneg == 1)  min = !min;
  if (alen != blen)  return  min  ?  (alen > blen) :  (blen > alen);

  for (i = 0; i < blen; i++)
    if (a[i] != b[i])
      return  min  ?  (a[i] > b[i])  :  (b[i] > a[i]);
  return 0; /* equal */
}

/* Absolute value.  out needs alen bytes. */
STRLEN strint_abs(char* out, const char* a, STRLEN alen)
{
  if (alen > 0 && a[0] == '-') { a++;  alen--; }
  memcpy(out, a, alen);
  return alen;
}

/* Negation.  out needs alen+1 bytes. */
STRLEN strint_neg(char* out, const char* a, STRLEN alen)
{
  int neg = strint_strip(&a, &alen);
  if (!neg && !(alen == 1 && a[0] == '0')) {
    out[0] = '-';
    memcpy(out+1, a, alen);
    return alen+1;
  }
  memcpy(out, a, alen);
  return alen;
}

/* Add or subtract two signed decimal integer strings.
 * If negate_b is false, computes a+b; if true, computes a-b (i.e. a+(-b)).
 * Write result to out; caller must ensure at least max(alen,blen)+2 bytes
 * (the +2 is tight for subtraction when a is negative and blen >= alen).
 * out may alias a (in-place update is safe); out must not alias b.
 * Returns the result length (no NULL terminator written). */
STRLEN strint_add_s(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen, bool negate_b)
{
  b9_t A, B;
  STRLEN rlen;
  if (blen == 0 || (blen == 1 && b[0] == '0')) { memcpy(out, a, alen); return alen; }
  if (alen == 0 || (alen == 1 && a[0] == '0')) {
    if (!negate_b) { memcpy(out, b, blen); return blen; }
    /* 0 - b = -b */
    if (b[0] == '-') { memcpy(out, b + 1, blen - 1); return blen - 1; }
    out[0] = '-'; memcpy(out + 1, b, blen); return blen + 1;
  }
  /* Fast +1 / -1 for cases with no carry */
  if (blen == 1 && b[0] == '1') {
    int aneg = a[0] == '-', negb = !!negate_b;
    int same_sign = (aneg == negb);
    STRLEN digits = alen - (aneg ? 1 : 0);

    if ( (same_sign && a[alen-1] != '9') ||
         (!same_sign && digits > 1 && a[alen-1] != '0') ) {
      memcpy(out, a, alen-1);
      out[alen-1] = a[alen-1] + (same_sign ? 1 : -1);
      return alen;
    }
  }
  b9_init_set_str(&A, a, alen);
  b9_init_set_str(&B, b, blen);
  if (negate_b) b9_neg(&B);
  b9_add(&A, &A, &B);
  rlen = b9_get_str(out, &A);
  b9_free(&A);  b9_free(&B);
  return rlen;
}

char* strint_vecsum(const char* const* a, const STRLEN* alen, size_t n, STRLEN* rlen)
{
  b9_t sum, tmp;
  size_t i;
  STRLEN len;
  char *out;

  b9_init(&sum);
  b9_init(&tmp);
  for (i = 0; i < n; i++) {
    if (alen[i] == 0 || (alen[i] == 1 && a[i][0] == '0'))
      continue;
    b9_set_str(&tmp, a[i], alen[i]);
    b9_add(&sum, &sum, &tmp);
  }

  len = b9_length(&sum);
  out = (char*) malloc((size_t)len + 1);
  if (out != 0) {
    len = b9_get_str(out, &sum);
    out[len] = '\0';
  }
  if (rlen) *rlen = len;

  b9_free(&sum);
  b9_free(&tmp);
  return out;
}

char* strint_vecprod(const char* const* a, const STRLEN* alen, size_t n, STRLEN* rlen)
{
  b9_t prod, tmp, T[4], *A;
  size_t i, k, nrem, ibase, prodn;
  STRLEN outlen;
  char *out;


  b9_init(&tmp);
  b9_init_set_uv(&prod, 1);

  /* Small number of items to multiply */
  if (n < 4) {
    for (i = 0; i < n; i++) {
      b9_set_str(&tmp, a[i], alen[i]);
      b9_mul(&prod, &prod, &tmp);
    }
    goto return_prod;
  }

  /* Fast check for a zero */
  for (i = 0; i < n; i++)
    if (alen[i] == 1 && a[i][0] == '0')
      break;
  if (i < n) {
    b9_set_uv(&prod, 0);
    goto return_prod;
  }

  /* Process 1/4 of the items, product tree those, next quarter, etc.
   *
   * There is no performance benefit for the quarter split, but it saves
   * us high water memory as we allocate only 1/4 of the b9 objects.
   */

  prodn = (n+3)/4;
  A = (b9_t*) malloc(prodn * sizeof(b9_t));
  for (i = 0; i < prodn; i++)
    b9_init(&A[i]);
  for (i = 0; i < 4; i++)
    b9_init(&T[i]);

  nrem = n % 4;
  for (k = 0, ibase = 0; k < 4; k++) {
    size_t qlen = n/4 + (k+1 <= nrem);  /* length of this window */
    for (i = 0; i < qlen; i++)
      b9_set_str(&A[i], a[ibase+i], alen[ibase+i]);
    b9_product(A, 0, qlen-1);
    b9_move(&T[k], &A[0]);
    ibase += qlen;
  }

  for (i = 0; i < prodn; i++)
    b9_free(&A[i]);
  free(A);

  b9_mul(&T[0],&T[0],&T[1]);
  b9_mul(&T[2],&T[2],&T[3]);
  b9_mul(&prod,&T[0],&T[2]);

  for (i = 0; i < 4; i++)
    b9_free(&T[i]);

return_prod:
  out = (char*) malloc((size_t)b9_length(&prod) + 1);
  outlen = 0;
  if (out != 0) {
    outlen = b9_get_str(out, &prod);
    out[outlen] = '\0';
  }
  if (rlen) *rlen = outlen;

  b9_free(&prod);
  b9_free(&tmp);
  return out;
}

/******************************************************************************/
/*                                 BINOMIAL                                   */
/******************************************************************************/


extern uint32_t range_prime_sieve_32(uint32_t** list, uint32_t n, uint32_t os);
extern void free_prime_sieve_32(uint32_t* list);

static uint32_t binomial_prime_exp_32(uint32_t n, uint32_t k, uint32_t p)
{
  uint32_t e = 0, nk = n - k;

  while (n > 0) {
    n  /= p;
    k  /= p;
    nk /= p;
    e  += n - k - nk;
  }
  return e;
}

static bool b9_binomial_window_u32(b9_t *out, uint32_t n, uint32_t k)
{
  uint32_t i, nstart, pidx, nprimes, *primes = 0, *A = 0;

  if (k > n || k > n-k || k == 0)
    croak("internal: b9_binomial_window_u32 bad k");

  if ((size_t)k > (size_t)MAX_SIZET / sizeof(uint32_t))
    return 0;
  A = (uint32_t *) malloc((size_t) k * sizeof(uint32_t));
  if (A == 0)
    return 0;

  if (k <= primes_small[NPRIMES_SMALL-1]) {
    nprimes = NPRIMES_SMALL-1;
  } else {
    nprimes = range_prime_sieve_32(&primes, k, 0);
    if (nprimes == 0 || primes == 0) {
      free(A);
      return 0;
    }
  }

  nstart = n-k+1;
  for (i = 0; i < k; i++)
    A[i] = nstart + i;

  for (pidx = 0; pidx < nprimes; pidx++) {
    uint32_t ppow, p = primes ? primes[pidx] : primes_small[pidx+1];
    if (p > k) break;
    for (ppow = p; ppow <= k; ) {
      uint32_t c, count = k / ppow;
      uint32_t rem = nstart % ppow;
      uint32_t off = (rem == 0) ? 0 : ppow - rem;
      for (c = 0; c < count; c++, off += ppow)
        A[off] /= p;
      if (ppow > k / p) break;
      ppow *= p;
    }
  }
  b9_product_u32(out, A, k);

  if (primes != 0)
    free_prime_sieve_32(primes);
  free(A);
  return 1;
}


char* strint_binomial_u32_u32(uint32_t n, uint32_t k, STRLEN* rlen)
{
  b9_t result;
  b9_t *A = 0;
  char *out;
  STRLEN outlen;
  uint32_t pidx, nprimes, nprod = 0, *primes = 0;

  if (rlen) *rlen = 0;
  b9_init_set_uv(&result, 0);
  if (k > n)
    goto return_result;
  if (k > n-k)
    k = n-k;
  b9_set_uv(&result, 1);
  if (k == 0)
    goto return_result;

  if (k <= 65536 || k < n/20) {
    if (!b9_binomial_window_u32(&result, n, k))
      goto return_null;
    goto return_result;
  }

  if (n <= primes_small[NPRIMES_SMALL-1]) {
    nprimes = NPRIMES_SMALL-1;
  } else {
    nprimes = range_prime_sieve_32(&primes, n, 0);
    if (nprimes == 0 || primes == 0)
      goto return_null;
  }

  if ((size_t)nprimes > (size_t)MAX_SIZET / sizeof(b9_t))
    goto return_null;
  A = (b9_t*) malloc((size_t)nprimes * sizeof(b9_t));
  if (A == 0)
    goto return_null;

  for (pidx = 0; pidx < nprimes; pidx++) {
    uint32_t e, p = primes ? primes[pidx] : primes_small[pidx+1];
    if (p > n) break;
    e = binomial_prime_exp_32(n, k, p);
    if (e > 0) {
      if (e == 1) {
        b9_init_set_uv(&A[nprod], p);
      } else if (e == 2 && p <= 65535) {
        b9_init_set_uv(&A[nprod], p*p);
      } else if (e == 3 && p <= 1625) {
        b9_init_set_uv(&A[nprod], p*p*p);
      } else {
        b9_t base;
        b9_init_set_uv(&base, p);
        b9_init(&A[nprod]);
        b9_pow(&A[nprod], &base, e);
        b9_free(&base);
      }
      nprod++;
    }
  }

  if (nprod != 0) {
    b9_product(A, 0, (size_t)nprod-1);
    b9_move(&result, &A[0]);
  }

return_result:
  out = (char*) malloc((size_t)b9_length(&result) + 1);
  outlen = 0;
  if (out) {
    outlen = b9_get_str(out, &result);
    out[outlen] = '\0';
  }
  goto cleanup;

return_null:
  out = 0;
  outlen = 0;

cleanup:
  if (rlen) *rlen = outlen;
  if (A != 0) {
    for (pidx = 0; pidx < nprod; pidx++)
      b9_free(&A[pidx]);
    free(A);
  }
  if (primes != 0)
    free_prime_sieve_32(primes);
  b9_free(&result);
  return out;
}


static bool b9_binomial_product_32(uint32_t k, const b9_t *N, b9_t *A)
{
  uint32_t pidx, nprimes, *primes = 0;

  if (k <= primes_small[NPRIMES_SMALL-1]) {
    nprimes = NPRIMES_SMALL-1;
  } else {
    nprimes = range_prime_sieve_32(&primes, k, 0);
    if (nprimes == 0 || primes == 0)
      return 0;
  }

  for (pidx = 0; pidx < nprimes; pidx++) {
    uint32_t ppow, p = primes ? primes[pidx] : primes_small[pidx+1];
    if (p > k) break;
    for (ppow = p; ppow <= k; ) {
      uint32_t c, count = k / ppow;
      uint32_t rem = b9_mod_u32(N, ppow);
      uint32_t off = (rem == 0) ? 0 : ppow - rem;
      for (c = 0; c < count; c++, off += ppow)
        b9_divexact_u32(&A[off], p);
      if (ppow > k / p) break;
      ppow *= p;
    }
  }
  if (primes != 0)
    free_prime_sieve_32(primes);
  b9_product(A, 0, (size_t)k-1);
  return 1;
}

char* strint_binomial_u32(const char* nstr, STRLEN nlen, uint32_t k, STRLEN* rlen)
{
  b9_t start, term, result;
  b9_t *A = 0;
  char *out;
  STRLEN outlen;
  UV nsmall;
  int negn, negout;
  uint32_t i;

  if (rlen) *rlen = 0;

  negn = (nlen > 0 && *nstr == '-');
  if (nlen > 0 && (*nstr == '+' || *nstr == '-')) {
    nstr++;
    nlen--;
  }
  while (nlen > 1 && *nstr == '0') {
    nstr++;
    nlen--;
  }

  b9_init(&start);
  b9_init_set_uv(&result, 1);

  if (k == 0)
    goto return_result;

  if (nlen == 0 || (nlen == 1 && nstr[0] == '0')) {
    b9_set_uv(&result, 0);
    goto return_result;
  }

  b9_set_str(&start, nstr, nlen);  /* absolute value, or positive n */
  negout = negn && (k & 1);

  if (!negn) {
    if (str_to_uv_s(nstr, nlen, &nsmall)) {
      if (k > nsmall) {
        b9_set_uv(&result, 0);
        goto return_result;
      }
      if (k > nsmall-k)
        k = (uint32_t)(nsmall-k);
      if (k == 0)
        goto return_result;
    }
    if (k > 1) {
      b9_t sub;
      b9_init_set_uv(&sub, (UV)k - 1);
      b9_neg(&sub);
      b9_add(&start, &start, &sub);  /* start = n-k+1 */
      b9_free(&sub);
    }
  }

  if ((size_t)k > (size_t)MAX_SIZET / sizeof(b9_t))
    goto return_null;
  A = (b9_t*) malloc((size_t)k * sizeof(b9_t));
  if (A == 0)
    goto return_null;

  {
    b9_init_set(&term, &start);
    for (i = 0; i < k; i++) {
      b9_init_set(&A[i], &term);
      b9_add_u32(&term, &term, 1);
    }
    b9_free(&term);
  }

  if (!b9_binomial_product_32(k, &start, A))
    goto return_null;

  /* Result is in the first element of A. */
  b9_move(&result, &A[0]);

  if (negout)
    b9_neg(&result);

return_result:
  out = (char*) malloc((size_t)b9_length(&result) + 1);
  outlen = 0;
  if (out) {
    outlen = b9_get_str(out, &result);
    out[outlen] = '\0';
  }
  goto cleanup;

return_null:
  out = 0;
  outlen = 0;

cleanup:
  if (rlen) *rlen = outlen;

  if (A != 0) {
    for (i = 0; i < k; i++)
      b9_free(&A[i]);
    free(A);
  }
  b9_free(&start);
  b9_free(&result);
  return out;
}


/******************************************************************************/


static int strint_digit_value(unsigned char c)
{
  if (c >= '0' && c <= '9') return c - '0';
  if (c >= 'A' && c <= 'Z') return c - 'A' + 10;
  if (c >= 'a' && c <= 'z') return c - 'a' + 10;
  return -1;
}

static char* b9_from_little_digits(b9_t *A, size_t len, UV base, STRLEN* rlen)
{
  b9_t b;
  char *out;
  size_t i, k = len;
  STRLEN outlen;

  b9_init_set_uv(&b, base);

  while (k > 1) {
    for (i = 0; i < k-1; i += 2) {
      b9_mul(&A[i+1], &A[i+1], &b);
      b9_add(&A[i>>1], &A[i], &A[i+1]);
    }
    if (k & 1)
      b9_move(&A[k>>1], &A[k-1]);
    k = (k+1) >> 1;
    if (k > 1)
      b9_mul(&b, &b, &b);
  }

  out = (char*) malloc((size_t)b9_length(&A[0]) + 1);
  outlen = 0;
  if (out) {
    outlen = b9_get_str(out, &A[0]);
    out[outlen] = '\0';
  }
  if (rlen) *rlen = outlen;

  for (i = 0; i < len; i++)
    b9_free(&A[i]);
  free(A);
  b9_free(&b);
  return out;
}

char* strint_fromdigits(const UV* d, size_t len, UV base, STRLEN* rlen)
{
  b9_t *A;
  char *out;
  size_t i;

  if (len == 0) {
    out = (char*) malloc(2);
    if (!out) return 0;
    out[0] = '0';
    out[1] = '\0';
    if (rlen) *rlen = 1;
    return out;
  }
  if (len == 1) {
    char buf[24];
    STRLEN len1 = uv_to_str(buf, d[0]);
    out = (char*) malloc((size_t)len1 + 1);
    if (!out) return 0;
    memcpy(out, buf, (size_t)len1);
    out[len1] = '\0';
    if (rlen) *rlen = len1;
    return out;
  }
  if (len > (size_t)MAX_SIZET / sizeof(b9_t))
    return 0;

  A = (b9_t*) malloc(len * sizeof(b9_t));
  if (!A) return 0;
  for (i = 0; i < len; i++)
    b9_init_set_uv(&A[i], d[len-1-i]);

  return b9_from_little_digits(A, len, base, rlen);
}

int strint_fromdigitstring(UV* rn, char** rstr, STRLEN* rlen, const char* s, STRLEN len, UV base)
{
  b9_t *A;
  UV n = 0;
  STRLEN i;
  int overflow = 0;

  while (len > 0 && *s == '0') { s++;  len--; }
  if (len == 0) {
    *rn = 0;
    return 1;
  }

  for (i = 0; i < len; i++) {
    int d = strint_digit_value((unsigned char)s[i]);
    if (d < 0 || (UV)d >= base)
      croak("fromdigits: invalid digit for base %"UVuf, base);
    if (!overflow) {
      if (n > (UV_MAX - (UV)d) / base) {
        overflow = 1;
      } else {
        n = n * base + (UV)d;
      }
    }
  }

  if (!overflow) {
    *rn = n;
    return 1;
  }

  if (len > (STRLEN)((size_t)MAX_SIZET / sizeof(b9_t)))
    croak("fromdigits: input too large");
  A = (b9_t*) malloc((size_t)len * sizeof(b9_t));
  if (!A)
    croak("fromdigits: allocation failed");
  for (i = 0; i < len; i++) {
    int d = strint_digit_value((unsigned char)s[len-1-i]);
    b9_init_set_uv(&A[i], (UV)d);
  }

  *rstr = b9_from_little_digits(A, (size_t)len, base, rlen);
  if (!*rstr)
    croak("fromdigits: allocation failed");
  return 2;
}

static void b9_mul_2pow_chunk(b9_t *out, const b9_t *a, unsigned bits)
{
  if (bits == 0) {
    if (out != a) {
      b9_ensure(out, a->n);
      if (a->n) memcpy(out->d, a->d, a->n * sizeof(b9limb_t));
      out->n = a->n;  out->neg = a->neg;
    }
  } else if (bits == 32) {
    b9_mul_u32(out, a, UINT32_C(65536));
    b9_mul_u32(out, out, UINT32_C(65536));
  } else {
    b9_mul_u32(out, a, UINT32_C(1) << bits);
  }
}

int strint_radix_to_int(UV* rn, char** rstr, STRLEN* rlen, const char* s, STRLEN len, UV base)
{
  b9_t R;
  UV n = 0;
  STRLEN i, start;

  if (rlen) *rlen = 0;
  if (rstr) *rstr = 0;
  while (len > 0 && *s == '0') { s++;  len--; }
  if (base != 2 && base != 16) return 0;

  if (len <= BITS_PER_WORD / (base == 16 ? 4 : 1)) {
    if (base == 16) {
      for (i = 0; i < len; i++)
        n = n * base + (UV)strint_digit_value((unsigned char)s[i]);
    } else {
      for (i = 0; i < len; i++)
        n = n * base + (UV)(s[i] - '0');
    }
    *rn = n;
    return 1;
  }

  b9_init_set_uv(&R, 0);
  start = 0;
  while (start < len) {
    STRLEN chunk_len = len - start;
    uint32_t chunk = 0;
    unsigned bits;

    if (base == 16 && chunk_len > 8)
      chunk_len = 8;
    else if (chunk_len > 32)
      chunk_len = 32;

    bits = (unsigned)chunk_len * ((base == 16) ? 4 : 1);
    for (i = 0; i < chunk_len; i++) {
      if (base == 16)
        chunk = (chunk << 4) | (uint32_t)strint_digit_value((unsigned char)s[start+i]);
      else
        chunk = (chunk << 1) | (uint32_t)(s[start+i] - '0');
    }
    b9_mul_2pow_chunk(&R, &R, bits);
    b9_add_u32(&R, &R, chunk);
    start += chunk_len;
  }

  *rstr = (char*) malloc((size_t)b9_length(&R) + 1);
  if (*rstr == 0) {
    b9_free(&R);
    croak("toint: allocation failed");
  }
  *rlen = b9_get_str(*rstr, &R);
  (*rstr)[*rlen] = '\0';
  b9_free(&R);
  return 2;
}

/******************************************************************************/
/*                              MULTIPLICATION                                */
/******************************************************************************/

/* Multiply two signed decimal integer strings.
 * Write result to out (caller ensures at least alen+blen bytes).
 * Returns the result length (no NULL terminator written). */
STRLEN strint_mul(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen)
{
  b9_t A, B, R;
  STRLEN rlen;

  b9_init_set_str(&A, a, alen);
  b9_init_set_str(&B, b, blen);
  b9_init(&R);
  b9_mul(&R, &A, &B);

  rlen = b9_get_str(out, &R);
  b9_free(&A);  b9_free(&B);  b9_free(&R);

  return rlen;
}

STRLEN strint_muladd_s(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen, const char* c, STRLEN clen, bool negate_c)
{
  b9_t A, B, C;
  STRLEN rlen;

  b9_init_set_str(&A, a, alen);
  b9_init_set_str(&B, b, blen);
  b9_init_set_str(&C, c, clen);
  if (negate_c) b9_neg(&C);

  b9_mul(&A, &A, &B);
  b9_add(&A, &A, &C);

  rlen = b9_get_str(out, &A);
  b9_free(&A);  b9_free(&B);  b9_free(&C);

  return rlen;
}

STRLEN strint_muladdmod_s(char* out,
                          const char* a, STRLEN alen,
                          const char* b, STRLEN blen,
                          const char* c, STRLEN clen, bool negate_c,
                          const char* m, STRLEN mlen)
{
  b9_t A, B, C, M, R;
  STRLEN rlen;

  b9_init_set_str(&A, a, alen);
  b9_init_set_str(&B, b, blen);
  b9_init_set_str(&C, c, clen);
  b9_init_set_str(&M, m, mlen);
  if (M.n == 0) {
    b9_free(&A);  b9_free(&B);  b9_free(&C);  b9_free(&M);
    return 0;
  }
  if (M.neg) M.neg = 0;         /* work with |M| */
  b9_init(&R);

  b9_mulmod(&A, &R, &A, &B, &M);
  if (C.n != 0) {
    if (negate_c)  b9_submod(&A, &R, &A, &C, &M);
    else           b9_addmod(&A, &R, &A, &C, &M);
  }

  rlen = b9_get_str(out, &A);
  b9_free(&A);  b9_free(&B);  b9_free(&C);  b9_free(&M);  b9_free(&R);

  return rlen;
}

STRLEN strint_powmod(char* out,
                     const char* a, STRLEN alen,
                     const char* e, STRLEN elen,
                     const char* m, STRLEN mlen)
{
  b9_t A, E, M, R, T;
  STRLEN rlen;

  b9_init_set_str(&A, a, alen);
  b9_init_set_str(&E, e, elen);
  b9_init_set_str(&M, m, mlen);
  if (M.n == 0 || E.neg) {
    b9_free(&A);  b9_free(&E);  b9_free(&M);
    return 0;
  }
  if (M.neg) M.neg = 0;         /* work with |M| */
  if (M.n == 1 && M.d[0] == 1) {
    out[0] = '0';
    b9_free(&A);  b9_free(&E);  b9_free(&M);
    return 1;
  }

  b9_init(&R);  b9_init(&T);
  b9_powmod(&R, &T, &A, &E, &M);
  rlen = b9_get_str(out, &R);
  b9_free(&A);  b9_free(&E);  b9_free(&M);  b9_free(&R);  b9_free(&T);

  return rlen;
}

/******************************************************************************/
/*                             EXPONENTIATION                                 */
/******************************************************************************/

/* Raise a to non-negative integer power exp.
 * out must have at least limit bytes.
 * Returns 0 if result would not fit in limit characters. */
STRLEN strint_pow(char* out, const char* a, STRLEN alen, UV exp, STRLEN limit)
{
  b9_t base, result;
  STRLEN rlen;

  /* exp = 0: a^0 = 1 for all a (including 0) */
  if (exp == 0) {
    out[0] = '1';
    return 1;
  }
  if (exp == 1) {
    if (alen > limit) return 0;
    memcpy(out, a, alen);
    return alen;
  }

  b9_init_set_str(&base, a, alen);
  b9_init(&result);
  b9_pow(&result, &base, exp);
  rlen = b9_length(&result) > limit  ?  0  :  b9_get_str(out, &result);
  b9_free(&base);  b9_free(&result);
  return rlen;
}


/******************************************************************************/
/*                                 DIVISION                                   */
/******************************************************************************/


/* Signed floor division and remainder of decimal integer strings.
 * Either qout or rout (and its length pointer) may be NULL.
 * qout needs alen+1 bytes; rout needs blen bytes.
 * Remainder has the sign of b (floor convention).
 * Returns false if b = 0, true otherwise. */
bool strint_fdivrem(char* qout, STRLEN* qlen, char* rout, STRLEN* rlen,
                   const char* a, STRLEN alen,
                   const char* b, STRLEN blen)
{
  b9_t ba, bb, bq, br;
  b9_init_set_str(&ba, a, alen);
  b9_init_set_str(&bb, b, blen);
  if (bb.n == 0) { b9_free(&ba);  b9_free(&bb);  return 0; }
  b9_init(&bq);  b9_init(&br);
  b9_fdivrem(&bq, &br, &ba, &bb);
  if (qout) { STRLEN l = b9_get_str(qout, &bq);  if (qlen) *qlen = l; }
  if (rout) { STRLEN l = b9_get_str(rout, &br);  if (rlen) *rlen = l; }
  b9_free(&ba);  b9_free(&bb);  b9_free(&bq);  b9_free(&br);
  return 1;
}

STRLEN strint_divint(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen)
{
  STRLEN len = 0;
  if (!strint_fdivrem(out, &len, 0, 0, a,alen, b,blen)) return 0;
  return len;
}

STRLEN strint_cdivint(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen)
{
  STRLEN len;
  b9_t ba, bb, bq, br;
  b9_init_set_str(&ba, a, alen);
  b9_init_set_str(&bb, b, blen);
  if (bb.n == 0) { b9_free(&ba);  b9_free(&bb);  return 0; }
  b9_init(&bq);  b9_init(&br);
  b9_fdivrem(&bq, &br, &ba, &bb);
  if (br.n > 0)
    b9_add_uv(&bq, &bq, 1);
  len = b9_get_str(out, &bq);
  b9_free(&ba);  b9_free(&bb);  b9_free(&bq);  b9_free(&br);
  return len;
}

STRLEN strint_modint(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen)
{
  STRLEN len = 0;
  if (!strint_fdivrem(0, 0, out, &len, a,alen, b,blen)) return 0;
  return len;
}

UV strint_moduv(const char* a, STRLEN alen, UV b)
{
  UV r = 0;
  STRLEN i;
  if (b <= 1) return 0;
  strint_strip(&a, &alen);
  for (i = 0; i < alen; i++)
    r = muladdmod(r, 10, (UV)(a[i] - '0'), b);
  return r;
}

/******************************************************************************/
/*                        SMALL FACTOR EXTRACTION                             */
/******************************************************************************/

STRLEN strint_remove_small_factors(char* str_out, UV* uv_out,
                                   UV* out_f, int* nf,
                                   const char* a, STRLEN alen)
{
#if BITS_PER_WORD == 64
  static const UV P = UVCONST(614889782588491410); /* primorial(47) */
  static const uint32_t sprimes[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47};
  static const int nsp = 15;
#else
  static const UV P = UVCONST(223092870);           /* primorial(23) */
  static const uint32_t sprimes[] = {2,3,5,7,11,13,17,19,23};
  static const int nsp = 9;
#endif
  UV g;
  int pi;

  g = gcduv(strint_moduv(a, alen, P), P);

  if (g > 1) {
    b9_t n;
    b9_init_set_str(&n, a, alen);

    for (pi = 0; pi < nsp && g > 1; pi++) {
      uint32_t p = sprimes[pi];
      if (g % p != 0) continue;

      while (b9_mod_u32(&n, p) == 0) {
        b9_divexact_u32(&n, p);
        out_f[(*nf)++] = p;

        /* Early UV detection: switch to UV arithmetic for the remaining tail */
        if ((STRLEN)n.n * B9_DIGS < sizeof(UV) * 3) {
          UV v = b9_get_uv(&n);
          b9_free(&n);
          while (v % p == 0) { v /= p; out_f[(*nf)++] = p; }
          while (g % p == 0) g /= p;
          for (pi++; pi < nsp && g > 1; pi++) {
            uint32_t q = sprimes[pi];
            if (g % q == 0) {
              while (v % q == 0) { v /= q; out_f[(*nf)++] = q; }
              while (g % q == 0) g /= q;
            }
          }
          *uv_out = v;
          return 0;
        }
      }
      while (g % p == 0) g /= p;
    }
    {
      STRLEN rlen = b9_get_str(str_out, &n);
      b9_free(&n);
      if (str_to_uv_s(str_out, rlen, uv_out)) return 0;
      return rlen;
    }
  }

  /* g == 1: no small factors */
  strint_strip(&a, &alen);
  if (str_out != a) memmove(str_out, a, alen);
  if (str_to_uv_s(str_out, alen, uv_out)) return 0;
  return alen;
}

STRLEN strint_trial_factor(char* str_out, UV* uv_out,
                           UV* out_f, int* nf,
                           const char* a, STRLEN alen,
                           const uint32_t* primes, uint32_t nprimes)
{
  b9_t n;
  int have_b9 = 0;
  uint32_t pi = 0;

  /* Process primes in batches whose product fits in a single integer.
   * Before b9 init: use UV-sized batches + strint_moduv (avoids b9 init
   * entirely for the common no-factor case).
   * After b9 init: cap batches at UINT32_MAX so b9_mod_u32 can be used. */
  while (pi < nprimes) {
    UV batch_cap  = have_b9 ? (UV)UINT32_MAX : UV_MAX;
    UV batch_prod = 1;
    uint32_t j, batch_start = pi;

    while (pi < nprimes && batch_prod <= batch_cap / primes[pi])
      batch_prod *= (UV)primes[pi++];

    /* Quick filter: skip batch if gcd(n mod batch_prod, batch_prod) == 1 */
    UV rem  = have_b9 ? (UV)b9_mod_u32(&n, (uint32_t)batch_prod)
                      : strint_moduv(a, alen, batch_prod);
    UV g    = gcduv(rem, batch_prod);
    if (g == 1) continue;

    /* At least one prime in this batch divides n */
    if (!have_b9) {
      b9_init_set_str(&n, a, alen);
      n.neg = 0;
      have_b9 = 1;
    }

    for (j = batch_start; j < pi; j++) {
      uint32_t p = primes[j];
      if (g % p != 0) continue;  /* not in gcd, skip */

      /* p divides n (guaranteed by gcd); divide it out completely */
      do {
        b9_divexact_u32(&n, p);
        out_f[(*nf)++] = p;
        if ((STRLEN)n.n * B9_DIGS < sizeof(UV) * 3) {
          UV v = b9_get_uv(&n);
          b9_free(&n);
          /* p may still divide v (e.g. p^2 | original n) */
          while (v % p == 0) { v /= p;  out_f[(*nf)++] = p; }
          /* Finish current batch in UV */
          for (j++; j < pi; j++) {
            UV q = (UV)primes[j];
            if (g % q != 0) continue;
            while (v % q == 0) { v /= q;  out_f[(*nf)++] = q; }
          }
          /* Remaining batches: plain UV loop with q*q > v termination */
          for (; pi < nprimes; pi++) {
            UV q = (UV)primes[pi];
            if (q > v / q) break;
            while (v % q == 0) { v /= q;  out_f[(*nf)++] = q; }
          }
          *uv_out = v;
          return 0;
        }
      } while (b9_mod_u32(&n, p) == 0);

      while (g % p == 0) g /= p;  /* consumed this prime from g */
    }
  }

  if (!have_b9) {
    /* No factors found; copy the stripped input as the cofactor */
    strint_strip(&a, &alen);
    if (str_out != a) memmove(str_out, a, alen);
    if (str_to_uv_s(str_out, alen, uv_out)) return 0;
    return alen;
  }

  /* All primes tried; write out the remaining cofactor */
  {
    STRLEN rlen = b9_get_str(str_out, &n);
    b9_free(&n);
    if (str_to_uv_s(str_out, rlen, uv_out)) return 0;
    return rlen;
  }
}


/******************************************************************************/
/*                           INTEGER LOG AND ROOT                             */
/******************************************************************************/

/* Returns floor(log_base(a)), or UV_MAX on error (a < 1 or base < 2).
 * Uses a floating-point estimate to land within 1-2 of the answer,
 * then verifies/adjusts with b9 arithmetic. */
UV strint_logint(const char* a, STRLEN alen, UV base)
{
  STRLEN nd, i;
  UV k;
  double approx, log10_a, k_est_f;
  b9_t bn, bbase, pow, tmp;

  if (base < 2) return UV_MAX;

  if (strint_strip(&a, &alen) || (alen == 1 && a[0] == '0'))
    return UV_MAX;

  /* Floating-point estimate: log10(a) via leading digits + digit count */
  nd = (alen < 15) ? alen : 15;
  approx = 0.0;
  for (i = 0; i < nd; i++)
    approx = approx * 10.0 + (double)(a[i] - '0');
  log10_a = (double)(alen - nd) + log10(approx);
  k_est_f = log10_a / log10((double)base);
  k = (k_est_f <= 1.0) ? 0 : (UV)k_est_f - 1;  /* start one below estimate */

  b9_init_set_str(&bn, a, alen);
  b9_init_set_uv(&bbase, base);
  b9_init(&pow);
  b9_pow(&pow, &bbase, k);
  b9_init(&tmp);

  /* Float overshot (very rare): step down by dividing until base^k ≤ a */
  while (k > 0 && b9_cmp(&pow, &bn) > 0) {
    b9_fdivrem(&tmp, NULL, &pow, &bbase);
    b9_swap(&pow, &tmp);
    k--;
  }

  /* Step up while base^(k+1) ≤ a */
  for (;;) {
    b9_mul(&tmp, &pow, &bbase);
    if (b9_cmp(&tmp, &bn) > 0) break;
    b9_swap(&pow, &tmp);
    k++;
  }

  b9_free(&bn);  b9_free(&bbase);  b9_free(&pow);  b9_free(&tmp);
  return k;
}

/* floor(n^(1/k)) via Newton's method, starting from above.
 * r must be initialized.  n and its canonical decimal string must represent
 * the same value. */
static void b9_root_ui(b9_t *r, const b9_t *n,
                       const char* a, STRLEN alen, UV k)
{
  STRLEN nd, i;
  double approx, log10_a, r_est_f;
  b9_t bkm1, bk, pk1, q, tmp;

  /* Float estimate of floor(log10(n^(1/k))) */
  nd = (alen < 15) ? alen : 15;
  approx = 0.0;
  for (i = 0; i < nd; i++)
      approx = approx * 10.0 + (double)(a[i] - '0');
  log10_a = (double)(alen - nd) + log10(approx);
  r_est_f = log10_a / (double)k;

  b9_init_set_uv(&bkm1, k - 1);
  b9_init_set_uv(&bk,   k);
  b9_init(&pk1);  b9_init(&q);  b9_init(&tmp);

  /* Initial r: float within ~0.02% of true root, biased above.
   * r_est_f = q_exp + frac.  r_mantissa ≈ 10^frac * 1.0002e8 (9 sig figs).
   * For q_exp >= mlen-1: r = r_mantissa * 10^(q_exp+1-mlen) (char temp, one malloc).
   * For q_exp < mlen-1:  r_est_f < 9, r fits in a UV. */
  {
    UV q_exp = (r_est_f >= 1.0) ? (UV)r_est_f : 0;
    double frac = r_est_f - (double)q_exp;
    UV r_mantissa = (UV)(pow(10.0, frac) * 1.0002e8) + 1;
    char mant_buf[12];
    STRLEN mlen = uv_to_str(mant_buf, r_mantissa);
    if (q_exp + 1 >= mlen) {
      char *r_char = (char*) malloc(q_exp + 4);
      memcpy(r_char, mant_buf, mlen);
      memset(r_char + mlen, '0', q_exp + 1 - mlen);
      b9_set_str(r, r_char, q_exp + 1);
      free(r_char);
    } else {
      b9_set_uv(r, (UV)(pow(10.0, r_est_f) * 1.0002) + 1);
    }
  }

  /* Newton: r_new = floor(((k-1)*r + floor(n/r^(k-1))) / k).
   * Converges from above; terminate when r_new >= r. */
  for (;;) {
    /* pk1 = r^(k-1) */
    b9_pow(&pk1, r, k - 1);

    /* q = floor(n / r^(k-1)) */
    b9_fdivrem(&q, NULL, n, &pk1);

    /* tmp = (k-1)*r + q */
    b9_mul(&tmp, &bkm1, r);
    b9_add(&tmp, &tmp, &q);

    /* pk1 = floor(tmp / k)  (reuse as r_new) */
    b9_fdivrem(&pk1, NULL, &tmp, &bk);

    if (b9_cmp(&pk1, r) >= 0) break;

    /* r = r_new */
    b9_swap(r, &pk1);
  }

  b9_free(&bkm1);  b9_free(&bk);  b9_free(&pk1);
  b9_free(&q);     b9_free(&tmp);
}

/* floor(n^(1/k)) via Newton's method, starting from above.
 * out must have at least alen+2 bytes.
 * Returns the length of the result, or 0 on error (k==0, n<0). */
STRLEN strint_rootint(char* out, const char* a, STRLEN alen, UV k)
{
  STRLEN rlen;
  b9_t n, r;

  if (k == 0) return 0;
  if (strint_strip(&a, &alen)) return 0;
  if (k == 1) { memcpy(out, a, alen); return alen; }
  if (alen == 1 && (a[0] == '0' || a[0] == '1')) {
    out[0] = a[0];
    return 1;
  }

  b9_init_set_str(&n, a, alen);
  b9_init(&r);
  b9_root_ui(&r, &n, a, alen, k);
  rlen = b9_get_str(out, &r);
  b9_free(&n);  b9_free(&r);
  return rlen;
}

bool strint_is_perfect_square(const char* a, STRLEN alen)
{
  uint32_t m, r99, r105;
  b9_t n, r, square;
  bool ret = 0;

  if (strint_strip(&a, &alen)) return 0;
  if (alen == 1 && (a[0] == '0' || a[0] == '1')) return 1;

  b9_init_set_str(&n, a, alen);

  /* One b9 remainder supplies filters modulo multiple small integers.
   * The three filters remove over 98% of candidates. */
  m = b9_mod_u32(&n, 64 * 105 * 33);

#if HAVE_UINT64
  if ((UINT64_C(1) << (m & 63)) & UINT64_C(0xfdfdfdedfdfcfdec)) goto done;
#else
  if ((UINT32_C(1) << (m & 31)) & UINT32_C(0xfdfcfdec)) goto done;
#endif
  r105 = m % 105;
  if ((r105*0xd24554cd) & (r105*0x0929579a) & 0x38020141) goto done;
  r99  = m % 99;
  if ((r99*0x5411171d) & (r99*0xe41dd1c7) & 0x80028a80) goto done;

  b9_init(&r);  b9_init(&square);
  b9_root_ui(&r, &n, a, alen, 2);
  b9_mul(&square, &r, &r);
  ret = b9_cmp(&square, &n) == 0;
  b9_free(&r);  b9_free(&square);
done:
  b9_free(&n);
  return ret;
}

STRLEN strint_lshiftint(char* out, const char* a, STRLEN alen, UV k)
{
  b9_t n, pow2;
  STRLEN rlen;
  if (k == 0) { memcpy(out, a, alen); return alen; }
  b9_init_set_str(&n, a, alen);
  b9_init_set_pow2(&pow2, k);
  b9_mul(&n, &n, &pow2);
  rlen = b9_get_str(out, &n);
  b9_free(&n);  b9_free(&pow2);
  return rlen;
}

STRLEN strint_rshiftint(char* out, const char* a, STRLEN alen, UV k)
{
  b9_t n, pow2, q;
  STRLEN rlen;
  int neg;
  if (k == 0) { memcpy(out, a, alen); return alen; }
  b9_init_set_str(&n, a, alen);
  b9_init_set_pow2(&pow2, k);
  neg = n.neg;  n.neg = 0;
  b9_init(&q);
  b9_fdivrem(&q, NULL, &n, &pow2);
  q.neg = (q.n != 0) ? neg : 0;
  rlen = b9_get_str(out, &q);
  b9_free(&n);  b9_free(&pow2);  b9_free(&q);
  return rlen;
}

STRLEN strint_rashiftint(char* out, const char* a, STRLEN alen, UV k)
{
  b9_t n, pow2, q;
  STRLEN rlen;
  if (k == 0) { memcpy(out, a, alen); return alen; }
  b9_init_set_str(&n, a, alen);
  b9_init_set_pow2(&pow2, k);
  b9_init(&q);
  b9_fdivrem(&q, NULL, &n, &pow2);
  rlen = b9_get_str(out, &q);
  b9_free(&n);  b9_free(&pow2);  b9_free(&q);
  return rlen;
}
