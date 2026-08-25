#include <stdlib.h>
#include <string.h>

#include "b9.h"

/******************************************************************************/
/* Internal allocation helpers. */
/******************************************************************************/

/* Internal allocations cannot propagate failure through the b9 interface. */
static void* b9_xmalloc(size_t count, size_t size)
{
  void *p;
  if (size != 0 && count > (size_t)MAX_SIZET / size)
    croak("internal: b9 allocation too large");
  p = malloc(count * size);
  if (p == 0 && count != 0 && size != 0)
    croak("internal: b9 allocation failed");
  return p;
}

static void* b9_xcalloc(size_t count, size_t size)
{
  void *p;
  if (size != 0 && count > (size_t)MAX_SIZET / size)
    croak("internal: b9 allocation too large");
  p = calloc(count, size);
  if (p == 0 && count != 0 && size != 0)
    croak("internal: b9 allocation failed");
  return p;
}

static void* b9_xrealloc(void *ptr, size_t count, size_t size)
{
  void *p;
  if (size != 0 && count > (size_t)MAX_SIZET / size)
    croak("internal: b9 allocation too large");
  p = realloc(ptr, count * size);
  if (p == 0 && count != 0 && size != 0)
    croak("internal: b9 allocation failed");
  return p;
}

void b9_free(b9_t *x)
{
  if (x->d != (b9limb_t*)x->d_small) free(x->d);
  b9_init(x);
}

void b9_neg(b9_t *x)
{
  x->neg = x->n == 0 ? 0 : !x->neg;
}


/******************************************************************************/
/* Base-10^N arithmetic. */
/******************************************************************************/

void b9_ensure(b9_t *x, uint32_t need)
{
  b9limb_t *newd;
  if (x->alloc >= need) return;
  if (x->d == (b9limb_t*)x->d_small) {
    newd = (b9limb_t*)b9_xmalloc((size_t)need, sizeof(b9limb_t));
    if (x->n) memcpy(newd, x->d_small, x->n * sizeof(b9limb_t));
  } else {
    newd = (b9limb_t*)b9_xrealloc(x->d, (size_t)need,
                                  sizeof(b9limb_t));
  }
  x->d = newd;
  x->alloc = need;
}

/* Exchange two b9_t values, correctly handling inline storage. */
void b9_swap(b9_t *a, b9_t *b)
{
  int a_inline = (a->d == (b9limb_t*)a->d_small);
  int b_inline = (b->d == (b9limb_t*)b->d_small);
  b9limb_t small[B9_INLINE_LIMBS];
  b9limb_t *d;
  b9_t *t;
  uint32_t an, bn, alloc, n;
  int aneg, bneg, neg;

  if (a_inline && b_inline) {
    an = a->n;  aneg = a->neg;
    bn = b->n;  bneg = b->neg;
    if (an) memcpy(small, a->d, (size_t)an * sizeof(b9limb_t));
    if (bn) memcpy(a->d, b->d, (size_t)bn * sizeof(b9limb_t));
    if (an) memcpy(b->d, small, (size_t)an * sizeof(b9limb_t));
    a->n = bn;  a->neg = bneg;
    b->n = an;  b->neg = aneg;
    return;
  }

  if (!a_inline && !b_inline) {
    d = a->d;          a->d = b->d;          b->d = d;
    alloc = a->alloc;  a->alloc = b->alloc;  b->alloc = alloc;
    n = a->n;          a->n = b->n;          b->n = n;
    neg = a->neg;      a->neg = b->neg;      b->neg = neg;
    return;
  }

  /* Put the inline value in a and the heap value in b. */
  if (!a_inline) { t = a;  a = b;  b = t; }
  an = a->n;  aneg = a->neg;
  if (an) memcpy(small, a->d, (size_t)an * sizeof(b9limb_t));
  a->d = b->d;  a->alloc = b->alloc;  a->n = b->n;  a->neg = b->neg;
  b->d = (b9limb_t*)b->d_small;
  b->alloc = B9_INLINE_LIMBS;  b->n = an;  b->neg = aneg;
  if (an) memcpy(b->d, small, (size_t)an * sizeof(b9limb_t));
}

/* Free dst, move src into dst.  src is blanked. */
void b9_move(b9_t *dst, b9_t *src)
{
  b9_free(dst);
  if (src->d == (b9limb_t*)src->d_small) {
    if (src->n)
      memcpy(dst->d, src->d, (size_t)src->n * sizeof(b9limb_t));
  } else {
    dst->d = src->d;
    dst->alloc = src->alloc;
  }
  dst->n = src->n;
  dst->neg = src->neg;
  b9_init(src);
}

/* Copy src into an initialized dst.  src is unchanged. */
void b9_set(b9_t *dst, const b9_t *src)
{
  if (dst != src) {
    b9_ensure(dst, src->n);
    if (src->n) memcpy(dst->d, src->d, src->n * sizeof(b9limb_t));
    dst->n = src->n;
    dst->neg = src->neg;
  }
}

/* Exact length in decimal digits including possible sign character */
STRLEN b9_length(const b9_t *x)
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
void b9_set_str(b9_t *x, const char *s, STRLEN len)
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

void b9_init_set_str(b9_t *x, const char *s, STRLEN len)
{
  b9_init(x);
  b9_set_str(x, s, len);
}

/* Write x to buf as a signed decimal string (no NUL).
 * buf must have at least x->n * B9_DIGS + 2 bytes.
 * Returns string length. */
STRLEN b9_get_str(char *buf, const b9_t *x)
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
void b9_set_uv(b9_t *x, UV v)
{
  uint32_t n = 0;
  b9_ensure(x, (uint32_t)B9_NLIMBS(B9_UV_DEC_DIGS));
  x->neg = 0;
  if (v == 0) { x->n = 0;  return; }
  while (v > 0) { x->d[n++] = (b9limb_t)((UV)(v % B9_BASE)); v /= (UV)B9_BASE; }
  x->n = n;
}

void b9_init_set_uv(b9_t *x, UV v)
{
  b9_init(x);
  b9_set_uv(x, v);
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
UV b9_get_uv(const b9_t* x) {
  UV v = 0;
  uint32_t i;
  for (i = x->n; i-- > 0; )
    v = v * (UV)B9_BASE + (UV)x->d[i];
  return v;
}

/* Simpler helpers */
static void b9_init_set_acc(b9_t *x, b9acc_t v)
  { b9_init(x);  b9_set_acc(x,v); }
void b9_init_set(b9_t *x, const b9_t *y)
{
  b9_init(x);
  b9_set(x, y);
}


/* Signed comparison.  Returns -1, 0, or 1. */
int b9_cmp(const b9_t *a, const b9_t *b)
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
void b9_add(b9_t *out, const b9_t *a, const b9_t *b)
{
  int cmp, use_tmp;
  b9limb_t *dst;
  uint32_t need, rn;

  if (a->n == 0) {
    b9_set(out, b);
    return;
  }
  if (b->n == 0) {
    b9_set(out, a);
    return;
  }

  use_tmp = (out->d == a->d || out->d == b->d);
  if (a->neg == b->neg) {
    need = (a->n > b->n ? a->n : b->n) + 1;
    if (use_tmp) {
      dst = (b9limb_t*)b9_xmalloc(need, sizeof(b9limb_t));
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
      dst = (b9limb_t*)b9_xmalloc(need, sizeof(b9limb_t));
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
  b9_t nb;
  nb.d = b->d;  nb.alloc = b->alloc;  nb.n = b->n;  nb.neg = b->neg;
  b9_neg(&nb);
  b9_add(out, a, &nb);
}
#endif

void b9_add_u32(b9_t *out, const b9_t *a, uint32_t v)
{
  b9acc_t carry;
  uint32_t i, need;

  if (v == 0) {
    b9_set(out, a);
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

void b9_add_uv(b9_t *out, const b9_t *a, UV v)
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

void b9_mul(b9_t *out, const b9_t *a, const b9_t *b)
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
    acc = (b9acc_t*)b9_xcalloc((size_t)rn + 1, sizeof(b9acc_t));
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

void b9_mul_u32(b9_t *out, const b9_t *a, uint32_t v)
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
void b9_pow(b9_t *out, const b9_t *a, UV exp)
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
void b9_init_set_pow2(b9_t *x, UV k)
{
  if (k < BITS_PER_WORD) {
    b9_init_set_uv(x, UVCONST(1) << k);
  } else {
    b9_init_set_uv(x, 2);
    b9_pow(x, x, k);
  }
}

/* out = a * 2^bits.  out may alias a. */
void b9_mul_2exp(b9_t *out, const b9_t *a, UV bits)
{
  if (bits == 0) {
    b9_set(out, a);
  } else if (b9_is_zero(a)) {
    out->n = 0;
    out->neg = 0;
  } else if (bits < 32) {
    b9_mul_u32(out, a, UINT32_C(1) << bits);
  } else if (bits == 32) {
    b9_mul_u32(out, a, UINT32_C(65536));
    b9_mul_u32(out, out, UINT32_C(65536));
  } else {
    b9_t pow2;
    b9_init_set_pow2(&pow2, bits);
    b9_mul(out, a, &pow2);
    b9_free(&pow2);
  }
}

/* Signed floor division and remainder (floor convention: rem has sign of b).
 * q or r may be NULL (that output is skipped).
 * b must be non-zero.  q and r must not alias a or b. */
void b9_fdivrem(b9_t *q, b9_t *r, const b9_t *a, const b9_t *b)
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
    b9_set(rp, a);
    b9_abs(rp);
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
    size_t ulen;
    uint32_t d;

    ulen = (size_t)na + 1;
    u = (b9limb_t*)b9_xmalloc(ulen + (size_t)nb, sizeof(b9limb_t));
    v = u + ulen;

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

    free(u);
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
    { b9limb_t *tmp = (b9limb_t*)b9_xmalloc((size_t)nb,
                                            sizeof(b9limb_t));
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
uint32_t b9_mod_u32(const b9_t* a, uint32_t p) {
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
void b9_divexact_u32(b9_t* a, uint32_t p) {
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

uint32_t b9_divrem_u32_inplace(b9_t* a, uint32_t p)
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


int b9_cmp_abs(const b9_t *a, const b9_t *b)
{
  return b9mag_cmp(a->d, a->n, b->d, b->n);
}

static void b9_reduce_mod(b9_t *out, b9_t *tmp, const b9_t *a, const b9_t *m)
{
  b9_fdivrem(NULL, tmp, a, m);
  b9_move(out, tmp);
}

void b9_addmod(b9_t *out, b9_t *tmp, const b9_t *a, const b9_t *b, const b9_t *m)
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

void b9_submod(b9_t *out, b9_t *tmp, const b9_t *a, const b9_t *b, const b9_t *m)
{
  b9_t nb;
  nb.d = b->d;  nb.alloc = b->alloc;  nb.n = b->n;  nb.neg = b->neg;
  b9_neg(&nb);
  b9_addmod(out, tmp, a, &nb, m);
}

static void b9_mulmod_reduced(b9_t *out, b9_t *tmp, const b9_t *a, const b9_t *b, const b9_t *m)
{
  b9_mul(out, a, b);
  b9_reduce_mod(out, tmp, out, m);
}

void b9_mulmod(b9_t *out, b9_t *tmp, const b9_t *a, const b9_t *b, const b9_t *m)
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

void b9_powmod(b9_t *out, b9_t *tmp, const b9_t *a, const b9_t *e, const b9_t *m)
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

void b9_product(b9_t A[], size_t a, size_t b) {
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

void b9_product_u32(b9_t *out, const uint32_t A[], size_t len)
{
  b9_t *B = 0;
  b9acc_t prod = 1;
  size_t i, max_chunks, nprod = 0;
  const size_t n_u32_per_acc = sizeof(b9acc_t) / sizeof(uint32_t);

  if (len == 0) {
    b9_set_uv(out, 1);
    return;
  }

  max_chunks = len / n_u32_per_acc + (len % n_u32_per_acc != 0);
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
