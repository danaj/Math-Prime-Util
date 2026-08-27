#include <stdlib.h>
#include <string.h>

#include "b9.h"

/* Use schoolbook multiplication below this many limbs (must be >= 3). */
#define B9_KARATSUBA_THRESHOLD 400

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

#if HAVE_UINT64
static bool b9_get_u32_checked(uint32_t *result, const b9_t *x)
{
  uint64_t v = 0;
  uint32_t i;

  if (x->n > (uint32_t)B9_NLIMBS(10)) return 0;
  for (i = x->n; i-- > 0; )
    v = v * B9_BASE + x->d[i];
  if (v > UINT32_MAX) return 0;
  *result = (uint32_t)v;
  return 1;
}
#endif

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

/* out = a * b (signed).  Every multiplication path permits out to alias
 * either input. */

/* Keep common schoolbook products and Karatsuba leaves off the heap without
 * making a 4 KiB stack array when b9acc_t is uint64_t. */
#define B9_MUL_STACK_LIMBS 448

/* Multiply a magnitude by one limb.  out may alias a. */
static void b9_mul_limb_abs(b9_t *out, const b9_t *a, b9limb_t v)
{
  b9acc_t carry = 0;
  uint32_t i;

  if (v == 0 || a->n == 0) {
    out->n = 0;
    out->neg = 0;
    return;
  }
  if (v == 1) {
    b9_set(out, a);
    out->neg = 0;
    return;
  }
  if (a->n == UINT32_MAX)
    croak("internal: b9 multiplication too large");
  b9_ensure(out, a->n + 1);
  for (i = 0; i < a->n; i++) {
    b9acc_t product = (b9acc_t)a->d[i] * v + carry;
    out->d[i] = (b9limb_t)(product % B9_BASE);
    carry = product / B9_BASE;
  }
  out->n = a->n;
  if (carry != 0) out->d[out->n++] = (b9limb_t)carry;
  out->neg = 0;
}

/* All input limbs are consumed before out is resized or written. */
static void b9_mul_schoolbook(b9_t *out, const b9_t *a, const b9_t *b)
{
  uint32_t i, j, rn;
  b9acc_t  stack_acc[B9_MUL_STACK_LIMBS + 1];
  b9acc_t *acc;

  if (a->n == 0 || b->n == 0) { out->n = 0;  out->neg = 0;  return; }

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

  while (rn > 0 && out->d[rn-1] == 0) rn--;
  out->n = rn;
  out->neg = 0;
}

/* Make a read-only, non-owning view of a limb range.  The view must not be
 * passed as an output or to b9_free. */
static void b9_init_slice_view(b9_t *view, const b9_t *a,
                               uint32_t start, uint32_t len)
{
  b9_init(view);
  if (start >= a->n) return;
  if (len > a->n-start) len = a->n-start;
  view->d = a->d + start;
  view->alloc = len;
  view->n = len;
  while (view->n > 0 && view->d[view->n-1] == 0) view->n--;
}

/* Initialize a temporary using its fixed-capacity portion of a work buffer. */
static void b9_init_work(b9_t *x, b9limb_t **next, uint32_t capacity)
{
  x->d = *next;
  x->alloc = capacity;
  x->n = 0;
  x->neg = 0;
  *next += capacity;
}

static void b9_add_work_size(size_t *total, uint32_t capacity)
{
  if (*total > (size_t)MAX_SIZET - (size_t)capacity)
    croak("internal: b9 Karatsuba workspace too large");
  *total += capacity;
}

/* Subtract |b| from nonnegative a in place. */
static void b9_sub_abs_inplace(b9_t *a, const b9_t *b)
{
  MPUassert(!a->neg && !b->neg &&
            b9mag_cmp(a->d, a->n, b->d, b->n) >= 0,
            "b9 Karatsuba subtraction underflow");
  a->n = b9mag_sub(a->d, a->d, a->n, b->d, b->n);
  if (a->n == 1 && a->d[0] == 0) a->n = 0;
  a->neg = 0;
}

/* Add x at a whole-limb offset into an allocated, zero-padded result. */
static void b9_add_shifted(b9_t *out, uint32_t outn,
                           const b9_t *x, uint32_t offset)
{
  uint32_t carry = 0, i, pos;

  MPUassert(offset <= outn && x->n <= outn-offset,
            "b9 Karatsuba result overflow");
  for (i = 0; i < x->n; i++) {
    b9acc_t sum = (b9acc_t)out->d[offset+i] + x->d[i] + carry;
    if (sum >= B9_BASE) {
      out->d[offset+i] = (b9limb_t)(sum - B9_BASE);
      carry = 1;
    } else {
      out->d[offset+i] = (b9limb_t)sum;
      carry = 0;
    }
  }
  pos = offset + x->n;
  while (carry != 0) {
    b9acc_t sum;
    MPUassert(pos < outn, "b9 Karatsuba carry overflow");
    sum = (b9acc_t)out->d[pos] + carry;
    if (sum >= B9_BASE) {
      out->d[pos++] = (b9limb_t)(sum - B9_BASE);
    } else {
      out->d[pos++] = (b9limb_t)sum;
      carry = 0;
    }
  }
}

static void b9_mul_abs(b9_t *out, const b9_t *a, const b9_t *b);

static void b9_mul_karatsuba(b9_t *out, const b9_t *a, const b9_t *b)
{
  b9_t a0, a1, b0, b1, asum, bsum, z0, z1, z2;
  uint32_t an = a->n, bn = b->n;
  uint32_t maxn = an > bn ? an : bn;
  uint32_t split = maxn >> 1;
  uint32_t rn = an + bn;
  uint32_t a0cap = an < split ? an : split;
  uint32_t b0cap = bn < split ? bn : split;
  uint32_t a1cap = an-a0cap;
  uint32_t b1cap = bn-b0cap;
  uint32_t asumcap = (a0cap > a1cap ? a0cap : a1cap) + 1;
  uint32_t bsumcap = (b0cap > b1cap ? b0cap : b1cap) + 1;
  uint32_t z0cap = a0cap + b0cap;
  uint32_t z1cap = asumcap + bsumcap;
  uint32_t z2cap = a1cap + b1cap;
  size_t worklimbs = 0;
  b9limb_t *work, *next;

  b9_add_work_size(&worklimbs, asumcap);
  b9_add_work_size(&worklimbs, bsumcap);
  b9_add_work_size(&worklimbs, z0cap);
  b9_add_work_size(&worklimbs, z1cap);
  b9_add_work_size(&worklimbs, z2cap);
  work = (b9limb_t*)b9_xmalloc(worklimbs, sizeof(b9limb_t));
  next = work;

  b9_init_work(&asum, &next, asumcap);
  b9_init_work(&bsum, &next, bsumcap);
  b9_init_work(&z0, &next, z0cap);
  b9_init_work(&z1, &next, z1cap);
  b9_init_work(&z2, &next, z2cap);
  MPUassert(next == work + worklimbs,
            "b9 Karatsuba workspace mismatch");

  b9_init_slice_view(&a0, a, 0, split);
  b9_init_slice_view(&a1, a, split, an-split);
  b9_init_slice_view(&b0, b, 0, split);
  b9_init_slice_view(&b1, b, split, bn-split);

  b9_mul_abs(&z0, &a0, &b0);
  b9_mul_abs(&z2, &a1, &b1);
  b9_add(&asum, &a0, &a1);
  b9_add(&bsum, &b0, &b1);
  b9_mul_abs(&z1, &asum, &bsum);
  b9_sub_abs_inplace(&z1, &z0);
  b9_sub_abs_inplace(&z1, &z2);

  b9_ensure(out, rn);
  memset(out->d, 0, (size_t)rn * sizeof(b9limb_t));
  b9_add_shifted(out, rn, &z0, 0);
  b9_add_shifted(out, rn, &z1, split);
  b9_add_shifted(out, rn, &z2, split + split);
  out->n = rn;
  out->neg = 0;
  while (out->n > 0 && out->d[out->n-1] == 0) out->n--;

  free(work);
}

static void b9_mul_abs(b9_t *out, const b9_t *a, const b9_t *b)
{
  uint32_t minn = a->n < b->n ? a->n : b->n;
  uint32_t maxn = a->n > b->n ? a->n : b->n;

  /* Very unequal operands do not recover the extra additions and copies. */
  if (minn >= B9_KARATSUBA_THRESHOLD && maxn-minn <= minn/2)
    b9_mul_karatsuba(out, a, b);
  else
    b9_mul_schoolbook(out, a, b);
}

void b9_mul(b9_t *out, const b9_t *a, const b9_t *b)
{
  const b9_t *large;
  b9limb_t small;
  int neg;
#if HAVE_UINT64
  uint32_t small32;
#endif

  if (a->n == 0 || b->n == 0) {
    out->n = 0;
    out->neg = 0;
    return;
  }
  neg = a->neg != b->neg;
  if (a->n == 1 || b->n == 1) {
    large = a->n == 1 ? b : a;
    small = a->n == 1 ? a->d[0] : b->d[0];
    b9_mul_limb_abs(out, large, small);
    out->neg = neg && out->n != 0;
    return;
  }
#if HAVE_UINT64
  if (b9_get_u32_checked(&small32, a)) {
    b9_mul_u32(out, b, small32);
    out->neg = neg && out->n != 0;
    return;
  }
  if (b9_get_u32_checked(&small32, b)) {
    b9_mul_u32(out, a, small32);
    out->neg = neg && out->n != 0;
    return;
  }
#endif
  if (a->n > UINT32_MAX-b->n)
    croak("internal: b9 multiplication too large");
  b9_mul_abs(out, a, b);
  out->neg = neg && out->n != 0;
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
  uint32_t i, need, extra = (uint32_t)B9_NLIMBS(10) + 1;
  int neg = a->neg;

  if (v == 0 || a->n == 0) {
    out->n = 0;
    out->neg = 0;
    return;
  }

  if (a->n > UINT32_MAX-extra)
    croak("internal: b9 multiplication too large");
  need = a->n + extra;
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

/* Test bit k of |x| without modifying x. */
bool b9_testbit(const b9_t *x, uint32_t k)
{
  b9_t q;
  bool ret;

  if (b9_is_zero(x)) return 0;
  if (k < 4) return ((x->d[0] >> k) & 1) != 0;

  b9_init(&q);
  b9_tdiv_2exp(&q, x, (UV)k);
  b9_abs(&q);
  ret = !b9_is_even(&q);
  b9_free(&q);
  return ret;
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

#if HAVE_UINT64 /* 2^44 * 10^6 fits in the 64-bit accumulator. */
  #define B9_TDIV_2EXP_CHUNK_BITS  44
#else           /* 2^18 * 10^4 fits in the 32-bit accumulator. */
  #define B9_TDIV_2EXP_CHUNK_BITS  18
#endif

/* Divide by 2^bits using one pass over the decimal limbs.  The caller limits
 * bits so carry * B9_BASE fits in b9acc_t. */
static b9acc_t b9_tdiv_2exp_chunk(b9_t *a, uint32_t bits)
{
  b9acc_t carry = 0;
  b9acc_t mask = ((b9acc_t)1 << bits) - 1;
  uint32_t i;

  for (i = a->n; i-- > 0; ) {
    b9acc_t cur = carry * B9_BASE + a->d[i];
    a->d[i] = (b9limb_t)(cur >> bits);
    carry = cur & mask;
  }
  while (a->n > 0 && a->d[a->n-1] == 0) a->n--;
  if (a->n == 0) a->neg = 0;
  return carry;
}

/* Divide in place by 2^bits and report whether any discarded bit was set. */
static bool b9_tdiv_2exp_inplace(b9_t *a, UV bits)
{
  bool discarded = 0;
  if (bits == 1) {
    discarded = !b9_is_even(a);
    b9_tdiv2(a);
    return discarded;
  }
  if (bits == 4) {
    discarded = !b9_is_zero(a) && (a->d[0] & 15) != 0;
    b9_tdiv16(a);
    return discarded;
  }
  while (bits > 0 && !b9_is_zero(a)) {
    uint32_t chunk = bits > B9_TDIV_2EXP_CHUNK_BITS
                   ? B9_TDIV_2EXP_CHUNK_BITS : (uint32_t)bits;
    if (b9_tdiv_2exp_chunk(a, chunk) != 0) discarded = 1;
    bits -= chunk;
  }
  return discarded;
}

/* out = trunc(a / 2^bits).  out may alias a. */
void b9_tdiv_2exp(b9_t *out, const b9_t *a, UV bits)
{
  if (out != a) b9_set(out, a);
  (void)b9_tdiv_2exp_inplace(out, bits);
}

/* out = floor(a / 2^bits).  out may alias a. */
void b9_fdiv_2exp(b9_t *out, const b9_t *a, UV bits)
{
  bool neg = b9_is_negative(a);
  bool adjust;
  if (out != a) b9_set(out, a);
  adjust = b9_tdiv_2exp_inplace(out, bits);
  if (neg && adjust) {
    b9_abs(out);
    b9_add_u32(out, out, 1);
    b9_set_negative(out, 1);
  }
}

static uint32_t b9_ctz_acc(b9acc_t v)
{
  uint32_t count = 0;
  MPUassert(v != 0, "b9_ctz_acc called with zero");
  while (!(v & 1)) { v >>= 1; count++; }
  return count;
}

/* Find the first set bit, destructively scanning x from bit zero. */
static bool b9_scan1_from_zero(uint32_t *result, b9_t *x, uint32_t limit)
{
  uint32_t count = 0;

  while (!b9_is_zero(x)) {
    b9acc_t rem = b9_tdiv_2exp_chunk(x, B9_TDIV_2EXP_CHUNK_BITS);
    if (rem != 0) {
      uint32_t skip = b9_ctz_acc(rem);
      if (skip > limit-count) return 0;
      *result = count + skip;
      return 1;
    }
    if (B9_TDIV_2EXP_CHUNK_BITS > limit-count) return 0;
    count += B9_TDIV_2EXP_CHUNK_BITS;
  }
  return 0;
}

static bool b9_bitscan(uint32_t *result, const b9_t *x, uint32_t start,
                       bool want_one)
{
  b9_t q;
  uint32_t skip;
  bool found;

  b9_init_set(&q, x);
  b9_abs(&q);
  if (start != 0) b9_tdiv_2exp(&q, &q, (UV)start);
  if (!want_one) b9_add_u32(&q, &q, 1);
  found = b9_scan1_from_zero(&skip, &q, UINT32_MAX-start);
  if (found) *result = start + skip;
  b9_free(&q);
  return found;
}

bool b9_bitscan0(uint32_t *result, const b9_t *x, uint32_t start)
{
  return b9_bitscan(result, x, start, 0);
}

bool b9_bitscan1(uint32_t *result, const b9_t *x, uint32_t start)
{
  return b9_bitscan(result, x, start, 1);
}

#if B9_DIGS == 9
  #define B9_BINARY_LIMB_BITS  30
#elif B9_DIGS == 6
  #define B9_BINARY_LIMB_BITS  20
#else
  #define B9_BINARY_LIMB_BITS  14
#endif
#if B9_BINARY_WORD_BITS == 32
  #define B9_BINARY_WORD_MASK  UINT32_MAX
#else
  #define B9_BINARY_WORD_MASK  UINT32_C(65535)
#endif

/* Upper bound on the binary words needed by a. */
static size_t b9_binary_word_bound(const b9_t *a)
{
  size_t bits;
  if ((size_t)a->n > ((size_t)MAX_SIZET - B9_BINARY_WORD_BITS + 1) /
                     B9_BINARY_LIMB_BITS)
    croak("internal: b9 binary value too large");
  bits = (size_t)a->n * B9_BINARY_LIMB_BITS;
  return (bits + B9_BINARY_WORD_BITS - 1) / B9_BINARY_WORD_BITS;
}

#define B9_BINARY_DC_THRESHOLD  16

static unsigned b9_size_log2(size_t n)
{
  unsigned log = 0;
  while (n >>= 1) log++;
  return log;
}

static void b9_set_binary_horner(b9_t *out, const b9binary_t *words,
                                 size_t n)
{
  b9_set_uv(out, 0);
  while (n-- > 0) {
    b9_mul_2exp(out, out, B9_BINARY_WORD_BITS);
    b9_add_u32(out, out, words[n]);
  }
}

/* Convert binary words using powers[k] = 2^(word_bits * 2^k). */
static void b9_set_binary_dc(b9_t *out, const b9binary_t *words, size_t n,
                             const b9_t powers[])
{
  size_t low_n;
  unsigned power_index;
  b9_t hi, lo, tmp;

  if (n < B9_BINARY_DC_THRESHOLD) {
    b9_set_binary_horner(out, words, n);
    return;
  }
  power_index = b9_size_log2(n / 2);
  low_n = (size_t)1 << power_index;
  b9_init(&hi);
  b9_init(&lo);
  b9_init(&tmp);
  b9_set_binary_dc(&hi, words+low_n, n-low_n, powers);
  b9_set_binary_dc(&lo, words, low_n, powers);
  b9_mul(&tmp, &hi, &powers[power_index]);
  b9_add(out, &tmp, &lo);
  b9_free(&hi);
  b9_free(&lo);
  b9_free(&tmp);
}

void b9_set_binary(b9_t *out, const b9binary_t *words, size_t n)
{
  b9_t *powers;
  size_t i, npowers;

  while (n > 0 && words[n-1] == 0) n--;
  if (n < B9_BINARY_DC_THRESHOLD) {
    b9_set_binary_horner(out, words, n);
    return;
  }

  npowers = (size_t)b9_size_log2(n / 2) + 1;
  powers = (b9_t*)b9_xmalloc(npowers, sizeof(b9_t));
  for (i = 0; i < npowers; i++) b9_init(&powers[i]);
  b9_init_set_pow2(&powers[0], B9_BINARY_WORD_BITS);
  for (i = 1; i < npowers; i++)
    b9_mul(&powers[i], &powers[i-1], &powers[i-1]);
  b9_set_binary_dc(out, words, n, powers);
  for (i = 0; i < npowers; i++) b9_free(&powers[i]);
  free(powers);
}

enum {
  B9_BITOP_AND,
  B9_BITOP_OR,
  B9_BITOP_XOR,
  B9_BITOP_ANDNOT
};

static void b9_bitop(b9_t *out, const b9_t *a, const b9_t *b, int op)
{
  b9_t A, B;
  b9binary_t *words;
  size_t awords, bwords, alloc, n = 0;

  if (op == B9_BITOP_AND && (b9_is_zero(a) || b9_is_zero(b))) {
    b9_set_uv(out, 0);
    return;
  }
  if (op == B9_BITOP_ANDNOT && b9_is_zero(a)) {
    b9_set_uv(out, 0);
    return;
  }
  if ((op == B9_BITOP_OR || op == B9_BITOP_XOR) && b9_is_zero(a)) {
    b9_set(out, b);  b9_abs(out);
    return;
  }
  if ((op == B9_BITOP_OR || op == B9_BITOP_XOR ||
       op == B9_BITOP_ANDNOT) && b9_is_zero(b)) {
    b9_set(out, a);  b9_abs(out);
    return;
  }

  awords = b9_binary_word_bound(a);
  bwords = b9_binary_word_bound(b);
  if (op == B9_BITOP_AND)
    alloc = awords < bwords ? awords : bwords;
  else if (op == B9_BITOP_ANDNOT)
    alloc = awords;
  else
    alloc = awords > bwords ? awords : bwords;
  words = (b9binary_t*)b9_xcalloc(alloc, sizeof(b9binary_t));
  b9_init_set(&A, a);  b9_abs(&A);
  b9_init_set(&B, b);  b9_abs(&B);

  while (n < alloc) {
    uint32_t av, bv;
    if (op == B9_BITOP_AND && (b9_is_zero(&A) || b9_is_zero(&B)))
      break;
    if (op == B9_BITOP_ANDNOT && b9_is_zero(&A))
      break;
    if (op != B9_BITOP_AND && op != B9_BITOP_ANDNOT &&
        b9_is_zero(&A) && b9_is_zero(&B))
      break;
    av = b9_is_zero(&A) ? 0 :
         (uint32_t)b9_tdiv_2exp_chunk(&A, B9_BINARY_WORD_BITS);
    bv = b9_is_zero(&B) ? 0 :
         (uint32_t)b9_tdiv_2exp_chunk(&B, B9_BINARY_WORD_BITS);
    switch (op) {
      case B9_BITOP_AND:    words[n] = av & bv;  break;
      case B9_BITOP_OR:     words[n] = av | bv;  break;
      case B9_BITOP_XOR:    words[n] = av ^ bv;  break;
      default:              words[n] = av & ~bv; break;
    }
    n++;
  }

  b9_set_binary(out, words, n);
  b9_free(&A);  b9_free(&B);
  free(words);
}

void b9_bitand(b9_t *out, const b9_t *a, const b9_t *b)
{
  b9_bitop(out, a, b, B9_BITOP_AND);
}

void b9_bitor(b9_t *out, const b9_t *a, const b9_t *b)
{
  b9_bitop(out, a, b, B9_BITOP_OR);
}

void b9_bitxor(b9_t *out, const b9_t *a, const b9_t *b)
{
  b9_bitop(out, a, b, B9_BITOP_XOR);
}

void b9_bitandnot(b9_t *out, const b9_t *a, const b9_t *b)
{
  b9_bitop(out, a, b, B9_BITOP_ANDNOT);
}

void b9_bitnot(b9_t *out, const b9_t *a, bool fixed_width, uint32_t width)
{
  b9_t A;
  b9binary_t *words;
  uint32_t topword = 0;
  size_t alloc, n = 0;

  if (!fixed_width && b9_is_zero(a)) {
    b9_set_uv(out, 1);
    return;
  }
  if (fixed_width && width == 0) {
    b9_set_uv(out, 0);
    return;
  }

  alloc = fixed_width
        ? (size_t)(width / B9_BINARY_WORD_BITS) +
          (width % B9_BINARY_WORD_BITS != 0)
        : b9_binary_word_bound(a);
  words = (b9binary_t*)b9_xcalloc(alloc, sizeof(b9binary_t));
  b9_init_set(&A, a);  b9_abs(&A);

  while (n < alloc) {
    uint32_t av = b9_is_zero(&A) ? 0 :
      (uint32_t)b9_tdiv_2exp_chunk(&A, B9_BINARY_WORD_BITS);
    words[n++] = (~av) & B9_BINARY_WORD_MASK;
    topword = av;
    if (!fixed_width && b9_is_zero(&A)) break;
  }

  if (fixed_width) {
    uint32_t topbits = width % B9_BINARY_WORD_BITS;
    if (topbits != 0)
      words[n-1] &= (UINT32_C(1) << topbits) - 1;
  } else {
    uint32_t topbits = 0;
    while (topword != 0) { topword >>= 1; topbits++; }
    if (topbits < B9_BINARY_WORD_BITS)
      words[n-1] &= (UINT32_C(1) << topbits) - 1;
  }

  b9_set_binary(out, words, n);
  b9_free(&A);
  free(words);
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

/* Divide nonnegative a in place by p and return the remainder. */
uint32_t b9_divrem_u32_inplace(b9_t* a, uint32_t p)
{
  if (p == 2 || p == 4 || p == 8 || p == 16) {
    uint32_t bits = (p == 2) ? 1 : (p == 4) ? 2 : (p == 8) ? 3 : 4;
    uint32_t rem = (a->n == 0) ? 0 : (a->d[0] & (p-1));
    b9_tdiv_2exp(a, a, (UV)bits);
    return rem;
  }

#if HAVE_UINT64
  uint32_t rem = 0;
  uint32_t i;

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
