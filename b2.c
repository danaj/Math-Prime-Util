#include <stdlib.h>
#include <string.h>

#include "b2.h"
#include "util_bits.h"

/* Use schoolbook multiplication below this many limbs (must be >= 3). */
#define B2_KARATSUBA_THRESHOLD     28
#define B2_SQR_KARATSUBA_THRESHOLD 56

/******************************************************************************/
/* Allocation and ownership. */
/******************************************************************************/

static void* b2_xmalloc(size_t count, size_t size)
{
  void *p;
  if (size != 0 && count > (size_t)MAX_SIZET / size)
    croak("internal: b2 allocation too large");
  p = malloc(count * size);
  if (p == 0 && count != 0 && size != 0)
    croak("internal: b2 allocation failed");
  return p;
}

static void* b2_xrealloc(void *ptr, size_t count, size_t size)
{
  void *p;
  if (size != 0 && count > (size_t)MAX_SIZET / size)
    croak("internal: b2 allocation too large");
  p = realloc(ptr, count * size);
  if (p == 0 && count != 0 && size != 0)
    croak("internal: b2 allocation failed");
  return p;
}

static void b2_trim(b2_t *x)
{
  while (x->n > 0 && x->d[x->n-1] == 0) x->n--;
  if (x->n == 0) x->neg = 0;
}

void b2_free(b2_t *x)
{
  if (x->d != x->d_small) free(x->d);
  b2_init(x);
}

void b2_neg(b2_t *x)
{
  x->neg = x->n != 0 && !x->neg;
}

void b2_ensure(b2_t *x, uint32_t need)
{
  b2limb_t *newd;
  if (x->alloc >= need) return;
  if (x->d == x->d_small) {
    newd = (b2limb_t*)b2_xmalloc((size_t)need, sizeof(b2limb_t));
    if (x->n != 0)
      memcpy(newd, x->d_small, (size_t)x->n * sizeof(b2limb_t));
  } else {
    newd = (b2limb_t*)b2_xrealloc(x->d, (size_t)need,
                                  sizeof(b2limb_t));
  }
  x->d = newd;
  x->alloc = need;
}

void b2_set(b2_t *dst, const b2_t *src)
{
  if (dst == src) return;
  b2_ensure(dst, src->n);
  if (src->n != 0)
    memcpy(dst->d, src->d, (size_t)src->n * sizeof(b2limb_t));
  dst->n = src->n;
  dst->neg = src->neg;
}

void b2_move(b2_t *dst, b2_t *src)
{
  if (dst == src) return;
  b2_free(dst);
  if (src->d == src->d_small) {
    if (src->n != 0)
      memcpy(dst->d, src->d, (size_t)src->n * sizeof(b2limb_t));
  } else {
    dst->d = src->d;
    dst->alloc = src->alloc;
  }
  dst->n = src->n;
  dst->neg = src->neg;
  b2_init(src);
}

void b2_swap(b2_t *a, b2_t *b)
{
  bool a_inline, b_inline, neg;
  b2limb_t small[B2_INLINE_LIMBS];
  b2limb_t *d;
  b2_t *t;
  uint32_t an, bn, alloc, n;

  if (a == b) return;
  a_inline = a->d == a->d_small;
  b_inline = b->d == b->d_small;
  if (a_inline && b_inline) {
    an = a->n;  bn = b->n;
    if (an != 0)
      memcpy(small, a->d, (size_t)an * sizeof(b2limb_t));
    if (bn != 0)
      memcpy(a->d, b->d, (size_t)bn * sizeof(b2limb_t));
    if (an != 0)
      memcpy(b->d, small, (size_t)an * sizeof(b2limb_t));
    neg = a->neg;  a->neg = b->neg;  b->neg = neg;
    a->n = bn;  b->n = an;
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
  an = a->n;
  if (an != 0)
    memcpy(small, a->d, (size_t)an * sizeof(b2limb_t));
  a->d = b->d;  a->alloc = b->alloc;  a->n = b->n;
  neg = a->neg;  a->neg = b->neg;
  b->d = b->d_small;  b->alloc = B2_INLINE_LIMBS;  b->n = an;
  b->neg = neg;
  if (an != 0)
    memcpy(b->d, small, (size_t)an * sizeof(b2limb_t));
}

/******************************************************************************/
/* Native and decimal conversion. */
/******************************************************************************/

void b2_set_uv(b2_t *x, UV v)
{
  uint32_t n = 0;
  b2_ensure(x, (uint32_t)((BITS_PER_WORD + B2_LIMB_BITS - 1) /
                          B2_LIMB_BITS));
  x->neg = 0;
  while (v != 0) {
    x->d[n++] = (b2limb_t)v;
#if B2_LIMB_BITS < BITS_PER_WORD
    v >>= B2_LIMB_BITS;
#else
    v = 0;
#endif
  }
  x->n = n;
}

UV b2_get_uv(const b2_t *x)
{
  UV v = 0;
#if B2_LIMB_BITS == BITS_PER_WORD
  if (x->n != 0) v = (UV)x->d[0];
#else
  uint32_t i;
  for (i = x->n; i-- > 0; )
    v = (v << B2_LIMB_BITS) | (UV)x->d[i];
#endif
  return v;
}

#ifndef B2_WIDE_INPUT
  #define B2_WIDE_INPUT  (HAVE_UINT128 && B2_LIMB_BITS == 32)
#endif

#if B2_WIDE_INPUT
  typedef uint64_t b2input_t;
  #define B2_INPUT_BASE    UINT64_C(10000000000000000000)
  #define B2_INPUT_DIGITS  19

static void b2_muladd_input_inplace(b2_t *x, b2input_t mul, b2input_t add)
{
  uint128_t carry = add;
  uint32_t i;
  if (x->n > UINT32_MAX-2)
    croak("internal: b2 decimal input too large");
  b2_ensure(x, x->n + 2);
  for (i = 0; i < x->n; i++) {
    uint128_t t = (uint128_t)x->d[i] * mul + carry;
    x->d[i] = (b2limb_t)t;
    carry = t >> B2_LIMB_BITS;
  }
  while (carry != 0) {
    x->d[x->n++] = (b2limb_t)carry;
    carry >>= B2_LIMB_BITS;
  }
}
#else
  typedef b2limb_t b2input_t;
  #define B2_INPUT_BASE    B2_DEC_BASE
  #define B2_INPUT_DIGITS  B2_DEC_DIGITS

static void b2_mul_limb_inplace(b2_t *x, b2limb_t v)
{
  b2acc_t carry = 0;
  uint32_t i;
  if (v == 0 || x->n == 0) {
    x->n = 0;
    x->neg = 0;
    return;
  }
  if (x->n == UINT32_MAX)
    croak("internal: b2 multiplication too large");
  b2_ensure(x, x->n + 1);
  for (i = 0; i < x->n; i++) {
    b2acc_t t = (b2acc_t)x->d[i] * v + carry;
    x->d[i] = (b2limb_t)t;
    carry = t >> B2_LIMB_BITS;
  }
  if (carry != 0) x->d[x->n++] = (b2limb_t)carry;
}

static void b2_add_limb_inplace(b2_t *x, b2limb_t v)
{
  b2acc_t carry = v;
  uint32_t i = 0;
  if (carry == 0) return;
  if (x->n == UINT32_MAX)
    croak("internal: b2 addition too large");
  b2_ensure(x, x->n + 1);
  while (carry != 0 && i < x->n) {
    b2acc_t t = (b2acc_t)x->d[i] + carry;
    x->d[i++] = (b2limb_t)t;
    carry = t >> B2_LIMB_BITS;
  }
  if (carry != 0) x->d[x->n++] = (b2limb_t)carry;
}
#endif

/* Parse unsigned decimal digits using the linear Horner method. */
static void b2_set_decimal_horner(b2_t *x, const char *s, STRLEN len)
{
  size_t estimate;
  STRLEN first, pos;

  estimate = (size_t)len / (B2_LIMB_BITS / 4) +
             ((size_t)len % (B2_LIMB_BITS / 4) != 0);
  if (estimate > UINT32_MAX)
    croak("internal: b2 decimal input too large");
  b2_ensure(x, (uint32_t)estimate);
  x->n = 0;
  x->neg = 0;

  first = len % B2_INPUT_DIGITS;
  if (first == 0) first = B2_INPUT_DIGITS;
  pos = 0;
  while (pos < len) {
    STRLEN end = pos + (pos == 0 ? first : B2_INPUT_DIGITS);
    b2input_t v = 0;
    while (pos < end) v = v * 10 + (uint32_t)(s[pos++] - '0');
#if B2_WIDE_INPUT
    b2_muladd_input_inplace(x, B2_INPUT_BASE, v);
#else
    if (x->n != 0)
      b2_mul_limb_inplace(x, (b2limb_t)B2_DEC_BASE);
    b2_add_limb_inplace(x, (b2limb_t)v);
#endif
  }
}

void b2_set_str(b2_t *x, const char *s, STRLEN len)
{
  bool neg;

  neg = len > 0 && s[0] == '-';
  if (len > 0 && (s[0] == '-' || s[0] == '+')) { s++; len--; }
  while (len > 1 && s[0] == '0') { s++; len--; }
  if (len == 0 || (len == 1 && s[0] == '0')) {
    x->n = 0;
    x->neg = 0;
    return;
  }
  b2_set_decimal_horner(x, s, len);
  x->neg = neg && x->n != 0;
}

static uint32_t b2_radix_digit(unsigned char c)
{
  if (c >= '0' && c <= '9') return (uint32_t)(c - '0');
  if (c >= 'A' && c <= 'F') return (uint32_t)(c - 'A' + 10);
  if (c >= 'a' && c <= 'f') return (uint32_t)(c - 'a' + 10);
  return UINT32_MAX;
}

void b2_set_radix(b2_t *x, const char *s, STRLEN len, uint32_t base)
{
  size_t nbits, nlimbs_size;
  STRLEN i;
  uint32_t bit = 0, shift, word = 0;

  MPUassert(base == 2 || base == 8 || base == 16,
            "b2_set_radix bad base");
  while (len > 0 && *s == '0') { s++; len--; }
  if (len == 0) {
    x->n = 0;
    x->neg = 0;
    return;
  }
  shift = base == 16 ? 4 : base == 8 ? 3 : 1;
  if ((size_t)len > (size_t)MAX_SIZET / shift)
    croak("internal: b2 radix input too large");
  nbits = (size_t)len * shift;
  nlimbs_size = nbits / B2_LIMB_BITS +
                (nbits % B2_LIMB_BITS != 0);
  if (nlimbs_size > UINT32_MAX)
    croak("internal: b2 radix input too large");
  b2_ensure(x, (uint32_t)nlimbs_size);
  memset(x->d, 0, nlimbs_size * sizeof(b2limb_t));

  for (i = len; i-- > 0; ) {
    uint32_t d = b2_radix_digit((unsigned char)s[i]);
    if (d >= base)
      croak("internal: b2 invalid radix digit");
    x->d[word] |= (b2limb_t)((b2limb_t)d << bit);
    if (bit + shift > B2_LIMB_BITS)
      x->d[word+1] |= (b2limb_t)(d >> (B2_LIMB_BITS-bit));
    bit += shift;
    if (bit >= B2_LIMB_BITS) {
      bit -= B2_LIMB_BITS;
      word++;
    }
  }
  x->n = (uint32_t)nlimbs_size;
  x->neg = 0;
  b2_trim(x);
}

STRLEN b2_string_size(const b2_t *x)
{
  size_t len;
  if (x->n == 0) return 1;
  if ((size_t)x->n > ((size_t)MAX_SIZET - 2) /
                     B2_DEC_BOUND_DIGITS)
    croak("internal: b2 decimal output too large");
  len = (size_t)x->n * B2_DEC_BOUND_DIGITS + (x->neg ? 1 : 0);
  return (STRLEN)len;
}

static STRLEN b2_write_decimal(char *buf, uint32_t v, STRLEN width)
{
  char tmp[12];
  STRLEN i = width;
  do { tmp[--i] = '0' + (char)(v % 10); v /= 10; } while (v != 0);
  while (i > 0) tmp[--i] = '0';
  memcpy(buf, tmp, (size_t)width);
  return width;
}

static STRLEN b2_get_str_impl(char *buf, const b2_t *x,
                              bool limited, STRLEN limit)
{
  b2limb_t small_chunks[16];
  b2limb_t *chunks;
  size_t actual, maxchars, maxchunks, nchunks = 0, topdigits;
  uint32_t i;
  STRLEN pos = 0;

  if (x->n == 0) {
    if (limited && limit < 1) return 0;
    buf[0] = '0';
    return 1;
  }
  maxchars = (size_t)b2_string_size(x);
  maxchunks = maxchars / B2_DEC_DIGITS +
              (maxchars % B2_DEC_DIGITS != 0);
  chunks = maxchunks <= sizeof(small_chunks)/sizeof(*small_chunks)
         ? small_chunks
         : (b2limb_t*)b2_xmalloc(maxchunks, sizeof(b2limb_t));

  /* Convert from base 2^B2_LIMB_BITS to base B2_DEC_BASE.  Unlike repeated
   * division of the entire binary value, each pass uses only division by a
   * compile-time constant. */
  for (i = x->n; i-- > 0; ) {
    b2acc_t carry = x->d[i];
    size_t j;
    for (j = 0; j < nchunks; j++) {
      b2acc_t cur = ((b2acc_t)chunks[j] << B2_LIMB_BITS) + carry;
      carry = cur / B2_DEC_BASE;
      chunks[j] = (b2limb_t)(cur - carry * B2_DEC_BASE);
    }
    while (carry != 0) {
      b2acc_t next = carry / B2_DEC_BASE;
      MPUassert(nchunks < maxchunks, "b2_get_str chunk overflow");
      chunks[nchunks++] = (b2limb_t)(carry - next * B2_DEC_BASE);
      carry = next;
    }
  }

  {
    uint32_t v = (uint32_t)chunks[nchunks-1];
    for (topdigits = 1; v >= 10; topdigits++) v /= 10;
  }
  actual = (x->neg ? 1 : 0) + topdigits;
  if (nchunks-1 > ((size_t)MAX_SIZET-actual) / B2_DEC_DIGITS) {
    if (chunks != small_chunks) free(chunks);
    croak("internal: b2 decimal output too large");
  }
  actual += (nchunks-1) * B2_DEC_DIGITS;
  if (limited && actual > (size_t)limit) {
    if (chunks != small_chunks) free(chunks);
    return 0;
  }

  if (x->neg) buf[pos++] = '-';
  {
    uint32_t v = (uint32_t)chunks[--nchunks];
    char tmp[12];
    STRLEN i = 0, j;
    do { tmp[i++] = '0' + (char)(v % 10); v /= 10; } while (v != 0);
    for (j = i; j-- > 0; ) buf[pos++] = tmp[j];
  }
  while (nchunks > 0) {
    pos += b2_write_decimal(buf + pos, (uint32_t)chunks[--nchunks],
                            B2_DEC_DIGITS);
  }
  if (chunks != small_chunks) free(chunks);
  return pos;
}

STRLEN b2_get_str(char *buf, const b2_t *x)
{
  return b2_get_str_impl(buf, x, 0, 0);
}

STRLEN b2_get_str_limit(char *buf, const b2_t *x, STRLEN limit)
{
  return b2_get_str_impl(buf, x, 1, limit);
}

void b2_init_set(b2_t *x, const b2_t *y)
{
  b2_init(x);
  b2_set(x, y);
}

void b2_init_set_uv(b2_t *x, UV v)
{
  b2_init(x);
  b2_set_uv(x, v);
}

void b2_init_set_str(b2_t *x, const char *s, STRLEN len)
{
  b2_init(x);
  b2_set_str(x, s, len);
}

void b2_init_set_radix(b2_t *x, const char *s, STRLEN len, uint32_t base)
{
  b2_init(x);
  b2_set_radix(x, s, len, base);
}

void b2_init_set_pow2(b2_t *x, uint32_t k)
{
  b2_init(x);
  b2_setbit(x, k);
}

/******************************************************************************/
/* Comparison and shifts. */
/******************************************************************************/

int b2_cmp_abs(const b2_t *a, const b2_t *b)
{
  uint32_t i;
  if (a->n != b->n) return a->n < b->n ? -1 : 1;
  for (i = a->n; i-- > 0; )
    if (a->d[i] != b->d[i]) return a->d[i] < b->d[i] ? -1 : 1;
  return 0;
}

int b2_cmp(const b2_t *a, const b2_t *b)
{
  int cmp;
  if (a->neg != b->neg) return a->neg ? -1 : 1;
  cmp = b2_cmp_abs(a, b);
  return a->neg ? -cmp : cmp;
}

static uint32_t b2mag_add(b2limb_t *out,
                          const b2limb_t *a, uint32_t an,
                          const b2limb_t *b, uint32_t bn)
{
  b2acc_t carry = 0;
  uint32_t i, rn = an > bn ? an : bn;
  for (i = 0; i < rn; i++) {
    b2acc_t av = i < an ? a[i] : 0;
    b2acc_t bv = i < bn ? b[i] : 0;
    b2acc_t sum = av + bv + carry;
    out[i] = (b2limb_t)sum;
    carry = sum >> B2_LIMB_BITS;
  }
  if (carry != 0) out[rn++] = (b2limb_t)carry;
  return rn;
}

/* |a| >= |b|. */
static uint32_t b2mag_sub(b2limb_t *out,
                          const b2limb_t *a, uint32_t an,
                          const b2limb_t *b, uint32_t bn)
{
  uint32_t borrow = 0, i;
  for (i = 0; i < an; i++) {
    b2limb_t av = a[i];
    b2limb_t bv = i < bn ? b[i] : 0;
    out[i] = (b2limb_t)(av - bv - borrow);
    borrow = av < bv || (borrow != 0 && av == bv);
  }
  while (an > 0 && out[an-1] == 0) an--;
  return an;
}

void b2_add(b2_t *out, const b2_t *a, const b2_t *b)
{
  const b2_t *hi, *lo;
  uint32_t need, rn;
  bool neg;
  int cmp;

  if (a->n == 0) { b2_set(out, b); return; }
  if (b->n == 0) { b2_set(out, a); return; }
  if (a->neg == b->neg) {
    need = a->n > b->n ? a->n : b->n;
    if (need == UINT32_MAX)
      croak("internal: b2 addition too large");
    neg = a->neg;
    b2_ensure(out, need + 1);
    rn = b2mag_add(out->d, a->d, a->n, b->d, b->n);
    out->n = rn;
    out->neg = neg;
    return;
  }

  cmp = b2_cmp_abs(a, b);
  if (cmp == 0) {
    out->n = 0;
    out->neg = 0;
    return;
  }
  hi = cmp > 0 ? a : b;
  lo = cmp > 0 ? b : a;
  neg = hi->neg;
  b2_ensure(out, hi->n);
  rn = b2mag_sub(out->d, hi->d, hi->n, lo->d, lo->n);
  out->n = rn;
  out->neg = neg;
}

/* Multiply a magnitude by one limb.  out may alias a. */
static void b2_mul_limb_abs(b2_t *out, const b2_t *a, b2limb_t v)
{
  b2acc_t carry = 0;
  uint32_t i;

  if (v == 0 || a->n == 0) {
    out->n = 0;
    out->neg = 0;
    return;
  }
  if (v == 1) {
    b2_set(out, a);
    out->neg = 0;
    return;
  }
  if (a->n == UINT32_MAX)
    croak("internal: b2 multiplication too large");
  b2_ensure(out, a->n + 1);
  for (i = 0; i < a->n; i++) {
    b2acc_t product = (b2acc_t)a->d[i] * v + carry;
    out->d[i] = (b2limb_t)product;
    carry = product >> B2_LIMB_BITS;
  }
  out->n = a->n;
  if (carry != 0) out->d[out->n++] = (b2limb_t)carry;
  out->neg = 0;
}

/* out must not alias either input. */
static void b2_mul_schoolbook(b2_t *out, const b2_t *a, const b2_t *b)
{
  const b2_t *tmp;
  b2acc_t carry;
  uint32_t i, j, rn;

  if (a->n == 0 || b->n == 0) {
    out->n = 0;
    out->neg = 0;
    return;
  }
  if (a->n < b->n) {
    tmp = a;
    a = b;
    b = tmp;
  }
  rn = a->n + b->n;
  b2_ensure(out, rn);

  carry = 0;
  for (j = 0; j < b->n; j++) {
    b2acc_t cur = (b2acc_t)a->d[0] * b->d[j] + carry;
    out->d[j] = (b2limb_t)cur;
    carry = cur >> B2_LIMB_BITS;
  }
  out->d[b->n] = (b2limb_t)carry;

  for (i = 1; i < a->n; i++) {
    if (a->d[i] == 0) {
      out->d[i+b->n] = 0;
      continue;
    }
    carry = 0;
    for (j = 0; j < b->n; j++) {
      b2acc_t cur = (b2acc_t)a->d[i] * b->d[j] +
                    out->d[i+j] + carry;
      out->d[i+j] = (b2limb_t)cur;
      carry = cur >> B2_LIMB_BITS;
    }
    out->d[i+b->n] = (b2limb_t)carry;
  }
  out->n = rn;
  out->neg = 0;
  b2_trim(out);
}

/* Add a double-limb product at a whole-limb offset. */
static void b2_add_product_at(b2_t *out, uint32_t outn,
                              uint32_t pos, b2acc_t product)
{
  b2acc_t sum, carry;

  MPUassert(pos < outn-1, "b2 square product overflow");
  sum = (b2acc_t)out->d[pos] + (b2limb_t)product;
  out->d[pos++] = (b2limb_t)sum;
  carry = sum >> B2_LIMB_BITS;

  sum = (b2acc_t)out->d[pos] +
        (b2limb_t)(product >> B2_LIMB_BITS) + carry;
  out->d[pos++] = (b2limb_t)sum;
  carry = sum >> B2_LIMB_BITS;
  while (carry != 0 && pos < outn) {
    sum = (b2acc_t)out->d[pos] + carry;
    out->d[pos++] = (b2limb_t)sum;
    carry = sum >> B2_LIMB_BITS;
  }
  MPUassert(carry == 0, "b2 square carry overflow");
}

/* out must not alias a. */
static void b2_sqr_schoolbook(b2_t *out, const b2_t *a)
{
  b2acc_t carry;
  uint32_t i, j, rn = a->n + a->n;

  b2_ensure(out, rn);
  memset(out->d, 0, (size_t)rn * sizeof(b2limb_t));

  /* Accumulate each off-diagonal product once. */
  for (i = 0; i < a->n; i++) {
    carry = 0;
    for (j = i+1; j < a->n; j++) {
      b2acc_t cur = (b2acc_t)a->d[i] * a->d[j] +
                    out->d[i+j] + carry;
      out->d[i+j] = (b2limb_t)cur;
      carry = cur >> B2_LIMB_BITS;
    }
    out->d[i+a->n] = (b2limb_t)carry;
  }

  /* Double the off-diagonal sum, then add each diagonal product. */
  carry = 0;
  for (i = 0; i < rn; i++) {
    b2acc_t cur = ((b2acc_t)out->d[i] << 1) + carry;
    out->d[i] = (b2limb_t)cur;
    carry = cur >> B2_LIMB_BITS;
  }
  MPUassert(carry == 0, "b2 square shift overflow");
  for (i = 0; i < a->n; i++)
    b2_add_product_at(out, rn, i+i,
                      (b2acc_t)a->d[i] * a->d[i]);

  out->n = rn;
  out->neg = 0;
  b2_trim(out);
}

/* Make a read-only, non-owning view of a limb range.  The view must not be
 * passed as an output or to b2_free. */
static void b2_init_slice_view(b2_t *view, const b2_t *a,
                               uint32_t start, uint32_t len)
{
  b2_init(view);
  if (start >= a->n) return;
  if (len > a->n-start) len = a->n-start;
  view->d = a->d + start;
  view->alloc = len;
  view->n = len;
  while (view->n > 0 && view->d[view->n-1] == 0) view->n--;
}

/* Initialize a temporary using its fixed-capacity portion of a work buffer. */
static void b2_init_work(b2_t *x, b2limb_t **next, uint32_t capacity)
{
  x->d = *next;
  x->alloc = capacity;
  x->n = 0;
  x->neg = 0;
  *next += capacity;
}

static void b2_add_work_size(size_t *total, uint32_t capacity)
{
  if (*total > (size_t)MAX_SIZET - (size_t)capacity)
    croak("internal: b2 Karatsuba workspace too large");
  *total += capacity;
}

/* Subtract |b| from the nonnegative a in place. */
static void b2_sub_abs_inplace(b2_t *a, const b2_t *b)
{
  MPUassert(!a->neg && !b->neg && b2_cmp_abs(a, b) >= 0,
            "b2 Karatsuba subtraction underflow");
  a->n = b2mag_sub(a->d, a->d, a->n, b->d, b->n);
  a->neg = 0;
}

/* Add x at a whole-limb offset into an allocated, zero-padded result. */
static void b2_add_shifted(b2_t *out, uint32_t outn,
                           const b2_t *x, uint32_t offset)
{
  b2acc_t carry = 0;
  uint32_t i, pos;

  MPUassert(offset <= outn && x->n <= outn-offset,
            "b2 Karatsuba result overflow");
  for (i = 0; i < x->n; i++) {
    b2acc_t sum = (b2acc_t)out->d[offset+i] + x->d[i] + carry;
    out->d[offset+i] = (b2limb_t)sum;
    carry = sum >> B2_LIMB_BITS;
  }
  pos = offset + x->n;
  while (carry != 0) {
    b2acc_t sum;
    MPUassert(pos < outn, "b2 Karatsuba carry overflow");
    sum = (b2acc_t)out->d[pos] + carry;
    out->d[pos++] = (b2limb_t)sum;
    carry = sum >> B2_LIMB_BITS;
  }
}

static void b2_mul_abs(b2_t *out, const b2_t *a, const b2_t *b);
static void b2_sqr_abs(b2_t *out, const b2_t *a);

static void b2_mul_karatsuba(b2_t *out, const b2_t *a, const b2_t *b)
{
  b2_t a0, a1, b0, b1, asum, bsum, z0, z1, z2;
  uint32_t an = a->n, bn = b->n;
  uint32_t maxn = an > bn ? an : bn;
  uint32_t split = maxn >> 1;
  uint32_t rn = an + bn;
  uint32_t a0cap = an < split ? an : split;
  uint32_t b0cap = bn < split ? bn : split;
  uint32_t a1cap = an - a0cap;
  uint32_t b1cap = bn - b0cap;
  uint32_t asumcap = (a0cap > a1cap ? a0cap : a1cap) + 1;
  uint32_t bsumcap = (b0cap > b1cap ? b0cap : b1cap) + 1;
  uint32_t z0cap = a0cap + b0cap;
  uint32_t z2cap = a1cap + b1cap;
  uint32_t z1cap = asumcap + bsumcap;
  size_t worklimbs = 0;
  b2limb_t *work, *next;

  b2_add_work_size(&worklimbs, asumcap);
  b2_add_work_size(&worklimbs, bsumcap);
  b2_add_work_size(&worklimbs, z0cap);
  b2_add_work_size(&worklimbs, z1cap);
  b2_add_work_size(&worklimbs, z2cap);
  work = (b2limb_t*)b2_xmalloc(worklimbs, sizeof(b2limb_t));
  next = work;

  b2_init_work(&asum, &next, asumcap);
  b2_init_work(&bsum, &next, bsumcap);
  b2_init_work(&z0, &next, z0cap);
  b2_init_work(&z1, &next, z1cap);
  b2_init_work(&z2, &next, z2cap);
  MPUassert(next == work + worklimbs,
            "b2 Karatsuba workspace mismatch");

  b2_init_slice_view(&a0, a, 0, split);
  b2_init_slice_view(&a1, a, split, an-split);
  b2_init_slice_view(&b0, b, 0, split);
  b2_init_slice_view(&b1, b, split, bn-split);

  b2_mul_abs(&z0, &a0, &b0);
  b2_mul_abs(&z2, &a1, &b1);
  b2_add(&asum, &a0, &a1);
  b2_add(&bsum, &b0, &b1);
  b2_mul_abs(&z1, &asum, &bsum);
  b2_sub_abs_inplace(&z1, &z0);
  b2_sub_abs_inplace(&z1, &z2);

  b2_ensure(out, rn);
  memset(out->d, 0, (size_t)rn * sizeof(b2limb_t));
  b2_add_shifted(out, rn, &z0, 0);
  b2_add_shifted(out, rn, &z1, split);
  b2_add_shifted(out, rn, &z2, split + split);
  out->n = rn;
  out->neg = 0;
  b2_trim(out);

  free(work);
}

static void b2_sqr_karatsuba(b2_t *out, const b2_t *a)
{
  b2_t a0, a1, z0, z1, z2;
  uint32_t an = a->n;
  uint32_t split = an >> 1;
  uint32_t rn = an + an;
  uint32_t a0cap = split;
  uint32_t a1cap = an-split;
  uint32_t z0cap = a0cap + a0cap;
  uint32_t z1cap = a0cap + a1cap;
  uint32_t z2cap = a1cap + a1cap;
  size_t worklimbs = 0;
  b2limb_t *work, *next;

  b2_add_work_size(&worklimbs, z0cap);
  b2_add_work_size(&worklimbs, z1cap);
  b2_add_work_size(&worklimbs, z2cap);
  work = (b2limb_t*)b2_xmalloc(worklimbs, sizeof(b2limb_t));
  next = work;

  b2_init_work(&z0, &next, z0cap);
  b2_init_work(&z1, &next, z1cap);
  b2_init_work(&z2, &next, z2cap);
  MPUassert(next == work + worklimbs,
            "b2 Karatsuba square workspace mismatch");

  b2_init_slice_view(&a0, a, 0, split);
  b2_init_slice_view(&a1, a, split, an-split);
  b2_sqr_abs(&z0, &a0);
  b2_mul_abs(&z1, &a0, &a1);
  b2_sqr_abs(&z2, &a1);

  b2_ensure(out, rn);
  memset(out->d, 0, (size_t)rn * sizeof(b2limb_t));
  b2_add_shifted(out, rn, &z0, 0);
  b2_add_shifted(out, rn, &z1, split);
  b2_add_shifted(out, rn, &z1, split);
  b2_add_shifted(out, rn, &z2, split + split);
  out->n = rn;
  out->neg = 0;
  b2_trim(out);

  free(work);
}

static void b2_mul_abs(b2_t *out, const b2_t *a, const b2_t *b)
{
  uint32_t minn = a->n < b->n ? a->n : b->n;
  uint32_t maxn = a->n > b->n ? a->n : b->n;

  /* Very unequal operands do not recover the extra additions and copies. */
  if (minn >= B2_KARATSUBA_THRESHOLD && maxn-minn <= minn/2)
    b2_mul_karatsuba(out, a, b);
  else
    b2_mul_schoolbook(out, a, b);
}

static void b2_sqr_abs(b2_t *out, const b2_t *a)
{
  if (a->n >= B2_SQR_KARATSUBA_THRESHOLD)
    b2_sqr_karatsuba(out, a);
  else
    b2_sqr_schoolbook(out, a);
}

void b2_mul(b2_t *out, const b2_t *a, const b2_t *b)
{
  const b2_t *large;
  b2limb_t small;
  b2_t tmp, *r = out;
  bool alias = out == a || out == b;
  bool neg;

  if (a->n == 0 || b->n == 0) {
    out->n = 0;
    out->neg = 0;
    return;
  }
  neg = a->neg != b->neg;
  if (a->n == 1 || b->n == 1) {
    large = a->n == 1 ? b : a;
    small = a->n == 1 ? a->d[0] : b->d[0];
    b2_mul_limb_abs(out, large, small);
    out->neg = neg && out->n != 0;
    return;
  }
  if (a == b) {
    b2_sqr(out, a);
    return;
  }
  if (a->n > UINT32_MAX-b->n)
    croak("internal: b2 multiplication too large");
  if (alias) {
    b2_init(&tmp);
    r = &tmp;
  }
  b2_mul_abs(r, a, b);
  r->neg = neg && r->n != 0;
  if (alias) b2_move(out, &tmp);
}

void b2_sqr(b2_t *out, const b2_t *a)
{
  b2_t tmp, *r = out;

  if (a->n == 0) {
    out->n = 0;
    out->neg = 0;
    return;
  }
  if (a->n == 1) {
    b2acc_t product = (b2acc_t)a->d[0] * a->d[0];
    b2_ensure(out, 2);
    out->d[0] = (b2limb_t)product;
    out->d[1] = (b2limb_t)(product >> B2_LIMB_BITS);
    out->n = out->d[1] != 0 ? 2 : 1;
    out->neg = 0;
    return;
  }
  if (a->n > UINT32_MAX/2)
    croak("internal: b2 square too large");
  if (out == a) {
    b2_init(&tmp);
    r = &tmp;
  }
  b2_sqr_abs(r, a);
  r->neg = 0;
  if (out == a) b2_move(out, &tmp);
}

/******************************************************************************/
/* Division by one native divisor. */
/******************************************************************************/

/* Return |a| modulo nonzero p. */
uint32_t b2_mod_u32(const b2_t *a, uint32_t p)
{
  uint32_t i, rem = 0;

  MPUassert(p != 0, "b2_mod_u32 divisor is zero");
  if (p == 1 || a->n == 0) return 0;
  if ((p & (p-1)) == 0) {
    rem = (uint32_t)a->d[0];
#if B2_LIMB_BITS < 32
    if (a->n > 1) rem |= (uint32_t)a->d[1] << B2_LIMB_BITS;
#endif
    return rem & (p-1);
  }
#if HAVE_UINT64
  for (i = a->n; i-- > 0; ) {
    b2acc_t cur = ((b2acc_t)rem << B2_LIMB_BITS) | a->d[i];
    rem = (uint32_t)(cur % p);
  }
#else
  if (p <= B2_LIMB_MAX) {
    for (i = a->n; i-- > 0; ) {
      b2acc_t cur = ((b2acc_t)rem << B2_LIMB_BITS) | a->d[i];
      rem = (uint32_t)(cur % p);
    }
  } else {
    for (i = a->n; i-- > 0; ) {
      b2limb_t limb = a->d[i];
      uint32_t bit;
      for (bit = B2_LIMB_BITS; bit-- > 0; ) {
        if (rem >= p-rem) rem -= p-rem;
        else              rem += rem;
        if ((limb >> bit) & 1) {
          if (++rem == p) rem = 0;
        }
      }
    }
  }
#endif
  return rem;
}

/* Divide a by nonzero p in place, preserving its sign. */
void b2_divexact_u32(b2_t *a, uint32_t p)
{
  uint32_t i, rem = 0;

  MPUassert(p != 0, "b2_divexact_u32 divisor is zero");
  if (p == 1 || a->n == 0) return;
  if ((p & (p-1)) == 0) {
    MPUassert(b2_mod_u32(a, p) == 0,
              "b2_divexact_u32 division is not exact");
    b2_tdiv_2exp(a, a, (UV)ctz((UV)p));
    return;
  }
#if HAVE_UINT64
  for (i = a->n; i-- > 0; ) {
    b2acc_t cur = ((b2acc_t)rem << B2_LIMB_BITS) | a->d[i];
    a->d[i] = (b2limb_t)(cur / p);
    rem = (uint32_t)(cur % p);
  }
#else
  if (p <= B2_LIMB_MAX) {
    for (i = a->n; i-- > 0; ) {
      b2acc_t cur = ((b2acc_t)rem << B2_LIMB_BITS) | a->d[i];
      a->d[i] = (b2limb_t)(cur / p);
      rem = (uint32_t)(cur % p);
    }
  } else {
    for (i = a->n; i-- > 0; ) {
      b2limb_t input = a->d[i], quotient = 0;
      uint32_t bit;
      for (bit = B2_LIMB_BITS; bit-- > 0; ) {
        bool qbit;
        if (rem >= p-rem) {
          rem -= p-rem;
          qbit = 1;
        } else {
          rem += rem;
          qbit = 0;
        }
        if ((input >> bit) & 1) {
          if (++rem == p) {
            rem = 0;
            qbit = 1;
          }
        }
        if (qbit) quotient |= (b2limb_t)1 << bit;
      }
      a->d[i] = quotient;
    }
  }
#endif
  MPUassert(rem == 0, "b2_divexact_u32 division is not exact");
  b2_trim(a);
}

void b2_pow(b2_t *out, const b2_t *a, UV exp)
{
  b2_t sq, tmp;
  b2_init_set(&sq, a);
  b2_set_uv(out, 1);
  b2_init(&tmp);
  while (exp > 0) {
    if (exp & 1) {
      b2_mul(&tmp, out, &sq);
      b2_move(out, &tmp);
    }
    exp >>= 1;
    if (exp > 0) {
      b2_sqr(&tmp, &sq);
      b2_move(&sq, &tmp);
    }
  }
  b2_free(&sq);
  b2_free(&tmp);
}

void b2_product(b2_t A[], size_t a, size_t b)
{
  if (b <= a) {
    /* A[a] already contains the result. */
  } else if (b == a+1) {
    b2_mul(&A[a], &A[a], &A[b]);
  } else if (b == a+2) {
    b2_mul(&A[a+1], &A[a+1], &A[a+2]);
    b2_mul(&A[a], &A[a], &A[a+1]);
  } else {
    size_t c = a + (b-a+1)/2;
    b2_product(A, a, c-1);
    b2_product(A, c, b);
    b2_mul(&A[a], &A[a], &A[c]);
  }
}

void b2_product_u32(b2_t *out, const uint32_t A[], size_t len)
{
  const unsigned per_word = BITS_PER_WORD / 32;
  b2_t *B;
  UV prod = 1;
  size_t i, max_chunks, nprod = 0;
  unsigned used = 0;

  if (len == 0) {
    b2_set_uv(out, 1);
    return;
  }
  max_chunks = len / per_word + (len % per_word != 0);
  B = (b2_t*)b2_xmalloc(max_chunks, sizeof(b2_t));
  for (i = 0; i < len; i++) {
    if (A[i] == 0) { prod = 0; break; }
    if (A[i] == 1) continue;
    if (used == per_word) {
      b2_init_set_uv(&B[nprod++], prod);
      prod = A[i];
      used = 1;
    } else {
      prod *= (UV)A[i];
      used++;
    }
  }
  if (prod == 0) {
    b2_set_uv(out, 0);
  } else {
    if (used != 0) b2_init_set_uv(&B[nprod++], prod);
    if (nprod == 0) {
      b2_set_uv(out, 1);
    } else {
      b2_product(B, 0, nprod-1);
      b2_move(out, &B[0]);
    }
  }
  for (i = 0; i < nprod; i++) b2_free(&B[i]);
  free(B);
}

size_t b2_bit_length(const b2_t *x)
{
  size_t bits;
  if (x->n == 0) return 0;
  if ((size_t)x->n > (size_t)MAX_SIZET / B2_LIMB_BITS)
    croak("internal: b2 bit length overflow");
  bits = (size_t)(x->n - 1) * B2_LIMB_BITS;
  return bits + (size_t)log2floor((UV)x->d[x->n-1]) + 1;
}

void b2_mul_2exp(b2_t *out, const b2_t *a, UV bits)
{
  UV word_uv;
  uint32_t oldn, word, newn, i, extra;
  uint32_t rem;
  b2acc_t carry = 0;

  if (out != a) b2_set(out, a);
  if (bits == 0 || out->n == 0) return;
  word_uv = bits / B2_LIMB_BITS;
  rem = (uint32_t)(bits % B2_LIMB_BITS);
  oldn = out->n;
  extra = rem != 0;
  if (oldn > UINT32_MAX-extra ||
      word_uv > (UV)(UINT32_MAX-oldn-extra))
    croak("internal: b2 left shift too large");
  word = (uint32_t)word_uv;
  b2_ensure(out, oldn + word + extra);
  if (word != 0) {
    memmove(out->d + word, out->d, (size_t)oldn * sizeof(b2limb_t));
    memset(out->d, 0, (size_t)word * sizeof(b2limb_t));
  }
  newn = oldn + word;
  if (rem != 0) {
    for (i = word; i < newn; i++) {
      b2acc_t cur = ((b2acc_t)out->d[i] << rem) | carry;
      out->d[i] = (b2limb_t)cur;
      carry = cur >> B2_LIMB_BITS;
    }
    if (carry != 0) out->d[newn++] = (b2limb_t)carry;
  }
  out->n = newn;
}

/* Divide in place by 2^bits and report whether any discarded bit was set. */
static bool b2_tdiv_2exp_inplace(b2_t *x, UV bits)
{
  UV word_uv;
  uint32_t word, rem, i, newn;
  bool discarded = 0;
  b2limb_t carry = 0, mask;

  if (bits == 0 || x->n == 0) return 0;
  word_uv = bits / B2_LIMB_BITS;
  if (word_uv >= x->n) {
    x->n = 0;
    x->neg = 0;
    return 1;
  }
  word = (uint32_t)word_uv;
  rem = (uint32_t)(bits % B2_LIMB_BITS);
  for (i = 0; i < word && !discarded; i++)
    discarded = x->d[i] != 0;
  if (rem != 0) {
    mask = (b2limb_t)(((b2limb_t)1 << rem) - 1);
    if ((x->d[word] & mask) != 0) discarded = 1;
  }

  newn = x->n - word;
  if (word != 0)
    memmove(x->d, x->d + word, (size_t)newn * sizeof(b2limb_t));
  x->n = newn;
  if (rem == 0) return discarded;

  for (i = x->n; i-- > 0; ) {
    b2limb_t limb = x->d[i];
    x->d[i] = (b2limb_t)((limb >> rem) |
                         (carry << (B2_LIMB_BITS-rem)));
    carry = limb & mask;
  }
  b2_trim(x);
  return discarded;
}

void b2_tdiv_2exp(b2_t *out, const b2_t *a, UV bits)
{
  if (out != a) b2_set(out, a);
  (void)b2_tdiv_2exp_inplace(out, bits);
}

static void b2_add_one(b2_t *x)
{
  uint32_t i = 0;
  if (x->n == UINT32_MAX)
    croak("internal: b2 addition too large");
  b2_ensure(x, x->n + 1);
  while (i < x->n && x->d[i] == B2_LIMB_MAX) x->d[i++] = 0;
  if (i < x->n)
    x->d[i]++;
  else
    x->d[x->n++] = 1;
}

void b2_fdiv_2exp(b2_t *out, const b2_t *a, UV bits)
{
  bool neg = b2_is_negative(a);
  bool adjust;
  if (out != a) b2_set(out, a);
  adjust = b2_tdiv_2exp_inplace(out, bits);
  if (neg && adjust) {
    b2_abs(out);
    b2_add_one(out);
    b2_set_negative(out, 1);
  }
}

/******************************************************************************/
/* Magnitude-oriented bit operations. */
/******************************************************************************/

UV b2_popcount(const b2_t *x)
{
  UV count = 0;
  uint32_t i;
  for (i = 0; i < x->n; i++) {
    UV add = popcnt((UV)x->d[i]);
    if (count > UV_MAX-add)
      croak("internal: b2 popcount overflow");
    count += add;
  }
  return count;
}

bool b2_testbit(const b2_t *x, uint32_t k)
{
  uint32_t word = k / B2_LIMB_BITS;
  uint32_t bit = k % B2_LIMB_BITS;
  return word < x->n && ((x->d[word] >> bit) & 1) != 0;
}

void b2_setbit(b2_t *x, uint32_t k)
{
  uint32_t word = k / B2_LIMB_BITS;
  uint32_t bit = k % B2_LIMB_BITS;
  uint32_t oldn = x->n;
  b2_ensure(x, word + 1);
  if (oldn <= word) {
    memset(x->d + oldn, 0, (size_t)(word + 1-oldn) * sizeof(b2limb_t));
    x->n = word + 1;
  }
  x->d[word] |= (b2limb_t)((b2limb_t)1 << bit);
  x->neg = 0;
}

void b2_clearbit(b2_t *x, uint32_t k)
{
  uint32_t word = k / B2_LIMB_BITS;
  uint32_t bit = k % B2_LIMB_BITS;
  x->neg = 0;
  if (word >= x->n) return;
  x->d[word] &= (b2limb_t)~((b2limb_t)1 << bit);
  b2_trim(x);
}

void b2_flipbit(b2_t *x, uint32_t k)
{
  uint32_t word = k / B2_LIMB_BITS;
  uint32_t bit = k % B2_LIMB_BITS;
  uint32_t oldn = x->n;
  b2_ensure(x, word + 1);
  if (oldn <= word) {
    memset(x->d + oldn, 0, (size_t)(word + 1-oldn) * sizeof(b2limb_t));
    x->n = word + 1;
  }
  x->d[word] ^= (b2limb_t)((b2limb_t)1 << bit);
  x->neg = 0;
  b2_trim(x);
}

static bool b2_make_index(uint32_t *result, uint32_t word, uint32_t bit)
{
  uint32_t base;
  if (word > UINT32_MAX / B2_LIMB_BITS) return 0;
  base = word * B2_LIMB_BITS;
  if (bit > UINT32_MAX-base) return 0;
  *result = base + bit;
  return 1;
}

bool b2_bitscan1(uint32_t *result, const b2_t *x, uint32_t start)
{
  uint32_t word = start / B2_LIMB_BITS;
  uint32_t bit = start % B2_LIMB_BITS;
  b2limb_t v, mask;

  if (word >= x->n) return 0;
  mask = bit == 0 ? B2_LIMB_MAX :
         (b2limb_t)(B2_LIMB_MAX << bit);
  v = x->d[word] & mask;
  while (1) {
    if (v != 0)
      return b2_make_index(result, word, (uint32_t)ctz((UV)v));
    if (++word >= x->n) return 0;
    v = x->d[word];
  }
}

bool b2_bitscan0(uint32_t *result, const b2_t *x, uint32_t start)
{
  uint32_t word = start / B2_LIMB_BITS;
  uint32_t bit = start % B2_LIMB_BITS;
  b2limb_t v, mask;

  if (word >= x->n) { *result = start; return 1; }
  mask = bit == 0 ? B2_LIMB_MAX :
         (b2limb_t)(B2_LIMB_MAX << bit);
  v = (b2limb_t)~x->d[word] & mask;
  while (1) {
    if (v != 0)
      return b2_make_index(result, word, (uint32_t)ctz((UV)v));
    if (++word >= x->n) return b2_make_index(result, word, 0);
    v = (b2limb_t)~x->d[word];
  }
}

enum {
  B2_BITOP_AND,
  B2_BITOP_OR,
  B2_BITOP_XOR,
  B2_BITOP_ANDNOT
};

static void b2_bitop(b2_t *out, const b2_t *a, const b2_t *b, int op)
{
  uint32_t an = a->n, bn = b->n, rn, i;
  if (op == B2_BITOP_AND)
    rn = an < bn ? an : bn;
  else if (op == B2_BITOP_ANDNOT)
    rn = an;
  else
    rn = an > bn ? an : bn;
  b2_ensure(out, rn);
  for (i = 0; i < rn; i++) {
    b2limb_t av = i < an ? a->d[i] : 0;
    b2limb_t bv = i < bn ? b->d[i] : 0;
    switch (op) {
      case B2_BITOP_AND:    out->d[i] = av & bv;  break;
      case B2_BITOP_OR:     out->d[i] = av | bv;  break;
      case B2_BITOP_XOR:    out->d[i] = av ^ bv;  break;
      default:              out->d[i] = av & (b2limb_t)~bv;  break;
    }
  }
  out->n = rn;
  out->neg = 0;
  b2_trim(out);
}

void b2_bitand(b2_t *out, const b2_t *a, const b2_t *b)
{
  b2_bitop(out, a, b, B2_BITOP_AND);
}

void b2_bitor(b2_t *out, const b2_t *a, const b2_t *b)
{
  b2_bitop(out, a, b, B2_BITOP_OR);
}

void b2_bitxor(b2_t *out, const b2_t *a, const b2_t *b)
{
  b2_bitop(out, a, b, B2_BITOP_XOR);
}

void b2_bitandnot(b2_t *out, const b2_t *a, const b2_t *b)
{
  b2_bitop(out, a, b, B2_BITOP_ANDNOT);
}

void b2_bitnot(b2_t *out, const b2_t *a, bool fixed_width, uint32_t width)
{
  size_t nbits, rn_size;
  uint32_t an = a->n, rn, i, topbits;

  if (!fixed_width && a->n == 0) {
    b2_set_uv(out, 1);
    return;
  }
  if (fixed_width && width == 0) {
    b2_set_uv(out, 0);
    return;
  }
  nbits = fixed_width ? (size_t)width : b2_bit_length(a);
  rn_size = nbits / B2_LIMB_BITS + (nbits % B2_LIMB_BITS != 0);
  if (rn_size > UINT32_MAX)
    croak("internal: b2 complement width too large");
  rn = (uint32_t)rn_size;
  b2_ensure(out, rn);
  for (i = 0; i < rn; i++) {
    b2limb_t av = i < an ? a->d[i] : 0;
    out->d[i] = (b2limb_t)~av;
  }
  topbits = (uint32_t)(nbits % B2_LIMB_BITS);
  if (topbits != 0)
    out->d[rn-1] &= (b2limb_t)(((b2limb_t)1 << topbits) - 1);
  out->n = rn;
  out->neg = 0;
  b2_trim(out);
}
