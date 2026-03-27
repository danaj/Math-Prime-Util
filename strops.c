#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "strops.h"

/******************************************************************************/
/*               BASE-B9 INTERNAL REPRESENTATION FOR FAST ARITHMETIC         */
/*                                                                            */
/* Numbers are stored internally as little-endian arrays of b9limb_t, each   */
/* holding B9_DIGS decimal digits (one "digit" in base B9_BASE).  The        */
/* accumulator type b9acc_t is wide enough that the schoolbook inner-product  */
/* loop never overflows for any practical input length.                       */
/*                                                                            */
/* Tier selection is based on compile-time capability flags from ptypes.h:   */
/*   HAVE_UINT128: base 10^9, 128-bit accumulator (unlimited safe length)    */
/*   HAVE_UINT64:  base 10^6,  64-bit accumulator (per-row flush, any length) */
/*   else:         base 10^4,  32-bit accumulator (uses per-row carry flush) */
/******************************************************************************/

#if HAVE_UINT128
  #define B9_DIGS  9
  #define B9_BASE  UINT32_C(1000000000)
  typedef uint32_t  b9limb_t;
  typedef uint128_t b9acc_t;
#elif HAVE_UINT64
  #define B9_DIGS  6
  #define B9_BASE  UINT32_C(1000000)
  typedef uint32_t  b9limb_t;
  typedef uint64_t  b9acc_t;
#else
  #define B9_DIGS  4
  #define B9_BASE  UINT32_C(10000)
  typedef uint32_t  b9limb_t;
  typedef uint32_t  b9acc_t;
#endif

/* Limbs needed for a decimal string of len digits. */
#define B9_NLIMBS(len)  (((len) + B9_DIGS - 1) / B9_DIGS)

/******************************************************************************/
/*                       B9 ARITHMETIC  (GMP-STYLE INTERFACE)                */
/*                                                                            */
/* b9_t holds a signed big integer as a managed little-endian limb array.    */
/* All b9_* functions are static (file-scope only).                          */
/******************************************************************************/

typedef struct {
    b9limb_t *d;      /* limb array, little-endian */
    uint32_t  alloc;  /* allocated limb slots */
    uint32_t  n;      /* significant limbs (0 means value is zero) */
    int       neg;    /* 1 if negative, 0 otherwise */
} b9_t;

static void b9_init(b9_t *x)
  { x->d = NULL;  x->alloc = 0;  x->n = 0;  x->neg = 0; }

static void b9_free(b9_t *x)
  { free(x->d);  x->d = NULL;  x->alloc = 0;  x->n = 0;  x->neg = 0; }

static void b9_ensure(b9_t *x, uint32_t need)
{
  if (x->alloc >= need) return;
  x->d     = (b9limb_t *) realloc(x->d, (size_t)need * sizeof(b9limb_t));
  x->alloc = need;
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
  b9_ensure(x, (20 + B9_DIGS - 1) / B9_DIGS);
  x->neg = 0;
  if (v == 0) { x->n = 0;  return; }
  while (v > 0) { x->d[n++] = (b9limb_t)((UV)(v % B9_BASE)); v /= (UV)B9_BASE; }
  x->n = n;
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

/* out = a - b (signed).  out may alias a or b. */
static void b9_sub(b9_t *out, const b9_t *a, const b9_t *b)
{
  b9_t nb = *b;   /* shallow copy: nb.d aliases b->d, only neg changes */
  b9_neg(&nb);
  b9_add(out, a, &nb);
}

static void b9_add_uv(b9_t *out, const b9_t *a, UV v)
{
  b9_t b;
  b9_init(&b);
  b9_set_uv(&b, v);
  b9_add(out, a, &b);
  b9_free(&b);
}

/* out = a * b (signed).
 * Reads all of a->d and b->d before writing to out->d, so out may alias
 * a or b (b9_ensure is called only after the multiply loop). */
static void b9_mul(b9_t *out, const b9_t *a, const b9_t *b)
{
  uint32_t i, j, rn;
  b9acc_t *acc;
  int neg;

  if (a->n == 0 || b->n == 0) { out->n = 0;  out->neg = 0;  return; }

  neg = (a->neg != b->neg) ? 1 : 0;
  rn  = a->n + b->n;
  acc = (b9acc_t*) calloc((size_t)(rn + 1), sizeof(b9acc_t));

  for (i = 0; i < a->n; i++) {
    if (a->d[i] == 0) continue;
    for (j = 0; j < b->n; j++)
      acc[i + j] += (b9acc_t)a->d[i] * b->d[j];
#if !HAVE_UINT128
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
  free(acc);

  while (rn > 1 && out->d[rn-1] == 0) rn--;
  out->n   = rn;
  out->neg = (out->n == 1 && out->d[0] == 0) ? 0 : neg;
}

/* out = a^exp.  out must not alias a; a is not modified. */
static void b9_pow(b9_t *out, const b9_t *a, UV exp)
{
    b9_t sq, tmp;
    b9_init(&sq);
    /* sq = copy of a */
    b9_ensure(&sq, a->n);
    if (a->n) memcpy(sq.d, a->d, a->n * sizeof(b9limb_t));
    sq.n = a->n;  sq.neg = a->neg;
    b9_set_uv(out, 1);
    b9_init(&tmp);
    while (exp > 0) {
        if (exp & 1) {
            b9_mul(&tmp, out, &sq);
            b9_free(out);  *out = tmp;  b9_init(&tmp);
        }
        exp >>= 1;
        if (exp > 0) {
            b9_mul(&tmp, &sq, &sq);
            b9_free(&sq);  sq = tmp;  b9_init(&tmp);
        }
    }
    b9_free(&sq);  b9_free(&tmp);
}

/* Signed floor division and remainder (floor convention: rem has sign of b).
 * q or r may be NULL (that output is skipped).
 * b must be non-zero.  q and r must not alias a or b. */
static void b9_divmod(b9_t *q, b9_t *r, const b9_t *a, const b9_t *b)
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


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


/* Parse a signed decimal integer string: strip optional sign and leading
 * zeros, validate that all remaining characters are digits, and update *sp
 * and *slen to point to the canonical digit sequence.  Returns 1 if negative,
 * 0 otherwise.  Croaks on invalid input. */
static int strint_parse(const char **sp, STRLEN *slen)
{
  const char* s = *sp;
  STRLEN i = 0, len = *slen;
  int neg = 0;

  if (s != 0 && len > 0) {
    neg = (s[0] == '-');
    if (s[0] == '-' || s[0] == '+') { s++; len--; }
    while (len > 0 && *s == '0') { s++; len--; }
    if (len == 0) { s--; len = 1; neg = 0; }  /* value is 0 */
    for (i = 0; i < len; i++)
      if (!isDIGIT(s[i]))
        break;
  }
  if (s == 0 || len == 0 || i < len) croak("Parameter must be an integer");
  *sp = s;
  *slen = len;
  return neg;
}

int strint_cmp(const char* a, STRLEN alen, const char* b, STRLEN blen) {
  STRLEN i;
  int aneg = strint_parse(&a, &alen);
  int bneg = strint_parse(&b, &blen);
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
  bneg = strint_parse(&b, &blen);

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
#if 1
  if (blen == 1 && b[0] == '1') {
    int aneg = a[0] == '-', negb = !!negate_b;
    if ((aneg==negb && a[alen-1] != '9') || (aneg!=negb && a[alen-1] != '0')) {
      memcpy(out, a, alen-1);
      out[alen-1] = a[alen-1] + (aneg==negb ? 1 : -1);
      return alen;
    }
  }
#endif
  b9_init(&A);  b9_init(&B);
  b9_set_str(&A, a, alen);
  b9_set_str(&B, b, blen);
  if (negate_b) b9_neg(&B);
  b9_add(&A, &A, &B);
  rlen = b9_get_str(out, &A);
  b9_free(&A);  b9_free(&B);
  return rlen;
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

  b9_init(&A);  b9_init(&B);  b9_init(&R);
  b9_set_str(&A, a, alen);
  b9_set_str(&B, b, blen);
  b9_mul(&R, &A, &B);

  rlen = b9_get_str(out, &R);
  b9_free(&A);  b9_free(&B);  b9_free(&R);

  return rlen;
}

STRLEN strint_muladd_s(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen, const char* c, STRLEN clen, bool negate_c)
{
  b9_t A, B, C;
  STRLEN rlen;

  b9_init(&A);  b9_init(&B);  b9_init(&C);
  b9_set_str(&A, a, alen);
  b9_set_str(&B, b, blen);
  b9_set_str(&C, c, clen);
  if (negate_c) b9_neg(&C);

  b9_mul(&A, &A, &B);
  b9_add(&A, &A, &C);

  rlen = b9_get_str(out, &A);
  b9_free(&A);  b9_free(&B);  b9_free(&C);

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

  b9_init(&base);  b9_init(&result);
  b9_set_str(&base, a, alen);
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
 * qout needs alen+2 bytes; rout needs blen+1 bytes.
 * Remainder has the sign of b (floor convention).
 * Returns false if b = 0, true otherwise. */
bool strint_divmod(char* qout, STRLEN* qlen, char* rout, STRLEN* rlen,
                   const char* a, STRLEN alen,
                   const char* b, STRLEN blen)
{
  b9_t ba, bb, bq, br;
  b9_init(&ba);  b9_init(&bb);  b9_init(&bq);  b9_init(&br);
  b9_set_str(&ba, a, alen);
  b9_set_str(&bb, b, blen);
  if (bb.n == 0) { b9_free(&ba);  b9_free(&bb);  return 0; }
  b9_divmod(&bq, &br, &ba, &bb);
  if (qout) { STRLEN l = b9_get_str(qout, &bq);  if (qlen) *qlen = l; }
  if (rout) { STRLEN l = b9_get_str(rout, &br);  if (rlen) *rlen = l; }
  b9_free(&ba);  b9_free(&bb);  b9_free(&bq);  b9_free(&br);
  return 1;
}

STRLEN strint_divint(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen)
{
  STRLEN len;
  if (!strint_divmod(out, &len, 0, 0, a,alen, b,blen)) return 0;
  return len;
}

STRLEN strint_modint(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen)
{
  STRLEN len;
  if (!strint_divmod(0, 0, out, &len, a,alen, b,blen)) return 0;
  return len;
}

STRLEN strint_cdivint(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen)
{
  STRLEN len;
  b9_t ba, bb, bq, br;
  b9_init(&ba);  b9_init(&bb);  b9_init(&bq);  b9_init(&br);
  b9_set_str(&ba, a, alen);
  b9_set_str(&bb, b, blen);
  if (bb.n == 0) { b9_free(&ba);  b9_free(&bb);  return 0; }
  b9_divmod(&bq, &br, &ba, &bb);
  if (br.n > 0)
    b9_add_uv(&bq, &bq, 1);
  len = b9_get_str(out, &bq);
  b9_free(&ba);  b9_free(&bb);  b9_free(&bq);  b9_free(&br);
  return len;
}


/******************************************************************************/
/*                           INTEGER LOG AND ROOT                             */
/******************************************************************************/

/* Write UV v as a decimal string; return length.  buf must have >= 21 bytes. */
static STRLEN uv_to_str(char* buf, UV v)
{
  STRLEN i = 0, j;
  if (v == 0) { buf[0] = '0'; return 1; }
  while (v > 0) { buf[i++] = '0' + (char)(v % 10); v /= 10; }
  for (j = 0; j < i / 2; j++) {
    char c = buf[j]; buf[j] = buf[i-1-j]; buf[i-1-j] = c;
  }
  return i;
}

/* Returns floor(log_base(a)), or UV_MAX on error (a < 1 or base < 2).
 * Uses a floating-point estimate to land within 1-2 of the answer,
 * then verifies/adjusts with b9 arithmetic. */
UV strint_logint(const char* a, STRLEN alen, UV base)
{
    STRLEN nd, i;
    UV k;
    double approx, log10_a, k_est_f;
    int aneg;
    b9_t bn, bbase, pow, tmp;

    if (base < 2) return UV_MAX;

    aneg = strint_parse(&a, &alen);
    if (aneg || (alen == 1 && a[0] == '0')) return UV_MAX;

    /* Floating-point estimate: log10(a) via leading digits + digit count */
    nd = (alen < 15) ? alen : 15;
    approx = 0.0;
    for (i = 0; i < nd; i++)
        approx = approx * 10.0 + (double)(a[i] - '0');
    log10_a = (double)(alen - nd) + log10(approx);
    k_est_f = log10_a / log10((double)base);
    k = (k_est_f <= 1.0) ? 0 : (UV)k_est_f - 1;  /* start one below estimate */

    b9_init(&bn);    b9_set_str(&bn, a, alen);
    b9_init(&bbase); b9_set_uv(&bbase, base);
    b9_init(&pow);   b9_pow(&pow, &bbase, k);
    b9_init(&tmp);

    /* Float overshot (very rare): step down by dividing until base^k ≤ a */
    while (k > 0 && b9_cmp(&pow, &bn) > 0) {
        b9_divmod(&tmp, NULL, &pow, &bbase);
        { b9_t s = pow;  pow = tmp;  tmp = s; }
        k--;
    }

    /* Step up while base^(k+1) ≤ a */
    for (;;) {
        b9_mul(&tmp, &pow, &bbase);
        if (b9_cmp(&tmp, &bn) > 0) break;
        { b9_t s = pow;  pow = tmp;  tmp = s; }
        k++;
    }

    b9_free(&bn);  b9_free(&bbase);  b9_free(&pow);  b9_free(&tmp);
    return k;
}

/* floor(n^(1/k)) via Newton's method, starting from above.
 * out must have at least alen+2 bytes.
 * Returns the length of the result, or 0 on error (k==0, n<0). */
STRLEN strint_rootint(char* out, const char* a, STRLEN alen, UV k)
{
    STRLEN nd, i, rlen;
    double approx, log10_a, r_est_f;
    int aneg;
    b9_t bn, bkm1, bk, r, pk1, q, tmp;

    if (k == 0) return 0;

    aneg = strint_parse(&a, &alen);
    if (aneg) return 0;

    if (k == 1) { memcpy(out, a, alen); return alen; }

    if (alen == 1 && (a[0] == '0' || a[0] == '1')) { out[0] = a[0]; return 1; }

    /* Float estimate of floor(log10(n^(1/k))) */
    nd = (alen < 15) ? alen : 15;
    approx = 0.0;
    for (i = 0; i < nd; i++)
        approx = approx * 10.0 + (double)(a[i] - '0');
    log10_a = (double)(alen - nd) + log10(approx);
    r_est_f = log10_a / (double)k;

    b9_init(&bn);   b9_set_str(&bn, a, alen);
    b9_init(&bkm1); b9_set_uv(&bkm1, k - 1);
    b9_init(&bk);   b9_set_uv(&bk,   k);
    b9_init(&r);
    b9_init(&pk1);  b9_init(&q);  b9_init(&tmp);

    /* Better initial r: float within ~0.02% of true root, biased above.
     * r_est_f = q_exp + frac.  r_mantissa ≈ 10^frac * 1.0002e8 (9 sig figs).
     * For q_exp >= mlen-1: r = r_mantissa * 10^(q_exp+1-mlen) (char temp, one malloc).
     * For q_exp < mlen-1:  r_est_f < 9, direct UV. */
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
            b9_set_str(&r, r_char, q_exp + 1);
            free(r_char);
        } else {
            b9_set_uv(&r, (UV)(pow(10.0, r_est_f) * 1.0002) + 1);
        }
    }

    /* Newton: r_new = floor(((k-1)*r + floor(n/r^(k-1))) / k).
     * Converges from above; terminate when r_new >= r. */
    for (;;) {
        /* pk1 = r^(k-1) */
        b9_pow(&pk1, &r, k - 1);

        /* q = floor(n / r^(k-1)) */
        b9_divmod(&q, NULL, &bn, &pk1);

        /* tmp = (k-1)*r + q */
        b9_mul(&tmp, &bkm1, &r);
        b9_add(&tmp, &tmp, &q);

        /* pk1 = floor(tmp / k)  (reuse as r_new) */
        b9_divmod(&pk1, NULL, &tmp, &bk);

        if (b9_cmp(&pk1, &r) >= 0) break;

        /* r = r_new */
        { b9_t s = r;  r = pk1;  pk1 = s; }
    }

    rlen = b9_get_str(out, &r);
    b9_free(&bn);  b9_free(&bkm1);  b9_free(&bk);
    b9_free(&r);   b9_free(&pk1);   b9_free(&q);  b9_free(&tmp);
    return rlen;
}
