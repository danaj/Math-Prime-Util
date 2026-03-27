#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "strops.h"

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

/******************************************************************************/

/* Increment the decimal digits in buf[0..len-1] in place.
 * Returns len+1 on carry-out (all 9s), len otherwise. */
static STRLEN mag_incr(char* buf, STRLEN len)
{
  STRLEN i = len - 1;
  while (buf[i] == '9') {
    buf[i] = '0';
    if (i == 0) {
      memmove(buf + 1, buf, len);
      buf[0] = '1';
      return len + 1;
    }
    i--;
  }
  buf[i]++;
  return len;
}

/* Decrement the decimal digits in buf[0..len-1] in place.
 * Caller must ensure value > 1 so the result stays positive.
 * Returns len-1 if a leading zero is shed, len otherwise. */
static STRLEN mag_decr(char* buf, STRLEN len)
{
  STRLEN i = len - 1;
  while (buf[i] == '0') { buf[i] = '9'; i--; }
  buf[i]--;
  if (len > 1 && buf[0] == '0') {
    memmove(buf, buf + 1, len - 1);
    return len - 1;
  }
  return len;
}

/* Compare two unsigned decimal digit strings (no signs).
 * Returns -1, 0, or 1. */
static int mag_cmp(const char* a, STRLEN alen, const char* b, STRLEN blen)
{
  STRLEN i;
  if (alen != blen) return (alen > blen) ? 1 : -1;
  for (i = 0; i < alen; i++)
    if (a[i] != b[i]) return (a[i] > b[i]) ? 1 : -1;
  return 0;
}

/* Add two unsigned decimal digit strings.
 * out needs max(alen,blen)+1 bytes.  Returns result length. */
static STRLEN mag_add(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen)
{
  STRLEN rmax, k;
  int carry = 0;
  if (alen < blen) { /* ensure a is the longer string */
    const char* t = a;  STRLEN tl = alen;
    a = b;  alen = blen;  b = t;  blen = tl;
  }
  rmax = alen + 1;
  k = rmax;
  while (alen > 0 || carry) {
    int d = carry + (alen > 0 ? a[--alen] - '0' : 0)
                  + (blen > 0 ? b[--blen] - '0' : 0);
    carry = (d >= 10);
    out[--k] = '0' + (carry ? d - 10 : d);
  }
  if (k > 0) memmove(out, out + k, rmax - k);
  return rmax - k;
}

/* Subtract b from a (a >= b, caller ensures this).
 * out needs alen bytes.  Returns result length (leading zeros stripped). */
static STRLEN mag_sub(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen)
{
  STRLEN k = alen, ai = alen, bi = blen, start;
  int borrow = 0;
  while (ai > 0) {
    int d = (a[--ai] - '0') - borrow - (bi > 0 ? b[--bi] - '0' : 0);
    if (d < 0) { borrow = 1; d += 10; } else borrow = 0;
    out[--k] = '0' + d;
  }
  /* strip leading zeros, keeping at least one digit */
  for (start = 0; start < alen - 1 && out[start] == '0'; start++) ;
  if (start > 0) memmove(out, out + start, alen - start);
  return alen - start;
}

/* Core signed addition with signs already parsed and stripped. */
static STRLEN str_add_impl(char* out,
                            const char* a, STRLEN alen, int aneg,
                            const char* b, STRLEN blen, int bneg)
{
  STRLEN off, rlen;
  int cmp, neg_result;
  if (aneg == bneg) {
    off = aneg ? 1 : 0;
    rlen = mag_add(out + off, a, alen, b, blen);
    if (aneg) out[0] = '-';
    return rlen + off;
  }
  cmp = mag_cmp(a, alen, b, blen);
  if (cmp == 0) { out[0] = '0'; return 1; }
  neg_result = (cmp > 0) ? aneg : bneg;
  off = neg_result ? 1 : 0;
  rlen = (cmp > 0) ? mag_sub(out + off, a, alen, b, blen)
                   : mag_sub(out + off, b, blen, a, alen);
  if (neg_result) out[0] = '-';
  return rlen + off;
}

/* Add or subtract two signed decimal integer strings.
 * If negate_b is false, computes a+b; if true, computes a-b (i.e. a+(-b)).
 * Write result to out; caller must ensure at least max(alen,blen)+2 bytes
 * (the +2 is tight for subtraction when a is negative and blen >= alen).
 * out may alias a (in-place update is safe); out must not alias b.
 * Returns the result length (no NULL terminator written). */
STRLEN strint_add_s(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen, bool negate_b)
{
  int aneg, bneg;
  if (blen == 0 || (blen == 1 && b[0] == '0')) { memcpy(out, a, alen); return alen; }
  if (alen == 0 || (alen == 1 && a[0] == '0')) {
    if (!negate_b) { memcpy(out, b, blen); return blen; }
    /* 0 - b = -b */
    if (b[0] == '-') { memcpy(out, b + 1, blen - 1); return blen - 1; }
    out[0] = '-'; memcpy(out + 1, b, blen); return blen + 1;
  }
  aneg = (a[0] == '-'); if (aneg) { a++; alen--; }
  bneg = (b[0] == '-'); if (bneg) { b++; blen--; }
  return str_add_impl(out, a, alen, aneg, b, blen, bneg ^ negate_b);
}

/******************************************************************************/
/*                              MULTIPLICATION                                */
/******************************************************************************/

/* Multiply two unsigned decimal digit strings.
 * out needs alen+blen bytes.  Returns result length (no leading zeros). */
static STRLEN mag_mul(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen)
{
  STRLEN i, j, rlen = alen + blen, start;
  uint32_t carry;
  /* Accumulate products into LSB-first uint32_t array, deferring carry.
   * This moves the division-by-10 from the O(m*n) inner loop to a
   * single O(m+n) pass, and simplifies the index arithmetic. */
  uint32_t *acc = (uint32_t*)calloc(rlen, sizeof(uint32_t));
  for (i = 0; i < blen; i++) {
    uint32_t bd = b[blen - 1 - i] - '0';
    if (bd == 0) continue;
    for (j = 0; j < alen; j++)
      acc[i + j] += bd * (uint32_t)(a[alen - 1 - j] - '0');
  }
  carry = 0;
  for (i = 0; i < rlen; i++) {
    uint32_t d = acc[i] + carry;
    carry = d / 10;
    out[rlen - 1 - i] = '0' + d % 10;
  }
  free(acc);
  for (start = 0; start < rlen - 1 && out[start] == '0'; start++) ;
  if (start > 0) memmove(out, out + start, rlen - start);
  return rlen - start;
}

/* Multiply two signed decimal integer strings.
 * Write result to out (caller ensures at least alen+blen bytes).
 * Returns the result length (no NULL terminator written). */
STRLEN strint_mul(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen)
{
  int aneg, bneg;
  STRLEN off, rlen;
  if (alen == 0 || (alen == 1 && a[0] == '0') ||
      blen == 0 || (blen == 1 && b[0] == '0'))
    { out[0] = '0'; return 1; }
  aneg = (a[0] == '-'); if (aneg) { a++; alen--; }
  bneg = (b[0] == '-'); if (bneg) { b++; blen--; }
  off = (aneg != bneg) ? 1 : 0;
  rlen = mag_mul(out + off, a, alen, b, blen);
  if (off) out[0] = '-';
  return rlen + off;
}

STRLEN strint_muladd_s(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen, const char* c, STRLEN clen, bool negate_c)
{
  /* strint_add_s can work in place */
  STRLEN rlen = strint_mul(out, a, alen, b, blen);
  rlen = strint_add_s(out, out, rlen, c, clen, negate_c);
  return rlen;
}

/******************************************************************************/
/*                             EXPONENTIATION                                 */
/******************************************************************************/

/* a^exp (unsigned magnitudes) using binary exponentiation.
 * out must have at least limit+1 bytes.
 * Returns 0 if result would exceed limit digits; otherwise result length. */
static STRLEN mag_pow(char* out, STRLEN limit, const char* a, STRLEN alen, UV exp)
{
  char *base, *acc, *tmp;
  STRLEN baselen, acclen;

  if (exp == 1) {
    if (alen > limit) return 0;
    memcpy(out, a, alen);
    return alen;
  }
  /* a^exp < 10^(exp*alen), so result has at most exp*alen digits */
  if (exp > limit / alen) return 0;

  base = (char*) malloc(limit + 2);
  acc  = (char*) malloc(limit + 2);
  tmp  = (char*) malloc(2 * limit + 4);  /* mag_mul needs acclen+baselen bytes */

  memcpy(base, a, alen);  baselen = alen;
  acc[0] = '1';           acclen  = 1;

  while (exp > 0) {
    if (exp & 1) {
      acclen = mag_mul(tmp, acc, acclen, base, baselen);
      memcpy(acc, tmp, acclen);
    }
    exp >>= 1;
    if (exp > 0) {
      baselen = mag_mul(tmp, base, baselen, base, baselen);
      memcpy(base, tmp, baselen);
    }
  }

  memcpy(out, acc, acclen);
  free(tmp);  free(base);  free(acc);
  return acclen;
}

/* Raise a to non-negative integer power exp.
 * out must have at least limit+1 bytes (limit for magnitude + 1 for sign).
 * Returns 0 if result magnitude would exceed limit digits. */
STRLEN strint_pow(char* out, const char* a, STRLEN alen, UV exp, STRLEN limit)
{
  STRLEN off, rlen;
  int aneg;

  /* exp = 0: a^0 = 1 for all a (including 0) */
  if (exp == 0) { out[0] = '1'; return 1; }

  aneg = strint_parse(&a, &alen);

  /* 0^exp = 0 for exp > 0 */
  if (alen == 1 && a[0] == '0') { out[0] = '0'; return 1; }

  /* 1^exp = 1; (-1)^exp = +/-1 */
  if (alen == 1 && a[0] == '1') {
    if (aneg && (exp & 1)) { out[0] = '-'; out[1] = '1'; return 2; }
    out[0] = '1'; return 1;
  }

  off  = (aneg && (exp & 1)) ? 1 : 0;
  rlen = mag_pow(out + off, limit, a, alen, exp);
  if (rlen == 0) return 0;
  if (off) out[0] = '-';
  return rlen + off;
}


/******************************************************************************/
/*                                 DIVISION                                   */
/******************************************************************************/


/* Multiply unsigned decimal digit string a by single digit d (0-9).
 * out needs alen+1 bytes.  Returns result length. */
static STRLEN mag_scale(char* out, const char* a, STRLEN alen, int d)
{
  int i, carry = 0;
  if (d == 0) { out[0] = '0'; return 1; }
  for (i = (int)alen - 1; i >= 0; i--) {
    int v = (a[i] - '0') * d + carry;
    carry = v / 10;
    out[i] = '0' + v % 10;
  }
  if (!carry) return alen;
  memmove(out + 1, out, alen);
  out[0] = '0' + carry;
  return alen + 1;
}

/* Unsigned long division: a / b (both non-zero).
 * qout needs alen bytes, rout needs blen bytes.
 * Returns quotient length; sets *rlenp to remainder length. */
static STRLEN mag_div(char* qout, char* rout, STRLEN* rlenp,
                      const char* a, STRLEN alen,
                      const char* b, STRLEN blen)
{
  STRLEN i, qi, qstart, cur_len;
  int k;
  char *mult_buf, *cur, *mult[10];
  STRLEN mult_len[10];

  cur = (char*) malloc(11 * (blen+2));
  mult_buf = cur + blen+2;

  /* Precompute 0*b .. 9*b (each needs at most blen+1 digits) */
  for (k = 0; k < 10; k++) mult[k] = mult_buf + k * (blen + 2);
  mult[0][0] = '0';  mult_len[0] = 1;
  for (k = 1; k <= 9; k++) mult_len[k] = mag_scale(mult[k], b, blen, k);

  /* Sliding-window partial dividend; at most blen+1 digits */
  cur[0] = '0';  cur_len = 1;

  for (i = 0, qi = 0; i < alen; i++) {
    /* current = current * 10 + a[i] */
    if (cur_len == 1 && cur[0] == '0') cur[0] = a[i];
    else                               cur[cur_len++] = a[i];

    /* find largest k s.t. k*b <= current */
    for (k = 9; k > 0 && mag_cmp(mult[k], mult_len[k], cur, cur_len) > 0; k--)
      ;

    qout[qi++] = '0' + k;
    if (k > 0)
      cur_len = mag_sub(cur, cur, cur_len, mult[k], mult_len[k]);
  }

  /* strip leading zeros from quotient */
  for (qstart = 0; qstart < qi - 1 && qout[qstart] == '0'; qstart++) ;
  if (qstart > 0) { memmove(qout, qout + qstart, qi - qstart); qi -= qstart; }

  memcpy(rout, cur, cur_len);
  *rlenp = cur_len;
  free(cur);
  return qi;
}

/* Signed floor division and remainder of decimal integer strings.
 * Either qout or rout (and its length pointer) may be NULL.
 * qout needs alen+1 bytes; rout needs blen bytes.
 * Remainder has the sign of b (floor convention).
 * Returns false if b = 0, true otherwise. */
bool strint_divmod(char* qout, STRLEN* qlen, char* rout, STRLEN* rlen,
                   const char* a, STRLEN alen,
                   const char* b, STRLEN blen)
{
  int aneg, bneg, neg_q;
  STRLEN uqlen, urlen, off;
  char *uq, *ur;
  bool nonzero_rem;

  aneg = (alen > 0 && a[0] == '-');  if (aneg) { a++; alen--; }
  bneg = (blen > 0 && b[0] == '-');  if (bneg) { b++; blen--; }

  if (blen == 0 || (blen == 1 && b[0] == '0')) return 0;  /* b = 0 */

  neg_q = (aneg != bneg);

  /* a = 0: both results are 0 */
  if (alen == 0 || (alen == 1 && a[0] == '0')) {
    if (qout) { qout[0] = '0';  if (qlen) *qlen = 1; }
    if (rout) { rout[0] = '0';  if (rlen) *rlen = 1; }
    return 1;
  }

  if (!qout && !rout) return 1;

  /* uq offset by 1 so mag_incr carry-out has room (memmoves within buffer) */
  uq = (char*) malloc(alen+2 + blen+2);
  ur = uq + alen+2;

  uqlen = mag_div(uq + 1, ur, &urlen, a, alen, b, blen);

  /* Floor adjustment: when signs differ and remainder != 0,
   * |q_floor| = |q_trunc| + 1  and  |r_floor| = |b| - |r_trunc| */
  nonzero_rem = !(urlen == 1 && ur[0] == '0');
  if (neg_q && nonzero_rem) {
    uqlen = mag_incr(uq + 1, uqlen);
    urlen = mag_sub(ur, b, blen, ur, urlen);
  }

  if (qout) {
    off = neg_q ? 1 : 0;
    if (off) qout[0] = '-';
    memcpy(qout + off, uq + 1, uqlen);
    if (qlen) *qlen = uqlen + off;
  }

  if (rout) {
    nonzero_rem = !(urlen == 1 && ur[0] == '0');
    off = (bneg && nonzero_rem) ? 1 : 0;
    if (off) rout[0] = '-';
    memcpy(rout + off, ur, urlen);
    if (rlen) *rlen = urlen + off;
  }

  free(uq);
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
  char *rem;
  STRLEN qlen, rlen;
  rem = (char*) malloc(blen + 2);
  if (!strint_divmod(out, &qlen, rem, &rlen, a,alen, b,blen)) {
    free(rem);
    return 0;
  }
  /* ceil(a/b) = floor(a/b) + 1 when remainder is non-zero */
  if (rlen > 1 || (rlen == 1 && rem[0] != '0'))
    qlen = strint_add_s(out, out, qlen, "1", 1, 0);
  free(rem);
  return qlen;
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
 * then verifies/adjusts with strint_mul and strint_divint. */
UV strint_logint(const char* a, STRLEN alen, UV base)
{
  char base_str[24];
  STRLEN base_len, nd, i, pow_limit;
  UV k;
  char *pow_buf, *next_buf;
  STRLEN pow_len, next_len;
  int aneg;
  double approx, log10_a, k_est_f;

  if (base < 2) return UV_MAX;

  aneg = strint_parse(&a, &alen);
  if (aneg || (alen == 1 && a[0] == '0')) return UV_MAX;

  base_len = uv_to_str(base_str, base);

  /* Floating-point estimate: log10(a) via leading digits + digit count */
  nd = (alen < 15) ? alen : 15;
  approx = 0.0;
  for (i = 0; i < nd; i++)
    approx = approx * 10.0 + (double)(a[i] - '0');
  log10_a = (double)(alen - nd) + log10(approx);
  k_est_f = log10_a / log10((double)base);
  k = (k_est_f <= 1.0) ? 0 : (UV)k_est_f - 1;  /* start one below estimate */

  /* Buffers large enough for base^k (≤ a) and base^(k+1).
   * pow_limit must satisfy k ≤ pow_limit/base_len so mag_pow's conservative
   * check doesn't trigger.  (k+2)*base_len gives comfortable headroom. */
  pow_limit = (STRLEN)(k + 2) * base_len + alen + 4;
  pow_buf   = (char*) malloc(pow_limit + 2);
  next_buf  = (char*) malloc(pow_limit + base_len + 2);

  pow_len = strint_pow(pow_buf, base_str, base_len, k, pow_limit);

  /* Float overshot (very rare): step down by dividing until base^k ≤ a */
  while (k > 0 && strint_cmp(pow_buf, pow_len, a, alen) > 0) {
    pow_len = strint_divint(next_buf, pow_buf, pow_len, base_str, base_len);
    memcpy(pow_buf, next_buf, pow_len);
    k--;
  }

  /* Step up while base^(k+1) ≤ a */
  for (;;) {
    next_len = strint_mul(next_buf, pow_buf, pow_len, base_str, base_len);
    if (strint_cmp(next_buf, next_len, a, alen) > 0) break;
    memcpy(pow_buf, next_buf, next_len);
    pow_len = next_len;
    k++;
  }

  free(pow_buf);
  free(next_buf);
  return k;
}

/* floor(n^(1/k)) via Newton's method, starting from above.
 * out must have at least alen+2 bytes.
 * Returns the length of the result, or 0 on error (k==0, n<0). */
STRLEN strint_rootint(char* out, const char* a, STRLEN alen, UV k)
{
  char k_str[24], km1_str[24];
  STRLEN k_len, km1_len;
  char *r_buf, *pow_buf, *tmp_buf;
  STRLEN r_len, pow_len, tmp_len, pk1_limit, buf_size;
  STRLEN nd, i;
  double approx, log10_a, r_est_f;
  int aneg;
  UV r_digits;

  if (k == 0) return 0;

  aneg = strint_parse(&a, &alen);
  if (aneg) return 0;

  /* k==1: result is n */
  if (k == 1) { memcpy(out, a, alen); return alen; }

  /* n==0 or n==1: result is n */
  if (alen == 1 && (a[0] == '0' || a[0] == '1')) {
    out[0] = a[0]; return 1;
  }

  /* Float estimate of floor(log10(n^(1/k))) */
  nd = (alen < 15) ? alen : 15;
  approx = 0.0;
  for (i = 0; i < nd; i++)
    approx = approx * 10.0 + (double)(a[i] - '0');
  log10_a = (double)(alen - nd) + log10(approx);
  r_est_f = log10_a / (double)k;

  /* Start at 10^(floor(r_est_f)+2): guaranteed strictly above the true root.
   * The +2 absorbs floating-point error in r_est_f. */
  r_digits = (r_est_f < 0.0) ? 1 : (UV)r_est_f + 2;

  /* pk1_limit: upper bound on digits of r^(k-1).
   * r has r_digits+1 magnitude digits; r^(k-1) has at most (k-1)*(r_digits+1) digits.
   * We need pk1_limit >= (k-1)*(r_digits+1) so mag_pow's check doesn't trigger. */
  pk1_limit = (k - 1) * (r_digits + 2) + 4;
  buf_size  = pk1_limit + alen + 8;

  tmp_buf = (char*) malloc(2*buf_size + r_digits+4);
  pow_buf = tmp_buf + buf_size;
  r_buf   = tmp_buf + 2*buf_size;

  k_len   = uv_to_str(k_str,   k);
  km1_len = uv_to_str(km1_str, k - 1);

  /* Better initial r: use float to get within ~0.01% of the true root.
   * r_est_f = q + f, q integer, f in [0,1).  True root ≈ 10^f * 10^q.
   * Build r_init = r_mantissa * 10^(q-mlen+1), where r_mantissa is a
   * ceil(10^f * 1.0002) rounded up to mlen significant digits.
   * The 0.02% upward bias guarantees r_init >= true root despite float
   * error in r_est_f (actual error is < 7e-15 relative).
   * 1.0002e8 fits safely in a 32-bit UV (< 4.29e9). */
  {
    UV q = (r_est_f >= 1.0) ? (UV)r_est_f : 0;
    double frac = r_est_f - (double)q;
    UV r_mantissa = (UV)(pow(10.0, frac) * 1.0002e8) + 1;
    char mant_buf[12];   /* scratch: holds up to 10-digit mantissa safely */
    STRLEN mlen = uv_to_str(mant_buf, r_mantissa);
    if (q + 1 >= mlen) {
      /* Large result: copy mantissa into r_buf and pad with trailing zeros. */
      memcpy(r_buf, mant_buf, mlen);
      memset(r_buf + mlen, '0', q + 1 - mlen);
      r_len = q + 1;
    } else {
      /* Small result (r_est_f < 9): compute initial r directly from float.
       * pow(10, r_est_f) < 10^9 < UV_MAX on both 32- and 64-bit. */
      UV r_direct = (UV)(pow(10.0, r_est_f) * 1.0002) + 1;
      r_len = uv_to_str(r_buf, r_direct);
    }
  }

  /* Newton: r_new = floor(((k-1)*r + floor(n/r^(k-1))) / k).
   * Converges from above: each r_new < r until r = floor(n^(1/k)). */
  for (;;) {
    STRLEN new_len;

    /* r^(k-1) */
    pow_len = strint_pow(pow_buf, r_buf, r_len, k - 1, pk1_limit);
    if (pow_len == 0) break;  /* should not happen with correct pk1_limit */

    /* q = floor(n / r^(k-1)) */
    tmp_len = strint_divint(tmp_buf, a, alen, pow_buf, pow_len);

    /* numerator = (k-1)*r + q */
    pow_len = strint_mul(pow_buf, r_buf, r_len, km1_str, km1_len);
    pow_len = strint_add_s(pow_buf, pow_buf, pow_len, tmp_buf, tmp_len, 0);

    /* r_new = floor(numerator / k) */
    new_len = strint_divint(tmp_buf, pow_buf, pow_len, k_str, k_len);

    /* Converged: r is floor(n^(1/k)) when Newton no longer decreases it. */
    if (strint_cmp(tmp_buf, new_len, r_buf, r_len) >= 0) break;

    memcpy(r_buf, tmp_buf, new_len);
    r_len = new_len;
  }

  memcpy(out, r_buf, r_len);
  free(tmp_buf);
  return r_len;
}
