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

/* Add 1 to the signed decimal integer string s/len.
 * Write result to out (caller ensures at least len+2 bytes).
 * Returns the result length (no NULL terminator written).
 * Input must be a canonical decimal integer (no leading zeros except "0"). */
STRLEN strint_incr(char* out, const char* s, STRLEN len)
{
  if (len == 0 || (len == 1 && s[0] == '0'))
    { out[0] = '1'; return 1; }
  if (s[0] != '-') {
    memcpy(out, s, len);
    return mag_incr(out, len);
  } else {
    s++; len--;
    if (len == 1 && s[0] == '1') { out[0] = '0'; return 1; }
    out[0] = '-';  memcpy(out + 1, s, len);
    return 1 + mag_decr(out + 1, len);
  }
}

/* Subtract 1 from the signed decimal integer string s/len.
 * Write result to out (caller ensures at least len+2 bytes).
 * Returns the result length (no NULL terminator written).
 * Input must be a canonical decimal integer (no leading zeros except "0"). */
STRLEN strint_decr(char* out, const char* s, STRLEN len)
{
  if (len == 0 || (len == 1 && s[0] == '0'))
    { out[0] = '-'; out[1] = '1'; return 2; }
  if (s[0] != '-') {
    if (len == 1 && s[0] == '1') { out[0] = '0'; return 1; }
    memcpy(out, s, len);
    return mag_decr(out, len);
  } else {
    s++; len--;
    out[0] = '-';  memcpy(out + 1, s, len);
    return 1 + mag_incr(out + 1, len);
  }
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
