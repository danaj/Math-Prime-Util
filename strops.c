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
 * Returns the result length (no NUL terminator written).
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
 * Returns the result length (no NUL terminator written).
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

/* Add two signed decimal integer strings.
 * Write result to out (caller ensures at least max(alen,blen)+1 bytes).
 * Returns the result length (no NUL terminator written). */
STRLEN strint_add(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen)
{
  int aneg, bneg;
  if (alen == 0 || (alen == 1 && a[0] == '0')) { memcpy(out, b, blen); return blen; }
  if (blen == 0 || (blen == 1 && b[0] == '0')) { memcpy(out, a, alen); return alen; }
  aneg = (a[0] == '-'); if (aneg) { a++; alen--; }
  bneg = (b[0] == '-'); if (bneg) { b++; blen--; }
  return str_add_impl(out, a, alen, aneg, b, blen, bneg);
}

/* Subtract b from a (signed decimal integer strings).
 * Write result to out (caller ensures at least max(alen,blen)+1 bytes).
 * Returns the result length (no NUL terminator written). */
STRLEN strint_sub(char* out, const char* a, STRLEN alen, const char* b, STRLEN blen)
{
  int aneg, bneg;
  if (blen == 0 || (blen == 1 && b[0] == '0')) { memcpy(out, a, alen); return alen; }
  if (alen == 0 || (alen == 1 && a[0] == '0')) {
    /* 0 - b = -b */
    if (b[0] == '-') { memcpy(out, b + 1, blen - 1); return blen - 1; }
    out[0] = '-'; memcpy(out + 1, b, blen); return blen + 1;
  }
  aneg = (a[0] == '-'); if (aneg) { a++; alen--; }
  bneg = (b[0] == '-'); if (bneg) { b++; blen--; }
  return str_add_impl(out, a, alen, aneg, b, blen, !bneg);
}
