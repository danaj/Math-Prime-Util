#include <string.h>

#include "ptypes.h"
#include "strops.h"

/* Parse a signed decimal integer string: strip optional sign and leading
 * zeros, validate that all remaining characters are digits, and update *sp
 * and *slen to point to the canonical digit sequence.  Returns 1 if negative,
 * 0 otherwise.  Croaks on invalid input. */
static int strnum_parse(const char **sp, STRLEN *slen)
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

int strnum_cmp(const char* a, STRLEN alen, const char* b, STRLEN blen) {
  STRLEN i;
  int aneg = strnum_parse(&a, &alen);
  int bneg = strnum_parse(&b, &blen);
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
bool strnum_minmax(bool min, const char* a, STRLEN alen, const char* b, STRLEN blen)
{
  int aneg, bneg;
  STRLEN i;

  /* a is checked, process b */
  bneg = strnum_parse(&b, &blen);

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

/* Add 1 to the signed decimal integer string s/len.
 * Write result to out (caller ensures at least len+2 bytes).
 * Returns the result length (no NUL terminator written).
 * Input must be a canonical decimal integer (no leading zeros except "0"). */
STRLEN strincr(char* out, const char* s, STRLEN len)
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
STRLEN strdecr(char* out, const char* s, STRLEN len)
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
