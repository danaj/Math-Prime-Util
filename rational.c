#include <stdio.h>

#include "ptypes.h"
#include "rational.h"
#define FUNC_gcd_ui 1
#include "util.h"


int contfrac(UV** cfrac, UV *rem, UV num, UV den)
{
  UV *cf;
  int steps = 0;
  New(0, cf, 2 * BITS_PER_WORD, UV);  /* Upper limit for gcd steps */
  while (den > 0) {
    UV q = num/den;
    UV r = num - q*den;
    num = den;
    den = r;
    cf[steps++] = q;
  }
  *rem = num;
  *cfrac = cf;
  return steps;
}

int next_calkin_wilf(UV* num, UV* den)
{
  UV n, d;
  if (num == 0 || den == 0) return 0;
  n = *num;
  d = *den;
  if (n == 0 || d == 0 || gcd_ui(n,d) != 1) return 0;

  /* next = d / (n+d-2*(n%d)) = d / (2(n/d)+1)*d-n */
  if (n < d) {  /* n/d is 0 */
    *den = d-n;
  } else if (d == 1) {
    if (n == UV_MAX) return 0;
    *den = n + 1;
  } else {      /* n >= d and d >= 2 */
    UV nd = n % d;  /* nd is less than d and less than n */
    UV nr = n-nd,  dr = d-nd;
    if (nr > UV_MAX-dr)  return 0;
    *den = nr + dr;
  }
  *num = d;
  return 1;
}
int next_stern_brocot(UV* num, UV* den)
{
  UV n, d;
  if (num == 0 || den == 0) return 0;
  n = *num;
  d = *den;
  if (n == 0 || d == 0 || gcd_ui(n,d) != 1) return 0;

  /* Backhouse and Ferreira show how to do this *if* we had a 2x2 matrix
   * for the node.  We could also exploit that given a/b and the next c/d
   *     bc-ad=3 if they share a parent
   * but this doesn't give us enough information to solve for both c,d.
   */

  if (*den == 1) { /* At end of the row, go to the start of the next. */
    if (*num == UV_MAX) return 0;
    *den = *num+1;
    *num = 1;
    return 1;
  }
  /* Given the tree e.g. LLLRRLLRR, we can go up to the nearest ancestor,
   * then back down.  That is, from the right, invert all L/R from the end
   * to and including the right L.  This really isn't a huge savings over
   * doing the full process.  Doing nth(n(F)+1) is clean. */
  return nth_stern_brocot(num, den, 1+stern_brocot_n(*num, *den));
}


#if 0  /* A recursive version */
UV calkin_wilf_n(UV num, UV den)
{
  if (num == den) {
    return 1;
  } else if (num > den) {
    UV f = calkin_wilf_n(num-den, num);
    if (f == 0 || f == UV_MAX) return 0;
    return 1 + f;
  } else {
    UV f = calkin_wilf_n(num, den-num);
    if (f == 0 || f > (UV_MAX/2)) return 0;
    return 2 * f;
  }
}
#endif

UV calkin_wilf_n(UV num, UV den)
{
  UV *cf = 0, n = 0, rem;
  uint32_t bit, d = 1, shift = 0;
  int i, steps = contfrac(&cf, &rem, num, den);

  if (rem != 1)   croak("Rational must be reduced");
  if (steps == 0) return 0;

  cf[steps-1]--;
  for (i = 0; i < steps; i++) {
    if ((shift+cf[i]) >= BITS_PER_WORD)
      break;
    if (d)
      for (bit = 0; bit < cf[i]; bit++)
        n |= UVCONST(1) << (shift+bit);
    shift += cf[i];
    d ^= 1;  /* d = 1-d; */
  }
  Safefree(cf);
  if (i < steps)  return 0;
  n |= UVCONST(1) << shift;
  return n;
}
UV stern_brocot_n(UV num, UV den)
{
  /* Reverse bits in the Calkin-Wilf n */
  UV n, M = calkin_wilf_n(num,den);
  if (M == 0) return 0;
  for (n = 1; M > 1; M >>= 1)
    n = (n << 1) | (M & 1);
  return n;
}


int nth_calkin_wilf(UV* num, UV* den, UV n)
{
  uint32_t b = 1;
  UV p = 0, q = 1;    /* p odd  q even */
  { UV v = n; while (v >>= 1) b++; }
  while (b--) {
    if ((n >> b) & 1)  p += q;
    else               q += p;
  }
  *num = p;
  *den = q;
  return 1;
}
int nth_stern_brocot(UV* num, UV* den, UV n)
{
  UV p = 1, q = 1;    /* p odd  q even */
  while (n > 1) {
    if (n & 1)  p += q;
    else        q += p;
    n >>= 1;
  }
  *num = p;
  *den = q;
  return 1;
}

UV nth_stern_diatomic(UV n)
{
  UV p = 0, q = 1;
  while (n) {
    if (n & 1)  p += q;
    else        q += p;
    n >>= 1;
  }
  return p;
}
