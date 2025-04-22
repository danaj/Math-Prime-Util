#include <stdio.h>

#include "ptypes.h"
#include "rational.h"
#define FUNC_gcd_ui 1
#include "util.h"
#include "totients.h"


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



UV farey_length(uint32_t n)
{
  UV t = sumtotient(n);
  return (t == 0)  ?  0  :  1 + sumtotient(n);
}

int next_farey(uint32_t n, uint32_t* p, uint32_t* q)
{
  IV ivu, ivg;
  UV u, uvp, uvq;

  if (n == 0 || p == 0 || q == 0 || *p >= *q) return 0;

  ivg = gcdext( (IV)*p, (IV)*q, &ivu, 0, 0, 0);

  u = ivu;
  uvp = *p / ivg;
  uvq = *q / ivg;

  *q = ((n+u) / uvq) * uvq - u;
  *p = (*q * uvp + 1) / uvq;
  return 1;
}

UV farey_array(uint32_t n, uint32_t **rnum, uint32_t **rden)
{
  uint32_t *num, *den;
  UV i, j, p0 = 0, q0 = 1, p1 = 1, q1 = n, p2, q2;
  UV len = farey_length(n);

  if (n < 1 || len < 2 || len >= UVCONST(4294967295))
    return 0;

  New(0, num, len, uint32_t);
  New(0, den, len, uint32_t);

  for (i = 0; i < len; i++) {
    num[i] = p0;
    den[i] = q0;
    /* Haros (1802), gives p/q using two previous terms */
    j = (q0 + n) / q1;
    p2 = j * p1 - p0;
    q2 = j * q1 - q0;
    p0 = p1; q0 = q1; p1 = p2; q1 = q2;
  }
  *rnum = num;
  *rden = den;
  return len;
}



/*
 * See:
 *      Pătraşcu and Pătraşcu (ANTS 2004)
 *      https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=d8882e782674d5cd312129823287768e123674e1
 *
 *      Pawlewicz (2007)
 *      https://www.mimuw.edu.pl/~pan/papers/farey-esa.pdf
 *
 *      Pawlewicz and Pătraşcu (2008)
 *      https://www.researchgate.net/publication/225715205_Order_Statistics_in_the_Farey_Sequences_in_Sublinear_Time_and_Counting_Primitive_Lattice_Points_in_Polygons
 *
 * For the rank, we're using a very simple but fast version.
 * TODO: Use the method from Pawlewicz 2007 (see page 7).
 *
 * For the kth member, binary search on rank.
 */

UV farey_rank(uint32_t n, uint32_t p, uint32_t q)
{
  uint32_t *count, i, g;
  UV sum;

  if (n == 0 || q == 0 || p == 0) return 0;

  g = gcd_ui(p,q);
  if (g != 1) { p /= g;  q /= g; }

  New(0, count, n+1, uint32_t);

  for (i = 2; i <= n; i++)
    count[i] = ((UV)i * p - 1) / q;
  sum = 1;
  for (i = 2; i <= n; i++) {
    uint32_t j, icount = count[i];
    for (j = i; j <= n-i; j += i)
      count[j+i] -= icount;
    sum += icount;
  }
  Safefree(count);
  return sum;
}


#if 0 /* Naive method. */
int kth_farey(uint32_t n, UV k, uint32_t* p, uint32_t* q)
{
  UV i, j, p0 = 0, q0 = 1, p1 = 1, q1 = n, p2, q2;
  UV len = farey_length(n);

  if (n  > 0 && len < 2)   return -1;  /* overflow  */
  if (n == 0 || k >= len)  return 0;   /* undefined */

  if (k > len/2) {  /* Exploit symmetry about 1/2, iterate backwards */
    p0 = 1;
    p1 = n-1;
    k = (len-1)-k;
  }
  for (i = 0; i < k; i++) {
    j = (q0 + n) / q1;
    p2 = j * p1 - p0;
    q2 = j * q1 - q0;
    p0 = p1; q0 = q1; p1 = p2; q1 = q2;
  }
  *p = p0;  *q = q0;
  return 1;
}
#else
static int _walk_to_k(uint32_t a, uint32_t n, uint32_t k, uint32_t* p, uint32_t* q)
{
  uint32_t g, j, p0, q0, p1, q1, p2, q2;

  g = gcd_ui(a,n);
  p0 = a/g;
  q0 = n/g;

  if (k == 0) { *p = p0;  *q = q0;  return 1; }

  /* From the single point, use extgcd to get the exact next fraction */
  p1 = p0; q1 = q0;
  next_farey(n, &p1, &q1);

  /* Now we have two fractions, so quick step through */
  while (--k) {
    j = (q0 + n) / q1;
    p2 = j * p1 - p0;
    q2 = j * q1 - q0;
    p0 = p1; p1 = p2;  q0 = q1; q1 = q2;
  }
  *p = p1;
  *q = q1;
  return 1;
}
int kth_farey(uint32_t n, UV k, uint32_t* p, uint32_t* q)
{
  uint32_t lo = 1, hi = n;
  UV cnt = 1;

  if (k < 2) {
    if (k == 0) { *p = 0;  *q = 1; }
    else        { *p = 1;  *q = n; }
    return 1;
  }
  if (n < 2) return 0;

  /* For a substantial performance benefit, we will estimate the position
   * and get its rank.  Then look a small distance in the other direction.
   * For small n this often completely brackets the value after only one
   * or two calls.  For large n, we can often do 2-5 times fewer calls.
   *
   * The downside is this is ugly, but it makes this call 2-4x faster.
   */
  if (n >= 5) {
    uint32_t const ginc = ((UV)n+8191)>>13;
    double   const dlen  = 1 + (0.304*(double)n*n + .29*(double)n + 0.95);
    uint32_t guess = k * ((double)n/dlen);
    UV gcnt = 0;
    if (guess <= lo) guess = lo+1; else if (guess >= hi) guess = hi-1;

    if (lo < hi) {
      gcnt = farey_rank(n, guess, n);
      if (gcnt <= k) { lo = guess; cnt = gcnt; } else { hi = guess-1; }
    }

    /* Look a small distance in the other direction.  We want it to be
     * far enough that we bracket the value, but not so far that we make
     * too many calls getting back. */
    if (gcnt <= k) { guess = (hi-ginc   < guess) ? hi : guess+ginc; }
    else           { guess = (lo+ginc+1 > guess) ? lo : guess-1-ginc; }

    if (lo < hi && guess > lo && guess < hi) {
      gcnt = farey_rank(n,guess,n);
      if (gcnt <= k) { lo = guess; cnt = gcnt; } else { hi = guess-1; }
    }
  }

  /* Now the binary search. */

  while (lo < hi) {
    uint32_t mid = lo + ((hi-lo+1)>>1);
    UV midcnt = farey_rank(n, mid, n);
    if (midcnt <= k)  { lo = mid;  cnt = midcnt; }
    else              { hi = mid-1; }
  }
  if (lo == n) {
    *p = (cnt == k);
    *q = 1;
    return k <= cnt;
  }
  return _walk_to_k(lo, n, k-cnt, p, q);
}
#endif
