#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * A nice discussion of congruent numbers can be found in:
 *
 * https://pub.math.leidenuniv.nl/~stevenhagenp/ANTproc/19yui.pdf
 *
 */

#include "ptypes.h"
#include "congruent_numbers.h"
#include "util.h"
#include "util_bits.h"
#include "util_math.h"
#include "factor.h"
#include "rootmod.h"

/******************************************************************************/

/* We only look at the non-square portion of n. */
static void remove_square_part(factored_t *nf)   /* Turn n*c^2 into n */
{
  if (nf->n > 3) {
    uint16_t i, j;
    for (i = 0; i < nf->nfactors; i++)
      if (nf->e[i] > 1)
        break;
    if (i < nf->nfactors) {
      UV N = 1;
      for (i = 0, j = 0; i < nf->nfactors; i++)
        if (nf->e[i] & 1) {
          N *= nf->f[i];
          nf->e[j] = 1;
          nf->f[j++] = nf->f[i];
        }
      nf->n = N;
      nf->nfactors = j;
    }
  }
  /* factoredp_validate(nf); */
}

/******************************************************************************/


/* Tunnell's method, counting integer solutions to ternary quadratics.
 * Assumes the weak BSD conjecture.
 * Weak BSD holds for n <   42553.   (Nemenzo 1998)
 * Weak BSD holds for n <   71474.   (Wado 2005)
 * Weak BSD holds for n <  300000.   (Matsuno 2005]
 * Weak BSD holds for n < 1000000.   (Matsuno 2006 + Elkies 2002)
 *
 * The final answer.
 * Most of the rest of this file is trying to avoid calling this if possible.
 */
static bool _is_congruent_number_tunnell(UV n)
{
  UV x, y, z, limz, limy, limx, n8z, zsols, sols[2] = {0,0};

  /* The input MUST be square-free or the result will not be correct. */

  if (n&1) {
    for (z = 0, limz = isqrt(n/8);  z <= limz;  z++) {
      zsols = 0;
      n8z = n - 8*z*z;
      for (y = 0, limy = isqrt(n8z/2);  y <= limy;  y++) {
        x = n8z - 2*y*y;   /* n odd => n8z odd => x odd */
        if (is_perfect_square(x))
          zsols += 1 << (1+(y>0)+(z>0));
      }
      sols[z&1] += zsols;
    }
  } else {
    for (z = 0, limz = isqrt((n/2)/8);  z <= limz;  z++) {
      zsols = 0;
      n8z = n/2 - 8*z*z;
      for (x = 1, limx = isqrt(n8z);  x <= limx;  x += 2) {
        y = n8z - x*x;
        if (y == 0 || is_perfect_square(y))
          zsols += 1 << (1+(y>0)+(z>0));
      }
      sols[z&1] += zsols;
    }
  }
  return (sols[0] == sols[1]);
}

/******************************************************************************/

#define SWAP4(x,y) { UV t; t=x; x=y; y=t;  t=x##8; x##8=y##8; y##8=t; }
#define KPQ kronecker_uu(p,q)
#define KQP kronecker_uu(q,p)
/* For 3 factors */
#define KPR kronecker_uu(p,r)
#define KQR kronecker_uu(q,r)
#define KRP kronecker_uu(r,p)
#define KRQ kronecker_uu(r,q)
/* For 4 factors */
#define KPS kronecker_uu(p,s)
#define KQS kronecker_uu(q,s)
#define KRS kronecker_uu(r,s)
#define KSP kronecker_uu(s,p)
#define KSQ kronecker_uu(s,q)
#define KSR kronecker_uu(s,r)

#define LAGRANGE_COND1 (KPQ==KPR || KQR==KQP || KRP==KRQ)
#define LAGRANGE_COND2 ((KPQ==1 && KPR==1) || (KQR==1 && KQP==1) || (KRP==1 && KRQ==1))
#define LAGRANGE_COND3 ((KPQ==-1 && KPR==-1) || (KQR==-1 && KQP==-1) || (KRP==-1 && KRQ==-1))

static bool _can_order_kronecker2(UV p, UV q) {
  return KPQ==-1 || KQP == -1;
}
static bool _can_order_kronecker3(UV p, UV q, UV r) {
  return (KPR==-1 && KQR==-1 && _can_order_kronecker2(p,q)) ||
         (KPQ==-1 && KRQ==-1 && _can_order_kronecker2(p,r)) ||
         (KQP==-1 && KRP==-1 && _can_order_kronecker2(q,r));
}
static bool _can_order_kronecker4(UV p, UV q, UV r, UV s) {
  return (KPS==-1 && KQS==-1 && KRS==-1 && _can_order_kronecker3(p,q,r)) ||
         (KPR==-1 && KQR==-1 && KSR==-1 && _can_order_kronecker3(p,q,s)) ||
         (KPQ==-1 && KRQ==-1 && KSQ==-1 && _can_order_kronecker3(p,r,s)) ||
         (KQP==-1 && KRP==-1 && KSP==-1 && _can_order_kronecker3(q,r,s));
}

/* Returns -1 if not known, 0 or 1 indicate definite results. */
static int _is_congruent_number_filter1(const factored_t *nf) {
  const UV n         = nf->n;
  const bool isodd   = n & 1;
  const bool iseven  = !isodd;

  UV p,  q,  r,  s;   /* Four odd factors in mod 8 order */
  UV p8, q8, r8, s8;  /* values mod 8 */
  UV fac[MPU_MAX_DFACTORS] = {0};  /* The odd factors, in mod 8 order */
  uint32_t nfac = 0;

  MPUassert(n >= 13, "n too small in icn_filter");

  /* The ACK conjecture (Alter, Curtz, and Kubota 1972):
   *     n = {5,6,7} mod 8   =>  n is a congruent number
   * also follows from the weak BSD conjecture.
   */
  if ((n % 8) >= 5) return 1;

  /* No filter here handles more than 4 odd factors */
  if (nf->nfactors-iseven > 4)  return -1;

  { /* Sort odd factors into fac array by mod 8 */
    uint32_t i;
    for (i=0; i<nf->nfactors; i++) if (nf->f[i] % 8 == 1) fac[nfac++] = nf->f[i];
    for (i=0; i<nf->nfactors; i++) if (nf->f[i] % 8 == 3) fac[nfac++] = nf->f[i];
    for (i=0; i<nf->nfactors; i++) if (nf->f[i] % 8 == 5) fac[nfac++] = nf->f[i];
    for (i=0; i<nf->nfactors; i++) if (nf->f[i] % 8 == 7) fac[nfac++] = nf->f[i];
  }
  p = fac[0];
  q = (nfac > 1) ? fac[1] : 0;
  r = (nfac > 2) ? fac[2] : 0;
  s = (nfac > 3) ? fac[3] : 0;
  p8 = p % 8;  q8 = q % 8;  r8 = r % 8;  s8 = s % 8;

  /* Evink 2021    https://arxiv.org/pdf/2105.01450.pdf
   * Feng 1996     http://matwbn.icm.edu.pl/ksiazki/aa/aa75/aa7513.pdf
   * Monsky 1990   https://gdz.sub.uni-goettingen.de/id/PPN266833020_0204
   * Lagrange 1974 https://www.numdam.org/item/SDPP_1974-1975__16_1_A11_0.pdf
   */

  if (isodd && nfac == 1) {          /* n = p */

    UV root;
    if (p8 == 3) return 0;  /* Genocchi 1855 */

    /* https://arxiv.org/pdf/2105.01450.pdf, Prop 2.1.2 */
    if (sqrtmodp(&root, 2, p) && kronecker_uu(1+root, p) == -1)
      return 0;
#if 0   /* Requires mulmod.h */
    { /* Evink 2021 shows these are equivalent to the sqrt test above */
      UV a,b;
      if (1 && cornacchia(&a, &b, 1, p)) {
        if (p != (a*a+b*b)) croak("bad corn for %lu\n",p);
        if (sqrmod(a+b,16) != 1)
          { printf("ret\n"); return 0; }
      }
      if (1 && cornacchia(&a, &b, 4, p))
        if (kronecker_uu(a+2*b, p) == -1)
          { printf("ret 2\n"); return 0; }
    }
#endif

  } else if (iseven && nfac == 1) {  /* n = 2p */

    if (p8 == 3 || p8 == 7)     return 1;  /* we already returned 1 earlier */
    if (p8 == 5)                return 0;  /* Genocchi 1855 */
    if (p % 16 == 9)            return 0;  /* Bastien 1915 */

  } else if (isodd && nfac == 2) {  /* n = pq */

    if (p8 == 3 && q8 == 3) return 0;  /* Genocchi 1855 */
#if 0  /* Monsky, all produce n mod 8 = 5 or 7: we already returned 1 */
    if (p8 == 3 && q8 == 7) return 1;
    if (p8 == 3 && q8 == 5) return 1;
    if (p8 == 1 && q8 == 5 && KPQ == -1) return 1;
    if (p8 == 1 && q8 == 7 && KPQ == -1) return 1;
#endif
    /* Lagrange 1974 */
    if (p8 == 1 && q8 == 3 && KPQ == -1) return 0;
    if (p8 == 5 && q8 == 7 && KPQ == -1) return 0;

  } else if (iseven && nfac == 2) {  /* n = 2pq */

    if (p8 == 5 && q8 == 5) return 0; /* Genocchi 1855 */
#if 0  /* Monsky, all produce n mod 8 = 6: we already returned 1 */
    if (p8 == 3 && q8 == 5) return 1;
    if (p8 == 5 && q8 == 7) return 1;
    if (p8 == 1 && q8 == 7 && KPQ == -1) return 1;
    if (p8 == 1 && q8 == 3 && KPQ == -1) return 1;
#endif
    /* Lagrange 1974 */
    if (p8 == 3 && q8 == 3) return 0;
    if (p8 == 1 && q8 == 5 && KPQ == -1) return 0;
    if (p8 == 3 && q8 == 7 && KPQ == -1) return 0;
    if (p8 == 7 && q8 == 7 && ((q % 16 == 7 && KPQ == 1) ||
                               (p % 16 == 7 && KQP == 1) )) return 0;
    if (p8 == 1 && q8 == 1 && KPQ == -1 && (p*q) % 16 == 9) return 0;

  } else if (isodd && nfac == 3) {  /* n = pqr */

#if 0 /* Serf 1991, all produce n mod 8 = 5 or 7: we already returned 1 */
    if (p8 == 3 && q8 == 3 && r8 == 5) return 1;
    if (p8 == 3 && q8 == 3 && r8 == 7) return 1;
    if (p8 == 7 && q8 == 7 && r8 == 7 && KPQ == -KPR && KPQ == KQR) return 1;
#endif
    /* Lagrange 1974 */
    if (p8 == 1 && q8 == 3 && r8 == 3 && KPQ == -KPR) return 0;
    if (p8 == 3 && q8 == 5 && r8 == 7 && KQR == -1) return 0;
    if (p8 == 3 && q8 == 7 && r8 == 7 && KPQ == -KPR && KPQ == KQR) return 0;
    if (p8 == 1 && q8 == 1 && r8 == 3 && LAGRANGE_COND3) return 0;
    if (p8 == 1 && q8 == 5 && r8 == 7 && LAGRANGE_COND3) return 0;
    if (p8 == 3 && q8 == 5 && r8 == 5 && LAGRANGE_COND3) return 0;
    if (p8 == 3 && q8 == 3 && r8 == 3 && LAGRANGE_COND1) return 0;
    if (p8 == 1 && q8 == 1 && r8 == 1 && LAGRANGE_COND3) {
      UV c,d;
      if (cornacchia(&c, &d, 8, n) && d&1)
        return 0;
    }

  } else if (iseven && nfac == 3) {  /* n = 2pqr */

#if 0 /* Serf 1991, all produce n mod 8 = 6: we already returned 1 */
    if (p8 == 3 && q8 == 3 && r8 == 7) return 1;
    if (p8 == 3 && q8 == 5 && r8 == 5) return 1;
    if (p8 == 5 && q8 == 5 && r8 == 7) return 1;
    if (p8 == 7 && q8 == 7 && r8 == 7 && KPQ == -KPR && KPQ == KQR) return 1;
#endif
    /* Lagrange 1974 */
    if (p8 == 1 && q8 == 3 && r8 == 3 && KPQ == -KPR) return 0;
    if (p8 == 1 && q8 == 5 && r8 == 5 && KPQ == -KPR) return 0;
    if (p8 == 3 && q8 == 5 && r8 == 7 && KRP == KRQ) return 0;
    if (p8 == 1 && q8 == 1 && r8 == 1 && LAGRANGE_COND3 && (p*q*r) % 16 == 9) return 0;
    if (p8 == 5 && q8 == 7 && r8 == 7 && KQP == KQR && KQP == -KRP) return 0;
    if (p8 == 1 && q8 == 1 && r8 == 5 && LAGRANGE_COND3) return 0;
    /* Lagrange's 1 3 7 seems to be incorrect.
     * 13706 = 2*7*11*89 = 2*89*7*11, so p = 1, q = -1, r = 3 mod 8.
     *       cond3 (q|r)= (q|p) = -1.
     * but 13706 is a congruent number.
     * Cheng/Guo 2018 show this case with a third congruency.
     */
    if (p8 == 3 && q8 == 3 && r8 == 5 && LAGRANGE_COND1) return 0;
    if (p8 == 5 && q8 == 5 && r8 == 5 && LAGRANGE_COND2) return 0;
    /* Cheng/Guo 2018, Theorem 1.3.6 */
    if (p8 == 1 && q8 == 3 && r8 == 7 && KPQ == -1 && KPR == -1 && KQR == -1) return 0;

  } else if (isodd && nfac == 4) {  /* n = pqrs */

    /* Serf 1991 */
    if (p8 == 5 && q8 == 5 && r8 == 7 && s8 == 7 &&
        ( (KPR ==  1 && KQR == -1 && KPS == -1) ||
          (KPR == -1 && KPS ==  1 && KQS == -1) ||
          (KPR == -1 && KPS == -1 && KQR == -KQS)))
      return 0;
    /* Cheng/Guo 2018, Theorem 1.3.1 */
    if (p8 == 3 && q8 == 3 && r8 == 5 && s8 == 7 &&
        (KPQ == -1 || KQP == -1) &&
        KPR == -1 && KQR == -1 && KPS == -1 && KQS == -1 && KRS == -1)
      return 0;
    /* Cheng/Guo 2018, Theorem 1.2.5 */
    if (p8 == 1 && q8 == 3 && r8 == 3 && s8 == 3 &&
        KPQ == -1 && KPR == -1 && KPS == -1 &&
        _can_order_kronecker3(q,r,s))
      return 0;
    /* Iskra 1996 (also Cheng/Guo 2018, Theorem 1.1.2) */
    if (p8 == 3 && q8 == 3 && r8 == 3 && s8 == 3 &&
        _can_order_kronecker4(p,q,r,s))
      return 0;
    /* Das 2020 4 factors */
    if (p8 == 1 && q8 == 3 && r8 == 5 && s8 == 7 &&
        (KQS == -1 || KSQ == -1) &&
         KPR == 1 && KRP == 1 && KPS == 1 && KRQ == 1 && KPQ == -1 && KRS == -1)
      return 0;
    /* Das 2020 4 factors */
    if (p8 == 1 && q8 == 1 && r8 == 3 && s8 == 3 &&
        (KRS == -1 || KSR == -1) && KPQ == 1 && KQP == 1 &&
        ( (KPR == -1 && KQS == -1 && KPS == 1 && KQR == 1) ||
          (KQR == -1 && KPS == -1 && KQS == 1 && KPR == 1) ))
      return 0;

  } else if (iseven && nfac == 4) {  /* n = 2pqrs */

    /* Serf 1991 */
    if (p8 == 1 && q8 == 1 && r8 == 3 && s8 == 3 &&
        ( (KPQ ==  1 && KPR == -KPS && KQR == -KQS) ||
          (KPQ == -1 && KPR ==  KPS && KQR == -KQS) ||
          (KPQ == -1 && KPR == -KPS)))
      return 0;
    /* Cheng/Guo 2018, Theorem 1.3.2 */
    if (p8 == 3 && q8 == 5 && r8 == 5 && s8 == 7 &&
        (KQR == -1 || KRQ == -1) &&
        KPQ == -1 && KPR == -1 && KPS == -1 && KQS == -1 && KRS == -1)
      return 0;
    /* Cheng/Guo 2018, Theorem 1.3.8 */
    if (p8 == 1 && q8 == 3 && r8 == 3 && s8 == 5 &&
        (KQR == -1 || KRQ == -1) &&
        KPQ == -1 && KPR == -1 && KPS == -1 && KQS == -1 && KRS == -1)
      return 0;
    /* Cheng/Guo 2018, Theorem 1.2.3 */
    if (p8 == 3 && q8 == 3 && r8 == 3 && s8 == 7 &&
        KPS == -1 && KQS == -1 && KRS == -1 &&
        _can_order_kronecker3(p,q,r))
      return 0;
    /* Cheng/Guo 2018, Theorem 1.2.6 */
    if (p8 == 1 && q8 == 5 && r8 == 5 && s8 == 5 &&
        KPQ == -1 && KPR == -1 && KPS == -1 &&
        _can_order_kronecker3(q,r,s))
      return 0;
    /* Cheng/Guo 2018, Theorem 1.1.1 (even analog of Iskra 1996) */
    if (p8 == 3 && q8 == 3 && r8 == 3 && s8 == 3 &&
        _can_order_kronecker4(p,q,r,s))
      return 0;
    /* Cheng/Guo 2018, Theorem 1.1.3 */
    if (p8 == 5 && q8 == 5 && r8 == 5 && s8 == 5 &&
        _can_order_kronecker4(p,q,r,s))
      return 0;
    /* Cheng/Guo 2018, Theorem 1.2.1 */
    if (p8 == 3 && q8 == 3 && r8 == 5 && s8 == 5 &&
        (KPQ == -1 || KQP == -1) && (KRS == -1 || KSR == -1) &&
        KPR == -1 && KQR == -1 && KPS == -1 && KQS == -1)
      return 0;

  }

  return -1;
}

/******************************************************************************/

/* This has more complicated filters that take arbitrary numbers of factors,
 * and have to handle permutations.
 */

/* One bit per distinct odd prime factor (at most 15).  The diagonal is unused. */
typedef struct {
  const UV *fac;
  uint16_t known[MPU_MAX_DFACTORS];
  uint16_t neg[MPU_MAX_DFACTORS];
} cn_symbols;

static bool _cn_negative(cn_symbols *K, unsigned i, unsigned j) {
  const uint16_t jbit = (uint16_t)(1U << j);
  if (!(K->known[i] & jbit)) {
    const uint16_t ibit = (uint16_t)(1U << i);
    UV p = K->fac[i], q = K->fac[j];
    bool neg = kronecker_uu(p,q) == -1;
    if (neg) K->neg[i] |= jbit;
    /* Reciprocity reverses the sign only when both primes are 3 mod 4. */
    if (neg != ((p & q & 2) != 0)) K->neg[j] |= ibit;
    K->known[i] |= jbit;
    K->known[j] |= ibit;
  }
  return (K->neg[i] & jbit) != 0;
}

static unsigned _cn_outdegree(cn_symbols *K, unsigned group, unsigned i) {
  unsigned missing = group & ~(K->known[i] | (1U << i));
  for (; missing; missing &= missing-1)
    _cn_negative(K, i, ctz(missing));
  return popcnt(K->neg[i] & group);
}

/* For primes all 3 mod 4, each pair forces one direction.  The resulting
 * tournament can be ordered iff its outdegrees are all distinct. */
static bool _cn_can_order(cn_symbols *K, unsigned group) {
  unsigned left, degrees = 0;
  if (popcnt(group) <= 2) return TRUE;
  for (left = group; left; left &= left-1) {
    unsigned bit = 1U << _cn_outdegree(K, group, ctz(left));
    if (degrees & bit) return FALSE;
    degrees |= bit;
  }
  return TRUE;
}

/* Can one reversal give an order with an odd gap between its endpoints?
 * Only the two endpoint degrees change, and in an order the degrees give
 * the positions counted from the end. */
static bool _cn_can_order_one_flip(cn_symbols *K, unsigned group) {
  unsigned left, edges, rest;
  unsigned degree[MPU_MAX_DFACTORS];

  for (left = group; left; left &= left-1) {
    unsigned i = ctz(left);
    degree[i] = _cn_outdegree(K, group, i);
  }
  for (left = group; left; left &= left-1) {
    unsigned i = ctz(left);
    for (edges = K->neg[i] & group; edges; edges &= edges-1) {
      unsigned j = ctz(edges);
      unsigned di = degree[i]-1, dj = degree[j]+1;
      unsigned used;
      if (dj <= di || !((dj-di) & 1)) continue;
      used = (1U << di) | (1U << dj);
      for (rest = group & ~((1U << i) | (1U << j)); rest; rest &= rest-1) {
        unsigned bit = 1U << degree[ctz(rest)];
        if (used & bit) break;
        used |= bit;
      }
      if (!rest) return TRUE;
    }
  }
  return FALSE;
}

/* Cheng / Guo, "The non-congruent numbers via Monsky's formula".
 * Test the residue counts before computing any Legendre symbols. */
static bool _cheng_guo_shape(unsigned eps, const unsigned g[8]) {
  /* Theorem 1.1 */
  if ( (g[1] == 0 && g[5] == 0 && g[7] == 0 && eps == 2 && g[3] % 2 == 0)
    || (g[1] == 0 && g[5] == 0 && g[7] == 0 && eps == 1)
    || (g[1] == 0 && g[3] == 0 && g[7] == 0 && eps == 2 && g[5] % 2 == 0))
    return TRUE;

  /* Theorem 1.2 */
  if ( (g[1] == 0 && g[7] == 0 && eps == 1 && (g[3] % 2) == 1 && g[5] >= 1 && g[5] % 2 == 0)
    || (g[1] == 0 && g[7] == 0 && eps == 2 && g[3] >= 1 && g[5] >= 1 && g[3] % 2 == 0 && g[5] % 2 == 0)
    || (g[1] == 0 && g[7] == 0 && eps == 2 && g[3] >= 1 && g[5] == 1 && g[3] % 2 == 0)
    || (g[1] == 0 && g[5] == 0 && eps == 2 && g[7] == 1 && g[3] % 2 == 1)
    || (g[1] == 0 && g[3] == 0 && eps == 1 && g[5] == 1 && g[7] == 1)
    || (g[5] == 0 && g[7] == 0 && eps == 1 && g[1] > 1 && g[1] % 2 == 0 && g[3] % 2 == 1)
    || (g[5] == 0 && g[7] == 0 && eps == 1 && g[1] == 1 && g[3] % 2 == 1)
    || (g[3] == 0 && g[7] == 0 && eps == 2 && g[1] == 1 && g[5] % 2 == 1)
    || (g[3] == 0 && g[7] == 0 && eps == 2 && g[5] == 1 && g[1] > 1 && g[1] % 2 == 0) )
    return TRUE;

  /* Theorem 1.3 */
  if ( (g[1] == 0 && eps == 1 && g[5] == 1 && g[7] == 1 && g[3] >= 1)
    || (g[1] == 0 && eps == 1 && g[5] >= 2 && g[7] == 1 && g[3] % 2 == 1 && g[5] % 2 == 1)
    || (g[1] == 0 && eps == 2 && g[5] >= 1 && g[7] == 1 && g[3] % 2 == 1 && g[5] % 2 == 0)
    || (g[3] == 0 && eps == 1 && g[7] == 1 && g[1] % 2 == 1 && g[5] % 2 == 1)
  /* No examples 1.3.4 */
    || (g[3] == 0 && eps == 2 && g[5] == 1 && g[7] == 1 && g[1] % 2 == 1)
  /* No examples 1.3.5 */
    || (g[5] == 0 && eps == 1 && g[1] >= 1 && g[3] >= 1 && g[7] == 1 && g[1] % 2 == 0 && g[3] % 2 == 0)
    || (g[5] == 0 && eps == 2 && g[1] >= 1 && g[3] >= 1 && g[7] == 1 && g[1] % 2 == 1 && g[3] % 2 == 1)
    || (g[7] == 0 && eps == 1 && g[1] >= 1 && g[5] >= 1 && g[1] % 2 == 0 && g[5] % 2 == 0 && g[3] % 2 == 1)
    || (g[7] == 0 && eps == 2 && g[1] == 1 && g[3] >= 1 && g[3] % 2 == 0 && g[5] % 2 == 1)
    || (g[7] == 0 && eps == 2 && g[1] >= 1 && g[3] >= 1 && g[5] == 1 && g[1] % 2 == 0 && g[3] % 2 == 0) )
    return TRUE;

  /* Theorem 1.4 */
  if ( (eps == 1 && g[1] >= 1 && g[7] == 1 && g[1] % 2 == 0 && g[3] % 2 == 1 && g[5] % 2 == 1)
    || (eps == 1 && g[3] >= 1 && g[7] == 1 && g[1] % 2 == 1 && g[5] % 2 == 1 && g[3] % 2 == 0) )
    return TRUE;

  return FALSE;
}

static bool _cheng_guo_symbols(cn_symbols *K, unsigned n, const unsigned groups[8]) {
  unsigned i, j;
  /* Across residue groups the direction is fixed.  Within the 1 and 5
   * groups reciprocity makes the direction irrelevant. */
  for (i = 0; i < n; i++) {
    unsigned ri = K->fac[i] & 7;
    for (j = i+1; j < n; j++) {
      unsigned rj = K->fac[j] & 7;
      if (ri == rj && (ri & 2)) continue;
      if (!(ri <= rj ? _cn_negative(K, i, j) : _cn_negative(K, j, i)))
        return FALSE;
    }
  }
  return _cn_can_order(K, groups[3]) && _cn_can_order(K, groups[7]);
}

/* Das / Saikia 2020: https://math.colgate.edu/~integers/u55/u55.pdf */
static bool _das_saikia_symbols(cn_symbols *K, const unsigned groups[8]) {
  const unsigned qs = groups[3] | groups[7];
  unsigned ps, left, used = 0;

  for (ps = groups[1] | groups[5]; ps; ps &= ps-1) {
    unsigned i = ctz(ps), match = 0;
    unsigned compatible = groups[(K->fac[i] & 7) + 2];

    /* All p-p symbols must be positive, in either direction. */
    for (left = ps & (ps-1); left; left &= left-1)
      if (_cn_negative(K, i, ctz(left))) return FALSE;

    /* One negative q per p, with the right residue and no shared matches. */
    for (left = qs; left; left &= left-1) {
      unsigned j = ctz(left), bit = 1U << j;
      if (!_cn_negative(K, i, j)) continue;
      if (match || !(compatible & bit) || (used & bit)) return FALSE;
      match = bit;
    }
    if (!match) return FALSE;
    used |= match;
  }
  /* Order all q's together; the two residue classes may interleave. */
  return _cn_can_order(K, qs);
}

/* Returns -1 if not known, 0 or 1 indicate definite results. */
static int _is_congruent_number_filter2(const factored_t *nf) {
  const unsigned iseven = !(nf->n & 1);
  const unsigned noddfac = nf->nfactors - iseven;
  const UV *oddfac = nf->f + iseven;
  unsigned i, g[8] = {0}, groups[8] = {0};
  bool allmod3, cheng, das;
  cn_symbols K;

  for (i = 0; i < noddfac; i++) {
    unsigned r = oddfac[i] & 7;
    g[r]++;
    groups[r] |= 1U << i;
  }
  allmod3 = g[3] == noddfac;
  cheng = !allmod3 && _cheng_guo_shape(1 + iseven, g);
  das = !iseven && noddfac >= 4 && g[1] == g[3] && g[5] == g[7];
  if (!allmod3 && !cheng && !das) return -1;
  memset(&K, 0, sizeof K);
  K.fac = oddfac;

  if (allmod3) {
    /* Iskra 1996 (odd); Cheng/Guo 2019 Theorem 1.1.1 (even). */
    if (_cn_can_order(&K, groups[3])) return 0;
    /* Cheng/Guo 2019 and Reinholz/Spearman/Yang 2013.  Reversing the
     * entire order interchanges positive and negative forward symbols
     * but preserves the gap.  Reinholz's first-position rule is a
     * special case of the odd-gap rule.
     * https://www.sciencedirect.com/science/article/pii/S0022314X18302774
     * https://central.bac-lac.gc.ca/.item?id=TC-BVAU-44941&op=pdf
     */
    return _cn_can_order_one_flip(&K, groups[3]) ? 0 : -1;
  }

  /* Cheng/Guo 2018   http://maths.nju.edu.cn/~guoxj/articles/IJNT2019.pdf */
  if (cheng && _cheng_guo_symbols(&K, noddfac, groups)) return 0;

  /* Das/Saikia 2020  https://math.colgate.edu/~integers/u55/u55.pdf */
  if (das && _das_saikia_symbols(&K, groups)) return 0;

  return -1;
}

/******************************************************************************/

/* Allow testing the filters and the counting functions separately */

int is_congruent_number_filter(UV n) {
  int res;
  factored_t nf;
  factorintp(&nf, n);
  remove_square_part(&nf);
  if (nf.n < 13) return (nf.n >= 5 && nf.n <= 7);

  res = _is_congruent_number_filter1(&nf);
  if (res != -1) return res;
  return _is_congruent_number_filter2(&nf);
}
bool is_congruent_number_tunnell(UV n) {
  factored_t nf;
  factorintp(&nf, n);
  remove_square_part(&nf);
  if (nf.n < 13) return (nf.n >= 5 && nf.n <= 7);

  return _is_congruent_number_tunnell(nf.n);
}

/******************************************************************************/

/* is_congruent_number(n).  OEIS A003273. */
bool is_congruent_number(UV n)
{
  int res;
  factored_t nf;
  factorintp(&nf, n);
  remove_square_part(&nf);
  if (nf.n < 13) return (nf.n >= 5 && nf.n <= 7);

  /* Relatively simple filters.  Order doesn't matter. */
  res = _is_congruent_number_filter1(&nf);
  if (res != -1) return res;
  /* More complicated filters.  Permutation is handled. */
  res = _is_congruent_number_filter2(&nf);
  if (res != -1) return res;

  return _is_congruent_number_tunnell(nf.n);
}
