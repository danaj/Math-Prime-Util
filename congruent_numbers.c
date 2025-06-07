#include <stdio.h>

/*
 * A nice discussion of congruent numbers can be found in:
 *
 * https://pub.math.leidenuniv.nl/~stevenhagenp/ANTproc/19yui.pdf
 *
 */

#include "ptypes.h"
#include "congruent_numbers.h"
#define FUNC_isqrt 1
#define FUNC_is_perfect_square 1
#include "util.h"
#include "factor.h"
#include "rootmod.h"
/* This is for randomly permuting the factors.  It will go away */
#include "csprng.h"

/******************************************************************************/

/******************************************************************************/

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

static factored_t permute_odd_factors(factored_t nf, UV seed)
{
  void* rng = prng_new((seed>>16)>>16, seed & 0xFFFFFFFF, 0, 0);
  uint16_t i, j;

  /* for (i=0;i<nf.nfactors i++) printf(" %lu",nf.f[i]);  printf(" -> "); */
  for (i = 1-(nf.n&1); i < nf.nfactors-1; i++) {
    j = (uint16_t) prng_next(rng) % (nf.nfactors-i);
    { UV t      = nf.f[i];  nf.f[i] = nf.f[i+j];  nf.f[i+j] = t; }
    { uint8_t t = nf.e[i];  nf.e[i] = nf.e[i+j];  nf.e[i+j] = t; }
  }
  /* for (i=0;i<nf.nfactors i++) printf(" %lu",nf.f[i]);  printf("\n"); */
  Safefree(rng);
  return nf;
}

/******************************************************************************/


/* Tunnell's method, counting integer solutions to ternary quadratics.
 * Assumes the weak BSD conjecture.
 * Weak BSD holds for n <   42553.   (Nemenzo 1998)
 * Weak BSD holds for n <   71474.   (Wado 2005)
 * Weak BSD holds for n <  300000.   (Matsuno 2005]
 * Weak BSD holds for n < 1000000.   (Matsuno 2006 + Elkies 2002)
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
#define KPR kronecker_uu(p,r)
#define KQR kronecker_uu(q,r)
/* Only used for a few cases: */
#define KQP kronecker_uu(q,p)
#define KRP kronecker_uu(r,p)
#define KRQ kronecker_uu(r,q)
#define LAGRANGE_COND1 (KPQ==KPR || KQR==KQP || KRP==KRQ)
#define LAGRANGE_COND2 ((KPQ==1 && KPR==1) || (KQR==1 && KQP==1) || (KRP==1 && KRQ==1))
#define LAGRANGE_COND3 ((KPQ==-1 && KPR==-1) || (KQR==-1 && KQP==-1) || (KRP==-1 && KRQ==-1))

/* Returns -1 if not known, 0 or 1 indicate definite results. */
int _is_congruent_number_filter1(const factored_t nf) {
  const UV *fac      = nf.f;
  const UV n         = nf.n;
  const int nfactors = nf.nfactors;

  MPUassert(n >= 13, "n too small in icn_filter");

  /* The ACK conjecture (Alter, Curtz, and Kubota 1972):
   *     n = {5,6,7} mod 8   =>  n is a congruent number
   * also follows from the weak BSD conjecture.
   */
  if (n % 8 == 5 || n % 8 == 6 || n % 8 == 7)  return 1;

  /* Evink 2021    https://arxiv.org/pdf/2105.01450.pdf
   * Feng 1996     http://matwbn.icm.edu.pl/ksiazki/aa/aa75/aa7513.pdf
   * Monsky 1990   https://gdz.sub.uni-goettingen.de/id/PPN266833020_0204
   * Lagrange 1974 https://www.numdam.org/item/SDPP_1974-1975__16_1_A11_0.pdf
   */

  if (nfactors == 1) {                   /* n = p */

    UV r, p = n, p8 = p % 8;
    if (p8 == 3) return 0;  /* Genocchi 1855 */

    /* https://arxiv.org/pdf/2105.01450.pdf, Prop 2.1.2 */
    if (sqrtmodp(&r, 2, p) && kronecker_uu(1+r, p) == -1)
      return 0;
#if 0
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

  } else if (!(n&1) && nfactors == 2) {  /* n = 2p */

    UV p = n >> 1, p8 = p % 8;
    if (p8 == 3 || p8 == 7)     return 1;  /* we already returned 1 earlier */
    if (p8 == 5)                return 0;  /* Genocchi 1855 */
    if (p % 16 == 9)            return 0;  /* Bastien 1915 */

  } else if ( (n&1) && nfactors == 2) {  /* n = pq */

    UV p = fac[0], q = fac[1],  p8 = p % 8, q8 = q % 8;
    if (p8 > q8) SWAP4(p,q);
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

  } else if (!(n&1) && nfactors == 3) {  /* n = 2pq */

    UV p = fac[1], q = fac[2],  p8 = p % 8, q8 = q % 8;
    if (p8 > q8) SWAP4(p,q);
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
    if (p8 == 7 && q8 == 7 && KPQ == 1 && q % 16 == 7) return 0;
    if (p8 == 1 && q8 == 1 && KPQ == -1 && (p*q) % 16 == 9) return 0;

  } else if ( (n&1) && nfactors == 3) {  /* n = pqr */

    UV p=fac[0], q=fac[1], r=fac[2], p8=fac[0]%8, q8=fac[1]%8, r8=fac[2]%8;
    if (q8 < p8) SWAP4(p,q);
    if (r8 < q8) SWAP4(q,r);
    if (q8 < p8) SWAP4(p,q);
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

  } else if (!(n&1) && nfactors == 4) {  /* n = 2pqr */

    UV p=fac[1], q=fac[2], r=fac[3], p8=fac[1]%8, q8=fac[2]%8, r8=fac[3]%8;
    if (q8 < p8) SWAP4(p,q);
    if (r8 < q8) SWAP4(q,r);
    if (q8 < p8) SWAP4(p,q);
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
    /* 1 3 7 is not checked.
     * 13706 = 2*7*11*89 = 2*89*7*11, so p = 1, q = -1, r = 3 mod 8.
     *       cond3 (q|r)= (q|p) = -1.
     * but 13706 is a congruent number.  So this seems to be incorrect.
     */
    if (p8 == 3 && q8 == 3 && r8 == 5 && LAGRANGE_COND1) return 0;
    if (p8 == 5 && q8 == 5 && r8 == 5 && LAGRANGE_COND2) return 0;

  } else if ( (n&1) && nfactors == 4) {  /* n = pqrs */

    UV p=fac[0], q=fac[1], r=fac[2], s=fac[3];
    UV p8=fac[0]%8, q8=fac[1]%8, r8=fac[2]%8, s8=fac[3]%8;
    if (r8 < p8) SWAP4(p,r);
    if (s8 < q8) SWAP4(q,s);
    if (q8 < p8) SWAP4(p,q);
    if (s8 < r8) SWAP4(r,s);
    if (r8 < q8) SWAP4(q,r);
    /* Serf 1991 */
    if (p8 == 5 && q8 == 5 && r8 == 7 && s8 == 7 &&
        ( (KPR ==  1 && KQR == -1 && kronecker_uu(p,s) == -1) ||
          (KPR == -1 && kronecker_uu(p,s) ==  1 && kronecker_uu(q,s) == -1) ||
          (KPR == -1 && kronecker_uu(p,s) == -1 && KQR == -kronecker_uu(q,s))))
      return 0;
  }

  return -1;
}

/******************************************************************************/

/* Returns -1 if not known, 0 or 1 indicate definite results. */
int _is_congruent_number_filter2(const factored_t nf) {
  const UV *fac      = nf.f;
  const UV n         = nf.n;
  const int nfactors = nf.nfactors;
  int i, j;

  /* All of these tests seem complicated and tedious.  The goal is to avoid
   * the Tunnell counting loop at the end if at all possible.  While that
   * loop is fast and simple compared to the method of descent, it is still
   * very, very time consuming for any reasonable size input.  So we make
   * an effort to identify known families.
   *
   * These filters all require looking at multiple orderings of the factors.
   * The conditions are true if *any* ordering meets the conditions.
   *
   * TODO: currently none of these look at different orderings.  So they will
   *       answer "I don't know" instead of flagging some numbers as non-cong.
   */

  /* Reinholz 2013 https://central.bac-lac.gc.ca/.item?id=TC-BVAU-44941&op=pdf
   * Cheng 2018   http://maths.nju.edu.cn/~guoxj/articles/IJNT2019.pdf
   * Cheng 2019   https://www.sciencedirect.com/science/article/pii/S0022314X18302774
   * Das 2020     https://math.colgate.edu/~integers/u55/u55.pdf
   */

  {
    const int noddfactors  = (n&1)  ?  nfactors  :  nfactors-1;
    const UV* oddfac       = (n&1)  ?  fac       :  fac+1;
    int k, l, allmod3 = 1;

    for (i = 1;  allmod3 && i <= noddfactors;  i++)
      if ((oddfac[i-1] % 8) != 3)
        allmod3 = 0;

    if (allmod3) { /* Iskra 1996 (odd) ; Cheng,Guo 2019 (even) */
      int iskra = 1;
      for (i = 2; iskra && i <= noddfactors; i++)
        for (j = 1; iskra && j < i; j++)
          iskra &= kronecker_uu(oddfac[j-1],oddfac[i-1]) == -1;
      if (iskra) return 0;
    }

    /* Reinholz, Spearman, Yang 2013 */
    if (allmod3 && (n&1)) {
      int m;
      for (m = 2; m <= nfactors; m += 2) {
        int reinholz = 1;
        for (i = 1; reinholz && i < nfactors; i++)
          for (j = 0; reinholz && j < i; j++)
            if (j == 0 && i == m-1)
              reinholz &= kronecker_uu(fac[j],fac[i]) ==  1;
            else
              reinholz &= kronecker_uu(fac[j],fac[i]) == -1;
        if (reinholz) return 0;
      }
    }

    /* Cheng/Guo 2019 "Some new families of non-congruent numbers" */
    if (allmod3) {
      for (k = 2;  k <= noddfactors;  k++) {
        for (l = 1;  l < k; l++) {
          int cheng = 1;
          if (!((k - l) & 1)) continue;
          for (i = 2; cheng && i <= noddfactors; i++)
            for (j = 1; cheng && j < i; j++)
              if (i == k && j == l)
                cheng &= kronecker_uu(oddfac[j-1],oddfac[i-1]) == -1;
              else
                cheng &= kronecker_uu(oddfac[j-1],oddfac[i-1]) ==  1;
          if (cheng) return 0;
        }
      }
    }

    /* Cheng / Guo 2018 "The non-congruent numbers via Monsky’s formula" */
    if (1) {
      bool quad;
      int g[8] = {0};  /* The number in each mod */
      UV P[MPU_MAX_DFACTORS], Q[MPU_MAX_DFACTORS], R[MPU_MAX_DFACTORS], S[MPU_MAX_DFACTORS];
      const int eps = (n&1) ? 1 : 2;
      for (i = 0; i < noddfactors; i++)  {
        UV m = oddfac[i] % 8;
        if (m == 1) P[ g[m]++ ] = oddfac[i];
        if (m == 3) Q[ g[m]++ ] = oddfac[i];
        if (m == 5) R[ g[m]++ ] = oddfac[i];
        if (m == 7) S[ g[m]++ ] = oddfac[i];
      }
      quad = 1;
      for (i = 2; quad && i <= g[1]; i++)
        for (j = 1; j < i; j++)
          quad &= kronecker_uu(P[j-1],P[i-1]) == -1;
      for (i = 2; quad && i <= g[3]; i++)
        for (j = 1; j < i; j++)
          quad &= kronecker_uu(Q[j-1],Q[i-1]) == -1;
      for (i = 2; quad && i <= g[5]; i++)
        for (j = 1; j < i; j++)
          quad &= kronecker_uu(R[j-1],R[i-1]) == -1;
      for (i = 2; quad && i <= g[7]; i++)
        for (j = 1; j < i; j++)
          quad &= kronecker_uu(S[j-1],S[i-1]) == -1;
      for (i = 1; quad && i <= g[3]; i++)
        for (j = 1; j <= g[1]; j++)
          quad &= kronecker_uu(P[j-1],Q[i-1]) == -1;
      for (i = 1; quad && i <= g[5]; i++)
        for (j = 1; j <= g[1]; j++)
          quad &= kronecker_uu(P[j-1],R[i-1]) == -1;
      for (i = 1; quad && i <= g[7]; i++)
        for (j = 1; j <= g[1]; j++)
          quad &= kronecker_uu(P[j-1],S[i-1]) == -1;
      for (i = 1; quad && i <= g[5]; i++)
        for (j = 1; j <= g[3]; j++)
          quad &= kronecker_uu(Q[j-1],R[i-1]) == -1;
      for (i = 1; quad && i <= g[7]; i++)
        for (j = 1; j <= g[3]; j++)
          quad &= kronecker_uu(Q[j-1],S[i-1]) == -1;
      for (i = 1; quad && i <= g[7]; i++)
        for (j = 1; j <= g[5]; j++)
          quad &= kronecker_uu(R[j-1],S[i-1]) == -1;
      if (quad) {
#if 1  /* Theorem 1.1 */
        if ( (g[1] == 0 && g[5] == 0 && g[7] == 0 && eps == 2 && g[3] % 2 == 0)
          || (g[1] == 0 && g[5] == 0 && g[7] == 0 && eps == 1)
          || (g[1] == 0 && g[3] == 0 && g[7] == 0 && eps == 2 && g[5] % 2 == 0))
          return 0;
#endif
#if 1  /* Theorem 1.2 */
        if ( (g[1] == 0 && g[7] == 0 && eps == 1 && (g[3] % 2) == 1 && g[5] >= 1 && g[5] % 2 == 0)
          || (g[1] == 0 && g[7] == 0 && eps == 2 && g[3] >= 1 && g[5] >= 1 && g[3] % 2 == 0 && g[5] % 2 == 0)
          || (g[1] == 0 && g[7] == 0 && eps == 2 && g[3] >= 1 && g[5] == 1 && g[3] % 2 == 0)
          || (g[1] == 0 && g[5] == 0 && eps == 2 && g[7] == 1 && g[3] % 2 == 1)
          || (g[1] == 0 && g[3] == 0 && eps == 1 && g[5] == 1 && g[7] == 1)
          || (g[5] == 0 && g[7] == 0 && eps == 1 && g[1] > 1 && g[1] % 2 == 0 && g[3] % 2 == 1)
          || (g[5] == 0 && g[7] == 0 && eps == 1 && g[1] == 1 && g[3] % 2 == 1)
          || (g[3] == 0 && g[7] == 0 && eps == 2 && g[1] == 1 && g[5] % 2 == 1)
          || (g[3] == 0 && g[7] == 0 && eps == 2 && g[5] == 1 && g[1] > 1 && g[1] % 2 == 0) )
          return 0;
#endif
#if 1  /* Theorem 1.3 */
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
          return 0;
#endif
#if 1  /* Theorem 1.4 */
        if ( (eps == 1 && g[1] >= 1 && g[7] == 1 && g[1] % 2 == 0 && g[3] % 2 == 1 && g[5] % 2 == 1)
          || (eps == 1 && g[3] >= 1 && g[7] == 1 && g[1] % 2 == 1 && g[5] % 2 == 1 && g[3] % 2 == 0) )
          return 0;
#endif
      }
    }
  }

  /**************************************************************************/

  /* Das / Saikia 2020, extending Lagrange 1974 and Serf 1989 */
  if ((n&1) && nfactors % 2 == 0 && nfactors >= 4 && nfactors <= 20) {
    int cntmod[8] = {0};
    for (i = 0; i < nfactors; i++)  {
      int m = fac[i] % 8;
      cntmod[m]++;
    }
    if (cntmod[1] == cntmod[3] && cntmod[5] == cntmod[7]) {
      /* We can separate all factors into (1,3) and (5,7) pairs. */
      UV pf[10], qf[10];
      int pindexbymod[8], qindexbymod[8];
      const int npairs = nfactors >> 1;
      bool das;

      pindexbymod[1] = qindexbymod[3] = 0;
      pindexbymod[5] = qindexbymod[7] = cntmod[1];
      for (i = 0; i < nfactors; i++) {
        int m = fac[i] % 8;
        if (m == 1 || m == 5) pf[pindexbymod[m]++] = fac[i];
        else                  qf[qindexbymod[m]++] = fac[i];
      }

      /* See if these conditions hold for all pairs */
      das = TRUE;
      for (i = 0; i < npairs; i++)
        das &= kronecker_uu(pf[i],qf[i]) == -1;
      for (i = 0; das && i < npairs; i++) {
        for (j = 0; j < npairs; j++) {
          if (i  > j && kronecker_uu(qf[j],qf[i]) != -1) das = FALSE;
          if (i != j && kronecker_uu(pf[i],pf[j]) !=  1) das = FALSE;
          if (i != j && kronecker_uu(pf[i],qf[j]) !=  1) das = FALSE;
        }
      }
      if (das) return 0;
    }
  }

  return -1;
}

/******************************************************************************/
/* Allow testing the filters and the counting functions separately */

int is_congruent_number_filter(UV n) {
  int res;
  factored_t nf = factorint(n);
  remove_square_part(&nf);
  if (nf.n < 13) return (nf.n >= 5 && nf.n <= 7);

  res = _is_congruent_number_filter1(nf);
  if (res != -1) return res;
  res = _is_congruent_number_filter2(nf);
  return res;
}
bool is_congruent_number_tunnell(UV n) {
  factored_t nf = factorint(n);
  remove_square_part(&nf);
  if (nf.n < 13) return (nf.n >= 5 && nf.n <= 7);

  return _is_congruent_number_tunnell(nf.n);
}

/******************************************************************************/

/* is_congruent_number(n).  OEIS A003273. */
bool is_congruent_number(UV n)
{
  int res;
  factored_t nf = factorint(n);
  remove_square_part(&nf);
  if (nf.n < 13) return (nf.n >= 5 && nf.n <= 7);

  res = _is_congruent_number_filter1(nf);
  if (res != -1) return res;
  res = _is_congruent_number_filter2(nf);
  if (res != -1) return res;

  /* Hacky fun time, randomly permute factors a few times */
  if (0) {
    int i;
    for (i = 1; res == -1 && i <= 10; i++) {
      nf = permute_odd_factors(nf, i);
      res = _is_congruent_number_filter2(nf);
    }
    if (res != -1) return res;
  }

  return _is_congruent_number_tunnell(nf.n);
}
