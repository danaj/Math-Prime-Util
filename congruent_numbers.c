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

#define SWAP4(x,y) { UV t; t=x; x=y; y=t;  t=x##8; x##8=y##8; y##8=t; }

/* is_congruent_number(n).  OEIS A003273. */
int is_congruent_number(UV n) {
  UV N, fac[MPU_MAX_FACTORS+1], exp[MPU_MAX_FACTORS+1];
  UV m8 = n % 8;
  int i, j, nfactors;

  if (n < 13)   return (n >= 5 && n <= 7);
  /* ACK conjecture, n in CN if n = {5,6,7} mod 8.  (BSD implies ACK) */
  if (m8 == 5 || m8 == 6 || m8 == 7)  return 1;

  nfactors = factor_exp(n, fac, exp);
  for (i = 0, j = 0, N = 1; i < nfactors; i++) {
    if (exp[i] & 1) {
      fac[j++] = fac[i];
      N *= fac[i];
    }
  }
  if (n != N) {    /* n wasn't square free, replace with N. */
    nfactors = j;
    n = N;
    m8 = n % 8;
    if (n < 13)   return (n >= 5 && n <= 7);
    if (m8 == 5 || m8 == 6 || m8 == 7)  return 1;
  }
  /* n is now square free and we can ignore exp[] */

  /*
   * Evink 2021   https://arxiv.org/pdf/2105.01450.pdf
   * Feng 1996    http://matwbn.icm.edu.pl/ksiazki/aa/aa75/aa7513.pdf
   * Monsky 1990  https://gdz.sub.uni-goettingen.de/id/PPN266833020_0204
   * Reinholz 2013 https://central.bac-lac.gc.ca/.item?id=TC-BVAU-44941&op=pdf
   * Das 2020     https://math.colgate.edu/~integers/u55/u55.pdf
   * Cheng 2018   http://maths.nju.edu.cn/~guoxj/articles/IJNT2019.pdf
   * Cheng 2019   https://www.sciencedirect.com/science/article/pii/S0022314X18302774
   */

  if (nfactors == 1) {                   /* n = p */

    UV r, p = n;
    if (m8 == 3)  return 0;
    if (m8 == 5 || m8 == 7)  return 1;  /* not seen here, handled earlier */

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

    UV p = n >> 1, m8 = p % 8;
    if (m8 == 3 || m8 == 7)     return 1;  /* we don't see these here */
    if (m8 == 5 || (p%16) == 9) return 0;  /* Bastien 1915 */

  } else if ( (n&1) && nfactors == 2) {  /* n = pq */

    UV p = fac[0], q = fac[1],  p8 = p % 8, q8 = q % 8;
    if (p8 > q8) SWAP4(p,q);
    if (p8 == 3 && q8 == 3) return 0;
#if 0  /* Monsky, all special cases of ACK */
    if (p8 == 3 && q8 == 7) return 1;
    if (p8 == 3 && q8 == 5) return 1;
    if (p8 == 1 && q8 == 5 && kronecker_uu(p,q) == -1) return 1;
    if (p8 == 1 && q8 == 7 && kronecker_uu(p,q) == -1) return 1;
#endif
    if (p8 == 1 && q8 == 3 && kronecker_uu(p,q) == -1) return 0;
    if (p8 == 5 && q8 == 7 && kronecker_uu(p,q) == -1) return 0;

  } else if (!(n&1) && nfactors == 3) {  /* n = 2pq */

    UV p = fac[1], q = fac[2],  p8 = p % 8, q8 = q % 8;
    if (p8 > q8) SWAP4(p,q);
    if (p8 == 5 && q8 == 5) return 0;
    if (p8 == 3 && q8 == 3) return 0; /* Lagrange 1974 */
#if 0  /* Monsky, all special cases of ACK */
    if (p8 == 3 && q8 == 5) return 1;
    if (p8 == 5 && q8 == 7) return 1;
    if (p8 == 1 && q8 == 7 && kronecker_uu(p,q) == -1) return 1;
    if (p8 == 1 && q8 == 3 && kronecker_uu(p,q) == -1) return 1;
#endif
    if (p8 == 1 && q8 == 5 && kronecker_uu(p,q) == -1) return 0;
    if (p8 == 3 && q8 == 7 && kronecker_uu(p,q) == -1) return 0;

  } else if ( (n&1) && nfactors == 3) {  /* n = pqr */
    UV p=fac[0], q=fac[1], r=fac[2], p8=fac[0]%8, q8=fac[1]%8, r8=fac[2]%8;
    if (p8 == 3 && q8 == 1) SWAP4(p,q);
    if (q8 == 1 && r8 == 3) SWAP4(q,r);
    if (p8 == 1 && q8 == 3 && r8 == 1 && kronecker_uu(p,q) == -1 && kronecker_uu(q,r) == -1) return 0;

    if (q8 < p8) SWAP4(p,q);
    if (r8 < q8) SWAP4(q,r);
    if (q8 < p8) SWAP4(p,q);
    /* Serf 1991 */
    if (p8 == 3 && q8 == 3 && r8 == 5) return 1;
    if (p8 == 3 && q8 == 3 && r8 == 7) return 1;
    if (p8 == 7 && q8 == 7 && r8 == 7 && kronecker_uu(p,q) == -kronecker_uu(p,r) == kronecker_uu(q,r)) return 1;
    if (p8 == 1 && q8 == 3 && r8 == 3 && kronecker_uu(p,q) == -kronecker_uu(p,r)) return 0;
    if (p8 == 3 && q8 == 5 && r8 == 7 && kronecker_uu(q,r) == -1) return 0;
    if (p8 == 3 && q8 == 7 && r8 == 7 && kronecker_uu(p,q) == -kronecker_uu(p,r) && kronecker_uu(p,q) == kronecker_uu(q,r)) return 0;

  } else if (!(n&1) && nfactors == 4) {  /* n = 2pqr */
    UV p=fac[1], q=fac[2], r=fac[3], p8=fac[1]%8, q8=fac[2]%8, r8=fac[3]%8;
    if (p8 == 5 && q8 == 1) SWAP4(p,q);
    if (q8 == 1 && r8 == 5) SWAP4(q,r);
    if (p8 == 1 && q8 == 5 && r8 == 1 && kronecker_uu(p,q) == -1 && kronecker_uu(q,r) == -1) return 0;

    if (q8 < p8) SWAP4(p,q);
    if (r8 < q8) SWAP4(q,r);
    if (q8 < p8) SWAP4(p,q);
    /* Serf 1991 */
    if (p8 == 3 && q8 == 3 && r8 == 7) return 1;
    if (p8 == 3 && q8 == 5 && r8 == 5) return 1;
    if (p8 == 5 && q8 == 5 && r8 == 7) return 1;
    if (p8 == 7 && q8 == 7 && r8 == 7 && kronecker_uu(p,q) == -kronecker_uu(p,r) && kronecker_uu(p,q) == kronecker_uu(q,r)) return 1;
    /* Lagrange 1974 */
    if (p8 == 1 && q8 == 3 && r8 == 3 && kronecker_uu(p,q) == -kronecker_uu(p,r)) return 0;
    if (p8 == 1 && q8 == 5 && r8 == 5 && kronecker_uu(p,q) == -kronecker_uu(p,r)) return 0;
    if (p8 == 3 && q8 == 5 && r8 == 7 && kronecker_uu(p,r) == -kronecker_uu(q,r)) return 0;
    if (p8 == 5 && q8 == 7 && r8 == 7 && kronecker_uu(p,q) == -kronecker_uu(p,r) && kronecker_uu(p,q) == kronecker_uu(q,r)) return 0;

  }

  {
    int noddfactors = (n&1)  ?  nfactors  :  nfactors-1;
    UV* oddfac      = (n&1)  ?  fac       :  fac+1;
    int i, j, k, l, allmod3 = 1;

    for (i = 1;  allmod3 && i <= noddfactors;  i++)
      if ((oddfac[i-1] % 8) != 3)
        allmod3 = 0;

    if (allmod3) { /* Iskra 1996 (odd) ; Cheng,Guo 2019 (even) */
      int iskra = 1;
      for (i = 2; iskra && i <= noddfactors; i++)
        for (j = 1; iskra && j < i; j++)
          if (kronecker_uu(oddfac[j-1],oddfac[i-1]) != -1)
            iskra = 0;
      if (iskra) return 0;
    }

    /* Reinholz, Spearman, Yang 2013 */
    if (allmod3 && (n&1)) {
      int m;
      for (m = 2; m <= nfactors; m += 2) {
        int reinholz = 1;
        for (i = 1; reinholz && i < nfactors; i++)
          for (j = 0; reinholz && j < i; j++)
            if (j == 0 && i == m-1) {
              if (kronecker_uu(fac[j],fac[i]) !=  1)  reinholz = 0;
            } else {
              if (kronecker_uu(fac[j],fac[i]) != -1)  reinholz = 0;
            }
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
              if (i == k && j == l) {
                if (kronecker_uu(oddfac[j-1],oddfac[i-1]) != -1)  cheng = 0;
              } else {
                if (kronecker_uu(oddfac[j-1],oddfac[i-1]) !=  1)  cheng = 0;
              }
          if (cheng) return 0;
        }
      }
    }

    /* Cheng / Guo 2018 "The non-congruent numbers via Monsky’s formula" */
    if (1) {
      int quad = 1;
      int mod8[MPU_MAX_FACTORS+1];
      int g[8] = {0};  /* The number in each mod */
      UV P[MPU_MAX_FACTORS+1], Q[MPU_MAX_FACTORS+1], R[MPU_MAX_FACTORS+1], S[MPU_MAX_FACTORS+1];
      int eps = (n&1) ? 1 : 2;
      for (i = 0; i < noddfactors; i++)  {
        UV m = oddfac[i] % 8;
        mod8[i] = m;
        if (m == 1) P[ g[m]++ ] = oddfac[i];
        if (m == 3) Q[ g[m]++ ] = oddfac[i];
        if (m == 5) R[ g[m]++ ] = oddfac[i];
        if (m == 7) S[ g[m]++ ] = oddfac[i];
      }
      for (i = 2; quad && i <= g[1]; i++)
        for (j = 1; j < i; j++)
          if (kronecker_uu(P[j-1],P[i-1]) != -1)
            quad = 0;
      for (i = 2; quad && i <= g[3]; i++)
        for (j = 1; j < i; j++)
          if (kronecker_uu(Q[j-1],Q[i-1]) != -1)
            quad = 0;
      for (i = 2; quad && i <= g[5]; i++)
        for (j = 1; j < i; j++)
          if (kronecker_uu(R[j-1],R[i-1]) != -1)
            quad = 0;
      for (i = 2; quad && i <= g[7]; i++)
        for (j = 1; j < i; j++)
          if (kronecker_uu(S[j-1],S[i-1]) != -1)
            quad = 0;
      for (i = 1; quad && i <= g[3]; i++)
        for (j = 1; j <= g[1]; j++)
          if (kronecker_uu(P[j-1],Q[i-1]) != -1)
            quad = 0;
      for (i = 1; quad && i <= g[5]; i++)
        for (j = 1; j <= g[1]; j++)
          if (kronecker_uu(P[j-1],R[i-1]) != -1)
            quad = 0;
      for (i = 1; quad && i <= g[7]; i++)
        for (j = 1; j <= g[1]; j++)
          if (kronecker_uu(P[j-1],S[i-1]) != -1)
            quad = 0;
      for (i = 1; quad && i <= g[5]; i++)
        for (j = 1; j <= g[3]; j++)
          if (kronecker_uu(Q[j-1],R[i-1]) != -1)
            quad = 0;
      for (i = 1; quad && i <= g[7]; i++)
        for (j = 1; j <= g[3]; j++)
          if (kronecker_uu(Q[j-1],S[i-1]) != -1)
            quad = 0;
      for (i = 1; quad && i <= g[7]; i++)
        for (j = 1; j <= g[5]; j++)
          if (kronecker_uu(R[j-1],S[i-1]) != -1)
            quad = 0;
      if (quad) {
        int cheng = 0;
#if 1  /* Theorem 1.1 */
        if (g[1] == 0 && g[5] == 0 && g[7] == 0 && eps == 2 && (g[3] % 2) == 0) cheng = 1;
        if (g[1] == 0 && g[5] == 0 && g[7] == 0 && eps == 1) cheng = 1;
        if (g[1] == 0 && g[3] == 0 && g[7] == 0 && eps == 2 && (g[5] % 2) == 0) cheng = 1;
#endif
#if 1  /* Theorem 1.2 */
        if (g[1] == 0 && g[7] == 0 && eps == 1 && (g[3] % 2) == 1 && g[5] >= 1 && (g[5] % 2) == 0) cheng = 1;
        if (g[1] == 0 && g[7] == 0 && eps == 2 && g[3] >= 1 && g[5] >= 1 && (g[3] % 2) == 0 && (g[5] % 2) == 0) cheng = 1;
        if (g[1] == 0 && g[7] == 0 && eps == 2 && g[3] >= 1 && g[5] == 1 && (g[3] % 2) == 0) cheng = 1;
        if (g[1] == 0 && g[5] == 0 && eps == 2 && g[7] == 1 && (g[3] % 2) == 1) cheng = 1;
        if (g[1] == 0 && g[3] == 0 && eps == 1 && g[5] == 1 && g[7] == 1) cheng = 1;
        if (g[5] == 0 && g[7] == 0 && eps == 1 && g[1] > 1 && (g[1] % 2) == 0 && (g[3] % 2) == 1) cheng = 1;
        if (g[5] == 0 && g[7] == 0 && eps == 1 && g[1] == 1 && (g[3] % 2) == 1) cheng = 1;
        if (g[3] == 0 && g[7] == 0 && eps == 2 && g[1] == 1 && (g[5] % 2) == 1) cheng = 1;
        if (g[3] == 0 && g[7] == 0 && eps == 2 && g[5] == 1 && g[1] > 1 && (g[1] % 2) == 0) cheng = 1;
#endif
#if 1  /* Theorem 1.3 */
        if (g[1] == 0 && eps == 1 && g[5] == 1 && g[7] == 1 && g[3] >= 1) cheng = 1;
        if (g[1] == 0 && eps == 1 && g[5] >= 2 && g[7] == 1 && (g[3] % 2) == 1 && (g[5] % 2) == 1) cheng = 1;
        if (g[1] == 0 && eps == 2 && g[5] >= 1 && g[7] == 1 && (g[3] % 2) == 1 && (g[5] % 2) == 0) cheng = 1;
        if (g[3] == 0 && eps == 1 && g[7] == 1 && (g[1] % 2) == 1 && (g[5] % 2) == 1) cheng = 1;
        /* No examples 1.3.4 */
        if (g[3] == 0 && eps == 2 && g[5] == 1 && g[7] == 1 && (g[1] % 2) == 1) cheng = 1;
        /* No examples 1.3.5 */
        if (g[5] == 0 && eps == 1 && g[1] >= 1 && g[3] >= 1 && g[7] == 1 && (g[1] % 2) == 0 && (g[3] % 2) == 0) cheng = 1;
        if (g[5] == 0 && eps == 2 && g[1] >= 1 && g[3] >= 1 && g[7] == 1 && (g[1] % 2) == 1 && (g[3] % 2) == 1) cheng = 1;
        if (g[7] == 0 && eps == 1 && g[1] >= 1 && g[5] >= 1 && (g[1] % 2) == 0 && (g[5] % 2) == 0 && (g[3] % 2) == 1) cheng = 1;
        if (g[7] == 0 && eps == 2 && g[1] == 1 && g[3] >= 1 && (g[3] % 2) == 0 && (g[5] % 2) == 1) cheng = 1;
        if (g[7] == 0 && eps == 2 && g[1] >= 1 && g[3] >= 1 && g[5] == 1 && (g[1] % 2) == 0 && (g[3] % 2) == 0) cheng = 1;
#endif
#if 1  /* Theorem 1.4 */
        if (eps == 1 && g[1] >= 1 && g[7] == 1 && (g[1] % 2) == 0 && (g[3] % 2) == 1 && (g[5] % 2) == 1) cheng = 1;
        if (eps == 1 && g[3] >= 1 && g[7] == 1 && (g[1] % 2) == 1 && (g[5] % 2) == 1 && (g[3] % 2) == 0) cheng = 1;
#endif
        if (cheng) return 0;
      }
    }
  }

  /* Das / Saikia 2020, extending Lagrange 1974 and Serf 1989 */
  if ((n&1) && (nfactors == 4 || nfactors == 6 || nfactors == 8)) {
    int mod8[10], npairs, epairs, das;
    UV pf[5],qf[5], totmod[8] = {0};
    for (i = 0; i < nfactors; i++)  {
      UV m = fac[i] % 8;
      mod8[i] = m;
      totmod[m]++;
    }
    if (totmod[1] == totmod[3] && totmod[5] == totmod[7]) {
      epairs = nfactors >> 1;
      for (npairs = 0; npairs <epairs; npairs++) {
        for (i = 0; i < nfactors; i++)
          if (mod8[i] == 1)
            break;
        if (i < nfactors) {
          for (j = 0; j < nfactors; j++)
            if (mod8[j] == 3)
              break;
          if (j >= nfactors) break;
          mod8[i] = mod8[j] = 0;
          pf[npairs] = fac[i];
          qf[npairs] = fac[j];
          continue;
        }
        for (i = 0; i < nfactors; i++)
          if (mod8[i] == 5)
            break;
        if (i < nfactors) {
          for (j = 0; j < nfactors; j++)
            if (mod8[j] == 7)
              break;
          if (j >= nfactors) break;
          mod8[i] = mod8[j] = 0;
          pf[npairs] = fac[i];
          qf[npairs] = fac[j];
          continue;
        }
        break;
      }
      if (npairs == epairs) {
        das = 1;
        for (i = 0; i < npairs; i++) {
          for (j = 0; j < npairs; j++) {
            if (i  > j && kronecker_uu(qf[j],qf[i]) != -1) das = 0;
            if (i != j && kronecker_uu(pf[i],pf[j]) !=  1) das = 0;
            if (i != j && kronecker_uu(pf[i],qf[j]) !=  1) das = 0;
            if (i == j && kronecker_uu(pf[i],qf[j]) != -1) das = 0;
          }
        }
        if (das) return 0;
      }
    }
  }

  /* Tunnell's method, counting integer solutions to ternary quadratics. */
  /* Assumes the weak BSD conjecture. */
  {
    UV x, y, z, limz, limy, limx, n8z, zsols, sols[2] = {0,0};

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

  return -1;
}
