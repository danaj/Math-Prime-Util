#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "ptypes.h"
#include "factor.h"
#define FUNC_isqrt 1
#include "util.h"
#include "cache.h"

/* Compile with:
 *  gcc -O3 -fomit-frame-pointer -march=native -Wall -DSTANDALONE -I. examples/findomegaseq.c factor.c util.c primality.c cache.c sieve.c chacha.c csprng.c prime_nth_count.c lmo.c -lm
 */

/* A simple example to use the ranged factoring API to find <VSEQ> sequential
 * integers whose omega (count of distinct prime factors) sums to <VSUM>.
 * Testing is from 1 to <VLEN> with results with indices less than <VPRINT>
 * ignored. */

int main(int argc, char *argv[])
{
  UV VSEQ, VSUM, VLEN, VPRINT;

  if (argc != 5) {
    printf("Usage: findseq  SEQLEN  SEQSUM  N  PRINTN\n");
    printf("E.g. findseq 11 23 10000000 10000\n");
    exit(-1);
  }

  VSEQ = strtoul(argv[1], 0, 10);
  VSUM = strtoul(argv[2], 0, 10);
  VLEN = strtoul(argv[3], 0, 10);
  VPRINT = strtoul(argv[4], 0, 10);

  /* ... add input error / range checking here ... */

#if 0
  {
    factor_range_context_t fctx;
    UV V[VSEQ], k, n, nf, T=0, N = VSEQ+VLEN;
    for (k=0; k<VSEQ; k++)  V[k]=0;
    fctx = factor_range_init(1, N, 0);
    for (n = 1; n <= N; n++) {
      int nft = factor_range_next(&fctx);
      /* We have the number of factors!  ... with multiplicity :(. */
      UV *factors = fctx.factors;
      for (nf=1, k=1; k<nft; k++)
        if (factors[k] != factors[k-1])
          nf++;
      /* Now we have the number of distinct factors. */
      T += nf - V[n%VSEQ];
      V[n%VSEQ] = nf;
      if (T == VSUM && n > VPRINT) printf("%"UVuf"\n", n-VSEQ+1);
    }
    factor_range_destroy(&fctx);
  }
#else
  {
    UV V[VSEQ], k, n, nf, T=0, N = VSEQ+VLEN;
    for (k=0; k<VSEQ; k++)  V[k]=0;
#if 0 /* Monolithic.  Very inefficient. */
    unsigned char* nfs = range_omega_sieve(0, N);
    for (n = 1; n <= N; n++) {
      nf = nfs[n];
      T += nf - V[n%VSEQ];
      V[n%VSEQ] = nf;
      if (T == VSUM && n > VPRINT) printf("%"UVuf"\n", n-VSEQ+1);
    }
    Safefree(nfs);
#endif
    prime_precalc(isqrt(N));
    n = 1;
    while (n <= N) {
      UV seg_beg = n;
      UV seg_end = N;
      if ((seg_end - seg_beg + 1) > 200000) seg_end = seg_beg + 200000 - 1;
      unsigned char* nfs = range_nfactor_sieve(seg_beg, seg_end, 0);
      while (n <= seg_end) {
        nf = nfs[n-seg_beg];
        T += nf - V[n%VSEQ];
        V[n%VSEQ] = nf;
        if (T == VSUM && n > VPRINT) printf("%"UVuf"\n", n-VSEQ+1);
        n++;
      }
      Safefree(nfs);
      //n = seg_end+1;
    }
  }
#endif
  return(0);
}
