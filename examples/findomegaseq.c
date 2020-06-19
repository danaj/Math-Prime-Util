#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "ptypes.h"
#include "factor.h"

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
  return(0);
}
