#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "threadlock.h"
#include "intcomplexity.h"

/* TODO: n/2 space optimization */

static INLINE UV _calcg(UV N) {   /* A000792 */
  UV res = 1;
  if (N >= ((BITS_PER_WORD == 64) ? 122 : 61))  return UV_MAX;
  while (N >= 5 || N == 3) { res *= 3;  N -= 3; }
  return res << (N/2);
}

#if 0  /* Simple no caching */

UV integer_complexity(UV n) { /* A005425 based on Martin N. Fuller's solution */
  uint8_t *A;
  UV G[128];
  UV i, j, ij, res;

  New(0, A, n+1, uint8_t);
  A[0]=0;  A[1]=1;
  for (i = 2; i <= n; i++)  A[i] = 127;
  for (i = 0; i <= 127; i++)   G[i] = _calcg(i);

  for (i = 2; i <= n; i++) {
    UV t = A[i-1], k = 1;
    if (t+1 < A[i])  A[i] = t+1;
    while (k < t/2 && G[k+1]+G[t-k-1] >= i) k++;

    /* Additions */
    for (j = 6; j <= G[k]; j++)
      if (A[j]+A[i-j] < A[i])
        A[i] = A[j] + A[i-j];
    
    /* Multiplications */
    for (j = 2, ij = 2*i;  j <= i && ij <= n;  j++, ij+=i)
      if (A[i]+A[j] < A[ij])
        A[ij] = A[i] + A[j];
  }
  res = A[n];
  Safefree(A);
  return res;
}

#else

static int mutex_init = 0;
READ_WRITE_LOCK_DECL(iccache);
static uint8_t* ic_cache_data = 0;
static UV       ic_cache_size = 0;

/* This runs during a write lock, so modify iccache without concern */
void expand_ic(UV to) {
  uint8_t *A;
  UV G[128];
  UV i, j, ij, first;

  /* Build more than they asked for, but avoid UV wrap. */
  to = (to < UV_MAX - 32767) ? to + 32767 : UV_MAX - 1;

  if (ic_cache_size >= to)
    return;

  Renew(ic_cache_data, to+1, uint8_t);
  A = ic_cache_data;  /* For programmer ease */

  if (ic_cache_size < 2) {
    A[0]=0; A[1]=1;
    first = 2;
  } else {
    first = ic_cache_size+1;
  }

  for (i = first; i <= to; i++)  A[i] = 127;
  for (i = 0; i <= 127; i++)   G[i] = _calcg(i);

  /* Multiplication loop for new data */
  for (i = 2; i < first; i++) {
    for (j = 2, ij = 2*i;  j <= i && ij <= to;  j++, ij+=i)
      if (ij >= first && A[i]+A[j] < A[ij])
        A[ij] = A[i] + A[j];
  }

  /* Fill in rest */
  for (i = first; i <= to; i++) {
    UV t = A[i-1], k = 1;
    if (t+1 < A[i])  A[i] = t+1;
    while (k < t/2 && G[k+1]+G[t-k-1] >= i) k++;

    /* Additions */
    for (j = 6; j <= G[k]; j++)
      if (A[j]+A[i-j] < A[i])
        A[i] = A[j] + A[i-j];
    
    /* Multiplications */
    for (j = 2, ij = 2*i;  j <= i && ij <= to;  j++, ij+=i)
      if (A[i]+A[j] < A[ij])
        A[ij] = A[i] + A[j];
  }

  ic_cache_size = to;
}

UV integer_complexity(UV n) { /* A005425 based on Martin N. Fuller's solution */
  UV res;
  if (!mutex_init) {  /* This needs to be done at boot */
    MUTEX_INIT(&iccache_mutex);
    COND_INIT(&iccache_turn);
    mutex_init = 1;
  }
  if (n == 0) {
    WRITE_LOCK_START(iccache);
    if (ic_cache_data != 0 || ic_cache_size != 0) {
      Safefree(ic_cache_data);
      ic_cache_data = 0;
      ic_cache_size = 0;
    }
    WRITE_LOCK_END(iccache);
    return 0;
  }
  if (n == UV_MAX)
    croak("integer_complexity: input too large");
  READ_LOCK_START(iccache);
  if (ic_cache_size >= n) {
    res = ic_cache_data[n];
    READ_LOCK_END(iccache);
  } else {
    READ_LOCK_END(iccache);
    WRITE_LOCK_START(iccache);
      if (ic_cache_size < n)
        expand_ic(n);
      res = ic_cache_data[n];
    WRITE_LOCK_END(iccache);
  }
  return res;
}
#endif
