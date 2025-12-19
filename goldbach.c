#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "goldbach.h"
#define FUNC_is_prime_in_sieve 1
#include "sieve.h"
#include "cache.h"
#include "primality.h"
#include "util.h"

/* TODO: Consider adding  Waring-Goldbach(n,k,t) */


UV minimal_goldbach_pair(UV n)
{
  UV p;
  if (n < 4) return 0;
  if (n & 1 || n == 4)
    return (is_prime(n-2)) ? 2 : 0;
  /* Maybe this could be faster using a sieve.  Max p < 4*10^18 is 9781 */
  for (p=3; p <= n/2; p = next_prime(p))
    if (is_prime(n-p))
      return p;
  return 0;
}

#if 0
  // Some ways for finding Goldbach pairs

  // 1. Simple
  START_DO_FOR_EACH_PRIME(3,n/2) {
    if (is_prob_prime(n-p))
      L[s++] = p;
  } END_DO_FOR_EACH_PRIME
  Renew(L, s, UV);  *size=s;  return L;

  // 2. Get a full list then walk from the edges.  Not bad for small sizes.
  if (n >= 22 && n < 4294967295U) {
    uint32_t *pr;
    UV nprimes = range_prime_sieve_32(&pr, n, 0);  /* pr[0]=2, pr[1]=3, */
    UV i = 4, j = nprimes-1;
    while (i <= j) {
      UV sum = pr[i] + pr[j];
      if (sum > n) {
        j--;
      } else {
        if (sum == n)
          L[s++] = pr[i];
        i++;
      }
    }
    Safefree(pr);
  }

  // 3. Single sieve
  if (n >= 22) {
    unsigned char* segment;
    UV seg_base, seg_low, seg_high;
    void* ctx = start_segment_primes(n/2, n-11, &segment);
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
        UV q = n-p;
        if (is_prob_prime(q))
          L[s++] = q;
      END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
    sort_uv_array(L, s);
  }

  // The double sieve, low as inner, comes out fastest.
#endif

static UV sieve_pairs(UV* L, UV n) {
  size_t s = 0;
  if (n >= 22) {
    unsigned char* segment;
    UV seg_base, seg_low, seg_high;
    void* ctx = start_segment_primes(n/2, n-11, &segment);
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      size_t qbeg  = n-seg_high,  qend  = n-seg_low;
      UV     qdbeg = qbeg/30,     qdend = (qend+29)/30;
      unsigned char* lowsieve;
      New(0, lowsieve, qdend-qdbeg+1, unsigned char);
      sieve_segment(lowsieve, qdbeg, qdend);
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
        UV q = n-p;
        if (L) {
          if (is_prime_in_sieve(lowsieve, q-qdbeg*30))
            L[s++] = q;
        } else {
          if (is_prime_in_sieve(lowsieve, q-qdbeg*30))
            s++;
        }
      END_DO_FOR_EACH_SIEVE_PRIME
      Safefree(lowsieve);
    }
    end_segment_primes(ctx);
    if (L && s > 1) { /* Reverse the list */
      size_t i = 0, j = s-1;
      while (i < j) { UV t=L[i]; L[i]=L[j]; L[j]=t; i++; j--; }
    }
  }
  return s;
}

UV* goldbach_pairs(size_t *size, UV n)
{
  UV *L;
  size_t s = 0;
  if (n < 4) return 0;
  if (n & 1 || n == 4) {
    if (!is_prime(n-2))
      return 0;
    New(0, L, 1, UV);
    L[0] = 2;
    *size = 1;
    return L;
  }
  /* Overestimate */
  New(0, L, max_nprimes(n/2) >> (n > 30030 ? 1 : 0), UV);

  if (n >=  6 && is_prime(n-3))  L[s++] = 3;
  if (n >= 10 && is_prime(n-5))  L[s++] = 5;
  if (n >= 14 && is_prime(n-7))  L[s++] = 7;

  s += sieve_pairs(L+s, n);

  Renew(L, s, UV);  /* Possibly reduce storage */
  *size = s;
  return L;
}

UV goldbach_pair_count(UV n)
{
  size_t s = 0;
  if (n < 4) return 0;
  if (n & 1 || n == 4)
    return (is_prime(n-2)) ? 1 : 0;

  if (n >=  6 && is_prime(n-3))  s++;
  if (n >= 10 && is_prime(n-5))  s++;
  if (n >= 14 && is_prime(n-7))  s++;

  s += sieve_pairs(0, n);

  return s;
}
