#include <stdio.h>
#include <stdlib.h>
#include "csprng.h"
#include "primality.h"
#include "util.h"
#include "prime_counts.h"
#include "mulmod.h"
#include "constants.h"
#include "random_prime.h"

UV random_nbit_prime(void* ctx, UV b)
{
  uint32_t start = 0, range;
  UV n, p;
  switch (b) {
    case 0:
    case 1:  return 0;
    case 2:  return urandomb(ctx,1) ?  2 :  3;
    case 3:  return urandomb(ctx,1) ?  5 :  7;
    case 4:  return urandomb(ctx,1) ? 11 : 13;
    case 5:  start =   7;  range =   5;  break;
    case 6:  start =  12;  range =   7;  break;
    case 7:  start =  19;  range =  13;  break;
    case 8:  start =  32;  range =  23;  break;
    case 9:  start =  55;  range =  43;  break;
    default: break;
  }

  if (start)
    return nth_prime(start + urandomm32(ctx,range));

  if (b > BITS_PER_WORD)
    return 0;

  /* Trivial method */
  p = (UVCONST(1) << (b-1)) + 1;
  while (1) {
    n = p + (urandomb(ctx,b-2) << 1);
    if (is_prob_prime(n))
      return n;
  }
}

UV random_ndigit_prime(void* ctx, UV d)
{
  UV lo, hi;
  if ( (d == 0) || (BITS_PER_WORD == 32 && d >= 10) || (BITS_PER_WORD == 64 && d >= 20) ) return 0;
  if (d == 1) return nth_prime(1 + urandomm32(ctx,4));
  if (d == 2) return nth_prime(5 + urandomm32(ctx,21));
  lo = powmod(10,d-1,UV_MAX)+1;
  hi = 10*lo-11;
  while (1) {
    UV n = (lo + urandomm64(ctx,hi-lo+1)) | 1;
    if (is_prob_prime(n))
      return n;
  }
}

UV random_prime(void* ctx, UV lo, UV hi)
{
  UV n, oddrange;
  if (lo > hi) return 0;
  /* Pull edges in to nearest primes */
  lo = (lo <= 2) ? 2 : next_prime(lo-1);
  hi = (hi >= MPU_MAX_PRIME) ? MPU_MAX_PRIME : prev_prime(hi+1);
  if (lo > hi) return 0;
  /* There must be at least one prime in the range */
  if (!(lo&1)) lo--;             /* treat 2 as 1 */
  oddrange = ((hi-lo)>>1) + 1;   /* look for odds */
  while (1) {
    n = lo + 2 * urandomm64(ctx, oddrange);
    if (n == 1 || is_prob_prime(n))
      return (n == 1) ? 2 : n;
  }
}

/* Note that 7 chosen bases or the first 12 prime bases are enough
 * to guarantee sucess.  We could choose to limit to those. */
int is_mr_random(void* ctx, UV n, UV k) {
  if (k >= 3*(n/4))
    return is_prob_prime(n);

  /* TODO: do 16 at a time */
  while (k--) {
    UV base = 2 + urandomm64(ctx, n-2);
    if (!miller_rabin(n, &base, 1))
      return 0;
  }
  return 1;
}

UV random_semiprime(void* ctx, UV b) { /* Even split of bits */
  static const uint16_t small_semi[] = {35,35,49,65,77,91,143,143,169,299,319,341,377,403};
  UV min, max, n, L, N;

  if (b < 4 || b > BITS_PER_WORD)
    return 0;

  switch (b) {
    case 4:  return 9;
    case 5:  return 21;
    case 6:  return small_semi[  0 + urandomm32(ctx,3) ];
    case 7:  return small_semi[  3 + urandomm32(ctx,3) ];
    case 8:  return small_semi[  6 + urandomm32(ctx,3) ];
    case 9:  return small_semi[  9 + urandomm32(ctx,5) ];
    default: break;
  }

  min = UVCONST(1) << (b-1);
  max = min + (min-1);
  L = b / 2;
  N = b - L;

  do {
    n = random_nbit_prime(ctx,L) * random_nbit_prime(ctx,N);
  } while (n < min || n > max);
  return n;
}

UV random_unrestricted_semiprime(void* ctx, UV b) { /* generic semiprime */
  static const unsigned char small_semi[] = {4,6,9,10,14,15,21,22,25,26,33,34,35,38,39,46,49,51,55,57,58,62,65,69,74,77,82,85,86,87,91,93,94,95,106,111,115,118,119,121,122,123};
  UV min, n;

  if (b < 3 || b > BITS_PER_WORD)
    return 0;

  switch (b) {
    case 3:  return small_semi[  0 + urandomm32(ctx, 2) ];
    case 4:  return small_semi[  2 + urandomm32(ctx, 4) ];
    case 5:  return small_semi[  6 + urandomm32(ctx, 4) ];
    case 6:  return small_semi[ 10 + urandomm32(ctx,12) ];
    case 7:  return small_semi[ 22 + urandomm32(ctx,20) ];
    default: break;
  }
  /* There are faster ways to generate if we could be lax on distribution.
   * Picking a random prime followed by a second that makes a semiprime in
   * the range seems obvious and is fast, but the distribution is wrong.
   * With that method, some semiprimes are much more likely than others. */
  min = UVCONST(1) << (b-1);
  do {
    n = min + urandomb(ctx,b-1);
  } while (!is_semiprime(n));
  return n;
}

UV random_safe_prime(void* ctx, UV bits)
{
  static const unsigned char small_safe[] = {5,7,11,23,47,59,83,107};
  const uint16_t p15mask = 14079;
  UV p, q, B;

  if (bits < 3 || bits > BITS_PER_WORD)
    return 0;

  switch (bits) {
    case 3:  return small_safe[  0 + urandomm32(ctx, 2) ];
    case 4:  return 11;
    case 5:  return 23;
    case 6:  return small_safe[  4 + urandomm32(ctx, 2) ];
    case 7:  return small_safe[  6 + urandomm32(ctx, 2) ];
    default: break;
  }

  /* do { q = 2 * random_nbit_prime(ctx, bits-1) + 1; } while (!is_prob_prime(q)); */

  /* Alternately we could construct p with last 2 bits set, then q = p >> 1. */
  B = (UVCONST(1) << (bits-2)) + 1;
  do {
    p = B + (urandomb(ctx, bits-3) << 1);
    q = 2*p+1;
  } while ( (p15mask & (1U << (p%15))) ||
            (p%7) == 0 || (p%7) == 3 ||
            !is_prob_prime(p) || !is_prob_prime(q) );

  return q;
}
