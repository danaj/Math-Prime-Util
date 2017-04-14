#include <stdio.h>
#include <stdlib.h>
#include "csprng.h"
#include "primality.h"
#include "util.h"
#include "prime_nth_count.h"
#include "lmo.h"
#include "mulmod.h"
#include "constants.h"
#include "random_prime.h"

UV random_nbit_prime(UV b)
{
  uint32_t start = 0, range;
  UV n, p;
  switch (b) {
    case 0:
    case 1:  return 0;
    case 2:  return urandomb(1) ?  2 :  3;
    case 3:  return urandomb(1) ?  5 :  7;
    case 4:  return urandomb(1) ? 11 : 13;
    case 5:  start =   7;  range =   5;  break;
    case 6:  start =  12;  range =   7;  break;
    case 7:  start =  19;  range =  13;  break;
    case 8:  start =  32;  range =  23;  break;
    case 9:  start =  55;  range =  43;  break;
    default: break;
  }

  if (start)
    return nth_prime(start + urandomm32(range));

  if (b > BITS_PER_WORD)
    return 0;

  /* Trivial method */
  p = (UVCONST(1) << (b-1)) + 1;
  while (1) {
    n = p + (urandomb(b-2) << 1);
    if (is_prob_prime(n))
      return n;
  }
}

UV random_ndigit_prime(UV d)
{
  UV lo, hi;
  if ( (d == 0) || (BITS_PER_WORD == 32 && d >= 10) || (BITS_PER_WORD == 64 && d >= 20) ) return 0;
  if (d == 1) return nth_prime(1 + urandomm32(4));
  if (d == 2) return nth_prime(5 + urandomm32(21));
  lo = powmod(10,d-1,UV_MAX)+1;
  hi = 10*lo-11;
  while (1) {
    UV n = (lo + urandomm64(hi-lo+1)) | 1;
    if (is_prob_prime(n))
      return n;
  }
}

UV random_prime(UV lo, UV hi)
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
    n = lo + 2 * urandomm64(oddrange);
    if (n == 1 || is_prob_prime(n))
      return (n == 1) ? 2 : n;
  }
}

/* Note that 7 chosen bases or the first 12 prime bases are enough
 * to guarantee sucess.  We could choose to limit to those. */
int is_mr_random(UV n, UV k) {
  if (k >= 3*(n/4))
    return is_prob_prime(n);

  /* TODO: do 16 at a time */
  while (k--) {
    UV base = 2 + urandomm64(n-2);
    if (!miller_rabin(n, &base, 1))
      return 0;
  }
  return 1;
}
