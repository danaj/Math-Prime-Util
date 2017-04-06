#include <stdio.h>
#include <stdlib.h>
#include "isaac.h"
#include "primality.h"
#include "util.h"
#include "lmo.h"
#include "mulmod.h"
#include "random_prime.h"

UV random_nbit_prime(UV b)
{
  uint32_t start = 0, range;
  UV n, p;
  switch (b) {
    case 0:
    case 1:  return 0;
    case 2:  return irandb(1) ?  2 :  3;
    case 3:  return irandb(1) ?  5 :  7;
    case 4:  return irandb(1) ? 11 : 13;
    case 5:  start =   7;  range =   5;  break;
    case 6:  start =  12;  range =   7;  break;
    case 7:  start =  19;  range =  13;  break;
    case 8:  start =  32;  range =  23;  break;
    case 9:  start =  55;  range =  43;  break;
    default: break;
  }

  if (start)
    return nth_prime(start + isaac_rand(range));

  if (b > BITS_PER_WORD)
    return 0;

  /* Trivial method */
  p = (UVCONST(1) << (b-1)) + 1;
  while (1) {
    n = p + (irandb(b-2) << 1);
    if (is_prob_prime(n))
      return n;
  }
}

UV random_ndigit_prime(UV d)
{
  UV lo, hi;
  if ( (d == 0) || (BITS_PER_WORD == 32 && d >= 10) || (BITS_PER_WORD == 64 && d >= 20) ) return 0;
  if (d == 1) return nth_prime(1 + isaac_rand(4));
  if (d == 2) return nth_prime(5 + isaac_rand(21));
  lo = powmod(10,d-1,UV_MAX)+1;
  hi = 10*lo+9;
  while (1) {
    UV n = (lo + isaac_rand64(hi-lo+1)) | 1;
    if (is_prob_prime(n))
      return n;
  }
}
