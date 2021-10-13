#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "constants.h"
#include "prime_powers.h"
#define FUNC_ctz 1
#define FUNC_log2floor 1
#define FUNC_ipow 1
#include "util.h"
#include "cache.h"
#include "sieve.h"
#include "primality.h"
#include "prime_counts.h"
#include "inverse_interpolate.h"

/******************************************************************************/
/*                               PRIME POWERS                                 */
/******************************************************************************/

int prime_power(UV n, UV* prime)
{
  int power = 0;
  if (n < 2) return 0;
  /* Check for small divisors */
  if (!(n&1)) {
    if (n & (n-1)) return 0;
    if (prime) *prime = 2;
    return ctz(n);
  }
  if ((n%3) == 0) {
    /* if (UVCONST(12157665459056928801) % n) return 0; */
    do { n /= 3; power++; } while (n > 1 && (n%3) == 0);
    if (n != 1) return 0;
    if (prime) *prime = 3;
    return power;
  }
  if ((n%5) == 0) {
    do { n /= 5; power++; } while (n > 1 && (n%5) == 0);
    if (n != 1) return 0;
    if (prime) *prime = 5;
    return power;
  }
  if ((n%7) == 0) {
    do { n /= 7; power++; } while (n > 1 && (n%7) == 0);
    if (n != 1) return 0;
    if (prime) *prime = 7;
    return power;
  }
  if (is_prob_prime(n))
    { if (prime) *prime = n; return 1; }
  /* Composite.  Test for perfect power with prime root. */
  power = powerof(n);
  if (power == 1) power = 0;
  if (power) {
    UV root = rootint(n, (UV)power);
    if (is_prob_prime(root))
      { if (prime) *prime = root; }
    else
      power = 0;
  }
  return power;
}


UV next_prime_power(UV n)
{
  UV i, bit;

  if (n < 2) return 2;
  if (n >= MPU_MAX_PRIME) return 0; /* Overflow (max power = max prime) */

#if 0
  /* Straightforward loop */
  for (i = n+1;  !is_prime_power(i);  i++)
    ;
  return i;
#else
  /* Skip evens */
  bit = UVCONST(1) << log2floor(n);
  for (i = n+1+(n&1);  i & bit;  i += 2)
    if (is_prime_power(i))
      return i;
  return i-1;  /* We went past a power of two */
#endif
}

UV prev_prime_power(UV n)
{
  UV i, bit;
  if (n <= 2) return 0;

#if 0
  for (i = n-1;  !is_prime_power(i);  i--)
    ;
  return i;
#else
  n--;
  bit = UVCONST(1) << log2floor(n);
  for (i = n-!(n&1);  i & bit;  i -= 2)
    if (is_prime_power(i))
      return i;
  return i+1;  /* We went past a power of two */
#endif
}

/* The prime powers without the primes */
UV prime_power_sieve2(UV** list, UV lo, UV hi) {
  UV k, log2n, *powers, np = 0, npmax = 0;

  if (hi < 2 || lo > hi) { *list = 0; return 0; }

  /* Estimate how many powers we'll have */
  log2n = log2floor(hi);
  for (k = 2; k <= log2n; k++) {
    npmax += prime_count_upper(rootint(hi,k));
    if (lo > 2) npmax -= prime_count_lower(rootint(lo-1,k));
  }

  New(0, powers, npmax, UV);

  /* Find all powers and add to list */
  for (k = 2; k <= log2n; k++) {
    START_DO_FOR_EACH_PRIME(2, rootint(hi,k)) {
      UV pk = ipow(p,k);
      //if (pk >= lo)  powers[np++] = pk;
      if (pk >= lo)  { if (np >= npmax) croak("overflow"); powers[np++] = pk; }
    } END_DO_FOR_EACH_PRIME
  }

  /* Sort them and return */
  qsort(powers, np, sizeof(UV), _numcmp);
  *list = powers;
  return np;
}

/* The prime powers with the primes */
UV prime_power_sieve(UV** list, UV lo, UV hi) {
  UV npower, nprime, ipower, iprime, ntotal, i, *powers, *primes, *tot;

  if (hi < 2 || lo > hi) { *list = 0; return 0; }

  /* For better performance / memory:
   *   1) realloc primes, use reverse merge to add powers in with one pass
   *   2) sieve the primes here and merge the powers in.
   */

  npower = prime_power_sieve2(&powers, lo, hi);
  nprime = range_prime_sieve(&primes, lo, hi);

  /* The powers get sparse, so this isn't impossible. */
  if (npower == 0) { Safefree(powers); *list = primes; return nprime; }

  ipower = 0;
  iprime = 0;
  ntotal = nprime + npower;
  New(0, tot, ntotal, UV);
  for (i = 0; i < ntotal; i++) {
    if      (ipower == npower)  tot[i] = primes[iprime++];
    else if (iprime == nprime)  tot[i] = powers[ipower++];
    else                        tot[i] = (primes[iprime] < powers[ipower]) ?  primes[iprime++]  :  powers[ipower++];
  }
  Safefree(powers);
  Safefree(primes);
  *list = tot;
  return ntotal;
}


UV prime_power_count_range(UV lo, UV hi) {
  if (hi < 2 || hi < lo) return 0;
  return prime_power_count(hi)  -  ((lo <= 2) ? 0 : prime_power_count(lo-1));
}

/* n A025528; 10^n A267712 */
UV prime_power_count(UV n) {
  uint32_t k, log2n;
  UV sum;

  if (n <= 5) return (n==0) ? 0 : n-1;

  sum = prime_count(n);
  log2n = log2floor(n);
  for (k = 2; k <= log2n; k++)
    sum += prime_count(rootint(n,k));
  return sum;
}

UV prime_power_count_lower(UV n) {
  uint32_t k, log2n;
  UV sum;

  if (n <= 5) return (n==0) ? 0 : n-1;

  sum = prime_count_lower(n);
  log2n = log2floor(n);
  for (k = 2; k <= log2n; k++)
    sum += prime_count_lower(rootint(n,k));
  return sum;
}
UV prime_power_count_upper(UV n) {
  uint32_t k, log2n;
  UV sum;

  if (n <= 5) return (n==0) ? 0 : n-1;

  sum = prime_count_upper(n);
  log2n = log2floor(n);
  for (k = 2; k <= log2n; k++)
    sum += prime_count_upper(rootint(n,k));
  return sum;
}
UV prime_power_count_approx(UV n) {
  uint32_t k, log2n;
  UV sum;

  if (n <= 5) return (n==0) ? 0 : n-1;

  sum = prime_count_approx(n);
  log2n = log2floor(n);
  for (k = 2; k <= log2n; k++)
    sum += prime_count_approx(rootint(n,k));
  return sum;
}

static UV _simple_nth_prime_power_lower(UV n) {
  if (n <= 100) return n+1;
  return (0.98 * nth_prime_lower(n)) - 400;
}
static UV _simple_nth_prime_power_upper(UV n) {
  return nth_prime_upper(n);
}

UV nth_prime_power_lower(UV n) {
  UV lo, hi;
  if (n <= 7) return (n==0) ? 0 : n+1+(n/5);
  lo = _simple_nth_prime_power_lower(n);
  hi = _simple_nth_prime_power_upper(n);
  return inverse_interpolate(lo, hi, n, &prime_power_count_upper, 0);
}
UV nth_prime_power_upper(UV n) {
  UV lo, hi;
  if (n <= 7) return (n==0) ? 0 : n+1+(n/5);
  lo = _simple_nth_prime_power_lower(n);
  hi = _simple_nth_prime_power_upper(n);
  return inverse_interpolate(lo, hi, n, &prime_power_count_lower, 0);
}
UV nth_prime_power_approx(UV n) {
  UV lo, hi;
  if (n <= 7) return (n==0) ? 0 : n+1+(n/5);
  lo = _simple_nth_prime_power_lower(n);
  hi = _simple_nth_prime_power_upper(n);
  return inverse_interpolate(lo, hi, n, &prime_power_count_approx, 0);
}
UV nth_prime_power(UV n) {
  if (n <= 7) return (n==0) ? 0 : n+1+(n/5);
  if (n >= MPU_MAX_PRIME_IDX) return MPU_MAX_PRIME;

#if 0    /* Bilinear interpolation.  Not bad, but not great. */
  UV lo, hi, pp;
  if (n <= 7) return (n==0) ? 0 : n+1+(n/5);

  lo = nth_prime_power_lower(n);
  hi = nth_prime_power_upper(n);
  pp = inverse_interpolate(lo, hi, n, &prime_power_count, 10000);
  return prev_prime_power(pp+1);
#endif

#if 0    /* Approximating interpolation.  Very good, but prefer simpler. */
  UV g, count;
  g = interpolate_with_approx(n, &count, 500,
                              &nth_prime_power_approx, &prime_power_count,
                              0);
  if (g > MPU_MAX_PRIME)
    g = MPU_MAX_PRIME;

  if (count >= n) {
    for (g = prev_prime_power(g+1);  count > n;  count--)
      g = prev_prime_power(g);
  } else {
    for (; count < n; count++)
      g = next_prime_power(g);
  }
  return g;
#endif

  /* Interpolation using functions for approximate nth and exact count.
   * This works quite well, and uses the is_prime_power() function to get
   * the exact result.  Our next/prev functions save negligible time. */
  return interpolate_with_approx(n, 0, 800,
                                 &nth_prime_power_approx, &prime_power_count,
                                 &is_prime_power);
}
