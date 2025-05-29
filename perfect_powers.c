#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "constants.h"
#include "perfect_powers.h"
#define FUNC_log2floor 1
#define FUNC_ipow 1
#include "util.h"
#include "inverse_interpolate.h"

/******************************************************************************/
/*                             PERFECT POWERS                                 */
/******************************************************************************/

bool is_perfect_power(UV n) {
  return (n == 1 || powerof(n) > 1);
}
bool is_perfect_power_neg(UV n) {
  uint32_t k = powerof(n);
  /* An exponent other than 0,1,2,4,8,16,... is ok */
  return (n == 1 || (k > 2 && (k & (k-1)) != 0));
}


UV next_perfect_power(UV n)
{
  uint32_t k, log2n;
  UV best = MPU_MAX_PERFECT_POW;

  if (n <= 1) return (n) ? 4 : 1;
  if (n >= MPU_MAX_PERFECT_POW) return 0; /* Overflow */

  log2n = log2floor(n);
  for (k = 2; k <= 1+log2n; k++) {
    UV c = ipow( rootint(n,k)+1, k);
    if (c < best && c > n) best = c;
  }
  return best;
}

UV prev_perfect_power(UV n)
{
  uint32_t k, log2n;
  UV best = 4;

  if (n <= 4) return (n > 1);
  if (n > MPU_MAX_PERFECT_POW) return MPU_MAX_PERFECT_POW;

  log2n = log2floor(n);
  for (k = 2; k <= log2n; k++) {
    UV c, r = rootint(n,k);
    c = ipow(r,k);
    if (c >= n) c = ipow(r-1,k);
    if (c > best && c < n) best = c;
  }
  return best;
}

/* Should we have a generator / sieve?   This is a common exercise using PQs. */

UV perfect_power_count_range(UV lo, UV hi) {
  if (hi < 1 || hi < lo) return 0;
  return perfect_power_count(hi) - ((lo <= 1) ? 0 : perfect_power_count(lo-1));
}

static const signed char _moebius[65] = {0,1,-1,-1,0,-1,1,-1,0,0,1,-1,0,-1,1,1,0,-1,0,-1,0,1,1,-1,0,0,1,0,0,-1,-1,-1,0,1,1,1,0,-1,1,1,0,-1,-1,-1,0,0,1,-1,0,0,0,1,0,-1,0,1,0,1,1,-1,0,-1,1,0,0};

/* n A069623; 10^n A070428 */
UV perfect_power_count(UV n) {
  uint32_t k, log2n;
  UV sum;

  if (n < 8) return 0+(n>=1)+(n>=4);

  log2n = log2floor(n);
  for (sum = 1, k = 2; k <= log2n; k++)
    if (_moebius[k])
      sum -= _moebius[k] * (rootint(n, k) - 1);
  return sum;
}

/* About 50 ns per call for exact, so not really worth truncation. */

UV perfect_power_count_lower(UV n) {  return perfect_power_count(n);  }

UV perfect_power_count_upper(UV n) {  return perfect_power_count(n);  }

UV perfect_power_count_approx(UV n) {  return perfect_power_count(n);  }


UV nth_perfect_power_lower(UV n) {
  double pp;
  if (n <= 1) return n;
  if (n >= MPU_MAX_PERFECT_POW_IDX) return MPU_MAX_PERFECT_POW;

  pp = pow(n,2.)  +  (13./3.)*pow(n,4./3.)  +  (32./15.)*pow(n,16./15.);
  pp += -2*pow(n, 5./ 3.) - 2*pow(n, 7./ 5.) - 2*pow(n, 9./ 7.) + 2*pow(n,12./10.);
  pp += -2*pow(n,13./11.) - 2*pow(n,15./13.);
  pp += 5.5;
  if (pp >= UV_MAX) return UV_MAX;
  return (UV)pp;
}
UV nth_perfect_power_upper(UV n) {
  double pp;
  if (n <= 1) return n;
  if (n >= MPU_MAX_PERFECT_POW_IDX) return MPU_MAX_PERFECT_POW;

  pp = pow(n,2.)  +  (13./3.)*pow(n,4./3.)  +  (32./15.)*pow(n,16./15.);
  pp += -2*pow(n, 5./ 3.) - 2*pow(n, 7./ 5.) - 2*pow(n, 9./ 7.) + 2*pow(n,12./10.);
  pp += /* skip 11 and 13 */ 2*pow(n,16./14.);
  pp -= 3.5;
  if (pp >= UV_MAX) return UV_MAX;
  return (UV)pp;
}
UV nth_perfect_power_approx(UV n) {
  double pp;
  if (n <= 1) return n;
  if (n >= MPU_MAX_PERFECT_POW_IDX) return MPU_MAX_PERFECT_POW;

  pp = pow(n,2.)  +  (13./3.)*pow(n,4./3.)  +  (32./15.)*pow(n,16./15.);

#if 0
  uint32_t q;
  for (q = 3; q <= 26; q++) {
    int m = moebius(q);
    if (m == 0 || q == 2 || q == 6 || q == 30) continue;
    pp += m * 2.0 * pow(n, (double)(q+2)/(double)q);
  }
#endif

  pp += -2*pow(n, 5./ 3.) - 2*pow(n, 7./ 5.) - 2*pow(n, 9./ 7.) + 2*pow(n,12./10.);
  pp += -2*pow(n,13./11.) - 2*pow(n,15./13.) + 2*pow(n,16./14.) + 2*pow(n,17./15.);
  pp -= 0.48 * pow(n,19.0/17.0);
  pp -= 1.5;
  if (pp >= UV_MAX) return UV_MAX;
  return (UV)pp;
}

UV nth_perfect_power(UV n) {
  UV g, count;

  if (n <= 1) return n;  /* 1,4,8,9,16,25,... */
  if (n >= MPU_MAX_PERFECT_POW_IDX) return MPU_MAX_PERFECT_POW;

  g = interpolate_with_approx(n, &count, 1000,
                              &nth_perfect_power_approx, &perfect_power_count,
                              0);
  if (g > MPU_MAX_PERFECT_POW)
    g = MPU_MAX_PERFECT_POW;

  if (count >= n) {
    for (g = prev_perfect_power(g+1);  count > n;  count--)
      g = prev_perfect_power(g);
  } else {
    for (; count < n; count++)
      g = next_perfect_power(g);
  }
  return g;
}
