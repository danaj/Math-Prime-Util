#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "constants.h"
#include "lucky_numbers.h"
//#include "util.h"
#include "inverse_interpolate.h"

/******************************************************************************/
/*                             LUCKY NUMBERS                                  */
/******************************************************************************/

static const char _lmask63[63+2] = {1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0,0,1,1,0,1,1,0,0,0,0,1,1,0,1,1,0,1,1,0,0,0,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0,0,0,0,0,1,1};
static const unsigned char _small_lucky[48] = {1,3,7,9,13,15,21,25,31,33,37,43,49,51,63,67,69,73,75,79,87,93,99,105,111,115,127,129,133,135,141,151,159,163,169,171,189,193,195,201,205,211,219,223,231,235,237,241};
static const unsigned char _small_lucky_count[48] = {0,1,1,2,2,2,2,3,3,4,4,4,4,5,5,6,6,6,6,6,6,7,7,7,7,8,8,8,8,8,8,9,9,10,10,10,10,11,11,11,11,11,11,12,12,12,12,12};

/* Lucky Number sieve for 32-bit inputs.
 * Pre-sieve for first 4-7 levels, then in-place deletion using memmove,
 * plus an optimization for a single pass for all single skips.
 * On x86 at least, faster than the other sieves, but uses more memory.
 */
uint32_t* lucky_sieve32(UV *size, UV n) {
  UV i, m, l15, l21, lsize = 0, level, init_level;
  uint32_t *lucky;

  if (n == 0) { *size = 0; return 0; }
  MPUassert(n < 0xFFFFFFFFUL, "lucky_sieve32 n too large");

  /* make initial list, with more culling if it seems worth doing */
  if (n < 300000) {
    const UV fsize = 96*(n+377)/378;  /* 378 = 2*3*7*9 */
    New(0, lucky, 1 + fsize, uint32_t);
    /* Cut 2,3 using mod 6 wheel, 7,9 using a mod 63 mask */
    for (i = 1, m = 1; i <= n; i += 6) {
      if (_lmask63[m  ])  lucky[lsize++] = i;
      if (_lmask63[m+2])  lucky[lsize++] = i+2;
      if ((m += 6) >= 63) m -= 63;
    }
    init_level = 4;
  } else {
    char mask819[819+2];
    const uint16_t v13[16] = {45,96,147,198,252,300,351,402,453,507,555,609,660,708,762,810};
    const UV fsize = (n < 1000000) ? 16128*(n+73709)/73710
                                   : 322560*(n+1547909)/1547910;
    New(0, lucky, 1 + fsize, uint32_t);
    /* Create the mod 819 mask from the smaller one */
    for (i = 0; i < 13; i++) memcpy(mask819+63*i,_lmask63,65);
    for (i = 0; i < 16; i++) mask819[v13[i]] = mask819[v13[i]+1] = 0;
    /* Use the mask and additionally two counters for another two levels */
    for (i = 1, m = 1, l15 = 0, l21 = 0; i <= n; i += 6) {
      if (mask819[m  ] && ++l15 != 15 && ++l21 != 21)  lucky[lsize++] = i;
      if (mask819[m+2] && ++l15 != 15 && ++l21 != 21)  lucky[lsize++] = i+2;
      if ((m += 6) >= 819) m -= 819;
      if (l15 >= 15) l15 -= 15;
      if (l21 >= 21) l21 -= 21;
    }
    init_level = 7;
  }
  if (lucky[lsize-1] > n) lsize--;   /* The +2 could have gone past N */

  /* After the fill-in, we'll start deleting at 13 or 25 */
  for (level = init_level; level < lsize && lucky[level]-1 < lsize; level++) {
    UV skip = lucky[level]-1, nlsize = skip;
    if (2*(skip+1) > lsize) break;  /* Only single skips left */
    for (i = skip+1; i < lsize; i += skip+1) {
      UV ncopy = (skip <= (lsize-i)) ? skip : (lsize-i);
      memmove( lucky + nlsize, lucky + i, ncopy * sizeof(uint32_t) );
      nlsize += ncopy;
    }
    lsize = nlsize;
  }
  /* Now we just have single skips.  Process them all in one pass. */
  if (level < lsize && lucky[level]-1 < lsize) {
    UV skip = lucky[level], nlsize = skip-1;
    while (skip < lsize) {
      UV ncopy = lucky[level+1] - lucky[level];
      if (ncopy > lsize-skip)  ncopy = lsize - skip;
      memmove(lucky + nlsize, lucky + skip, ncopy * sizeof(uint32_t));
      nlsize += ncopy;
      skip += ncopy + 1;
      level++;
    }
    lsize = nlsize;
  }
  *size = lsize;
  return lucky;
}
/* Lucky Number sieve for 64-bit inputs.
 * Uses running counters to skip entries while we add them.
 * Based substantially on Hugo van der Sanden's cgen_lucky.c.
 */
UV* lucky_sieve(UV *size, UV n) {
  UV i, j, c3, lsize, lmax, lindex, *lucky, *count;

  if (n == 0) { *size = 0; return 0; }

  /* Init */
  lmax = (n < 1000) ? 153 : 100 + n/log(n);
  New(0, lucky, lmax, UV);
  New(0, count, lmax, UV);
  lucky[0] = 1;
  lucky[1] = 3;
  lucky[2] = 7;
  lindex = 2;
  lsize = 1;
  c3 = 2;

  for (i = 3; i <= n; i += 2) {
    if (!--c3) { c3 = 3; continue; }  /* Shortcut count[1] */
    for (j = 2; j < lindex; j++) {
      if (--count[j] == 0) {
        count[j] = lucky[j];
        break;
      }
    }
    if (j < lindex) continue;

    if (lsize >= lmax) {  /* Given the estimate, we probably never do this. */
      lmax = 1 + lsize * 1.2;
      Renew(lucky, lmax, UV);
      Renew(count, lmax, UV);
    }
    lucky[lsize] = count[lsize] = i;
    lsize++;

    if (lucky[lindex] == lsize) {
      lindex++;  lsize--;  /* Discard immediately */
    }
  }
  Safefree(count);
  *size = lsize;
  return lucky;
}

/* static UV lucky_count_approx(UV n) { return 0.5 + 0.970 * n / log(n); } */
/* static UV lucky_count_upper(UV n) { return 200 + lucky_count_approx(n) * 1.025; } */

static UV _cb_nlu(UV mid, UV k) { return nth_lucky_upper(mid); }
static UV _cb_nll(UV mid, UV k) { return nth_lucky_lower(mid); }
static UV _cb_nla(UV mid, UV k) { return nth_lucky_approx(mid); }

static UV _simple_lucky_count_approx(UV n) {
  return   (n <        7)  ?  (n > 0) + (n > 2)
         : (n <= 1000000)  ?  0.9957 * n/log(n)
                           : (1.03670 - log(n)/299) * n/log(n);
}
static UV _simple_lucky_count_upper(UV n) {
  return (n <= 10000) ?  10 + _simple_lucky_count_approx(n) * 1.1
                      : 140 + _simple_lucky_count_approx(n) * 1.004;
}
static UV _simple_lucky_count_lower(UV n) {
  return (n <= 10000) ? _simple_lucky_count_approx(n) * 0.9
                      : _simple_lucky_count_approx(n) * 0.99;
}

UV lucky_count_approx(UV n) {
  UV lo, hi;
  if (n < 48) return _small_lucky_count[n];
  /* return _simple_lucky_count_approx(n); */
  lo = _simple_lucky_count_lower(n);
  hi = _simple_lucky_count_upper(n);
  return inverse_interpolate(lo, hi, n, 0, &_cb_nla, 0);
}
UV lucky_count_upper(UV n) {   /* Holds under 1e9 */
  UV lo, hi;
  if (n < 48) return _small_lucky_count[n];
  lo = _simple_lucky_count_lower(n);
  hi = 1.01 * _simple_lucky_count_upper(n);
  return inverse_interpolate(lo, hi, n, 0, &_cb_nll, 0);
}
UV lucky_count_lower(UV n) {   /* Holds under 1e9 */
  UV lo, hi;
  if (n < 48) return _small_lucky_count[n];
  lo = _simple_lucky_count_lower(n);
  hi = _simple_lucky_count_upper(n);
  return inverse_interpolate(lo, hi, n, 0, &_cb_nlu, 0);
}
UV lucky_count(UV n) {
  UV nlucky;
  if (n < 48) {
    nlucky = _small_lucky_count[n];
  } else if (n <= UVCONST(2000000000)) {
    uint32_t *lucky32 = lucky_sieve32(&nlucky, n);
    Safefree(lucky32);
  } else {
    UV *lucky64 = lucky_sieve(&nlucky, n);
    Safefree(lucky64);
  }
  return nlucky;
}

UV nth_lucky_approx(UV n) {
  double corr, fn = n, logn = log(n), loglogn = log(logn);
  if (n <= 48)  return (n == 0) ? 0 : _small_lucky[n-1];
  if (n <= 80000) {
    corr = (n <= 10000)  ?  0.2502  :  0.2581;
    return (UV)(fn*(logn + 0.5*loglogn + corr*loglogn*loglogn) + 0.5);
  } else {
    corr = (n <=    10000) ? -0.0173 :
           (n <=   100000) ? -0.0318 :
           (n <=  1000000) ? -0.0384 :
           (n <= 10000000) ? -0.0422 :
                             -0.0440 ;
    /* Hawkins and Briggs (1958), attributed to S. Chowla. */
    return (UV)(fn*(logn + 0.5*loglogn*loglogn + corr*loglogn*loglogn) + 0.5);
  }
}
UV nth_lucky_upper(UV n) {
  double corr;
  if (n <= 48)  return (n == 0) ? 0 : _small_lucky[n-1];
  corr = (n <= 100) ? 1.05  :
         (n <= 300) ? 1.03  :
         (n <= 800) ? 1.01  :
                      1.0033;   /* verified to 1e9 */
  return (UV)(corr * nth_lucky_approx(n) + 0.5);
}
UV nth_lucky_lower(UV n) {
  double corr;
  if (n <= 48)  return (n == 0) ? 0 : _small_lucky[n-1];
  corr = (n <=  130) ? 0.985 :
         (n <= 2000) ? 0.992 :
                       0.996 ;   /* verified to 1e9 */
  return (UV)(corr * nth_lucky_approx(n));
}

UV nth_lucky(UV n) {
  UV i, k, nlucky;

  if (n <= 48)  return (n == 0) ? 0 : _small_lucky[n-1];

  /* Apply the backward sieve, ala Wilson, for entry n */
  if (n <= UVCONST(2000000000)) {
    uint32_t *lucky32 = lucky_sieve32(&nlucky, n);
    for (i = nlucky-1, k = n-1; i >= 1; i--)
      k += k/(lucky32[i]-1);
    Safefree(lucky32);
  } else {
    UV *lucky64 = lucky_sieve(&nlucky, n);
    for (i = nlucky-1, k = n-1; i >= 1; i--)
      k += k/(lucky64[i]-1);
    Safefree(lucky64);
  }
  return (2 * k + 1);
}



int is_lucky(UV n) {
  uint32_t *lucky32;
  UV i, l, quo, pos, nlucky, lsize;

  /* Simple pre-tests */
  if ( !(n & 1) || (n%6) == 5 || !_lmask63[n % 63]) return 0;
  if (n < 45) return 1;

  /* Check valid position using the static list */
  pos = (n+1) >> 1;  /* Initial position in odds */
  for (i = 1; i < 48; i++) {
    l = _small_lucky[i];
    if (pos < l) return 1;
    quo = pos / l;
    if (pos == quo*l) return 0;
    pos -= quo;
  }

  /* Check more small values */
  if (n >= 1000000U) {
    lucky32 = lucky_sieve32(&nlucky, lsize = lucky_count_upper(n)/25);
    while (i < nlucky) {
      l = lucky32[i++];
      if (pos < l) break;
      quo = pos / l;
      if (pos == quo*l) { Safefree(lucky32); return 0; }
      pos -= quo;
    }
    Safefree(lucky32);
    if (pos < l) return 1;
  }

  /* Generate all needed values and continue checking from where we left off. */
  lucky32 = lucky_sieve32(&nlucky, lsize = lucky_count_upper(n));
  while (1) {
    if (i >= nlucky) { Safefree(lucky32); lucky32 = lucky_sieve32(&nlucky, lsize *= 1.02); }
    l = lucky32[i++];
    if (pos < l)  break;
    quo = pos / l;
    if (pos == quo*l) { Safefree(lucky32); return 0; }
    pos -= quo;
  }
  Safefree(lucky32);
  return 1;
}
