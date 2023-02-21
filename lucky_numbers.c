#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "constants.h"
#include "lucky_numbers.h"
#include "inverse_interpolate.h"
#include "ds_bitmask126.h"

static const int _verbose = 0;

/******************************************************************************/
/*                             LUCKY NUMBERS                                  */
/******************************************************************************/

static const unsigned char _small_lucky[48] = {1,3,7,9,13,15,21,25,31,33,37,43,49,51,63,67,69,73,75,79,87,93,99,105,111,115,127,129,133,135,141,151,159,163,169,171,189,193,195,201,205,211,219,223,231,235,237,241};
static const unsigned char _small_lucky_count[48] = {0,1,1,2,2,2,2,3,3,4,4,4,4,5,5,6,6,6,6,6,6,7,7,7,7,8,8,8,8,8,8,9,9,10,10,10,10,11,11,11,11,11,11,12,12,12,12,12};
/* True for any position where (n % 7*9) could be a lucky number */
static const char _lmask63[63+2] = {1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0,0,1,1,0,1,1,0,0,0,0,1,1,0,1,1,0,1,1,0,0,0,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0,0,0,0,0,1,1};
/* mpufile '$n++; chomp; $v=$_; next unless $v > 10000; $m[ ($v>>1) % 4095 ]++; END { for (0..4094) { next unless $m[$_]; $b[$_ >> 5] |= (1 << ($_%32)); } say join ",",@b; }' ~/misc/ntheory/lucky_1e8.txt */
/* A large bitmask for ((n>>1) % 3*7*3*13) (819).  Covers 2,3,7,9,13. */
static const uint32_t _lmask5[26] = {2334495963,2261929142,1169344621,2204739155,2727961910,1639207725,3513561243,2430232978,1754683725,3630970059,3025873062,1278646881,3658323539,3055177010,1830209833,3406669457,3054200212,1837519692,1531293898,650340770,757258597,2606838995,2530306226,1169218145,3408442969,11572};


/* Lucky Number sieves.
 *
 * Mask presieving for the first 5 levels, followed by pre-sieving with a small
 * number of initial values.
 *
 * For fairly small sieves, e.g. 250k or so values, we use a simplied pagelist.
 * Unlike the full pagelist method, this does not use an index tree.
 *
 * For sieving of non-small sizes, a bitmask (32 bits per 126 integers) is
 * used, with an index tree allowing log(n) time index lookups.  This is much
 * faster and uses substantially less memory than the other methods.  Memory
 * use grows linearly with the sieve size n.
 *
 * Generate first 10M lucky numbers (from 1 to 196502733) on 2020 M1 Mac:
 *           1.8s  bitmask126        memory:  n/25                (  8MB)
 *           3.1s  pagelist_sieve32  memory:  4 * count * ~2.5    (100MB)
 *           4.2s  pagelist_sieve64  memory:  8 * count * ~2.3    (190MB)
 *        1356s    lucky_cgen        memory:  8 * count * 2       (160MB)
 *        8950s    Wilson            memory:  8 * count * 1       ( 80MB)
 *
 * pagelist:
 * nth_lucky(1<<31):   55291335127    47 sec using lucky_sieve32  930MB
 * nth_lucky(1<<32):  113924214621   140 sec using lucky_sieve64  3.2GB
 * nth_lucky(1<<33):  234516370291   312 sec using lucky_sieve64  6.3GB
 * nth_lucky(1<<34):  482339741617   733 sec using lucky_sieve64 12.1GB
 *
 * bitmask:
 * nth_lucky(1<<31):   55291335127    24 sec using lucky_sieve32   89MB
 * nth_lucky(1<<32):  113924214621    51 sec using lucky_sieve64  173MB
 * nth_lucky(1<<33):  234516370291   107 sec using lucky_sieve64  341MB
 * nth_lucky(1<<34):  482339741617   225 sec using lucky_sieve64  677MB
 * nth_lucky(1<<35):  991238156013   471 sec using lucky_sieve64  1.3GB
 * nth_lucky(1<<36): 2035487409679   994 sec using lucky_sieve64  2.6GB
 * nth_lucky(1<<37): 4176793875529  2063 sec using lucky_sieve64  5.3GB
 *
 * A Graviton3 r7g takes about 1.6x more CPU time.
 * nth_lucky(1<<39)  17551419620869 in 264min on Graviton3 r7g, 21GB.
 * nth_lucky(1<<40)  35944896074391 in 582min on Graviton3 r7g, 42GB.
 * nth_lucky(1<<41)  73571139180453 in 1205min on Graviton3 r7g, 84GB.
 * nth_lucky(1<<42) 150499648533909 in 2530min on Graviton3 r7g, 168GB.
 */


/* Simple 32-bit pagelist: fast for small (less than 10M or so) inputs.
 * Simple filtering, then sieve a big block using memmove.
 * This is memory intensive and has poor performance with large n.
 */
uint32_t* _small_lucky_sieve32(UV *size, uint32_t n) {
  uint32_t i, m, c13, level, init_level, fsize, lsize, *lucky;

  if (n < 259) {
    if (n == 0) { *size = 0; return 0; }
    New(0, lucky, 5+n/5, uint32_t);
    for (lsize = 0; lsize < 48 && _small_lucky[lsize] <= n; lsize++)
      lucky[lsize] = _small_lucky[lsize];
    *size = lsize;
    return lucky;
  }

  /* @l=(2,3,7,9,13); $n=vecprod(@l); $n -= divint($n,$_) for @l; say $n */
  fsize = 1152*(n+4913)/4914;
  New(0, lucky, 1 + fsize, uint32_t);
  lsize = c13 = 0;

  /* Create initial list, filtering out 3,7,9,13 */
  for (i = 1, m = 1; i <= n; i += 6) {
    if (_lmask63[m  ]) {
      if (++c13 == 13) c13 = 0;   else lucky[lsize++] = i;
    }
    if (_lmask63[m+2] && (i+2) <= n) {
      if (++c13 == 13) c13 = 0;   else lucky[lsize++] = i+2;
    }
    if ((m += 6) >= 63) m -= 63;
  }
  init_level = 5;

  /* After the fill-in, we'll start deleting at 15 */
  for (level = init_level; level < lsize && lucky[level]-1 < lsize; level++) {
    uint32_t skip = lucky[level]-1, nlsize = skip;
    if (2*(skip+1) > lsize) break;  /* Only single skips left */
    for (i = skip+1; i < lsize; i += skip+1) {
      uint32_t ncopy = (skip <= (lsize-i)) ? skip : (lsize-i);
      memmove( lucky + nlsize, lucky + i, ncopy * sizeof(uint32_t) );
      nlsize += ncopy;
    }
    lsize = nlsize;
  }
  /* Now we just have single skips.  Process them all in one pass. */
  if (level < lsize && lucky[level]-1 < lsize) {
    uint32_t skip = lucky[level], nlsize = skip-1;
    while (skip < lsize) {
      uint32_t ncopy = lucky[level+1] - lucky[level];
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

#if 0  /* No longer used */
#include "ds_pagelist32.h"
uint32_t* _pagelist_lucky_sieve32(UV *size, uint32_t n) {
  uint32_t i, m, lsize, level, init_level, *lucky;
  pagelist32_t *pl;

  if (n > 4294967275U)  n = 4294967275U;  /* Max 32-bit lucky number */

  if (n <= 280000) return _small_lucky_sieve32(size, n);

  pl = pagelist32_create(n);

  /* make initial list using filters for small lucky numbers. */
  {
    UV slsize;
    uint32_t sln, ln, lbeg, lend, *count, *slucky;

    /* Decide how much additional filtering we'll do. */
    sln = (n <= 1000000U)  ?  133  :  (n <= 100000000)  ?  241  :  925;
    slucky = _small_lucky_sieve32(&slsize, sln);
    Newz(0, count, slsize, uint32_t);
    lbeg = 5;
    lend = slsize-1;

    if (1) {
      uint32_t ntarget = (2.4 * (double)n/log(n));
      uint32_t ninit = n/2;
      for (i = 1; i < slsize && ninit > ntarget; i++)
        ninit -= ninit/slucky[i];
      if (i < slsize) lend = i;
      if (lend < lbeg) lend = lbeg;
    }

    if (_verbose) printf("lucky_sieve32 pre-sieve using %u lucky numbers up to %u\n", lend, slucky[lend]);

    /* Construct the initial list */
    for (i = 1, m = 0; i <= n; i += 2, m += 1) {
      if (m >= 819) m -= 819;  /* m = (i>>1) % 819 */
      if (_lmask5[m >> 5] & (1U << (m & 0x1F))) {
        for (ln = lbeg; ln <= lend; ln++) {
          if (++count[ln] == slucky[ln]) {
            count[ln] = 0;
            break;
          }
        }
        if (ln > lend)
          pagelist32_append(pl,i);
      }
    }
    init_level = lend+1;
    Safefree(slucky);
    Safefree(count);
  }

  lsize = pl->nelems;
  if (_verbose) printf("lucky_sieve32 done inserting.  values:  %u   pages: %u\n", lsize, pl->npages[0]);

  if (init_level < lsize) {
    /* Use an iterator rather than calling pagelist32_val(pl,level) */
    pagelist32_iter_t iter = pagelist32_iterator_create(pl, init_level);
    for (level = init_level; level < lsize; level++) {
      uint32_t skip = pagelist32_iterator_next(&iter) - 1;
      if (skip >= lsize) break;
      for (i = skip; i < lsize; i += skip) {
        pagelist32_delete(pl, i);
        lsize--;
      }
    }
    if (_verbose) printf("lucky_sieve32 done sieving.    values:  %u   pages: %u\n", lsize, pl->npages[0]);
  }

  lucky = pagelist32_to_array(size, pl);
  if (*size != lsize) croak("bad sizes in lucky sieve 32");
  if (_verbose) printf("lucky_sieve32 done copying.\n");
  pagelist32_destroy(pl);
  return lucky;
}
#endif

bitmask126_t* bitmask126_sieve(UV* size, UV n) {
  UV i, lsize, level, init_level;
  bitmask126_t *pl;

  pl = bitmask126_create(n);

  /* make initial list using filters for small lucky numbers. */
  {
    UV slsize;
    uint32_t m, sln, ln, lbeg, lend, *count, *slucky;

    /* Decide how much additional filtering we'll do. */
    sln =  (n <=  200000000)  ?   21  :
           (n <= 0xFFFFFFFF)  ?   25  :   87;
    slucky = _small_lucky_sieve32(&slsize, sln);
    Newz(0, count, slsize, uint32_t);
    lbeg = 5;
    lend = slsize-1;

    if (_verbose) printf("bitmask lucky pre-sieve using %u lucky numbers up to %u\n", lend, slucky[lend]);

    /* Construct the initial list */
    for (i = 1, m = 0; i <= n; i += 2, m += 1) {
      if (m >= 819) m -= 819;  /* m = (i>>1) % 819 */
      if (_lmask5[m >> 5] & (1U << (m & 0x1F))) {
        for (ln = lbeg; ln <= lend; ln++) {
          if (++count[ln] == slucky[ln]) {
            count[ln] = 0;
            break;
          }
        }
        if (ln > lend)
          bitmask126_append(pl,i);
      }
    }
    init_level = lend+1;
    Safefree(slucky);
    Safefree(count);
  }

  lsize = pl->nelems;
  if (_verbose) printf("bitmask lucky done inserting.  values:  %lu\n",lsize);

  if (init_level < lsize) {
    bitmask126_iter_t iter = bitmask126_iterator_create(pl, init_level);
    for (level = init_level; level < lsize; level++) {
      UV skip = bitmask126_iterator_next(&iter) - 1;
      if (skip >= lsize) break;
      for (i = skip; i < lsize; i += skip) {
        bitmask126_delete(pl, i);
        lsize--;
      }
    }
    if (_verbose) printf("bitmask lucky done sieving.    values:  %lu\n",lsize);
  }
  *size = lsize;
  return pl;
}

uint32_t* lucky_sieve32(UV *size, uint32_t n) {
  uint32_t *lucky;
  bitmask126_t *pl;

  if (n == 0) { *size = 0; return 0; }
  if (n > 4294967275U)  n = 4294967275U;  /* Max 32-bit lucky number */

  if (n <=   280000U) return _small_lucky_sieve32(size, n);

  pl = bitmask126_sieve(size, n);

  lucky = bitmask126_to_array32(size, pl);
  if (_verbose) printf("lucky_sieve32 done copying.\n");
  bitmask126_destroy(pl);
  return lucky;
}

UV* lucky_sieve64(UV *size, UV n) {
  UV *lucky;
  bitmask126_t *pl;

  if (n == 0) { *size = 0; return 0; }

  pl = bitmask126_sieve(size, n);

  lucky = bitmask126_to_array(size, pl);
  if (_verbose) printf("lucky_sieve64 done copying.\n");
  bitmask126_destroy(pl);
  return lucky;
}


/* Lucky Number sieve for 64-bit inputs.
 * Uses running counters to skip entries while we add them.
 * Based substantially on Hugo van der Sanden's cgen_lucky.c.
 */
UV* lucky_sieve_cgen(UV *size, UV n) {
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

/******************************************************************************/

/* static UV lucky_count_approx(UV n) { return 0.5 + 0.970 * n / log(n); } */
/* static UV lucky_count_upper(UV n) { return 200 + lucky_count_approx(n) * 1.025; } */

static UV _simple_lucky_count_approx(UV n) {
  double logn = log(n);
  return   (n <           7)  ?  (n > 0) + (n > 2)
         : (n <=      10000)  ?  1.03591 * n/logn
         : (n <=    1000000)  ?  0.99575 * n/logn
         : (n <=   10000000)  ? (1.03523 - logn/305) * n/logn
         : (n <=  100000000)  ? (1.03432 - logn/304) * n/logn
         : (n <= 4000000000U) ? (1.03613 - logn/(100*log(logn))) * n/logn
         /* fit 1e9 to 1e10 */
                              : (1.03654 - logn/(100*log(logn))) * n/logn;
}
static UV _simple_lucky_count_upper(UV n) {
  return   (n <= 10000) ?   5 + _simple_lucky_count_approx(n) * 1.01
         : (n <= 1e10)  ?  50 + _simple_lucky_count_approx(n) * 1.011
                        :  60 + _simple_lucky_count_approx(n) * 1.0360;
}
static UV _simple_lucky_count_lower(UV n) {
  return   (n <= 10000) ?  _simple_lucky_count_approx(n) * 0.98
                        :  _simple_lucky_count_approx(n) * 0.99;
}

UV lucky_count_approx(UV n) {
  UV lo, hi;
  if (n < 48) return _small_lucky_count[n];
  /* return _simple_lucky_count_approx(n); */
  lo = _simple_lucky_count_lower(n);
  hi = _simple_lucky_count_upper(n);
  return inverse_interpolate(lo, hi, n, &nth_lucky_approx, 0);
}
UV lucky_count_upper(UV n) {   /* Holds under 1e9 */
  UV lo, hi;
  if (n < 48) return _small_lucky_count[n];
    return _simple_lucky_count_upper(n);
#if 0 && BITS_PER_WORD == 64
  if (n > UVCONST(18409850581000000000))
    return _simple_lucky_count_upper(n);
#endif
  lo = _simple_lucky_count_lower(n);
  hi = _simple_lucky_count_upper(n);
  return inverse_interpolate(lo, hi, n, &nth_lucky_lower, 0);
}
UV lucky_count_lower(UV n) {   /* Holds under 1e9 */
  UV lo, hi;
  if (n < 48) return _small_lucky_count[n];
  lo = _simple_lucky_count_lower(n);
  hi = _simple_lucky_count_upper(n);
  return inverse_interpolate(lo, hi, n, &nth_lucky_upper, 0);
}
UV lucky_count_range(UV lo, UV hi) {
  UV nlucky, lsize;

  if (hi < lo)
    return 0;
  if (hi < 48)
    return _small_lucky_count[hi] - (lo == 0 ? 0 : _small_lucky_count[lo-1]);

  /*
   * Analogous to how nth_lucky works, we sieve enough lucky numbers to
   * ensure we cover everything up to 'hi'.  We can then get an exact
   * count by determining exactly how many values will be removed.
   */

  if ((lo & 1)) lo--;    /* Both lo and hi will be even */
  if ((hi & 1)) hi++;
  lsize = 1+lucky_count_upper(hi);

  if (hi <= UVCONST(2000000000)) {
    uint32_t i, hicount = hi/2, locount = lo/2;
    uint32_t *lucky32 = lucky_sieve32(&nlucky, lsize);
    for (i = 1; i < nlucky && lucky32[i] <= lo; i++) {
      locount -= locount/lucky32[i];
      hicount -= hicount/lucky32[i];
    }
    for ( ; i < nlucky && lucky32[i] <= hicount; i++)
      hicount -= hicount/lucky32[i];
    Safefree(lucky32);
    return hicount - locount;
  } else {
    /* We use the iterator here to cut down on memory use. */
    UV i, hicount = hi/2, locount = lo/2;
    bitmask126_t* pl = bitmask126_sieve(&nlucky, lsize);
    bitmask126_iter_t iter = bitmask126_iterator_create(pl, 1);
    for (i = 1; i < nlucky; i++) {
      UV l = bitmask126_iterator_next(&iter);
      if (l <= lo)  locount -= locount/l;
      if (l > hicount)  break;
      hicount -= hicount/l;
    }
    bitmask126_destroy(pl);
    return hicount - locount;
  }
}
UV lucky_count(UV n) {
  return lucky_count_range(0,n);
}

UV nth_lucky_approx(UV n) {
  double est, corr, fn, logn, loglogn, loglogn2;
  if (n <= 48)  return (n == 0) ? 0 : _small_lucky[n-1];
  fn = n;  logn = log(fn);  loglogn = log(logn);  loglogn2 = loglogn * loglogn;

  /* Use interpolation so we have monotonic growth, as well as good results.
   * We use one formula for small values, and another for larger. */

  /* p1=1<<14; e1=199123;  p2=1<<16; e2=904225;
   * x1=log(log(p1))^2;  x2=log(log(p2))^2;  y1=(e1/p1-log(p1)-0.5*log(log(p1)))/x1;  y2=(e2/p2-log(p2)-0.5*log(log(p2)))/x2;  m=(y2-y1)/(x2-x1);  printf("      corr = %13.11f + %.11f * (loglogn2 - %.11f);\n", y1, m, x1);
   */
  if (n <= 65536) {
    if      (n >= 16384)             /* 16384 -- 65536 */
      corr = 0.25427076035 + 0.00883698771 * (loglogn2 - 5.16445809103);
    else if (n >= 2048)              /*  2048 -- 16384 */
      corr = 0.24513311782 + 0.00880360023 * (loglogn2 - 4.12651426090);
    else if (n >= 256)               /*   256 --   2048 */
      corr = 0.25585213066 - 0.00898952075 * (loglogn2 - 2.93412446098);
    else                             /*    49 --   256 */
      corr = 0.38691439589 - 0.12050840608 * (loglogn2 - 1.84654667704);
    est = fn * (logn + 0.5*loglogn + corr*loglogn2) + 0.5;
  } else {
    /* p1=1<<32; e1=113924214621;   p2=1<<37; e2=4176793875529;
     * x1=log(log(p1))^2;  x2=log(log(p2))^2;  y1=(e1/p1-log(p1)-0.5*x1)/x1;  y2=(e2/p2-log(p2)-0.5*x2)/x2;  m=(y2-y1)/(x2-x1);  printf("      corr = %13.11f + %.11f * (loglogn2 - %.11f);\n", y1, m, x1);
     */
    if    (fn >= 1099511627776.0)    /* 2^40 -- 2^42 */
      corr = -0.05012215934 - 0.00141424235 * (loglogn2 - 11.03811938314);
    else if (fn >= 68719476736.0)    /* 2^36 -- 2^40 */
      corr = -0.04904974983 - 0.00155649126 * (loglogn2 - 10.34912771904);
    else if (fn >= 4294967296.0)     /* 2^32 -- 2^36 */
      corr = -0.04770894029 - 0.00180229750 * (loglogn2 -  9.60518309351);
    else if (fn >=   67108864)       /* 2^26 -- 2^32 */
      corr = -0.04484819198 - 0.00229977135 * (loglogn2 -  8.36125581665);
    else if (fn >=    1048576)       /* 2^20 -- 2^26 */
      corr = -0.03971615189 - 0.00354309756 * (loglogn2 -  6.91279440604);
    else if (n >=      65536)       /* 2^16 -- 2^20 */
      corr = -0.03240114452 - 0.00651036735 * (loglogn2 -  5.78920076332);
    else if (n >=        512)       /* 2^9  -- 2^16 */
      corr =  0.00990254026 - 0.01735396532 * (loglogn2 -  3.35150517018);
    else                             /* 2^6  -- 2^9 */
      corr =  0.13714087150 - 0.09637971899 * (loglogn2 -  2.03132772443);
    /* Hawkins and Briggs (1958), attributed to S. Chowla. */
    est = fn * (logn + (0.5+corr)*loglogn2) + 0.5;
  }
  if (est >= MPU_MAX_LUCKY) return MPU_MAX_LUCKY;
  return (UV)est;
}
UV nth_lucky_upper(UV n) {
  double est, corr;
  if (n <= 48)  return (n == 0) ? 0 : _small_lucky[n-1];
  corr = (n <=  1000) ? 1.01   :
         (n <=  8200) ? 1.005   :
                        1.001;   /* verified to n=3e9 / v=1e11 */
  est = corr * nth_lucky_approx(n) + 0.5;
  if (est >= MPU_MAX_LUCKY) return MPU_MAX_LUCKY;
  return (UV)est;
}
UV nth_lucky_lower(UV n) {
  double est, corr;
  if (n <= 48)  return (n == 0) ? 0 : _small_lucky[n-1];
  est = nth_lucky_approx(n);
  corr = (n <=        122) ? 0.95  :
         (n <=       4096) ? 0.97  :
         (n <=     115000) ? 0.998 :
                             0.999 ;    /* verified to n=3e9 / v=1e11 */
  est = corr * nth_lucky_approx(n);
  return (UV)est;
}

UV nth_lucky(UV n) {
  UV i, k, nlucky;

  if (n <= 48)  return (n == 0) ? 0 : _small_lucky[n-1];

  /* Apply the backward sieve, ala Wilson, for entry n */
  if (n <= UVCONST(100000000)) {
    uint32_t *lucky32 = lucky_sieve32(&nlucky, n);
    for (i = nlucky-1, k = n-1; i >= 1; i--)
      k += k/(lucky32[i]-1);
    Safefree(lucky32);
  } else { /* Iterate backwards through the sieve directly to save memory. */
    bitmask126_t* pl = bitmask126_sieve(&nlucky, n);
    bitmask126_iter_t iter = bitmask126_iterator_create(pl, nlucky-1);
    for (i = nlucky-1, k = n-1; i >= 1; i--)
      k += k / (bitmask126_iterator_prev(&iter) - 1);
    bitmask126_destroy(pl);
  }
  return (2 * k + 1);
}


static int test_lucky_to(UV lsize, UV *beg, UV *end) {
  UV i = *beg, pos = *end, l, quo, nlucky;
  int ret = -1;

  if (lsize <= 700000000U) {
    uint32_t *lucky32 = lucky_sieve32(&nlucky, lsize);
    while (i < nlucky) {
      l = lucky32[i++];
      if (pos < l)      { ret = 1; break; }
      quo = pos / l;
      if (pos == quo*l) { ret = 0; break; }
      pos -= quo;
    }
    Safefree(lucky32);
  } else {
    /* For 64-bit, iterate directly through the bit-mask to save memory. */
    bitmask126_t* pl = bitmask126_sieve(&nlucky, lsize);
    if (i < nlucky) {
      bitmask126_iter_t iter = bitmask126_iterator_create(pl, i);
      while (i < nlucky) {
        l = bitmask126_iterator_next(&iter);
        i++;
        if (pos < l)      { ret = 1; break; }
        quo = pos / l;
        if (pos == quo*l) { ret = 0; break; }
        pos -= quo;
      }
    }
    bitmask126_destroy(pl);
  }
  /* printf("tested lsize = %lu  from %lu to %lu\n", lsize, *beg, i-1); */
  *beg = i;
  *end = pos;
  return ret;
}

int is_lucky(UV n) {
  UV i, l, quo, pos, lsize;
  int res;

  /* Simple pre-tests */
  if ( !(n & 1) || (n%6) == 5 || !_lmask63[n % 63]) return 0;
  if (n < 45) return 1;
  if (n > MPU_MAX_LUCKY) return 0;

  /* Check valid position using the static list */
  pos = (n+1) >> 1;  /* Initial position in odds */

  for (i = 1; i < 48; i++) {
    l = _small_lucky[i];
    if (pos < l) return 1;
    quo = pos / l;
    if (pos == quo*l) return 0;
    pos -= quo;
  }

  lsize = 1+lucky_count_upper(n);

  { /* Check more small values */
    UV psize = 800, gfac = 6;
    while (psize < lsize/3) {
      res = test_lucky_to(psize, &i, &pos);
      if (res != -1) return res;
      psize *= gfac;
      gfac += 2;
    }
  }
  res = test_lucky_to(lsize, &i, &pos);
  return (res == 0) ? 0 : 1;
}
