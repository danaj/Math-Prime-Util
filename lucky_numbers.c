#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#include "constants.h"
#include "lucky_numbers.h"
#include "inverse_interpolate.h"
#include "ds_ull.h"   /* The unrolled linked list that we use for sieving */

static const int _verbose = 0;

/******************************************************************************/
/*                             LUCKY NUMBERS                                  */
/******************************************************************************/

static const unsigned char _small_lucky[48] = {1,3,7,9,13,15,21,25,31,33,37,43,49,51,63,67,69,73,75,79,87,93,99,105,111,115,127,129,133,135,141,151,159,163,169,171,189,193,195,201,205,211,219,223,231,235,237,241};
static const unsigned char _small_lucky_count[48] = {0,1,1,2,2,2,2,3,3,4,4,4,4,5,5,6,6,6,6,6,6,7,7,7,7,8,8,8,8,8,8,9,9,10,10,10,10,11,11,11,11,11,11,12,12,12,12,12};
/* True for any position where (n % 7*9) could be a lucky number */
static const char _lmask63[63+2] = {1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0,0,1,1,0,1,1,0,0,0,0,1,1,0,1,1,0,1,1,0,0,0,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0,0,0,0,0,1,1};
/* A large bitmask for ((n>>1) % 3*7*9*13) (2457) */
static const uint32_t _lmask5[77] = {2334495963,2261929142,1169344621,2204739155,2727961910,1639207725,3513561243,2430232978,1754683725,3630970059,3025873062,1278646881,3658323539,3055177010,1830209833,3406669457,3054200212,1837519692,1531293898,650340770,757258597,2606838995,2530306226,1169218145,3408442969,2799185204,1706318124,1131165330,848965014,699669325,2574062795,2228948386,1284934501,3396110035,2790999218,1161216865,3540886233,2996462004,1771491436,1229828306,881551542,749099117,2590355547,2253614130,1293081133,3408475675,2795071766,1167383305,3540825739,2999594374,1772509509,3378853569,3029544354,1823611745,3664352465,2790870432,1830079336,1529620184,916087988,630608740,453849818,382347702,229021805,2204472531,2257800374,1169312365,3278480475,2195285298,1637143341,3513528987,2967103638,1219909961,3628938443,3025808550,1815517285,3391985233,740658};

/* Lucky Number sieves.
 * Presieving for the first 4-11 levels, then standard sieving.
 *
 * 1) the presieving saves a lot of work but has diminishing returns.
 *
 * 2) using a hybrid array / linked list makes a good tradeoff, especially
 *    with a simple position cache so we don't walk from the start every time.
 *
 * This is not the most memory efficient, but is very fast.
 *
 * Generate first 10M lucky numbers (from 1 to 196502733) on 2020 M1 Mac:
 *          13.8s  lucky_sieve32   memory:  4 * count * ~2.5    (100MB)
 *          27.1s  lucky_sieve64   memory:  8 * count * ~2.3    (190MB)
 *        1356s    lucky_cgen      memory:  8 * count * 2       (160MB)
 *        8950s    Wilson          memory:  8 * count * 1       ( 80MB)
 *
 * nth_lucky(1<<31):   55291335127    25 min using lucky_sieve32, 930MB
 * nth_lucky(1<<32):  113924214621   107 min using lucky_sieve64  3.2GB
 *
 * lucky_sieve64 for the first example uses 398 el/page, second uses 2303.
 */


/* This is fast for small (less than 10M or so) inputs.
 * Simple filtering, then sieve a big block using memmove.
 * This memory intensive and has poor performance with large n.
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

uint32_t* lucky_sieve32(UV *size, uint32_t n) {
  uint32_t i, m, lsize, level, init_level, *lucky;
  ull32_t *pl;

  if (n <= 1000000) return _small_lucky_sieve32(size, n);

  pl = ull32_create();

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
      if (m >= 2457) m -= 2457;  /* m = (i>>1) % 2457 */
      if (_lmask5[m >> 5] & (1U << (m & 0x1F))) {
        for (ln = lbeg; ln <= lend; ln++) {
          if (++count[ln] == slucky[ln]) {
            count[ln] = 0;
            break;
          }
        }
        if (ln > lend)
          ull32_append(pl, i);
      }
    }
    init_level = lend+1;
    Safefree(slucky);
    Safefree(count);
  }

  lsize = pl->nelems;
  if (_verbose) printf("lucky_sieve32 done inserting.  values:  %u   pages: %u\n", lsize, ull32_npages(pl));

  if (init_level < lsize) {
    ull32_iter_t iter = ull32_iterator_create(pl, init_level);
    for (level = init_level; level < lsize; level++) {
      uint32_t skip = ull32_iterator_next(&iter) - 1;
      if (skip >= lsize) break;
      ull32_set_fastpage(pl, skip);
      for (i = skip; i < lsize; i += skip) {
        ull32_delete(pl, i);
        lsize--;
      }
    }
    if (_verbose) printf("lucky_sieve32 done sieving.    values:  %u   pages: %u\n", lsize, ull32_npages(pl));
  }

  lucky = ull32_toarray(size, pl);
  if (*size != lsize) croak("bad sizes in lucky sieve 32");
  if (_verbose) printf("lucky_sieve32 done copying.\n");
  ull32_destroy(pl);
  return lucky;
}

UV* lucky_sieve64(UV *size, UV n) {
  UV i, m, lsize, level, init_level, *lucky;
  ull_t *pl;

  if (n == 0) { *size = 0; return 0; }

  pl = ull_create();

  {
    UV slsize;
    uint32_t sln, ln, lbeg, lend, *count, *slucky;

    /* Decide how much additional filtering we'll do. */
    sln = !(n>>26) ? 391 : !(n>>32) ? 1879 : !(n>>36) ? 5385 : 60000;
    slucky = _small_lucky_sieve32(&slsize, sln);
    Newz(0, count, slsize, uint32_t);
    lbeg = 5;
    lend = slsize-1;

    if (1) {
      UV ntarget = (2.2 * (double)n/log(n));
      UV ninit = n/2;
      for (i = 1; i < slsize && ninit > ntarget; i++)
        ninit -= ninit/slucky[i];
      if (i < slsize) lend = i;
      if (lend < lbeg) lend = lbeg;
    }

    if (_verbose) printf("lucky_sieve64 pre-sieve using %u lucky numbers up to %u\n", lend, slucky[lend]);

    /* Construct the initial list */
    for (i = 1, m = 0; i <= n; i += 2, m += 1) {
      if (m >= 2457) m -= 2457;  /* m = (i>>1) % 2457 */
      if (_lmask5[m >> 5] & (1U << (m & 0x1F))) {
        for (ln = lbeg; ln <= lend; ln++) {
          if (++count[ln] == slucky[ln]) {
            count[ln] = 0;
            break;
          }
        }
        if (ln > lend)
          ull_append(pl, i);
      }
    }
    init_level = lend+1;
    Safefree(slucky);
    Safefree(count);
  }

  lsize = pl->nelems;
  if (_verbose) printf("lucky_sieve64 done inserting.  values:  %lu   pages: %lu\n", lsize, ull_npages(pl));

  if (init_level < lsize) {
    ull_iter_t iter = ull_iterator_create(pl, init_level);
    for (level = init_level; level < lsize; level++) {
      UV skip = ull_iterator_next(&iter) - 1;
      if (skip >= lsize) break;
      ull_set_fastpage(pl, skip);
      for (i = skip; i < lsize; i += skip) {
        ull_delete(pl, i);
        lsize--;
      }
    }
    if (_verbose) printf("lucky_sieve64 done sieving.    values:  %lu   pages: %lu\n", lsize, ull_npages(pl));
  }

  lucky = ull_toarray(size, pl);
  if (*size != lsize) croak("bad sizes in lucky sieve 64");
  if (_verbose) printf("lucky_sieve64 done copying.\n");
  ull_destroy(pl);
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

static UV _cb_nlu(UV mid, UV k) { return nth_lucky_upper(mid); }
static UV _cb_nll(UV mid, UV k) { return nth_lucky_lower(mid); }
static UV _cb_nla(UV mid, UV k) { return nth_lucky_approx(mid); }

static UV _simple_lucky_count_approx(UV n) {
  return   (n <           7)  ?  (n > 0) + (n > 2)
         : (n <=    1000000)  ?  0.9957 * n/log(n)
         : (n <= 1000000000)  ? (1.03670 - log(n)/299) * n/log(n)
                              : (1.03670 - log(n)/(120*log(log(n)))) * n/log(n);
}
static UV _simple_lucky_count_upper(UV n) {
  return (n <= 10000) ?  10 + _simple_lucky_count_approx(n) * 1.1
                      : 140 + _simple_lucky_count_approx(n) * 1.05;
}
static UV _simple_lucky_count_lower(UV n) {
  return (n <= 10000) ? _simple_lucky_count_approx(n) * 0.9
                      : _simple_lucky_count_approx(n) * 0.98;
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
  hi = _simple_lucky_count_upper(n);
  return inverse_interpolate(lo, hi, n, 0, &_cb_nll, 0);
}
UV lucky_count_lower(UV n) {   /* Holds under 1e9 */
  UV lo, hi;
  if (n < 48) return _small_lucky_count[n];
  lo = _simple_lucky_count_lower(n);
  hi = _simple_lucky_count_upper(n);
  return inverse_interpolate(lo, hi, n, 0, &_cb_nlu, 0);
}
UV lucky_count_range(UV lo, UV hi) {
  UV nlo = 0, nlucky, lsize;

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

  if (hi <= UVCONST(4000000000)) {
    uint32_t i, hicount = hi/2, locount = lo/2;
    uint32_t *lucky32 = lucky_sieve32(&nlucky, lsize);
    for (i = 1; i < nlucky && lucky32[i] <= lo; i++) {
      locount -= locount/lucky32[i];
      hicount -= hicount/lucky32[i];
    }
    for ( ; i < nlucky && lucky32[i] <= hi; i++)
      hicount -= hicount/lucky32[i];
    Safefree(lucky32);
    return hicount - locount;
  } else {
    UV i, hicount = hi/2, locount = lo/2;
    UV *lucky64 = lucky_sieve64(&nlucky, lsize);
    for (i = 1; i < nlucky && lucky64[i] <= lo; i++) {
      locount -= locount/lucky64[i];
      hicount -= hicount/lucky64[i];
    }
    for ( ; i < nlucky && lucky64[i] <= hi; i++)
      hicount -= hicount/lucky64[i];
    Safefree(lucky64);
    return hicount - locount;
  }
}
UV lucky_count(UV n) {
  return lucky_count_range(0,n);
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
    return (UV)(fn*(logn + (0.5+corr)*loglogn*loglogn) + 0.5);
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
  if (n <= UVCONST(4000000000)) {
    uint32_t *lucky32 = lucky_sieve32(&nlucky, n);
    for (i = nlucky-1, k = n-1; i >= 1; i--)
      k += k/(lucky32[i]-1);
    Safefree(lucky32);
  } else {
    UV *lucky64 = lucky_sieve64(&nlucky, n);
    for (i = nlucky-1, k = n-1; i >= 1; i--)
      k += k/(lucky64[i]-1);
    Safefree(lucky64);
  }
  return (2 * k + 1);
}


static int test_lucky_to(UV lsize, UV *beg, UV *end) {
  UV i = *beg, pos = *end, l, quo, nlucky;
  int ret = -1;

  if (lsize <= 0xFFFFFFFFU) {
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
    UV* lucky64 = lucky_sieve64(&nlucky, lsize);
    while (i < nlucky) {
      l = lucky64[i++];
      if (pos < l)      { ret = 1; break; }
      quo = pos / l;
      if (pos == quo*l) { ret = 0; break; }
      pos -= quo;
    }
    Safefree(lucky64);
  }
  *beg = i;
  *end = pos;
  return ret;
}

int is_lucky(UV n) {
  UV i, l, quo, pos, nlucky, lsize;
  int res;

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
  if (n > 100000) {
    lsize = n/(1.3*log(n)*log(n));
    if (lsize > 100000) {
      res = test_lucky_to(lsize/100, &i, &pos);
      if (res != -1) return res;
    }
    res = test_lucky_to(lsize, &i, &pos);
    if (res != -1) return res;
  }

  /* Generate all needed values and continue checking from where we left off. */
  lsize = 1+lucky_count_upper(n);
  res = test_lucky_to(lsize, &i, &pos);
  return (res == 0) ? 0 : 1;
}
