#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ptypes.h"
#define FUNC_isqrt 1
#define FUNC_icbrt 1
#define FUNC_ctz 1
#define FUNC_ipow 1
#include "util.h"
#include "totients.h"
#include "sieve.h"
#include "factor.h"
#include "keyval.h"


UV* range_totient(UV lo, UV hi) {
  UV* totients;
  UV i, seg_base, seg_low, seg_high, count = hi-lo+1;
  unsigned char* segment;
  void* ctx;

  if (hi < lo) croak("range_totient error hi %"UVuf" < lo %"UVuf"\n", hi, lo);
  New(0, totients, count, UV);

  /* Do via factoring if very small or if we have a small range */
  if (hi < 100 || count <= 10 || hi/count > 1000) {
    for (i = 0; i < count; i++)
      totients[i] = totient(lo+i);
    return totients;
  }

  if (hi == UV_MAX) {
    totients[--count] = totient(UV_MAX);
    hi--;
  }

  /* If doing a full sieve, do it monolithic.  Faster. */
  if (lo == 0) {
    UV* prime;
    double loghi = log(hi);
    UV max_index = (hi < 67)     ? 18
                 : (hi < 355991) ? 15+(hi/(loghi-1.09))
                 : (hi/loghi) * (1.0+1.0/loghi+2.51/(loghi*loghi));
    UV j, index, nprimes = 0;

    New(0, prime, max_index, UV);  /* could use prime_count_upper(hi) */
    memset(totients, 0, count * sizeof(UV));
    for (i = 2; i <= hi/2; i++) {
      index = 2*i;
      if ( !(i&1) ) {
        if (i == 2) { totients[2] = 1; prime[nprimes++] = 2; }
        totients[index] = totients[i]*2;
      } else {
        if (totients[i] == 0) {
          totients[i] = i-1;
          prime[nprimes++] = i;
        }
        for (j=0; j < nprimes && index <= hi; index = i*prime[++j]) {
          if (i % prime[j] == 0) {
            totients[index] = totients[i]*prime[j];
            break;
          } else {
            totients[index] = totients[i]*(prime[j]-1);
          }
        }
      }
    }
    Safefree(prime);
    /* All totient values have been filled in except the primes.  Mark them. */
    for (i = ((hi/2) + 1) | 1; i <= hi; i += 2)
      if (totients[i] == 0)
        totients[i] = i-1;
    totients[1] = 1;
    totients[0] = 0;
    return totients;
  }

  for (i = 0; i < count; i++) {
    UV v = lo+i, nv = v;
    if (v % 2 == 0)  nv -= nv/2;
    if (v % 3 == 0)  nv -= nv/3;
    if (v % 5 == 0)  nv -= nv/5;
    totients[i] = nv;
  }

  ctx = start_segment_primes(7, hi/2, &segment);
  while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
    START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high ) {
      for (i = P_GT_LO(2*p,p,lo); i >= lo && i <= hi; i += p)
        totients[i-lo] -= totients[i-lo]/p;
    } END_DO_FOR_EACH_SIEVE_PRIME
  }
  end_segment_primes(ctx);

  /* Fill in all primes */
  for (i = (lo | 1) - lo; i < count; i += 2)
    if (totients[i] == i+lo)
      totients[i]--;
  if (lo <= 1) totients[1-lo] = 1;

  return totients;
}

/******************************************************************************/

#define HAVE_SUMTOTIENT_128 (BITS_PER_WORD == 64 && HAVE_UINT128)
#define MAX_TOTSUM ( (BITS_PER_WORD == 64) ? UVCONST(7790208950) : 118868 )

typedef struct {
  UV    hsize;
  UV   *nhash;  /* n value */
  UV   *shash;  /* sum for n */
} sumt_hash_t;
#define _CACHED_SUMT(x) \
  (((x)<csize)  ?  cdata[x]  :  _sumt((x), cdata, csize, thash))
static UV _sumt(UV n, const UV *cdata, UV csize, sumt_hash_t thash) {
  UV sum, s, k, lim, hn = n % thash.hsize;

  if (n < csize) return cdata[n];
  if (thash.nhash[hn] == n)  return thash.shash[hn];

  sum = (n & 1)  ?  n*((n+1)>>1)  :  (n>>1)*(n+1);
  s = isqrt(n);
  lim = n/(s+1);

  sum -= (n - n/2) * _CACHED_SUMT(1);
  for (k = 2; k <= lim; k++) {
    sum -= _CACHED_SUMT(n/k);
    sum -= ((n/k) - (n/(k+1))) * _CACHED_SUMT(k);
  }
  if (s > lim)
    sum -= ((n/s) - (n/(s+1))) * _CACHED_SUMT(s);

  thash.nhash[hn] = n;
  thash.shash[hn] = sum;
  return sum;
}
UV sumtotient(UV n) {
  UV sum, i, csize, *sumcache;
  sumt_hash_t thash;

  if (n <= 2)  return n;
  if (n > MAX_TOTSUM)  return 0;

  if (n < 500) {    /* For very small values, do a simple sum */
    UV *phi = range_totient(0,n);
    for (sum = 0, i = 1; i <= n; i++)
      sum += phi[i];
    Safefree(phi);
    return sum;
  }

  csize = ((n < 1900) ? 2 : 1)  * icbrt(n) * icbrt(n);

  sumcache = range_totient(0, csize-1);
  for (i = 2; i < csize; i++)
    sumcache[i] += sumcache[i-1];

  thash.hsize = next_prime(10 + csize/16);
  Newz(0, thash.nhash, thash.hsize, UV);
  New( 0, thash.shash, thash.hsize, UV);

  sum = _sumt(n, sumcache, csize, thash);
  Safefree(thash.nhash);
  Safefree(thash.shash);
  Safefree(sumcache);
  return sum;
}

#if HAVE_SUMTOTIENT_128
#define _CACHED_SUMT128(x) \
  (((x)<csize)  ?  (uint128_t)cdata[x]  :  _sumt128((x), cdata, csize, thash))
typedef struct {
  UV         hsize;
  UV        *nhash;  /* n value */
  uint128_t *shash;  /* sum for n */
} sumt_hash_128_t;
static uint128_t _sumt128(UV n, const UV *cdata, UV csize, sumt_hash_128_t thash) {
  uint128_t sum;
  UV s, k, lim, hn = n % thash.hsize;

  if (n < csize) return cdata[n];
  if (thash.nhash[hn] == n)  return thash.shash[hn];

  sum = (n & 1)  ?  (uint128_t)n * ((n+1)>>1)  :  (uint128_t)(n+1) * (n>>1);
  s = isqrt(n);
  lim = n/(s+1);
  sum -= (n - n/2);
  for (k = 2; k <= lim; k++) {
    sum -= _CACHED_SUMT128(n/k);
    sum -= ((n/k) - (n/(k+1))) * _CACHED_SUMT128(k);
  }
  if (s > lim)
    sum -= ((n/s) - (n/(s+1))) * _CACHED_SUMT128(s);

  thash.nhash[hn] = n;
  thash.shash[hn] = sum;

  return sum;
}

int sumtotient128(UV n, UV *hi_sum, UV *lo_sum) {
  UV i, cbrtn, csize, *sumcache;
  uint128_t sum;
  sumt_hash_128_t thash;

  if (n <= 2)  { *hi_sum = 0;  *lo_sum = n; return 1; }

  cbrtn = icbrt(n);
  csize = 1 * cbrtn * cbrtn;
  if (csize > 400000000U)  csize = 400000000;  /* Limit to 3GB */

  sumcache = range_totient(0, csize-1);
  for (i = 2; i < csize; i++)
    sumcache[i] += sumcache[i-1];

#if 0
  thash.hsize = (n <= UVCONST( 10000000000000))  ?    500009
              : (n <= UVCONST(100000000000000))  ?   8000009
                                                 :  24000001;
#else
  thash.hsize = next_prime(10 + (cbrtn * cbrtn)/128);
  /* if (thash.hsize > 290000000U) thash.hsize=290000000U; */ /* Limit to 7GB */
#endif
  /* printf("hsize %lu  csize %lu  sumcache[%lu] = %lu\n", thash.hsize, csize, csize-1, sumcache[csize-1]); */
  Newz(0, thash.nhash, thash.hsize, UV);
  New( 0, thash.shash, thash.hsize, uint128_t);

  sum = _sumt128(n, sumcache, csize, thash);
  *hi_sum = (sum >> 64) & UV_MAX;
  *lo_sum = (sum      ) & UV_MAX;
  Safefree(thash.nhash);
  Safefree(thash.shash);
  Safefree(sumcache);
  return 1;
}
#else
int sumtotient128(UV n, UV *hi_sum, UV *lo_sum) {
  return 0;
}
#endif



UV totient_factored(UV n, UV nfacs, UV* fac, UV* exp) {
  UV i, totient = 1;
  if (n <= 1) return n;
  /* while ((n & 0x3) == 0) { n >>= 1; totient <<= 1; }         */
  /* if    ((n & 0x1) == 0) { n >>= 1; nfacs--; fac++; exp++; } */
  if (fac[0] == 2) {  /* Handle factors of 2 more efficiently */
    n       >>= exp[0];
    totient <<= (exp[0]-1);
    nfacs--; fac++; exp++;
  }
  for (i = 0; i < nfacs; i++) {
    UV f = fac[i], e = exp[i];
    totient *= f-1;
    while (e > 1) { totient *= f; e--; }
  }
  return totient;
}
#if 0
UV totient(UV n) {
  UV nfacs, fac[MPU_MAX_FACTORS+1], exp[MPU_MAX_FACTORS+1];
  if (n <= 1) return n;
  nfacs = factor_exp(n, fac, exp); /* factor and calculate totient */
  return totient_factored(n, nfacs, fac, exp);
}
#else
UV totient(UV n) {
  UV i, nfacs, totient, lastf, facs[MPU_MAX_FACTORS+1];
  if (n <= 1) return n;
  totient = 1;
  /* phi(2m) = 2phi(m) if m even, phi(m) if m odd */
  while ((n & 0x3) == 0) { n >>= 1; totient <<= 1; }
  if ((n & 0x1) == 0) { n >>= 1; }
  /* factor and calculate totient */
  nfacs = factor(n, facs);
  lastf = 0;
  for (i = 0; i < nfacs; i++) {
    UV f = facs[i];
    if (f == lastf) { totient *= f;               }
    else            { totient *= f-1;  lastf = f; }
  }
  return totient;
}
#endif

static const UV jordan_overflow[5] =
#if BITS_PER_WORD == 64
  {UVCONST(4294967311), 2642249, 65537, 7133, 1627};
#else
  {UVCONST(     65537),    1627,   257,   85,   41};
#endif
UV jordan_totient(UV k, UV n) {
  UV factors[MPU_MAX_FACTORS+1];
  int nfac, i;
  UV totient;
  if (k == 0 || n <= 1) return (n == 1);
  if (k > 6 || (k > 1 && n >= jordan_overflow[k-2])) return 0;

  totient = 1;
  /* Similar to Euler totient, shortcut even inputs */
  while ((n & 0x3) == 0) { n >>= 1; totient *= (1<<k); }
  if ((n & 0x1) == 0) { n >>= 1; totient *= ((1<<k)-1); }
  nfac = factor(n,factors);
  for (i = 0; i < nfac; i++) {
    UV p = factors[i];
    UV pk = ipow(p,k);
    totient *= (pk-1);
    while (i+1 < nfac && p == factors[i+1]) {
      i++;
      totient *= pk;
    }
  }
  return totient;
}



static int _totpred(UV n, UV maxd) {
  UV i, ndivisors, *divs;
  int res;

  if (n & 1) return 0;
  if ((n & (n-1)) == 0) return 1;
  n >>= 1;
  if (n == 1) return 1;
  if (n < maxd && is_prime(2*n+1)) return 1;

  divs = _divisor_list(n, &ndivisors);
  for (i = 0, res = 0; i < ndivisors && divs[i] < maxd && res == 0; i++) {
    UV r, d = divs[i], p = 2*d+1;
    if (!is_prime(p)) continue;
    r = n/d;
    while (1) {
      if (r == p || _totpred(r, d)) { res = 1; break; }
      if (r % p) break;
      r /= p;
    }
  }
  Safefree(divs);
  return res;
}

int is_totient(UV n) {
  return (n == 0 || (n & 1))  ?  (n==1)  :  _totpred(n,n);
}


UV inverse_totient_count(UV n) {
  set_t set, sumset;
  keyval_t keyval;
  UV res, i, ndivisors, *divs;

  if (n == 1) return 2;
  if (n < 1 || n & 1) return 0;
  if (is_prime(n >> 1)) { /* Coleman Remark 3.3 (Thm 3.1) and Prop 6.2 */
    if (!is_prime(n+1)) return 0;
    if (n >= 10)        return 2;
  }

  divs = _divisor_list(n, &ndivisors);

  init_set(&set, 2*ndivisors);
  keyval.key = 1;  keyval.val = 1;
  set_addsum(&set, keyval);

  for (i = 0; i < ndivisors; i++) {
    UV d = divs[i],  p = d+1;
    if (is_prime(p)) {
      UV j,  np = d,  v = valuation(n, p);
      init_set(&sumset, ndivisors/2);
      for (j = 0; j <= v; j++) {
        UV k, ndiv = n/np;  /* Loop over divisors of n/np */
        if (np == 1) {
          keyval_t kv;  kv.key = 1; kv.val = 1;
          set_addsum(&sumset, kv);
        } else {
          for (k = 0; k < ndivisors && divs[k] <= ndiv; k++) {
            UV val, d2 = divs[k];
            if ((ndiv % d2) != 0) continue;
            val = set_getval(set, d2);
            if (val > 0) {
              keyval_t kv;  kv.key = d2*np; kv.val = val;
              set_addsum(&sumset, kv);
            }
          }
        }
        /* if (j < v && np > UV_MAX/p) croak("overflow np  d %lu", d); */
        np *= p;
      }
      set_merge(&set, sumset);
      free_set(&sumset);
    }
  }
  Safefree(divs);
  res = set_getval(set, n);
  free_set(&set);
  return res;
}

UV* inverse_totient_list(UV *ntotients, UV n) {
  set_list_t setlist, divlist;
  UV i, ndivisors, *divs, *tlist;
  UV *totlist = 0;

  MPUassert(n <= UV_MAX/7.5, "inverse_totient_list n too large");

  if (n == 1) {
    New(0, totlist, 2, UV);
    totlist[0] = 1;  totlist[1] = 2;
    *ntotients = 2;
    return totlist;
  }
  if (n < 1 || n & 1) {
    *ntotients = 0;
    return totlist;
  }
  if (is_prime(n >> 1)) { /* Coleman Remark 3.3 (Thm 3.1) and Prop 6.2 */
    if (!is_prime(n+1)) {
      *ntotients = 0;
      return totlist;
    }
    if (n >= 10) {
      New(0, totlist, 2, UV);
      totlist[0] = n+1;  totlist[1] = 2*n+2;
      *ntotients = 2;
      return totlist;
    }
  }

  divs = _divisor_list(n, &ndivisors);

  init_setlist(&setlist, 2*ndivisors);
  setlist_addval(&setlist, 1, 1);   /* Add 1 => [1] */

  for (i = 0; i < ndivisors; i++) {
    UV d = divs[i],  p = d+1;
    if (is_prime(p)) {
      UV j,  dp = d,  pp = p,  v = valuation(n, p);
      init_setlist(&divlist, ndivisors/2);
      for (j = 0; j <= v; j++) {
        UV k, ndiv = n/dp;  /* Loop over divisors of n/dp */
        if (dp == 1) {
          setlist_addval(&divlist, 1, 2);   /* Add 1 => [2] */
        } else {
          for (k = 0; k < ndivisors && divs[k] <= ndiv; k++) {
            UV nvals, *vals, d2 = divs[k];
            if ((ndiv % d2) != 0) continue;
            vals = setlist_getlist(&nvals, setlist, d2);
            if (vals != 0)
              setlist_addlist(&divlist, d2 * dp, nvals, vals, pp);
          }
        }
        dp *= p;
        pp *= p;
      }
      setlist_merge(&setlist, divlist);
      free_setlist(&divlist);
    }
  }
  Safefree(divs);
  tlist = setlist_getlist(ntotients, setlist, n);
  if (tlist != 0 && *ntotients > 0) {
    New(0, totlist, *ntotients, UV);
    memcpy(totlist, tlist, *ntotients * sizeof(UV));
    qsort(totlist, *ntotients, sizeof(UV), _numcmp);
  }
  free_setlist(&setlist);
  return totlist;
}
