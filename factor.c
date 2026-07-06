#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define FUNC_isqrt  1
#define FUNC_gcd_ui 1
#define FUNC_is_perfect_square 1
#define FUNC_clz 1
#define FUNC_log2floor 1
#include "ptypes.h"
#include "constants.h"
#include "factor.h"
#include "sieve.h"
#include "util.h"
#include "sort.h"
#include "mulmod.h"
#include "cache.h"
#include "primality.h"
#include "lucas_seq.h"
#include "montmath.h"
static int holf32(uint32_t n, UV *factors, uint32_t rounds);

/*
 * You need to remember to use UV for unsigned and IV for signed types that
 * are large enough to hold our data.
 *   If you use int, that's 32-bit on LP64 and LLP64 machines.  You lose.
 *   If you use long, that's 32-bit on LLP64 machines.  You lose.
 *   If you use long long, you may be too large which isn't so bad, but some
 *                         compilers may not understand the type at all.
 * perl.h already figured all this out, and provided us with these types which
 * match the native integer type used inside our Perl, so just use those.
 */


/* For doing trial division loops over the small primes.
 * Returns either 1 or the new unfactored n.
 * Puts any factors in place and increments *nfactors.
 * Assumes n has no factors smaller than primes_small[sp].
 * Will check primes_small[sp] .. primes_small[endsp] inclusive.
 * endsp will be clamped to NPRIMES_SMALL-1.
 */
static uint32_t _trial32(uint32_t n, UV *factors, int *nfactors, uint32_t sp, uint32_t endsp) {
  uint32_t f;
  if (sp < 1) sp = 1;
  if (endsp > NPRIMES_SMALL-1) endsp = NPRIMES_SMALL-1;
  if (sp > endsp || n == 1) return n;

  do {
    f = primes_small[sp];
    if (f*f > n) break;
    while (n % f == 0) {
      factors[(*nfactors)++] = f;
      n /= f;
    }
  } while (++sp <= endsp);

  if (f*f > n && n != 1) {
    factors[(*nfactors)++] = n;
    n = 1;
  }
  return n;
}
static UV _trialuv(UV n, UV *factors, int *nfactors, uint32_t sp, uint32_t endsp) {
  uint32_t f;
  if (sp < 1) sp = 1;
  if (endsp > NPRIMES_SMALL-1) endsp = NPRIMES_SMALL-1;
  if (sp > endsp || n == 1) return n;

  do {
    f = primes_small[sp];
    if (f*f > n) break;
    while (n % f == 0) {
      factors[(*nfactors)++] = f;
      n /= f;
    }
  } while (++sp <= endsp);

  if (f*f > n && n != 1) {
    factors[(*nfactors)++] = n;
    n = 1;
  }
  return n;
}

static int _small_trial_factor(UV n, UV *factors, UV *newn, uint32_t *lastf)
{
  int nfactors = 0;
  uint32_t const endsp = 82;
  uint32_t sp = 4;
  uint32_t f = 7;

  if (n > 1) {
    while ( (n & 1) == 0 ) { factors[nfactors++] = 2; n /= 2; }
    while ( (n % 3) == 0 ) { factors[nfactors++] = 3; n /= 3; }
    while ( (n % 5) == 0 ) { factors[nfactors++] = 5; n /= 5; }
  }

  /* Trial primes 7 to 421 */
  n = (n <= 4294967295U) ? _trial32(n, factors, &nfactors, sp, endsp)
                         : _trialuv(n, factors, &nfactors, sp, endsp);
  sp = endsp+1; /* 83 */
  f = primes_small[sp];  /* 431 */

  if (n < 2017*2017 && f*f <= n) {  /* Trial division from 431 to 2011 */
    uint32_t const lastsp = NPRIMES_SMALL-1;
    n = _trial32(n, factors, &nfactors, sp, lastsp);
    f = 2017;
  }

  if (f*f > n && n != 1) {
    factors[nfactors++] = n;
    n = 1;
  }
  if (newn) *newn = n;
  if (lastf) *lastf = f;
  return nfactors;
}

static int _power_factor(UV n, UV *factors)
{
  uint32_t root;
  int nfactors, i, j, k;
  if (n > 3 && (k = powerof_ret(n, &root))) {
    nfactors = factor(root, factors);
    for (i = nfactors; i >= 0; i--)
      for (j = 0; j < k; j++)
        factors[k*i+j] = factors[i];
    return k * nfactors;
  }
  factors[0] = n;
  return 1;
}

/* Find one factor of an input n. */
int factor_one(UV n, UV *factors, bool primality, bool trial)
{
  int nfactors;
  if (n < 4) {
    factors[0] = n;
    return (n == 1) ? 0 : 1;
  }
  /* TODO: deal with small n */
  if (trial) {
    uint32_t sp, f;
    if (!(n&1)) { factors[0] = 2; factors[1] = n >> 1; return 2; }
    if (!(n%3)) { factors[0] = 3; factors[1] = n / 3; return 2; }
    if (!(n%5)) { factors[0] = 5; factors[1] = n / 5; return 2; }
    for (sp = 4; (f = primes_small[sp]) < 2011; sp++) {
      if ( (n % f) == 0 )
        { factors[0] = f; factors[1] = n/f; return 2; }
    }
    if (n < f*f) { factors[0] = n; return 1; }
  }
  if (primality && is_prime(n)) {
    factors[0] = n;
    return 1;
  }
#if 0  /* Simple solution, just fine on x86_64 */
  nfactors = (n < 1073741824UL) ? holf32(n, factors, 1000000)
                                : pbrent_factor(n, factors, 500000, 1);
  if (nfactors < 2) croak("factor_one failed on %lu\n", n);
#endif

#if BITS_PER_WORD == 64
    /* For small semiprimes the fastest solution is HOLF under 32, then
     * Lehman (no trial) under 38.  On random inputs, HOLF is best somewhere
     * between 28 and 32 bits.  Adding Lehman is always slower. */
    if (n <= 0xFFFFFFFFU) {
      nfactors = holf32(n, factors, 10000);  /* 2400 is enough */
      if (nfactors > 1) return nfactors;
    }
#endif

  {
    /* Adjust the number of rounds based on the number size and speed */
    UV const nbits = BITS_PER_WORD - clz(n);
#if USE_MONTMATH
    UV const br_rounds = 8000 + (9000 * ((nbits <= 45) ? 0 : (nbits-45)));
    UV const sq_rounds = 200000;
#elif MULMODS_ARE_FAST
    UV const br_rounds =  100 + ( 100 * ((nbits <= 45) ? 0 : (nbits-45)));
    UV const sq_rounds = 100000;
#else
    UV const br_rounds = (nbits >= 63) ? 120000 : (nbits >= 58) ? 500 : 0;
    UV const sq_rounds = 200000;
#endif
    /* Almost all inputs are factored here */
    if (br_rounds > 0) {
      nfactors = pbrent_factor(n, factors, br_rounds, 1);
      if (nfactors > 1) return nfactors;
    }
#if USE_MONTMATH
    nfactors = pbrent_factor(n, factors, 2*br_rounds, 3);
    if (nfactors > 1) return nfactors;
#endif
#if HAS_ECM64
    /* SQUFOF is still good for smaller inputs, better to skip for larger */
    if (nbits <= 48) {
      nfactors = squfof_factor(n, factors, sq_rounds);
      if (nfactors > 1) return nfactors;
    }

    nfactors = tinyecm64_factor(n, factors, 500, 5000, 20, 0);
    if (nfactors > 1) return nfactors;
    nfactors = tinyecm64_factor(n, factors, 1000, 10000, 20, 20);
    if (nfactors > 1) return nfactors;
    /* Essentially all 64-bit inputs have been found by this point. */
    nfactors = tinyecm64_factor(n, factors, 2000, 20000, 40, 40);
    if (nfactors > 1) return nfactors;
    nfactors = tinyecm64_factor(n, factors, 4000, 40000, 20, 80);
    if (nfactors > 1) return nfactors;
    nfactors = tinyecm64_factor(n, factors, 8000, 80000, 10, 100);
    if (nfactors > 1) return nfactors;
#else
    /* SQUFOF is useful in all cases before p-1 */
    if (nbits <= 62) {
      nfactors = squfof_factor(n, factors, sq_rounds);
      if (nfactors > 1) return nfactors;
    }

    nfactors = pminus1_factor(n, factors, 8000, 120000);
    if (nfactors > 1) return nfactors;
#endif

    /* Get the stragglers */
    nfactors = pbrent_factor(n, factors, 500000, 5);
    if (nfactors > 1) return nfactors;
    nfactors = prho_factor(n, factors, 180000);
    if (nfactors > 1) return nfactors;
    nfactors = cheb_factor(n, factors, 1000000, 0);
    if (nfactors > 1) return nfactors;
    croak("factor_one failed on %lu\n", n);
  }
  return nfactors;
}

/******************************************************************************/
/*                              Main factor loop                              */
/*                                                                            */
/* Puts factors in factors[] and returns the number found.                    */
/******************************************************************************/
int factor(UV n, UV *factors)
{
  UV tofac_stack[MPU_MAX_FACTORS+1];
  int nsmallfactors, npowerfactors, nfactors, i, j, ntofac = 0;
  uint32_t f;

  nfactors = _small_trial_factor(n, factors, &n, &f);
  if (n == 1) return nfactors;

#if BITS_PER_WORD == 64
  /* For small values less than f^3, use simple factor to split semiprime */
  if (n < 100000000 && n < f*f*f) {
    if (MR32(n)) factors[nfactors++] = n;
    else         nfactors += holf32(n, factors+nfactors, 10000);
    return nfactors;
  }
#endif

  nsmallfactors = nfactors;

  /* Perfect powers.  Factor root only once. */
  npowerfactors = _power_factor(n, factors+nsmallfactors);
  if (npowerfactors > 1) return nsmallfactors + npowerfactors;

  /* loop over each remaining factor, until ntofac == 0 */
  do {
    while ( (n >= f*f) && (!is_def_prime(n)) ) {
      int split_success = factor_one(n, tofac_stack+ntofac, 0, 0) - 1;
      if (split_success != 1 || tofac_stack[ntofac] == 1 || tofac_stack[ntofac] == n)
        croak("internal: factor_one failed to factor %"UVuf"\n", n);
      ntofac++; /* Leave one on the to-be-factored stack */
      n = tofac_stack[ntofac];  /* Set n to the other one */
    }
    /* n is now prime (or 1), so add to already-factored stack */
    if (n != 1)  factors[nfactors++] = n;
    /* Pop the next number off the to-factor stack */
    if (ntofac > 0)  n = tofac_stack[ntofac-1];
  } while (ntofac-- > 0);

  /* Sort the non-small factors */
  for (i = nsmallfactors+1; i < nfactors; i++) {
    UV fi = factors[i];
    for (j = i; j > 0 && factors[j-1] > fi; j--)
      factors[j] = factors[j-1];
    factors[j] = fi;
  }
  return nfactors;
}

void factorintp(factored_t *nf, UV n)
{
  UV fac[MPU_MAX_FACTORS], *f = nf->f;
  uint8_t *e = nf->e;
  uint32_t nfactors, i, j;

  nf->n = n;
  if (n < 4) {
    f[0] = n;
    e[0] = 1;
    nf->nfactors = 1 - (n==1);
    return;
  }
  nfactors = factor(n, fac);
  f[0] = fac[0];
  e[0] = 1;
  for (i = 1, j = 0; i < nfactors; i++) {
    if (fac[i] == fac[i-1])
      e[j]++;
    else
      f[++j] = fac[i], e[j] = 1;
  }
  nf->nfactors = (uint16_t)j+1;
}

void factoredp_validate(const factored_t *nf)
{
  if (nf->n == 0) {
    MPUassert(nf->nfactors == 1, "factored_t n=0  =>  nfactors = 1");
    MPUassert(nf->f[0] == 0 && nf->e[0] == 1, "factored_t n=0  =>  vecprod = n");
  } else if (nf->n == 1) {
    MPUassert(nf->nfactors == 0, "factored_t n=1  =>  nfactors = 0");
  } else {
    UV lf = 0, N = 1, t;
    uint32_t i;
    MPUassert(nf->nfactors <= MPU_MAX_DFACTORS, "factored_t n has too many factors");
    for (i = 0; i < nf->nfactors; i++) {
      MPUassert(is_prime(nf->f[i]), "factored_t n has non-prime factor");
      MPUassert(lf < nf->f[i], "factored_t factors not in order");
      lf = nf->f[i];
      MPUassert(nf->e[i] < BITS_PER_WORD, "factored_t exponent k too high");
      MPUassert(nf->e[i] > 0, "factored_t exponent k too low");
      if (nf->e[i] == 1) {
        N *= nf->f[i];
      } else {
        t = ipowsafe(nf->f[i], nf->e[i]);
        MPUassert(t != UV_MAX, "factored_t f^e overflows")
        N *= t;
      }
    }
    MPUassert(N == nf->n, "factored_t n is not equal to f^e * f^e ...");
  }
}
uint32_t factoredp_total_factors(const factored_t *nf) {
  uint32_t i, nfacs = 0;
  for (i = 0; i < nf->nfactors; i++)
    nfacs += nf->e[i];
  return nfacs;
}
bool factoredp_is_square_free(const factored_t *nf) {
  uint32_t i;
  for (i = 0; i < nf->nfactors; i++)
    if (nf->e[i] > 1)
      break;
  return i >= nf->nfactors;
}
signed char factoredp_moebius(const factored_t *nf) {
#if 0
  return !factoredp_is_square_free(nf) ? 0 : nf->nfactors % 2 ? -1 : 1;
#else
  uint32_t i;
  for (i = 0; i < nf->nfactors; i++)
    if (nf->e[i] > 1)
      return 0;
  return nf->nfactors % 2 ? -1 : 1;
#endif
}
uint32_t factoredp_linear_factors(UV fac[], const factored_t *nf) {
  uint32_t i, nfac = 0;
  for (i = 0; i < nf->nfactors; i++) {
    UV f = nf->f[i], e = nf->e[i];
    while (e--)
      fac[nfac++] = f;
  }
  return nfac;
}

/******************************************************************************/


int prime_bigomega(UV n)
{
  UV factors[MPU_MAX_FACTORS+1];
  return factor(n, factors);
}
int prime_omega(UV n)
{
  if (n <= 1) return (n==0);
  return factorint(n).nfactors;
}


int trial_factor(UV n, UV *factors, UV f, UV last)
{
  int sp, nfactors = 0;

  if (f < 2) f = 2;
  if (last == 0 || last*last > n) last = UV_MAX;

  if (n < 4 || last < f) {
    factors[0] = n;
    return (n == 1) ? 0 : 1;
  }

  /* possibly do uint32_t specific code here */

  if (f < primes_small[NPRIMES_SMALL-1]) {
    while ( (n & 1) == 0 ) { factors[nfactors++] = 2; n >>= 1; }
    if (3<=last) while ( (n % 3) == 0 ) { factors[nfactors++] = 3; n /= 3; }
    if (5<=last) while ( (n % 5) == 0 ) { factors[nfactors++] = 5; n /= 5; }
    for (sp = 4; sp < (int)NPRIMES_SMALL; sp++) {
      f = primes_small[sp];
      if (f*f > n || f > last) break;
      while ( (n%f) == 0 ) {
        factors[nfactors++] = f;
        n /= f;
      }
    }
  }
  /* Trial division using a mod-30 wheel for larger values */
  if (f*f <= n && f <= last) {
    UV m, newlimit, limit = isqrt(n);
    if (limit > last) limit = last;
    m = f % 30;
    while (f <= limit) {
      if ( (n%f) == 0 ) {
        do {
          factors[nfactors++] = f;
          n /= f;
        } while ( (n%f) == 0 );
        newlimit = isqrt(n);
        if (newlimit < limit)  limit = newlimit;
      }
      f += wheeladvance30[m];
      m = nextwheel30[m];
    }
  }
  /* All done! */
  if (n != 1)
    factors[nfactors++] = n;
  return nfactors;
}


static UV _divisors_from_factors(UV* res, factored_t nf, UV k) {
  UV count;
  uint32_t i;

  res[0] = count = 1;
  for (i = 0; i < nf.nfactors; i++) {
    UV s, scount = count, p = nf.f[i], mult = 1;
    uint32_t j, e = nf.e[i];
    for (j = 0; j < e; j++) {
      mult *= p;
      for (s = 0; s < scount; s++) {
        UV t = res[s] * mult;
        if (t <= k)
          res[count++] = t;
      }
    }
  }
  return count;
}

UV* divisor_list(UV n, UV *num_divisors, UV maxd)
{
  factored_t nf;
  UV ndivisors, *divs;
  uint32_t i;

  if (n == 0 || maxd == 0) {
    *num_divisors = 0;
    return 0;
  } else if (n == 1 || maxd == 1) {
    New(0, divs, 1, UV);
    divs[0] = 1;
    *num_divisors = 1;
    return divs;
  }

  if (maxd > n) maxd = n;

  /* Factor and convert to factor/exponent pair */
  nf = factorint(n);
  /* Calculate number of divisors, allocate space, fill with divisors */
  ndivisors = nf.e[0] + 1;
  for (i = 1; i < nf.nfactors; i++)
    ndivisors *= (nf.e[i] + 1);
  New(0, divs, ndivisors, UV);
  ndivisors = _divisors_from_factors(divs, nf, maxd);
  /* Sort divisors (numeric ascending) */
  sort_uv_array(divs, ndivisors);
  /* Return number of divisors and list */
  *num_divisors = ndivisors;
  return divs;
}


/* The usual method, on OEIS for instance, is:
 *    (p^(k*(e+1))-1) / (p^k-1)
 * but that overflows quicky.  Instead we rearrange as:
 *    1 + p^k + p^k^2 + ... p^k^e
 * Return 0 if the result overflowed.
 */
static const UV sigma_overflow[11] =
#if BITS_PER_WORD == 64
         {UVCONST(3000000000000000000),UVCONST(3000000000),2487240,64260,7026,
         1622, 566, 256, 139, 85, 57};
#else
         {UVCONST(845404560), 52560, 1548, 252, 84, 41, 24, 16, 12, 10, 8};
#endif
UV divisor_sum(UV n, UV k)
{
  factored_t nf;
  UV product;
  uint32_t i, j;

  if (k > 11 || (k > 0 && n >= sigma_overflow[k-1])) return 0;
  /*   divisors(0) = []   divisors(1) = [1]  */
  if (n <= 1)  return n;
  nf = factorint(n);
  product = 1;
  if (k == 0) {
    for (i = 0; i < nf.nfactors; i++)
      product *= (nf.e[i]+1);
  } else if (k == 1) {
    for (i = 0; i < nf.nfactors; i++) {
      UV       f = nf.f[i];
      uint16_t e = nf.e[i];
      UV pke = f, fmult = 1 + f;
      while (e-- > 1) {
        pke *= f;
        fmult += pke;
      }
      product *= fmult;
    }
  } else {
    for (i = 0; i < nf.nfactors; i++) {
      UV       f = nf.f[i];
      uint16_t e = nf.e[i];
      UV fmult, pke, pk = f;
      for (j = 1; j < k; j++)  pk *= f;
      fmult = 1 + pk;
      pke = pk;
      while (e-- > 1) {
        pke *= pk;
        fmult += pke;
      }
      product *= fmult;
    }
  }
  return product;
}


UV aliquot_sum(UV n) {
  factored_t nf;
  UV product;
  uint32_t i;

  if (n <= 1) return 0;

  nf = factorint(n);

  /* Prime power: aliquot(p^e) = 1 + p + ... + p^(e-1), always < p^e = n */
  if (nf.nfactors == 1) {
    UV f = nf.f[0];
    uint16_t e = nf.e[0];
    UV pke = 1, fmult = 1;
    while (e-- > 1) {
      pke *= f;
      fmult += pke;
    }
    return fmult;
  }

  /* Multiple prime factors: compute σ(n) with overflow detection, subtract n */
  product = 1;
  for (i = 0; i < nf.nfactors; i++) {
    UV       f = nf.f[i];
    uint16_t e = nf.e[i];
    UV pke = f, fmult = 1 + f;
    while (e-- > 1) {
      pke *= f;
      if (fmult > UV_MAX - pke) return 0;
      fmult += pke;
    }
    if (product > UV_MAX / fmult) return 0;
    product *= fmult;
  }
  return product - n;
}


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


static int found_factor(UV n, UV f, UV* factors)
{
  UV g = n/f;
  if (f == 1 || f == n) {
    factors[0] = n;
    return 1;
  }
  factors[f >= g] = f;
  factors[f <  g] = g;
  MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
  return 2;
}
static int no_factor(UV n, UV* factors)
{
  factors[0] = n;
  return 1;
}

/* Knuth volume 2, algorithm C.
 * Can't compete with HOLF, SQUFOF, pbrent, etc.
 */
int fermat_factor(UV n, UV *factors, UV rounds)
{
  IV sqn, x, y, r;
  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in fermat_factor");
  sqn = isqrt(n);
  x = 2 * sqn + 1;
  y = 1;
  r = (sqn*sqn) - n;

  while (r != 0) {
    if (rounds-- == 0) return no_factor(n,factors);
    r += x;
    x += 2;
    do {
      r -= y;
      y += 2;
    } while (r > 0);
  }
  r = (x-y)/2;
  return found_factor(n, r, factors);
}

/* Hart's One Line Factorization. */
int holf_factor(UV n, UV *factors, UV rounds)
{
  UV i, s, m, f;
  uint32_t root;

  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in holf_factor");

  /* We skip the perfect-square test for s in the loop, so we
   * will never succeed if n is a perfect square.  Test that now. */
  if (is_perfect_square_ret(n,&root))
    return found_factor(n, root, factors);

  if (n <= (UV_MAX >> 6)) {    /* Try with premultiplier first */
    UV npre = n * ( (n <= (UV_MAX >> 13)) ? 720 :
                    (n <= (UV_MAX >> 11)) ? 480 :
                    (n <= (UV_MAX >> 10)) ? 360 :
                    (n <= (UV_MAX >>  8)) ?  60 : 30 );
    UV ni = npre;
    while (rounds--) {
      s = 1 + (UV)isqrt(ni);
      m = (s*s) - ni;
      if (is_perfect_square_ret(m, &root)) {
        f = gcd_ui(n, s - root);
        if (f > 1 && f < n)
          return found_factor(n, f, factors);
      }
      if (ni >= (ni+npre)) break;
      ni += npre;
    }
    if (rounds == (UV) -1)
      return no_factor(n,factors);
  }

  for (i = 1; i <= rounds; i++) {
    s = (UV) sqrt( (double)n * (double)i );
    /* Assume s^2 isn't a perfect square.  We're rapidly losing precision
     * so we won't be able to accurately detect it anyway. */
    s++;    /* s = ceil(sqrt(n*i)) */
    m = sqrmod(s, n);
    if (is_perfect_square_ret(m, &root)) {
      f = gcd_ui( (s>root) ? s-root : root-s, n);
      /* This should always succeed, but with overflow concerns.... */
      return found_factor(n, f, factors);
    }
  }
  return no_factor(n,factors);
}
static int holf32(uint32_t n, UV *factors, uint32_t rounds) {
  UV npre, ni;    /* These should be 64-bit */
  uint32_t s, m, f;

  if (n < 3) return no_factor(n,factors);
  if (!(n&1)) { factors[0] = 2; factors[1] = n/2; return 2; }
  if (is_perfect_square_ret(n,&f)) { factors[0] = factors[1] = f; return 2; }

  ni = npre = (UV) n * ((BITS_PER_WORD == 64) ? 5040 : 1);
  while (rounds--) {
    s = 1 + isqrt(ni);
    m = ((UV)s*(UV)s) - ni;
    if (is_perfect_square_ret(m, &f)) {
      f = gcd_ui(n, s - f);
      if (f > 1 && f < n)
        return found_factor(n, f, factors);
    }
    if (ni >= (ni+npre)) break; /* We've overflowed */
    ni += npre;
  }
  return no_factor(n,factors);
}


#define ABSDIFF(x,y)  (x>y) ? x-y : y-x
#if USE_MONTMATH
/* Pollard Rho with Brent's updates, using Montgomery reduction. */
int pbrent_factor(UV n, UV *factors, UV rounds, UV a)
{
  UV const nbits = BITS_PER_WORD - clz(n);
  const UV inner = (nbits <= 31) ? 32 : (nbits <= 35) ? 64 : (nbits <= 40) ? 160 : (nbits <= 52) ? 256 : 320;
  UV f, m, r, rleft, Xi, Xm, Xs;
  int irounds, fails = 6;
  const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);

  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in pbrent_factor");
  r = f = 1;
  Xi = Xm = Xs = mont1;
  a = mont_geta(a,n);

  while (rounds > 0) {
    rleft = (r > rounds) ? rounds : r;
    Xm = Xi;
    /* Do rleft rounds, inner at a time */
    while (rleft > 0) {
      irounds = (rleft > (UV)inner) ? inner : rleft;
      rleft -= irounds;
      rounds -= irounds;
      Xs = Xi;
      if (n < (1ULL << 63)) {
        Xi = mont_mulmod63(Xi,Xi+a,n);
        m = ABSDIFF(Xi,Xm);
        while (--irounds > 0) {
          Xi = mont_mulmod63(Xi,Xi+a,n);
          f = ABSDIFF(Xi,Xm);
          m = mont_mulmod63(m, f, n);
        }
      } else if (a == mont1) {
        Xi = mont_mulmod64(Xi,Xi+a,n);
        m = ABSDIFF(Xi,Xm);
        while (--irounds > 0) {
          Xi = mont_mulmod64(Xi,Xi+a,n);
          f = ABSDIFF(Xi,Xm);
          m = mont_mulmod64(m, f, n);
        }
      } else {
        Xi = addmod(mont_mulmod64(Xi,Xi,n), a, n);
        m = ABSDIFF(Xi,Xm);
        while (--irounds > 0) {
          Xi = addmod(mont_mulmod64(Xi,Xi,n), a, n);
          f = ABSDIFF(Xi,Xm);
          m = mont_mulmod64(m, f, n);
        }
      }
      f = gcd_ui(m, n);
      if (f != 1)
        break;
    }
    /* If f == 1, then we didn't find a factor.  Move on. */
    if (f == 1) {
      r *= 2;
      continue;
    }
    if (f == n) {  /* back up, with safety */
      Xi = Xs;
      do {
        if (n < (1ULL << 63) || a == mont1)
          Xi = mont_mulmod(Xi,Xi+a,n);
        else
          Xi = addmod(mont_mulmod(Xi,Xi,n),a,n);
        m = ABSDIFF(Xi,Xm);
        f = gcd_ui(m, n);
      } while (f == 1 && r-- != 0);
    }
    if (f == 0 || f == n) {
      if (fails-- <= 0) break;
      Xi = Xm = mont1;
      a = addmod(a, mont_geta(11,n), n);
      continue;
    }
    return found_factor(n, f, factors);
  }
  return no_factor(n,factors);
}
#else
/* Pollard Rho with Brent's updates. */
int pbrent_factor(UV n, UV *factors, UV rounds, UV a)
{
  UV f, m, r, Xi, Xm;
  const UV inner = (n <= 4000000000UL) ? 32 : 160;
  int fails = 6;

  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in pbrent_factor");

  r = f = Xi = Xm = 1;
  while (rounds > 0) {
    UV rleft = (r > rounds) ? rounds : r;
    UV saveXi = Xi;
    /* Do rleft rounds, inner at a time */
    while (rleft > 0) {
      UV dorounds = (rleft > inner) ? inner : rleft;
      saveXi = Xi;
      rleft -= dorounds;
      rounds -= dorounds;
      Xi = sqraddmod(Xi, a, n);        /* First iteration, no mulmod needed */
      m = ABSDIFF(Xi,Xm);
      while (--dorounds > 0) {         /* Now do inner-1=63 more iterations */
        Xi = sqraddmod(Xi, a, n);
        f = ABSDIFF(Xi,Xm);
        m = mulmod(m, f, n);
      }
      f = gcd_ui(m, n);
      if (f != 1)
        break;
    }
    /* If f == 1, then we didn't find a factor.  Move on. */
    if (f == 1) {
      r *= 2;
      Xm = Xi;
      continue;
    }
    if (f == n) {  /* back up, with safety */
      Xi = saveXi;
      do {
        Xi = sqraddmod(Xi, a, n);
        f = gcd_ui( ABSDIFF(Xi,Xm), n);
      } while (f == 1 && r-- != 0);
    }
    if (f == 0 || f == n) {
      if (fails-- <= 0) break;
      Xm = addmod(Xm, 11, n);
      Xi = Xm;
      a++;
      continue;
    }
    return found_factor(n, f, factors);
  }
  return no_factor(n,factors);
}
#endif

/* Pollard's Rho. */
int prho_factor(UV n, UV *factors, UV rounds)
{
  UV f, i, m, oldU, oldV;
  const UV inner = 64;
  UV U = 7;
  UV V = 7;
  UV a = 1;
  int fails = 3;

  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in prho_factor");

  rounds = (rounds + inner - 1) / inner;

  while (rounds-- > 0) {
    m = 1; oldU = U; oldV = V;
    for (i = 0; i < inner; i++) {
      U = sqraddmod(U, a, n);
      V = sqraddmod(V, a, n);
      V = sqraddmod(V, a, n);
      f = (U > V) ? U-V : V-U;
      m = mulmod(m, f, n);
    }
    f = gcd_ui(m, n);
    if (f == 1)
      continue;
    if (f == n) {  /* back up to find a factor*/
      U = oldU; V = oldV;
      i = inner;
      do {
        U = sqraddmod(U, a, n);
        V = sqraddmod(V, a, n);
        V = sqraddmod(V, a, n);
        f = gcd_ui( (U > V) ? U-V : V-U, n);
      } while (f == 1 && i-- != 0);
    }
    if (f == 0 || f == n) {
      if (fails-- <= 0) break;
      U = addmod(U,2,n);
      V = U;
      a += 2;
      continue;
    }
    return found_factor(n, f, factors);
  }
  return no_factor(n,factors);
}


/******************************************************************************/
/* Tiny ECM -- elliptic curve factoring, 64-bit prototype                     */
/******************************************************************************/

#if HAS_ECM64

static INLINE uint64_t tecm64_addmod(uint64_t a, uint64_t b, uint64_t n)
{
  uint64_t r = a + b;
  if (r < a || r >= n) r -= n;
  return r;
}

static INLINE uint64_t tecm64_submod(uint64_t a, uint64_t b, uint64_t n)
{
  return (a >= b) ? a - b : n - (b - a);
}

static INLINE uint64_t tecm64_halfmod(uint64_t a, uint64_t n)
{
  return (a & 1)  ?  (a >> 1) + (n >> 1) + 1  :  a >> 1;
}

static uint64_t tecm64_mulmod(uint64_t a, uint64_t b, uint64_t n)
{
#if HAVE_UINT128
  return (uint64_t)(((uint128_t)a * b) % n);
#else
  uint64_t r = 0;
  if (a >= n) a %= n;
  if (b >= n) b %= n;
  if (a < b) { uint64_t t = a; a = b; b = t; }
  while (b) {
    if (b & 1) r = tecm64_addmod(r, a, n);
    b >>= 1;
    if (b) a = tecm64_addmod(a, a, n);
  }
  return r;
#endif
}

#define tecm64_sqrmod(a,n) tecm64_mulmod(a,a,n)

static uint64_t tecm64_powmod(uint64_t a, UV k, uint64_t n)
{
  uint64_t r = 1;
  if (a >= n) a %= n;
  while (k) {
    if (k & 1) r = tecm64_mulmod(r, a, n);
    k >>= 1;
    if (k) a = tecm64_sqrmod(a, n);
  }
  return r;
}

static uint64_t tecm64_gcd(uint64_t a, uint64_t b)
{
  while (b) { uint64_t t = b; b = a % b; a = t; }
  return a;
}

/* Binary extended modular inverse for odd n.  Returns 0 if gcd(a,n) > 1. */
static uint64_t tecm64_modinv(uint64_t a, uint64_t n)
{
  uint64_t u, v, x1, x2;
  if (a == 0) return 0;
  u = a % n;  v = n;
  x1 = 1;     x2 = 0;
  while (u != 1 && v != 1) {
    while ((u & 1) == 0) {
      u >>= 1;
      x1 = (x1 & 1) ? tecm64_halfmod(x1, n) : x1 >> 1;
    }
    while ((v & 1) == 0) {
      v >>= 1;
      x2 = (x2 & 1) ? tecm64_halfmod(x2, n) : x2 >> 1;
    }
    if (u >= v) {
      u -= v;
      x1 = tecm64_submod(x1, x2, n);
    } else {
      v -= u;
      x2 = tecm64_submod(x2, x1, n);
    }
    if (u == 0 || v == 0) return 0;
  }
  return (u == 1) ? x1 : x2;
}

typedef struct {
  uint64_t n;
#if HAVE_UINT128
  uint64_t ninv, r2;
#endif
} mont64_t;

static void tecm64_mont_setup(mont64_t *ctx, uint64_t n)
{
  ctx->n = n;
#if HAVE_UINT128
  {
    uint64_t x = (3*n)^2;
    x *= (uint64_t)2 - n * x;
    x *= (uint64_t)2 - n * x;
    x *= (uint64_t)2 - n * x;
    x *= (uint64_t)2 - n * x;
    ctx->ninv = (uint64_t)0 - x;
    ctx->r2 = tecm64_powmod(2, 128, n);
  }
#endif
}

static INLINE uint64_t tecm64_mont_mulmod(uint64_t a, uint64_t b,
                                          const mont64_t *ctx)
{
#if HAVE_UINT128
  /* REDC needs the carry above 2^128 when n > 2^63. */
  uint128_t ab = (uint128_t)a * b;
  uint64_t lo = (uint64_t)ab;
  uint64_t hi = (uint64_t)(ab >> 64);
  uint64_t m = lo * ctx->ninv;
  uint64_t mn_hi = (uint64_t)(((uint128_t)m * ctx->n) >> 64);
  uint64_t u = hi + mn_hi;
  int ov = (u < hi);
  u += (lo != 0);
  ov += (u < (uint64_t)(lo != 0));
  if (ov || u >= ctx->n) u -= ctx->n;
  return u;
#else
  return tecm64_mulmod(a, b, ctx->n);
#endif
}

#define tecm64_mont_sqrmod(a,ctx) tecm64_mont_mulmod(a,a,ctx)

static INLINE uint64_t tecm64_mont_enter(uint64_t a, const mont64_t *ctx)
{
#if HAVE_UINT128
  return tecm64_mont_mulmod(a % ctx->n, ctx->r2, ctx);
#else
  return a % ctx->n;
#endif
}

static INLINE uint64_t tecm64_mont_exit(uint64_t a, const mont64_t *ctx)
{
#if HAVE_UINT128
  return tecm64_mont_mulmod(a, 1, ctx);
#else
  return a;
#endif
}

typedef struct { uint64_t X, Z; } ecpt64_t;

static INLINE void ecm_double64(ecpt64_t *R, const ecpt64_t *P,
                                uint64_t A24, const mont64_t *ctx)
{
  uint64_t n = ctx->n;
  uint64_t u = tecm64_mont_sqrmod(tecm64_submod(P->X, P->Z, n), ctx);
  uint64_t v = tecm64_mont_sqrmod(tecm64_addmod(P->X, P->Z, n), ctx);
  uint64_t w = tecm64_submod(v, u, n);
  R->X = tecm64_mont_mulmod(u, v, ctx);
  R->Z = tecm64_mont_mulmod(w, tecm64_addmod(u, tecm64_mont_mulmod(A24, w, ctx), n), ctx);
}

static INLINE void ecm_dadd64(ecpt64_t *R,
                              const ecpt64_t *P, const ecpt64_t *Q,
                              const ecpt64_t *Pm, const mont64_t *ctx)
{
  uint64_t n = ctx->n;
  uint64_t u = tecm64_mont_mulmod(tecm64_submod(P->X, P->Z, n),
                                  tecm64_addmod(Q->X, Q->Z, n), ctx);
  uint64_t v = tecm64_mont_mulmod(tecm64_addmod(P->X, P->Z, n),
                                  tecm64_submod(Q->X, Q->Z, n), ctx);
  uint64_t s = tecm64_addmod(u, v, n);
  uint64_t d = tecm64_submod(u, v, n);
  R->X = tecm64_mont_mulmod(tecm64_mont_sqrmod(s, ctx), Pm->Z, ctx);
  R->Z = tecm64_mont_mulmod(tecm64_mont_sqrmod(d, ctx), Pm->X, ctx);
}

static void ecm_mul64(ecpt64_t *R, const ecpt64_t *P, UV k,
                      uint64_t A24, const mont64_t *ctx)
{
  ecpt64_t R0, R1;
  UV bit;
  if (k == 1) { *R = *P; return; }
  R0 = *P;
  ecm_double64(&R1, &R0, A24, ctx);
  if (k == 2) { *R = R1; return; }
  bit = (UV)1 << (8*sizeof(UV) - 1);
  while (!(bit & k)) bit >>= 1;
  for (bit >>= 1; bit; bit >>= 1) {
    if (k & bit) {
      ecm_dadd64(&R0, &R0, &R1, P, ctx);
      ecm_double64(&R1, &R1, A24, ctx);
    } else {
      ecm_dadd64(&R1, &R0, &R1, P, ctx);
      ecm_double64(&R0, &R0, A24, ctx);
    }
  }
  *R = R0;
}

/* Same fixed Suyama sigma sequence as tinyecm128. */
static const uint16_t ecm64_sigmas[] = {
   11,   13,   17,   19,   23,   29,   31,   37,   41,   43,
   47,   53,   59,   61,   67,   71,   73,   79,   83,   89,
  103,  127,  139,  149,  151,  157,  163,  181,  191,  197,
  199,  211,  223,  227,  233,  239,  257,  271,  293,  313,
  331,  337,  347,  359,  379,  389,  397,  401,  409,  421,
  443,  449,  457,  479,  487,  509,  521,  523,  547,  557,
  587,  641,  653,  659,  673,  677,  683,  691,  719,  727,
  739,  751,  769,  797,  809,  853,  919,  929,  941,  997,
 1049, 1051, 1063, 1091, 1093, 1109, 1117, 1129, 1153, 1201,
 1217, 1229, 1283, 1327, 1361, 1381, 1427, 1447, 1459, 1471,
 1481, 1489, 1543, 1549, 1571, 1621, 1709, 1723, 1753, 1759,
 1801, 1811, 1867, 1987, 2039, 2099, 2113, 2131, 2251, 2309,
 2347, 2381, 2399, 2447, 2473, 2551, 2557, 2663, 2677, 2689,
 2713, 2719, 2749, 2857, 2879, 2887, 2939, 3001, 3061, 3067,
 3121, 3137, 3187, 3251, 3259, 3271, 3307, 3359, 3371, 3373,
 3467, 3593, 3607, 3623, 3643, 3709, 3733, 3793, 3851, 3923,
 3989, 4019, 4049, 4129, 4231, 4253, 4283, 4339, 4349, 4441,
 4523, 4649, 4787, 4987, 4999, 5171, 5237, 5273, 5297, 5333,
 5387, 5471, 5479, 5647, 5749, 5791, 6101, 6163, 6257, 6299,
 6337, 6451, 6491, 6659, 6793, 6823, 6967, 7013, 7229, 7253,
 7333, 7369, 7477, 7621, 7793, 7817, 8059, 8167, 8209, 8263,
 8311, 8377, 8573, 8641, 8741, 8837, 8863, 8963, 9001, 9151,
 9203, 9433, 9697, 9743, 9781, 9883,10007,10069,10099,10139,
10163,10193,10267,10429,10457,10487,10691,10837,10949,11087,
11243,11321,11411,11681,11813,11903,12011,12263,12277,12401,
12409,12437,12479,12569,12619,12739,12911,13331,13367,13537,
13721,13789,13841,13873,14051,14149,14221,14419,14431,14827,
14887,15077,15289,15467,15511,15649,15773,15797,15859,15901,
16057,16141,16217,16529,16547,16553,16619,17299,17393,17419,
17449,17737,17921,18049,18223,19073,19183,19477,20021,20323,
20347,20759,20929,21023,21157,21587,21611,21613,21673,21751,
21799,21821,22109,22469,22651,22943,23327,23459,23567,23767,
23911,23957,24001,24197,24281,24407,24799,24851,25147,25183,
25469,25679,25703,26561,26683,26701,26821,27073,27191,27271,
27277,27427,27487,27539,27617,27647,27673,27749,27983,28319,
28789,28843,29017,29123,29209,29669,29803,29921,30323,30809,
30851,30911,30983,31397,31541,31963,32369,32561,32771,32969,
33029,33083,33487,33637,33757,34057,34381,34513,34613,34807,
35083,35171,35311,35381,36013,36251,36493,36529,36551,36913,
36919,37363,37517,37699,37907,38047,38177,38273,38749,38903
};
#define NECM64_SIGMAS ((int)(sizeof(ecm64_sigmas)/sizeof(ecm64_sigmas[0])))

static int ecm_batch_normalize_x64(uint64_t *xout, uint64_t *fout,
                                   const ecpt64_t *P, UV npoints,
                                   const mont64_t *ctx)
{
  uint64_t n = ctx->n, acc = 1, inv, g;
  uint64_t *zv, *prefix;
  UV i;

  if (npoints == 0) { *fout = 0; return 1; }

  New(0, zv,     npoints, uint64_t);
  New(0, prefix, npoints, uint64_t);

  for (i = 0; i < npoints; i++) {
    zv[i] = tecm64_mont_exit(P[i].Z, ctx);
    g = tecm64_gcd(zv[i], n);
    if (g > 1) {
      *fout = (g < n) ? g : 0;
      Safefree(prefix);
      Safefree(zv);
      return 0;
    }
    acc = tecm64_mulmod(acc, zv[i], n);
    prefix[i] = acc;
  }

  g = tecm64_gcd(acc, n);
  if (g > 1) {
    *fout = (g < n) ? g : 0;
    Safefree(prefix);
    Safefree(zv);
    return 0;
  }

  inv = tecm64_modinv(acc, n);
  if (inv == 0) {
    *fout = 0;
    Safefree(prefix);
    Safefree(zv);
    return 0;
  }

  for (i = npoints; i > 0; i--) {
    UV j = i - 1;
    uint64_t prev = (j == 0) ? 1 : prefix[j-1];
    uint64_t zinv = tecm64_mulmod(inv, prev, n);
    inv = tecm64_mulmod(inv, zv[j], n);
    xout[j] = tecm64_mont_mulmod(P[j].X, tecm64_mont_enter(zinv, ctx), ctx);
  }

  *fout = 0;
  Safefree(prefix);
  Safefree(zv);
  return 1;
}

static uint64_t tinyecm64_stage2(const ecpt64_t *Q, uint64_t A24,
                                 const mont64_t *ctx, UV B1, UV B2)
{
  uint64_t n = ctx->n, f, *nqx, *Sx;
  ecpt64_t *nq, *S;
  UV D, twoD, m, mend, nwindows, w, i;
  uint64_t gprod;

  if (B2 <= B1) return 0;
  D = isqrt(B2 >> 1);
  if (D & 1) D++;
  if (D == 0 || D > (UV_MAX-1)/2) return 0;
  twoD = 2*D;

  mend = (B2 > UV_MAX-D) ? UV_MAX : B2 + D;
  for (m = 1, nwindows = 0; m < mend; nwindows++) {
    if (m > UV_MAX-twoD) break;
    m += twoD;
  }
  if (nwindows == 0) return 0;

  New(0, nq,  twoD+1, ecpt64_t);
  New(0, nqx, D+1, uint64_t);

  nq[1] = *Q;
  for (i = 2; i <= twoD; i++) {
    if (i & 1) {
      ecm_dadd64(&nq[i], &nq[(i-1)/2], &nq[(i+1)/2], Q, ctx);
    } else {
      ecm_double64(&nq[i], &nq[i/2], A24, ctx);
    }
  }

  nqx[0] = 0;
  if (!ecm_batch_normalize_x64(nqx+1, &f, nq+1, D, ctx)) {
    Safefree(nqx);
    Safefree(nq);
    return f;
  }

  New(0, S,  nwindows, ecpt64_t);
  New(0, Sx, nwindows, uint64_t);

  S[0] = *Q;
  {
    ecpt64_t Xm = nq[twoD-1];
    for (w = 1; w < nwindows; w++) {
      ecpt64_t oldS = S[w-1];
      ecm_dadd64(&S[w], &nq[twoD], &S[w-1], &Xm, ctx);
      Xm = oldS;
    }
  }

  if (!ecm_batch_normalize_x64(Sx, &f, S, nwindows, ctx)) {
    Safefree(Sx);
    Safefree(S);
    Safefree(nqx);
    Safefree(nq);
    return f;
  }

  gprod = tecm64_mont_enter(1, ctx);
  m = 1;
  for (w = 0; w < nwindows; w++) {
    UV hi = (m > UV_MAX-D) ? UV_MAX : m + D;
    UV lo = (m > D) ? m - D : 0;

    if (hi > B2) hi = B2;
    if (hi > B1) {
      if (lo <= B1) lo = B1 + 1;
      if (lo < 2) lo = 2;
      if (lo <= hi) {
        START_DO_FOR_EACH_PRIME(lo,hi) {
          UV idx;
          if (p < m) {
            idx = m - p;
          } else if (p > m) {
            if (m <= UV_MAX/2) {
              UV mm = 2*m;
              if (p <= mm) {
                UV mirror = mm - p;
                if (mirror > B1 && mirror >= lo && is_prime(mirror))
                  continue;
              }
            }
            idx = p - m;
          } else {
            continue;
          }
          if (idx <= D) {
            uint64_t diff = tecm64_submod(Sx[w], nqx[idx], n);
            gprod = tecm64_mont_mulmod(gprod, diff, ctx);
          }
        } END_DO_FOR_EACH_PRIME
        f = tecm64_gcd(tecm64_mont_exit(gprod, ctx), n);
        if (f > 1) {
          Safefree(Sx);
          Safefree(S);
          Safefree(nqx);
          Safefree(nq);
          return (f < n) ? f : 0;
        }
      }
    }
    if (m > UV_MAX-twoD) break;
    m += twoD;
  }

  f = tecm64_gcd(tecm64_mont_exit(gprod, ctx), n);
  Safefree(Sx);
  Safefree(S);
  Safefree(nqx);
  Safefree(nq);
  return (f > 1 && f < n) ? f : 0;
}

static uint64_t tinyecm64(uint64_t n, UV B1, UV B2,
                          uint32_t ncurves, uint32_t sigma_offset)
{
  mont64_t ctx;
  UV sqrtB1, j;
  uint32_t ci;

  if (n < 3 || (n & 1) == 0) return 0;
  tecm64_mont_setup(&ctx, n);
  sqrtB1 = (UV)sqrt((double)B1);

  for (ci = 0; ci < ncurves && sigma_offset+ci < NECM64_SIGMAS; ci++) {
    uint64_t sigma, ui, vi, umv, t3uv;
    uint64_t u2, u3, v2, v3, umv2, umv3;
    uint64_t abs_num_r, num, den_r, den_inv, A24;
    ecpt64_t P;

    sigma = ecm64_sigmas[sigma_offset+ci];
    ui = sigma * sigma - 5;
    vi = 4 * sigma;
    umv = ui - vi;
    t3uv = 3 * ui + vi;

    u2 = tecm64_mulmod(ui % n, ui % n, n);
    u3 = tecm64_mulmod(u2, ui % n, n);
    v2 = tecm64_mulmod(vi % n, vi % n, n);
    v3 = tecm64_mulmod(v2, vi % n, n);
    umv2 = tecm64_mulmod(umv % n, umv % n, n);
    umv3 = tecm64_mulmod(umv2, umv % n, n);

    abs_num_r = tecm64_mulmod(umv3, t3uv % n, n);
    num = abs_num_r ? n - abs_num_r : 0;
    den_r = tecm64_mulmod(tecm64_mulmod(u3, vi % n, n), 16 % n, n);
    den_inv = tecm64_modinv(den_r, n);
    if (den_inv == 0) continue;
    A24 = tecm64_mont_mulmod(tecm64_mont_enter(num, &ctx),
                             tecm64_mont_enter(den_inv, &ctx), &ctx);

    P.X = tecm64_mont_enter(u3, &ctx);
    P.Z = tecm64_mont_enter(v3, &ctx);

    j = 0;
    START_DO_FOR_EACH_PRIME(2,B1) {
      UV k = p;
      if (p <= sqrtB1) { UV pm = B1 / p; while (k <= pm) k *= p; }
      ecm_mul64(&P, &P, k, A24, &ctx);
      if ((j++ % 64) == 0) {
        uint64_t g = tecm64_gcd(tecm64_mont_exit(P.Z, &ctx), n);
        if (g > 1 && g < n) RETURN_FROM_EACH_PRIME(return g);
        if (g == n)         RETURN_FROM_EACH_PRIME(goto next_curve);
      }
    } END_DO_FOR_EACH_PRIME
    {
      uint64_t g = tecm64_gcd(tecm64_mont_exit(P.Z, &ctx), n);
      if (g > 1 && g < n) return g;
    }

    if (B2 > B1) {
      uint64_t g = tinyecm64_stage2(&P, A24, &ctx, B1, B2);
      if (g > 1 && g < n) return g;
    }

next_curve:;
  }
  return 0;
}

int tinyecm64_factor(UV n, UV *factors, UV B1, UV B2, UV ncurves, UV sigma_offset)
{
  uint64_t f;
  if (B2 == 0) B2 = (B1 > UV_MAX/20) ? UV_MAX : 20*B1;
  if (ncurves      > NECM64_SIGMAS) ncurves      = NECM64_SIGMAS;
  if (sigma_offset > NECM64_SIGMAS) sigma_offset = NECM64_SIGMAS;
  f = tinyecm64((uint64_t)n, B1, B2, ncurves, sigma_offset);
  return (f > 1 && f < n) ? found_factor(n, (UV)f, factors)
                          : no_factor(n, factors);
}

#endif  /* HAS_ECM64 */


/* Pollard's P-1 */
int pminus1_factor(UV n, UV *factors, UV B1, UV B2)
{
  UV f, k, kmin;
  UV a     = 2, q     = 2;
  UV savea = 2, saveq = 2;
  UV j = 1;
  UV sqrtB1 = isqrt(B1);
#if USE_MONTMATH
  const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);
  UV ma = mont_geta(a,n);
  #define PMINUS1_APPLY_POWER  ma = mont_powmod(ma, k, n)
  #define PMINUS1_RECOVER_A    a = mont_recover(ma,n)
#else
  #define PMINUS1_APPLY_POWER  a = powmod(a, k, n)
  #define PMINUS1_RECOVER_A
#endif
  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in pminus1_factor");

  if (B1 <= primes_small[NPRIMES_SMALL-2]) {
    UV i;
    for (i = 1; primes_small[i] <= B1; i++) {
      q = k = primes_small[i];
      if (q <= sqrtB1) {
        k = q*q;  kmin = B1/q;
        while (k <= kmin)  k *= q;
      }
      PMINUS1_APPLY_POWER;
      if ( (j++ % 32) == 0) {
        PMINUS1_RECOVER_A;
        if (a == 0 || gcd_ui(a-1, n) != 1)
          break;
        savea = a;  saveq = q;
      }
    }
    PMINUS1_RECOVER_A;
  } else {
    START_DO_FOR_EACH_PRIME(2, B1) {
      q = k = p;
      if (q <= sqrtB1) {
        k = q*q;  kmin = B1/q;
        while (k <= kmin)  k *= q;
      }
      PMINUS1_APPLY_POWER;
      if ( (j++ % 32) == 0) {
        PMINUS1_RECOVER_A;
        if (a == 0 || gcd_ui(a-1, n) != 1)
          break;
        savea = a;  saveq = q;
      }
    } END_DO_FOR_EACH_PRIME
    PMINUS1_RECOVER_A;
  }
  if (a == 0) return no_factor(n,factors);
  f = gcd_ui(a-1, n);

  /* If we found more than one factor in stage 1, backup and single step */
  if (f == n) {
    a = savea;
    START_DO_FOR_EACH_PRIME(saveq, B1) {
      k = p;  kmin = B1/p;
      while (k <= kmin)  k *= p;
      a = powmod(a, k, n);
      f = gcd_ui(a-1, n);
      q = p;
      if (f != 1)
        break;
    } END_DO_FOR_EACH_PRIME
    /* If f == n again, we could do:
     * for (savea = 3; f == n && savea < 100; savea = next_prime(savea)) {
     *   a = savea;
     *   for (q = 2; q <= B1; q = next_prime(q)) {
     *     ...
     *   }
     * }
     * but this could be a huge time sink if B1 is large, so just fail.
     */
  }

  /* STAGE 2 */
  if (f == 1 && B2 > B1) {
    UV bm = a;
    UV b = 1;
    UV bmdiff;
    UV precomp_bm[111] = {0};    /* Enough for B2 = 189M */

    /* calculate (a^q)^2, (a^q)^4, etc. */
    bmdiff = sqrmod(bm, n);
    precomp_bm[0] = bmdiff;
    for (j = 1; j < 20; j++) {
      bmdiff = mulmod(bmdiff,bm,n);
      bmdiff = mulmod(bmdiff,bm,n);
      precomp_bm[j] = bmdiff;
    }

    a = powmod(a, q, n);
    j = 1;
    START_DO_FOR_EACH_PRIME( q+1, B2 ) {
      UV lastq = q;
      UV qdiff;
      q = p;
      /* compute a^q = a^lastq * a^(q-lastq) */
      qdiff = (q - lastq) / 2 - 1;
      if (qdiff >= 111) {
        bmdiff = powmod(bm, q-lastq, n);  /* Big gap */
      } else {
        bmdiff = precomp_bm[qdiff];
        if (bmdiff == 0) {
          if (precomp_bm[qdiff-1] != 0)
            bmdiff = mulmod(mulmod(precomp_bm[qdiff-1],bm,n),bm,n);
          else
            bmdiff = powmod(bm, q-lastq, n);
          precomp_bm[qdiff] = bmdiff;
        }
      }
      a = mulmod(a, bmdiff, n);
      if (a == 0) break;
      b = mulmod(b, a-1, n);   /* if b == 0, we found multiple factors */
      if ( (j++ % 64) == 0 ) {
        f = gcd_ui(b, n);
        if (f != 1)
          break;
      }
    } END_DO_FOR_EACH_PRIME
    f = gcd_ui(b, n);
  }
  return found_factor(n, f, factors);
}

/* Simple Williams p+1 */
static void pp1_pow(UV *cX, UV exp, UV n)
{
  UV X0 = *cX;
  UV X  = *cX;
  UV Y = mulsubmod(X, X, 2, n);
  UV bit = UVCONST(1) << (clz(exp)-1);
  while (bit) {
    UV T = mulsubmod(X, Y, X0, n);
    if ( exp & bit ) {
      X = T;
      Y = mulsubmod(Y, Y, 2, n);
    } else {
      Y = T;
      X = mulsubmod(X, X, 2, n);
    }
    bit >>= 1;
  }
  *cX = X;
}
int pplus1_factor(UV n, UV *factors, UV B1)
{
  UV X1, X2, f;
  UV sqrtB1 = isqrt(B1);
  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in pplus1_factor");

  X1 =  7 % n;
  X2 = 11 % n;
  f = 1;
  START_DO_FOR_EACH_PRIME(2, B1) {
    UV k = p;
    if (p < sqrtB1) {
      UV kmin = B1/p;
      while (k <= kmin)
        k *= p;
    }
    pp1_pow(&X1, k, n);
    if (X1 != 2) {
      f = gcd_ui( submod(X1, 2, n), n);
      if (f != 1 && f != n) break;
    }
    pp1_pow(&X2, k, n);
    if (X2 != 2) {
      f = gcd_ui( submod(X2, 2, n), n);
      if (f != 1 && f != n) break;
    }
  } END_DO_FOR_EACH_PRIME

  return found_factor(n, f, factors);
}


/* SQUFOF, based on Ben Buhrow's racing version. */
#if 1
  /* limit to 62-bit inputs, use 32-bit types, faster */
  #define SQUFOF_TYPE uint32_t
  #define SQUFOF_MAX  (UV_MAX >> 2)
#else
  /* All 64-bit inputs possible, though we severely limit multipliers */
  #define SQUFOF_TYPE UV
  #define SQUFOF_MAX  UV_MAX
#endif
typedef struct
{
  int valid;
  SQUFOF_TYPE P;
  SQUFOF_TYPE bn;
  SQUFOF_TYPE Qn;
  SQUFOF_TYPE Q0;
  SQUFOF_TYPE b0;
  SQUFOF_TYPE it;
  SQUFOF_TYPE imax;
  SQUFOF_TYPE mult;
} mult_t;

/* N < 2^63 (or 2^31).  Returns 0 or a factor */
static UV squfof_unit(UV n, mult_t* mult_save)
{
  SQUFOF_TYPE imax,i,Q0,Qn,bn,b0,P,bbn,Ro,S,So,t1,t2;
  uint32_t root;

  P = mult_save->P;
  bn = mult_save->bn;
  Qn = mult_save->Qn;
  Q0 = mult_save->Q0;
  b0 = mult_save->b0;
  i  = mult_save->it;
  imax = i + mult_save->imax;

#define SQUARE_SEARCH_ITERATION \
      t1 = P; \
      P = bn*Qn - P; \
      t2 = Qn; \
      Qn = Q0 + bn*(t1-P); \
      Q0 = t2; \
      bn = (b0 + P) / Qn; \
      i++;

  while (1) {
    int j = 0;
    if (i & 0x1) {
      SQUARE_SEARCH_ITERATION;
    }
    /* i is now even */
    while (1) {
      /* We need to know P, bn, Qn, Q0, iteration count, i  from prev */
      if (i >= imax) {
        /* save state and try another multiplier. */
        mult_save->P = P;
        mult_save->bn = bn;
        mult_save->Qn = Qn;
        mult_save->Q0 = Q0;
        mult_save->it = i;
        return 0;
      }

      SQUARE_SEARCH_ITERATION;

      /* Even iteration.  Check for square: Qn = S*S */
      if (is_perfect_square_ret(Qn,&root))
        break;

      /* Odd iteration. */
      SQUARE_SEARCH_ITERATION;
    }
    S = root; /* isqrt(Qn); */
    mult_save->it = i;

    /* Reduce to G0 */
    Ro = P + S*((b0 - P)/S);
    So = (n - (UV)Ro*(UV)Ro)/(UV)S;
    bbn = (b0+Ro)/So;

    /* Search for symmetry point */
#define SYMMETRY_POINT_ITERATION \
      t1 = Ro; \
      Ro = bbn*So - Ro; \
      t2 = So; \
      So = S + bbn*(t1-Ro); \
      S = t2; \
      bbn = (b0+Ro)/So; \
      if (Ro == t1) break;

    j = 0;
    while (1) {
      SYMMETRY_POINT_ITERATION;
      SYMMETRY_POINT_ITERATION;
      SYMMETRY_POINT_ITERATION;
      SYMMETRY_POINT_ITERATION;
      if (j++ > 2000000) {
         mult_save->valid = 0;
         return 0;
      }
    }

    t1 = gcd_ui(Ro, n);
    if (t1 > 1)
      return t1;
  }
}

/* Gower and Wagstaff 2008:
 *    http://www.ams.org/journals/mcom/2008-77-261/S0025-5718-07-02010-8/
 * Section 5.3.  I've added some with 13,17,19.  Sorted by F(). */
static const UV squfof_multipliers[] =
  /* { 3*5*7*11, 3*5*7, 3*5*11, 3*5, 3*7*11, 3*7, 5*7*11, 5*7,
       3*11,     3,     5*11,   5,   7*11,   7,   11,     1   }; */
  { 3*5*7*11, 3*5*7,  3*5*7*11*13, 3*5*7*13, 3*5*7*11*17, 3*5*11,
    3*5*7*17, 3*5,    3*5*7*11*19, 3*5*11*13,3*5*7*19,    3*5*7*13*17,
    3*5*13,   3*7*11, 3*7,         5*7*11,   3*7*13,      5*7,
    3*5*17,   5*7*13, 3*5*19,      3*11,     3*7*17,      3,
    3*11*13,  5*11,   3*7*19,      3*13,     5,           5*11*13,
    5*7*19,   5*13,   7*11,        7,        3*17,        7*13,
    11,       1 };
#define NSQUFOF_MULT (sizeof(squfof_multipliers)/sizeof(squfof_multipliers[0]))

int squfof_factor(UV n, UV *factors, UV rounds)
{
  mult_t mult_save[NSQUFOF_MULT];
  UV i, nn64, sqrtnn64, mult, f64,rounds_done = 0;
  int mults_racing = NSQUFOF_MULT;

  /* Caller should have handled these trivial cases */
  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in squfof_factor");

  /* Too big */
  if (n > SQUFOF_MAX)
    return no_factor(n,factors);

  for (i = 0; i < NSQUFOF_MULT; i++) {
    mult_save[i].valid = -1;
    mult_save[i].it = 0;
  }

  /* Race each multiplier for a bit (20-20k rounds) */
  while (mults_racing > 0 && rounds_done < rounds) {
    for (i = 0; i < NSQUFOF_MULT && rounds_done < rounds; i++) {
      if (mult_save[i].valid == 0)  continue;
      mult = squfof_multipliers[i];
      nn64 = n * mult;
      if (mult_save[i].valid == -1) {
        if ((SQUFOF_MAX / mult) < n) {
          mult_save[i].valid = 0; /* This multiplier would overflow 64-bit */
          mults_racing--;
          continue;
        }
        sqrtnn64 = isqrt(nn64);
        mult_save[i].valid = 1;
        mult_save[i].Q0    = 1;
        mult_save[i].b0    = sqrtnn64;
        mult_save[i].P     = sqrtnn64;
        mult_save[i].Qn    = (SQUFOF_TYPE)(nn64 - sqrtnn64 * sqrtnn64);
        if (mult_save[i].Qn == 0)
          return found_factor(n, sqrtnn64, factors);
        mult_save[i].bn    = (2 * sqrtnn64) / (UV)mult_save[i].Qn;
        mult_save[i].it    = 0;
        mult_save[i].mult  = mult;
        mult_save[i].imax  = (UV) (sqrt(sqrtnn64) / 16);
        if (mult_save[i].imax < 20)     mult_save[i].imax = 20;
        if (mult_save[i].imax > rounds) mult_save[i].imax = rounds;
      }
      if (mults_racing == 1)  /* Do all rounds if only one multiplier left */
        mult_save[i].imax = (rounds - rounds_done);
      f64 = squfof_unit(nn64, &mult_save[i]);
      if (f64 > 1) {
        UV f64red = f64 / gcd_ui(f64, mult);
        if (f64red > 1) {
          /* unsigned long totiter = 0;
             {int K; for (K = 0; K < NSQUFOF_MULT; K++) totiter += mult_save[K].it; }
             printf("  n %lu mult %lu it %lu (%lu)\n",n,mult,totiter,(UV)mult_save[i].it); */
          return found_factor(n, f64red, factors);
        }
        /* Found trivial factor.  Quit working with this multiplier. */
        mult_save[i].valid = 0;
      }
      if (mult_save[i].valid == 0)
        mults_racing--;
      rounds_done += mult_save[i].imax;   /* Assume we did all rounds */
    }
  }
  return no_factor(n,factors);
}

#define SQR_TAB_SIZE 512
static int sqr_tab_init = 0;
static double sqr_tab[SQR_TAB_SIZE];
static void make_sqr_tab(void) {
  int i;
  for (i = 0; i < SQR_TAB_SIZE; i++)
    sqr_tab[i] = sqrt((double)i);
  sqr_tab_init = 1;
}

/* Lehman written and tuned by Warren D. Smith.
 * Revised by Ben Buhrow and Dana Jacobsen. */
int lehman_factor(UV n, UV *factors, bool do_trial) {
  const double Tune = ((n >> 31) >> 5) ? 3.5 : 5.0;
  double x, sqrtn;
  UV a,c,kN,kN4,B2;
  uint32_t b,p,k,r,B,U,Bred,inc,ip=2;

  if (!(n&1)) return found_factor(n, 2, factors);

  B = Tune * (1+icbrt(n));

  if (do_trial) {
    uint32_t FirstCut = 0.1 * B;
    if (FirstCut < 84) FirstCut = 84;
    if (FirstCut > 65535) FirstCut = 65535;
    for (ip = 2;  ip < NPRIMES_SMALL; ip++) {
      p = primes_small[ip];
      if (p >= FirstCut)
        break;
      if (n % p == 0)
        return found_factor(n, p, factors);
    }
  }

#if BITS_PER_WORD == 64
  if (n >= UVCONST(8796393022207)) return no_factor(n,factors);
#endif
  Bred = B / (Tune * Tune * Tune);
  B2 = B*B;
  kN = 0;

  if (!sqr_tab_init) make_sqr_tab();
  sqrtn = sqrt(n);

  for (k = 1; k <= Bred; k++) {
    if (k&1) { inc = 4; r = (k+n) % 4; }
    else     { inc = 2; r = 1; }
    kN += n;
#if BITS_PER_WORD == 64
    if (kN >= UVCONST(1152921504606846976)) return no_factor(n,factors);
#endif
    kN4 = kN*4;

    x = (k < SQR_TAB_SIZE) ? sqrtn * sqr_tab[k] : sqrt((double)kN);
    a = x;
    if ((UV)a * (UV)a == kN)
      return found_factor(n, gcd_ui(a,n), factors);
    x *= 2;
    a = x + 0.9999999665;  /* Magic constant */
    b = a % inc;
    b = a + (inc+r-b) % inc;
    c = (UV)b*(UV)b - kN4;
    U = x + B2/(2*x);
    for (a = b;  a <= U;  c += inc*(a+a+inc), a += inc) {
      /* Check for perfect square */
      if (is_perfect_square_ret(c,&b)) {
        B2 = gcd_ui(a+b, n);
        return found_factor(n, B2, factors);
      }
    }
  }
  if (do_trial) {
    if (B > 65535) B = 65535;
    /* trial divide from primes[ip] to B.  We could:
     *   1) use table of 6542 shorts for the primes.
     *   2) use a wheel
     *   3) let trial_factor handle it
     */
    if (ip >= NPRIMES_SMALL)  ip = NPRIMES_SMALL-1;
    return trial_factor(n, factors, primes_small[ip], B);
  }
  return no_factor(n,factors);
}

/* Chebyshev polynomials of the first kind T_n(x) = V_n(2x,1) / 2. */
/* Basic algorithm from Daniel "Trizen" Șuteu */
int cheb_factor(UV n, UV *factors, UV B, UV initx)
{
  UV sqrtB, inv, x, f, i;

  if (B == 0) { B = log2floor(n);  B = 8*B*B; }
  if (B > isqrt(n)) B = isqrt(n);
  sqrtB = isqrt(B);
  inv = modinverse(2,n);   /* multiplying by this will divide by two */
  x = (initx == 0) ? 72 : initx;
  f = 1;

  START_DO_FOR_EACH_PRIME(2, B) {
    if (p <= sqrtB) {
      UV lgp = logint(B, p);
      UV plgp = ipowsafe(p, lgp);
      if (plgp < UV_MAX) {
        x = mulmod(lucasvmod(addmod(x,x,n), 1, plgp, n), inv, n);
      } else {
        for (i = 1; i <= lgp; i++)
          x = mulmod(lucasvmod(addmod(x,x,n), 1, p, n), inv, n);
      }
    } else {
      x = mulmod(lucasvmod(addmod(x,x,n), 1, p, n), inv, n);
    }
    f = gcd_ui(x-1, n);  if (f > 1)  break;
  } END_DO_FOR_EACH_PRIME

  if (f > 1 && f < n)
    return found_factor(n, f, factors);
  return no_factor(n,factors);
}


/******************************************************************************/
/******************************************************************************/


static const uint32_t _fr_chunk = 256*1024;

/* Help performance by doing a cube root sieve for small ranges */
static bool _fr_full_sieve(UV sqrtn, UV range)  /* range = hi-lo */
{
  if (sqrtn <   10000000U) return 1;               /* Below 10^14 */
  if (sqrtn <   35000000U) return (range >   900); /* Below 10^15 */
  if (sqrtn <  100000000U) return (range >  1700); /* Below 10^16 */
  if (sqrtn <  350000000U) return (range >  3400); /* Below 10^17 */
  if (sqrtn < 1000000000U) return (range >  5500); /* Below 10^18 */
  if (sqrtn < 3500000000U) return (range > 17000); /* Below 10^19 */
  return (range > 19000);
}

static void _vec_factor(UV lo, UV hi, UV *nfactors, UV *farray, UV noffset, bool square_free)
{
  UV *N, j, n, sqrthi, sievelim;
  sqrthi = isqrt(hi);
  n = hi-lo+1;
  New(0, N, hi-lo+1, UV);
  for (j = 0; j < n; j++) {
    N[j] = 1;
    nfactors[j] = 0;
  }
  sievelim = _fr_full_sieve(sqrthi, hi-lo)  ?  sqrthi  :  icbrt(hi);
  START_DO_FOR_EACH_PRIME(2, sievelim) {
    UV q, t, A;
    if (square_free == 0) {
      UV kmin = hi / p;
      for (q = p; q <= kmin; q *= p) {
        t = lo / q, A = t * q;
        if (A < lo) A += q;
        for (j = A-lo; j < n; j += q) {
          farray[ j*noffset + nfactors[j]++ ] = p;
          N[j] *= p;
        }
      }
    } else {
      q = p*p, t = lo / q, A = t * q;
      if (A < lo) A += q;
      for (j = A-lo; j < n; j += q) {
        N[j] = 0;
        nfactors[j] = 0;
      }
      q = p, t = lo / q, A = t * q;
      if (A < lo) A += q;
      for (j = A-lo; j < n; j += q) {
        if (N[j] > 0) {
          farray[ j*noffset + nfactors[j]++ ] = p;
          N[j] *= p;
        }
      }
    }
  } END_DO_FOR_EACH_PRIME

  if (sievelim == sqrthi) {
    /* Handle the unsieved results, which are prime */
    for (j = 0; j < n; j++) {
      if (N[j] == 1)
        farray[ j*noffset + nfactors[j]++ ] = j+lo;
      else if (N[j] > 0 && N[j] != j+lo)
        farray[ j*noffset + nfactors[j]++ ] = (j+lo) / N[j];
    }
  } else {
    /* Handle the unsieved results, which are prime or semi-prime */
    for (j = 0; j < n; j++) {
      UV rem = j+lo;
      if (N[j] > 0 && N[j] != rem) {
        if (N[j] != 1)
          rem /= N[j];
        if (square_free && is_perfect_square(rem)) {
          nfactors[j] = 0;
        } else {
          UV* f = farray + j*noffset + nfactors[j];
          nfactors[j] += factor_one(rem, f, 1, 0);
        }
      }
    }
  }
  Safefree(N);
}

factor_range_context_t factor_range_init(UV lo, UV hi, bool square_free) {
  factor_range_context_t ctx;
  UV span = hi-lo;

  MPUassert(hi >= lo, "factor_range_init: hi < lo");

  ctx.lo = lo;
  ctx.hi = hi;
  ctx.n = lo-1;
  ctx.is_square_free = square_free;

  if (span >= 100) {          /* Sieve in chunks */
    if (square_free) ctx._noffset = (hi <= 42949672965UL) ? 10 : 15;
    else             ctx._noffset = BITS_PER_WORD - clz(hi);
    ctx._coffset = _fr_chunk;
    New(0, ctx._nfactors, _fr_chunk, UV);
    New(0, ctx._farray, _fr_chunk * ctx._noffset, UV);
    { /* Prealloc all the sieving primes now. */
      UV t = isqrt(hi);
      if (!_fr_full_sieve(t, span))  t = icbrt(hi);
      get_prime_cache(t, 0);
    }
  } else {                    /* factor each number */
    New(0, ctx.factors, square_free ? 15 : 63, UV);
    ctx._nfactors = 0;
    ctx._farray = ctx.factors;
    ctx._noffset = 0;
  }
  return ctx;
}

int factor_range_next(factor_range_context_t *ctx) {
  int j, nfactors;
  UV n;
  if (ctx->n >= ctx->hi)
    return -1;
  n = ++(ctx->n);
  if (ctx->_nfactors) {
    if (ctx->_coffset >= _fr_chunk) {
      UV clo = n;
      UV chi = n + _fr_chunk - 1;
      if (chi > ctx->hi || chi < clo) chi = ctx->hi;
      _vec_factor(clo, chi, ctx->_nfactors, ctx->_farray, ctx->_noffset, ctx->is_square_free);
      ctx->_coffset = 0;
    }
    nfactors = ctx->_nfactors[ctx->_coffset];
    ctx->factors = ctx->_farray + ctx->_coffset * ctx->_noffset;
    ctx->_coffset++;
  } else {
    if (ctx->is_square_free && n >= 49 && (!(n% 4) || !(n% 9) || !(n%25) || !(n%49)))
      return 0;
    nfactors = factor(n, ctx->factors);
    if (ctx->is_square_free) {
      for (j = 1; j < nfactors; j++)
        if (ctx->factors[j] == ctx->factors[j-1])
          break;
      if (j < nfactors) return 0;
    }
  }
  return nfactors;
}

void factor_range_destroy(factor_range_context_t *ctx) {
  if (ctx->_farray != 0) Safefree(ctx->_farray);
  if (ctx->_nfactors != 0) Safefree(ctx->_nfactors);
  ctx->_farray = ctx->_nfactors = ctx->factors = 0;
}


/******************************************************************************/
/* Find number of factors for all values in a range */
/******************************************************************************/

unsigned char* range_nfactor_sieve(UV lo, UV hi, bool with_multiplicity) {
  unsigned char* nf;
  UV *N, i, range, sqrtn = isqrt(hi);

  if (hi < lo || hi-lo == UV_MAX) croak("Invalid range in range_nfactor_sieve");
  range = hi - lo + 1;
  Newz(0, nf, range, unsigned char);
  New(0, N, range, UV);

  /* We could set to 1 and sieve from 2, or do this initialization */
  for (i = 0; i < range; i++) {
    UV v = lo+i;
    N[i] = 1;
    if (!(v&1) && v >= 2) {
      UV k = v >> 1;
      unsigned char nz = 1;
      while (!(k&1)) { nz++; k >>= 1; }
      nf[i] = (with_multiplicity) ? nz : 1;
      N[i] = UVCONST(1) << nz;
    }
  }

  START_DO_FOR_EACH_PRIME(3, sqrtn) {
    UV pk, maxpk = UV_MAX/p; \
    for (i = P_GT_LO_0(p,p,lo); i < range; i += p)
      { N[i] *= p;  nf[i]++; }
    for (pk = p*p; pk <= hi; pk *= p) {
      for (i = P_GT_LO_0(pk,pk,lo); i < range; i += pk)
        { N[i] *= p;  if (with_multiplicity) nf[i]++; }
      if (pk >= maxpk) break;  /* Overflow protection */
    }
  } END_DO_FOR_EACH_PRIME

  for (i = 0; i < range; i++)
    if (N[i] < (lo+i))
      nf[i]++;
  Safefree(N);
  if (lo == 0) nf[0] = 1;
  return nf;
}


UV sopfr(UV n)
{
  factored_t nf;
  UV sum;
  uint32_t i;

  if (n <= 5) return n - (n==1);

  nf = factorint(n);
  sum = 0;
  for (i = 0; i < nf.nfactors; i++)
    sum += (nf.f[i] * (UV)nf.e[i]);
  return sum;
}

UV sopf(UV n)
{
  factored_t nf;
  UV sum;
  uint32_t i;

  if (n <= 3) return n - (n==1);

  nf = factorint(n);
  sum = 0;
  for (i = 0; i < nf.nfactors; i++)
    sum += nf.f[i];
  return sum;
}

/******************************************************************************/

/* Compile with:
 *  gcc -O3 -fomit-frame-pointer -march=native -Wall -DSTANDALONE -DFACTOR_STANDALONE factor.c util.c primality.c cache.c sieve.c chacha.c csprng.c prime_counts.c prime_count_cache.c lmo.c legendre_phi.c real.c inverse_interpolate.c rootmod.c lucas_seq.c prime_powers.c sort.c -lm
 */
#ifdef FACTOR_STANDALONE
#include <errno.h>
int main(int argc, char *argv[])
{
  UV n;
  UV factors[MPU_MAX_FACTORS+1];
  int nfactors, i, a;

  if (argc <= 1) {
    char line[1024];
    while (1) {
      if (!fgets(line,sizeof(line),stdin)) break;
      n = strtoull(line, 0, 10);
      nfactors = factor(n, factors);
      if (nfactors == 1) {
        printf("%"UVuf": %"UVuf"\n",n,n);
      } else if (nfactors == 2) {
        printf("%"UVuf": %"UVuf" %"UVuf"\n",n,factors[0],factors[1]);
      } else if (nfactors == 3) {
        printf("%"UVuf": %"UVuf" %"UVuf" %"UVuf"\n",n,factors[0],factors[1],factors[2]);
      } else {
        printf("%"UVuf": %"UVuf" %"UVuf" %"UVuf" %"UVuf"",n,factors[0],factors[1],factors[2],factors[3]);
        for (i = 4; i < nfactors; i++) printf(" %"UVuf"", factors[i]);
        printf("\n");
      }
    }
    exit(0);
  }

  for (a = 1; a < argc; a++) {
    n = strtoul(argv[a], 0, 10);
    if (n == ULONG_MAX && errno == ERANGE) { printf("Argument larger than ULONG_MAX\n"); return(-1); }
    nfactors = factor(n, factors);
    printf("%"UVuf":", n);
    for (i = 0; i < nfactors; i++)
      printf(" %"UVuf"", factors[i]);
    printf("\n");
  }

  return(0);
}
#endif
