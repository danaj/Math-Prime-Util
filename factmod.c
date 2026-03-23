#define FUNC_isqrt 1
#define FUNC_ipow 1
#include "ptypes.h"
#include "sieve.h"
#include "util.h"
#include "factor.h"
#include "mulmod.h"
#include "montmath.h"
#include "constants.h"
#include "factmod.h"


/******************************************************************************/
/*                                  N! MOD M                                  */
/******************************************************************************/

static UV _powersin(UV p, UV d) {
  UV td = d/p, e = td;
  do { td/=p; e += td; } while (td > 0);
  return e;
}

static UV _facmod(UV n, UV m) {
  UV i, res = 1;

  if (n < 1000) {

    for (i = 2; i <= n && res != 0; i++)
      res = mulmod(res,i,m);

  } else {

    unsigned char* segment;
    UV seg_base, seg_low, seg_high;
    UV sqn = isqrt(n), nsqn = n/sqn, j = sqn, nlo = 0, nhi = 0, s1 = 1;
    void* ctx = start_segment_primes(7, n, &segment);

    for (i = 1; i <= 3; i++) {  /* Handle 2,3,5 assume n>=25*/
      UV p = primes_small[i];
      res = mulmod(res, powmod(p,_powersin(p, n),m), m);
    }
    while (res!=0 && next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
        if (p <= nsqn) {
          res = mulmod(res, powmod(p,_powersin(p,n),m), m);
        } else {
          while (p > nhi) {
            res = mulmod(res, powmod(s1,j,m), m);
            s1 = 1;
            j--;
            nlo = n/(j+1)+1;
            nhi = n/j;
          }
          if (p >= nlo)
            s1 = mulmod(s1, p, m);
        }
        if (res == 0) break;
      END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
    res = mulmod(res, s1, m);

  }

  return res;
}
#if USE_MONTMATH
static UV _facmod_mont(UV n, UV m) {
  const uint64_t npi = mont_inverse(m),  mont1 = mont_get1(m);
  uint64_t monti = mont1;
  UV i, res = mont1;

  if (n < 1000) {

    for (i = 2; i <= n && res != 0; i++) {
      monti = addmod(monti,mont1,m);
      res = mont_mulmod(res,monti,m);
    }

  } else {

    unsigned char* segment;
    UV seg_base, seg_low, seg_high;
    UV sqn = isqrt(n), nsqn = n/sqn, j = sqn, nlo = 0, nhi = 0;
    UV s1 = mont1;
    void* ctx = start_segment_primes(7, n, &segment);

    for (i = 1; i <= 3; i++) {  /* Handle 2,3,5 assume n>=25*/
      UV p = primes_small[i];
      UV mp = mont_geta(p,m);
      res = mont_mulmod(res, mont_powmod(mp,_powersin(p,n),m), m);
    }
    while (res!=0 && next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
        UV mp = mont_geta(p,m);
        if (p <= nsqn) {
          res = mont_mulmod(res, mont_powmod(mp,_powersin(p,n),m), m);
        } else {
          while (p > nhi) {
            res = mont_mulmod(res, mont_powmod(s1,j,m), m);
            s1 = mont1;
            j--;
            nlo = n/(j+1)+1;
            nhi = n/j;
          }
          if (p >= nlo)
            s1 = mont_mulmod(s1, mp, m);
        }
        if (res == 0) break;
      END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
    res = mont_mulmod(res, s1, m);

  }

  res = mont_recover(res, m);
  return res;
}
#endif

UV factorialmod(UV n, UV m) {  /*  n! mod m */
  UV d = n, res = 1;
  uint32_t i;
  bool m_prime;

  if (n >= m || m == 1) return 0;
  if (n <= 1 || m == 2) return (n <= 1);

  if (n <= 10) { /* Keep things simple for small n */
    for (i = 2; i <= n && res != 0; i++)
      res = (res * i) % m;
    return res;
  }

  m_prime = is_prime(m);
  if (n > m/2 && m_prime)    /* Check if we can go backwards */
    d = m-n-1;
  if (d < 2)
    return (d == 0) ? m-1 : 1;   /* Wilson's Theorem: n = m-1 and n = m-2 */

  if (d > 100 && !m_prime) {   /* Check for composite m that leads to 0 */
    factored_t mf = factorint(m);
    UV maxpk = 0;
    for (i = 0; i < mf.nfactors; i++) {
      UV t = mf.f[i] * mf.e[i];   /* Possibly too high if exp[j] > fac[j] */
      if (t > maxpk)
        maxpk = t;
    }
    /* Maxpk is >= S(m), the Kempner number A002034 */
    if (n >= maxpk)
      return 0;
  }

#if USE_MONTMATH
  if (m & 1) {
    res = _facmod_mont(d, m);
  } else
#endif
  {
    res = _facmod(d, m);
  }

  if (d != n && res != 0) {      /* Handle backwards case */
    if (!(d&1)) res = submod(m,res,m);
    res = modinverse(res,m);
  }

  return res;
}


/******************************************************************************/
/*                          BINOMIAL(N,K) MOD M                               */
/******************************************************************************/

static UV _factorial_valuation(UV n, UV p) {
  UV k = 0;
  while (n >= p) {
    n /= p;
    k += n;
  }
  return k;
}
static int _binoval(UV n, UV k, UV m) {
  return _factorial_valuation(n,m) - _factorial_valuation(k,m) - _factorial_valuation(n-k,m);
}
static UV _factorialmod_without_prime(UV n, UV p, UV m) {
  UV i, pmod, r = 1;
  MPUassert(p >= 2 && m >= p && (m % p) == 0, "_factorialmod called with wrong args");
  if (n <= 1) return 1;

  if (n >= m) {
    /* Note with p=2 the behaviour is different */
    if ( ((n/m) & 1) && (p > 2 || m == 4) )  r = m-1;
    n %= m;
  }

#if USE_MONTMATH
  if (m & 1) {
    const uint64_t npi = mont_inverse(m),  mont1 = mont_get1(m);
    uint64_t mi = mont1;
    r = mont_geta(r, m);
    for (i = pmod = 2; i <= n; i++) {
      mi = addmod(mi, mont1, m);
      if (pmod++ == p) pmod = 1;
      else             r = mont_mulmod(r, mi, m);
    }
    r = mont_recover(r, m);
  } else
#endif
  {
    for (i = pmod = 2; i <= n; i++) {
      if (pmod++ == p) pmod = 1;
      else             r = mulmod(r, i, m);
    }
  }
  return r;
}
static UV _factorialmod_without_prime_powers(UV n, UV p, UV m) {
  UV ip, r = 1;

  for (ip = n; ip > 1; ip /= p)
    r = mulmod(r, _factorialmod_without_prime(ip, p, m), m);

  return r;
}
static UV _binomial_mod_prime_power(UV n, UV k, UV p, UV e) {
  UV r, b, m, i, num, den, ip, ires;

  if (k > n) return 0;
  if (k == 0 || k == n) return 1;
  if (k > n/2)  k = n-k;

  b = _binoval(n,k,p);
  if (e <= b) return 0;
  m = ipow(p,e);

  if (k == 1) return n % m;

  /* Both methods work fine -- choose based on performance. */
  den  = _factorialmod_without_prime_powers(k, p, m);
  if (k >= m) {
    num  = _factorialmod_without_prime_powers(n, p, m);
    ip   = _factorialmod_without_prime_powers(n-k, p, m);
    den = mulmod(den, ip, m);
  } else {
#if USE_MONTMATH
    if (m & 1) {
      const uint64_t npi = mont_inverse(m),  mont1 = mont_get1(m);
      num = mont1;
      for (i = n-k+1, ires = (i-1)%p; i <= n; i++) {
        ip = i;
        if (++ires == p) { ires = 0; do { ip /= p; } while ((ip % p) == 0); }
        num = mont_mulmod(num, mont_geta(ip, m), m);
      }
      num = mont_recover(num, m);
    } else
#endif
    {
      num = 1;
      for (i = n-k+1, ires = (i-1) % p; i <= n; i++) {
        ip = i;
        if (++ires == p) { ires = 0; do { ip /= p; } while ((ip % p) == 0); }
        num = mulmod(num, ip, m);
      }
    }
  }

  r = divmod(num, den, m);
  if (b > 0) r = mulmod(r, ipow(p,b), m);
  return r;
}

static UV _binomial_lucas_mod_prime(UV n, UV k, UV p) {
  UV res, t, vn[BITS_PER_WORD], vk[BITS_PER_WORD];
  int i, ln, lk;

  if (p < 2) return 0;
  if (p == 2) return !(~n & k);

  for (t = n, ln = 0; t > 0; t /= p)
    vn[ln++] = t % p;
  for (t = k, lk = 0; t > 0; t /= p)
    vk[lk++] = t % p;

  res = 1;
  for (i = ln-1; i >= 0; i--) {
    UV ni = vn[i];
    UV ki = (i < lk) ? vk[i] : 0;
    res = mulmod(res, _binomial_mod_prime_power(ni, ki, p, 1), p);
  }
  return res;
}

/* Based on Granville's paper on the generalization of Lucas's theorem to
 * prime powers: https://www.dms.umontreal.ca/~andrew/Binomial/genlucas.html
 * and Max Alekseyev's binomod.gp program. */
static UV _binomial_lucas_mod_prime_power(UV n, UV k, UV p, UV q) {
  UV N[BITS_PER_WORD], K[BITS_PER_WORD], R[BITS_PER_WORD], e[BITS_PER_WORD];
  UV i, d, m, n1, k1, r1, m1, res;

  MPUassert(q < BITS_PER_WORD, "bad exponent in binomialmod generalized lucas");
  m = ipow(p, q);

  /* Construct the digits for N, K, and N-K (R). */
  n1 = n;   k1 = k;  r1 = n-k;
  for (d = 0; n1 > 0; d++) {
    N[d] = n1 % p;  n1 /= p;
    K[d] = k1 % p;  k1 /= p;
    R[d] = r1 % p;  r1 /= p;
  }
  /* Compute the number of carries. */
  for (i = 0; i < d; i++)
    e[i] = (N[i] < (K[i] + ((i > 0) ? e[i-1] : 0)));
  /* Turn the carries into a cumulative count. */
  for (i = d-1; i >= 1; i--)
    e[i-1] += e[i];

  if (e[0] >= q) return 0;
  q -= e[0];
  m1 = ipow(p, q);

  /* Now make the digits for the reduced N, K, N-K */
  n1 = n;   k1 = k;  r1 = n-k;
  for (d = 0; n1 > 0; d++) {
    N[d] = n1 % m1;  n1 /= p;
    K[d] = k1 % m1;  k1 /= p;
    R[d] = r1 % m1;  r1 /= p;
  }

  /* Theorem 1 from Granville indicates the +/- 1.  */
  res = ((p > 2 || q < 3) && q < d && e[q-1] % 2)  ?  m-1  :  1;
  res = mulmod(res, powmod(p, e[0], m), m);

  /* Compute the individual binomials (again, theorem 1) */
  for (i = 0; i < d; i++) {
    UV ni = _factorialmod_without_prime(N[i], p, m);
    UV ki = _factorialmod_without_prime(K[i], p, m);
    UV ri = _factorialmod_without_prime(R[i], p, m);
    UV r = divmod(ni, mulmod(ki, ri, m), m);
    res = mulmod(res, r, m);
  }
  return res;
}

bool binomialmod(UV *res, UV n, UV k, UV m) {

  if (m <= 1)           { *res = 0; return 1; }
  if (k == 0 || k >= n) { *res = (k == 0 || k == n); return 1; }

  if (m == 2) { *res = !(~n & k); return 1; }

#if 0
    if ( (*res = binomial(n,k)) )
      { *res %= m; return 1; }
#endif

  if (is_prime(m)) {
    *res = _binomial_lucas_mod_prime(n, k, m);
    return 1;
  }
  {
    UV bin[MPU_MAX_DFACTORS], mod[MPU_MAX_DFACTORS];
    uint32_t i;
    factored_t mf = factorint(m);

    for (i = 0; i < mf.nfactors; i++) {
      if (mf.e[i] == 1) {
        bin[i] = _binomial_lucas_mod_prime(n, k, mf.f[i]);
        mod[i] = mf.f[i];
      } else {
        /* bin[i] = _binomial_mod_prime_power(n, k, mf.f[i], mf.e[i]); */
        /* Use generalized Lucas */
        bin[i] = _binomial_lucas_mod_prime_power(n, k, mf.f[i], mf.e[i]);
        mod[i] = ipow(mf.f[i], mf.e[i]);
      }
    }
    /* chinese with p^e as modulos, so should never get -1 back */
    return chinese(res, 0, bin, mod, mf.nfactors) == 1;
  }
}
