/******************************************************************************/
/*                              MODULAR ROOTS                                 */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ptypes.h"
#define FUNC_isqrt 1
#define FUNC_is_perfect_square 1
#define FUNC_gcd_ui 1
#define FUNC_ipow 1
#include "util.h"
#include "mulmod.h"
#include "factor.h"
#include "rootmod.h"

/* Pick one or both */
#define USE_ROOTMOD_SPLITK 1    /* enables rootmod_composite1 */
#define USE_ROOTMOD_SPLITN 1    /* enables rootmod_composite2 */


/******************************************************************************/
/*                               SQRT(N) MOD M                                */
/******************************************************************************/

/* _sqrtmod_prime       assumes  1 < a < p,  n > 1,  p > 2,  p prime.
 * _sqrtmod_prime_power assumes  1 < a < p,  n > 1,  p > 2,  p prime.
 * If any of these are not true, the result is undefined.
 *
 * _sqrtmod_composite takes care of the edge conditions and factors n.
 *
 * _sqrtmod_composite and _sqrtmod_prime_power always return UV_MAX
 * if no root exists, while any other return value will be a valid root.
 *
 * The exported functions sqrtmod(a,n) and rootmod(a,2,n) further:
 *   - verify the result and return success / fail in a separate int.
 *   - always returns the smaller of the two roots.
 *
 * sqrtmodp / rootmodp does the same except n is assumed prime.
 */

#if !USE_MONTMATH
static UV _sqrtmod_prime(UV a, UV p) {
  if ((p % 4) == 3) {
    return powmod(a, (p+1)>>2, p);
  }
  if ((p % 8) == 5) { /* Atkin's algorithm.  Faster than Legendre. */
    UV a2, alpha, beta, b;
    a2 = addmod(a,a,p);
    alpha = powmod(a2,(p-5)>>3,p);
    beta  = mulmod(a2,sqrmod(alpha,p),p);
    b     = mulmod(alpha, mulmod(a, (beta ? beta-1 : p-1), p), p);
    return b;
  }
  if ((p % 16) == 9) { /* Müller's algorithm extending Atkin */
    UV a2, alpha, beta, b, d = 1;
    a2 = addmod(a,a,p);
    alpha = powmod(a2, (p-9)>>4, p);
    beta  = mulmod(a2, sqrmod(alpha,p), p);
    if (sqrmod(beta,p) != p-1) {
      do { d += 2; } while (kronecker_uu(d,p) != -1 && d < p);
      alpha = mulmod(alpha, powmod(d,(p-9)>>3,p), p);
      beta  = mulmod(a2, mulmod(sqrmod(d,p),sqrmod(alpha,p),p), p);
    }
    b = mulmod(alpha, mulmod(a, mulmod(d,(beta ? beta-1 : p-1),p),p),p);
    return b;
  }

  /* Verify Euler condition for odd p */
  if ((p & 1) && powmod(a,(p-1)>>1,p) != 1) return 0;

  /* Algorithm 1.5.1 from Cohen.  Tonelli/Shanks. */
  {
    UV x, q, e, t, z, r, m, b;
    q = p-1;
    e = valuation_remainder(q, 2, &q);
    t = 3;
    while (kronecker_uu(t, p) != -1) {
      t += 2;
      if (t == 201) {           /* exit if p looks like a composite */
        if ((p % 2) == 0 || powmod(2, p-1, p) != 1 || powmod(3, p-1, p) != 1)
          return 0;
      } else if (t >= 20000) {  /* should never happen */
        return 0;
      }
    }
    z = powmod(t, q, p);
    b = powmod(a, q, p);
    r = e;
    q = (q+1) >> 1;
    x = powmod(a, q, p);
    while (b != 1) {
      t = b;
      for (m = 0; m < r && t != 1; m++)
        t = sqrmod(t, p);
      if (m >= r) break;
      t = powmod(z, UVCONST(1) << (r-m-1), p);
      x = mulmod(x, t, p);
      z = mulmod(t, t, p);
      b = mulmod(b, z, p);
      r = m;
    }
    return x;
  }
  return 0;
}
#else
static UV _sqrtmod_prime(UV a, UV p) {
  const uint64_t npi = mont_inverse(p),  mont1 = mont_get1(p);
  a = mont_geta(a,p);

  if ((p % 4) == 3) {
    UV b = mont_powmod(a, (p+1)>>2, p);
    return mont_recover(b, p);
  }

  if ((p % 8) == 5) { /* Atkin's algorithm.  Faster than Legendre. */
    UV a2, alpha, beta, b;
    a2 = addmod(a,a,p);
    alpha = mont_powmod(a2,(p-5)>>3,p);
    beta  = mont_mulmod(a2,mont_sqrmod(alpha,p),p);
    beta  = submod(beta, mont1, p);
    b     = mont_mulmod(alpha, mont_mulmod(a, beta, p), p);
    return mont_recover(b, p);
  }
  if ((p % 16) == 9) { /* Müller's algorithm extending Atkin */
    UV a2, alpha, beta, b, d = 1;
    a2 = addmod(a,a,p);
    alpha = mont_powmod(a2, (p-9)>>4, p);
    beta  = mont_mulmod(a2, mont_sqrmod(alpha,p), p);
    if (mont_sqrmod(beta,p) != submod(0,mont1,p)) {
      do { d += 2; } while (kronecker_uu(d,p) != -1 && d < p);
      d = mont_geta(d,p);
      alpha = mont_mulmod(alpha, mont_powmod(d,(p-9)>>3,p), p);
      beta  = mont_mulmod(a2, mont_mulmod(mont_sqrmod(d,p),mont_sqrmod(alpha,p),p), p);
      beta  = mont_mulmod(submod(beta,mont1,p), d, p);
    } else {
      beta  = submod(beta, mont1, p);
    }
    b = mont_mulmod(alpha, mont_mulmod(a, beta, p), p);
    return mont_recover(b, p);
  }

  /* Verify Euler condition for odd p */
  if ((p & 1) && mont_powmod(a,(p-1)>>1,p) != mont1) return 0;

  /* Algorithm 1.5.1 from Cohen.  Tonelli/Shanks. */
  {
    UV x, q, e, t, z, r, m, b;
    q = p-1;
    e = valuation_remainder(q, 2, &q);
    t = 3;
    while (kronecker_uu(t, p) != -1) {
      t += 2;
      if (t == 201) {           /* exit if p looks like a composite */
        if ((p % 2) == 0 || powmod(2, p-1, p) != 1 || powmod(3, p-1, p) != 1)
          return 0;
      } else if (t >= 20000) {  /* should never happen */
        return 0;
      }
    }
    t = mont_geta(t, p);
    z = mont_powmod(t, q, p);
    b = mont_powmod(a, q, p);
    r = e;
    q = (q+1) >> 1;
    x = mont_powmod(a, q, p);
    while (b != mont1) {
      t = b;
      for (m = 0; m < r && t != mont1; m++)
        t = mont_sqrmod(t, p);
      if (m >= r) break;
      t = mont_powmod(z, UVCONST(1) << (r-m-1), p);
      x = mont_mulmod(x, t, p);
      z = mont_mulmod(t, t, p);
      b = mont_mulmod(b, z, p);
      r = m;
    }
    return mont_recover(x, p);
  }
  return 0;
}
#endif

static UV _sqrtmod_prime_power(UV a, UV p, UV e) {
  UV r, s, n, pk, apk, ered, np;

  if (e == 1) {
    if (a >= p) a %= p;
    if (p == 2 || a == 0)  return a;
    r = _sqrtmod_prime(a,p);
    return (sqrmod(r,p) == a) ? r : UV_MAX;
  }

  n  = ipow(p,e);
  pk = p*p;

  if ((a % n) == 0)
    return 0;

  if ((a % pk) == 0) {
    apk = a / pk;
    s = _sqrtmod_prime_power(apk, p, e-2);
    if (s == UV_MAX) return UV_MAX;
    return s * p;
  }

  if ((a % p) == 0)
    return UV_MAX;

  ered = (p > 2 || e < 5)  ?  (e+1)>>1  :  (e+3)>>1;
  s = _sqrtmod_prime_power(a, p, ered);
  if (s == UV_MAX) return UV_MAX;

  np = (p != 2 || (n > (UV_MAX/p)))  ?  n  :  n * p;
  r = addmod(s, gcddivmod(submod(a,sqrmod(s,np),np), addmod(s,s,np), n), n);
  if (sqrmod(r,n) != (a % n)) return UV_MAX;
  return r;
}

static UV _sqrtmod_composite(UV a, UV n) {
  UV r, s, t, i, fe, N, inv;
  UV fac[MPU_MAX_FACTORS+1];
  UV exp[MPU_MAX_FACTORS+1];
  int nfactors;

  if (n == 0) return UV_MAX;
  if (a >= n) a %= n;
  if (n <= 2 || a <= 1)      return a;
  if (is_perfect_square(a))  return isqrt(a);

  nfactors = factor_exp(n, fac, exp);
  N = ipow(fac[0], exp[0]);
  r = _sqrtmod_prime_power(a, fac[0], exp[0]);
  if (r == UV_MAX) return UV_MAX;
  for (i = 1; i < (UV) nfactors; i++) {
    fe = ipow(fac[i], exp[i]);
    s = _sqrtmod_prime_power(a, fac[i], exp[i]);
    if (s == UV_MAX) return UV_MAX;
    inv = modinverse(N, fe);
    t = mulmod(inv, submod(s % fe,r % fe,fe), fe);
    r = addmod(r, mulmod(N,t,n), n);
    N *= fe;
  }
  return r;
}

static int _sqrtmod_return(UV r, UV *s, UV a, UV p) {
  if (p-r < r)  r = p-r;
  if (mulmod(r, r, p) != (a % p)) return 0;
  if (s != 0) *s = r;
  return 1;
}
int sqrtmodp(UV *s, UV a, UV p) {
  if (p == 0) return 0;
  if (a >= p) a %= p;
  if (p <= 2 || a <= 1) return _sqrtmod_return(a, s, a, p);
  return _sqrtmod_return(_sqrtmod_prime(a,p), s, a, p);
}

int sqrtmod(UV *s, UV a, UV n) {
  /* return rootmod(s, a, 2, n); */
  return _sqrtmod_return(_sqrtmod_composite(a,n), s, a, n);
}




/******************************************************************************/
/*                          K-TH ROOT OF N MOD M                              */
/******************************************************************************/

static int _rootmod_return(UV r, UV *s, UV a, UV k, UV p) {
  if (k == 2 && p-r < r)  r = p-r;
  if (powmod(r, k, p) != (a % p)) return 0;
  if (s != 0) *s = r;
  return 1;
}


/* Generalized Tonelli-Shanks for k-th root mod a prime, with k prime */
static UV _ts_prime(UV a, UV k, UV p, UV *z) {
  UV A, B, y, x, e, r, T, ke, t;

  /* Assume:  k > 1,  1 < a < p,  p > 2,  k prime,  p prime */

  for (e = 0, r = p-1; !(r % k); r /= k) e++;
  /*  p-1 = r * k^e   =>   ke = ipow(k,e) = (p-1)/r  */
  ke = (p-1)/r;

  x = powmod(a, modinverse(k % r, r), p);
  B = mulmod(powmod(x, k, p), modinverse(a, p), p);

  for (T = 2, y = 1;  y == 1;  T++) {
    t = powmod(T, r, p);
    y = powmod(t, ke/k, p);
  }

  while (ke != k) {
    ke = ke/k;
    T = t;
    t = powmod(t, k, p);
    A = powmod(B, ke/k, p);
    while (A != 1) {
      x = mulmod(x, T, p);
      B = mulmod(B, t, p);
      A = mulmod(A, y, p);
    }
  }
  if (z)  *z = t;
  return x;
}

#if USE_ROOTMOD_SPLITK
/* Alternate, taking prime p but composite k. */
/* k-th root using Tonelli-Shanks for prime k and p */
/* This works much better for me than AMM (Holt 2003 or Cao/Sha/Fan 2011). */
/* See Algorithm 3.3 of van de Woestijne (2006). */
/* https://www.opt.math.tugraz.at/~cvdwoest/maths/dissertatie.pdf */
/* Also see Pari's Tonelli-Shanks by Bill Allombert, 2014,2017, which seems */
/* to be the same algorithm. */

/* Algorithm 3.3, step 2 "Find generator" */
static void _find_ts_generator(UV *py, UV *pm, UV a, UV k, UV p) {
  UV e, r, y, m, x, ke1;
  /* Assume:  k > 2,  1 < a < p,  p > 2,  k prime,  p prime */
  /* e = valuation_remainder(p-1,k,&r); */
  for (e = 0, r = p-1; !(r % k); r /= k)  e++;
  ke1 = ipow(k, e-1);
  for (x = 2, m = 1; m == 1; x++) {
    y = powmod(x, r, p);
    if (y != 1)
      m = powmod(y, ke1, p);
    MPUassert(x < p, "bad Tonelli-Shanks input\n");
  }
  *py = y;
  *pm = m;
}

static UV _ts_rootmod(UV a, UV k, UV p, UV y, UV m) {
  UV e, r, A, x,  l, T, z, kz;

  /* Assume:  k > 2,  1 < a < p,  p > 2,  k prime,  p prime */
  /* It is not expected to work with prime powers. */

  /* e = valuation_remainder(p-1,k,&r); */
  for (e = 0, r = p-1; !(r % k); r /= k) e++;
  /* p-1 = r * k^e */
  x = powmod(a, modinverse(k % r, r), p);
  A = (a == 0)  ?  0  :  mulmod(powmod(x, k, p), modinverse(a, p), p);

  if (y == 0 && A != 1)
    _find_ts_generator(&y, &m, a, k, p);

  while (A != 1) {
    for (l = 1, T = A;  T != 1;  l++) {
      if (l >= e) return 0;
      z = T;
      T = powmod(T, k, p);
    }
    kz = negmod( znlog_solve(z, m, p, k), k);  /* k = znorder(m,p) */
    m = powmod(m, kz, p);
    T = powmod(y, kz * ipow(k, e-l), p);
    /* In the loop we always end with l < e, so e always gets smaller */
    e = l-1;
    x = mulmod(x, T, p);
    y = powmod(T, k, p);
    if (y <= 1) return 0;  /* In theory this will never be hit. */
    A = mulmod(A, y, p);
  }
  return x;
}

static UV _compute_generator(UV l, UV e, UV r, UV p) {
  UV x, y, m = 1;
  UV lem1 = ipow(l, e-1);
  for (x = 2; m == 1; x++) {
    y = powmod(x, r, p);
    if (y == 1) continue;
    m = powmod(y, lem1, p);
  }
  return y;   /* We might want to also return m */
}

/* Following Pari, we calculate a root of unity to allow finding other roots */
static UV _rootmod_prime_splitk(UV a, UV k, UV p, UV *zeta) {
  UV g;

  if (zeta) *zeta = 1;
  if (a >= p) a %= p;
  if (a == 0 || (a == 1 && !zeta)) return a;

  /* Assume:  k >= 2,  1 < a < p,  p > 2, p prime */

  if (k == 2) {
    if (zeta) *zeta = p-1;
    return _sqrtmod_prime(a, p);
  }

  /* See Algorithm 2.1 of van de Woestijne (2006), or Lindhurst (1997) */
  /* The latter's proposition 7 generalizes to composite p */

  g = gcd_ui(k, p-1);

  if (g != 1) {
    UV fac[MPU_MAX_FACTORS+1];
    UV exp[MPU_MAX_FACTORS+1];
    int i, nfactors = factor_exp(g, fac, exp);
    for (i = 0; a != 0 && i < nfactors; i++) {
      UV y, m, F = fac[i], E = exp[i];
      if (zeta) {
        UV REM, V, Y;
        V = valuation_remainder(p-1, F, &REM);
        Y = _compute_generator(F, V, REM, p);
        *zeta = mulmod(*zeta, powmod(Y, ipow(F, V-E), p), p);
      }
      _find_ts_generator(&y, &m,  a, F, p);
      while (E-- > 0)
        a = _ts_rootmod(a, F, p,  y, m);
    }
  }
  if (g != k) {
    UV kg = k/g,  pg = (p-1)/g;
    a = powmod(a, modinverse(kg % pg, pg), p);
  }
  return a;
}
#endif


#if 0 /* For testing purposes only. */
static UV _trial_rootmod(UV a, UV k, UV n) {
  UV r;
  if (n == 0) return 0;
  if (a >= n) a %= n;
  if (a <= 1) return a;
  for (r = 2; r < n; r++)
    if (powmod(r, k, n) == a)
      return r;
  return 0;
}
static UV* _trial_allsqrtmod(UV* nroots, UV a, UV n) {
  UV i, *roots, numr = 0, allocr = 16;

  if (n == 0) return 0;
  if (a >= n) a %= n;

  New(0, roots, allocr, UV);
  for (i = 0; i <= n/2; i++) {
    if (mulmod(i,i,n) == a) {
      if (numr >= allocr-1)  Renew(roots, allocr += 256, UV);
      roots[numr++] = i;
      if (i != 0 && 2*i != n)
        roots[numr++] = n-i;
    }
  }
  qsort(roots, numr, sizeof(UV), _numcmp);
  *nroots = numr;
  return roots;
}
static UV* _trial_allrootmod(UV* nroots, UV a, UV g, UV n) {
  UV i, *roots, numr = 0, allocr = 16;

  if (n == 0) return 0;
  if (a >= n) a %= n;

  New(0, roots, allocr, UV);
  for (i = 0; i < n; i++) {
    if (powmod(i,g,n) == a) {
      if (numr >= allocr-1)  Renew(roots, allocr += 256, UV);
      roots[numr++] = i;
    }
  }
  *nroots = numr;
  return roots;
}
#endif


/******************************************************************************/
/*                          K-TH ROOT OF N MOD M (splitk)                     */
/******************************************************************************/

#if USE_ROOTMOD_SPLITK
/* Given a solution to r^k = a mod p^(e-1), return r^k = a mod p^e */
static int _hensel_lift(UV *re, UV r, UV a, UV k, UV pe) {
  UV f, fp, d;

  /* UV pe = ipow(p, e); */
  if (a >= pe) a %= pe;
  f = submod(powmod(r, k, pe), a, pe);
  if (f == 0) { *re = r; return 1; }
  fp = mulmod(k, powmod(r, k-1, pe), pe);
  d = divmod(f, fp, pe);
  if (d == 0) return 0;  /* We need a different base root */
  *re = submod(r, d, pe);
  return 1;
}

static UV _rootmod_composite1(UV a, UV k, UV n) {
  UV f, g, e, r;
  UV fac[MPU_MAX_FACTORS+1];
  UV exp[MPU_MAX_FACTORS+1];
  int i, nfactors;

  /* Assume:  k >= 2,  1 < a < n,  n > 2, n composite */

#if 0
  /* For square roots of p^k with gcd(a,p)==1, this is straightforward. */
  if (k == 2 && (i = primepower(n, &f)) && (a % f) > 1) {
    UV x = _sqrtmod_prime(a % f, f);
    UV r = n/f;
    UV j = powmod(x, r, n);
    UV k = powmod(a, (n - r - r + 1) >> 1, n);
    return mulmod(j, k, n);
  }
#endif

  nfactors = factor_exp(n, fac, exp);
  for (i = 0; i < nfactors; i++) {
    f = fac[i];
    /* Find a root mod f.  If none exists, there is no root for n. */
    r = _rootmod_prime_splitk(a%f, k, f, 0);
    if (powmod(r, k, f) != (a%f)) return 0;
    /* If we have a prime power, use Hensel lifting to solve for p^e */
    if (exp[i] > 1) {
      UV fe = f;
      for (e = 2; e <= exp[i]; e++) {
        fe *= f;
        /* We aren't guaranteed a solution, though we usually get one. */
        if (!_hensel_lift(&r, r, a, k, fe)) {
          /* Search for a different base root */
          UV t, m = fe / (f*f);
          for (t = 1; t < f; t++) {
            if (_hensel_lift(&r, r + t*m, a, k, fe))
              break;
          }
          /* That didn't work, do a stronger but time consuming search. */
          if (t >= f) {
            UV afe = a % fe;
            for (r = (a % f); r < fe; r += f)
              if (powmod(r, k, fe) == afe)
                break;
            if (r >= fe) return 0;
          }
        }
      }
      fac[i] = fe;
    }
    exp[i] = r;
  }
  if (chinese(&g, exp, fac, nfactors) != 1) return 0;
  return g;
}
#endif

/******************************************************************************/
/*                          K-TH ROOT OF N MOD M (splitn)                     */
/******************************************************************************/

/* _rootmod_composite2 factors k and combines:
 *    _rootmod_kprime takes prime k along with factored n:
 *       _rootmod_prime_power splits p^e into primes (prime k):
 *          _rootmod_prime finds a root (prime p and prime k)
 *             _sqrtmod_prime  (if k==2)
 *             _ts_prime
 */

#if USE_ROOTMOD_SPLITN && !USE_ROOTMOD_SPLITK
static UV _rootmod_prime(UV a, UV k, UV p) {
  UV r, g;

  /* Assume: p is prime, k is prime */

  if (a >= p) a %= p;
  if (p == 2 || a == 0)  return a;
  if (k == 2) {
    r = _sqrtmod_prime(a,p);
    return (sqrmod(r,p) == a)  ?  r  :  UV_MAX;
  }

  /* If co-prime, we have one root */
  g = gcd_ui(k, p-1);
  if (g == 1)
    return powmod(a, modinverse(k % (p-1), p-1), p);

  /* Check generalized Euler's criterion */
  if (powmod(a, (p-1)/g, p) != 1)
    return UV_MAX;

  return _ts_prime(a, k, p, 0);
}

static UV _rootmod_prime_power(UV a, UV k, UV p, UV e) {
  UV r, s, t, n, np, pk, apk, ered;

  /* Assume: p is prime, k is prime, e >= 1 */

  if (k == 2) return _sqrtmod_prime_power(a, p, e);
  if (e == 1) return _rootmod_prime(a, k, p);

  n  = ipow(p,e);
  pk = ipow(p,k);
  /* Note: a is not modded */

  if ((a % n) == 0)
    return 0;

  if ((a % pk) == 0) {
    apk = a / pk;
    s = _rootmod_prime_power(apk, k, p, e-k);
    if (s == UV_MAX) return UV_MAX;
    return s * p;
  }

  if ((a % p) == 0)
    return UV_MAX;

  ered = (p > 2 || e < 5)  ?  (e+1)>>1  :  (e+3)>>1;
  s = _rootmod_prime_power(a, k, p, ered);
  if (s == UV_MAX) return UV_MAX;

  np = (p != k || (n > (UV_MAX/p)))  ?  n  :  n * p;
  t = powmod(s, k-1, np);
  r = addmod(s, gcddivmod(submod(a,mulmod(t,s,np),np), mulmod(k,t,np), n), n);
  if (powmod(r, k, n) != (a % n)) return UV_MAX;
  return r;
}

static UV _rootmod_kprime(UV a, UV k, UV n, int nfactors, UV *fac, UV *exp) {
  UV i, N, fe, r, s, t, inv;

  /* Assume: k is prime */

  N = ipow(fac[0], exp[0]);
  r = _rootmod_prime_power(a, k, fac[0], exp[0]);
  if (r == UV_MAX) return UV_MAX;
  for (i = 1; i < (UV) nfactors; i++) {
    fe = ipow(fac[i], exp[i]);
    s = _rootmod_prime_power(a, k, fac[i], exp[i]);
    if (s == UV_MAX) return UV_MAX;
    inv = modinverse(N, fe);
    t = mulmod(inv, submod(s % fe,r % fe,fe), fe);
    r = addmod(r, mulmod(N,t,n), n);
    N *= fe;
  }
  return r;

}

static UV _rootmod_composite2(UV a, UV k, UV n) {
  UV i, r;
  UV kfac[MPU_MAX_FACTORS+1];
  UV nfac[MPU_MAX_FACTORS+1];
  UV nexp[MPU_MAX_FACTORS+1];
  int nfactors, kfactors;

  if (n == 0) return 0;
  if (a >= n) a %= n;

  if (n <= 2 || a <= 1) return a;
  if (k <= 1)           return (k == 0) ? 1 : a;

  /* Factor n */
  nfactors = factor_exp(n, nfac, nexp);

  if (is_prime(k))
    return _rootmod_kprime(a, k, n, nfactors, nfac, nexp);

  kfactors = factor(k, kfac);
  r = a;
  for (i = 0; i < (UV)kfactors; i++) {   /* for each prime k */
    r = _rootmod_kprime(r, kfac[i], n, nfactors, nfac, nexp);
    if (r == UV_MAX) { /* Bad path.  We have to use a fallback method. */
#if USE_ROOTMOD_SPLITK
      r = _rootmod_composite1(a,k,n);
#else
      UV *roots, numr;
      roots = allrootmod(&numr,a,k,n);
      r = (numr > 0)  ?  roots[0]  :  UV_MAX;
      Safefree(roots);
#endif
      break;
    }
  }
  return r;
}
#endif

int rootmodp(UV *s, UV a, UV k, UV p) {
  UV r;
  if (p == 0) return 0;
  if (a >= p) a %= p;

  /* return _rootmod_return(_trial_rootmod(a,k,n), s, a, k, p); */

  if      (p <= 2 || a <= 1)  r = a;
  else if (k <= 1)            r = (k == 0) ? 1 : a;
  else if (is_power(a,k))     r = rootint(a,k);
#if USE_ROOTMOD_SPLITK
  else                        r = _rootmod_prime_splitk(a,k,p,0);
#else
  else                        r = _rootmod_composite2(a,k,p);
#endif
  return _rootmod_return(r, s, a, k, p);
}

int rootmod(UV *s, UV a, UV k, UV n) {
  UV r;
  if (n == 0) return 0;
  if (a >= n) a %= n;

  /* return _rootmod_return(_trial_rootmod(a,k,n), s, a, k, n); */

  if      (n <= 2 || a <= 1)  r = a;
  else if (k <= 1)            r = (k == 0) ? 1 : a;
  else if (is_power(a,k))     r = rootint(a,k);
#if USE_ROOTMOD_SPLITK
  else if (is_prime(n))       r = _rootmod_prime_splitk(a,k,n,0);
  else                        r = _rootmod_composite1(a,k,n);
#else
  else                        r = _rootmod_composite2(a,k,n);
#endif
  return _rootmod_return(r, s, a, k, n);
}




/******************************************************************************/
/*               SQRTMOD AND ROOTMOD RETURNING ALL RESULTS                    */
/******************************************************************************/


/* We could alternately just let the allocation fail */
#define MAX_ROOTS_RETURNED 600000000

/* Combine roots using Cartesian product CRT */
static UV* _rootmod_cprod(UV* nroots,
                          UV nr1, UV *roots1, UV p1,
                          UV nr2, UV *roots2, UV p2) {
  UV i, j,  nr, *roots, inv;

  nr = nr1 * nr2;
  if (nr > MAX_ROOTS_RETURNED) croak("Maximum returned roots exceeded");
  New(0, roots, nr, UV);

  inv = modinverse(p1, p2);
  for (i = 0; i < nr1; i++) {
    UV r1 = roots1[i];
    for (j = 0; j < nr2; j++) {
      UV r2 = roots2[j];
#if 0
      UV ca[2], cn[2];
      ca[0] = r1;  cn[0] = p1;
      ca[1] = r2;  cn[1] = p2;
      if (chinese(roots + i * nr2 + j, ca, cn, 2) != 1)
        croak("chinese fail in allrootmod");
#else
      UV t = mulmod(inv, submod(r2 % p2,r1 % p2,p2), p2);
      roots[i * nr2 + j] = addmod(r1, mulmod(p1,t,p1*p2), p1*p2);
#endif
    }
  }
  Safefree(roots1);
  Safefree(roots2);
  *nroots = nr;
  return roots;
}

static UV* _one_root(UV* nroots, UV r) {
  UV *roots;
  New(0, roots, 1, UV);
  roots[0] = r;
  *nroots = 1;
  return roots;
}
static UV* _two_roots(UV* nroots, UV r, UV s) {
  UV *roots;
  New(0, roots, 2, UV);
  roots[0] = r;  roots[1] = s;
  *nroots = 2;
  return roots;
}


/* allsqrtmod algorithm from Hugo van der Sanden, 2021 */

static UV* _allsqrtmodpk(UV *nroots, UV a, UV p, UV k) {
  UV *roots, *roots2, nr2 = 0;
  UV i, j, pk, pj, q, q2, a2;

  pk = ipow(p,k);
  *nroots = 0;

  if ((a % p) == 0) {
    if ((a % pk) == 0) {
      UV low = ipow(p, k >> 1);
      UV high = (k & 1) ? low * p : low;
      if (low > MAX_ROOTS_RETURNED) croak("Maximum returned roots exceeded");
      New(0, roots, low, UV);
      for (i = 0; i < low; i++)
        roots[i] = high * i;
      *nroots = low;
      return roots;
    }
    a2 = a / p;
    if ((a2 % p) != 0)
      return 0;
    pj = pk / p;
    roots2 = _allsqrtmodpk(&nr2, a2/p, p, k-2);
    if (roots2 == 0) return 0;
    *nroots = nr2 * p;
    if (*nroots > MAX_ROOTS_RETURNED) croak("Maximum returned roots exceeded");
    New(0, roots, *nroots, UV);
    for (i = 0; i < nr2; i++)
      for (j = 0; j < p; j++)
        roots[i*p+j] = roots2[i] * p  +  j * pj;
    Safefree(roots2);
    return roots;
  }

  q = _sqrtmod_prime_power(a,p,k);
  if (q == UV_MAX) return 0;

  New(0, roots, 4, UV);
  roots[0] = q;  roots[1] = pk - q;
  if      (p != 2) { *nroots = 2; }
  else if (k == 1) { *nroots = 1; }
  else if (k == 2) { *nroots = 2; }
  else {
    pj = pk / p;
    q2 = mulmod(q, pj-1, pk);
    roots[2] = q2;  roots[3] = pk - q2;
    *nroots = 4;
  }
  return roots;
}

static UV* _allsqrtmodfact(UV *nroots, UV a, UV n, int nf, UV *fac, UV *exp) {
  UV *roots, *roots1, *roots2, nr, nr1, nr2, p, k, pk, n2;

  MPUassert(nf > 0, "empty factor list in _allsqrtmodfact");

  p = fac[0], k = exp[0];
  *nroots = 0;

  /* nr1,roots1 are roots of p^k -- the first prime power */
  roots1 = _allsqrtmodpk(&nr1, a, p, k);
  if (roots1 == 0) return 0;
  if (nf == 1) {
    *nroots = nr1;
    return roots1;
  }
  pk = ipow(p, k);
  n2 = n / pk;

  /* nr2,roots2 are roots of all the rest, found recursively */
  roots2 = _allsqrtmodfact(&nr2, a, n2, nf-1, fac+1, exp+1);
  if (roots2 == 0) return 0;

  roots = _rootmod_cprod(&nr,  nr1, roots1, pk,  nr2, roots2, n2);

  *nroots = nr;
  return roots;
}

UV* allsqrtmod(UV* nroots, UV a, UV n) {
  UV *roots, numr = 0;
  UV fac[MPU_MAX_FACTORS+1];
  UV exp[MPU_MAX_FACTORS+1];
  int nfactors;

  if (n == 0) return 0;
  if (a >= n) a %= n;

  /* return _trial_allsqrtmod(nroots, a, n); */

  if (n <= 2)  return _one_root(nroots, a);

  nfactors = factor_exp(n, fac, exp);
  roots = _allsqrtmodfact(&numr, a, n, nfactors, fac, exp);
  if (numr > 0) qsort(roots, numr, sizeof(UV), _numcmp);
  *nroots = numr;
  return roots;
}


/* allrootmod factors k and combines:
 *    _allrootmod_kprime takes prime k and factored n:
 *       _allrootmod_prime_power splits p^e into primes:
 *          _allrootmod_prime finds all the roots for prime p and prime k
 *             _ts_prime (could alternately call _rootmod_prime_splitk)
 */

static UV* _allrootmod_prime(UV* nroots, UV a, UV k, UV p) {
  UV r, g, z, r2, *roots, numr = 0;

  *nroots = 0;
  if (a >= p) a %= p;

  /* Assume: p is prime, k is prime */

  /* simple case */
  if (p == 2 || a == 0)  return _one_root(nroots, a);

  /* If co-prime, we have one root */
  g = gcd_ui(k, p-1);
  if (g == 1) {
    UV r = powmod(a, modinverse(k % (p-1), p-1), p);
    return _one_root(nroots, r);
  }
  /* At this point k < p.  (k is a prime so if k>=p, g=1) */

  /*  Check generalized Euler's criterion:
   *  r^k = a mod p has a solution iff a^((p-1)/gcd(p-1,k)) = 1 mod p  */
  if (powmod(a, (p-1)/g, p) != 1)
    return 0;

  /* Special case p=3 for performance */
  if (p == 3)  return (k == 2 && a == 1)  ?  _two_roots(nroots, 1, 2)  :  0;

  /* functionally identical:   r = _rootmod_prime_splitk(a, k, p, &z); */
  r = _ts_prime(a, k, p, &z);
  if (powmod(r,k,p) != a || z == 0) croak("allrootmod: failed to find root");

  New(0, roots, k, UV);
  roots[numr++] = r;
  for (r2 = mulmod(r, z, p); r2 != r && numr < k; r2 = mulmod(r2, z, p) )
    roots[numr++] = r2;
  if (r2 != r) croak("allrootmod: excess roots found");

  *nroots = numr;
  return roots;
}


static UV* _allrootmod_prime_power(UV* nroots, UV a, UV k, UV p, UV e) {
  UV n, i, j, pk, s, t, r, numr = 0, *roots = 0, nr2 = 0, *roots2 = 0;

#if 0
  MPUassert(p >= 2, "_allrootmod_prime_power must be given a prime modulus");
  MPUassert(e >= 1, "_allrootmod_prime_power must be given exponent >= 1");
  MPUassert(k >= 2, "_allrootmod_prime_power must be given k >= 2");
  MPUassert(is_prime(k), "_allrootmod_prime_power must be given prime k");
  MPUassert(is_prime(p), "_allrootmod_prime_power must be given prime p");
#endif

  if (e == 1) return _allrootmod_prime(nroots, a, k, p);

  n = ipow(p,e);
  pk = ipow(p, k);
  /* Note: a is not modded */

  if ((a % n) == 0) {

    t = ((e-1) / k) + 1;
    s = ipow(p,t);
    numr  = ipow(p,e-t);
    New(0, roots, numr, UV);
    for (i = 0; i < numr; i++)
      roots[i] = mulmod(i, s, n);

  } else if ((a % pk) == 0) {

    UV apk = a / pk;
    UV pe1 = ipow(p, k-1);
    UV pek = ipow(p, e-k+1);
    roots2 = _allrootmod_prime_power(&nr2, apk, k, p, e-k);
    numr = pe1 * nr2;
    New(0, roots, numr, UV);
    for (i = 0; i < nr2; i++)
      for (j = 0; j < pe1; j++)
        roots[i*pe1+j] = addmod(mulmod(roots2[i],p,n), mulmod(j,pek,n), n);
    Safefree(roots2);

  } else if ((a % p) != 0) {

    UV np = (n > (UV_MAX/p)) ? n : n*p;
    UV ered = (p > 2 || e < 5)  ?  (e+1)>>1  :  (e+3)>>1;
    roots2 = _allrootmod_prime_power(&nr2, a, k, p, ered);

    if (k != p) {

      for (j = 0; j < nr2; j++) {
        s = roots2[j];
        t = powmod(s, k-1, n);
        r = addmod(s,gcddivmod(submod(a,mulmod(t,s,n),n),mulmod(k,t,n),n),n);
        roots2[j] = r;
      }
      roots = roots2;
      numr = nr2;

    } else {

      /* Step 1, transform roots, eliding any that aren't valid */
      for (j = 0; j < nr2; j++) {
        s = roots2[j];
        t = powmod(s, k-1, np);
        r = addmod(s,gcddivmod(submod(a,mulmod(t,s,np),np),mulmod(k,t,np),n),n);
        if (powmod(r, k, n) == (a % n))
          roots2[numr++] = r;
      }
      nr2 = numr;

      /* Step 2, Expand out by k */
      if (nr2 > 0) {
        numr = nr2 * k;
        New(0, roots, numr, UV);
        for (j = 0; j < nr2; j++) {
          r = roots2[j];
          for (i = 0; i < k; i++)
            roots[j*k+i] = mulmod(r, addmod( mulmod(i,n/p,n), 1, n), n);
        }
      }
      Safefree(roots2);

      /* Step 3, Remove any duplicates */
      if (numr == 2 && roots[0] == roots[1])
        numr = 1;
      if (numr > 2) {
        qsort(roots, numr, sizeof(UV), _numcmp);
        for (j = 0, i = 1; i < numr; i++)
          if (roots[j] != roots[i])
            roots[++j] = roots[i];
        numr = j+1;
      }

    }
  }
  *nroots = numr;
  return roots;
}

static UV* _allrootmod_kprime(UV* nroots, UV a, UV k, UV n, int nfactors, UV *fac, UV *exp) {
  UV i, fe, N, *roots = 0, *roots2, numr = 0, nr2;

  if (k == 2) return _allsqrtmodfact(nroots, a, n, nfactors, fac, exp);

  *nroots = 0;
  N = ipow(fac[0], exp[0]);
  roots = _allrootmod_prime_power(&numr, a, k, fac[0], exp[0]);
  if (numr == 0) { Safefree(roots); return 0; }
  for (i = 1; i < (UV) nfactors; i++) {
    fe = ipow(fac[i], exp[i]);
    roots2 = _allrootmod_prime_power(&nr2, a, k, fac[i], exp[i]);
    if (nr2 == 0) { Safefree(roots); Safefree(roots2); return 0; }
    /* Cartesian product using CRT.  roots and roots2 are freed. */
    roots = _rootmod_cprod(&numr,  numr, roots, N,  nr2, roots2, fe);
    N *= fe;
  }
  MPUassert(N == n, "allrootmod: Incorrect factoring");

  *nroots = numr;
  return roots;
}

UV* allrootmod(UV* nroots, UV a, UV k, UV n) {
  UV i, numr = 0, *roots = 0;
  UV kfac[MPU_MAX_FACTORS+1];
  UV nfac[MPU_MAX_FACTORS+1];
  UV nexp[MPU_MAX_FACTORS+1];
  int nfactors, kfactors;

  /* return _trial_allrootmod(nroots, a, k, n); */

  *nroots = 0;
  if (n == 0) return 0;
  if (a >= n) a %= n;

  if (n <= 2 || k == 1)
    return _one_root(nroots, a);   /* n=1 => [0],  n=2 => [0] or [1] */

  if (k == 0) {
    if (a != 1) return 0;
    New(0, roots, n, UV);
    for (i = 0; i < n; i++)
      roots[i] = i;
    *nroots = n;
    return roots;
  }

  /* Factor n */
  nfactors = factor_exp(n, nfac, nexp);

  if (is_prime(k)) {

    roots = _allrootmod_kprime(&numr, a, k, n, nfactors, nfac, nexp);

  } else { /* Split k into primes */

    kfactors = factor(k, kfac);
    roots = _allrootmod_kprime(&numr, a, kfac[0], n, nfactors, nfac, nexp);

    for (i = 1; numr > 0 && i < (UV)kfactors; i++) {   /* for each prime k */
      UV j, t, allocr = numr, primek = kfac[i];
      UV *roots2 = 0, nr2 = 0,  *roots3 = 0, nr3 = 0;
      New(0, roots3, allocr, UV);
      for (j = 0; j < numr; j++) {         /* get a list from each root */
        roots2 = _allrootmod_kprime(&nr2, roots[j], primek, n, nfactors, nfac, nexp);
        if (nr2 == 0) continue;
        /* Append to roots3 */
        if (nr3 + nr2 > MAX_ROOTS_RETURNED) croak("Maximum returned roots exceeded");
        if (nr3 + nr2 >= allocr)  Renew(roots3, allocr += nr2, UV);
        for (t = 0; t < nr2; t++)
          roots3[nr3++] = roots2[t];
        Safefree(roots2);
      }
      /* We've walked through all the roots combining to roots3 */
      Safefree(roots);
      roots = roots3;
      numr = nr3;
    }

  }
  if (numr > 1)
    qsort(roots, numr, sizeof(UV), _numcmp);
  *nroots = numr;
  return roots;
}
