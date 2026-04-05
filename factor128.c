#include <math.h>
#include "ptypes.h"
#include "constants.h"
#include "factor.h"
#include "factor128.h"
#include "cache.h"
#include "sieve.h"
#include "util.h"       /* for verbose and next_prime */

#if BITS_PER_WORD == 64 && HAVE_UINT128

/* Decimal string for a uint128_t, for debugging.
 * printf("  found factor %s of %s\n", u128_str(f), u128_str(t)); */
static const char *u128_str(uint128_t n) {
  static char bufs[4][42];
  static int  idx = 0;
  char *buf = bufs[idx++ & 3];
  char *p   = buf + 41;
  *p = '\0';
  if (n == 0) { *--p = '0'; return p; }
  while (n) { *--p = '0' + (int)(n % 10); n /= 10; }
  return p;
}

/*****************************************************************************
 * 128-bit modular arithmetic
 * All functions assume 0 <= a, b < n and n > 1.
 *****************************************************************************/

static INLINE uint128_t addmod128(uint128_t a, uint128_t b, uint128_t n) {
  uint128_t r = a + b;
  /* Handle both unsigned overflow (r < a) and r >= n */
  if (r < a || r >= n) r -= n;
  return r;
}

static INLINE uint128_t submod128(uint128_t a, uint128_t b, uint128_t n) {
  return (a >= b) ? a - b : n - (b - a);
}

/* (a + n) / 2 mod n, for odd n */
static INLINE uint128_t half_mod128(uint128_t a, uint128_t n) {
  if (a & 1) {
    uint128_t t = a + n;
    /* If overflow occurred, the true value is >= 2^128, so after /2
     * the top bit is set.  Reconstruct with (1<<127) | (t>>1). */
    return (t < a) ? ((uint128_t)1 << 127) | (t >> 1) : t >> 1;
  }
  return a >> 1;
}

/* Binary double-and-add mulmod.  O(128) addmod128 calls.
 * Used only for primality testing (variable modulus). */
static uint128_t mulmod128(uint128_t a, uint128_t b, uint128_t n) {
  uint128_t r = 0;
  a %= n;
  while (b > 0) {
    if (b & 1) r = addmod128(r, a, n);
    a = addmod128(a, a, n);
    b >>= 1;
  }
  return r;
}

#define sqrmod128(a, n)  mulmod128(a, a, n)

static uint128_t powmod128(uint128_t a, uint128_t k, uint128_t n) {
  uint128_t r = 1;
  a %= n;
  while (k > 0) {
    if (k & 1) r = mulmod128(r, a, n);
    k >>= 1;
    if (k) a = sqrmod128(a, n);
  }
  return r;
}

static uint128_t gcd128(uint128_t a, uint128_t b) {
  while (b) { uint128_t t = b; b = a % b; a = t; }
  return a;
}


/*****************************************************************************
 * Montgomery multiplication for fixed modulus.
 *
 * R = 2^128.  Given odd n, we keep values in Montgomery form: â = a*R mod n.
 * mont_mulmod(â, b̂) = â*b̂*R^{-1} mod n  (REDC algorithm).
 * Correct for all odd n < 2^128, including n > 2^127.
 *
 * Setup: ninv = -n^{-1} mod 2^128 (via Hensel lifting).
 *        r2   = R^2 mod n          (for converting to Montgomery form).
 *****************************************************************************/

typedef struct { uint128_t n, ninv, r2; } mont128_t;

/* Low 128 bits of the full 256-bit product a*b. */
static INLINE uint128_t lo128(uint128_t a, uint128_t b) {
  uint64_t a0 = (uint64_t)a,  a1 = (uint64_t)(a >> 64);
  uint64_t b0 = (uint64_t)b,  b1 = (uint64_t)(b >> 64);
  uint128_t lo = (uint128_t)a0 * b0;
  lo += (uint128_t)a0 * b1 << 64;
  lo += (uint128_t)a1 * b0 << 64;
  return lo;
}

/* Full 256-bit product a*b, split into lo/hi 128-bit halves. */
static INLINE void mul256(uint128_t a, uint128_t b,
                          uint128_t *lo_out, uint128_t *hi_out) {
  uint64_t a0 = (uint64_t)a,  a1 = (uint64_t)(a >> 64);
  uint64_t b0 = (uint64_t)b,  b1 = (uint64_t)(b >> 64);

  uint128_t p00 = (uint128_t)a0 * b0;
  uint128_t p01 = (uint128_t)a0 * b1;
  uint128_t p10 = (uint128_t)a1 * b0;
  uint128_t p11 = (uint128_t)a1 * b1;

  /* lo = p00 + ((p01 + p10) << 64), tracked with carry into hi */
  uint128_t mid  = p01 + p10;
  int       mc   = (mid < p01);          /* carry from mid addition */
  uint128_t lo   = p00 + (mid << 64);
  int       lc   = (lo < p00);           /* carry from lo addition */

  *lo_out = lo;
  *hi_out = p11 + (mid >> 64) + ((uint128_t)mc << 64) + lc;
}

/* High 128 bits of the full 256-bit product a*b (lo half discarded). */
static INLINE uint128_t mulhi128(uint128_t a, uint128_t b) {
  uint64_t a0 = (uint64_t)a,  a1 = (uint64_t)(a >> 64);
  uint64_t b0 = (uint64_t)b,  b1 = (uint64_t)(b >> 64);

  uint128_t p00 = (uint128_t)a0 * b0;
  uint128_t p01 = (uint128_t)a0 * b1;
  uint128_t p10 = (uint128_t)a1 * b0;
  uint128_t p11 = (uint128_t)a1 * b1;

  uint128_t mid = p01 + p10;
  int       mc  = (mid < p01);
  uint128_t lo  = p00 + (mid << 64);     /* only needed for carry lc */
  int       lc  = (lo < p00);

  return p11 + (mid >> 64) + ((uint128_t)mc << 64) + lc;
}

/* Montgomery REDC: returns a*b*R^{-1} mod n.
 * Inputs must already be in [0, n). */
static INLINE uint128_t mont_mulmod128(uint128_t a, uint128_t b,
                                       const mont128_t *ctx) {
  uint128_t lo, hi;
  mul256(a, b, &lo, &hi);

  /* m = lo * ninv  (only low 128 bits needed) */
  uint128_t m = lo128(lo, ctx->ninv);

  /* m·n ≡ -lo (mod R), so lo + mn_lo = R when lo≠0, 0 otherwise.
   * carry into the high half is simply (lo != 0); no need to compute mn_lo. */
  int carry = (lo != 0);
  uint128_t mn_hi = mulhi128(m, ctx->n);

  /* hi + mn_hi + carry may exceed 2^128 when n > 2^127 (true result still
   * in [0, 2n), so one conditional subtract suffices).  Track overflow bit. */
  uint128_t t = hi + mn_hi;
  int ov = (t < hi);
  t += (uint128_t)carry;
  ov += (t < (uint128_t)carry);

  /* When ov=1 the true value is t + 2^128 >= n; uint128_t wraparound in
   * t - ctx->n correctly gives (t + 2^128 - n) mod 2^128. */
  if (ov || t >= ctx->n) return t - ctx->n;
  return t;
}

#define mont_sqrmod128(a, ctx)  mont_mulmod128(a, a, ctx)

/* Convert a (normal) value into Montgomery form: a*R mod n. */
static INLINE uint128_t mont_enter128(uint128_t a, const mont128_t *ctx) {
  return mont_mulmod128(a % ctx->n, ctx->r2, ctx);
}

/* Convert from Montgomery form back to normal: â*R^{-1} mod n. */
static INLINE uint128_t mont_exit128(uint128_t a, const mont128_t *ctx) {
  return mont_mulmod128(a, 1, ctx);
}

/* Initialise a Montgomery context for modulus n (must be odd). */
static void mont_setup128(mont128_t *ctx, uint128_t n) {
  ctx->n = n;

  /* ninv = n^{-1} mod 2^128 via Hensel lifting.
   * Each step doubles the number of correct bits: x = x*(2 - n*x). */
  uint128_t x = (3*n)^2;                 /*   4 bits (n odd) */
  x *= 2 - n * x;                        /*   8 bits */
  x *= 2 - n * x;                        /*  16 bits */
  x *= 2 - n * x;                        /*  32 bits */
  x *= 2 - n * x;                        /*  64 bits */
  x *= 2 - n * x;                        /* 128 bits */
  /* REDC add-form needs m = T*(-n^{-1}) mod R; store the negation */
  ctx->ninv = -x;

  /* r2 = R^2 mod n = (2^128)^2 mod n = 2^256 mod n. */
  ctx->r2 = powmod128(2,256,n);
}

/*****************************************************************************
 * Montgomery powmod — exponent k is a plain UV, a/return in Montgomery form.
 *****************************************************************************/

static uint128_t mont_powmod128(uint128_t a, UV k, const mont128_t *ctx) {
  uint128_t r = mont_enter128(1, ctx);
  while (k > 0) {
    if (k & 1) r = mont_mulmod128(r, a, ctx);
    k >>= 1;
    if (k) a = mont_sqrmod128(a, ctx);
  }
  return r;
}

/*****************************************************************************
 * Primality — BPSW (Miller-Rabin base 2 + strong Lucas, Selfridge params)
 *****************************************************************************/

/* Returns 1 if n is a strong pseudoprime base `base`, 0 otherwise.
 * Assumes n > 2, n odd. */
static int miller_rabin128(uint128_t n, uint128_t base) {
  uint128_t d = n - 1;
  int s = 0;

  while (!(d & 1)) { d >>= 1; s++; }

  uint128_t x = powmod128(base % n, d, n);
  if (x == 1 || x == n - 1) return 1;
  while (--s > 0) {
    x = sqrmod128(x, n);
    if (x == n - 1) return 1;
  }
  return 0;
}

/* Jacobi symbol (a/n), n odd > 0.  a may be any signed 128-bit value. */
static int jacobi128(int128_t a_in, uint128_t n) {
  uint128_t a;
  int j = 1;

  if (n == 1) return 1;

  if (a_in < 0) {
    a = (uint128_t)(-a_in) % n;
    /* (-1/n) = (-1)^((n-1)/2); n odd so n%4 in {1,3} */
    if ((n & 3) == 3) j = -j;
  } else {
    a = (uint128_t)a_in % n;
  }

  while (a != 0) {
    while (!(a & 1)) {
      a >>= 1;
      /* (2/n) = (-1)^((n^2-1)/8) */
      if ((n & 7) == 3 || (n & 7) == 5) j = -j;
    }
    /* Quadratic reciprocity */
    uint128_t t = a; a = n; n = t;
    if ((a & 3) == 3 && (n & 3) == 3) j = -j;
    a %= n;
  }
  return (n == 1) ? j : 0;
}

/* Lucas sequence mod n: given (D, P=1, Q), index k, computes (U_k, V_k, Q^k).
 * Uses left-to-right binary method. */
static void lucas_seq128(uint128_t *U, uint128_t *V, uint128_t *Qk,
                         int128_t D_s, int64_t Q_s,
                         uint128_t k, uint128_t n) {
  /* Work in [0,n) using addmod/submod for signed values */
  uint128_t D  = (D_s >= 0) ? (uint128_t)D_s % n : n - (uint128_t)(-D_s) % n;
  uint128_t Q  = (Q_s >= 0) ? (uint128_t)Q_s % n : n - (uint128_t)(-Q_s) % n;
  uint128_t P  = 1;

  uint128_t u = 1, v = P, q = Q;
  uint128_t u2, v2, q2;
  uint128_t bit;

  if (k == 0) { *U = 0; *V = 2 % n; *Qk = 1; return; }

  /* Find highest bit below the leading 1 */
  bit = (uint128_t)1 << 126;
  while (bit > k) bit >>= 1;
  bit >>= 1;

  while (bit > 0) {
    /* Double: U_{2m}, V_{2m}, Q^{2m} */
    u2 = mulmod128(u, v, n);
    v2 = submod128(mulmod128(v, v, n), addmod128(q, q, n), n);
    q2 = sqrmod128(q, n);

    if (k & bit) {
      /* Advance by 1: U_{2m+1}, V_{2m+1}, Q^{2m+1} */
      /* U_{2m+1} = (P*U_{2m} + V_{2m}) / 2
         V_{2m+1} = (D*U_{2m} + P*V_{2m}) / 2  (P=1)
         Q^{2m+1} = Q^{2m} * Q                              */
      uint128_t pu2_v2 = addmod128(u2, v2, n); /* P*u2 + v2, P=1 */
      uint128_t du2_v2 = addmod128(mulmod128(D, u2, n), v2, n);
      u = half_mod128(pu2_v2, n);
      v = half_mod128(du2_v2, n);
      q = mulmod128(q2, Q, n);
    } else {
      u = u2; v = v2; q = q2;
    }
    bit >>= 1;
  }

  *U = u; *V = v; *Qk = q;
}

/* Strong Lucas probable prime test with Selfridge parameters.
 * n must be odd, > 2, not a perfect square.
 * D is chosen from sequence 5,-7,9,-11,13,... until jacobi(D,n)=-1.
 * P=1, Q=(1-D)/4. */
static int is_strong_lucas_pp128(uint128_t n) {
  int64_t D;
  int64_t sign = 1;
  int abs_D;

  for (abs_D = 5; ; abs_D += 2, sign = -sign) {
    D = sign * abs_D;
    int j = jacobi128((int128_t)D, n);
    if (j == -1) break;
    if (j == 0)  return 0; /* gcd(D,n) > 1 => composite */
    if (abs_D == 21) {
      /* Check for perfect square before going further */
      /* (simple check: if n is a perfect square, jacobi never hits -1) */
      /* We rely on caller having excluded perfect squares or just continue */
    }
  }

  int64_t Q_s = (1 - D) / 4;
  uint128_t d = n + 1;
  int s = 0;
  while (!(d & 1)) { d >>= 1; s++; }

  uint128_t U, V, Qk;
  lucas_seq128(&U, &V, &Qk, D, Q_s, d, n);

  if (U == 0) return 1;
  for (int r = 0; r < s; r++) {
    if (V == 0) return 1;
    /* V_{2m} = V_m^2 - 2*Q^m  (mod n) */
    V = submod128(sqrmod128(V, n), addmod128(Qk, Qk, n), n);
    Qk = sqrmod128(Qk, n);
  }
  return 0;
}

/* BPSW primality test, no trial division.  Assumes n > 2011 and odd.
 * Called from the factoring loop where small factors are already removed. */
static int is_bpsw128(uint128_t n) {
  if (!miller_rabin128(n, 2)) return 0;
  if (!is_strong_lucas_pp128(n)) return 0;
  return 1;
}

/* Returns 1 if n is (probably) prime, 0 if composite.
 * Uses trial division up to 2011, then BPSW. */
static int is_prime128(uint128_t n) {
  if (n < 2) return 0;
  if (n < 4) return 1;
  if (!(n & 1)) return 0;

  /* Trial division using primes_small[] */
  for (int i = 2; i < NPRIMES_SMALL; i++) {
    UV p = primes_small[i];
    if ((uint128_t)p * p > n) return 1;
    if (n % p == 0) return 0;
  }

  return is_bpsw128(n);
}

/*****************************************************************************
 * SQUFOF (Square Form Factorization) for 128-bit inputs.
 *
 * Adapted from squfof126.c (Dana Jacobsen / Ben Buhrow racing SQUFOF).
 * Internal state stays in uint64_t (values are O(sqrt(N*mult)) < 2^63).
 * Only the two places that need N itself use uint128_t:
 *   • isqrt128()  — initial sqrt
 *   • So = (N - Ro^2)/S  — one division during the reduction phase
 *   • gcd(Ro, N)  — gcd at the end of each inner cycle
 * Limit: N*mult must fit in 126 bits (same bound as the original).
 *****************************************************************************/

/* Integer square root of a 128-bit value: returns floor(sqrt(n)). */
static uint64_t isqrt128(uint128_t n) {
  if (n == 0) return 0;
  /* Initial estimate: r = 2^ceil(bits/2), an overestimate of sqrt(n) */
  uint128_t r = 1, tmp = n;
  while (tmp > 3) { tmp >>= 2; r <<= 1; }
  r <<= 1;   /* ensure we start from above so Newton descends to floor */
  /* Newton's method descends from above: r = (r + n/r) / 2 */
  while (1) {
    uint128_t rn = (r + n / r) >> 1;
    if (rn >= r) break;
    r = rn;
  }
  /* r may still be one too high; adjust down */
  while (r * r > n) r--;
  return (uint64_t)r;
}

typedef struct {
  int      valid;
  uint64_t P, bn, Qn, Q0, b0, it, imax, mult;
} squfof128_mult_t;

/* One racing "slice" for a single multiplier.
 * Returns a 64-bit factor (possibly divisible by mult) or 0 to continue. */
static uint64_t squfof_unit128(uint128_t nn64, squfof128_mult_t *ms) {
  uint64_t P  = ms->P,  bn = ms->bn, Qn = ms->Qn;
  uint64_t Q0 = ms->Q0, b0 = ms->b0;
  uint64_t i  = ms->it, imax = i + ms->imax;
  uint64_t t1, t2;

#define SQUFOF128_ITER \
    t1 = P;            \
    P  = bn*Qn - P;    \
    t2 = Qn;           \
    Qn = Q0 + bn*(t1-P); \
    Q0 = t2;           \
    bn = (b0 + P) / Qn; \
    i++;

  /* If iteration count is odd, do one step to make it even */
  if (i & 1) { SQUFOF128_ITER }

  while (1) {
    if (i >= imax) {
      ms->P = P; ms->bn = bn; ms->Qn = Qn; ms->Q0 = Q0; ms->it = i;
      return 0;
    }

    SQUFOF128_ITER  /* even iteration */

    /* Check if Qn is a perfect square (quick residue filter) */
    t2 = Qn & 127;
    if (!((t2*0x8bc40d7d) & (t2*0xa1e2f5d1) & 0x14020a)) {
      t1 = (uint32_t) sqrt((double)Qn);
      if (Qn == t1*t1) break;      /* perfect square found */
    }

    SQUFOF128_ITER  /* odd iteration */
  }

  uint64_t S = t1;
  ms->it = i;

  /* Reduce to G0 using the square root S */
  uint64_t Ro  = P + S * ((b0 - P) / S);
  uint64_t So  = (uint64_t)((nn64 - (uint128_t)Ro * Ro) / S);
  uint64_t bbn = (b0 + Ro) / So;

#define SQUFOF128_SYM \
    t1 = Ro;                      \
    Ro = bbn*So - Ro;             \
    if (Ro == t1) break;          \
    t2 = So;                      \
    So = S + bbn*(t1 - Ro);      \
    S  = t2;                      \
    bbn = (b0 + Ro) / So;

  /* Walk to the symmetry point */
  uint64_t j = 0;
  while (1) {
    SQUFOF128_SYM SQUFOF128_SYM SQUFOF128_SYM SQUFOF128_SYM
    if (j++ > imax) { ms->valid = 0; return 0; }
  }

  /* gcd(Ro, nn64) */
  t1 = (uint64_t) gcd128((uint128_t)Ro, nn64);
  return (t1 > 1) ? t1 : 0;

#undef SQUFOF128_ITER
#undef SQUFOF128_SYM
}

/* Multiplier table from Gower & Wagstaff 2008, plus 1680-series. */
static const uint64_t squfof128_multipliers[] = {
  33*1680, 11*1680, 66*1680,  3*1680,  2*1680,  6*1680, 22*1680, 78*1680,
   1*1680, 26*1680, 39*1680, 13*1680,102*1680, 30*1680, 34*1680, 10*1680,
  15*1680, 51*1680,  5*1680, 57*1680, 17*1680, 19*1680,
  3*5*7*11, 3*5*7,  3*5*7*11*13, 3*5*7*13, 3*5*7*11*17, 3*5*11,
  3*5*7*17, 3*5,    3*5*7*11*19, 3*5*11*13,3*5*7*19,    3*5*7*13*17,
  3*5*13,   3*7*11, 3*7,         5*7*11,   3*7*13,      5*7,
  3*5*17,   5*7*13, 3*5*19,      3*11,     3*7*17,      3,
  3*11*13,  5*11,   3*7*19,      3*13,     5,           5*11*13,
  5*7*19,   5*13,   7*11,        7,        3*17,        7*13,
  11,       1
};
#define NSQUFOF128_MULT \
    (int)(sizeof(squfof128_multipliers)/sizeof(squfof128_multipliers[0]))

/* Returns a nontrivial factor of n (up to 126-bit composites), or 0.
 * rounds: max total inner iterations across all racing multipliers. */
static uint128_t squfof128(uint128_t n, UV rounds) {
  /* 126-bit limit: internal state uses uint64_t which holds values < 2^63 */
  static const uint128_t MAX_NN64 = ((uint128_t)1 << 126) - 1;

  squfof128_mult_t ms[NSQUFOF128_MULT];
  int i, mults_racing = NSQUFOF128_MULT;
  UV rounds_done = 0;

  for (i = 0; i < NSQUFOF128_MULT; i++) {
    ms[i].valid = -1;
    ms[i].it    = 0;
  }

  while (mults_racing > 0 && rounds_done < rounds) {
    for (i = 0; i < NSQUFOF128_MULT && rounds_done < rounds; i++) {
      if (ms[i].valid == 0) continue;

      uint64_t mult = squfof128_multipliers[i];
      uint128_t nn64 = n * mult;

      if (ms[i].valid == -1) {
        if (nn64 > MAX_NN64) {
          ms[i].valid = 0; mults_racing--; continue;
        }
        uint64_t sqrtnn64 = isqrt128(nn64);
        uint128_t rem     = nn64 - (uint128_t)sqrtnn64 * sqrtnn64;
        if (rem == 0) {
          /* nn64 is a perfect square — find factor via gcd */
          uint64_t f = (uint64_t)gcd128((uint128_t)sqrtnn64, n);
          if (f > 1 && f < (uint64_t)-1) return (uint128_t)f;
        }
        /* imax ≈ 0.5 * nn64^(1/5): compute via Newton's method on 5th root */
        uint64_t imax_v;
        {
          /* one Newton step: x = (4*x + nn64/(x^4)) / 5 would work but
           * is complex; simpler: use floating-point with known precision */
          double fbits = 0;
          { uint128_t tmp = nn64; while (tmp >>= 1) fbits += 1.0; }
          imax_v = (uint64_t)(0.5 * pow(2.0, fbits / 5.0));
          if (imax_v < 20) imax_v = 20;
        }
        ms[i].valid = 1;
        ms[i].Q0    = 1;
        ms[i].b0    = sqrtnn64;
        ms[i].P     = sqrtnn64;
        ms[i].Qn    = (uint64_t)rem;
        ms[i].bn    = (2 * sqrtnn64) / ms[i].Qn;
        ms[i].mult  = mult;
        ms[i].imax  = imax_v;
      }

      /* Give this multiplier at most (rounds - rounds_done) iterations */
      if (mults_racing == 1 || ms[i].imax > (rounds - rounds_done))
        ms[i].imax = rounds - rounds_done;

      uint64_t f64 = squfof_unit128(nn64, &ms[i]);
      if (f64 > 1) {
        /* Divide out any factor shared with the multiplier */
        uint64_t g = (uint64_t)gcd128((uint128_t)f64, (uint128_t)mult);
        uint64_t f64red = f64 / g;
        if (f64red > 1) {
          /* Verify f64red actually divides n (guards against spurious results
           * from degenerate SQUFOF states, e.g. when Ro computation overflows) */
          uint128_t factor = gcd128((uint128_t)f64red, n);
          if (factor > 1 && factor < n) return factor;
        }
        /* Trivial or spurious factor — retire this multiplier */
        ms[i].valid = 0;
      }
      if (ms[i].valid == 0) mults_racing--;
      rounds_done += ms[i].imax;
    }
  }
  return 0;
}

/*****************************************************************************
 * Pollard P-1 factoring — Stage 1 + Stage 2 using Montgomery arithmetic.
 *
 * Stage 1: a = 2^E mod n  where E = ∏ q^e, q prime, q^e ≤ B1.
 *          Check gcd(a−1, n).  GCD batched every 32 primes; backtrack on
 *          gcd == n to isolate the exact factor.
 *
 * Stage 2: standard prime continuation.  bm = a after stage 1.
 *          For each prime p ∈ (B1, B2]: a = a · bm^(p−p_prev) mod n.
 *          Accumulate b = ∏(a−1) mod n; GCD every 64 primes.
 *          Small prime-gap powers bm^{2k} (k=1..111) are cached lazily.
 *****************************************************************************/

static uint128_t pminus1_128(uint128_t n, UV B1, UV B2) {
  if (B1 < 7) return 0;

  mont128_t ctx;
  mont_setup128(&ctx, n);

  /*--- Stage 1 -----------------------------------------------------------*/
  uint128_t a     = mont_enter128(2, &ctx);
  uint128_t savea = a;
  UV q = 2, saveq = 2;
  UV sqrtB1 = (UV) sqrt((double)B1);
  UV j = 15;   /* checkpoint counter, start offset like GMP */

  START_DO_FOR_EACH_PRIME(2,B1) {
    UV k = p;
    if (p <= sqrtB1) { UV km = B1/p; while (k <= km) k *= p; }
    a = mont_powmod128(a, k, &ctx);
    if ((j++ % 32) == 0) {
      uint128_t an = mont_exit128(a, &ctx);
      uint128_t g  = gcd128(an > 0 ? an - 1 : n - 1, n);
      if (g == n)  RETURN_FROM_EACH_PRIME(goto stage1_backtrack);
      if (g > 1)   RETURN_FROM_EACH_PRIME(return g);
      saveq = p;
      savea = a;
    }
  } END_DO_FOR_EACH_PRIME

  /* Final stage-1 GCD */
  {
    uint128_t an = mont_exit128(a, &ctx);
    uint128_t g  = gcd128(an > 0 ? an - 1 : n - 1, n);
    if (g == n)  goto stage1_backtrack;
    if (g > 1)   return g;
  }
  goto stage2;

stage1_backtrack:
  /* Multiple factors found -- redo from last checkpoint one prime at a time */
  a = savea;
  for (q = saveq; q <= B1; q = next_prime(q)) {
    UV k = q;
    if (q <= sqrtB1) { UV km = B1/q; while (k <= km) k *= q; }
    a = mont_powmod128(a, k, &ctx);
    uint128_t an = mont_exit128(a, &ctx);
    uint128_t g  = gcd128(an > 0 ? an - 1 : n - 1, n);
    if (g == n)  return 0;
    if (g > 1)   return g;
  }
  return 0;

stage2:
  if (B2 <= B1) return 0;

  /*--- Stage 2 -----------------------------------------------------------*/
  /* bm = a (end-of-stage-1 value).  We precompute bm^2, bm^4, …, bm^44
   * eagerly, and lazily cache up to bm^222.  The index is:
   *   qdiff = (prime_gap)/2 − 1,  bm_pow[qdiff] = bm^(prime_gap).       */
  uint128_t bm      = a;
  uint128_t bm_sq2  = mont_sqrmod128(bm, &ctx);      /* bm^2 */
  uint128_t bm_pow[111];
  int       bm_init[111];
  int       ii;
  for (ii = 0; ii < 111; ii++) bm_init[ii] = 0;

  bm_pow[0]  = bm_sq2;  bm_init[0]  = 1;         /* bm^2  */
  for (ii = 1; ii < 22; ii++) {                  /* bm^4 … bm^44 */
    bm_pow[ii] = mont_mulmod128(bm_pow[ii-1], bm_sq2, &ctx);
    bm_init[ii] = 1;
  }

  q = next_prime(B1);  /* q is now the first prime > B1 */
  /* Advance a to the first stage-2 prime (= q, first prime > B1) */
  a = mont_powmod128(a, q, &ctx);

  uint128_t R_modn = mont_enter128(1, &ctx);   /* Mont. form of 1 */
  uint128_t b = R_modn;                        /* accumulator */
  j = 31;

  START_DO_FOR_EACH_PRIME(q,B2) {
    uint128_t step;
    UV gap   = p - q;  /* q is the previous prime */
    UV pdiff = gap / 2 - 1;
    if (pdiff < 111) {
      if (!bm_init[pdiff]) {
        /* lazily compute bm^gap = bm^(2*(qdiff+1)) */
        bm_pow[pdiff] = mont_powmod128(bm, gap, &ctx);
        bm_init[pdiff] = 1;
      }
      step = bm_pow[pdiff];
    } else {
      /* rare large gap — compute on the fly */
      step = mont_powmod128(bm, gap, &ctx);
    }
    a = mont_mulmod128(a, step, &ctx);
    /* (a − 1) in Montgomery form: a_mont − R_modn */
    uint128_t am1 = submod128(a, R_modn, n);
    b = mont_mulmod128(b, am1, &ctx);

    if ((j++ % 64) == 0) {
      uint128_t g = gcd128(mont_exit128(b, &ctx), n);
      if (g > 1 && g < n)  RETURN_FROM_EACH_PRIME(return g);
      if (g == n)          break;
    }
    q = p;
  } END_DO_FOR_EACH_PRIME
  /* Final stage-2 GCD */
  {
    uint128_t g = gcd128(mont_exit128(b, &ctx), n);
    if (g > 1 && g < n)  return g;
  }
  return 0;
}

/*****************************************************************************
 * Tiny ECM — elliptic curve factoring.
 *
 * Targets factors in the 40–55 bit range that P-1 and SQUFOF miss.
 * Uses Suyama's parameterization, the Montgomery ladder, and the existing
 * mont_mulmod128 infrastructure.
 *****************************************************************************/

/* Returns a^{-1} mod n, or 0 if gcd(a,n) > 1. */
static uint128_t modinv128(uint128_t a, uint128_t n) {
  uint128_t r0 = n,  r1 = a;   /* unsigned so n > 2^127 works correctly */
  int128_t  s0 = 0, s1 = 1;
  while (r1) {
    uint128_t q = r0 / r1, t;
    t = r0 - q * r1;  r0 = r1;  r1 = t;
    int128_t ts = s0 - (int128_t)q * s1;  s0 = s1;  s1 = ts;
  }
  if (r0 != 1) return 0;
  return s0 < 0 ? n - (uint128_t)(-s0) : (uint128_t)s0;
}

/* Projective (X:Z) point on a Montgomery curve.  Values in Montgomery form. */
typedef struct { uint128_t X, Z; } ecpt128_t;

/* Point doubling: R = 2P on By²=x³+Ax²+x, A24 = (A+2)/4. */
static INLINE void ecm_double128(ecpt128_t *R, const ecpt128_t *P,
                                 uint128_t A24, const mont128_t *ctx) {
  uint128_t u = mont_sqrmod128(submod128(P->X, P->Z, ctx->n), ctx); /* (X-Z)² */
  uint128_t v = mont_sqrmod128(addmod128(P->X, P->Z, ctx->n), ctx); /* (X+Z)² */
  uint128_t w = submod128(v, u, ctx->n);                             /* 4XZ   */
  R->X = mont_mulmod128(u, v, ctx);
  R->Z = mont_mulmod128(w, addmod128(u, mont_mulmod128(A24, w, ctx), ctx->n), ctx);
}

/* Differential addition: R = P+Q given P-Q = Pm.  R may alias P or Q. */
static INLINE void ecm_dadd128(ecpt128_t *R,
                               const ecpt128_t *P, const ecpt128_t *Q,
                               const ecpt128_t *Pm, const mont128_t *ctx) {
  uint128_t u = mont_mulmod128(submod128(P->X, P->Z, ctx->n),
                               addmod128(Q->X, Q->Z, ctx->n), ctx);
  uint128_t v = mont_mulmod128(addmod128(P->X, P->Z, ctx->n),
                               submod128(Q->X, Q->Z, ctx->n), ctx);
  uint128_t s = addmod128(u, v, ctx->n);
  uint128_t d = submod128(u, v, ctx->n);
  R->X = mont_mulmod128(mont_sqrmod128(s, ctx), Pm->Z, ctx);
  R->Z = mont_mulmod128(mont_sqrmod128(d, ctx), Pm->X, ctx);
}

/* Montgomery-ladder scalar multiply: R = k*P.
 * R1-R0 = P is maintained, so the original P is always the Pm arg.
 * Safe to call with R aliasing P: R0 is a local copy made before any writes. */
static void ecm_mul128(ecpt128_t *R, const ecpt128_t *P, UV k,
                       uint128_t A24, const mont128_t *ctx) {
  if (k == 1) { *R = *P; return; }
  ecpt128_t R0 = *P, R1;
  ecm_double128(&R1, &R0, A24, ctx);
  if (k == 2) { *R = R1; return; }
  UV bit = (UV)1 << (8*sizeof(UV) - 1);
  while (!(bit & k)) bit >>= 1;
  for (bit >>= 1; bit; bit >>= 1) {
    if (k & bit) {
      ecm_dadd128(&R0, &R0, &R1, P, ctx);
      ecm_double128(&R1, &R1, A24, ctx);
    } else {
      ecm_dadd128(&R1, &R0, &R1, P, ctx);
      ecm_double128(&R0, &R0, A24, ctx);
    }
  }
  *R = R0;
}

/* Fixed σ values for Suyama's parameterization.
 * All σ ≥ 11 so u = σ²-5 > v = 4σ (non-degenerate).
 * All prime factors of 16u³v are ≤ 2011 (trial-division limit),
 * so gcd(denominator, n) = 1 is guaranteed after our trial-division step. */
static const uint16_t ecm_sigmas[] = {
   11,   13,   17,   19,   23,   29,   31,   37,   41,   43,
   47,   53,   59,   61,   67,   71,   73,   79,   83,   89,
  103,  127,  139,  149,  151,  157,  163,  181,  191,  197,
  199,  211,  223,  227,  233,  239,  257,  271,  293,  313, /*  40 */
  331,  337,  347,  359,  379,  389,  397,  401,  409,  421,
  443,  449,  457,  479,  487,  509,  521,  523,  547,  557,
  587,  641,  653,  659,  673,  677,  683,  691,  719,  727,
  739,  751,  769,  797,  809,  853,  919,  929,  941,  997,
 1049, 1051, 1063, 1091, 1093, 1109, 1117, 1129, 1153, 1201,
 1217, 1229, 1283, 1327, 1361, 1381, 1427, 1447, 1459, 1471, /* 100 */
 1481, 1489, 1543, 1549, 1571, 1621, 1709, 1723, 1753, 1759,
 1801, 1811, 1867, 1987, 2039, 2099, 2113, 2131, 2251, 2309,
 2347, 2381, 2399, 2447, 2473, 2551, 2557, 2663, 2677, 2689,
 2713, 2719, 2749, 2857, 2879, 2887, 2939, 3001, 3061, 3067,
 3121, 3137, 3187, 3251, 3259, 3271, 3307, 3359, 3371, 3373,
 3467, 3593, 3607, 3623, 3643, 3709, 3733, 3793, 3851, 3923, /* 160 */
 3989, 4019, 4049, 4129, 4231, 4253, 4283, 4339, 4349, 4441,
 4523, 4649, 4787, 4987, 4999, 5171, 5237, 5273, 5297, 5333,
 5387, 5471, 5479, 5647, 5749, 5791, 6101, 6163, 6257, 6299,
 6337, 6451, 6491, 6659, 6793, 6823, 6967, 7013, 7229, 7253, /* 200 */
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
36919,37363,37517,37699,37907,38047,38177,38273,38749,38903  /* 400 */
};
#define NECM128_SIGMAS ((int)(sizeof(ecm_sigmas)/sizeof(ecm_sigmas[0])))

/* Returns a non-trivial factor of n, or 0.
 * ncurves curves, Suyama parameterization, stage-1 bound B1. */
static uint128_t tinyecm128(uint128_t n, UV B1, int ncurves, int sigma_offset) {
  mont128_t ctx;
  mont_setup128(&ctx, n);
  UV sqrtB1 = (UV)sqrt((double)B1);

  for (int ci = 0; ci < ncurves && sigma_offset+ci < NECM128_SIGMAS; ci++) {
    UV sigma = ecm_sigmas[sigma_offset+ci];

    /* Suyama parameterization. */
    UV ui          = sigma * sigma - 5;       /* u = σ²-5  */
    UV vi          = 4 * sigma;               /* v = 4σ    */
    uint128_t u3   = (uint128_t)ui * ui * ui; /* u³ */
    uint128_t v3   = (uint128_t)vi * vi * vi; /* v³ */
    uint128_t umv  = ui - vi;                 /* u-v > 0 for σ ≥ 6 */
    uint128_t umv3 = umv * umv * umv;
    uint128_t t3uv = 3*ui + vi;

    /* A24 = (v-u)³(3u+v) / (16u³v) = -(u-v)³(3u+v) / (16u³v)  mod n */
    uint128_t abs_num = umv3 * t3uv;
    uint128_t abs_num_r = abs_num % n;
    uint128_t num = abs_num_r ? n - abs_num_r : 0;

    uint128_t den     = 16 * u3 * vi;
    uint128_t den_r   = den % n;
    uint128_t den_inv = modinv128(den_r, n);
    if (den_inv == 0) continue;   /* gcd(den,n) > 1 — shouldn't happen */

    uint128_t mnum  = mont_enter128(num,     &ctx);   /* num * R mod n */
    uint128_t mdinv = mont_enter128(den_inv, &ctx);   /* den_inv * R mod n */
    uint128_t mA24  = mont_mulmod128(mnum, mdinv, &ctx);

    ecpt128_t P    = { mont_enter128(u3, &ctx),
                       mont_enter128(v3, &ctx) };

    /* Stage 1: apply prime powers ≤ B1; check gcd(Z, n) every 64 primes. */
    uint32_t j = 0;
    START_DO_FOR_EACH_PRIME(2,B1) {
      UV k = p;
      if (p <= sqrtB1) { UV pm = B1 / p; while (k <= pm) k *= p; }
      ecm_mul128(&P, &P, k, mA24, &ctx);
      if ((j++ % 64) == 0) {
        uint128_t g = gcd128(P.Z, n);
        if (g > 1 && g < n) RETURN_FROM_EACH_PRIME(return g);
        if (g == n)         RETURN_FROM_EACH_PRIME(goto next_curve);
      }
    } END_DO_FOR_EACH_PRIME
    { uint128_t g = gcd128(P.Z, n);
      if (g > 1 && g < n) return g; }

next_curve:;
  }
  return 0;
}

#if 0
/*****************************************************************************
 * Brent's rho factoring
 *****************************************************************************/

/* Returns a non-trivial factor of n, or 0 on failure.
 * c is the polynomial constant x^2+c; try several values of c on failure.
 * Uses Montgomery multiplication for speed. */
static uint128_t brent128(uint128_t n, uint128_t c, UV maxr) {
  mont128_t ctx;
  mont_setup128(&ctx, n);

  UV m = 128;
  uint128_t r = 1, p;

  /* Convert to Montgomery form */
  uint128_t cy = mont_enter128(2, &ctx);   /* y in Montgomery form */
  uint128_t cc = mont_enter128(c, &ctx);   /* c in Montgomery form */
  uint128_t cq = mont_enter128(1, &ctx);   /* accumulator, starts at 1 */
  uint128_t cx, cys;

#define MBRENT_F(v)  addmod128(mont_sqrmod128(v, &ctx), cc, n)

  cx = cy;
  do {
    cx = cy;
    for (UV i = 0; i < r; i++)
      cy = MBRENT_F(cy);
    UV k = 0;
    do {
      cys = cy;
      for (UV i = 0; i < (m < r - k ? m : r - k); i++) {
        cy = MBRENT_F(cy);
        uint128_t diff = (cx > cy) ? cx - cy : cy - cx;
        cq = mont_mulmod128(cq, diff, &ctx);
      }
      p = gcd128(cq, n);   /* gcd(q_mont,n)=gcd(q_value,n) since gcd(R,n)=1 */
      k += m;
    } while (k < r && p == 1);
    if (r > maxr) return 0;
    r <<= 1;
  } while (p == 1);

  if (p == n) {
    do {
      cys = MBRENT_F(cys);
      uint128_t diff = (cx > cys) ? cx - cys : cys - cx;
      p = gcd128(diff, n);
    } while (p == 1);
  }
#undef MBRENT_F

  return (p == n) ? 0 : p;
}
#endif

/*****************************************************************************
 * Result management
 *****************************************************************************/

/* Insert a UV prime factor p into nf (sorted).  Increments exponent if
 * already present.  Caller must ensure nfactors < MPU_MAX_DFACTORS128. */
static void insert_factor128(factored128_t *nf, UV p) {
  uint16_t i, nf_ = nf->nfactors;
  for (i = 0; i < nf_; i++) {
    if (nf->f[i] == p) { nf->e[i]++; return; }
    if (nf->f[i] >  p) break;
  }
  /* Shift right and insert */
  for (uint16_t j = nf_; j > i; j--) {
    nf->f[j] = nf->f[j-1];
    nf->e[j] = nf->e[j-1];
  }
  nf->f[i] = p;
  nf->e[i] = 1;
  nf->nfactors++;
}

/*****************************************************************************
 * Public entry point
 *****************************************************************************/

bool factorintp128(factored128_t *nf, uint128_t n) {
  int const verbose = _XS_get_verbose();

  nf->n        = n;
  nf->nfactors = 0;
  nf->flarge   = 0;

  if (n <= 1) return 0;

  /* Trial division using primes_small[] (primes 2..2011) */
  for (int i = 1; i < NPRIMES_SMALL; i++) {
    UV p = primes_small[i];
    if ((uint128_t)p * p > n) break;
    while (n % p == 0) {
      insert_factor128(nf, p);
      n /= p;
    }
  }
  if (n == 1) return 1;

  /* n is now > 2011; remaining cofactor stack */
  uint128_t stack[128];
  int top = 0;
  stack[top++] = n;

  while (top > 0) {
    uint128_t t = stack[--top];
    int show = verbose;
    /* printf("  factor128 working on %s\n", u128_str(t)); */

    /* Shortcut: fits in UV */
    if (t <= (uint128_t)UV_MAX) {
      UV uv = (UV)t;
      UV factors[MPU_MAX_FACTORS];
      int nf2 = factor(uv, factors);
      for (int i = 0; i < nf2; i++)
        insert_factor128(nf, factors[i]);
      continue;
    }

    if (is_bpsw128(t)) {
      /* Allow for shortcut above to be adjusted to allow 64-bit inputs.
       * If 64-bit, insert the prime.  Otherwise it's the ONE large prime. */
      if (t <= (uint128_t)UV_MAX) insert_factor128(nf, t);
      else                        nf->flarge = t;
      continue;
    }

    /* Composite: find a factor */
    uint128_t f = 0;

    /* P-1 stage 1+2: fast first try */
    if (!f) f = pminus1_128(t, 20000, 200000);
    if (f && show) {show=0;printf("p-1/128 found factor %s\n",u128_str(f));}

    /* SQUFOF: fast for balanced composites.
     * rounds budget scales with n's bit size: each multiplier needs
     * ~0.5*(n*mult)^(1/5) iters, so we allow ~10x that per bit-width
     * to cover enough multipliers before giving up.
     * At 79+ bits, skipping this for ECM is faster. */
    if (!f && ((t >> 60) >> 19) == 0) {
      int b128 = 0;
      { uint128_t tmp = t; while (tmp >>= 1) b128++; }
      UV sq_rounds = (UV)(10.0 * pow(2.0, (double)b128 / 5.0));
      if (sq_rounds < 300000)  sq_rounds = 300000;
      if (sq_rounds > 5000000) sq_rounds = 5000000;
      f = squfof128(t, sq_rounds);
      if (f && show) {show=0;printf("squfof128 found factor %s\n",u128_str(f));}
    }

    /* ECM */
    if (!f) f = tinyecm128(t, 2000, 40, 0);
    if (f && show) {show=0;printf("tinyecm128 2k found factor %s\n",u128_str(f));}

    if (!f) f = tinyecm128(t, 10000, 40, 40);
    if (f && show) {show=0;printf("tinyecm128 10k found factor %s\n",u128_str(f));}

    if (!f) f = tinyecm128(t, 50000, 100, 80);
    if (f && show) {show=0;printf("tinyecm128 50k found factor %s\n",u128_str(f));}

    /* At this point the smallest factor is highly likely to be >= 54 bits */
    if (!f) f = tinyecm128(t, 800000, 140, 180);
    if (f && show) {show=0;printf("tinyecm128 800k found factor %s\n",u128_str(f));}

    /* It is unlikely a 128-bit semiprime will get here at all. */
    if (!f) f = tinyecm128(t, 160000000, 4, 320);
    if (f && show) {show=0;printf("tinyecm128 160M found factor %s\n",u128_str(f));}

#if 0
    /* Pollard/Brent rho: fallback when others failed */
    for (uint32_t c = 1; c < 20 && f == 0; c += 2) {
      if (show)printf("  brent128 try c=%u\n",c);
      f = brent128(t, c, 1U << 20);
    }
    if (f && show) {show=0;printf("brent128 found factor %s\n",u128_str(f));}
#endif

    if (f == 0) {
      if (verbose) printf("failed to find factor in factor128: %s\n",u128_str(n));
      return 0;
    }

    /* Push both halves (f and t/f) for further factoring */
    stack[top++] = f;
    stack[top++] = t / f;
  }
  if (verbose) fflush(stdout);
  return 1;
}

#endif /* BITS_PER_WORD == 64 && HAVE_UINT128 */
