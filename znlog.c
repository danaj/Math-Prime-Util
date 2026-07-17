#include <stdio.h>

#define FUNC_ipow   1
#define FUNC_isqrt  1
#define FUNC_gcd_ui 1
#include "ptypes.h"
#include "factor.h"
#include "util.h"
#include "mulmod.h"
#include "primality.h"
#include "znlog.h"

#if BITS_PER_WORD == 64 && HAVE_UINT64 && HAVE_UINT128
  #define USE_ZNLOG_MONTMATH 1

/* Pollard rho spends almost all its time multiplying with one fixed modulus.
 * Keep its state in Montgomery form to avoid a 128-bit remainder each step.
 * R is 2^64 and all inputs to znlog_mont_mul are in [0,n). */
typedef struct {
  uint64_t n, ninv, r2, one;
} znlog_mont_t;

static INLINE uint64_t znlog_mont_mul(uint64_t a, uint64_t b,
                                      const znlog_mont_t *ctx) {
  uint128_t t = (uint128_t)a * b;
  uint64_t lo = (uint64_t)t;
  uint64_t m = lo * ctx->ninv;
  uint128_t mn = (uint128_t)m * ctx->n;
  uint64_t th = (uint64_t)(t >> 64);
  uint64_t mh = (uint64_t)(mn >> 64);
  uint64_t r = th + mh;
  int overflow = (r < th);
  uint64_t oldr = r;

  /* lo + low(m*n) is either 0 or 2^64.  Track that carry and any
   * overflow from the high-half addition before the final subtraction. */
  r += (lo != 0);
  overflow |= (r < oldr);
  return (overflow || r >= ctx->n) ? r - ctx->n : r;
}

static void znlog_mont_setup(znlog_mont_t *ctx, uint64_t n) {
  uint64_t i, x = (3*n) ^ 2;

  ctx->n = n;
  /* Hensel lifting doubles the correct bits of n^-1 at each step. */
  x *= (uint64_t)2 - n*x;
  x *= (uint64_t)2 - n*x;
  x *= (uint64_t)2 - n*x;
  x *= (uint64_t)2 - n*x;
  ctx->ninv = (uint64_t)0 - x;

  /* r2 = R^2 mod n.  Setup is infrequent, so avoid division entirely. */
  ctx->r2 = 1;
  for (i = 0; i < 128; i++)
    ctx->r2 = addmod(ctx->r2, ctx->r2, n);
  ctx->one = znlog_mont_mul(1, ctx->r2, ctx);
}

static INLINE uint64_t znlog_mont_enter(uint64_t a,
                                        const znlog_mont_t *ctx) {
  return znlog_mont_mul(a, ctx->r2, ctx);
}

static uint64_t znlog_mont_pow(uint64_t a, uint64_t k,
                               const znlog_mont_t *ctx) {
  uint64_t r = ctx->one;
  while (k) {
    if (k & 1) r = znlog_mont_mul(r, a, ctx);
    k >>= 1;
    if (k) a = znlog_mont_mul(a, a, ctx);
  }
  return r;
}
#else
  #define USE_ZNLOG_MONTMATH 0
#endif

/******************************************************************************/
/* DLP */
/******************************************************************************/

static UV dlp_trial(UV a, UV g, UV p, UV maxrounds) {
  UV k, t;
  if (maxrounds >= p) maxrounds = p-1;

#if USE_ZNLOG_MONTMATH
  if ((p&1) && maxrounds > 256) {
    znlog_mont_t ctx;
    UV stop = 128;

    /* Avoid Montgomery setup when the log is found almost immediately. */
    for (t = g, k = 1; k < stop; k++) {
      if (t == a)
        return k;
      t = mulmod(t, g, p);
      if (t == g) return 0;   /* Stop at cycle */
    }

    znlog_mont_setup(&ctx, p);
    t = znlog_mont_enter(t, &ctx);
    g = znlog_mont_enter(g, &ctx);
    a = znlog_mont_enter(a, &ctx);
    for ( ; k < maxrounds; k++) {
      if (t == a)
        return k;
      t = znlog_mont_mul(t, g, &ctx);
      if (t == g) break;   /* Stop at cycle */
    }
  } else
#endif
  {
    for (t = g, k = 1; k < maxrounds; k++) {
      if (t == a)
        return k;
      t = mulmod(t, g, p);
      if (t == g) break;   /* Stop at cycle */
    }
  }
  return 0;
}

/******************************************************************************/
/* DLP - Pollard Rho */
/******************************************************************************/

/* Compare with Pomerance paper (dartmouth dtalk4):
 * Type I/II/III = our case 1, 0, 2.
 * x_i = u, a_i = v, b_i = w
 *
 * Also see Bai/Brent 2008 for many ideas to speed this up.
 * https://maths-people.anu.edu.au/~brent/pd/rpb231.pdf
 * E.g. Teske adding-walk, Brent's cycle algo, Teske modified cycle
 */
#define pollard_rho_cycle(u,v,w,p,n,a,g) \
    switch (u % 3) { \
      case 0: u = mulmod(u,u,p);  v = addmod(v,v,n);  w = addmod(w,w,n); break;\
      case 1: u = mulmod(u,a,p);  v = addmod(v,1,n);                     break;\
      case 2: u = mulmod(u,g,p);                      w = addmod(w,1,n); break;\
    }

#if USE_ZNLOG_MONTMATH
#define pollard_rho_cycle_mont(u,v,w,n,a,g,ctx) \
    switch (u % 3) { \
      case 0: u = znlog_mont_mul(u,u,ctx);  v = addmod(v,v,n);  w = addmod(w,w,n); break;\
      case 1: u = znlog_mont_mul(u,a,ctx);  v = addmod(v,1,n);                     break;\
      case 2: u = znlog_mont_mul(u,g,ctx);                      w = addmod(w,1,n); break;\
    }

static UV dlp_prho_mont(UV a, UV g, UV p, UV n, UV maxrounds) {
  znlog_mont_t ctx;
  int const verbose = _XS_get_verbose();
  UV k = 0, round;
  UV am, gm;
  int nrestarts;

  znlog_mont_setup(&ctx, p);
  am = znlog_mont_enter(a, &ctx);
  gm = znlog_mont_enter(g, &ctx);

  for (nrestarts = 0; nrestarts <= 4; nrestarts++) {
    UV u, v, w, U, V, W, power, lam;

    u = U = (nrestarts == 0) ? ctx.one
                             : znlog_mont_pow(gm, (UV)nrestarts, &ctx);
    v = V = 0;
    w = W = (UV)nrestarts;
    power = lam = 1;

    pollard_rho_cycle_mont(U,V,W,n,am,gm,&ctx);
    round = 1;

    while (u != U && round < maxrounds) {
      if (power == lam) {
        u = U;  v = V;  w = W;
        power <<= 1;
        lam = 0;
      }
      pollard_rho_cycle_mont(U,V,W,n,am,gm,&ctx);
      lam++;
      round++;
    }

    if (u == U) {
      UV r1, r2, G, G2;
      r1 = submod(v, V, n);
      if (r1 == 0) {
        if (nrestarts >= 4) {
          if (verbose) printf("DLP Rho failure, r=0 (no more restarts)\n");
          break;
        }
        if (verbose) printf("DLP Rho r1=0, restart %d\n", nrestarts+1);
        continue;
      }

      r2 = submod(W, w, n);
      G = gcd_ui(r1,n);
      G2 = gcd_ui(G,r2);
      k = divmod(r2/G2, r1/G2, n/G2);
      if (G > 1) {
        if (powmod(g,k,p) == a) {
          if (verbose > 2) printf("  common GCD %"UVuf"\n", G2);
        } else {
          UV m, l = divmod(r2, r1, n/G);
          for (m = 0; m < G; m++) {
            k = addmod(l, mulmod(m,(n/G),n), n);
            if (powmod(g,k,p) == a) break;
          }
          if (m<G && verbose > 2) printf("  GCD %"UVuf", found with m=%"UVuf"\n", G, m);
        }
      }

      if (powmod(g,k,p) == a)
        break;

      if (verbose > 2) printf("r1 = %"UVuf"  r2 = %"UVuf" k = %"UVuf"\n", r1, r2, k);
      if (verbose) printf("Incorrect DLP Rho solution: %"UVuf"\n", k);
      k = 0;
      break;
    }
  }

  if (verbose && k) printf("DLP Rho solution found after %"UVuf" steps\n", round);
  return k;
}
#endif

static UV dlp_prho(UV a, UV g, UV p, UV n, UV maxrounds) {
  int const verbose = _XS_get_verbose();
  UV k = 0, round;
  int nrestarts;

  if (maxrounds > n) maxrounds = n;

#if USE_ZNLOG_MONTMATH
  /* Setup costs more than the complete search for very small orders. */
  if ((p & 1) && maxrounds > 256)
    return dlp_prho_mont(a, g, p, n, maxrounds);
#endif

  for (nrestarts = 0; nrestarts <= 4; nrestarts++) {
    UV u, v, w, U, V, W, power, lam;

    u = U = (nrestarts == 0) ? 1 : powmod(g, (UV)nrestarts, p);
    v = V = 0;
    w = W = (UV)nrestarts;
    power = lam = 1;

    pollard_rho_cycle(U,V,W,p,n,a,g);
    round = 1;

    while (u != U && round < maxrounds) {
      if (power == lam) {
        u = U;  v = V;  w = W;
        power <<= 1;
        lam = 0;
      }
      pollard_rho_cycle(U,V,W,p,n,a,g);
      lam++;
      round++;
    }

    if (u == U) {
      UV r1, r2, G, G2;
      r1 = submod(v, V, n);
      if (r1 == 0) {
        if (nrestarts >= 4) {
          if (verbose) printf("DLP Rho failure, r=0 (no more restarts)\n");
          break;
        }
        if (verbose) printf("DLP Rho r1=0, restart %d\n", nrestarts+1);
        continue;
      }

      r2 = submod(W, w, n);
      G = gcd_ui(r1,n);
      G2 = gcd_ui(G,r2);
      k = divmod(r2/G2, r1/G2, n/G2);
      if (G > 1) {
        if (powmod(g,k,p) == a) {
          if (verbose > 2) printf("  common GCD %"UVuf"\n", G2);
        } else {
          UV m, l = divmod(r2, r1, n/G);
          for (m = 0; m < G; m++) {
            k = addmod(l, mulmod(m,(n/G),n), n);
            if (powmod(g,k,p) == a) break;
          }
          if (m<G && verbose > 2) printf("  GCD %"UVuf", found with m=%"UVuf"\n", G, m);
        }
      }

      if (powmod(g,k,p) == a)
        break;

      if (verbose > 2) printf("r1 = %"UVuf"  r2 = %"UVuf" k = %"UVuf"\n", r1, r2, k);
      if (verbose) printf("Incorrect DLP Rho solution: %"UVuf"\n", k);
      k = 0;
      break;
    }
  }

  if (verbose && k) printf("DLP Rho solution found after %"UVuf" steps\n", round);
  return k;
}


/******************************************************************************/
/* DLP - BSGS */
/******************************************************************************/

typedef struct bsgs_hash_t {
  UV M;    /* The baby step index */
  UV V;    /* The powmod value */
  struct bsgs_hash_t* next;
} bsgs_hash_t;

/****************************************/
/*  Simple and limited pool allocation  */
#define BSGS_ENTRIES_PER_PAGE 8000
typedef struct bsgs_page_top_t {
  struct bsgs_page_t* first;
  bsgs_hash_t** table;
  UV  size;
  int nused;
  int npages;
} bsgs_page_top_t;

typedef struct bsgs_page_t {
  bsgs_hash_t entries[BSGS_ENTRIES_PER_PAGE];
  struct bsgs_page_t* next;
} bsgs_page_t;

static bsgs_hash_t* get_entry(bsgs_page_top_t* top) {
  if (top->nused == 0 || top->nused >= BSGS_ENTRIES_PER_PAGE) {
    bsgs_page_t* newpage;
    Newz(0, newpage, 1, bsgs_page_t);
    newpage->next = top->first;
    top->first = newpage;
    top->nused = 0;
    top->npages++;
  }
  return top->first->entries + top->nused++;
}
static void destroy_pages(bsgs_page_top_t* top) {
  bsgs_page_t* head = top->first;
  while (head != 0) {
    bsgs_page_t* next = head->next;
    Safefree(head);
    head = next;
  }
  top->first = 0;
}
/****************************************/

static void bsgs_hash_put(bsgs_page_top_t* pagetop, UV v, UV i) {
  UV idx = v % pagetop->size;
  bsgs_hash_t** table = pagetop->table;
  bsgs_hash_t* entry = table[idx];

  while (entry && entry->V != v)
    entry = entry->next;

  if (!entry) {
    entry = get_entry(pagetop);
    entry->M = i;
    entry->V = v;
    entry->next = table[idx];
    table[idx] = entry;
  }
}

static UV bsgs_hash_get(bsgs_page_top_t* pagetop, UV v) {
  bsgs_hash_t* entry = pagetop->table[v % pagetop->size];
  while (entry && entry->V != v)
    entry = entry->next;
  return (entry) ? entry->M : 0;
}

static UV bsgs_hash_put_get(bsgs_page_top_t* pagetop, UV v, UV i) {
  UV idx = v % pagetop->size;
  bsgs_hash_t** table = pagetop->table;
  bsgs_hash_t* entry = table[idx];

  while (entry && entry->V != v)
    entry = entry->next;

  if (entry)
    return entry->M;

  entry = get_entry(pagetop);
  entry->M = i;
  entry->V = v;
  entry->next = table[idx];
  table[idx] = entry;
  return 0;
}

static UV dlp_bsgs(UV a, UV g, UV p, UV n, UV maxent) {
  bsgs_page_top_t PAGES;
  UV i, m, maxm, hashmap_count;
  UV aa, S, gm, T, gs_i, bs_i;
  UV result = 0;
  int const verbose = _XS_get_verbose();

  if (n <= 2) return 0;   /* Shouldn't be here with gorder this low */

  if (a == 0) return 0;  /* We don't handle this case */

  maxm = isqrt(n);
  m = (maxent > maxm) ? maxm : maxent;
  if (m > 40000000) m = 40000000;  /* cap to match hash table capacity */

  hashmap_count = (m < 65537) ? 65537 : next_prime(m);  /* Ave depth around 2 */

  /* Create table.  Size: 8*hashmap_count bytes. */
  PAGES.size = hashmap_count;
  PAGES.first = 0;
  PAGES.nused = 0;
  PAGES.npages = 0;
  Newz(0, PAGES.table, hashmap_count, bsgs_hash_t*);

  aa = mulmod(a,a,p);
  S = a;
  gm = powmod(g, m, p);
  T = gm;
  gs_i = 0;
  bs_i = 0;

  bsgs_hash_put(&PAGES, S, 0);   /* First baby step */
  S = mulmod(S, g, p);
  /* Interleaved Baby Step Giant Step */
  for (i = 1; i <= m; i++) {
    gs_i = bsgs_hash_put_get(&PAGES, S, i);
    if (gs_i) { bs_i = i; break; }
    S = mulmod(S, g, p);
    if (S == aa) {  /* S = a*g^(i+1).  If S == a^2 then g^(i+1) == a. */
      if (verbose) printf("  dlp bsgs: solution at BS step %"UVuf"\n", i+1);
      result = i+1;
      break;
    }
    bs_i = bsgs_hash_put_get(&PAGES, T, i);
    if (bs_i) { gs_i = i; break; }
    T = mulmod(T, gm, p);
  }

  if (!result) {
    /* Extend Giant Step search */
    if (!(gs_i || bs_i)) {
      UV b = n/m + ((n % m) != 0);
      if (m < maxm && b > 8*m) b = 8*m;
      for (i = m+1; i < b; i++) {
        bs_i = bsgs_hash_get(&PAGES, T);
        if (bs_i) { gs_i = i; break; }
        T = mulmod(T, gm, p);
      }
    }

    if (gs_i || bs_i) {
      result = submod(mulmod(gs_i, m, n), bs_i, n);
    }
  }
  if (verbose) printf("  dlp bsgs using %d pages (%.1fMB+%.1fMB) for hash\n", PAGES.npages, ((double)PAGES.npages * sizeof(bsgs_page_t)) / (1024*1024), ((double)hashmap_count * sizeof(bsgs_hash_t*)) / (1024*1024));

  destroy_pages(&PAGES);
  Safefree(PAGES.table);
  if (result != 0 && powmod(g,result,p) != a) {
    if (verbose) printf("Incorrect DLP BSGS solution: %"UVuf"\n", result);
    result = 0;
  }
  return result;
}


/* Find smallest k where a = g^k mod p */
#define DLP_TRIAL_NUM  20000
UV znlog_solve(UV a, UV g, UV p, UV n) {
  UV k;
  const int verbose = _XS_get_verbose();

  if (a >= p) a %= p;
  if (g >= p) g %= p;

  if (a == 1 || g == 0 || p <= 2)
    return 0;

  if (verbose > 1 && n != p-1) printf("  g=%"UVuf" p=%"UVuf", order %"UVuf"\n", g, p, n);

  /* n = 0 means p is composite AND g is not coprime to p.
   * There may still be solutions, e.g. znlog(12,98,100) = 19
   * But most of the algorithms for DLP don't work in this case.
   */
  if (n == 0) {
    k = dlp_trial(a, g, p, p-1);
    if (verbose) printf("  dlp full trial %s\n", k!=0 ? "success" : "failure");
    return k;
  }

  /* If n = znorder(g,p) is small, then do a simple search */
  if (n <= DLP_TRIAL_NUM) {
    k = dlp_prho(a, g, p, n, isqrt(n)*4);
    if (k != 0) {
      if (verbose) printf("  dlp prho small success\n");
      return k;
    }
    k = dlp_trial(a, g, p, DLP_TRIAL_NUM);
    if (verbose) printf("  dlp trial 10k %s\n", (k!=0 || p <= DLP_TRIAL_NUM) ? "success" : "failure");
    return k;
  }

  /* Existence check, n != 0 here, meaning g is coprime to p. */
  if (powmod(a, n /* n is gorder */, p) != 1) return 0;

  k = dlp_trial(a, g, p, DLP_TRIAL_NUM);
  if (k != 0) {
    if (verbose) printf("  dlp trial 20k success\n");
    return k;
  }

  {
    UV maxent = (isqrt(n)+1) / 2;
    if (maxent > 1000000) {
      k = dlp_prho(a, g, p, n, 2000000000U);
      if (verbose) printf("  dlp rho large %s\n", k!=0 ? "success" : "failure");
    } else {
      k = dlp_bsgs(a, g, p, n, maxent);
      if (verbose) printf("  dlp bsgs %"UVuf"k %s\n", maxent/1000, k!=0 ? "success" : "failure");
    }
    if (k != 0) return k;
  }

  if (verbose) printf("  dlp doing exhaustive trial\n");
  k = dlp_trial(a, g, p, p-1);
  return k;
}

/* Find k in [0,q^e) with gamma^k = delta mod n, where gamma has order q^e.
 * This is Pohlig-Hellman digit extraction: e subproblems of order q.
 */
static int znlog_prime_power(UV *r, UV delta, UV gamma, UV n, UV q, uint32_t e) {
  UV pe, pj, gamma0, k;
  uint32_t j;

  pe = ipow(q, e);
  gamma0 = powmod(gamma, ipow(q, e-1), n);
  k = 0;
  pj = 1;

  for (j = 0; j < e; j++) {
    UV t, h, xj;

    if (k == 0) {
      t = delta;
    } else {
      t = mulmod(delta, powmod(gamma, pe-k, n), n);
    }

    h = (j == e-1) ? t : powmod(t, ipow(q, e-1-j), n);

    if (h == 1) {
      xj = 0;
    } else {
      xj = znlog_solve(h, gamma0, n, q);
      if (xj == 0)
        return 0;
      if (xj >= q)
        xj %= q;
    }

    k += xj * pj;
    if (j != e-1)
      pj *= q;
  }

  if (powmod(gamma, k, n) != delta)
    return 0;

  *r = k;
  return 1;
}

/* Silver-Pohlig-Hellman */
static UV znlog_ph(UV a, UV g, UV p, UV p1) {
  factored_t pf;
  UV x, sol[MPU_MAX_DFACTORS], mod[MPU_MAX_DFACTORS];
  uint32_t i;

  if (p1 == 0) return 0;
  pf = factorint(p1);
  if (pf.nfactors == 1 && pf.e[0] == 1)
    return znlog_solve(a, g, p, p1);
  for (i = 0; i < pf.nfactors; i++) {
    UV pi = ipow(pf.f[i],pf.e[i]);
    UV delta = powmod(a, p1/pi, p);
    UV gamma = powmod(g, p1/pi, p);
    /* printf(" solving znlog(%"UVuf",%"UVuf",%"UVuf")\n", delta, gamma, p); */
    /* Does a solution for this sub-problem exist? */
    if (delta != 1 && powmod(delta, pi, p) != 1) return 0;
    if (delta == 1) {
      sol[i] = 0;
    } else if (!znlog_prime_power(&sol[i], delta, gamma, p, pf.f[i], pf.e[i])) {
      return 0;
    }
    mod[i] = pi;
  }
  if (chinese(&x, 0, sol, mod, pf.nfactors) == 1 && powmod(g, x, p) == a)
    return x;
  return 0;
}

/* Reduce non-coprime logs until the remaining base is coprime to the
 * remaining modulus.  Return 0 if no solution is found.
 */
static UV znlog_non_coprime(UV a, UV g, UV n) {
  UV d, nr, ar, gp, dacc, k, gorig;

  gorig = g;
  d = gcd_ui(g, n);
  if (d == 1)
    return 0;

  nr = n;
  ar = a;
  gp = 1;
  dacc = 1;
  k = 0;

  while (d > 1) {
    if (gp == a)
      return k;
    if (ar % d != 0)
      return 0;

    nr /= d;
    ar /= d;

    if (nr > 1)
      dacc = mulmod(dacc, gorig/d, nr);
    else
      dacc = 0;

    gp = mulmod(gp, gorig, n);
    k++;

    d = gcd_ui(g, nr);
  }

  if (gp == a)
    return k;
  if (nr == 1)
    return 0;

  d = modinverse(dacc, nr);
  if (d == 0)
    return 0;

  ar = mulmod(ar % nr, d, nr);
  g = gorig % nr;
  d = znlog(ar, g, nr);
  if (d == 0 && ar != 1)
    return 0;

  k += d;
  return (powmod(gorig, k, n) == a) ? k : 0;
}

/* Find smallest k where a = g^k mod p */
UV znlog(UV a, UV g, UV p) {
  UV k, gorder;
  const int verbose = _XS_get_verbose();

  if (a >= p) a %= p;
  if (g >= p) g %= p;

  if (a == 1 || g == 0 || p <= 2)
    return 0;

  gorder = znorder(g,p);
  if (gorder != 0 && powmod(a, gorder, p) != 1) return 0;

  if (gorder == 0) {
    k = znlog_non_coprime(a, g, p);
    if (verbose) printf("  dlp non-coprime %s\n", k!=0 ? "success" : "failure");
    return k;
  }

  if (gorder != 0 && !is_prob_prime(gorder)) {
    k = znlog_ph(a, g, p, gorder);
    if (verbose) printf("  dlp PH %s\n", k!=0 ? "success" : "failure");
    return k;  /* solution or 0 (no solution exists) */
  }

  return znlog_solve(a, g, p, gorder);
}
