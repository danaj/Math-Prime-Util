#include <math.h>
#include <string.h>
#include "ptypes.h"
#include "tinyqs128.h"
#include "util.h"
#include "util_bits.h"

/*******************************************************************************
 *
 *             tinyqs128
 *
 * uint128_t adaptation of Jason Papadopoulos's tinyqs (SIQS).
 *
 * Attempts to factor integers up to 128 bits.
 *
 * Ideas from the MPQS code of F. Bahr, and lattice sievers of Jens Franke.
 *
 ******************************************************************************/

#if HAVE_FACTOR128

typedef  uint8_t    u8;
typedef  int16_t   s16;
typedef uint16_t   u16;
typedef  int32_t   s32;
typedef uint32_t   u32;
typedef uint64_t   u64;

typedef struct {
  uint128_t mag;
  bool neg;
} tq_cell;

/* This is a fixed-width scalar and never allocates. */
typedef tq_cell tq_t[1];

#define TQ_U128_MAX ((uint128_t)-1)

static void tq_normalize(tq_t a)
{
  if (a->mag == 0) a->neg = 0;
}

static void tq_init(tq_t a)
{
  a->mag = 0;
  a->neg = 0;
}

static void tq_set(tq_t r, const tq_t a)
{
  *r = *a;
}

static void tq_set_ui(tq_t r, uint128_t a)
{
  r->mag = a;
  r->neg = 0;
}

static uint64_t tq_get_ui(const tq_t a)
{
  return (uint64_t)a->mag;
}

static int tq_sgn(const tq_t a)
{
  return a->mag == 0 ? 0 : (a->neg ? -1 : 1);
}

static void tq_abs(tq_t r, const tq_t a)
{
  r->mag = a->mag;
  r->neg = 0;
}

static int tq_cmp(const tq_t a, const tq_t b)
{
  if (a->neg != b->neg) return a->neg ? -1 : 1;
  if (a->mag == b->mag) return 0;
  if (a->neg) return a->mag > b->mag ? -1 : 1;
  return a->mag > b->mag ? 1 : -1;
}

static int tq_cmp_ui(const tq_t a, uint128_t b)
{
  if (a->neg) return -1;
  return a->mag == b ? 0 : (a->mag > b ? 1 : -1);
}

static void tq_swap(tq_t a, tq_t b)
{
  tq_cell t = *a;
  *a = *b;
  *b = t;
}

static void tq_add(tq_t r, const tq_t a, const tq_t b)
{
  tq_cell t;

  if (a->neg == b->neg) {
    MPUassert(a->mag <= TQ_U128_MAX - b->mag,
              "tinyqs128 signed addition overflow");
    t.mag = a->mag + b->mag;
    t.neg = a->neg;
  } else if (a->mag >= b->mag) {
    t.mag = a->mag - b->mag;
    t.neg = a->neg;
  } else {
    t.mag = b->mag - a->mag;
    t.neg = b->neg;
  }
  *r = t;
  tq_normalize(r);
}

static void tq_sub(tq_t r, const tq_t a, const tq_t b)
{
  tq_cell nb = *b;
  if (nb.mag) nb.neg = !nb.neg;
  {
    tq_t t;
    t[0] = nb;
    tq_add(r, a, t);
  }
}

static void tq_mul(tq_t r, const tq_t a, const tq_t b)
{
  tq_cell t;

  if (a->mag == 0 || b->mag == 0) {
    tq_set_ui(r, 0);
    return;
  }
  MPUassert(a->mag <= TQ_U128_MAX / b->mag,
            "tinyqs128 signed multiplication overflow");
  t.mag = a->mag * b->mag;
  t.neg = a->neg != b->neg;
  *r = t;
}

static void tq_mul_ui(tq_t r, const tq_t a, uint64_t b)
{
  tq_cell t;

  if (a->mag == 0 || b == 0) {
    tq_set_ui(r, 0);
    return;
  }
  MPUassert(a->mag <= TQ_U128_MAX / b,
            "tinyqs128 small multiplication overflow");
  t.mag = a->mag * b;
  t.neg = a->neg;
  *r = t;
}

static uint64_t tq_tdiv_ui(const tq_t a, uint64_t d)
{
  MPUassert(d != 0, "tinyqs128 division by zero");
  return (uint64_t)(a->mag % d);
}

static uint64_t tq_tdiv_q_ui(tq_t q, const tq_t a, uint64_t d)
{
  uint64_t rem;
  MPUassert(d != 0, "tinyqs128 division by zero");
  rem = (uint64_t)(a->mag % d);
  q->mag = a->mag / d;
  q->neg = a->neg;
  tq_normalize(q);
  return rem;
}

static void tq_divexact_ui(tq_t q, const tq_t a, uint64_t d)
{
  MPUassert(d != 0 && a->mag % d == 0,
            "tinyqs128 inexact small division");
  q->mag = a->mag / d;
  q->neg = a->neg;
  tq_normalize(q);
}

static void tq_divexact(tq_t q, const tq_t a, const tq_t d)
{
  if (d->mag == 0 || a->mag % d->mag != 0) {
    char as[40], ds[40], rs[40];
    u128_to_str(as, a->mag);
    u128_to_str(ds, d->mag);
    u128_to_str(rs, d->mag ? a->mag % d->mag : 0);
    croak("tinyqs128 inexact division: %s / %s remainder %s",
          as, ds, rs);
  }
  q->mag = a->mag / d->mag;
  q->neg = (q->mag != 0) && (a->neg != d->neg);
}

static unsigned int tq_sizeinbase(const tq_t a, int base)
{
  uint128_t n = a->mag;
  unsigned int bits = 0;
  MPUassert(base == 2, "tinyqs128 only supports binary size");
  do { bits++; n >>= 1; } while (n);
  return bits;
}

static unsigned int tq_scan1(const tq_t a, unsigned int start)
{
  uint128_t n;
  unsigned int bit = start;
  MPUassert(a->mag != 0, "tinyqs128 cannot scan zero");
  n = a->mag >> start;
  while (!(n & 1)) { n >>= 1; bit++; }
  return bit;
}

static void tq_tdiv_q_2exp(tq_t q, const tq_t a, unsigned int bits)
{
  q->mag = bits >= 128 ? 0 : a->mag >> bits;
  q->neg = a->neg;
  tq_normalize(q);
}

static uint128_t tq_mod_value(const tq_t a, uint128_t n)
{
  uint128_t r;
  MPUassert(n != 0, "tinyqs128 modulus by zero");
  r = a->mag % n;
  if (a->neg && r) r = n - r;
  return r;
}

static uint128_t tq_gcd_value(uint128_t a, uint128_t b)
{
  while (b) {
    uint128_t t = b;
    b = a % b;
    a = t;
  }
  return a;
}

static void tq_gcd(tq_t r, const tq_t a, const tq_t b)
{
  tq_set_ui(r, tq_gcd_value(a->mag, b->mag));
}

static uint128_t tq_addmod_value(uint128_t a, uint128_t b, uint128_t n)
{
  uint128_t nb = n - b;
  return a >= nb ? a - nb : a + b;
}

static uint128_t tq_mulmod_value(uint128_t a, uint128_t b, uint128_t n)
{
  uint128_t r = 0;

  a %= n;
  b %= n;
  if (a <= (uint128_t)UINT64_MAX && b <= (uint128_t)UINT64_MAX)
    return ((uint128_t)(uint64_t)a * (uint64_t)b) % n;

  while (b) {
    if (b & 1) r = tq_addmod_value(r, a, n);
    b >>= 1;
    if (b) a = tq_addmod_value(a, a, n);
  }
  return r;
}

static void tq_mulmod(tq_t r, const tq_t a, const tq_t b, const tq_t n)
{
  uint128_t av = tq_mod_value(a, n->mag);
  uint128_t bv = tq_mod_value(b, n->mag);
  tq_set_ui(r, tq_mulmod_value(av, bv, n->mag));
}

static void tq_mulmod_ui(tq_t r, const tq_t a, uint64_t b, const tq_t n)
{
  uint128_t av = tq_mod_value(a, n->mag);
  tq_set_ui(r, tq_mulmod_value(av, b % n->mag, n->mag));
}

static uint64_t tq_isqrt_value(uint128_t n)
{
  uint128_t r = 1, t = n;

  if (n == 0) return 0;
  while (t > 3) { t >>= 2; r <<= 1; }
  r <<= 1;
  while (1) {
    uint128_t rn = (r + n / r) >> 1;
    if (rn >= r) break;
    r = rn;
  }
  while (r * r > n) r--;
  return (uint64_t)r;
}

static int tq_root(tq_t r, const tq_t a, unsigned int k)
{
  uint64_t root;
  MPUassert(k == 2 && !a->neg, "tinyqs128 only supports square roots");
  root = tq_isqrt_value(a->mag);
  tq_set_ui(r, root);
  return (uint128_t)root * root == a->mag;
}

static void tq_sqrt(tq_t r, const tq_t a)
{
  MPUassert(!a->neg, "tinyqs128 square root of negative value");
  tq_set_ui(r, tq_isqrt_value(a->mag));
}

static double tq_get_d(const tq_t a)
{
  double d = (double)a->mag;
  return a->neg ? -d : d;
}

/******************************************************************************/

#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

#ifndef M_SQRT2
#define M_SQRT2	1.41421356237309504880
#endif

#define RAND_MULT 2131995753

static u32 get_rand(u32 *seed1, u32 *seed2)
{
  /* A multiply-with-carry generator by George Marsaglia.
     The period is about 2^63. */

  u64 temp = (u64)(*seed1) * (u64)RAND_MULT + (u64)(*seed2);
  *seed1 = (u32)temp;
  *seed2 = (u32)(temp >> 32);
  return (u32)temp;
}

/* masks for picking out individual bits of 64-bit
   words, used for the linear algebra */

#define B(x) ((u64)(1) << (x))

static const u64 bitmask[] = {
  B( 0), B( 1), B( 2), B( 3), B( 4), B( 5), B( 6), B( 7),
  B( 8), B( 9), B(10), B(11), B(12), B(13), B(14), B(15),
  B(16), B(17), B(18), B(19), B(20), B(21), B(22), B(23),
  B(24), B(25), B(26), B(27), B(28), B(29), B(30), B(31),
  B(32), B(33), B(34), B(35), B(36), B(37), B(38), B(39),
  B(40), B(41), B(42), B(43), B(44), B(45), B(46), B(47),
  B(48), B(49), B(50), B(51), B(52), B(53), B(54), B(55),
  B(56), B(57), B(58), B(59), B(60), B(61), B(62), B(63),
};

static const u8 primes_tiny[30] =  /* First 30 primes */
{  2,  3,  5,  7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
  61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113};

/******************************************************************************/

/* Start with a small dependency margin, extending it in equal steps when
   the generated dependencies are trivial. */
#define EXTRA_RELATION_STEP_TINY 32
#define MAX_EXTRA_RELATIONS_TINY 256

/* largest number of relations that can go into the
   linear algebra (includes relations combined from
   pairs of partial relations */
#define MAX_RELATIONS_TINY 960

/* the largest possible factor base */
#define MAX_FB_SIZE_TINY (MAX_RELATIONS_TINY - \
                          MAX_EXTRA_RELATIONS_TINY)

/* offset of the first valid factor base prime */
#define MIN_FB_OFFSET_TINY 1

/* offset of the first factor base prime
   actually contributing to the sieving */
#define MIN_FB_OFFSET_TO_SIEVE_TINY 7

/* number of primes used when testing multipliers */
#define NUM_TEST_PRIMES_TINY 30

/* fudge factor to the target sieve value to account
   for not sieving with the smallest factor base primes */
#define SMALL_PRIME_FUDGE_TINY 10

/* maximum number of SIQS polynomials to be computed */
#define MAX_POLY_TINY 2048

/* maximum number of FB primes that contribute to
   a single polynomial 'A' value */
#define MAX_POLY_FACTORS_TINY 5

/* the size of the sieve interval. Each polynomial will
   sieve over this many positive and negative values */
#define SIEVE_SIZE_TINY 16384

/* Value of the sieve root used when sieving is not performed for a given
   FB prime.  It lies outside both the valid root and sieve-index ranges. */
#define DO_NOT_SIEVE_TINY 65535

/* maximum number of factors a relation can have (the
   large prime is stored separately) */
#define MAX_FACTORS_TINY 20

/* partial relations are listed in the order
   in which they occur, and a hashtable matches
   up partial relations with the same large prime. */
#define LOG2_PARTIAL_TABLE_SIZE 12
#define LARGE_PRIME_HASH(x) (((u32)(x) * (u32)0x9e3779b1U) >> \
                                (32 - LOG2_PARTIAL_TABLE_SIZE))

/* number of slots for each hashtable entry */
#define LP_HASH_DEPTH_TINY 2

/* scale factor for all log values */
#define LOGPRIME_SCALE_TINY 2

/* maximum number of relations to be saved for
   resieving, used in place of trial factoring */
#define SIEVE_BATCH_SIZE_TINY 128

/* maximum size of the pool of FB primes that
   can appear in a polynomial 'A' value */
#define POLY_SELECT_BITS_TINY 12

#define TINY_POSITIVE 0
#define TINY_NEGATIVE 1

/* Restrictions on the above:
 * SIEVE_BATCH_SIZE_TINY <= 128; batch indices use the low 7 sieve-byte bits.
 * SIEVE_SIZE_TINY must be divisible by 32 and no larger than 32768.
 * 4 * MAX_RELATIONS_TINY must fit below the u16 hash sentinel 0xffff.
 * MAX_POLY_TINY must fit u16 poly_num.
 * POLY_SELECT_BITS_TINY <= 16.
 * MAX_FACTORS_TINY + 1 must fit u8 because that value is an overflow sentinel.
 */



/* structure describing a single relation */

typedef struct {
  u32 large_prime;      /* the large prime (may be 1) */
  s16 sieve_offset;     /* the sieve offset of the relation */
  u16 poly_num;         /* ID of the poly that produced the relation */
  u8 num_factors;       /* number of factors from the factor base
                           (duplicates count) */
  u16 fb_offsets[MAX_FACTORS_TINY]; /* offsets into FB of primes that
                                       divide this relation */
} tiny_relation;

/* structure describing a factor base entry */

typedef struct {
  u16 prime;    /* the factor base prime */
  u16 modsqrt;  /* x that solves x^2 = N mod p */
  u32 recip;    /* integer reciprocal of 'prime' */
  u8 logprime;  /* log value used in sieve */
  u16 roots[2]; /* the two sieve roots for 'prime' */
} tiny_fb;

/* structure describing one SIQS polynomial */

typedef struct {
  u16 a_fb_offsets[MAX_POLY_FACTORS_TINY];  /* factors of 'A' value */
  tq_t b;                                   /* B value */
} tiny_poly;

/* main structure controlling the factorization */

typedef struct {

  /* basic stuff */

  uint128_t input_n;                /* original number to be factored */
  tq_t n;                           /* multiplier * input_n */
  u16 multiplier_fb[2];             /* fb offsets of factors of multiplier */
  u32 rand_seed1, rand_seed2;

  /* polynomial selection stuff */

  double target_a;                  /* the optimal size of poly A values */
  s32 poly_num;                     /* ID of current polynomial */
  s32 num_a_factors;                /* # of factors in poly 'A' values */
  s32 poly_select_idx;              /* ID of the combination of primes
                                       that will make current A value */
  u16 poly_select_mask;             /* generated combinations after table */
  u16 poly_select_offsets[POLY_SELECT_BITS_TINY]; /* pool of primes for A */
  tq_t poly_b_aux[MAX_POLY_FACTORS_TINY];      /* scratch values for com-
                                                   puting poly B values */
  tiny_poly poly_list[MAX_POLY_TINY];      /* list of SIQS polynomials */
  tq_t poly_a, poly_b, poly_c;

  /* sieve stuff */

  u64 sieve_block[SIEVE_SIZE_TINY / sizeof(u64)];  /* the sieve interval */

  /* factor base stuff */

  s32 fb_size;                      /* number of FB primes */
  float test_prime_contrib[NUM_TEST_PRIMES_TINY]; /* scratch space used in
                                                     multiplier selection */
  tiny_fb factor_base[MAX_FB_SIZE_TINY];          /* the factor base */
  u16 root_aux[MAX_POLY_FACTORS_TINY *
               MAX_FB_SIZE_TINY];      /* scratch value for initializing
                                          sieve roots */
  /* relation stuff */

  s32 num_full_relations;   /* where next full relation will go */
  s32 partial_idx;          /* where next partial relation will go */
  s32 target_relations;     /* factor base plus dependency margin */
  s32 large_prime_max;      /* max value of a large prime */
  s32 error_bits;           /* value used for trial factoring cutoff */
  tiny_relation sieve_batch[SIEVE_BATCH_SIZE_TINY]; /* resieved relations */

  /* all relations that survive sieving are put in relation_list.
     Full relations (and partial relations whose large prime has
     occurred more than once) are stored in a list that grows up
     from the beginning of the list, while partial relations that
     have not been matched up yet are stored in a list growing down
     from the end of relation_list. num_full_relations is the index
     of the first free space for full relations, and partial_idx
     does the same for unmatched partial relations. */

  tiny_relation relation_list[4 * MAX_RELATIONS_TINY];

  /* a hashtable is used to match up partial relations, using the
     large prime as a hash key. The hashtable stores the index in
     relation_list of the partial relation that connects up all the
     other partial relations with the same large prime (those other
     relations are treated as full relations) */

  u16 partial_hash[1 << LOG2_PARTIAL_TABLE_SIZE][LP_HASH_DEPTH_TINY];

  /* linear algebra stuff */

  u16 null_vectors[MAX_RELATIONS_TINY];
  u64 matrix[MAX_FB_SIZE_TINY][(MAX_RELATIONS_TINY+63) / 64];
} tiny_qs_params;


/* Utility routines.  These are not too performance sensitive,
 * but because the inputs are only 16-bit, these can be fast. */

/* Compute the Legendre symbol (a/p) */
static s32 legendre_16(s32 a, s32 p) {
  u32 x = (u32)a, y = (u32)p, neg = 0;

  while (x) {
    u32 shift = (u32)ctz((UV)x);
    u32 tmp;
    x >>= shift;
    neg ^= (shift & 1) && ((y + 2) & 4);  /* Supplement for (2/y). */
    neg ^= (x & y & 2) != 0;              /* Quadratic reciprocity. */
    tmp = x;
    x = y % x;
    y = tmp;
  }
  return (y == 1) ? 1 - 2 * (s32)neg : 0;
}

/* Compute a^b mod n */
static s32 powm_16(s32 a, s32 b, s32 n) {
  s32 res = 1;
  while (b) {
    if (b & 1)
      res = res * a % n;
    a = a * a % n;
    b = b >> 1;
  }
  return res;
}

/* Modular inverse of a mod p (a > 0, p prime) */
static s32 modinv_16(s32 a, s32 p) {
  u32 dividend = (u32)p, divisor = (u32)a;
  u32 ps1 = 1, ps2 = 0, parity = 0;

  while (divisor > 1) {
    u32 q = dividend / divisor;
    u32 rem = dividend - q * divisor;
    u32 next = q * ps1 + ps2;
    dividend = divisor;
    divisor = rem;
    ps2 = ps1;
    ps1 = next;
    parity ^= 1;
  }
  return parity ? p - (s32)ps1 : (s32)ps1;
}

/* Deterministic Tonelli-Shanks for the small factor-base primes. */
static s32 sqrtModP_16(s32 a, s32 p)
{
  s32 q, s, z, c, x, t, m;

  if (a == 0) return 0;
  if ((p & 3) == 3) {
    return powm_16(a, (p+1)/4, p);
  }
  if ((p & 7) == 5) {
    u32 alpha = (u32)powm_16(2*a % p, (p-5)/8, p);
    u32 root = (u32)a * alpha % (u32)p;
    u32 beta = root * alpha % (u32)p;
    beta += beta;
    if (beta >= (u32)p) beta -= (u32)p;
    return (s32)(root * (beta-1) % (u32)p);
  }
  if ((p & 15) == 9) {  /* Mueller's extension of Atkin. */
    u32 v = (u32)powm_16(2*a % p, (p-9)/16, p);
    u32 root = (u32)a * v % (u32)p;
    u32 beta = root * v % (u32)p;
    beta += beta;
    if (beta >= (u32)p) beta -= (u32)p;

    if (beta * beta % (u32)p != (u32)p-1) {
      u32 d = 1;
      do { d += 2; }
      while (legendre_16((s32)d, p) != -1 && d < (u32)p);
      v = v * (u32)powm_16((s32)d, (p-1)/8, p) % (u32)p;
      root = (u32)a * v % (u32)p;
      beta = root * v % (u32)p;
      beta += beta;
      if (beta >= (u32)p) beta -= (u32)p;
    }
    return (s32)(root * (beta-1) % (u32)p);
  }

  s = ctz((UV)(p-1));
  q = (p-1) >> s;

  for (z = 3; legendre_16(z, p) != -1; z += 2)
    ;
  c = powm_16(z, q, p);
  t = powm_16(a, q >> 1, p);
  x = t * a % p;
  t = t * x % p;
  m = s;

  while (t != 1) {
    s32 i = 1, j, b;
    s32 t2 = t * t % p;
    while (i < m && t2 != 1) {
      t2 = t2 * t2 % p;
      i++;
    }
    MPUassert(i < m, "tinyqs128 Tonelli-Shanks failed");
    b = c;
    for (j = m-i-1; j > 0; j--)
      b = b * b % p;
    x = x * b % p;
    b = b * b % p;
    t = t * b % p;
    c = b;
    m = i;
  }
  return x;
}


static void init_tinyqs(tiny_qs_params *p) {
  s32 i;

  tq_init(p->n);
  tq_init(p->poly_a);
  tq_init(p->poly_b);
  tq_init(p->poly_c);
  /* Fold all input words and avoid the MWC generator's absorbing zero state. */
  p->rand_seed1 = 11111111U ^ (u32)p->input_n ^
                                 (u32)(p->input_n >> 64);
  p->rand_seed2 = 22222222U ^ (u32)(p->input_n >> 32) ^
                                 (u32)(p->input_n >> 96);
  if ((p->rand_seed1 | p->rand_seed2) == 0) {
    p->rand_seed1 = 11111111U;
    p->rand_seed2 = 22222222U;
  }

  /* init the scratch values for polynomial 'B'
     value computations */

  for (i = 0; i < MAX_POLY_FACTORS_TINY; i++) {
    tq_init(p->poly_b_aux[i]);
  }

  /* set up the list of sieve polynomials */

  for (i = 0; i < MAX_POLY_TINY; i++) {
    tq_init(p->poly_list[i].b);
  }

  /* see the next routine for an explanation of what
     these quantities are */

  for (i = 1; i < NUM_TEST_PRIMES_TINY; i++) {
    p->test_prime_contrib[i] = 2 * log((double)primes_tiny[i]) /
                               (primes_tiny[i] - 1) / M_LN2;
  }
}


/* Implementation of the modified Knuth-Schroeppel multiplier
   algorithm. This borrows ideas from at least four different
   sources, and seems to choose multipliers that are better on
   average than many of the other methods available.

   There are many misconceptions about what this algorithm is
   supposed to do. We want to multiply the input number n by a
   small odd squarefree constant k, chosen so that the factor base
   for k * n contains as many small primes as possible. Since small primes
   occur more often than big ones, this makes sieve values smaller
   on average and so more likely to be smooth. We quantify this
   by measuring the average contribution of the first NUM_TEST_PRIMES_TINY
   primes to sieve values. There are two constraints: first, larger
   multipliers mean a larger number to factor. Second, we can't spend
   all day testing multipliers, so the set of multipliers to test should
   be small.

   The list of available multipliers depends on the value of n mod
   8, 3, and 5; each row of the table below gives the multipliers
   to try, pre-sorted by how well they approximately optimize sieving
   (the routine below computes a better approximation). Note that a
   multiplier of 1 (i.e. no multiplier) is always possible. Experiments
   show that 90% of the time the optimal multiplier is in one of the
   first four columns of the table */

#define MAX_MULTIPLIERS 13                           /* for residue classes: */
static const u8 mult_list[32][MAX_MULTIPLIERS] = {   /* mod 8  mod 3  mod 5 */
{ 1, 19, 61, 31, 21, 13,  7,  3, 73, 41,  5, 33, 37 }, /*  1      1      1 */
{ 1, 13,  7,  3, 73, 33, 37, 17, 57, 43,  5, 19, 15 }, /*  1      1      2 */
{ 1, 13,  7,  3, 73, 33, 37, 17, 57, 43,  5, 19, 15 }, /*  1      1      3 */
{ 1, 19, 61, 31, 21, 13,  7,  3, 73, 41,  5, 33, 37 }, /*  1      1      4 */
{ 1, 41,  5, 17, 11, 89, 29, 65, 21,  3, 59, 33, 35 }, /*  1      2      1 */
{ 1, 17,  5,  3, 33, 65, 57, 23, 41, 53, 47, 11, 89 }, /*  1      2      2 */
{ 1, 17,  5,  3, 33, 65, 57, 23, 41, 53, 47, 11, 89 }, /*  1      2      3 */
{ 1, 41,  5, 17, 11, 89, 29, 65, 21,  3, 59, 33, 35 }, /*  1      2      4 */
{ 1, 19,  3, 11, 31,  7, 51, 43, 15, 39, 61, 55, 21 }, /*  3      1      1 */
{ 1,  3,  7, 43, 19, 13, 37, 15, 55, 11, 73, 31, 35 }, /*  3      1      2 */
{ 1,  3,  7, 43, 19, 13, 37, 15, 55, 11, 73, 31, 35 }, /*  3      1      3 */
{ 1, 19,  3, 11, 31,  7, 51, 43, 15, 39, 61, 55, 21 }, /*  3      1      4 */
{ 1, 11,  3, 59, 35,  5, 51, 19, 29, 41, 15, 23, 39 }, /*  3      2      1 */
{ 1,  3, 11, 35,  5, 23, 17, 47,  7, 59, 43, 15, 53 }, /*  3      2      2 */
{ 1,  3, 11, 35,  5, 23, 17, 47,  7, 59, 43, 15, 53 }, /*  3      2      3 */
{ 1, 11,  3, 59, 35,  5, 51, 19, 29, 41, 15, 23, 39 }, /*  3      2      4 */
{ 1, 61, 21, 13,  5, 19, 37, 31, 29,  7,  3, 11, 15 }, /*  5      1      1 */
{ 1, 13, 37,  7,  3,  5, 73, 61, 21, 43, 33, 53, 17 }, /*  5      1      2 */
{ 1, 13, 37,  7,  3,  5, 73, 61, 21, 43, 33, 53, 17 }, /*  5      1      3 */
{ 1, 61, 21, 13,  5, 19, 37, 31, 29,  7,  3, 11, 15 }, /*  5      1      4 */
{ 1,  5, 29, 21, 11, 41, 53, 17, 89,  3, 59, 61, 65 }, /*  5      2      1 */
{ 1,  5, 53, 17,  3, 13, 29, 23, 21, 37, 47, 33, 11 }, /*  5      2      2 */
{ 1,  5, 53, 17,  3, 13, 29, 23, 21, 37, 47, 33, 11 }, /*  5      2      3 */
{ 1,  5, 29, 21, 11, 41, 53, 17, 89,  3, 59, 61, 65 }, /*  5      2      4 */
{ 1, 31,  7, 19, 15, 39, 55,  3, 11, 61, 21, 13, 51 }, /*  7      1      1 */
{ 1,  7,  3, 15, 13, 55, 31, 43, 23, 37, 19, 47, 73 }, /*  7      1      2 */
{ 1,  7,  3, 15, 13, 55, 31, 43, 23, 37, 19, 47, 73 }, /*  7      1      3 */
{ 1, 31,  7, 19, 15, 39, 55,  3, 11, 61, 21, 13, 51 }, /*  7      1      4 */
{ 1, 11,  5, 15, 23, 39,  3, 29, 47, 59, 31, 35,  7 }, /*  7      2      1 */
{ 1, 23,  3, 47,  7,  5, 15, 17, 11, 35, 53, 39, 33 }, /*  7      2      2 */
{ 1, 23,  3, 47,  7,  5, 15, 17, 11, 35, 53, 39, 33 }, /*  7      2      3 */
{ 1, 11,  5, 15, 23, 39,  3, 29, 47, 59, 31, 35,  7 }, /*  7      2      4 */
};

static void find_multiplier_tiny(tiny_qs_params *params) {
  s32 i, j;
  u16 test_nmodp[NUM_TEST_PRIMES_TINY];
  s32 best_mult = 1;
  s32 nmod8 = tq_get_ui(params->n) % 8;
  float best_score;
  const u8 *mult_row;
  s32 num_tests;

  /* Precompute information that will be needed for all multipliers */

  for (i = 1; i < NUM_TEST_PRIMES_TINY; i++)
    test_nmodp[i] = tq_tdiv_ui(params->n, primes_tiny[i]);

  /* Find the row of the table that is appropriate for this value of n */

  mult_row = mult_list[ test_nmodp[2] - 1 +
                        4*(test_nmodp[1] - 1) +
		        8*(nmod8 / 2) ];

  /* test less than the whole row if n is small */

  num_tests = tq_sizeinbase(params->n, 2) / 10;
  if (num_tests > MAX_MULTIPLIERS)
    num_tests = MAX_MULTIPLIERS;

  best_score = 1000.0;
  for (i = 0; i < num_tests; i++) {
    s32 curr_mult = mult_row[i];
    s32 knmod8 = (nmod8 * curr_mult) % 8;
    float score;

    /* The GMP implementation can always form k*n.  Native tinyQS keeps
     * every polynomial value within uint128_t, so reject overflowing k. */
    if (params->n->mag > TQ_U128_MAX / (u32)curr_mult)
      continue;

    /* measure the contribution of 2 as a factor of sieve
       values. The multiplier itself must also be taken into
       account in the score. 'score' is the correction that
       is implicitly applied to the size of sieve values; a
       negative score makes sieve values smaller, and so is
       better. */

    if (knmod8 == 1)
      score = 0.5 * log((double)curr_mult) / M_LN2 - 2;
    else if (knmod8 == 5)
      score = 0.5 * log((double)curr_mult) / M_LN2 - 1;
    else
      score = 0.5 * log((double)curr_mult) / M_LN2 - 0.5;

    for (j = 1; j < NUM_TEST_PRIMES_TINY; j++) {
      s32 prime = primes_tiny[j];
      s32 knmodp = (s32)test_nmodp[j] * curr_mult % prime;

      /* if prime j is actually in the factor base
         for k * n ... */

      if (legendre_16(knmodp, prime) != -1) {

        /* ...add its contribution. A prime p con-
           tributes log(p) to 1 in p sieve values, plus
           log(p) to 1 in p^2 sieve values, etc. The
           average contribution of all multiples of p
           to a random sieve value is thus

           log(p) * (1/p + 1/p^2 + 1/p^3 + ...)
           = (log(p) / p) * 1 / (1 - (1/p))
           = log(p) / (p-1)

           This contribution occurs once for each
           square root used for sieving. There are two roots
           for each factor base prime. If the prime divides
           the squarefree multiplier, its one root modulo p
           does not lift to higher powers, so its contribution
           is instead log(p)/p. The scores are premultiplied by
           2.0, and logarithms are in base 2 (though any base
           will do) */

        if (knmodp == 0)
          score -= 0.5 * params->test_prime_contrib[j]
                     * (prime - 1) / prime;
        else
          score -= params->test_prime_contrib[j];
      }
    }
    if (score < best_score) {
      best_score = score;
      best_mult = curr_mult;
    }
  }
  /* from now on we will factor best_mult * n */
  tq_mul_ui(params->n, params->n, best_mult);
}


static s32 init_fb_tiny(tiny_qs_params *params, s32 fb_size,
                        uint128_t *factor) {
  s32 i, mult_idx;
  tiny_fb *factor_base = params->factor_base;
  s32 prime = 3;

  i = MIN_FB_OFFSET_TINY;
  mult_idx = 0;
  factor_base[i].prime = 2;
  params->multiplier_fb[0] = 0;
  params->multiplier_fb[1] = 0;

  /* Keep setting up factor base primes until enough
     are found or the pool of primes runs out */

  for (i++; i < fb_size && prime < 16384; prime = next_prime(prime)) {
    tiny_fb *fbptr = factor_base + i;
    s32 nmodp = tq_tdiv_ui(params->n, prime);

    if (legendre_16(nmodp, prime) != -1) {
      fbptr->prime = prime;
      fbptr->logprime = (u8)(LOGPRIME_SCALE_TINY *
                             log((double)prime) / M_LN2 + 0.5);
      fbptr->recip = (u32)(B(32) / (u64)prime);

      /* if the prime divides n, it is part of n's
         multiplier and is treated separately */

      if (nmodp != 0) {
        fbptr->modsqrt = (u16)sqrtModP_16(nmodp, prime);
#ifdef DEBUGGING
        MPUassert((u32)fbptr->modsqrt * fbptr->modsqrt % prime ==
                  (u32)nmodp, "tinyqs128 invalid modular square root");
#endif
      }
      else {
        if (params->input_n % (u32)prime == 0) {
          *factor = prime;
          return 0;
        }
        fbptr->modsqrt = DO_NOT_SIEVE_TINY;
        MPUassert(mult_idx < 2, "tinyqs128 multiplier has too many factors");
        params->multiplier_fb[mult_idx++] = i;
      }
      i++;
    }
  }
  params->fb_size = i;
  return i;
}


static INLINE void sieve_one_root_tiny(u8 *sieve_block, u32 root,
                                       u32 prime, u8 logprime)
{
  u32 prime4 = prime << 2;

  while (root + prime4 < SIEVE_SIZE_TINY) {
    sieve_block[root]             -= logprime;
    sieve_block[root + prime]     -= logprime;
    sieve_block[root + 2 * prime] -= logprime;
    sieve_block[root + 3 * prime] -= logprime;
    root += prime4;
  }
  while (root < SIEVE_SIZE_TINY) {
    sieve_block[root] -= logprime;
    root += prime;
  }
}

/* Alternating-root-gap sieve structure inspired by PARI/GP's MPQS sieve. */
static INLINE void sieve_two_roots_tiny(u8 *sieve_block, u32 root1,
                                        u32 root2, u32 prime, u8 logprime)
{
  u32 root, gap1, gap2;
  u32 prime4 = prime << 2;

  if (root1 > root2) {
    u32 tmp = root1;
    root1 = root2;
    root2 = tmp;
  }

  /* Visit the sorted roots using their alternating gaps.  This keeps one
     induction variable and one loop bound for both progressions. */
  root = root1;
  gap1 = root2 - root1;
  gap2 = prime - gap1;
  while (root + prime4 < SIEVE_SIZE_TINY) {
    sieve_block[root] -= logprime; root += gap1;
    sieve_block[root] -= logprime; root += gap2;
    sieve_block[root] -= logprime; root += gap1;
    sieve_block[root] -= logprime; root += gap2;
    sieve_block[root] -= logprime; root += gap1;
    sieve_block[root] -= logprime; root += gap2;
    sieve_block[root] -= logprime; root += gap1;
    sieve_block[root] -= logprime; root += gap2;
  }
  while (root < SIEVE_SIZE_TINY) {
    sieve_block[root] -= logprime;
    root += gap1;
    if (root >= SIEVE_SIZE_TINY)
      break;
    sieve_block[root] -= logprime;
    root += gap2;
  }
}

/* Core sieving routine */
static void fill_sieve_block_tiny(tiny_qs_params *params) {
  s32 i, fb_size = params->fb_size;
  u8 *sieve_block = (u8*) params->sieve_block;
  tiny_fb *factor_base = params->factor_base;

  /* The interval is small enough that it does not need to be segmented.
     Interleave the two root progressions, however: this removes one of the
     induction variables and exposes a compact sequence of sieve writes. */

  for (i = MIN_FB_OFFSET_TO_SIEVE_TINY; i < fb_size; i++) {
    tiny_fb *fbptr = factor_base + i;
    u32 prime = fbptr->prime;
    u8 logprime = fbptr->logprime;
    u32 root1 = fbptr->roots[0];
    u32 root2 = fbptr->roots[1];

    if (root1 == DO_NOT_SIEVE_TINY) {
      if (root2 != DO_NOT_SIEVE_TINY)
        sieve_one_root_tiny(sieve_block, root2, prime, logprime);
    }
    else if (root2 == DO_NOT_SIEVE_TINY) {
      sieve_one_root_tiny(sieve_block, root1, prime, logprime);
    }
    else {
      sieve_two_roots_tiny(sieve_block, root1, root2, prime, logprime);
    }
  }
}


#define PACKED_SIEVE_MASK ((u64)0x80808080 << 32 | 0x80808080)

/* Walk a sieve block, find offsets corresponding to useful relations. */
static s32 mark_sieve_block_tiny(tiny_qs_params *params) {
  s32 i, j, k;
  u8  *sieve_block   = (u8*) params->sieve_block;
  u64 *sieve_block64 = params->sieve_block;

  /* standard technique for testing sieve locations
     in parallel: initialize each byte to the target
     sieve value, and subtract logs of the factor base
     primes instead of adding them. Sieve offsets that
     accumulate enough log values become negative,
     and it's easy to simultaneously test for the top
     bit in several bytes being set */

  for (i = j = 0; i < SIEVE_SIZE_TINY / 8; i += 4) {

    /* handle 32 bytes at a time */
    u64 accum = sieve_block64[i] |
                sieve_block64[i+1] |
                sieve_block64[i+2] |
                sieve_block64[i+3];

    if ((accum & PACKED_SIEVE_MASK) == (u64)(0))
      continue;

    /* at least one byte is a hit; go back and search
       the list one at a time. We treat the sieve interval
       as a hashtable, and associate entry j in the list
       of relations to be resieved (params->sieve_batch[])
       with a byte that is negative. The high-order bit of
       the byte is set to indicate that the low-order bits
       mean something */

    for (k = 0; k < 32; k++) {
      s32 offset = 8 * i + k;
      if (sieve_block[offset] & 0x80) {
        tiny_relation *r = params->sieve_batch + j;
        r->sieve_offset = offset;
        r->num_factors = 0;
        sieve_block[offset] = j | 0x80;

        if (++j == SIEVE_BATCH_SIZE_TINY)
          return j;
      }
    }
  }

  return j;
}


/* Like fill_sieve_block_tiny, but sieving to avoid trial division */
static void resieve_tiny(tiny_qs_params *params, s32 num_surviving) {
  s32 i, fb_size = params->fb_size;
  s32 sieve_end = params->sieve_batch[num_surviving - 1].sieve_offset + 1;
  u8 *sieve_block = (u8*) params->sieve_block;
  tiny_fb *factor_base = params->factor_base;

  /* Even though this routine does only a little more work
     than fill_sieve_block_tiny(), it runs almost 3x slower */

  for (i = MIN_FB_OFFSET_TO_SIEVE_TINY; i < fb_size; i++) {
    tiny_fb *fbptr = factor_base + i;
    s32 prime = fbptr->prime;
    s32 root = fbptr->roots[0];

    while (root < sieve_end) {
      s32 val = sieve_block[root];
      if (val & 0x80) {
        tiny_relation *r = params->sieve_batch + (val & 0x7f);
        if (r->num_factors < MAX_FACTORS_TINY)
          r->fb_offsets[r->num_factors++] = i;
        else
          r->num_factors = MAX_FACTORS_TINY + 1;
      }
      root += prime;
    }

    root = fbptr->roots[1];
    while (root < sieve_end) {
      s32 val = sieve_block[root];
      if (val & 0x80) {
        tiny_relation *r = params->sieve_batch + (val & 0x7f);
        if (r->num_factors < MAX_FACTORS_TINY)
          r->fb_offsets[r->num_factors++] = i;
        else
          r->num_factors = MAX_FACTORS_TINY + 1;
      }
      root += prime;
    }
  }
}


/* Trial factor a relation that survived sieving */
static s32 check_sieve_val_tiny(tiny_qs_params *params,
                                 tq_t a, tq_t b, tq_t c,
                                 tiny_relation *r,
                                 s32 sign_of_index)
{
  s32 i, j;
  tiny_fb *factor_base = params->factor_base;
  s32 num_factors = 0;
  s32 sieve_offset = r->sieve_offset;
  tiny_relation *relation = params->relation_list +
                            params->num_full_relations;
  u16 *fb_offsets = relation->fb_offsets;
  tq_t res, res2;

  tq_init(res);
  tq_init(res2);

  /* resieve_tiny records a sentinel rather than overrun fb_offsets. */
  if (r->num_factors > MAX_FACTORS_TINY)
    return 0;

  /* form the polynomial value */

  tq_mul_ui(res, a, sieve_offset);
  if (sign_of_index == TINY_POSITIVE)
    tq_add(res, res, b);
  else
    tq_sub(res, res, b);
  tq_mul_ui(res, res, sieve_offset);
  tq_add(res, res, c);
  if (tq_sgn(res) < 0) {
    tq_abs(res, res);
    fb_offsets[num_factors++] = 0;
  }

  /* extract powers of two */

  i = tq_scan1(res, 0);
  if (i) {
    tq_tdiv_q_2exp(res, res, i);
    do {
      if (num_factors >= MAX_FACTORS_TINY)
        return 0;
      fb_offsets[num_factors++] = MIN_FB_OFFSET_TINY;
    } while (--i);
  }

  /* divide out the unsieved factor base primes */

  for (i = MIN_FB_OFFSET_TINY + 1;
             i < MIN_FB_OFFSET_TO_SIEVE_TINY; i++) {
    tiny_fb *fbptr = factor_base + i;
    s32 prime = fbptr->prime;
    s32 root1 = fbptr->roots[0];
    s32 root2 = fbptr->roots[1];
    u32 recip = fbptr->recip;

    if (root1 == DO_NOT_SIEVE_TINY)
      continue;

    j = (s32)(((u64)sieve_offset * (u64)recip) >> 32);
    j = sieve_offset - j * prime;
    if (j >= prime)
      j -= prime;

    if (j == root1 || j == root2) {
      while (tq_tdiv_q_ui(res2, res, prime) == 0) {
        if (num_factors >= MAX_FACTORS_TINY)
          return 0;

        fb_offsets[num_factors++] = i;
        tq_swap(res, res2);
      }
    }
  }

  /* divide out the factors of the multiplier, if any */

  for (i = 0; i < 2; i++) {
    if (params->multiplier_fb[i]) {
      s32 prime;
      j = params->multiplier_fb[i];
      prime = factor_base[j].prime;
      while (tq_tdiv_q_ui(res2, res, prime) == 0) {
        if (num_factors >= MAX_FACTORS_TINY)
          return 0;

        fb_offsets[num_factors++] = j;
        tq_swap(res, res2);
      }
    }
  }

  /* We should probably have been adding log values
     to the log of this relation in the previous loops,
     and testing that the complete log value now
     exceeds the trial factoring cutoff. However,
     resieving has already found the remaining factors,
     so we wouldn't save much time bailing out at
     this point */

  for (i = 0; i < r->num_factors; i++) {
    s32 prime;
    j = r->fb_offsets[i];
    prime = factor_base[j].prime;

    while (tq_tdiv_q_ui(res2, res, prime) == 0) {
      if (num_factors >= MAX_FACTORS_TINY)
        return 0;

      fb_offsets[num_factors++] = j;
      tq_swap(res, res2);
    }
  }

  /* start filling in the final relation */

  if (sign_of_index == TINY_NEGATIVE)
    sieve_offset = -sieve_offset;
  relation->sieve_offset = sieve_offset;
  relation->num_factors = num_factors;
  relation->poly_num = params->poly_num;

  if (tq_cmp_ui(res, 1) == 0) {

    /* full relation; we're done */

    relation->large_prime = 1;
    params->num_full_relations++;
  }
  else if (tq_cmp_ui(res, params->large_prime_max) < 0) {
    u32 lp = tq_get_ui(res);
    u32 table_idx = LARGE_PRIME_HASH(lp);
    s32 partial_idx;

    /* partial relation; see if it has occurred already */

    relation->large_prime = lp;
    for (i = 0; i < LP_HASH_DEPTH_TINY; i++) {
      partial_idx = params->partial_hash[table_idx][i];
      if (partial_idx == 0xffff ||
          lp == params->relation_list[partial_idx].large_prime)
        break;
    }

    if (i == LP_HASH_DEPTH_TINY)
      return 0; /* not found, and no room to store it */

    if (partial_idx == 0xffff) {

      /* not found, but the hashtable entry has
         room to keep it; transfer the relation to
         the partial list */
      params->relation_list[params->partial_idx] = *relation;
      params->partial_hash[table_idx][i] = params->partial_idx--;
    }
    else {

      /* large prime has matched, new relation can stay */
      params->num_full_relations++;
    }
  }

  /* make sure the 'heap' of full relations has not
     overflowed into the 'stack' of partial relations */
  if (params->num_full_relations >= params->partial_idx)
    return -1;
  return 0;
}


/* Initialize the subsystem for forming SIQS sieve polynomials */
static void init_siqs_tiny(tiny_qs_params *params) {
  u32 i, j, plus_idx, minus_idx, fb_size = params->fb_size;
  u32 num_factors = params->num_a_factors;
  tiny_fb *factor_base = params->factor_base;

  /* compute the optimal size of the factors of
     the polynomial 'A' value. We know how many
     primes it should have, and know the optimal
     A value that will minimize sieving time.
     Assume further that all factors are the
     same size. First compute the factor size,
     then locate the factor base offset where
     it approximately occurs */

  j = (u32)(exp(log(params->target_a) / num_factors) + 0.5);
  for (i = MIN_FB_OFFSET_TINY + 1; i < fb_size - 1; i++) {
    if (factor_base[i].prime > j)
      break;
  }
  if (i == MIN_FB_OFFSET_TINY + 1)
    i++;

  /* polynomial A values are built by selecting from
     a pool of primes. There are POLY_SELECT_BITS_TINY
     primes in the pool, evenly distributed above and
     below the optimal factor base offset */

  memset(params->poly_select_offsets, 0,
         sizeof(params->poly_select_offsets));
  plus_idx = i;
  minus_idx = i-1;
  i = 0;
  while (plus_idx < fb_size || minus_idx > MIN_FB_OFFSET_TINY + 1) {
    if (plus_idx < fb_size &&
        factor_base[plus_idx].modsqrt != DO_NOT_SIEVE_TINY) {
      params->poly_select_offsets[i] = plus_idx;
      if (++i == POLY_SELECT_BITS_TINY)
        break;
    }
    if (minus_idx > MIN_FB_OFFSET_TINY + 1 &&
        factor_base[minus_idx].modsqrt != DO_NOT_SIEVE_TINY) {
      params->poly_select_offsets[i] = minus_idx;
      if (++i == POLY_SELECT_BITS_TINY)
        break;
    }
    plus_idx++;
    if (minus_idx > MIN_FB_OFFSET_TINY + 1)
      minus_idx--;
  }

  /* polynomial selection will begin at offset
     zero of the tables below */

  params->poly_select_idx = 0;
}


/* Choosing the primes for polynomial A is a compromise. A should be
   close to its optimal value, successive choices should be different
   enough to avoid duplicate relations, and selection must be cheap.

   init_siqs_tiny() builds an ordered pool of 12 factor-base primes.
   Bit i in a mask selects pool entry i. Pool entries alternate above
   and below the target factor size, so lower bits generally represent
   primes closer to the target. A valid mask has num_a_factors bits set.

   a_choice[] lists the preferred masks in increasing numeric order.
   Each A value produces 2^(num_a_factors-1) B polynomials. The table
   provides 64, 240, and 256 polynomials for A values containing 3, 4,
   and 5 factors.

   If the preferred masks are exhausted, find_poly_a() enumerates the
   remaining masks of the required weight from the same 12-prime pool. */

static const u16 a_choice[] = {
       0x007, 0x00f, 0x019, 0x01f, 0x02a, 0x033,
       0x034, 0x03c, 0x04c, 0x052, 0x055, 0x05a,
       0x061, 0x066, 0x067, 0x069, 0x079, 0x096,
       0x099, 0x0a5, 0x0aa, 0x0ab, 0x0b5, 0x0c3,
       0x0cc, 0x0cd, 0x0d3, 0x0f0, 0x12d, 0x133,
       0x14b, 0x155, 0x181, 0x187, 0x199, 0x1e1,
       0x22e, 0x256, 0x282, 0x303, 0x304, 0x30c,
       0x330, 0x3c0, 0x484, 0x502, 0x505, 0x50a,
       0x550, 0x5a0, 0x601, 0x606, 0x609, 0x660,
       0x690, 0x888, 0x906, 0x909, 0x910, 0x960,
       0x990, 0xa05, 0xa0a, 0xa20, 0xa50, 0xc40,
};

/* Compute the next polynomial A value */
static s32 find_poly_a(tiny_qs_params *params, tq_t a) {
  u32 i, j, mask, num_a_factors = params->num_a_factors;
  const u32 num_a_choices = sizeof(a_choice) / sizeof(a_choice[0]);
  tiny_fb *factor_base = params->factor_base;
  tiny_poly *poly = params->poly_list + params->poly_num;

  /* choose the next bitfield representing primes to use */

  for (i = params->poly_select_idx; i < num_a_choices; i++) {
    if ((u32)popcnt(a_choice[i]) == num_a_factors)
      break;
  }
  if (i < num_a_choices) {
    mask = a_choice[i];
    params->poly_select_idx = i + 1;
  } else {
    /* Native tinyQS may need more than the hand-tuned table supplies near
     * 128 bits, so continue with unused combinations from the same pool. */
    u32 ti;
    params->poly_select_idx = num_a_choices;
    for (mask = params->poly_select_mask ? params->poly_select_mask : 1;
         mask < (1U << POLY_SELECT_BITS_TINY); mask++) {
      if ((u32)popcnt(mask) != num_a_factors)
        continue;
      for (ti = 0; ti < num_a_choices; ti++)
        if (mask == a_choice[ti]) break;
      if (ti == num_a_choices)
        break;
    }
    if (mask >= (1U << POLY_SELECT_BITS_TINY))
      return -1;
    params->poly_select_mask = mask + 1;
  }

  /* gather the chosen primes */

  for (i = j = 0; i < POLY_SELECT_BITS_TINY; i++) {
    if (!(mask & (1 << i)))
      continue;

    if (params->poly_select_offsets[i] == 0)
      return -2;

    poly->a_fb_offsets[j] = params->poly_select_offsets[i];
    if (++j == num_a_factors)
      break;
  }

  /* multiply them together */

  tq_set_ui(a, 1);
  for (i = 0; i < num_a_factors; i++) {
    j = poly->a_fb_offsets[i];
    tq_mul_ui(a, a, factor_base[j].prime);
  }

  return 0;
}


/* Compute the first of a list of polynomial B values */
static void find_first_poly_b(tiny_qs_params *params,
                              tq_t a, tq_t b, tq_t c) {
  u32 i, j, fb_size = params->fb_size;
  u32 num_a_factors = params->num_a_factors;
  tiny_fb *factor_base = params->factor_base;
  tiny_poly *poly = params->poly_list + params->poly_num;

  tq_set_ui(b, 0);

  /* fill in the auxiliary quantities needed to
     compute future B values */

  for (i = 0; i < num_a_factors; i++) {
    tiny_fb *fbptr = factor_base + poly->a_fb_offsets[i];
    s32 g, prime = fbptr->prime;

    tq_divexact_ui(params->poly_b_aux[i], a, prime);
    g = tq_tdiv_ui(params->poly_b_aux[i], prime);
    g = modinv_16(g, prime);
    g = (s32)g * fbptr->modsqrt % prime;
    if (g > prime/2)
      g = prime - g;
    tq_mul_ui(params->poly_b_aux[i],
               params->poly_b_aux[i], g);
    tq_add(b, b, params->poly_b_aux[i]);
    tq_add(params->poly_b_aux[i],
            params->poly_b_aux[i],
            params->poly_b_aux[i]);
  }
  /* This first B is the sum of the auxiliary
     quantities computed previously */

  tq_set(poly->b, b);

  /* Form C, a helper for computing the value
     of a polynomial before trial factoring */

#ifdef DEBUGGING
  for (i = 0; i < num_a_factors; i++) {
    u32 p = factor_base[poly->a_fb_offsets[i]].prime;
    u32 bm = tq_tdiv_ui(b, p);
    u32 nm = tq_tdiv_ui(params->n, p);
    MPUassert((u64)bm * bm % p == nm,
              "tinyqs128 first polynomial B is not a square root");
  }
#endif

  tq_mul(c, b, b);
  tq_sub(c, c, params->n);
  tq_divexact(c, c, a);

  /* initialize the factor base for sieving */

  for (i = MIN_FB_OFFSET_TINY + 1; i < fb_size; i++) {
    tiny_fb *fbptr = factor_base + i;
    s32 prime = fbptr->prime;
    s32 modsqrt = fbptr->modsqrt;
    s32 amodp = tq_tdiv_ui(a, prime);
    s32 bmodp = prime - tq_tdiv_ui(b, prime);

    if (fbptr->modsqrt == DO_NOT_SIEVE_TINY) {
      /* factors of the multiplier never contribute to sieving */
      fbptr->roots[0] = DO_NOT_SIEVE_TINY;
      fbptr->roots[1] = DO_NOT_SIEVE_TINY;
      continue;
    }

    if (amodp == 0) {
      /* factor base primes that divide the A value get one sieve root
         and not two. */
      amodp = prime - tq_tdiv_ui(c, prime); /* C<0  =>  C mod prime */
      fbptr->roots[0] = amodp * modinv_16(2 * bmodp % prime, prime) % prime;
      fbptr->roots[1] = DO_NOT_SIEVE_TINY;
    } else {
      /* handle all the other FB primes, including allowing the next
         2^(num_a_factors-1)-1 polynomials to initialize quickly */
      amodp = modinv_16(amodp, prime);
      fbptr->roots[0] = amodp * (bmodp + modsqrt) % prime;
      fbptr->roots[1] = amodp * (bmodp + prime - modsqrt) % prime;
      for (j = 0; j < num_a_factors; j++) {
        bmodp = tq_tdiv_ui(params->poly_b_aux[j], prime);
        params->root_aux[j * fb_size + i] = bmodp * amodp % prime;
      }
    }
  }
}


/* Initialize B values beyond the first */
static void find_next_poly_b(tiny_qs_params *params,
                             tq_t a, tq_t b, tq_t c) {
  s32 i, j, fb_size = params->fb_size;
  s32 num_a_factors = params->num_a_factors;
  tiny_fb *factor_base = params->factor_base;
  tiny_poly *poly = params->poly_list + params->poly_num;
  u32 mask = params->poly_num & ((1 << (num_a_factors-1)) - 1);
  u8 do_sub;
  u16 *row;

  /* current poly starts with the previous poly */

  tq_set(b, poly[-1].b);
  for (i = 0; i < num_a_factors; i++)
    poly[0].a_fb_offsets[i] = poly[-1].a_fb_offsets[i];

  /* determine the auxiliary B constant that comes next
     in Gray code order, and add to or subtract from
     the current B. This also determines which of the
     rows from the table of corrections are applied to
     the factor base */

  i = 0;
  while ((mask & (1 << i)) == 0)
    i++;

  row = params->root_aux + fb_size * i;

  do_sub = 0;
  if (mask & (1 << (i+1))) {
    tq_add(b, b, params->poly_b_aux[i]);
    do_sub = 1;
  } else {
    tq_sub(b, b, params->poly_b_aux[i]);
  }

  /* form the C helper value */

  tq_mul(c, b, b);
  tq_sub(c, c, params->n);
  tq_divexact(c, c, a);

  /* set up the factor base for the next B */

  for (j = MIN_FB_OFFSET_TINY + 1; j < fb_size; j++) {
    tiny_fb *fbptr = factor_base + j;
    s32 prime = fbptr->prime;
    s32 root1 = fbptr->roots[0];
    s32 root2 = fbptr->roots[1];

    /* apply the correction to each sieve root */

    if (root2 != DO_NOT_SIEVE_TINY) {

      /* ordinary FB prime. Note that the previous sieving
         operation negated the roots to use, so they have
         to be negated again before the correction is applied */

      if (root1)
        root1 = prime - root1;
      if (root2)
        root2 = prime - root2;

      if (do_sub) {
        root1 -= row[j];
        root2 -= row[j];
      }
      else {
        root1 += row[j] - prime;
        root2 += row[j] - prime;
      }

      if (root1 < 0)
        root1 += prime;
      if (root2 < 0)
        root2 += prime;
      fbptr->roots[0] = root1;
      fbptr->roots[1] = root2;

    } else if (root1 != DO_NOT_SIEVE_TINY) {

      /* sieving with root1 but not root 2 only
         happens if the prime divides 'A'. Compute
         the new sieve root manually */

      s32 cmodp = prime - tq_tdiv_ui(c, prime);
      s32 bmodp = tq_tdiv_ui(b, prime);
      if (tq_sgn(b) > 0)
        bmodp = prime - bmodp;
      fbptr->roots[0] = cmodp * modinv_16(2 * bmodp % prime, prime) % prime;
    }
  }
  tq_set(poly[0].b, b);
}


/* Do all the sieving for one polynomial */
static s32 sieve_next_poly_tiny(tiny_qs_params *params) {
  s32 i, fb_size = params->fb_size;
  u8 *sieve_block = (u8*) params->sieve_block;
  tiny_fb *factor_base = params->factor_base;
  s32 cutoff1, num_surviving;
  s32 poly_num = params->poly_num;
  s32 target_relations = params->target_relations;
  tq_cell *a = params->poly_a;
  tq_cell *b = params->poly_b;
  tq_cell *c = params->poly_c;

  /* generate the polynomial */

  if (!(poly_num & ((1 << (params->num_a_factors-1))-1))) {
    i = find_poly_a(params, a);
    if (i)
      return i;
    find_first_poly_b(params, a, b, c);
  } else {
    find_next_poly_b(params, a, b, c);
  }

  /* compute the cutoff beyond which trial factoring
     will be used on sieve values. */

  cutoff1 = LOGPRIME_SCALE_TINY * ((s32)tq_sizeinbase(c, 2) -
                  params->error_bits - SMALL_PRIME_FUDGE_TINY - 1);

  /* the trial factoring code wants 2*B and not B */

  tq_add(b, b, b);

  /* sieve over positive offsets, mark the most
     promising offsets, resieve to trial factor
     them all at once and then finish each in turn */

  memset(sieve_block, cutoff1 - 1, SIEVE_SIZE_TINY);
  fill_sieve_block_tiny(params);
  num_surviving = mark_sieve_block_tiny(params);
  if (num_surviving) {
    resieve_tiny(params, num_surviving);
    for (i = 0; i < num_surviving; i++) {
      if (check_sieve_val_tiny(params, a, b, c,
                               params->sieve_batch + i,
                               TINY_POSITIVE) != 0) {
        return -3;
      }
      if (params->num_full_relations >= target_relations)
        break;
    }
  }

  /* flip the sieve roots from positive to negative values */
  for (i = MIN_FB_OFFSET_TINY + 1; i < fb_size; i++) {
    tiny_fb *fbptr = factor_base + i;
    s32 prime = fbptr->prime;
    s32 root1 = fbptr->roots[0];
    s32 root2 = fbptr->roots[1];
    if (root1 != DO_NOT_SIEVE_TINY && root1)
      fbptr->roots[0] = prime - root1;
    if (root2 != DO_NOT_SIEVE_TINY && root2)
      fbptr->roots[1] = prime - root2;
  }

  if (params->num_full_relations >= target_relations)
    return 0;

  /* repeat the sieve procedure for negative sieve offsets */

  memset(sieve_block, cutoff1 - 1, SIEVE_SIZE_TINY);
  fill_sieve_block_tiny(params);
  sieve_block[0] = 0;  /* x=0 was already handled by the positive sieve */
  num_surviving = mark_sieve_block_tiny(params);
  if (num_surviving) {
    resieve_tiny(params, num_surviving);
    for (i = 0; i < num_surviving; i++) {
      if (check_sieve_val_tiny(params, a, b, c,
                               params->sieve_batch + i,
                               TINY_NEGATIVE) != 0) {
        return -3;
      }
      if (params->num_full_relations >= target_relations)
        break;
    }
  }

  return 0;
}


static s32 find_partial_relation(const tiny_qs_params *params,
                                 u32 large_prime) {
  s32 i;
  s32 hash_idx = LARGE_PRIME_HASH(large_prime);

  for (i = 0; i < LP_HASH_DEPTH_TINY; i++) {
    u16 idx = params->partial_hash[hash_idx][i];
    if (idx != 0xffff &&
        params->relation_list[idx].large_prime == large_prime)
      return idx;
  }
  return -1;
}


/* Find linear dependencies among a set of relations */
static s32 solve_linear_system_tiny(tiny_qs_params *params) {
  s32 i, j, k, start_row, row;
  s32 nrows = params->fb_size;
  s32 ncols = params->num_full_relations;
  s32 num_a_factors = params->num_a_factors;
  u16 rowperm[MAX_FB_SIZE_TINY];
  u16 pivot[MAX_FB_SIZE_TINY];

  memset(params->matrix, 0, sizeof(params->matrix));

  /* build the matrix; relations become columns, and
     pairs of matched partial relations fuse into
     columns as well */

  for (i = 0; i < ncols; i++) {
    tiny_relation *r;
    tiny_poly *poly;
    for (j = 0; j < 2; j++) {

      r = params->relation_list + i;
      if (j == 1) {
        s32 partial_idx = find_partial_relation(params, r->large_prime);
        if (partial_idx < 0)
          return 0;
        r = params->relation_list + partial_idx;
      }
      poly = params->poly_list + r->poly_num;

      for (k = 0; k < r->num_factors; k++) {
        row = r->fb_offsets[k];
        params->matrix[row][i / 64] ^= bitmask[i % 64];
      }

      /* the factors in the polynomial A value
         figure into the matrix as well */

      for (k = 0; k < num_a_factors; k++) {
        row = poly->a_fb_offsets[k];
        params->matrix[row][i / 64] ^= bitmask[i % 64];
      }
      if (r->large_prime == 1)
        break;
    }
  }
  for (i = 0; i < nrows; i++)
    rowperm[i] = i;

  /* begin with a random vector of dependencies */

  for (i = 0; i < ncols; i++)
    params->null_vectors[i] = (u16)get_rand(
                      &params->rand_seed1, &params->rand_seed2);

  /* perform the elimination */

  for (i = start_row = 0; start_row < nrows && i < ncols; i++) {

    /* find the next pivot */

    for (j = start_row; j < nrows; j++) {
      row = rowperm[j];
      if (params->matrix[row][i / 64] & bitmask[i % 64])
        break;
    }
    if (j == nrows)
      continue;

    rowperm[j] = rowperm[start_row];
    rowperm[start_row] = row;
    pivot[start_row++] = i;

    /* eliminate it from the other rows */

    for (j++; j < nrows; j++) {
      s32 row2 = rowperm[j];
      if (params->matrix[row2][i / 64] & bitmask[i % 64]) {
        for (k = i / 64; k < (ncols + 63) / 64; k++) {
          params->matrix[row2][k] ^= params->matrix[row][k];
        }
      }
    }
  }

  /* perform back substitution */

  for (i = start_row - 1; i >= 0; i--) {
    u16 accum;
    row = rowperm[i];

    for (j = pivot[i] + 1, accum = 0; j < ncols; j++) {
      if (params->matrix[row][j / 64] & bitmask[j & 63])
        accum ^= params->null_vectors[j];
    }
    params->null_vectors[pivot[i]] = accum;
  }
  return 1;
}


/* Perform MPQS square root phase */
static u32 find_factors_tiny(tiny_qs_params *params,
                             tq_t factor1, tq_t factor2) {
  s32 i, j, k;
  u16 mask, fb_counts[MAX_FB_SIZE_TINY];
  tiny_fb *factor_base = params->factor_base;
  tq_t x, y, t0, t1;

  tq_init(x); tq_init(y);
  tq_init(t0); tq_init(t1);

  /* for each dependency */

  for (mask = 1; mask; mask <<= 1) {

    memset(fb_counts, 0, sizeof(fb_counts));
    tq_set_ui(x, 1);
    tq_set_ui(y, 1);

    /* for each relation allowed in the dependency */

    for (i = 0; i < params->num_full_relations; i++) {

      if (!(params->null_vectors[i] & mask))
        continue;

      for (j = 0; j < 2; j++) {
        tiny_relation *r = params->relation_list + i;
        tiny_poly *poly;

        /* match up partials with the same large prime */

        if (j == 1) {
          s32 partial_idx = find_partial_relation(params, r->large_prime);
          if (partial_idx < 0)
            return 0;
          r = params->relation_list + partial_idx;
          tq_mulmod_ui(y, y, r->large_prime, params->n);
        }
        poly = params->poly_list + r->poly_num;

        /* add the factors of this relation to the table
           of factors. Include the factors of A as well */

        for (k = 0; k < r->num_factors; k++)
          fb_counts[r->fb_offsets[k]]++;

        tq_set_ui(t1, 1);
        for (k = 0; k < params->num_a_factors; k++) {
          s32 idx = poly->a_fb_offsets[k];
          fb_counts[idx]++;
          tq_mul_ui(t1, t1, factor_base[idx].prime);
        }

        /* multiply A * sieve_offset + B into the left
           side of the congruence */

        if (r->sieve_offset < 0) {
          tq_mul_ui(t1, t1, -(r->sieve_offset));
          tq_sub(t1, t1, poly->b);
        }
        else {
          tq_mul_ui(t1, t1, r->sieve_offset);
          tq_add(t1, t1, poly->b);
        }
        tq_mulmod(x, x, t1, params->n);

        if (r->large_prime == 1)
          break;
      }
    }

    /* Form the right side of the congruence; given its
       prime factorization, cut the exponent of each prime
       in half and perform a modular exponentiation */

    for (i = MIN_FB_OFFSET_TINY; i < params->fb_size; i++) {
      u16 mask2 = 0x8000;
      u16 exponent = fb_counts[i] / 2;
      u32 prime = params->factor_base[i].prime;

      if (exponent == 0)
        continue;

      tq_set_ui(t0, prime);
      while (!(exponent & mask2))
        mask2 >>= 1;

      for (mask2 >>= 1; mask2; mask2 >>= 1) {
        tq_mulmod(t0, t0, t0, params->n);
        if (exponent & mask2) {
          tq_mulmod_ui(t0, t0, prime, params->n);
        }
      }
      tq_mulmod(y, t0, y, params->n);
    }

    /* For x and y the halves of the congruence, compute gcd(x+-y, n) */

    for (i = 0; i < 2; i++) {
      if (i == 0)
        tq_set_ui(t0, tq_addmod_value(x->mag, y->mag, params->n->mag));
      else
        tq_sub(t0, x, y);

      tq_gcd(t1, t0, params->n);
      if (tq_cmp_ui(t1, 1) && tq_cmp(t1, params->n)) {

        /* We've possibly found a nontrivial factor of n.
           Divide any factors of the multiplier out from both factors. */

        u32 mult1 = 0;
        u32 mult2 = 0;

        if (params->multiplier_fb[0])
          mult1 = params->factor_base[params->multiplier_fb[0]].prime;
        if (params->multiplier_fb[1])
          mult2 = params->factor_base[params->multiplier_fb[1]].prime;

        tq_divexact(t0, params->n, t1);
        if (mult1) {
          if (tq_tdiv_ui(t0, mult1) == 0)
            tq_divexact_ui(t0, t0, mult1);
          if (tq_tdiv_ui(t1, mult1) == 0)
            tq_divexact_ui(t1, t1, mult1);
        }
        if (mult2) {
          if (tq_tdiv_ui(t0, mult2) == 0)
            tq_divexact_ui(t0, t0, mult2);
          if (tq_tdiv_ui(t1, mult2) == 0)
            tq_divexact_ui(t1, t1, mult2);
        }

        /* If both remaining factors exceed unity,
           we've factored n and can stop */
        if (tq_cmp_ui(t0, 1) && tq_cmp_ui(t1, 1)) {
          tq_set(factor1, t0);
          tq_set(factor2, t1);
          return 1;
        }
      }
    }

    /* otherwise try the next dependency */
  }
  return 0;
}



typedef struct {
  s32 fb_size;
  s32 num_poly_factors;
} tiny_qs_config;

/* Factor base sizes in four-bit buckets through 128-bit factorizations. */

static const tiny_qs_config static_config[] = {
 { 40, 3 },
 { 50, 3 },
 { 60, 3 },
 { 70, 3 },
 { 80, 3 },
 { 90, 3 },
 { 110, 3 },
 { 120, 3 },
 { 140, 3 },
 { 140, 3 },
 { 160, 3 },
 { 180, 4 },
 { 230, 4 },
 { 280, 4 },
 { 350, 4 },
 { 420, 4 },
 { 490, 5 },
 { 560, 5 },  /* 120 bits */
 { 630, 5 },  /* 124 bits */
 { 700, 5 },  /* 128 bits */
};

uint128_t tinyqs128(uint128_t n)
{
  tiny_qs_params *params;
  s32 bits, fb_size, max_relations, bound, large_prime_mult, status = 0;
  const tiny_qs_config *config;
  tq_t factor, other, tmp;
  uint128_t result = 0;

  if (n < 4) return n-(n==1);
  if (!(n & 1)) return 2;
  if (n % 3 == 0) return 3;
  if (n % 5 == 0) return 5;
  tq_init(factor);
  tq_init(other);
  tq_init(tmp);
  tq_set_ui(tmp, n);

  /* make sure the input isn't a perfect square.
     We may also want to add a test for a perfect
     cube, but that's so unlikely it's probably
     not worth worrying about */

  if (tq_root(factor, tmp, 2) != 0)
    return factor->mag;

  /* start the initialization */

  params = (tiny_qs_params *)mpu_aligned_alloc(1, sizeof(*params), 16);
  memset(params, 0, sizeof(*params));
  params->input_n = n;
  init_tinyqs(params);
  tq_set_ui(params->n, n);
  params->num_full_relations = 0;
  params->partial_idx = 4 * MAX_RELATIONS_TINY - 1;
  params->poly_num = 0;

  bits = tq_sizeinbase(params->n, 2);
  find_multiplier_tiny(params);

  /* determine the factor base size and the
     number of primes in a polynomial A value */

  if (bits < 50)
    bits = 50;
  config = static_config + ((bits - 50) / 4);
  fb_size = config->fb_size;
  params->num_a_factors = config->num_poly_factors;

  /* build the factor base */

  fb_size = init_fb_tiny(params, fb_size, &result);
  if (result)
    goto done;
  params->target_relations = fb_size + EXTRA_RELATION_STEP_TINY;
  max_relations = fb_size + MAX_EXTRA_RELATIONS_TINY;

  /* compute the optimal A value */

  tq_sqrt(tmp, params->n);
  params->target_a = tq_get_d(tmp) * M_SQRT2 / SIEVE_SIZE_TINY;
  init_siqs_tiny(params);

  /* compute the large prime cutoff and the
     size of the fudge factor needed to account
     for it in the sieving cutoff */

  large_prime_mult = 15;
  bound = params->factor_base[fb_size - 1].prime;
  bound *= large_prime_mult;
  params->large_prime_max = bound;
  params->error_bits = (u32)(log(bound) / M_LN2 + 1);

  /* empty out the hashtable for partial relations */

  memset(params->partial_hash, 0xff, sizeof(params->partial_hash));

  /* Collect an initial relation set, then extend it only when all generated
     dependencies are trivial.  This keeps easy inputs cheap while giving
     difficult relation sets more independent dependencies. */

  while (!status && params->target_relations <= max_relations) {
    s32 attempt, attempts, sieve_status = 0;
    s32 exhausted;

    while (params->poly_num < MAX_POLY_TINY &&
           params->num_full_relations < params->target_relations) {
      sieve_status = sieve_next_poly_tiny(params);
      if (sieve_status != 0)
        break;
      params->poly_num++;
    }

    if (sieve_status == -3)
      goto done;
    exhausted = sieve_status != 0 || params->poly_num >= MAX_POLY_TINY;

    /* On exhaustion, use any final matrix that has guaranteed dependencies. */
    if (params->num_full_relations <= fb_size)
      break;

    attempts = exhausted || params->target_relations == max_relations ? 4 : 1;
    for (attempt = 0; attempt < attempts && !status; attempt++) {
      if (!solve_linear_system_tiny(params))
        break;
      status = find_factors_tiny(params, factor, other);
    }
    if (status || exhausted || params->target_relations >= max_relations)
      break;
    params->target_relations += EXTRA_RELATION_STEP_TINY;
  }

  if (_XS_get_verbose() >= 2)
    printf("tinyqs128 status %d, relations %d/%d, polynomials %d\n",
           status, params->num_full_relations,
           params->target_relations, params->poly_num);

  /* Return smaller of two factors */
  if (status && tq_cmp(factor, other) > 0) tq_swap(factor, other);

  if (status) {
    result = tq_gcd_value(factor->mag, n);
    if (result <= 1 || result >= n) result = 0;
  }

done:
  mpu_aligned_free(params);
  return result;
}

#endif /* HAVE_FACTOR128 */
