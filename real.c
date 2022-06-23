#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define FUNC_log2floor 1
#include "ptypes.h"
#include "util.h"
#include "sieve.h"
#include "real.h"
#include "mathl.h"

#define KAHAN_INIT(s) \
  LNV s ## _y, s ## _t; \
  LNV s ## _c = 0.0; \
  LNV s = 0.0;

#define KAHAN_SUM(s, term) \
  do { \
    s ## _y = (term) - s ## _c; \
    s ## _t = s + s ## _y; \
    s ## _c = (s ## _t - s) - s ## _y; \
    s = s ## _t; \
  } while (0)

static const unsigned short primes_tiny[] =
  {0,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
   101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,
   193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,
   293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,
   409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503};
#define NPRIMES_TINY (sizeof(primes_tiny)/sizeof(primes_tiny[0]))


/******************************************************************************/
/*                       REAL FUNCTIONS (EI,LI,etc.)                          */
/******************************************************************************/

/*
 * See:
 *  "Multiple-Precision Exponential Integral and Related Functions"
 *      by David M. Smith
 *  "On the Evaluation of the Complex-Valued Exponential Integral"
 *      by Vincent Pegoraro and Philipp Slusallek
 *  "Numerical Recipes" 3rd edition
 *      by William H. Press et al.
 *  "Rational Chevyshev Approximations for the Exponential Integral E_1(x)"
 *      by W. J. Cody and Henry C. Thacher, Jr.
 *  "High-precision Computation of Uniform Asymptotic Expansions for Special Functions"
 *      by Guillermo Navas-Palencia (2019)
 *
 * Any mistakes here are mine.  This code has not been rigorously verified.
 * Alternates: Navas-Palencia, Boost, MPFR, Pari/GP, Arb.
 *
 * We are trying to get close to maximum precision for all x with double, long
 * double, and quadmath.  Hence the rational Chebyshev approximations should
 * not be used with quadmath (unless they are are modified).
 *
 * Performance, i7-6700HQ, 2.6GHz, 1e-9 to 1000 step 0.001
 * range x > 0:
 *   0.22 microseconds, NV = double         max rel error 1.4e-14
 *   0.19 microseconds, NV = long double    max rel error 4.3e-17
 *  18.97 microseconds, NV = quad           max rel error 4.4e-32
 * range x < 0:
 *   0.18 microseconds, NV = double         max rel error 1.4e-14
 *   0.15 microseconds, NV = long double    max rel error 1.2e-17
 *   9.31 microseconds, NV = quad           max rel error 1.7e-32
 *
 * The maximum error is near the root 0.3725074...
 * The relative error profile for double precision is essentially identical
 * to the Navas-Palencia expintei(x) function.
 * Using long double on x86 improves the results with no time penalty.
 * Using quadmath gives improved results at a substantial time penalty.
 */

static LNV const euler_mascheroni = LNVCONST(0.57721566490153286060651209008240243104215933593992);
static LNV const li2 = LNVCONST(1.045163780117492784844588889194613136522615578151);

/* Rational Chebyshev approximation (Cody, Thacher), good for -1 < x < 0 */
static LNV _ei_chebyshev_neg(const LNV x) {
  static const LNV C6p[7] = { LNVCONST(-148151.02102575750838086),
                              LNVCONST( 150260.59476436982420737),
                              LNVCONST(  89904.972007457256553251),
                              LNVCONST(  15924.175980637303639884),
                              LNVCONST(   2150.0672908092918123209),
                              LNVCONST(    116.69552669734461083368),
                              LNVCONST(      5.0196785185439843791020) };
  static const LNV C6q[7] = { LNVCONST( 256664.93484897117319268),
                              LNVCONST( 184340.70063353677359298),
                              LNVCONST(  52440.529172056355429883),
                              LNVCONST(   8125.8035174768735759866),
                              LNVCONST(    750.43163907103936624165),
                              LNVCONST(     40.205465640027706061433),
                              LNVCONST(      1.0000000000000000000000) };
  LNV sumn = C6p[0]-x*(C6p[1]-x*(C6p[2]-x*(C6p[3]-x*(C6p[4]-x*(C6p[5]-x*C6p[6])))));
  LNV sumd = C6q[0]-x*(C6q[1]-x*(C6q[2]-x*(C6q[3]-x*(C6q[4]-x*(C6q[5]-x*C6q[6])))));
  return loglnv(-x) - sumn/sumd;
}
/* Cody / Thacher rational Chebyshev for x > 24 */
static LNV _ei_chebyshev_pos24(const LNV x) {
  static const LNV P2[10] = {
      LNVCONST( 1.75338801265465972390E02),
      LNVCONST(-2.23127670777632409550E02),
      LNVCONST(-1.81949664929868906455E01),
      LNVCONST(-2.79798528624305389340E01),
      LNVCONST(-7.63147701620253630855E00),
      LNVCONST(-1.52856623636929636839E01),
      LNVCONST(-7.06810977895029358836E00),
      LNVCONST(-5.00006640413131002475E00),
      LNVCONST(-3.00000000320981265753E00),
      LNVCONST( 1.00000000000000485503E00) };
  static const LNV Q2[9] = {
      LNVCONST( 3.97845977167414720840E04),
      LNVCONST( 3.97277109100414518365E00),
      LNVCONST( 1.37790390235747998793E02),
      LNVCONST( 1.17179220502086455287E02),
      LNVCONST( 7.04831847180424675988E01),
      LNVCONST(-1.20187763547154743238E01),
      LNVCONST(-7.99243595776339741065E00),
      LNVCONST(-2.99999894040324959612E00),
      LNVCONST( 1.99999999999048104167E00) };
  LNV invx = LNV_ONE / x, frac = 0.0;
  uint32_t n;
  for (n = 0; n <= 8; n++)
    frac = Q2[n] / (P2[n] + x + frac);
  frac += P2[9];
  return explnv(x) * (invx + invx*invx*frac);
}
#if 0
/* Continued fraction, good for x < -1 */
static LNV _ei_cfrac_neg(const LNV x) {
  LNV lc = 0, ld = LNV_ONE / (LNV_ONE - x);
  LNV val = ld * (-explnv(x));
  uint32_t n;
  for (n = 1; n <= 20000; n++) {
    LNV old, t, n2 = n * n;
    t = (LNV)(2*n + 1) - x;
    lc = LNV_ONE / (t - n2 * lc);
    ld = LNV_ONE / (t - n2 * ld);
    old = val;
    val *= ld/lc;
    if ( fabslnv(val-old) <= LNV_EPSILON*fabslnv(val) )
      break;
  }
  return val;
}
#endif
/* eint_v using Laguerre series, Navas-Palencia (2019). */
static LNV _eintv_laguerre_series(const LNV v, const LNV x) {
  LNV L_k = 1.0, L_k1 = x + v;
  LNV q, r, u = LNV_ONE, d = LNV_ONE;
  uint32_t k;
  KAHAN_INIT(sum);
  KAHAN_SUM(sum, (LNV_ONE/L_k1));
  for (k = 1; k < 500; k++) {
    u *= v + k - 1;
    d *= 1 + k;
    q = L_k1 * (x + 2*k + v) / (k + 1)  -  L_k * (k + v - 1) / (k + 1);
    r = u / (d * (q * L_k1));
    KAHAN_SUM(sum, r);
    L_k = L_k1;
    L_k1 = q;
    if (fabslnv(r) < 0.1 * LNV_EPSILON)
      break;
  }
  return sum * explnv(-x);
}
/* Convergent series for small negative x through medium positive x */
static LNV _ei_series_convergent(LNV const x) {
  LNV term, fact_n = x;
  uint32_t n;
  KAHAN_INIT(sum);
  for (n = 2; n <= 400; n++) {
    LNV invn = LNV_ONE / n;
    fact_n *= (LNV)x * invn;
    term = fact_n * invn;
    KAHAN_SUM(sum, term);
    /* printf("C  after adding %.20Lf, val = %.20Lf\n", term, sum); */
    if (fabslnv(term) < LNV_EPSILON*fabslnv(sum)) break;
  }
  KAHAN_SUM(sum, euler_mascheroni);
  KAHAN_SUM(sum, loglnv(fabslnv(x)));
  KAHAN_SUM(sum, x);
  return sum;
}
/* Asymptotic divergent series, for large positive x */
static LNV _ei_series_divergent(LNV const x) {
  LNV invx = LNV_ONE / x, term = invx;
  unsigned int n;
  KAHAN_INIT(sum);
  for (n = 2; n <= 400; n++) {
    LNV last_term = term;
    term = term * ( (LNV)n * invx );
    if (term < LNV_EPSILON*sum) break;
    if (term < last_term) {
      KAHAN_SUM(sum, term);
      /* printf("A  after adding %.20llf, sum = %.20llf\n", term, sum); */
    } else {
      KAHAN_SUM(sum, (-last_term/1.07) );
      /* printf("A  after adding %.20llf, sum = %.20llf\n", -last_term/1.07, sum); */
      break;
    }
  }
  KAHAN_SUM(sum, invx);
  KAHAN_SUM(sum, LNV_ONE);
  return explnv(x) * sum * invx;
}

NV Ei(NV x) {
  if (x == 0) croak("Invalid input to ExponentialIntegral:  x must be != 0");
  /* Protect against messed up rounding modes */
  if (x >=  12000) return INFINITY;
  if (x <= -12000) return 0;

  if (x < 0) {
    if (x >= -1.0 && !LNV_IS_QUAD) return _ei_chebyshev_neg(x);
    else if (x < -0.80)            return -_eintv_laguerre_series(1, -x);
    else                           return _ei_series_convergent(x);
  } else {
    if (x < (-2 * loglnv(LNV_EPSILON)))         return _ei_series_convergent(x);
    if (x >= 24 && (!LNV_IS_QUAD || x <= 43.2)) return _ei_chebyshev_pos24(x);
    else                                        return _ei_series_divergent(x);
  }
}

NV Li(NV x) {
  if (x == 0) return 0;
  if (x == 1) return -INFINITY;
  if (x == 2) return li2;
  if (x < 0) croak("Invalid input to LogarithmicIntegral:  x must be >= 0");
  if (x >= NV_MAX) return INFINITY;

  /* Calculate directly using Ramanujan's series. */
  if (x > 1) {
    const LNV logx = loglnv(x);
    LNV sum = 0, inner_sum = 0, old_sum, factorial = 1, power2 = 1;
    LNV q, p = -1;
    int k = 0, n = 0;

    for (n = 1, k = 0; n < 200; n++) {
      factorial *= n;
      p *= -logx;
      q = factorial * power2;
      power2 *= 2;
      for (; k <= (n - 1) / 2; k++)
        inner_sum += LNV_ONE / (2 * k + 1);
      old_sum = sum;
      sum += (p / q) * inner_sum;
      if (fabslnv(sum - old_sum) <= LNV_EPSILON) break;
    }
    return euler_mascheroni + loglnv(logx) + sqrtlnv(x) * sum;
  }

  return Ei(loglnv(x));
}

static long double ld_inverse_li(long double lx) {
  int i;
  long double t, term, old_term = 0;
  /* Iterate Halley's method until error grows. */
  t = (lx <= 2)  ?  2  :  lx * logl(lx);
  for (i = 0; i < 4; i++) {
    long double dn = Li(t) - lx;
    term = dn*logl(t) / (1.0L + dn/(2*t));
    if (i > 0 && fabsl(term) >= fabsl(old_term)) { t -= term/4; break; }
    old_term = term;
    t -= term;
  }
  return t;
}

UV inverse_li(UV x) {
  UV r, i;
  long double lx = (long double) x;

  if (x <= 2) return x + (x > 0);
  r = (UV) ceill( ld_inverse_li(lx) );
  /* Meet our more stringent goal of an exact answer. */
  i = (x > 4e16) ? 2048 : 128;
  if (Li(r-1) >= lx) {
    while (Li(r-i) >= lx) r -= i;
    for (i = i/2; i > 0; i /= 2)
      if (Li(r-i) >= lx) r -= i;
  } else {
    while (Li(r+i-1) < lx) r += i;
    for (i = i/2; i > 0; i /= 2)
      if (Li(r+i-1) < lx) r += i;
  }
  return r;
}

static long double ld_inverse_R(long double lx) {
  int i;
  long double t, dn, term, old_term = 0;

  /* Rough estimate */
  if (lx <= 3.5) {
    t = lx + 2.24*(lx-1)/2;
  } else {
    t = lx * logl(lx);
    if      (lx <   50) { t *= 1.2; }
    else if (lx < 1000) { t *= 1.15; }
    else {   /* use inverse Li (one iteration) for first inverse R approx */
      dn = Li(t) - lx;
      term = dn * logl(t) / (1.0L + dn/(2*t));
      t -= term;
    }
  }
  /* Iterate 1-n rounds of Halley, usually only 3 needed. */
  for (i = 0; i < 100; i++) {
    dn = RiemannR(t, 1e-12) - lx;
#if 1  /* Use f(t) = li(t) for derivatives */
    term = dn * logl(t) / (1.0L + dn/(2*t));
#else  /* Use f(t) = li(t) - li(sqrt(t))/2 for derivatives */
    long double logt = logl(t);
    long double sqrtt = sqrtl(t);
    long double FA = 2 * sqrtt * logt;
    long double FB = 2 * sqrtt - 1;
    long double ifz = FA / FB;
    long double iffz = (logt - 2*FB) / (2 * sqrtt * FA * FA * FA * FA);
    term = dn * ifz * (1.0L - dn * iffz);
#endif
    if (i > 0 && fabsl(term) >= fabsl(old_term)) { t -= term/4; break; }
    old_term = term;
    t -= term;
  }
  return t;
}

UV inverse_R(UV x) {
  if (x < 2) return x + (x > 0);
  return (UV) ceill( ld_inverse_R( (long double) x) );
}


/*
 * Storing the first 10-20 Zeta values makes sense.  Past that it is purely
 * to avoid making the call to generate them ourselves.  We could cache the
 * calculated values. These all have 1 subtracted from them.  */
static const long double riemann_zeta_table[] = {
  0.6449340668482264364724151666460251892L,  /* zeta(2) */
  0.2020569031595942853997381615114499908L,
  0.0823232337111381915160036965411679028L,
  0.0369277551433699263313654864570341681L,
  0.0173430619844491397145179297909205279L,
  0.0083492773819228268397975498497967596L,
  0.0040773561979443393786852385086524653L,
  0.0020083928260822144178527692324120605L,
  0.0009945751278180853371459589003190170L,
  0.0004941886041194645587022825264699365L,
  0.0002460865533080482986379980477396710L,
  0.0001227133475784891467518365263573957L,
  0.0000612481350587048292585451051353337L,
  0.0000305882363070204935517285106450626L,
  0.0000152822594086518717325714876367220L,
  0.0000076371976378997622736002935630292L,  /* zeta(17)  Past here all we're */
  0.0000038172932649998398564616446219397L,  /* zeta(18)  getting is speed.   */
  0.0000019082127165539389256569577951013L,
  0.0000009539620338727961131520386834493L,
  0.0000004769329867878064631167196043730L,
  0.0000002384505027277329900036481867530L,
  0.0000001192199259653110730677887188823L,
  0.0000000596081890512594796124402079358L,
  0.0000000298035035146522801860637050694L,
  0.0000000149015548283650412346585066307L,
  0.0000000074507117898354294919810041706L,
  0.0000000037253340247884570548192040184L,
  0.0000000018626597235130490064039099454L,
  0.0000000009313274324196681828717647350L,
  0.0000000004656629065033784072989233251L,
  0.0000000002328311833676505492001455976L,
  0.0000000001164155017270051977592973835L,
  0.0000000000582077208790270088924368599L,
  0.0000000000291038504449709968692942523L,
  0.0000000000145519218910419842359296322L,
  0.0000000000072759598350574810145208690L,
  0.0000000000036379795473786511902372363L,
  0.0000000000018189896503070659475848321L,
  0.0000000000009094947840263889282533118L,
  0.0000000000004547473783042154026799112L,
  0.0000000000002273736845824652515226821L,
  0.0000000000001136868407680227849349105L,
  0.0000000000000568434198762758560927718L,
  0.0000000000000284217097688930185545507L,
  0.0000000000000142108548280316067698343L,
  0.00000000000000710542739521085271287735L,
  0.00000000000000355271369133711367329847L,
  0.00000000000000177635684357912032747335L,
  0.000000000000000888178421093081590309609L,
  0.000000000000000444089210314381336419777L,
  0.000000000000000222044605079804198399932L,
  0.000000000000000111022302514106613372055L,
  0.0000000000000000555111512484548124372374L,
  0.0000000000000000277555756213612417258163L,
  0.0000000000000000138777878097252327628391L,
};
#define NPRECALC_ZETA (sizeof(riemann_zeta_table)/sizeof(riemann_zeta_table[0]))

/* Riemann Zeta on the real line, with 1 subtracted.
 * Compare to Math::Cephes zetac.  Also zeta with q=1 and subtracting 1.
 *
 * The Cephes zeta function uses a series (2k)!/B_2k which converges rapidly
 * and has a very wide range of values.  We use it here for some values.
 *
 * Note: Calculations here are done on long doubles and we try to generate as
 *       much accuracy as possible.  They will get returned to Perl as an NV,
 *       which is typically a 64-bit double with 15 digits.
 *
 * For values 0.5 to 5, this code uses the rational Chebyshev approximation
 * from Cody and Thacher.  This method is extraordinarily fast and very
 * accurate over its range (slightly better than Cephes for most values).  If
 * we had quad floats, we could use the 9-term polynomial.
 */
long double ld_riemann_zeta(long double x) {
  int i;

  if (x < 0)  croak("Invalid input to RiemannZeta:  x must be >= 0");
  if (x == 1) return INFINITY;

  if (x == (unsigned int)x) {
    int k = x - 2;
    if ((k >= 0) && (k < (int)NPRECALC_ZETA))
      return riemann_zeta_table[k];
  }

  /* Cody / Thacher rational Chebyshev approximation for small values */
  if (x >= 0.5 && x <= 5.0) {
    static const long double C8p[9] = { 1.287168121482446392809e10L,
                                        1.375396932037025111825e10L,
                                        5.106655918364406103683e09L,
                                        8.561471002433314862469e08L,
                                        7.483618124380232984824e07L,
                                        4.860106585461882511535e06L,
                                        2.739574990221406087728e05L,
                                        4.631710843183427123061e03L,
                                        5.787581004096660659109e01L };
    static const long double C8q[9] = { 2.574336242964846244667e10L,
                                        5.938165648679590160003e09L,
                                        9.006330373261233439089e08L,
                                        8.042536634283289888587e07L,
                                        5.609711759541920062814e06L,
                                        2.247431202899137523543e05L,
                                        7.574578909341537560115e03L,
                                       -2.373835781373772623086e01L,
                                        1.000000000000000000000L    };
    long double sumn = C8p[0]+x*(C8p[1]+x*(C8p[2]+x*(C8p[3]+x*(C8p[4]+x*(C8p[5]+x*(C8p[6]+x*(C8p[7]+x*C8p[8])))))));
    long double sumd = C8q[0]+x*(C8q[1]+x*(C8q[2]+x*(C8q[3]+x*(C8q[4]+x*(C8q[5]+x*(C8q[6]+x*(C8q[7]+x*C8q[8])))))));
    long double sum = (sumn - (x-1)*sumd) / ((x-1)*sumd);
    return sum;
  }

  if (x > 17000.0)
    return 0.0;

#if 0
  {
    KAHAN_INIT(sum);
    /* Simple defining series, works well. */
    for (i = 5; i <= 1000000; i++) {
      long double term = powl(i, -x);
      KAHAN_SUM(sum, term);
      if (term < LDBL_EPSILON*sum) break;
    }
    KAHAN_SUM(sum, powl(4, -x) );
    KAHAN_SUM(sum, powl(3, -x) );
    KAHAN_SUM(sum, powl(2, -x) );
    return sum;
  }
#endif

  /* The 2n!/B_2k series used by the Cephes library. */
  {
    /* gp/pari:
     *   for(i=1,13,printf("%.38g\n",(2*i)!/bernreal(2*i)))
     * MPU:
     *   use bignum;
     *   say +(factorial(2*$_)/bernreal(2*$_))->bround(38) for 1..13;
     */
    static const long double A[] = {
      12.0L,
     -720.0L,
      30240.0L,
     -1209600.0L,
      47900160.0L,
     -1892437580.3183791606367583212735166425L,
      74724249600.0L,
     -2950130727918.1642244954382084600497650L,
      116467828143500.67248729113000661089201L,
     -4597978722407472.6105457273596737891656L,
      181521054019435467.73425331153534235290L,
     -7166165256175667011.3346447367083352775L,
      282908877253042996618.18640556532523927L,
    };
    long double a, b, s, t;
    const long double w = 10.0;
    s = 0.0;
    b = 0.0;
    for (i = 2; i < 11; i++) {
      b = powl( i, -x );
      s += b;
      if (fabsl(b) < fabsl(LDBL_EPSILON * s))
        return s;
    }
    s = s + b*w/(x-1.0) - 0.5 * b;
    a = 1.0;
    for (i = 0; i < 13; i++) {
      long double k = 2*i;
      a *= x + k;
      b /= w;
      t = a*b/A[i];
      s = s + t;
      if (fabsl(t) < fabsl(LDBL_EPSILON * s))
        break;
      a *= x + k + 1.0;
      b /= w;
    }
    return s;
  }
}

long double RiemannR(long double x, long double eps) {
  long double part_term, term, flogx, ki, old_sum;
  unsigned int k;
  KAHAN_INIT(sum);

  if (x <= 0) croak("Invalid input to RiemannR:  x must be > 0");
  if (eps < LDBL_EPSILON) eps = LDBL_EPSILON;

  if (x > 1e19) {
    const signed char* amob = range_moebius(0, 100);
    KAHAN_SUM(sum, Li(x));
    for (k = 2; k <= 100; k++) {
      if (amob[k] == 0) continue;
      ki = 1.0L / (long double) k;
      part_term = powl(x,ki);
      if (part_term > LDBL_MAX) return INFINITY;
      term = amob[k] * ki * Li(part_term);
      old_sum = sum;
      KAHAN_SUM(sum, term);
      if (fabsl(sum - old_sum) <= eps) break;
    }
    Safefree(amob);
    return sum;
  }

  KAHAN_SUM(sum, 1.0);
  flogx = logl(x);
  part_term = 1;

  for (k = 1; k <= 10000; k++) {
    ki = (k-1 < NPRECALC_ZETA) ? riemann_zeta_table[k-1] : ld_riemann_zeta(k+1);
    part_term *= flogx / k;
    term = part_term / (k + k * ki);
    old_sum = sum;
    KAHAN_SUM(sum, term);
    /* printf("R %5d after adding %.18Lg, sum = %.19Lg (%Lg)\n", k, term, sum, fabsl(sum-old_sum)); */
    if (fabsl(sum - old_sum) <= eps) break;
  }

  return sum;
}

/* Options for LambertW initial approximation:
 *
 * - Four regions, we used before:
 *      Pade(3,2), Winitzki 2003, Vargas 2013, Corless 1993
 *   Has issues near -1/e but ok around zero.
 *
 * - Iacono and Boyd (2017).  Very simple function over whole range.
 *   Doesn't work right very near -1/e and around zero.
 *
 * - Vazquez-Leal et al. (2019).  Divides into four regions, power
 *   series for each.  Great results.  Also has issues near -1/e and zero.
 *
 * We use known solutions for near -1/e and around zero.  See Fukushima (2013)
 * and Johannson (2017,2020) for lots of discussion and solutions.
 * Use Vazquez-Leal (PSEM Approximations) for the rest.
 */
static long double _lambertw_approx(long double x) {
  long double w, k1, k2, k3;

  if (x < -0.312) {
    /* Use Puiseux series, e.g. Verberic 2009, Boost, Johannson (2020). */
    /* Near the branch point.  See Fukushima (2013) section 2.5. */
    k2 = 2.0L * (1.0L + 2.7182818284590452353603L * x);
    if (k2 <= 0) return -1.0L + 1*LDBL_EPSILON;
    k1 = sqrtl(k2);
    w = -1.0L + (1.0L + (-1.0L/3.0L + (11.0L/72.0L + (-43.0L/540.0L + (769.0L/17280.0L + (-221.0L/8505.0L + (680863.0L/43545600.0L + (-1963.0L/204120.0L + 226287557.0L/37623398400.0L
    * k1) * k1) * k1) * k1) * k1) * k1) * k1) * k1) * k1;

  } else if (x > -0.14 && x < 0.085) {
    /* Around zero.  See Fukushima (2013) section 2.6. */
    w = (1.0L + (-1.0L + (3.0L/2.0L + (-8.0L/3.0L + (125.0L/24.0L + (-54.0L/5.0L + (16807.0L/720.0L + (-16384.0L/315.0L + 531441.0L/4480.0L
        * x) * x) * x) * x) * x) * x) * x) * x) * x;

  } else if (x < 1) {
    /* This and the rest from Vazquez-Leal et al. (2019). */
    k1 = sqrtl(1.0L + 2.7182818284590452353603L * x);
    k2 = 0.33333333333333333333333L + 0.70710678118654752440084L / k1 - 0.058925565098878960366737L * k1 +
         (x + 0.36787944117144L) * (0.050248489761611L + (0.11138904851051 + 0.040744556245195L * x) * x)
         /
         (1.0L + (2.7090878606183L + (1.5510922597820L + 0.095477712183841L * x) * x) * x);
    w = -(k2-1)/k2;

  } else if (x < 40) {
    k1 = 1.0L + (5.950065500550155L + (13.96586471370701L + (10.52192021050505L + (3.065294254265870L + 0.1204576876518760L * x) * x) * x) * x) * x;
    w = 0.1600049638651493L * logl(k1);

  } else if (x < 20000) {
    k1 = 1.0L + (-3.16866642511229e11L + (3.420439800038598e10L +
         (-1.501433652432257e9L + (3.44887729947585e7L + (-4.453783741137856e5L +
         (3257.926478908996L + (-10.82545259305382L + (0.6898058947898353e-1L +
         0.4703653406071575e-4L * x) * x) * x) * x) * x) * x) * x) * x) * x;
    w = 0.9898045358731312e-1L * logl(k1);

  } else {
    k1 = 1.0L / (1.0L + logl(1.0L + x));
    k2 = 1.0L / k1;
    k3 = logl(k2);
    w = k2-1-k3+(1+k3+(-1/2+(1/2)*k3*k3 +(-1/6+(-1+(-1/2+
        (1/3) * k3) * k3) * k3) * k1) * k1) * k1;
  }
  return w;
}

NV lambertw(NV x) {
  long double w;
  int i;

  if (x < -0.36787944117145L)
    croak("Invalid input to LambertW:  x must be >= -1/e");
  if (x == 0.0L) return 0.0L;

  /* Estimate initial value */
  w = _lambertw_approx(x);

  /* TODO: this section might not be best for quad precision */
  /* If input is too small, return .99999.... */
  /* if (w <= -1.0L) return -1.0L + LDBL_EPSILON; */
  /* For very small inputs, don't iterate, return approx directly. */
  if (x < -0.36768) return w;

#if 0  /* Halley */
  long double lastw = w;
  for (i = 0; i < 100; i++) {
    long double ew = expl(w);
    long double wew = w * ew;
    long double wewx = wew - x;
    long double w1 = w + 1;
    w = w - wewx / (ew * w1 - (w+2) * wewx/(2*w1));
    if (w != 0.0L && fabsl((w-lastw)/w) <= 8*LDBL_EPSILON) break;
    lastw = w;
  }
#else  /* Fritsch, see Veberic 2009.  1-2 iterations are enough. */
  for (i = 0; i < 6 && w != 0.0L; i++) {
    long double w1 = 1 + w;
    long double zn = logl((long double)x/w) - w;
    long double qn = 2 * w1 * (w1+(2.0L/3.0L)*zn);
    long double en = (zn/w1) * (qn-zn)/(qn-2.0L*zn);
    /* w *= 1.0L + en;  if (fabsl(en) <= 16*LDBL_EPSILON) break; */
    long double wen = w * en;
    if (isnan(wen)) return 0;
    w += wen;
    if (fabsl(wen) <= 64*LDBL_EPSILON) break;
  }
#endif

#if LNV_IS_QUAD /* For quadmath, one high precision correction */
  if (w != LNV_ZERO) {
    LNV lw = w;
    LNV w1 = LNV_ONE + lw;
    LNV zn = loglnv((LNV)x/lw) - lw;
    LNV qn = LNVCONST(2.0) * w1 * (w1+(LNVCONST(2.0)/LNVCONST(3.0))*zn);
    LNV en = (zn/w1) * (qn-zn)/(qn-LNVCONST(2.0)*zn);
    return lw + lw * en;
  }
#endif

  /* With long double = 64-bit double, we have 15 digits precision
   * near the branch point, and 16 over the rest of the range.
   * With long double = x86 extended precision, we have over 17 digits
   * over the entire range.
   * Correcting to the exact LDBL_EPSILON does not improve this. */

  return w;
}


/******************************************************************************/
/*                        Chebyshev PSI / THETA                               */
/******************************************************************************/

NV chebyshev_psi(UV n)
{
  UV k;
  KAHAN_INIT(sum);

  for (k = log2floor(n); k > 0; k--) {
    KAHAN_SUM(sum, chebyshev_theta(rootint(n,k)));
  }
  return sum;
}

#if BITS_PER_WORD == 64
typedef struct {
  UV n;
  LNV theta;
} cheby_theta_t;
static const cheby_theta_t _cheby_theta[] = { /* >= quad math precision */
  { UVCONST(      67108864),LNVCONST(    67100507.6357700963903836828562472350035880) },
  { UVCONST(     100000000),LNVCONST(    99987730.0180220043832124342600487053812729) },
  { UVCONST(     134217728),LNVCONST(   134204014.5735572091791081610859055728165544) },
  { UVCONST(     268435456),LNVCONST(   268419741.6134308193112682817754501071404173) },
  { UVCONST(     536870912),LNVCONST(   536842885.8045763840625719515011160692495056) },
  { UVCONST(    1000000000),LNVCONST(   999968978.5775661447991262386023331863364793) },
  { UVCONST(    1073741824),LNVCONST(  1073716064.8860663337617909073555831842945484) },
  { UVCONST(    2147483648),LNVCONST(  2147432200.2475857676814950053003448716360822) },
  { UVCONST(    4294967296),LNVCONST(  4294889489.1735446386752045191908417183337361) },
  { UVCONST(    8589934592),LNVCONST(  8589863179.5654263491545135406516173629373070) },
  { UVCONST(   10000000000),LNVCONST(  9999939830.6577573841592219954033850595228736) },
  { UVCONST(   12884901888),LNVCONST( 12884796620.4324254952601520445848183460347362) },
  { UVCONST(   17179869184),LNVCONST( 17179757715.9924077567777285147574707468995695) },
  { UVCONST(   21474836480),LNVCONST( 21474693322.0998273969188369449626287713082943) },
  { UVCONST(   25769803776),LNVCONST( 25769579799.3751535467593954636665656772211515) },
  { UVCONST(   30064771072),LNVCONST( 30064545001.2305211029215168703433831598544454) },
  { UVCONST(   34359738368),LNVCONST( 34359499180.0126643918259085362039638823175054) },
  { UVCONST(   51539607552),LNVCONST( 51539356394.9531019037592855639826469993402730) },
  { UVCONST(   68719476736),LNVCONST( 68719165213.6369838785284711480925219076501720) },
  { UVCONST(   85899345920),LNVCONST( 85899083852.3471545629838432726841470626910905) },
  { UVCONST(  100000000000),LNVCONST( 99999737653.1074446948519125729820679772770146) },
  { UVCONST(  103079215104),LNVCONST(103079022007.113299711630969211422868856259124) },
  { UVCONST(  120259084288),LNVCONST(120258614516.787336970535750737470005730125261) },
  { UVCONST(  137438953472),LNVCONST(137438579206.444595884982301543904849253294539) },
  { UVCONST(  171798691840),LNVCONST(171798276885.585945657918751085729734540334501) },
  { UVCONST(  206158430208),LNVCONST(206158003808.160276853604927822609009916573462) },
  { UVCONST(  240518168576),LNVCONST(240517893445.995868018331936763125264759516048) },
  { UVCONST(  274877906944),LNVCONST(274877354651.045354829956619821889825596300686) },
  { UVCONST(  309237645312),LNVCONST(309237050379.850690561796126460858271984023198) },
  { UVCONST(  343597383680),LNVCONST(343596855806.595496630500062749631211394707114) },
  { UVCONST(  377957122048),LNVCONST(377956498560.227794386327526022452943941537993) },
  { UVCONST(  412316860416),LNVCONST(412316008796.349553568121442261222464590518293) },
  { UVCONST(  446676598784),LNVCONST(446675972485.936512329625489223180824947531484) },
  { UVCONST(  481036337152),LNVCONST(481035608287.572961376833237046440177624505864) },
  { UVCONST(  515396075520),LNVCONST(515395302740.633513931333424447688399032397200) },
  { UVCONST(  549755813888),LNVCONST(549755185085.539613556787409928561107952681488) },
  { UVCONST(  584115552256),LNVCONST(584115015741.698143680148976236958207248900725) },
  { UVCONST(  618475290624),LNVCONST(618474400071.621528348965919774195984612254220) },
  { UVCONST(  652835028992),LNVCONST(652834230470.583317059774197550110194348469358) },
  { UVCONST(  687194767360),LNVCONST(687193697328.927006867624832386534836384752774) },
  { UVCONST(  721554505728),LNVCONST(721553211683.605313067593521060195071837766347) },
  { UVCONST(  755914244096),LNVCONST(755913502349.878525212441903698096011352015192) },
  { UVCONST(  790273982464),LNVCONST(790273042590.053075430445971969285969445183076) },
  { UVCONST(  824633720832),LNVCONST(824633080997.428352876758261549475609957696369) },
  { UVCONST(  858993459200),LNVCONST(858992716288.318498931165663742671579465316192) },
  { UVCONST(  893353197568),LNVCONST(893352235882.851072417721659027263613727927680) },
  { UVCONST(  927712935936),LNVCONST(927711881043.628817668337317445143018372892386) },
  { UVCONST(  962072674304),LNVCONST(962071726126.508938539006575212272731584070786) },
  { UVCONST(  996432412672),LNVCONST(996431411588.361462717402562171913706963939018) },
  { UVCONST( 1099511627776),LNVCONST(1099510565082.05800550569923209414874779035972) },
  { UVCONST( 1168231104512),LNVCONST(1168230478726.83399452743801182220790107593115) },
  { UVCONST( 1236950581248),LNVCONST(1236949680081.02610603189530371762093291521116) },
  { UVCONST( 1305670057984),LNVCONST(1305668780900.04255251887970870257110498423202) },
  { UVCONST( 1374389534720),LNVCONST(1374388383792.63751003694755359184583212193880) },
  { UVCONST( 1443109011456),LNVCONST(1443107961091.80955496949174183091839841371227) },
  { UVCONST( 1511828488192),LNVCONST(1511827317611.91227277802426032456922797572429) },
  { UVCONST( 1580547964928),LNVCONST(1580546753969.30607547506449941085747942395437) },
  { UVCONST( 1649267441664),LNVCONST(1649265973878.75361554498682516738256005501353) },
  { UVCONST( 1717986918400),LNVCONST(1717985403764.24562741452793071287954107946922) },
  { UVCONST( 1786706395136),LNVCONST(1786704769212.04241689416220650800274263053933) },
  { UVCONST( 1855425871872),LNVCONST(1855425013030.54920163513184322741954734357404) },
  { UVCONST( 1924145348608),LNVCONST(1924143701943.02957992419280264060220278182021) },
  { UVCONST( 1992864825344),LNVCONST(1992863373568.84039296068619447120308124302086) },
  { UVCONST( 2061584302080),LNVCONST(2061583632335.91985095534685076604018573279204) },
  { UVCONST( 2130303778816),LNVCONST(2113122935598.01727180199783433992649406589029) },
  { UVCONST( 2199023255552),LNVCONST(2199021399611.18488312543276191461914978761981) },
  { UVCONST( 2267742732288),LNVCONST(2267740947106.05038218811506263712808318234921) },
  { UVCONST( 2336462209024),LNVCONST(2336460081480.34962633829077377680844065198307) },
  { UVCONST( 2405181685760),LNVCONST(2405179969505.38642629423585641169740223940265) },
  { UVCONST( 2473901162496),LNVCONST(2473899311193.37872375168104562948639924654178) },
  { UVCONST( 2542620639232),LNVCONST(2542619362554.88893589220737167756411653816418) },
  { UVCONST( 2611340115968),LNVCONST(2611338370515.94936514022501267847930999670553) },
  { UVCONST( 2680059592704),LNVCONST(2680057722824.52981820001574883706268873541107) },
  { UVCONST( 2748779069440),LNVCONST(2748777610452.18903407570165081726781627254885) },
  { UVCONST( 2817498546176),LNVCONST(2817497017165.31924616507392971415494161401775) },
  { UVCONST( 2886218022912),LNVCONST(2886216579432.32232322707222172612181994322081) },
  { UVCONST( 2954937499648),LNVCONST(2954936100812.97301730406598982753121204977388) },
  { UVCONST( 3023656976384),LNVCONST(3023654789503.82041452274471455184651411931920) },
  { UVCONST( 3298534883328),LNVCONST(3298533215621.76606493931157388037915263658637) },
  { UVCONST( 3573412790272),LNVCONST(3573411344351.74163523704886736624674718378131) },
  { UVCONST( 3848290697216),LNVCONST(3848288415701.82534219216958446478503907262807) },
  { UVCONST( 4123168604160),LNVCONST(4123166102085.86116301709394219323327831487542) },
  { UVCONST( 4398046511104),LNVCONST(4398044965678.05143041707871320554940671182665) },
  { UVCONST( 4672924418048),LNVCONST(4672922414672.04998927945349278916525727295687) },
  { UVCONST( 4947802324992),LNVCONST(4947800056419.04384937181159608905993450182729) },
  { UVCONST( 5222680231936),LNVCONST(5222678728087.69487334278665824384732845008859) },
  { UVCONST( 5497558138880),LNVCONST(5497555766573.55159115560501595606332808978878) },
  { UVCONST( 5772436045824),LNVCONST(5772433560746.27053256770924553245647027548204) },
  { UVCONST( 6047313952768),LNVCONST(6047310750621.24497633828761530843255989494448) },
  { UVCONST( 6322191859712),LNVCONST(6322189275338.39747421237532473168802646234745) },
  { UVCONST( 6597069766656),LNVCONST(6579887620000.56226807898107616294821989189226) },
  { UVCONST( 6871947673600),LNVCONST(6871945430474.61791600096091374271286154432006) },
  { UVCONST( 7146825580544),LNVCONST(7146823258390.34361980709600216319269118247416) },
  { UVCONST( 7421703487488),LNVCONST(7421700443390.35536080251964387835425662360121) },
  { UVCONST( 7696581394432),LNVCONST(7696578975137.73249441643024336954233783264803) },
  { UVCONST( 7971459301376),LNVCONST(7971457197928.90863708984184849978605273042512) },
  { UVCONST( 8246337208320),LNVCONST(8246333982863.77146812177727648999195989358960) },
  { UVCONST( 8521215115264),LNVCONST(8529802085075.55635100929751669785228592926043) },
  { UVCONST( 8796093022208),LNVCONST(8796089836425.34909684634625258535266362465034) },
  { UVCONST( 9345848836096),LNVCONST(9345845828116.77456046925508587313) },
  { UVCONST( 9895604649984),LNVCONST(9895601077915.26821447819584407150) },
  { UVCONST(10000000000000),LNVCONST(9999996988293.03419965318214160284) },
  { UVCONST(15000000000000),LNVCONST(14999996482301.7098815115045166858) },
  { UVCONST(20000000000000),LNVCONST(19999995126082.2286880312461318496) },
  { UVCONST(25000000000000),LNVCONST(24999994219058.4086216020475916538) },
  { UVCONST(30000000000000),LNVCONST(29999995531389.8454274046657200568) },
  { UVCONST(35000000000000),LNVCONST(34999992921190.8049427456456479005) },
  { UVCONST(40000000000000),LNVCONST(39999993533724.3168289589273168844) },
  { UVCONST(45000000000000),LNVCONST(44999993567606.9795798378256194424) },
  { UVCONST(50000000000000),LNVCONST(49999992543194.2636545758235373677) },
  { UVCONST(55000000000000),LNVCONST(54999990847877.2435105757522625171) },
  { UVCONST(60000000000000),LNVCONST(59999990297033.6261976055811111726) },
  { UVCONST(65000000000000),LNVCONST(64999990861395.5522142429859245014) },
  { UVCONST(70000000000000),LNVCONST(69999994316409.8717306521862685981) },
  { UVCONST(75000000000000),LNVCONST(74999990126219.8344899338374090165) },
  { UVCONST(80000000000000),LNVCONST(79999990160858.3042387288372250950) },
  { UVCONST(85000000000000),LNVCONST(84999987096970.5915212896832780715) },
  { UVCONST(90000000000000),LNVCONST(89999989501395.0738966599857919767) },
  { UVCONST(95000000000000),LNVCONST(94999990785908.6672552042792168144) },
  { UVCONST(100000000000000),LNVCONST(99999990573246.9785384070303475639) },
};
#define NCHEBY_VALS (sizeof(_cheby_theta)/sizeof(_cheby_theta[0]))
#endif

NV chebyshev_theta(UV n)
{
  uint16_t i = 0;
  UV tp, startn, seg_base, seg_low, seg_high;
  unsigned char* segment;
  void* ctx;
  LNV initial_sum, prod = LNV_ONE;
  KAHAN_INIT(sum);

  if (n < 500) {
    for (i = 1;  (tp = primes_tiny[i]) <= n; i++) {
      KAHAN_SUM(sum, loglnv(tp));
    }
    return sum;
  }

#if defined NCHEBY_VALS
  if (n >= _cheby_theta[0].n) {
    for (i = 1; i < NCHEBY_VALS; i++)
      if (n < _cheby_theta[i].n)
        break;
    startn = _cheby_theta[i-1].n;
    initial_sum = _cheby_theta[i-1].theta;
  } else
#endif
  {
    KAHAN_SUM(sum, loglnv(2*3*5*7*11*13));
    startn = 17;
    initial_sum = 0;
  }

  ctx = start_segment_primes(startn, n, &segment);
#if 0
  while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
    START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high ) {
      KAHAN_SUM(sum, loglnv(p));
    } END_DO_FOR_EACH_SIEVE_PRIME
  }
#else
  while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
    START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high ) {
      prod *= (LNV) p;
      if (++i >= (LNV_IS_QUAD ? 64 : 8)) {
        KAHAN_SUM(sum, loglnv(prod));
        prod = LNV_ONE;
        i = 0;
      }
    } END_DO_FOR_EACH_SIEVE_PRIME
  }
  if (prod > 1.0) { KAHAN_SUM(sum, loglnv(prod));  prod = LNV_ONE; }
#endif
  end_segment_primes(ctx);

  if (initial_sum > 0) KAHAN_SUM(sum, initial_sum);
  return sum;
}


/******************************************************************************/
/*                                Other                                       */
/******************************************************************************/

#if 0
/* This is the de Bruijn approximation, not exact! */
static long double dickman_rho(long double u) {
  int i;
  long double zeta;

  if (u <= 1) return 1;
  if (u <= 2) return 1-logl(u);

  /* Also see Granville (2008) Smooth numbers */

  /* Calculate zeta.  See Bach and Sorenson (2013) page 10 */
  zeta = 2*(u-1);
  for (i = 0; i < 7; i++) {
    long double uz1 = 1 + u*zeta;
    zeta = zeta - ( (zeta-logl(uz1))*uz1 ) / (uz1-u);
  }
  /* Alternately:  zeta = -1/u - LambertW1(-exp(-1/u)/u) */

  return expl(-u*zeta+Ei(zeta)) / (zeta * sqrtl(2*3.1415926535*u));
}
#endif


