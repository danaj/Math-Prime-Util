#ifndef MPU_MATHL_H
#define MPU_MATHL_H

#include <float.h>
#include <math.h>

/* Use long double to get a little more precision when we're calculating the
 * math functions -- especially those calculated with a series.  Long double
 * is defined in C89 (ISO C).  Note that 'long double' on many platforms is
 * identical to 'double so it may buy us nothing.  But it's worth trying.
 *
 * While the type was in C89, math functions using it are in C99.  A few
 * systems lacked support for years (e.g. NetBSD and FreeBSD).
 */

#undef HAS_LDBL_FUNCS
#if defined(__FreeBSD__) && (__FreeBSD_version < 1000034)
  /* Added in 2013. */
#elif _MSC_VER || defined(__IBMC__) || defined(__IBMCPP__)
#  define HAS_LDBL_FUNCS 1
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#  define HAS_LDBL_FUNCS 1
#else
  /* We don't have them */
#endif

#ifdef HAS_LDBL_FUNCS
  /* math.h should give us these as functions or macros.
   *
   *  extern long double fabsl(long double);
   *  extern long double floorl(long double);
   *  extern long double ceill(long double);
   *  extern long double sqrtl(long double);
   *  extern long double powl(long double, long double);
   *  extern long double expl(long double);
   *  extern long double logl(long double);
   *  extern long double log1pl(long double);
   */
#else
  #define fabsl(x)    (long double) fabs( (double) (x) )
  #define floorl(x)   (long double) floor( (double) (x) )
  #define ceill(x)    (long double) ceil( (double) (x) )
  #define sqrtl(x)    (long double) sqrt( (double) (x) )
  #define powl(x, y)  (long double) pow( (double) (x), (double) (y) )
  #define expl(x)     (long double) exp( (double) (x) )
  #define logl(x)     (long double) log( (double) (x) )
  #define log1pl(x)   (long double) log1p( (double) (x) )
#endif

#ifdef LDBL_INFINITY
  #undef INFINITY
  #define INFINITY LDBL_INFINITY
#elif !defined(INFINITY)
  #define INFINITY (DBL_MAX + DBL_MAX)
#endif
#ifndef LDBL_EPSILON
  #define LDBL_EPSILON 1e-16
#endif
#ifndef LDBL_MAX
  #define LDBL_MAX DBL_MAX
#endif

#endif
