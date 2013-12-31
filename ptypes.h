#ifndef MPU_PTYPES_H
#define MPU_PTYPES_H

#ifdef _MSC_VER
 /* No stdint.h for MS C, so we lose the chance to possibly optimize
  * some operations on 64-bit machines running a 32-bit Perl.  It's probably
  * a rare enough case that we don't need to be too concerned.  If we do want,
  * see:  http://gauss.cs.ucsb.edu/~aydin/CombBLAS/html/stdint_8h_source.html
  * for some ideas.
  *
  *  Thanks to Sisyphus for bringing the MSC issue to my attention (and even
  *  submitting a working patch!).
  */
typedef unsigned __int8  uint8_t;
typedef unsigned __int16 uint16_t;
typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64_t;
typedef __int64 int64_t;
typedef __int32 int32_t;
typedef __int16 int16_t;
typedef __int8 int8_t;
#define inline __inline

#ifdef _M_X64
#  define __x86_64__
#  define __x86_64
#  define __amd64__
#  define __amd64
#endif
#ifdef _M_IX86
#  define __i386__
#  define __i386
#  define i386
#  define _X86_
#endif
#ifdef _M_IA64
#  define __ia64__
#  define __ia64
#  define __IA64__
#  define __itanium__
#endif

#elif defined(__sun) || defined(__sun__)
 /* stdint.h is only in Solaris 10+. */
 #if defined(__SunOS_5_10) || defined(__SunOS_5_11) || defined(__SunOS_5_12)
  #define __STDC_LIMIT_MACROS
  #include <stdint.h>
 #endif
#else
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#endif

#include "EXTERN.h"
#include "perl.h"

/* From perl.h, wrapped in PERL_CORE */
#ifndef U32_CONST
# if INTSIZE >= 4
#  define U32_CONST(x) ((U32TYPE)x##U)
# else
#  define U32_CONST(x) ((U32TYPE)x##UL)
# endif
#endif

/* From perl.h, wrapped in PERL_CORE */
#ifndef U64_CONST
# ifdef HAS_QUAD
#  if INTSIZE >= 8
#   define U64_CONST(x) ((U64TYPE)x##U)
#  elif LONGSIZE >= 8
#   define U64_CONST(x) ((U64TYPE)x##UL)
#  elif QUADKIND == QUAD_IS_LONG_LONG
#   define U64_CONST(x) ((U64TYPE)x##ULL)
#  else /* best guess we can make */
#   define U64_CONST(x) ((U64TYPE)x##UL)
#  endif
# endif
#endif


#ifdef HAS_QUAD
  #define BITS_PER_WORD  64
  #define UVCONST(x)     U64_CONST(x)
#else
  #define BITS_PER_WORD  32
  #define UVCONST(x)     U32_CONST(x)
#endif

/* Try to determine if we have 64-bit available via uint64_t */
#define HAVE_STD_U64 0
#if defined(UINT64_MAX) && defined(__UINT64_C)
  #if (UINT64_MAX >= __UINT64_C(18446744073709551615))
    #undef HAVE_STD_U64
    #define HAVE_STD_U64 1
  #endif
#endif

#define MAXBIT        (BITS_PER_WORD-1)
#define NWORDS(bits)  ( ((bits)+BITS_PER_WORD-1) / BITS_PER_WORD )
#define NBYTES(bits)  ( ((bits)+8-1) / 8 )

#define MPUassert(c,text) if (!(c)) { croak("Math::Prime::Util internal error: " text); }

#ifndef DEBUGGING
#  if __has_builtin(__builtin_unreachable) \
     || (__GNUC__ == 4 && __GNUC_MINOR__ >= 5 || __GNUC__ > 5) /* 4.5 -> */
#    define MPUASSUME(x) ((x) ? (void) 0 : __builtin_unreachable())
#  elif defined(_MSC_VER)
#    define MPUASSUME(x) __assume(x)
#  elif defined(__ARMCC_VERSION) /* untested */
#    define MPUASSUME(x) __promise(x)
#  else
/* a random compiler might define assert to its own special optimization token
   so pass it through to C lib as a last resort */
#    define MPUASSUME(x) assert(x)
#  endif
#else
#  define MPUASSUME(x) assert(x)
#endif

#define MPUNOT_REACHED MPUASSUME(0)

#endif
