#ifndef MPU_PTYPES_H
#define MPU_PTYPES_H

#ifdef _MSC_VER
 /* No stdint.h for MS C, but all the types can be defined.
  *
  * Thanks to Sisyphus and bulk88 for all the help with MSC,
  * including working patches.
  */
typedef unsigned __int8  uint8_t;
typedef unsigned __int16 uint16_t;
typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64_t;
typedef __int64 int64_t;
typedef __int32 int32_t;
typedef __int16 int16_t;
typedef __int8 int8_t;

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

#ifndef PERL_VERSION
  #include <limits.h>
  #include <stdio.h>
  #include <stdlib.h>
  #include <ctype.h>
  typedef unsigned long UV;
  typedef   signed long IV;
  typedef size_t        STRLEN;
  #define UV_MAX ULONG_MAX
  #define IV_MAX LONG_MAX
  #define UVCONST(x) ((unsigned long)x##UL)
  #define U32_CONST(x) ((unsigned int)x##U)
  #define UVuf "lu"
  #define IVdf "ld"
  #define croak(fmt,...)            { printf(fmt,##__VA_ARGS__); exit(3); }
  #define New(id, mem, size, type)  mem = (type*) malloc((size)*sizeof(type))
  #define Newz(id, mem, size, type) mem = (type*) calloc(size, sizeof(type))
  #define Renew(mem, size, type)    mem =(type*)realloc(mem,(size)*sizeof(type))
  #define Safefree(mem)             free((void*)mem)
  #define isDIGIT(x)                isdigit(x)
#if ULONG_MAX >> 31 == 1
  #define BITS_PER_WORD  32
#elif ULONG_MAX >> 63 == 1
  #define BITS_PER_WORD  64
#else
  #error Unsupported bits per word (must be 32 or 64)
#endif

#else

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


/* See:
 *   http://www.nntp.perl.org/group/perl.perl5.porters/2013/09/msg207524.html
 * for some discussion.
 */
#ifdef HAS_QUAD
  #define BITS_PER_WORD  64
  #define UVCONST(x)     U64_CONST(x)
#else
  #define BITS_PER_WORD  32
  #define UVCONST(x)     U32_CONST(x)
#endif

#endif   /* End of Perl specific section */

/* Try to determine if we have 64-bit available via uint64_t */
#if defined(UINT64_MAX) || defined(_UINT64_T) || defined(__UINT64_TYPE__)
  #define HAVE_STD_U64 1
#elif defined(_MSC_VER)   /* We set up the types earlier */
  #define HAVE_STD_U64 1
#else
  #define HAVE_STD_U64 0
#endif

#define MAXBIT        (BITS_PER_WORD-1)
#define NWORDS(bits)  ( ((bits)+BITS_PER_WORD-1) / BITS_PER_WORD )
#define NBYTES(bits)  ( ((bits)+8-1) / 8 )

#define MPUassert(c,text) if (!(c)) { croak("Math::Prime::Util internal error: " text); }

/* The ASSUME bits are from perl 5.19.6 perl.h */

#ifndef __has_builtin
#  define __has_builtin(x) 0 /* not a clang style compiler */
#endif

#ifndef DEBUGGING
#  if (__GNUC__ == 4 && __GNUC_MINOR__ >= 5 || __GNUC__ > 5) || __has_builtin(__builtin_unreachable)
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

#if (defined(__SIZEOF_INT128__) || __GNUC__ >= 4) && (defined(__x86_64__) || defined(__powerpc64__))
  #if __clang__ && (__clang_major__ > 4 || (__clang_major__ == 4 && __clang_minor__ >= 2))
    #define HAVE_UINT128 1
    typedef unsigned __int128 uint128_t;
  #elif __GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 4)
    #define HAVE_UINT128 0
  #elif __GNUC__ == 4 && __GNUC_MINOR__ >= 4 && __GNUC_MINOR__ < 6
    #define HAVE_UINT128 1
    typedef unsigned int uint128_t __attribute__ ((__mode__ (TI)));
  #else
    #define HAVE_UINT128 1
    typedef unsigned __int128 uint128_t;
  #endif
#else
  #define HAVE_UINT128 0
#endif

#if defined(__GNUC__) || defined(__clang__)
  #define INLINE inline
#elif defined(_MSC_VER)
  #define INLINE __inline
#else
  #define INLINE
#endif

#endif
