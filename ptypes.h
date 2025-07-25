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

#ifdef STANDALONE
  #include <limits.h>
  #include <stdio.h>
  #include <stdlib.h>
  #include <ctype.h>
  #include <stdbool.h>
  #define TRUE true
  #define FALSE false
  typedef unsigned long UV;
  typedef   signed long IV;
  typedef        double NV;
  typedef        size_t STRLEN;
  #define UV_MAX ULONG_MAX
  #define IV_MAX LONG_MAX
  #define NV_MAX DBL_MAX
  #define UVCONST(x) ((unsigned long)x##UL)
  #define U32_CONST(x) ((unsigned int)x##U)
  #define UVuf "lu"
  #define IVdf "ld"
  #define NVff "f"
  /* Technically this is sizeof(NV) but that's not valid for macros */
  #define NVSIZE 8
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
  #define HAVE_UINT64 1
#elif defined(_MSC_VER)   /* We set up the types earlier */
  #define HAVE_UINT64 1
#else
  #define HAVE_UINT64 0
#endif

#define MAXBIT        (BITS_PER_WORD-1)
#define NWORDS(bits)  ( ((bits)+BITS_PER_WORD-1) / BITS_PER_WORD )
#define NBYTES(bits)  ( ((bits)+8-1) / 8 )

#define MPUassert(c,text) if (!(c)) { croak("Math::Prime::Util internal error: " text); }
#define MPUverbose(level,fmt,...) \
  if (_XS_get_verbose() >= level) { printf(fmt,##__VA_ARGS__); fflush(stdout); }

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

#if defined(__SIZEOF_INT128__) && !defined(__CUDACC__)
  #define HAVE_UINT128 1
  typedef unsigned __int128 uint128_t;
#elif (__GNUC__ >= 4) && (defined(__x86_64__) || defined(__powerpc64__))
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
#elif defined(__BITINT_MAXWIDTH__) && __BITINT_MAXWIDTH__ >= 128
  /* Should have included <stdint.h> and <limits.h> already */
  #define HAVE_UINT128 1
  typedef unsigned _BitInt(128) uint128_t;
#else
  #define HAVE_UINT128 0
#endif

/* Perl 5.23.0 added the very helpful definition.  Without it, guess. */
#ifndef NVMANTBITS
  #if NVSIZE <= 8
    #define NVMANTBITS ((NVSIZE <= 4) ? 24 : 53)
  #elif defined(USE_QUADMATH)
    #define NVMANTBITS 112
  #elif defined(__LDBL_MANT_DIG__)
    #define NVMANTBITS __LDBL_MANT_DIG__
  #else
    #define NVMANTBITS 64
  #endif
#endif

#if defined(USE_QUADMATH)
  typedef __float128   LNV;
  #define LNV_ZERO     0.0Q
  #define LNV_ONE      1.0Q
  #define LNVCONST(x) ((__float128)x##Q)
  #define loglnv(x)    logq(x)
  #define explnv(x)    expq(x)
  #define sqrtlnv(x)   sqrtq(x)
  #define fabslnv(x)   fabsq(x)
  #define LNV_EPSILON  FLT128_EPSILON
  #define LNV_IS_QUAD  1
#else
  typedef long double  LNV;
  #define LNV_ZERO     0.0L
  #define LNV_ONE      1.0L
  #define LNVCONST(x) ((long double)x##L)
  #define loglnv(x)    logl(x)
  #define explnv(x)    expl(x)
  #define sqrtlnv(x)   sqrtl(x)
  #define fabslnv(x)   fabsl(x)
  #define LNV_EPSILON  LDBL_EPSILON
  #define LNV_IS_QUAD  0
#endif

#if (defined(__GNUC__) || defined(__clang__)) && __STDC_VERSION__ >= 199901L
  #define INLINE inline
#elif defined(_MSC_VER)
  #define INLINE __inline
#else
  #define INLINE
#endif

#if (defined(__GNUC__) || defined(__clang__)) && !defined(__INTEL_COMPILER)
  #define ISCONSTFUNC __attribute__((const))
  #define NOINLINE __attribute__((noinline))
#else
  #define ISCONSTFUNC
  #define NOINLINE
#endif

#if __BIG_ENDIAN__ || (defined(BYTEORDER) && (BYTEORDER == 0x4321 || BYTEORDER == 0x87654321))
#  if (__GNUC__ == 4 && __GNUC_MINOR__ >= 4) || __GNUC__ >= 5 || (__clang__ &&  __clang_major__ >= 4)
#    if BITS_PER_WORD == 64
#      define LEUV(x) __builtin_bswap64(x)
#    else
#      define LEUV(x) __builtin_bswap32(x)
#    endif
#  else
#    if BITS_PER_WORD == 64
       /* compare to 5 step interleave */
       static UV LEUV(UV x) {
         UV v = ((x & UVCONST(0xFF00000000000000)) >> 56) |
                ((x & UVCONST(0x00FF000000000000)) >> 40) |
                ((x & UVCONST(0x0000FF0000000000)) >> 24) |
                ((x & UVCONST(0x000000FF00000000)) >>  8) |
                ((x & UVCONST(0x00000000FF000000)) <<  8) |
                ((x & UVCONST(0x0000000000FF0000)) << 24) |
                ((x & UVCONST(0x000000000000FF00)) << 40) |
                ((x & UVCONST(0x00000000000000FF)) << 56);
         return v;
       }
#    else
       static UV LEUV(UV x) {
         UV v = ((x & 0xFF000000) >> 24) |
                ((x & 0x00FF0000) >>  8) |
                ((x & 0x0000FF00) <<  8) |
                ((x & 0x000000FF) << 24);
         return v;
       }
#    endif
#  endif
#else /* LE */
#  define LEUV(x) (x)
#endif

#endif
