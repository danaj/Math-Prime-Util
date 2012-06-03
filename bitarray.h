#ifndef MPU_BITARRAY_H
#define MPU_BITARRAY_H

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

#define MAXBIT        (BITS_PER_WORD-1)
#define NWORDS(bits)  ( ((bits)+BITS_PER_WORD-1) / BITS_PER_WORD )
#define NBYTES(bits)  ( ((bits)+8-1) / 8 )

/* The mod-30 wheel is byte-based so doesn't use these. */
#define SET_ARRAY_BIT(ar,n) \
   ar[(n)/BITS_PER_WORD]  |=  (UVCONST(1) << ((n)%BITS_PER_WORD))
#define XOR_ARRAY_BIT(ar,n) \
   ar[(n)/BITS_PER_WORD]  ^=  (UVCONST(1) << ((n)%BITS_PER_WORD))
#define CLR_ARRAY_BIT(ar,n) \
   ar[(n)/BITS_PER_WORD]  &=  ~(UVCONST(1) << ((n)%BITS_PER_WORD))
#define IS_SET_ARRAY_BIT(ar,n) \
   (ar[(n)/BITS_PER_WORD] & (UVCONST(1) << ((n)%BITS_PER_WORD)) )

#endif
