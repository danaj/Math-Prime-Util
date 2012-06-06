#ifndef MPU_BITARRAY_H
#define MPU_BITARRAY_H

#include "ptypes.h"

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
