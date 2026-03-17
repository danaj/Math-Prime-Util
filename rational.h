#ifndef MPU_RATIONAL_H
#define MPU_RATIONAL_H

#include "ptypes.h"

/* Returns number of steps (max 2 * BITS_PER_WORD) */
extern int contfrac(UV** cfrac, UV *rem, UV num, UV den);
/* Returns false if overflow */
extern bool convergents(UV** P, UV** Q, UV* cfrac, size_t len);

extern bool next_calkin_wilf(UV* num, UV* den);
extern bool next_stern_brocot(UV* num, UV* den);

extern UV calkin_wilf_n(UV num, UV den);
extern UV stern_brocot_n(UV num, UV den);

extern bool nth_calkin_wilf(UV* num, UV* den, UV n);
extern bool nth_stern_brocot(UV* num, UV* den, UV n);

extern UV nth_stern_diatomic(UV n);

extern UV farey_length(uint32_t n);
extern bool next_farey(uint32_t n, uint32_t* p, uint32_t* q);
extern bool kth_farey(uint32_t n, UV k, uint32_t* p, uint32_t* q);
extern UV farey_array(uint32_t n, uint32_t **num, uint32_t **den);

/* How many fractions are < p/q in F_n */
extern UV farey_rank(uint32_t n, uint32_t p, uint32_t q);


/* best rational approximation to x that has denominator <= dbound */
extern bool bestrational(UV* n, UV* d, NV x, UV dbound);

#endif
