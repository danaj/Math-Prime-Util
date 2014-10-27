#ifndef MPU_SIEVE_H
#define MPU_SIEVE_H

#include "ptypes.h"

extern unsigned char* sieve_erat30(UV end);
extern int sieve_segment(unsigned char* mem, UV startd, UV endd);
extern void* start_segment_primes(UV low, UV high, unsigned char** segmentmem);
extern int next_segment_primes(void* vctx, UV* base, UV* low, UV* high);
extern void end_segment_primes(void* vctx);


static const UV wheel30[] = {1, 7, 11, 13, 17, 19, 23, 29};
/* Used for moving between primes */
static const unsigned char nextwheel30[30] = {
    1,  7,  7,  7,  7,  7,  7, 11, 11, 11, 11, 13, 13, 17, 17,
   17, 17, 19, 19, 23, 23, 23, 23, 29, 29, 29, 29, 29, 29,  1 };
static const unsigned char prevwheel30[30] = {
   29, 29,  1,  1,  1,  1,  1,  1,  7,  7,  7,  7, 11, 11, 13,
   13, 13, 13, 17, 17, 19, 19, 19, 19, 23, 23, 23, 23, 23, 23 };
/* The bit mask within a byte */
static const unsigned char masktab30[30] = {
    0,  1,  0,  0,  0,  0,  0,  2,  0,  0,  0,  4,  0,  8,  0,
    0,  0, 16,  0, 32,  0,  0,  0, 64,  0,  0,  0,  0,  0,128  };
/* Inverse of masktab30 */
static const unsigned char imask30[129] = {
0,1,7,0,11,0,0,0,13,0,0,0,0,0,0,0,17,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,23,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,29};
/* Add this to a number and you'll ensure you're on a wheel location */
static const unsigned char distancewheel30[30] =
    {1,0,5,4,3,2,1,0,3,2,1,0,1,0,3,2,1,0,1,0,3,2,1,0,5,4,3,2,1,0};
/* add this to n to get to the next wheel location */
static const unsigned char wheeladvance30[30] =
    {1,6,5,4,3,2,1,4,3,2,1,2,1,4,3,2,1,2,1,4,3,2,1,6,5,4,3,2,1,2};
/* subtract this from n to get to the previous wheel location */
static const unsigned char wheelretreat[30] =
    {1,2,1,2,3,4,5,6,1,2,3,4,1,2,1,2,3,4,1,2,1,2,3,4,1,2,3,4,5,6};


#ifdef FUNC_is_prime_in_sieve
static int is_prime_in_sieve(const unsigned char* sieve, UV p) {
  UV d = p/30;
  UV m = p - d*30;
  /* If m isn't part of the wheel, we return 0 */
  return ( (masktab30[m] != 0) && ((sieve[d] & masktab30[m]) == 0) );
}
#endif

#ifdef FUNC_next_prime_in_sieve
/* Will return 0 if it goes past lastp */
static UV next_prime_in_sieve(const unsigned char* sieve, UV p, UV lastp) {
  UV d, m;
  if (p < 7)
    return (p < 2) ? 2 : (p < 3) ? 3 : (p < 5) ? 5 : 7;
  d = p/30;
  m = p - d*30;
  do {
    if (m != 29) {
      m = nextwheel30[m];
    } else {
      d++; m = 1;
      if (d*30 >= lastp) return 0; /* sieves have whole bytes filled */
    }
  } while (sieve[d] & masktab30[m]);
  return d*30+m;
}
#endif
#ifdef FUNC_prev_prime_in_sieve
static UV prev_prime_in_sieve(const unsigned char* sieve, UV p) {
  UV d, m;
  if (p <= 7)
    return (p <= 2) ? 0 : (p <= 3) ? 2 : (p <= 5) ? 3 : 5;
  d = p/30;
  m = p - d*30;
  do {
    m = prevwheel30[m];
    if (m==29) { if (d == 0) return 0;  d--; }
  } while (sieve[d] & masktab30[m]);
  return(d*30+m);
}
#endif

/* Useful macros for the wheel-30 sieve array */
#define START_DO_FOR_EACH_SIEVE_PRIME(sieve, a, b) \
  { \
    const unsigned char* sieve_ = sieve; \
    UV p = a; \
    UV l_ = b; \
    UV d_ = p/30; \
    UV lastd_ = l_/30; \
    UV mask_ = masktab30[ p-d_*30 + distancewheel30[ p-d_*30 ] ]; \
    UV s_ = sieve_[d_]; \
    while (d_ <= lastd_ && (s_ & mask_)) { \
      mask_ <<= 1;  if (mask_ > 128) { s_ = sieve_[++d_]; mask_ = 1; } \
    } \
    p = d_*30 + imask30[mask_]; \
    while ( p <= l_ ) { \

#define END_DO_FOR_EACH_SIEVE_PRIME \
      do { \
        mask_ <<= 1; \
        if (mask_ > 128) { \
          while (++d_ <= lastd_) { s_ = sieve_[d_]; if (s_ != 0xFF) break; } \
          if (d_ > lastd_) break; \
          mask_ = 1; \
        } \
      } while (s_ & mask_); \
      p = d_*30 + imask30[mask_]; \
    } \
  }

#define START_DO_FOR_EACH_PRIME(a, b) \
  { \
    const unsigned char* sieve_; \
    UV p  = a; \
    UV l_ = b; \
    UV d_ = p/30; \
    UV s_, mask_ = 2; \
    UV lastd_ = l_/30; \
    get_prime_cache(l_, &sieve_); \
    s_ = sieve_[d_]; \
    if (p <= 5) { \
      p = (p <= 2) ? 2 : (p <= 3) ? 3 : 5; \
    } else if (p != 7) { \
      mask_ = masktab30[ p-d_*30 + distancewheel30[ p-d_*30 ] ]; \
      while (d_ <= lastd_ && (s_ & mask_)) { \
        mask_ <<= 1;  if (mask_ > 128) { s_ = sieve_[++d_]; mask_ = 1; } \
      } \
      p = d_*30 + imask30[mask_]; \
    } \
    while ( p <= l_ ) {

#define RETURN_FROM_EACH_PRIME(retstmt) \
    do { release_prime_cache(sieve_); retstmt; } while (0)

#define END_DO_FOR_EACH_PRIME \
      if (p < 7) { \
        p += 1 + (p > 2); \
      } else { \
        do { \
          mask_ <<= 1; \
          if (mask_ > 128) { \
            if (++d_ > lastd_) break; \
            s_ = sieve_[d_]; \
            mask_ = 1; \
          } \
        } while (s_ & mask_); \
        p = d_*30 + imask30[mask_]; \
      } \
    } \
    release_prime_cache(sieve_); \
  }

#endif
