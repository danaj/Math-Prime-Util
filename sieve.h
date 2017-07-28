#ifndef MPU_SIEVE_H
#define MPU_SIEVE_H

#include "ptypes.h"
#define FUNC_ctz 1
#include "util.h"

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
static const unsigned char wheelretreat30[30] =
    {1,2,1,2,3,4,5,6,1,2,3,4,1,2,1,2,3,4,1,2,1,2,3,4,1,2,3,4,5,6};
/* Given a sieve byte, this indicates the first zero */
static const unsigned char nextzero30[256] =
  {0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,
   0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,6,0,1,0,2,0,1,0,3,0,1,0,2,
   0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,
   0,2,0,1,0,3,0,1,0,2,0,1,0,7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,
   0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,
   0,6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,
   0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,8};
/* At this m (p-30*(p/30)), OR with this to clear previous entries */
static const unsigned char clearprev30[30] =
    {  0,  0,  1,  1,  1,  1,  1,  1,  3,  3,  3,  3,  7,  7, 15,
      15, 15, 15, 31, 31, 63, 63, 63, 63,127,127,127,127,127,127};


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
  unsigned char s;
  if (p < 7)
    return (p < 2) ? 2 : (p < 3) ? 3 : (p < 5) ? 5 : 7;
  p++;
  if (p >= lastp) return 0;
  d = p/30;
  m = p - d*30;
  s = sieve[d] | clearprev30[m];
  while (s == 0xFF) {
    d++;
    if (d*30 >= lastp) return 0;
    s = sieve[d];
  }
  return d*30 + wheel30[nextzero30[s]];
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

#if 0
/* Useful macros for the wheel-30 sieve array */
#define START_DO_FOR_EACH_SIEVE_PRIME(sieve, base, a, b) \
  { \
    const unsigned char* sieve_ = sieve; \
    UV base_ = base; \
    UV p = a-base_; \
    UV l_ = b; \
    UV d_ = p/30; \
    UV lastd_ = (l_-base_)/30; \
    unsigned char bit_, s_ = sieve_[d_] | clearprev30[p-d_*30]; \
    base_ += d_*30; \
    while (1) { \
      if (s_ == 0xFF) { \
        do { \
          base_ += 30;  d_++; \
          if (d_ > lastd_) break; \
          s_ = sieve_[d_]; \
        } while (s_ == 0xFF); \
        if (d_ > lastd_) break; \
      } \
      bit_ = nextzero30[s_]; \
      s_ |= 1 << bit_; \
      p = base_ + wheel30[bit_]; \
      if (p > l_ || p < base_) break; /* handle overflow */ \
      {

#define END_DO_FOR_EACH_SIEVE_PRIME \
      } \
    } \
  }
#else
/* Extract word at a time, good suggestion from Kim Walisch */
static const unsigned char wheel240[] = {1,7,11,13,17,19,23,29,31,37,41,43,47,49,53,59,61,67,71,73,77,79,83,89,91,97,101,103,107,109,113,119,121,127,131,133,137,139,143,149,151,157,161,163,167,169,173,179,181,187,191,193,197,199,203,209,211,217,221,223,227,229,233,239};
#define START_DO_FOR_EACH_SIEVE_PRIME(sieve, base, a, b) \
  { \
    const UV* sieve_ = (const UV*)sieve;        /* word ptr to sieve */ \
    const UV nperw_ = 30*sizeof(UV); /* nums per word */     \
    UV base_ = base;                 /* start of sieve n */  \
    UV b_ = a;                       /* begin value n */     \
    UV f_ = b;                       /* final value n */     \
    UV begw_ = (b_-base_)/nperw_;    /* first word */        \
    UV endw_ = (f_-base_)/nperw_;    /* first word */        \
    UV sw_, tz_, p; \
    base_ += begw_*nperw_; \
    while (begw_ <= endw_) { \
      sw_ = ~ LEUV(sieve_[begw_]); \
      while (sw_ != 0) { \
        tz_ = ctz(sw_); \
        sw_ &= ~(UVCONST(1) << tz_); \
        p = base_ + wheel240[tz_]; \
        if (p > f_) break; \
        if (p >= b_) {

#define END_DO_FOR_EACH_SIEVE_PRIME \
        } \
      } \
      begw_++; \
      base_ += nperw_; \
    } \
  }
#endif

#define START_DO_FOR_EACH_PRIME(a, b) \
  { \
    const unsigned char* sieve_; \
    UV p  = a; \
    UV l_ = b; \
    UV d_ = p/30; \
    UV lastd_ = l_/30; \
    unsigned char s_, bit_; \
    get_prime_cache(l_, &sieve_); \
    if (p == 2) p = 1; \
    s_ = sieve_[d_] | clearprev30[p-d_*30]; \
    while (1) { \
      if (p < 5) { \
        p = (p < 2) ? 2 : (p < 3) ? 3 : 5; \
      } else { \
        if (s_ == 0xFF) { \
          do { \
            d_++; \
            if (d_ > lastd_) break; \
            s_ = sieve_[d_]; \
          } while (s_ == 0xFF); \
          if (d_ > lastd_) break; \
        } \
        bit_ = nextzero30[s_]; \
        s_ |= 1 << bit_; \
        p = d_*30 + wheel30[bit_]; \
        if (p < d_*30) break; \
      } \
      if (p > l_) break; \
      { \

#define RETURN_FROM_EACH_PRIME(retstmt) \
        do { release_prime_cache(sieve_); retstmt; } while (0)

#define END_DO_FOR_EACH_PRIME \
      } \
    } \
    release_prime_cache(sieve_); \
  }

#endif
