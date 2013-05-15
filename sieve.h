#ifndef MPU_SIEVE_H
#define MPU_SIEVE_H

#include "EXTERN.h"
#include "perl.h"

extern unsigned char* sieve_erat30(UV end);
extern int sieve_segment(unsigned char* mem, UV startd, UV endd);


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
/* Add this to a number and you'll ensure you're on a wheel location */
static const unsigned char distancewheel30[30] =
    {1,0,5,4,3,2,1,0,3,2,1,0,1,0,3,2,1,0,1,0,3,2,1,0,5,4,3,2,1,0};
/* Once on the wheel, add this to get to next spot.  In p space, not m. */
static const unsigned char wheeladvance30[30] =
    {0,6,0,0,0,0,0,4,0,0,0,2,0,4,0,0,0,2,0,4,0,0,0,6,0,0,0,0,0,2};

#ifdef FUNC_is_prime_in_sieve
static int is_prime_in_sieve(const unsigned char* sieve, UV p) {
  UV d = p/30;
  UV m = p - d*30;
  /* If m isn't part of the wheel, we return 0 */
  return ( (masktab30[m] != 0) && ((sieve[d] & masktab30[m]) == 0) );
}
#endif

#ifdef FUNC_next_prime_in_sieve
/* Warning -- can go off the end of the sieve */
static UV next_prime_in_sieve(const unsigned char* sieve, UV p) {
  UV d, m;
  if (p < 7)
    return (p < 2) ? 2 : (p < 3) ? 3 : (p < 5) ? 5 : 7;
  d = p/30;
  m = p - d*30;
  do {
    if (m==29) { d++; m = 1; while (sieve[d] == 0xFF) d++; }
    else       { m = nextwheel30[m]; }
  } while (sieve[d] & masktab30[m]);
  return(d*30+m);
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
    m = prevwheel30[m];  if (m==29) { if (d == 0) return 0;  d--; }
  } while (sieve[d] & masktab30[m]);
  return(d*30+m);
}
#endif

/* Useful macros for the wheel-30 sieve array */
#define START_DO_FOR_EACH_SIEVE_PRIME(sieve, a, b) \
  { \
    UV p = a; \
    UV l_ = b; \
    UV d_ = p/30; \
    UV m_ = p-d_*30; \
    m_ += distancewheel30[m_]; \
    p = d_*30 + m_; \
    while ( p <= l_ ) { \
      if ((sieve[d_] & masktab30[m_]) == 0)

#define END_DO_FOR_EACH_SIEVE_PRIME \
      m_ = nextwheel30[m_];  if (m_ == 1) { d_++; } \
      p = d_*30+m_; \
    } \
  }

#define START_DO_FOR_EACH_PRIME(a, b) \
  { \
    const unsigned char* sieve_; \
    UV d_ = 0; \
    UV m_ = 7; \
    UV p  = a; \
    UV l_ = b; \
    if (get_prime_cache(l_, &sieve_) < l_) { \
      release_prime_cache(sieve_); \
      croak("Could not generate sieve for %"UVuf, l_); \
    } \
    if (p <= 5) { \
      p = (p <= 2) ? 2 : (p <= 3) ? 3 : 5; \
    } else { \
      d_ = p/30; \
      m_ = p-d_*30; \
      m_ += distancewheel30[m_]; \
      p = d_*30 + m_; \
    } \
    while ( p <= l_ ) { \
      if (p < 7 || !(sieve_[d_] & masktab30[m_]))

#define RETURN_FROM_EACH_PRIME(x) \
    do { release_prime_cache(sieve_); return(x); } while (0)

#define END_DO_FOR_EACH_PRIME \
      if (p < 7) { \
        p += 1 + (p > 2); \
      } else { \
        m_ = nextwheel30[m_];  if (m_ == 1) { d_++; } \
        p = d_*30+m_; \
      } \
    } \
    release_prime_cache(sieve_); \
  }

#endif
