#ifndef MPU_FACTOR_H
#define MPU_FACTOR_H

#include "ptypes.h"

#define MPU_MAX_FACTORS 64

extern int trial_factor(UV n, UV *factors, UV maxtrial);

extern int fermat_factor(UV n, UV *factors, UV rounds);

extern int squfof_factor(UV n, UV *factors, UV rounds);

extern int pbrent_factor(UV n, UV *factors, UV maxrounds);

extern int prho_factor(UV n, UV *factors, UV maxrounds);

extern int pminus1_factor(UV n, UV *factors, UV maxrounds);

extern int miller_rabin(UV n, const UV *bases, UV nbases);
extern int is_prob_prime(UV n);

#if 0
/* Try to make a quick probable prime test */
static int is_perhaps_prime(UV n)
{
  static const UV bases2[2] = {31, 73};
  static const UV bases3[3] = {2,7,61};
  if ( (n == 2) || (n == 3) || (n == 5) || (n == 7) )
    return 2;
  if ( (n<2) || ((n% 2)==0) || ((n% 3)==0) || ((n% 5)==0) || ((n% 7)==0) )
    return 0;
  if (n < 121)
    return 2;
  if (n < UVCONST(9080191))
    return 2*miller_rabin(n, bases2, 2);
  else
    return miller_rabin(n, bases3, 3) * ((n < UVCONST(4759123141)) ? 2 : 1);
}

static int is_perhaps_prime7(UV n)
{
  static const UV bases2[2] = {31, 73};
  static const UV bases3[3] = {2,7,61};
  /* n must be >= 73 */
  if (n < UVCONST(9080191))
    return 2*miller_rabin(n, bases2, 2);
  else
    return miller_rabin(n, bases3, 3) * ((n < UVCONST(4759123141)) ? 2 : 1);
}
#endif

#endif
