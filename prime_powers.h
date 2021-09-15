#ifndef MPU_PRIME_POWERS_H
#define MPU_PRIME_POWERS_H

#include "ptypes.h"

extern int prime_power(UV n, UV* prime);
static INLINE int is_prime_power(UV n) { return !!prime_power(n,0); }

extern UV  next_prime_power(UV x);
extern UV  prev_prime_power(UV x);

extern UV prime_power_sieve2(UV** list, UV lo, UV hi);  /* p^e with e >= 2 */
extern UV prime_power_sieve(UV** list, UV lo, UV hi);   /* p^e with e >= 1 */

extern UV  prime_power_count_range(UV lo, UV hi);

extern UV  prime_power_count(UV n);
extern UV  prime_power_count_lower(UV n);
extern UV  prime_power_count_upper(UV n);
extern UV  prime_power_count_approx(UV n);

extern UV  nth_prime_power(UV n);
extern UV  nth_prime_power_lower(UV n);
extern UV  nth_prime_power_upper(UV n);
extern UV  nth_prime_power_approx(UV n);

#endif
