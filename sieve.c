#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sieve.h"
#include "ptypes.h"
#include "bitarray.h"



UV* sieve_erat(UV end)
{
  UV* mem;
  UV n, s;
  UV last = (end+1)/2;

  Newz(0, mem, NWORDS(last), UV );
  if (mem == 0) {
    croak("allocation failure in sieve_erat: could not alloc %"UVuf" bits", last);
    return 0;
  }

  n = 3;
  /* TODO: overflow */
  while ( (n*n) <= end) {
    for (s = n*n; s <= end; s += 2*n)
      SET_ARRAY_BIT(mem,s/2);
    do { n += 2; } while (IS_SET_ARRAY_BIT(mem,n/2));
  }

  SET_ARRAY_BIT(mem, 1/2);  /* 1 is composite */

  return mem;
}


/* Wheel 30 sieve.  Ideas from Terje Mathisen and Quesada / Van Pelt. */
unsigned char* sieve_erat30(UV end)
{
  unsigned char* mem;
  UV max_buf, buffer_words, limit;
  UV prime;

  max_buf = (end/30) + ((end%30) != 0);
  buffer_words = (max_buf + sizeof(UV) - 1) / sizeof(UV);
  Newz(0, mem, buffer_words*sizeof(UV), unsigned char );
  if (mem == 0) {
    croak("allocation failure in sieve_erat30: could not alloc %"UVuf" bytes", (buffer_words*sizeof(UV)));
    return 0;
  }

  /* Shortcut to mark 7.  Just an optimization. */
  if ( (7*7) <= end ) {
    UV d = 1;
    while ( (d+6) < max_buf) {
      mem[d+0] = 0x20;  mem[d+1] = 0x10;  mem[d+2] = 0x81;  mem[d+3] = 0x08;
      mem[d+4] = 0x04;  mem[d+5] = 0x40;  mem[d+6] = 0x02;  d += 7;
    }
    if ( d < max_buf )  mem[d++] = 0x20;
    if ( d < max_buf )  mem[d++] = 0x10;
    if ( d < max_buf )  mem[d++] = 0x81;
    if ( d < max_buf )  mem[d++] = 0x08;
    if ( d < max_buf )  mem[d++] = 0x04;
    if ( d < max_buf )  mem[d++] = 0x40;
    MPUassert( d >= max_buf, "sieve_erat30 incorrect 7 sieve");
  }
  limit = sqrt((double) end);  /* prime*prime can overflow */
  for (prime = 11; prime <= limit; prime = next_prime_in_sieve(mem,prime)) {
    UV d = (prime*prime)/30;
    UV m = (prime*prime) - d*30;
    UV dinc = (2*prime)/30;
    UV minc = (2*prime) - dinc*30;
    UV wdinc[8];
    unsigned char wmask[8];
    int i;

    /* Find the positions of the next composites we will mark */
    for (i = 1; i <= 8; i++) {
      UV dlast = d;
      do {
        d += dinc;
        m += minc;
        if (m >= 30) { d++; m -= 30; }
      } while ( masktab30[m] == 0 );
      wdinc[i-1] = d - dlast;
      wmask[i%8] = masktab30[m];
    }
    d -= prime;
#if 0
    assert(d == ((prime*prime)/30));
    assert(d < max_buf);
    assert(prime = (wdinc[0]+wdinc[1]+wdinc[2]+wdinc[3]+wdinc[4]+wdinc[5]+wdinc[6]+wdinc[7]));
#endif
    i = 0;        /* Mark the composites */
    do {
      mem[d] |= wmask[i];
      d += wdinc[i];
      i = (i+1) & 7;
    } while (d < max_buf);
  }

  mem[0] |= masktab30[1];  /* 1 is composite */

  return mem;
}



int sieve_segment(unsigned char* mem, UV startd, UV endd)
{
  const unsigned char* sieve;
  UV limit;
  UV pcsize;
  UV startp = 30*startd;
  UV endp = (endd >= (UV_MAX/30))  ?  UV_MAX-2  :  30*endd+29;
  UV ranged = endd - startd + 1;

  MPUassert( (mem != 0) && (endd >= startd) && (endp >= startp),
             "sieve_segment bad arguments");
  memset(mem, 0, ranged);

  limit = sqrt( (double) endp ) + 1;
  /* printf("segment sieve from %"UVuf" to %"UVuf" (aux sieve to %"UVuf")\n", startp, endp, limit); */
  pcsize = get_prime_cache(limit, &sieve);
  if (pcsize < limit) {
    croak("Couldn't generate small sieve for segment sieve");
    return 0;
  }

  START_DO_FOR_EACH_SIEVE_PRIME(sieve, 7, pcsize)
  {
    /* p increments from 7 to at least sqrt(endp) */
    UV p2 = p*p;   /* TODO: overflow */
    if (p2 > endp)  break;
    /* Find first multiple of p greater than p*p and larger than startp */
    if (p2 < startp) {
      p2 = (startp / p) * p;
      if (p2 < startp)  p2 += p;
    }
    /* Bump to next multiple that isn't divisible by 2, 3, or 5 */
    while (masktab30[p2%30] == 0) { p2 += p; }
    /* It is possible we've overflowed p2, so check for that */
    if ( (p2 <= endp) && (p2 >= startp) ) {
      /* Sieve from startd to endd starting at p2, stepping p */
#if 0
      /* Basic sieve */
      do {
        mem[(p2 - startp)/30] |= masktab30[p2%30];
        do { p2 += 2*p; } while (masktab30[p2%30] == 0);
      } while ( (p2 <= endp) && (p2 >= startp) );
#else
      UV d = (p2)/30;
      UV m = (p2) - d*30;
      UV dinc = (2*p)/30;
      UV minc = (2*p) - dinc*30;
      UV wdinc[8];
      unsigned char wmask[8];
      int i;

      /* Find the positions of the next composites we will mark */
      for (i = 1; i <= 8; i++) {
        UV dlast = d;
        do {
          d += dinc;
          m += minc;
          if (m >= 30) { d++; m -= 30; }
        } while ( masktab30[m] == 0 );
        wdinc[i-1] = d - dlast;
        wmask[i%8] = masktab30[m];
      }
      d -= p;
      i = 0;        /* Mark the composites */
      do {
        mem[d-startd] |= wmask[i];
        d += wdinc[i];
        i = (i+1) & 7;
      } while (d <= endd);
#endif
    }
  }
  END_DO_FOR_EACH_SIEVE_PRIME;

  if (startd == 0)
    mem[0] |= masktab30[1];  /* 1 is composite */

  return 1;
}
