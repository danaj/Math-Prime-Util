#include <stdio.h>
#include <stdlib.h>

#define FUNC_is_prime_in_sieve 1
#include "sieve.h"
#include "ptypes.h"
#include "util.h"
#include "primality.h"

typedef struct {
  uint32_t nmax;
  uint32_t nsize;
  UV* list;
} vlist;
#define INIT_VLIST(v) \
  v.nsize = 0; \
  v.nmax = 100; \
  New(0, v.list, v.nmax, UV);
#define PUSH_VLIST(v, n) \
  do { \
    if (v.nsize >= v.nmax) \
      Renew(v.list, v.nmax += 100, UV); \
    v.list[v.nsize++] = n; \
  } while (0)

/* Given p prime, is this a cluster? */
static int is_cluster(UV p, uint32_t nc, uint32_t* cl) {
  uint32_t c;
  for (c = 1; c < nc; c++)
    if (!_XS_is_prime(p+cl[c]))
      break;
  return (c == nc);
}

/* This is fine for small ranges.  Low overhead. */
UV* sieve_cluster_simple(UV beg, UV end, uint32_t nc, uint32_t* cl, UV* numret)
{
  vlist retlist;

  INIT_VLIST(retlist);
  if (beg <= 2 && end >= 2 && is_cluster(2, nc, cl)) PUSH_VLIST(retlist, 2);
  if (beg <= 3 && end >= 3 && is_cluster(3, nc, cl)) PUSH_VLIST(retlist, 3);
  if (beg <= 5 && end >= 5 && is_cluster(5, nc, cl)) PUSH_VLIST(retlist, 5);
  if (beg < 7)  beg = 7;
  if (beg <= end) {
    uint32_t c;
    unsigned char* segment;
    UV seg_base, seg_beg, seg_end;
    
    void* ctx = start_segment_primes(beg, end, &segment);
    while (next_segment_primes(ctx, &seg_base, &seg_beg, &seg_end)) {
      UV sp, last_sieve_cluster = (seg_end >= cl[nc-1]) ? seg_end-cl[nc-1] : 0;
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_beg, seg_end )
        if (p <= last_sieve_cluster) {
          sp = p - seg_base;
          for (c = 1; c < nc; c++)
            if (!is_prime_in_sieve(segment, sp+cl[c]))
              break;
          if (c == nc)
            PUSH_VLIST(retlist,p);
        } else {
          if (is_cluster(p, nc, cl))
            PUSH_VLIST(retlist, p);
        }
      END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
  }
  *numret = retlist.nsize;
  return retlist.list;
}
