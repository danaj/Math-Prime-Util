#include <stdio.h>
#include <stdlib.h>

#define FUNC_is_prime_in_sieve 1
#define FUNC_gcd_ui 1
#include "sieve.h"
#include "ptypes.h"
#include "util.h"
#include "primality.h"

#define NSMALLPRIMES 168
static const unsigned short sprimes[NSMALLPRIMES] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997};

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

#define ADDVAL32(v, n, max, val) \
  do { if (n >= max) Renew(v, max += 512, UV);  v[n++] = val; } while (0)
#define SWAPL32(l1, n1, m1,  l2, n2, m2) \
  { UV t_, *u_ = l1;  l1 = l2;  l2 = u_; \
            t_ = n1;  n1 = n2;  n2 = t_; \
            t_ = m1;  m1 = m2;  m2 = t_; }

static int is_admissible(uint32_t nc, uint32_t* cl) {
  uint32_t i, j, c;
  char rset[sprimes[NSMALLPRIMES-1]];

  if (nc > NSMALLPRIMES) return 1;    /* TODO */
  for (i = 0; i < nc; i++) {
    uint32_t p = sprimes[i];
    memset(rset, 0, p);
    for (c = 0; c < nc; c++)
      rset[cl[c] % p] = 1;
    for (j = 0; j < p; j++)
      if (rset[j] == 0)
        break;
    if (j == p) /* All values were 1 */
      return 0;
  }
  return 1;
}

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

  /* If not admissible, then don't keep looking. */
  if (!is_admissible(nc, cl) && end > sprimes[nc])
    end = sprimes[nc];

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


#define addmodded(r,a,b,n)  do { r = a + b; if (r >= n) r -= n; } while(0)

UV* sieve_cluster(UV low, UV high, uint32_t nc, uint32_t* cl, UV* numret)
{
  vlist retlist;
  UV i, ppr, nres, allocres;
  uint32_t const targres = 100000;
  UV *residues, *cres, num_mr = 0, num_lucas = 0;
  uint32_t pp_0, pp_1, pp_2, *resmod_0, *resmod_1, *resmod_2;
  uint32_t rem_0, rem_1, rem_2, remadd_0, remadd_1, remadd_2;
  uint32_t pi, startpi = 1, maxpi = 150;
  uint32_t lastspr = sprimes[maxpi-1];
  uint32_t c, smallnc;
  char crem_0[43*47], crem_1[53*59], crem_2[61*67], **VPrem;
  int const _verbose = _XS_get_verbose();

  if ((UV_MAX - cl[nc-1]) < high)  return 0;  /* Overflow */

  if (    ((high-low) < 10000)
       || (nc == 3 && ((high >> 32) >> 15) == 0) /* sieving large vals is slow */
       || (nc == 2 && ((high >> 32) >> 26) == 0)
       || (nc < 2) )
    return sieve_cluster_simple(low, high, nc, cl, numret);

  if (!(low&1))    low++;
  if (!(high&1))   high--;

  INIT_VLIST(retlist);

  if (low < lastspr) {
    UV t, chigh = (high > lastspr) ? lastspr : high;
    UV* s = sieve_cluster_simple(low, chigh, nc, cl, &t);
    for (i = 0; i < t; i++)
      PUSH_VLIST(retlist, s[i]);
    Safefree(s);
    low = chigh + 2;
  }
  if (low > high) {
    *numret = retlist.nsize;
    return retlist.list;
  }
  if (low&1) low--;

  /* Determine the primorial size and acceptable residues */
  New(0, residues, allocres = 1024, UV);
  {
    UV remr, *res2, allocres2, nres2, maxppr;
    /* Calculate residues for a small primorial */
    for (pi = 2, ppr = 1, i = 0;  i <= pi;  i++) ppr *= sprimes[i];
    remr = low % ppr;
    nres = 0;
    for (i = 1; i <= ppr; i += 2) {
      for (c = 0; c < nc; c++) {
        UV v = (remr + i + cl[c]) % ppr;
        if (gcd_ui(v, ppr) != 1) break;
      }
      if (c == nc)
        ADDVAL32(residues, nres, allocres, i);
    }
    /* Raise primorial size until we have plenty of residues */
    New(0, res2, allocres2 = 1024, UV);
    maxppr = high - low;
#if BITS_PER_WORD == 64
    while (pi++ < 12) {
#else
    while (pi++ < 8) {
#endif
      uint32_t j, p = sprimes[pi];
      UV r, newppr = ppr * p;
      if (nres == 0 || nres > targres/(p/2) || newppr > maxppr) break;
      if (_verbose > 1) printf("cluster sieve found %lu residues mod %lu\n", nres, ppr);
      remr = low % newppr;
      nres2 = 0;
      for (i = 0; i < p; i++) {
        for (j = 0; j < nres; j++) {
          r = i*ppr + residues[j];
          for (c = 0; c < nc; c++) {
            UV v = remr + r + cl[c];
            if ((v % p) == 0) break;
          }
          if (c == nc)
            ADDVAL32(res2, nres2, allocres2, r);
        }
      }
      ppr = newppr;
      SWAPL32(residues, nres, allocres,  res2, nres2, allocres2);
    }
    startpi = pi;
    Safefree(res2);
  }
  if (_verbose) printf("cluster sieve using %lu residues mod %lu\n", nres, ppr);

  /* Return if not admissible, maybe with a single small value */
  if (nres == 0) {
    Safefree(residues);
    *numret = retlist.nsize;
    return retlist.list;
  }

  /* Pre-mod the residues with first two primes for fewer modulos every chunk */
  {
    uint32_t p1 = sprimes[startpi+0], p2 = sprimes[startpi+1];
    uint32_t p3 = sprimes[startpi+2], p4 = sprimes[startpi+3];
    uint32_t p5 = sprimes[startpi+4], p6 = sprimes[startpi+5];
    pp_0 = p1*p2; pp_1 = p3*p4; pp_2 = p5*p6;
    memset(crem_0, 1, pp_0);
    memset(crem_1, 1, pp_1);
    memset(crem_2, 1, pp_2);
    /* Mark remainders that indicate a composite for this residue. */
    for (i = 0; i < p1; i++) { crem_0[i*p1]=0; crem_0[i*p2]=0; }
    for (     ; i < p2; i++) { crem_0[i*p1]=0;                }
    for (i = 0; i < p3; i++) { crem_1[i*p3]=0; crem_1[i*p4]=0; }
    for (     ; i < p4; i++) { crem_1[i*p3]=0;                }
    for (i = 0; i < p5; i++) { crem_2[i*p5]=0; crem_2[i*p6]=0; }
    for (     ; i < p6; i++) { crem_2[i*p5]=0;                }
    for (c = 1; c < nc; c++) {
      uint32_t c1=cl[c], c2=cl[c], c3=cl[c], c4=cl[c], c5=cl[c], c6=cl[c];
      if (c1 >= p1) c1 %= p1;
      if (c2 >= p2) c2 %= p2;
      for (i = 1; i <= p1; i++) { crem_0[i*p1-c1]=0; crem_0[i*p2-c2]=0; }
      for (     ; i <= p2; i++) { crem_0[i*p1-c1]=0;                   }
      if (c3 >= p3) c3 %= p3;
      if (c4 >= p4) c4 %= p4;
      for (i = 1; i <= p3; i++) { crem_1[i*p3-c3]=0; crem_1[i*p4-c4]=0; }
      for (     ; i <= p4; i++) { crem_1[i*p3-c3]=0;                   }
      if (c5 >= p5) c5 %= p5;
      if (c6 >= p6) c6 %= p6;
      for (i = 1; i <= p5; i++) { crem_2[i*p5-c5]=0; crem_2[i*p6-c6]=0; }
      for (     ; i <= p6; i++) { crem_2[i*p5-c5]=0;                   }
    }
    New(0, resmod_0, nres, uint32_t);
    New(0, resmod_1, nres, uint32_t);
    New(0, resmod_2, nres, uint32_t);
    for (i = 0; i < nres; i++) {
      resmod_0[i] = residues[i] % pp_0;
      resmod_1[i] = residues[i] % pp_1;
      resmod_2[i] = residues[i] % pp_2;
    }
  }

  /* Precalculate acceptable residues for more primes */
  New(0, VPrem, maxpi, char*);
  memset(VPrem, 0, maxpi);
  for (pi = startpi+6; pi < maxpi; pi++) {
    uint32_t p = sprimes[pi];
    New(0, VPrem[pi], p, char);
    memset(VPrem[pi], 1, p);
  }
  for (pi = startpi+6, smallnc = 0; pi < maxpi; pi++) {
    uint32_t p = sprimes[pi];
    char* prem = VPrem[pi];
    prem[0] = 0;
    while (smallnc < nc && cl[smallnc] < p)   smallnc++;
    for (c = 1; c < smallnc; c++) prem[p-cl[c]] = 0;
    for (     ; c <      nc; c++) prem[p-(cl[c]%p)] = 0;
  }

  New(0, cres, nres, UV);

  rem_0 = low % pp_0;  remadd_0 = ppr % pp_0;
  rem_1 = low % pp_1;  remadd_1 = ppr % pp_1;
  rem_2 = low % pp_2;  remadd_2 = ppr % pp_2;

  /* Loop over their range in chunks of size 'ppr' */
  while (low <= high) {
    uint32_t r, nr, remr, ncres;

    /* Reduce the allowed residues for this chunk using more primes */

    { /* Start making a list of this chunk's residues using three pairs */
      for (r = 0, ncres = 0; r < nres; r++) {
        addmodded(remr, rem_0, resmod_0[r], pp_0);
        if (crem_0[remr]) {
          addmodded(remr, rem_1, resmod_1[r], pp_1);
          if (crem_1[remr]) {
            addmodded(remr, rem_2, resmod_2[r], pp_2);
            if (crem_2[remr]) {
              cres[ncres++] = residues[r];
            }
          }
        }
      }
      addmodded(rem_0, rem_0, remadd_0, pp_0);
      addmodded(rem_1, rem_1, remadd_1, pp_1);
      addmodded(rem_2, rem_2, remadd_2, pp_2);
    }

    /* Sieve through more primes one at a time, removing residues. */
    for (pi = startpi+6; pi < maxpi && ncres > 0; pi++) {
      uint32_t p = sprimes[pi];
      uint32_t rem = low % p;
      char* prem = VPrem[pi];
      /* Check divisibility of each remaining residue with this p */
      /* If we extended prem we could remove the add in the loop below */
      if (startpi <= 9) {   /* Residues are 32-bit */
        for (r = 0, nr = 0; r < ncres; r++) {
          if (prem[ (rem+(uint32_t)cres[r]) % p ])
            cres[nr++] = cres[r];
        }
      } else {              /* Residues are 64-bit */
        for (r = 0, nr = 0; r < ncres; r++) {
          if (prem[ (rem+cres[r]) % p ])
            cres[nr++] = cres[r];
        }
      }
      ncres = nr;
    }
    if (_verbose > 2) printf("cluster sieve range has %u residues left\n", ncres);

    /* Now check each of the remaining residues for inclusion */
    for (r = 0; r < ncres; r++) {
      UV p = low + cres[r];
      if (p > high) break;
      /* PRP test.  Split to save time. */
      for (c = 0; c < nc; c++)
        if (num_mr++,!is_euler_plumb_pseudoprime(p+cl[c]))
          break;
      if (c < nc) continue;
      for (c = 0; c < nc; c++)
        if (num_lucas++,!is_almost_extra_strong_lucas_pseudoprime(p+cl[c], 1))
          break;
      if (c < nc) continue;
      PUSH_VLIST(retlist, p);
    }
    low += ppr;
    if (low < ppr) low = UV_MAX;
  }

  if (_verbose) printf("cluster sieve ran %lu MR and %lu Lucas tests\n", num_mr, num_lucas);
  for (pi = startpi+6; pi < maxpi; pi++)
    Safefree(VPrem[pi]);
  Safefree(VPrem);
  Safefree(resmod_0);
  Safefree(resmod_1);
  Safefree(resmod_2);
  Safefree(cres);
  Safefree(residues);
  *numret = retlist.nsize;
  return retlist.list;
}
