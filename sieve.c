#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define FUNC_isqrt 1
#define FUNC_next_prime_in_sieve
#include "sieve.h"
#include "ptypes.h"
#include "cache.h"
#include "util.h"
#include "primality.h"
#include "montmath.h"
#include "prime_counts.h"

/* Is it better to do a partial sieve + primality tests vs. full sieve? */
static bool do_partial_sieve(UV startp, UV endp) {
  UV range = endp - startp;
  if (USE_MONTMATH) range /= 8;  /* Fast primality tests */
#if BITS_PER_WORD == 64
  if ( (startp > UVCONST(     100000000000000) && range <     40000) ||
       (startp > UVCONST(    1000000000000000) && range <    150000) ||
       (startp > UVCONST(   10000000000000000) && range <    600000) ||
       (startp > UVCONST(  100000000000000000) && range <   2500000) ||
       (startp > UVCONST( 1000000000000000000) && range <  10000000) ||
       (startp > UVCONST(10000000000000000000) && range <  40000000) )
    return 1;
#endif
  return 0;
}

/* 1001 bytes of presieved mod-30 bytes.  If the area to be sieved is
 * appropriately filled with this data, then 7, 11, and 13 do not have
 * to be sieved.  It wraps, so multiple memcpy's can be used.  Do be
 * aware that if you start at 0, you'll have to correct the first byte.
 *
 * mpu '$g=7*11*13; @b=(0)x$g; for $d (0..$g-1) { $i=0; for $m (1,7,11,13,17,19,23,29) { $n=30*$d+$m; if (gcd($n,$g) != 1) { $b[$d] |= (1<<$i); } $i++; } } for (0..$#b) { printf "0x%02x,",$b[$_]; print "\n" unless ($_+1)%13; } print "\n"'
 */
#define PRESIEVE_SIZE (7*11*13)
static const unsigned char presieve13[PRESIEVE_SIZE] =
{ 0x0e,0x20,0x10,0x81,0x49,0x24,0xc2,0x06,0x2a,0x90,0xa1,0x0c,0x14,
  0x58,0x02,0x61,0x11,0xc3,0x28,0x0c,0x44,0x22,0xa4,0x10,0x91,0x18,
  0x4d,0x40,0x82,0x21,0x58,0xa1,0x28,0x04,0x42,0x92,0x20,0x51,0x91,
  0x8a,0x04,0x48,0x03,0x60,0x34,0x81,0x1c,0x06,0xc1,0x02,0xa2,0x10,
  0x89,0x08,0x24,0x45,0x42,0x30,0x10,0xc5,0x0a,0x86,0x40,0x0a,0x30,
  0x38,0x85,0x08,0x15,0x40,0x63,0x20,0x96,0x83,0x88,0x04,0x60,0x16,
  0x28,0x10,0x81,0x49,0x44,0xe2,0x02,0x2c,0x12,0xa1,0x0c,0x04,0x50,
  0x0a,0x61,0x10,0x83,0x48,0x2c,0x40,0x26,0x26,0x90,0x91,0x08,0x55,
  0x48,0x82,0x20,0x19,0xc1,0x28,0x04,0x44,0x12,0xa0,0x51,0x81,0x9a,
  0x0c,0x48,0x02,0x21,0x54,0xa1,0x18,0x04,0x43,0x82,0xa2,0x10,0x99,
  0x08,0x24,0x44,0x03,0x70,0x30,0xc1,0x0c,0x86,0xc0,0x0a,0x20,0x30,
  0x8d,0x08,0x14,0x41,0x43,0x20,0x92,0x85,0x0a,0x84,0x60,0x06,0x30,
  0x18,0x81,0x49,0x05,0xc2,0x22,0x28,0x14,0xa3,0x8c,0x04,0x50,0x12,
  0x69,0x10,0x83,0x09,0x4c,0x60,0x22,0x24,0x12,0x91,0x08,0x45,0x50,
  0x8a,0x20,0x18,0x81,0x68,0x24,0x40,0x16,0x22,0xd1,0x81,0x8a,0x14,
  0x48,0x02,0x20,0x15,0xc1,0x38,0x04,0x45,0x02,0xa2,0x10,0x89,0x18,
  0x2c,0x44,0x02,0x31,0x50,0xe1,0x08,0x86,0x42,0x8a,0x20,0x30,0x95,
  0x08,0x14,0x40,0x43,0x60,0xb2,0x81,0x0c,0x06,0xe0,0x06,0x20,0x10,
  0x89,0x49,0x04,0xc3,0x42,0x28,0x10,0xa5,0x0e,0x84,0x50,0x02,0x71,
  0x18,0x83,0x08,0x0d,0x40,0x22,0x24,0x14,0x93,0x88,0x45,0x40,0x92,
  0x28,0x18,0x81,0x29,0x44,0x60,0x12,0x24,0x53,0x81,0x8a,0x04,0x58,
  0x0a,0x20,0x14,0x81,0x58,0x24,0x41,0x06,0xa2,0x90,0x89,0x08,0x34,
  0x4c,0x02,0x30,0x11,0xc1,0x28,0x86,0x44,0x0a,0xa0,0x30,0x85,0x18,
  0x1c,0x40,0x43,0x21,0xd2,0xa1,0x08,0x04,0x62,0x86,0x20,0x10,0x91,
  0x49,0x04,0xc2,0x03,0x68,0x30,0xa1,0x0c,0x06,0xd0,0x02,0x61,0x10,
  0x8b,0x08,0x0c,0x41,0x62,0x24,0x10,0x95,0x0a,0xc5,0x40,0x82,0x30,
  0x18,0x81,0x28,0x05,0x40,0x32,0x20,0x55,0x83,0x8a,0x04,0x48,0x12,
  0x28,0x14,0x81,0x19,0x44,0x61,0x02,0xa6,0x12,0x89,0x08,0x24,0x54,
  0x0a,0x30,0x10,0xc1,0x48,0xa6,0x40,0x0e,0x22,0xb0,0x85,0x08,0x14,
  0x48,0x43,0x20,0x93,0xc1,0x28,0x04,0x64,0x06,0xa0,0x10,0x81,0x59,
  0x0c,0xc2,0x02,0x29,0x50,0xa1,0x0c,0x04,0x52,0x82,0x61,0x10,0x93,
  0x08,0x0c,0x40,0x23,0x64,0x30,0x91,0x0c,0x47,0xc0,0x82,0x20,0x18,
  0x89,0x28,0x04,0x41,0x52,0x20,0x51,0x85,0x8a,0x84,0x48,0x02,0x30,
  0x1c,0x81,0x18,0x05,0x41,0x22,0xa2,0x14,0x8b,0x88,0x24,0x44,0x12,
  0x38,0x10,0xc1,0x09,0xc6,0x60,0x0a,0x24,0x32,0x85,0x08,0x14,0x50,
  0x4b,0x20,0x92,0x81,0x48,0x24,0x60,0x06,0x22,0x90,0x81,0x49,0x14,
  0xca,0x02,0x28,0x11,0xe1,0x2c,0x04,0x54,0x02,0xe1,0x10,0x83,0x18,
  0x0c,0x40,0x22,0x25,0x50,0xb1,0x08,0x45,0x42,0x82,0x20,0x18,0x91,
  0x28,0x04,0x40,0x13,0x60,0x71,0x81,0x8e,0x06,0xc8,0x02,0x20,0x14,
  0x89,0x18,0x04,0x41,0x42,0xa2,0x10,0x8d,0x0a,0xa4,0x44,0x02,0x30,
  0x18,0xc1,0x08,0x87,0x40,0x2a,0x20,0x34,0x87,0x88,0x14,0x40,0x53,
  0x28,0x92,0x81,0x09,0x44,0x60,0x06,0x24,0x12,0x81,0x49,0x04,0xd2,
  0x0a,0x28,0x10,0xa1,0x4c,0x24,0x50,0x06,0x63,0x90,0x83,0x08,0x1c,
  0x48,0x22,0x24,0x11,0xd1,0x28,0x45,0x44,0x82,0xa0,0x18,0x81,0x38,
  0x0c,0x40,0x12,0x21,0x51,0xa1,0x8a,0x04,0x4a,0x82,0x20,0x14,0x91,
  0x18,0x04,0x41,0x03,0xe2,0x30,0x89,0x0c,0x26,0xc4,0x02,0x30,0x10,
  0xc9,0x08,0x86,0x41,0x4a,0x20,0x30,0x85,0x0a,0x94,0x40,0x43,0x30,
  0x9a,0x81,0x08,0x05,0x60,0x26,0x20,0x14,0x83,0xc9,0x04,0xc2,0x12,
  0x28,0x10,0xa1,0x0d,0x44,0x70,0x02,0x65,0x12,0x83,0x08,0x0c,0x50,
  0x2a,0x24,0x10,0x91,0x48,0x65,0x40,0x86,0x22,0x98,0x81,0x28,0x14,
  0x48,0x12,0x20,0x51,0xc1,0xaa,0x04,0x4c,0x02,0xa0,0x14,0x81,0x18,
  0x0c,0x41,0x02,0xa3,0x50,0xa9,0x08,0x24,0x46,0x82,0x30,0x10,0xd1,
  0x08,0x86,0x40,0x0b,0x60,0x30,0x85,0x0c,0x16,0xc0,0x43,0x20,0x92,
  0x89,0x08,0x04,0x61,0x46,0x20,0x10,0x85,0x4b,0x84,0xc2,0x02,0x38,
  0x18,0xa1,0x0c,0x05,0x50,0x22,0x61,0x14,0x83,0x88,0x0c,0x40,0x32,
  0x2c,0x10,0x91,0x09,0x45,0x60,0x82,0x24,0x1a,0x81,0x28,0x04,0x50,
  0x1a,0x20,0x51,0x81,0xca,0x24,0x48,0x06,0x22,0x94,0x81,0x18,0x14,
  0x49,0x02,0xa2,0x11,0xc9,0x28,0x24,0x44,0x02,0xb0,0x10,0xc1,0x18,
  0x8e,0x40,0x0a,0x21,0x70,0xa5,0x08,0x14,0x42,0xc3,0x20,0x92,0x91,
  0x08,0x04,0x60,0x07,0x60,0x30,0x81,0x4d,0x06,0xc2,0x02,0x28,0x10,
  0xa9,0x0c,0x04,0x51,0x42,0x61,0x10,0x87,0x0a,0x8c,0x40,0x22,0x34,
  0x18,0x91,0x08,0x45,0x40,0xa2,0x20,0x1c,0x83,0xa8,0x04,0x40,0x12,
  0x28,0x51,0x81,0x8b,0x44,0x68,0x02,0x24,0x16,0x81,0x18,0x04,0x51,
  0x0a,0xa2,0x10,0x89,0x48,0x24,0x44,0x06,0x32,0x90,0xc1,0x08,0x96,
  0x48,0x0a,0x20,0x31,0xc5,0x28,0x14,0x44,0x43,0xa0,0x92,0x81,0x18,
  0x0c,0x60,0x06,0x21,0x50,0xa1,0x49,0x04,0xc2,0x82,0x28,0x10,0xb1,
  0x0c,0x04,0x50,0x03,0x61,0x30,0x83,0x0c,0x0e,0xc0,0x22,0x24,0x10,
  0x99,0x08,0x45,0x41,0xc2,0x20,0x18,0x85,0x2a,0x84,0x40,0x12,0x30,
  0x59,0x81,0x8a,0x05,0x48,0x22,0x20,0x14,0x83,0x98,0x04,0x41,0x12,
  0xaa,0x10,0x89,0x09,0x64,0x64,0x02,0x34,0x12,0xc1,0x08,0x86,0x50,
  0x0a,0x20,0x30,0x85,0x48,0x34,0x40,0x47,0x22,0x92,0x81,0x08,0x14,
  0x68,0x06,0x20,0x11,0xc1,0x69,0x04,0xc6,0x02,0xa8,0x10,0xa1,0x1c,
  0x0c,0x50,0x02,0x61,0x50,0xa3,0x08,0x0c,0x42,0xa2,0x24,0x10,0x91,
  0x08,0x45,0x40,0x83,0x60,0x38,0x81,0x2c,0x06,0xc0,0x12,0x20,0x51,
  0x89,0x8a,0x04,0x49,0x42,0x20,0x14,0x85,0x1a,0x84,0x41,0x02,0xb2,
  0x18,0x89,0x08,0x25,0x44,0x22,0x30,0x14,0xc3,0x88,0x86,0x40,0x1a,
  0x28,0x30,0x85,0x09,0x54,0x60,0x43,0x24,0x92,0x81,0x08,0x04,0x70};

static const UV max_sieve_prime = (BITS_PER_WORD==64) ? 4294967291U : 65521U;

/* Tile <from> bytes from source to <to> bytes in dest */
static void memtile(unsigned char* dst, const unsigned char* src, size_t from, size_t to) {
  if (to < from)
    from = to;
  if (dst != src)
    memcpy(dst, src, from);
  while (from < to) {
    size_t bytes = (2*from > to) ? to-from : from;
    memcpy(dst+from, dst, bytes);
    from += bytes;
  }
}

static UV sieve_prefill(unsigned char* mem, UV startd, UV endd)
{
  UV vnext_prime = 17;
  UV nbytes = endd - startd + 1;
  MPUassert( (mem != 0) && (endd >= startd), "sieve_prefill bad arguments");

  if (startd != 0) {
    UV pstartd = startd % PRESIEVE_SIZE;
    UV tailbytes = PRESIEVE_SIZE - pstartd;
    if (tailbytes > nbytes) tailbytes = nbytes;
    memcpy(mem, presieve13 + pstartd, tailbytes); /* Copy tail to mem */
    mem += tailbytes;    /* Advance so mem points at the beginning */
    nbytes -= tailbytes;
  }
  if (nbytes > 0) {
    memcpy(mem, presieve13, (nbytes < PRESIEVE_SIZE) ? nbytes : PRESIEVE_SIZE);
    memtile(mem, mem, PRESIEVE_SIZE, nbytes);
    if (startd == 0) mem[0] = 0x01; /* Correct first byte */
  }
  /* Leaving option open to tile 17 out and sieve, then return 19 */
  return vnext_prime;
}

/* Marking primes is done the same way we used to do with tables, but
 * now uses heavily unrolled code based on Kim Walisch's mod-30 sieve.
 */
#define set_bit(s,n)  *(s) |= (1 << n);
static const unsigned char masknum30[30] =
    {0,0,0,0,0,0,0,1,0,0,0,2,0,3,0,0,0,4,0,5,0,0,0,6,0,0,0,0,0,7};
static const unsigned char qinit30[30] =
    {0,0,1,1,1,1,1,1,2,2,2,2,3,3,4,4,4,4,5,5,6,6,6,6,7,7,7,7,7,7};

typedef struct {
  uint32_t prime;
  UV       offset;
  uint8_t  index;
} wheel_t;

#define CROSS_INDEX(v, b0,b1,b2,b3,b4,b5,b6,b7,  i0,i1,i2,i3,i4,i5,i6,i7, it) \
  while (1) { \
    case (v+0): if(s>=send){w->index=v+0;break;} set_bit(s,b0); s += r*6+i0; \
    case (v+1): if(s>=send){w->index=v+1;break;} set_bit(s,b1); s += r*4+i1; \
    case (v+2): if(s>=send){w->index=v+2;break;} set_bit(s,b2); s += r*2+i2; \
    case (v+3): if(s>=send){w->index=v+3;break;} set_bit(s,b3); s += r*4+i3; \
    case (v+4): if(s>=send){w->index=v+4;break;} set_bit(s,b4); s += r*2+i4; \
    case (v+5): if(s>=send){w->index=v+5;break;} set_bit(s,b5); s += r*4+i5; \
    case (v+6): if(s>=send){w->index=v+6;break;} set_bit(s,b6); s += r*6+i6; \
    case (v+7): if(s>=send){w->index=v+7;break;} set_bit(s,b7); s += r*2+i7; \
    while (s + r*28 + it-1 < send) { \
      set_bit(s + r *  0 +  0, b0); \
      set_bit(s + r *  6 + i0, b1); \
      set_bit(s + r * 10 + i0+i1, b2); \
      set_bit(s + r * 12 + i0+i1+i2, b3); \
      set_bit(s + r * 16 + i0+i1+i2+i3, b4); \
      set_bit(s + r * 18 + i0+i1+i2+i3+i4, b5); \
      set_bit(s + r * 22 + i0+i1+i2+i3+i4+i5, b6); \
      set_bit(s + r * 28 + i0+i1+i2+i3+i4+i5+i6, b7); \
      s += r*30 + it; \
    } \
  }
static wheel_t create_wheel(UV startp, uint32_t prime)
{
  wheel_t w;
  UV q = prime;
  UV p2 = q*q;

  if (startp == 0) {
    wheel_t ws = { prime, p2/30, qinit30[q % 30]  +  8*masknum30[prime % 30] };
    return ws;
  }

  if (p2 < startp) {
    q = 1+(startp-1)/prime;
    q += distancewheel30[q % 30];
    p2 = prime * q;
    /* The offset if p2 overflows is still ok, or set to max_sieve_prime+1. */
    /* if (p2 < startp) p2 = max_sieve_prime+1; */
  }

  w.offset = (p2-startp) / 30;
  w.index = qinit30[q % 30]  +  8*masknum30[prime % 30];
  w.prime = prime;
  return w;
}

static void mark_primes(unsigned char* s, UV bytes, wheel_t* w)
{
  if (w->offset >= bytes) {
    w->offset -= bytes;
  } else {
    const unsigned char* send = s + bytes;
    uint32_t r = w->prime / 30;
    s += w->offset;
    switch (w->index) {
      CROSS_INDEX( 0, 0,1,2,3,4,5,6,7, 0,0,0,0,0,0,0,1,  1); break;
      CROSS_INDEX( 8, 1,5,4,0,7,3,2,6, 1,1,1,0,1,1,1,1,  7); break;
      CROSS_INDEX(16, 2,4,0,6,1,7,3,5, 2,2,0,2,0,2,2,1, 11); break;
      CROSS_INDEX(24, 3,0,6,5,2,1,7,4, 3,1,1,2,1,1,3,1, 13); break;
      CROSS_INDEX(32, 4,7,1,2,5,6,0,3, 3,3,1,2,1,3,3,1, 17); break;
      CROSS_INDEX(40, 5,3,7,1,6,0,4,2, 4,2,2,2,2,2,4,1, 19); break;
      CROSS_INDEX(48, 6,2,3,7,0,4,5,1, 5,3,1,4,1,3,5,1, 23); break;
      CROSS_INDEX(56, 7,6,5,4,3,2,1,0, 6,4,2,4,2,4,6,1, 29); break;
    }
    w->offset = s - send;
  }
}

/* Monolithic mod-30 wheel sieve */
unsigned char* sieve_erat30(UV end)
{
  unsigned char *mem;
  UV max_buf, limit, prime;

  max_buf = (end/30) + ((end%30) != 0);
  /* Round up to a word */
  max_buf = ((max_buf + sizeof(UV) - 1) / sizeof(UV)) * sizeof(UV);
  New(0, mem, max_buf, unsigned char );

  /* Fill buffer with marked 7, 11, and 13 */
  prime = sieve_prefill(mem, 0, max_buf-1);

  limit = isqrt(end);  /* prime*prime can overflow */
  for (  ; prime <= limit; prime = next_prime_in_sieve(mem,prime,end)) {
    wheel_t w = create_wheel(0, prime);
    mark_primes(mem, max_buf, &w);
  }
  return mem;
}

static void _primality_test_sieve(unsigned char* mem, UV startp, UV endp) {
  START_DO_FOR_EACH_SIEVE_PRIME(mem, 0, 0, endp-startp) {
    if (!BPSW(startp + p))           /* If the candidate is not prime, */
      mem[p/30] |= masktab30[p%30];  /* mark the sieve location.       */
  } END_DO_FOR_EACH_SIEVE_PRIME;
}
static void _sieve_range(unsigned char* mem, const unsigned char* sieve, UV startd, UV endd, UV limit) {
  UV startp = 30*startd;
  UV start_base_prime = sieve_prefill(mem, startd, endd);
  START_DO_FOR_EACH_SIEVE_PRIME(sieve, 0, start_base_prime, limit) { /* Sieve */
    wheel_t w = create_wheel(startp, p);
    mark_primes(mem, endd-startd+1, &w);
  } END_DO_FOR_EACH_SIEVE_PRIME;
}

bool sieve_segment_partial(unsigned char* mem, UV startd, UV endd, UV depth)
{
  const unsigned char* sieve;
  UV startp = 30*startd,  endp = (endd >= (UV_MAX/30)) ? UV_MAX-2 : 30*endd+29;
  UV limit = isqrt(endp);
  MPUassert(mem != 0 && endd >= startd && endp >= startp && depth >= 13,
            "sieve_segment_partial bad arguments");
  /* limit = min( sqrt(end), max-64-bit-prime, requested depth ) */
  if (limit > max_sieve_prime)  limit = max_sieve_prime;
  if (limit > depth)  limit = depth;
  get_prime_cache(limit, &sieve);                      /* Get sieving primes  */
  _sieve_range(mem, sieve, startd, endd, limit);
  release_prime_cache(sieve);
  return 1;
}

/* Segmented mod-30 wheel sieve */
bool sieve_segment(unsigned char* mem, UV startd, UV endd)
{
  const unsigned char* sieve;
  UV startp = 30*startd,  endp = (endd >= (UV_MAX/30)) ? UV_MAX-2 : 30*endd+29;
  UV sieve_size,  limit = isqrt(endp);
  int do_partial = do_partial_sieve(startp, endp);
  MPUassert(mem != 0 && endd >= startd && endp >= startp,
            "sieve_segment bad arguments");

  sieve_size = get_prime_cache(0, &sieve);

  if (sieve_size >= endp) {

    /* We can just use the primary cache */
    memcpy(mem, sieve+startd, endd-startd+1);
    release_prime_cache(sieve);

  } else if (!do_partial && sieve_size >= limit) {

    /* Full sieve and we have all sieving primes in hand */
    _sieve_range(mem, sieve, startd, endd, limit);
    release_prime_cache(sieve);

  } else {

    release_prime_cache(sieve);
    if (do_partial)
      limit >>= ((startp < (UV)1e16) ? 8 : 10);

    /* sieve_segment_partial(mem, startd, endd, limit); */
    get_prime_cache(limit, &sieve);
    _sieve_range(mem, sieve, startd, endd, limit);
    release_prime_cache(sieve);

    if (do_partial)
      _primality_test_sieve(mem, startp, endp);

  }
  return 1;
}

bool sieve_segment_wheel(unsigned char* mem, UV startd, UV endd, wheel_t *warray, uint32_t wsize)
{
  uint32_t i = 0, limit, start_base_prime;
  uint32_t segsize = endd - startd + 1;
  UV startp = 30*startd;
  UV endp = (endd >= (UV_MAX/30))  ?  UV_MAX-2  :  30*endd+29;
  MPUassert(mem != 0 && endd >= startd && endp >= startp,
            "sieve_segment bad arguments");

  /* possibly use primary cache directly */

  /* Fill buffer with marked 7, 11, and 13 */
  start_base_prime = sieve_prefill(mem, startd, endd);
  while (i < wsize && warray[i].prime < start_base_prime)
    i++;

  limit = isqrt(endp);
  if (limit > max_sieve_prime)  limit = max_sieve_prime;

  while (i < wsize && warray[i].prime <= limit) {
    if (warray[i].index >= 64)
      warray[i] = create_wheel(startp, warray[i].prime);
    mark_primes(mem, segsize, &(warray[i++]));
  }

  if (limit > warray[wsize-1].prime && warray[wsize-1].prime < max_sieve_prime)
    _primality_test_sieve(mem, startp, endp);

  return 1;
}

/**************************************************************************/

static UV simple_prime_count_upper(UV n) {
  double pc, logn = log(n);
  if (n < 5)                return 0 + (n>1) + (n>2);
  if (n < 355991)           pc = n / (logn-1.112);
  else if (n < 2953652287U) pc = n / logn * (1 + 1/logn + 2.51 / (logn*logn));
  else                      pc = n / logn * (1 + 1/logn + 2.334 / (logn*logn));
  return (UV) ceil(pc);
}

typedef struct {
  UV lod;
  UV hid;
  UV low;
  UV high;
  UV endp;
  UV segment_size;
  unsigned char* segment;
  unsigned char* base;
  wheel_t *warray;
  uint32_t wsize;
} segment_context_t;

/*
 * unsigned char* segment;
 * UV seg_base, seg_low, seg_high;
 * void* ctx = start_segment_primes(low, high, &segment);
 * while (beg < 7) {
 *   beg = (beg <= 2) ? 2 : (beg <= 3) ? 3 : 5;
 *   .... with beg ....
 *   beg += 1 + (beg > 2);
 * }
 * while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
 *   START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_low - seg_base, seg_high - seg_base )
 *     .... with seg_base + p ....
 *   END_DO_FOR_EACH_SIEVE_PRIME
 * }
 * end_segment_primes(ctx);
 */

void* start_segment_primes(UV low, UV high, unsigned char** segmentmem)
{
  segment_context_t* ctx;
  UV nsegments, range;

  MPUassert( high >= low, "start_segment_primes bad arguments");
  New(0, ctx, 1, segment_context_t);
  ctx->low = low;
  ctx->high = high;
  ctx->lod = low / 30;
  ctx->hid = high / 30;
  ctx->endp = (ctx->hid >= (UV_MAX/30))  ?  UV_MAX-2  :  30*ctx->hid+29;
  range = ctx->hid - ctx->lod + 1;   /* range in bytes */

#if BITS_PER_WORD == 64
  if (high > 1e10 && range > 32*1024-16) {
    UV size, div;
    /* Use larger segments */
    size = isqrt(32*isqrt(high)) * (logint(high,2)-2);
    if (size < 128*1024) size = 128*1024;
    /* Evenly split the range into segments */
    div = (range+size-1)/size;
    size = (div <= 1)  ?  range  :  (range+div-1)/div;
    ctx->segment_size = size;
    New(0, ctx->segment, size, unsigned char);
  } else
#endif
  ctx->segment = get_prime_segment( &(ctx->segment_size) );
  *segmentmem = ctx->segment;
  nsegments = (((high-low+29)/30)+ctx->segment_size-1) / ctx->segment_size;

  MPUverbose(3, "segment sieve: byte range %lu split into %lu segments of size %lu\n", (unsigned long)range, (unsigned long)nsegments, (unsigned long)ctx->segment_size);

  ctx->base = 0;
  ctx->warray = 0;
  ctx->wsize = 0;
#if 1
  { /* Generate wheel data for this segment sieve */
    const UV maxsieve = UVCONST(400000000);
    UV limit, nprimes;
    wheel_t *warray;
    wheel_t w = {0,0,128};
    uint32_t wsize = 0;
    /* Number of primes for a full sieve */
    limit = isqrt(ctx->endp);
    /* For small ranges a partial sieve is much faster */
    if (do_partial_sieve(low, high))
      limit >>= ((low < (UV)1e16) ? 8 : 10);
    if (limit <= maxsieve) {
      /* Bump to one more than needed. */
      limit = next_prime(limit);
      /* We'll make space for this many */
      nprimes = simple_prime_count_upper(limit);
      MPUverbose(4, "segment sieve %lu - %lu, primes to %lu (max %lu)\n", (unsigned long)low, (unsigned long)high, (unsigned long)limit, (unsigned long)nprimes);
      New(0, warray, nprimes, wheel_t);
      START_DO_FOR_EACH_PRIME(0,limit) {
        if (wsize >= nprimes) croak("segment bad upper count");
        w.prime = p;
        warray[wsize++] = w;
      } END_DO_FOR_EACH_PRIME;
      ctx->warray = warray;
      ctx->wsize = wsize;
    }
  }
#endif

  return (void*) ctx;
}

bool next_segment_primes(void* vctx, UV* base, UV* low, UV* high)
{
  UV seghigh_d, range_d;
  segment_context_t* ctx = (segment_context_t*) vctx;

  if (ctx->lod > ctx->hid) return 0;

  seghigh_d = ((ctx->hid - ctx->lod) < ctx->segment_size)
            ? ctx->hid
           : (ctx->lod + ctx->segment_size - 1);
  range_d = seghigh_d - ctx->lod + 1;
  *low = ctx->low;
  *high = (seghigh_d == ctx->hid) ? ctx->high : (seghigh_d*30 + 29);
  *base = ctx->lod * 30;

  MPUassert( seghigh_d >= ctx->lod, "next_segment_primes: highd < lowd");
  MPUassert( range_d <= ctx->segment_size, "next_segment_primes: range > segment size");

  if (ctx->warray != 0)
    sieve_segment_wheel(ctx->segment, ctx->lod, seghigh_d, ctx->warray, ctx->wsize);
  else
    sieve_segment(ctx->segment, ctx->lod, seghigh_d);

  ctx->lod += range_d;
  ctx->low = *high + 2;

  return 1;
}

void end_segment_primes(void* vctx)
{
  segment_context_t* ctx = (segment_context_t*) vctx;
  MPUassert(ctx != 0, "end_segment_primes given a null pointer");
  if (ctx->segment != 0) {
    release_prime_segment(ctx->segment);
    ctx->segment = 0;
  }
  if (ctx->base != 0) {
    Safefree(ctx->base);
    ctx->base = 0;
  }
  if (ctx->warray != 0) {
    Safefree(ctx->warray);
    ctx->warray = 0;
  }
  Safefree(ctx);
}

UV range_prime_sieve(UV**list, UV lo, UV hi)
{
  UV *P, Psize, i = 0;
  if (hi < lo) { *list = 0; return 0; }
  Psize = prime_count_upper(hi) - prime_count_lower(lo) + 1;
  New(0, P, Psize, UV);
  if (lo <= 2 && hi >= 2) P[i++] = 2;
  if (lo <= 3 && hi >= 3) P[i++] = 3;
  if (lo <= 5 && hi >= 5) P[i++] = 5;
  {
    unsigned char* segment;
    UV seg_base, seg_low, seg_high;
    void* ctx = start_segment_primes(lo, hi, &segment);
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
        P[i++] = p;
      END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
  }
  *list = P;
  return i;
}

uint32_t range_prime_sieve_32(uint32_t** list, uint32_t n, uint32_t offset)
{
  uint32_t *P, i = offset;

  if (n < 2) { *list = 0; return 0; }
  New(0, P, prime_count_upper(n) + offset + 3, uint32_t);   /* Allocate list */
  if (offset > 0)  memset(P, 0, offset * sizeof(uint32_t)); /* Zero to offset */
  P[i++] = 2;  P[i++] = 3;  P[i++] = 5;                     /* Fill in 2/3/5 */

  if (n >= 7) {
    unsigned char* segment;
    UV seg_base, seg_low, seg_high;
    void* ctx = start_segment_primes(7, n, &segment);
    while (next_segment_primes(ctx, &seg_base, &seg_low, &seg_high)) {
      START_DO_FOR_EACH_SIEVE_PRIME( segment, seg_base, seg_low, seg_high )
        P[i++] = p;
      END_DO_FOR_EACH_SIEVE_PRIME
    }
    end_segment_primes(ctx);
  }
  while (P[i-1] > n)  i--;  /* Truncate the count if necesssary. */
  *list = P;
  return i-offset;          /* Returns number of primes, excluding offset */
}
