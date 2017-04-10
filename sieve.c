#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define FUNC_next_prime_in_sieve
#include "sieve.h"
#include "ptypes.h"
#include "cache.h"
#define FUNC_isqrt 1
#include "util.h"
#include "primality.h"
#include "montmath.h"

/* Is it better to do a partial sieve + primality tests vs. full sieve? */
static int do_partial_sieve(UV startp, UV endp) {
  UV range = endp - startp;
  if (USE_MONTMATH) range /= 8;  /* Fast primality tests */
#if BITS_PER_WORD == 64
  if ( (startp > UVCONST(     100000000000000) && range <     20000) ||
       (startp > UVCONST(    1000000000000000) && range <    100000) ||
       (startp > UVCONST(   10000000000000000) && range <    200000) ||
       (startp > UVCONST(  100000000000000000) && range <   2000000) ||
       (startp > UVCONST( 1000000000000000000) && range <  10000000) ||
       (startp > UVCONST(10000000000000000000) && range <  20000000) )
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

static const unsigned char stepdata[8][8][8] = {
  { {96,65,34,67,36,69,102,47},
    {64,47,75,106,46,105,77,44},
    {32,86,33,87,115,45,114,84},
    {72,46,85,42,73,127,44,123},
    {96,35,100,79,33,66,37,78},
    {64,108,34,109,67,47,65,46},
    {32,76,117,33,118,74,35,87},
    {32,127,86,45,84,43,82,121},
  },
  { {65,34,67,36,69,102,47,96},
    {105,77,44,64,47,75,106,46},
    {33,87,115,45,114,84,32,86},
    {73,127,44,123,72,46,85,42},
    {33,66,37,78,96,35,100,79},
    {65,46,64,108,34,109,67,47},
    {33,118,74,35,87,32,76,117},
    {121,32,127,86,45,84,43,82},
  },
  { {34,67,36,69,102,47,96,65},
    {106,46,105,77,44,64,47,75},
    {114,84,32,86,33,87,115,45},
    {42,73,127,44,123,72,46,85},
    {66,37,78,96,35,100,79,33},
    {34,109,67,47,65,46,64,108},
    {74,35,87,32,76,117,33,118},
    {82,121,32,127,86,45,84,43},
  },
  { {67,36,69,102,47,96,65,34},
    {75,106,46,105,77,44,64,47},
    {115,45,114,84,32,86,33,87},
    {123,72,46,85,42,73,127,44},
    {35,100,79,33,66,37,78,96},
    {67,47,65,46,64,108,34,109},
    {35,87,32,76,117,33,118,74},
    {43,82,121,32,127,86,45,84},
  },
  { {36,69,102,47,96,65,34,67},
    {44,64,47,75,106,46,105,77},
    {84,32,86,33,87,115,45,114},
    {44,123,72,46,85,42,73,127},
    {100,79,33,66,37,78,96,35},
    {108,34,109,67,47,65,46,64},
    {76,117,33,118,74,35,87,32},
    {84,43,82,121,32,127,86,45},
  },
  { {69,102,47,96,65,34,67,36},
    {77,44,64,47,75,106,46,105},
    {45,114,84,32,86,33,87,115},
    {85,42,73,127,44,123,72,46},
    {37,78,96,35,100,79,33,66},
    {109,67,47,65,46,64,108,34},
    {117,33,118,74,35,87,32,76},
    {45,84,43,82,121,32,127,86},
  },
  { {102,47,96,65,34,67,36,69},
    {46,105,77,44,64,47,75,106},
    {86,33,87,115,45,114,84,32},
    {46,85,42,73,127,44,123,72},
    {78,96,35,100,79,33,66,37},
    {46,64,108,34,109,67,47,65},
    {118,74,35,87,32,76,117,33},
    {86,45,84,43,82,121,32,127},
  },
  { {47,96,65,34,67,36,69,102},
    {47,75,106,46,105,77,44,64},
    {87,115,45,114,84,32,86,33},
    {127,44,123,72,46,85,42,73},
    {79,33,66,37,78,96,35,100},
    {47,65,46,64,108,34,109,67},
    {87,32,76,117,33,118,74,35},
    {127,86,45,84,43,82,121,32},
  },
};

static const int wheelmap[30] =
  {0,0,0,0,0,0,0,1,0,0,0,2,0,3,0,0,0,4,0,5,0,0,0,6,0,0,0,0,0,7};
static const int wheel2xmap[30] =     /* (2*p)%30 => 2,14,22,26,4,8,16,28 */
  {0,0,0,0,4,0,0,0,5,0,0,0,0,0,1,0,6,0,0,0,0,0,2,0,0,0,3,0,7,0};
/*     2   4       8          14  16          22      26  28     (2*p)%30 */

#define FIND_COMPOSITE_POSITIONS(d, m, p) \
  do { \
    int v; \
    UV dinc = (2*p) / 30; \
    UV minc = (2*p) - dinc*30; \
    const unsigned char* steps = stepdata [wheelmap[m]] [wheel2xmap[minc]]; \
    v = steps[0]; wdinc[0] = dinc*(v>>5)+((v>>3)&0x3); wmask[0] = 1<<(v&0x7); \
    v = steps[1]; wdinc[1] = dinc*(v>>5)+((v>>3)&0x3); wmask[1] = 1<<(v&0x7); \
    v = steps[2]; wdinc[2] = dinc*(v>>5)+((v>>3)&0x3); wmask[2] = 1<<(v&0x7); \
    v = steps[3]; wdinc[3] = dinc*(v>>5)+((v>>3)&0x3); wmask[3] = 1<<(v&0x7); \
    v = steps[4]; wdinc[4] = dinc*(v>>5)+((v>>3)&0x3); wmask[4] = 1<<(v&0x7); \
    v = steps[5]; wdinc[5] = dinc*(v>>5)+((v>>3)&0x3); wmask[5] = 1<<(v&0x7); \
    v = steps[6]; wdinc[6] = dinc*(v>>5)+((v>>3)&0x3); wmask[6] = 1<<(v&0x7); \
    v = steps[7]; wdinc[7] = dinc*(v>>5)+((v>>3)&0x3); wmask[7] = 1<<(v&0x7); \
  } while (0)

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

/* Wheel 30 sieve.  Ideas from Terje Mathisen and Quesada / Van Pelt. */
unsigned char* sieve_erat30(UV end)
{
  unsigned char* mem;
  UV max_buf, limit, prime;

  max_buf = (end/30) + ((end%30) != 0);
  /* Round up to a word */
  max_buf = ((max_buf + sizeof(UV) - 1) / sizeof(UV)) * sizeof(UV);
  New(0, mem, max_buf, unsigned char );

  /* Fill buffer with marked 7, 11, and 13 */
  prime = sieve_prefill(mem, 0, max_buf-1);

  limit = isqrt(end);  /* prime*prime can overflow */
  for (  ; prime <= limit; prime = next_prime_in_sieve(mem,prime,end)) {
    UV p2 = prime*prime;
    UV d = p2 / 30;
    UV m = p2 - d*30;
    UV wdinc[8];
    unsigned char wmask[8];

    /* Find the positions of the next composites we will mark */
    FIND_COMPOSITE_POSITIONS(d, m, prime);
#if 0
    assert(d == ((prime*prime)/30));
    assert(d < max_buf);
    assert(prime = (wdinc[0]+wdinc[1]+wdinc[2]+wdinc[3]+wdinc[4]+wdinc[5]+wdinc[6]+wdinc[7]));
#endif
    /* Regular code to mark composites:
    *  i = 0;
    *  do {mem[d] |= wmask[i]; d += wdinc[i]; i = (i+1)&7;} while (d < max_buf);
    * Unrolled version: */
    while ( (d+prime) < max_buf ) {
      mem[d] |= wmask[0];  d += wdinc[0];
      mem[d] |= wmask[1];  d += wdinc[1];
      mem[d] |= wmask[2];  d += wdinc[2];
      mem[d] |= wmask[3];  d += wdinc[3];
      mem[d] |= wmask[4];  d += wdinc[4];
      mem[d] |= wmask[5];  d += wdinc[5];
      mem[d] |= wmask[6];  d += wdinc[6];
      mem[d] |= wmask[7];  d += wdinc[7];
    }
    while (1) {
      mem[d] |= wmask[0];  d += wdinc[0];  if (d >= max_buf) break;
      mem[d] |= wmask[1];  d += wdinc[1];  if (d >= max_buf) break;
      mem[d] |= wmask[2];  d += wdinc[2];  if (d >= max_buf) break;
      mem[d] |= wmask[3];  d += wdinc[3];  if (d >= max_buf) break;
      mem[d] |= wmask[4];  d += wdinc[4];  if (d >= max_buf) break;
      mem[d] |= wmask[5];  d += wdinc[5];  if (d >= max_buf) break;
      mem[d] |= wmask[6];  d += wdinc[6];  if (d >= max_buf) break;
      mem[d] |= wmask[7];  d += wdinc[7];  if (d >= max_buf) break;
    }
  }
  return mem;
}



int sieve_segment(unsigned char* mem, UV startd, UV endd)
{
  const unsigned char* sieve;
  UV limit, slimit, start_base_prime, sieve_size;
  UV startp = 30*startd;
  UV endp = (endd >= (UV_MAX/30))  ?  UV_MAX-2  :  30*endd+29;

  MPUassert( (mem != 0) && (endd >= startd) && (endp >= startp),
             "sieve_segment bad arguments");

  /* It's possible we can just use the primary cache */
  sieve_size = get_prime_cache(0, &sieve);
  if (sieve_size >= endp) {
    memcpy(mem, sieve+startd, endd-startd+1);
    release_prime_cache(sieve);
    return 1;
  }

  /* Fill buffer with marked 7, 11, and 13 */
  start_base_prime = sieve_prefill(mem, startd, endd);

  limit = isqrt(endp);  /* floor(sqrt(n)), will include p if p*p=endp */
  /* Don't use a sieve prime such that p*p > UV_MAX */
  if (limit > max_sieve_prime)  limit = max_sieve_prime;
  slimit = limit;
  if (do_partial_sieve(startp, endp))
    slimit >>= ((startp < (UV)1e16) ? 8 : 10);
  /* printf("segment sieve from %"UVuf" to %"UVuf" (aux sieve to %"UVuf")\n", startp, endp, slimit); */
  if (slimit > sieve_size) {
    release_prime_cache(sieve);
    get_prime_cache(slimit, &sieve);
  }

  START_DO_FOR_EACH_SIEVE_PRIME(sieve, 0, start_base_prime, slimit)
  {
    /* p increments from 17 to at most sqrt(endp).  Note on overflow:
     * 32-bit: limit=     65535, max p =      65521, p*p = ~0-1965854
     * 64-bit: limit=4294967295, max p = 4294967291, p*p = ~0-42949672934
     * No overflow here, but possible after the incrementing below. */
    UV p2 = p*p;
    if (p2 < startp) {
      UV f = 1+(startp-1)/p;
      p2 = p * (f + distancewheel30[f%30]);
    }
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
      UV d = p2 / 30;
      UV m = p2 - d*30;

      if ((p2 + 2*p) > endp) {
        /* There is only one composite to be marked in this segment */
        mem[d-startd] |= masktab30[m];
      } else {
        UV wdinc[8];
        unsigned char wmask[8];
        UV offset_endd = endd - startd;
        UV unrolls = (endd-d+1) / p;
        /* Find the positions of the next composites we will mark */
        FIND_COMPOSITE_POSITIONS(d, m, p);
        d -= startd;
        /* Unrolled inner loop for marking composites */
        while ( unrolls-- > 0) {
          mem[d] |= wmask[0];  d += wdinc[0];
          mem[d] |= wmask[1];  d += wdinc[1];
          mem[d] |= wmask[2];  d += wdinc[2];
          mem[d] |= wmask[3];  d += wdinc[3];
          mem[d] |= wmask[4];  d += wdinc[4];
          mem[d] |= wmask[5];  d += wdinc[5];
          mem[d] |= wmask[6];  d += wdinc[6];
          mem[d] |= wmask[7];  d += wdinc[7];
        }
        while (d <= offset_endd) {
          mem[d] |= wmask[0];  d += wdinc[0];  if (d > offset_endd) break;
          mem[d] |= wmask[1];  d += wdinc[1];  if (d > offset_endd) break;
          mem[d] |= wmask[2];  d += wdinc[2];  if (d > offset_endd) break;
          mem[d] |= wmask[3];  d += wdinc[3];  if (d > offset_endd) break;
          mem[d] |= wmask[4];  d += wdinc[4];  if (d > offset_endd) break;
          mem[d] |= wmask[5];  d += wdinc[5];  if (d > offset_endd) break;
          mem[d] |= wmask[6];  d += wdinc[6];  if (d > offset_endd) break;
          mem[d] |= wmask[7];  d += wdinc[7];
        }
      }
#endif
    }
  }
  END_DO_FOR_EACH_SIEVE_PRIME;
  release_prime_cache(sieve);

  if (limit > slimit) { /* We've sieved out most composites, but not all. */
    START_DO_FOR_EACH_SIEVE_PRIME(mem, 0, 0, endp-startp) {
      if (!BPSW(startp + p))        /* If the candidate is not prime, */
        mem[d_] |= 1 << bit_;       /* mark the sieve location.       */
    } END_DO_FOR_EACH_SIEVE_PRIME;
  }
  return 1;
}

/**************************************************************************/

typedef struct {
  UV lod;
  UV hid;
  UV low;
  UV high;
  UV endp;
  UV segment_size;
  unsigned char* segment;
  unsigned char* base;
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
  UV slimit;

  MPUassert( high >= low, "start_segment_primes bad arguments");
  New(0, ctx, 1, segment_context_t);
  ctx->low = low;
  ctx->high = high;
  ctx->lod = low / 30;
  ctx->hid = high / 30;
  ctx->endp = (ctx->hid >= (UV_MAX/30))  ?  UV_MAX-2  :  30*ctx->hid+29;

#if BITS_PER_WORD == 64
  if (high > 1e11 && high-low > 1e6) {
    UV range = (high-low+29)/30;
    /* Select what we think would be a good segment size */
    UV size = isqrt(isqrt(high)) * ((high < 1e15) ? 500 : 250);
    /* Evenly split the range into segments */
    UV div = (range+size-1)/size;
    size = (div <= 1)  ?  range  :  (range+div-1)/div;
    if (_XS_get_verbose() >= 2)
      printf("segment sieve: byte range %lu split into %lu segments of size %lu\n", (unsigned long)range, (unsigned long)div, (unsigned long)size);
    ctx->segment_size = size;
    New(0, ctx->segment, size, unsigned char);
  } else
#endif
  ctx->segment = get_prime_segment( &(ctx->segment_size) );
  *segmentmem = ctx->segment;

  ctx->base = 0;
  /* Expand primary cache so we won't regen each call */
  slimit = isqrt(ctx->endp)+1;
  if (do_partial_sieve(low, high))  slimit >>= 8;
  get_prime_cache( slimit, 0);

  return (void*) ctx;
}

int next_segment_primes(void* vctx, UV* base, UV* low, UV* high)
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
  Safefree(ctx);
}
