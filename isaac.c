/*
 * The ISAAC CSPRNG plus interface.
 * Slightly modified readable.c from Bob Jenkins 1996.
 */

#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include "ptypes.h"
#include "threadlock.h"
#include "isaac.h"

static int mutex_init = 0;
static int good_seed = 0;
MUTEX_DECL(state);

/*****************************************************************************/

static uint32_t randrsl[256];
static uint32_t randcnt;

/* internal state */
static uint32_t mm[256];
static uint32_t aa = 0, bb = 0, cc = 0;

static void isaac(void)
{
   uint32_t i,x,y;

   cc = cc + 1;    /* cc just gets incremented once per 256 results */
   bb = bb + cc;   /* then combined with bb */

   for (i=0; i<256; ++i)
   {
     x = mm[i];
     switch (i%4) {
       case 0: aa = aa^(aa<<13); break;
       case 1: aa = aa^(aa>>6); break;
       case 2: aa = aa^(aa<<2); break;
       case 3: aa = aa^(aa>>16); break;
     }
     aa              = mm[(i+128)%256] + aa;
     mm[i]      = y  = mm[(x>>2)%256] + aa + bb;
     randrsl[i] = bb = mm[(y>>10)%256] + x;
   }
   randcnt = 0;
}
#define mix(a,b,c,d,e,f,g,h) \
{ \
   a^=b<<11; d+=a; b+=c; \
   b^=c>>2;  e+=b; c+=d; \
   c^=d<<8;  f+=c; d+=e; \
   d^=e>>16; g+=d; e+=f; \
   e^=f<<10; h+=e; f+=g; \
   f^=g>>4;  a+=f; g+=h; \
   g^=h<<8;  b+=g; h+=a; \
   h^=a>>9;  c+=h; a+=b; \
}

static void randinit(void) {
   int i;
   uint32_t a,b,c,d,e,f,g,h;
   aa=bb=cc=0;
   a=b=c=d=e=f=g=h=0x9e3779b9;  /* the golden ratio */

   for (i=0; i<4; ++i)          /* scramble it */
   {
     mix(a,b,c,d,e,f,g,h);
   }

   for (i=0; i<256; i+=8)   /* fill in mm[] with messy stuff */
   {
     a+=randrsl[i  ]; b+=randrsl[i+1]; c+=randrsl[i+2]; d+=randrsl[i+3];
     e+=randrsl[i+4]; f+=randrsl[i+5]; g+=randrsl[i+6]; h+=randrsl[i+7];
     mix(a,b,c,d,e,f,g,h);
     mm[i  ]=a; mm[i+1]=b; mm[i+2]=c; mm[i+3]=d;
     mm[i+4]=e; mm[i+5]=f; mm[i+6]=g; mm[i+7]=h;
   }

   {        /* do a second pass to make all of the seed affect all of mm */
     for (i=0; i<256; i+=8)
     {
       a+=mm[i  ]; b+=mm[i+1]; c+=mm[i+2]; d+=mm[i+3];
       e+=mm[i+4]; f+=mm[i+5]; g+=mm[i+6]; h+=mm[i+7];
       mix(a,b,c,d,e,f,g,h);
       mm[i  ]=a; mm[i+1]=b; mm[i+2]=c; mm[i+3]=d;
       mm[i+4]=e; mm[i+5]=f; mm[i+6]=g; mm[i+7]=h;
     }
   }

   isaac();            /* fill in the first set of results */
   randcnt=256;        /* first use will run isaac() again */
}

/*****************************************************************************/

static void memtile(unsigned char* src, UV from, UV to) {
  while (from < to) {
    UV bytes = (2*from > to) ? to-from : from;
    memcpy(src+from, src, bytes);
    from += bytes;
  }
}

void isaac_seed(uint32_t bytes, const unsigned char* data)
{
  if (!mutex_init) {
    MUTEX_INIT(&state_mutex);
    mutex_init = 1;
  }
  MUTEX_LOCK(&state_mutex);
  memset(mm, 0, 4*256);
  memset(randrsl, 0, 4*256);
  if (bytes > 0 && data != 0) {
    unsigned char* rdata = (unsigned char*) randrsl;
    if (bytes > 1024) bytes = 1024;
    memcpy(rdata, data, bytes);
    /* Tile short seeds as recommended by Bob Jenkins. */
    memtile(rdata, bytes, 1024);
  }
  randinit();
  good_seed = (bytes >= 16);
  MUTEX_UNLOCK(&state_mutex);
}
int isaac_well_seeded(void) { return good_seed; }

/*****************************************************************************/

/* Unfortunately the mutexes are very expensive, even with only one thread */

void isaac_rand_bytes(uint32_t bytes, unsigned char* data)
{
  if (!mutex_init) croak("ISAAC used before init");
  MUTEX_LOCK(&state_mutex);
  if ( 4*(256-randcnt) >= bytes) {
    /* We have enough data, just copy it and leave */
    memcpy(data, (unsigned char*) (randrsl+randcnt), bytes);
    randcnt += (bytes+3)/4;
  } else {
    /* Loop copying up to 1024 bytes at a time */
    uint32_t n_rand_bytes, n_copy_bytes;
    while (bytes > 0) {
      if (randcnt > 255)
        isaac();
      n_rand_bytes = 4 * (256-randcnt);
      n_copy_bytes = (n_rand_bytes > bytes) ? bytes : n_rand_bytes;
      memcpy(data, (unsigned char*) (randrsl+randcnt), n_copy_bytes);
      data += n_copy_bytes;
      randcnt += (n_copy_bytes+3)/4;
      bytes -= n_copy_bytes;
    }
  }
  MUTEX_UNLOCK(&state_mutex);
}

uint32_t irand32(void)
{
#ifdef USE_ITHREADS
  uint32_t d;
  isaac_rand_bytes(4, (unsigned char*)&d);
  return d;
#else
  if (randcnt > 255)
    isaac();
  return randrsl[randcnt++];
#endif
}
UV irand64(void)
{
#if BITS_PER_WORD == 64
  UV d;
  isaac_rand_bytes(8, (unsigned char*)&d);
  return d;
#else
  croak("irand64 too many bits for UV");
#endif
}

static NV _tonv_32 = -1.0;
static NV _tonv_64 = -1.0;
NV drand64(void)
{
  if (_tonv_32 < 0 && mutex_init) {
    MUTEX_LOCK(&state_mutex);
    if (_tonv_32 < 0) {
      int i;
      NV t64, t32;
      for (t32 = 1.0, i = 0; i < 32; i++)
        t32 /= 2.0;
      for (t64 = t32, i = 0; i < 32; i++)
        t64 /= 2.0;
      _tonv_64 = t64;
      _tonv_32 = t32;
    }
    MUTEX_UNLOCK(&state_mutex);
  }
#if BITS_PER_WORD == 64
  return irand64() * _tonv_64;
#else
  return irand32() * _tonv_64 + irand32() * _tonv_32;
#endif
}

/*****************************************************************************/

/* Return rand 32-bit integer between 0 to n-1 inclusive */
uint32_t isaac_rand(uint32_t n)
{
  uint32_t r, rmax;
  if (n <= 1) return 0;
  if ((n & (n-1)) == 0) return (irand32() % n);
  rmax = (4294967295U / n) * n;
  do {
    r = irand32();
  } while (r >= rmax);
  return r % n;
}

UV isaac_rand64(UV n)
{
  if (n <= 1) return 0;
  if (n <= 4294967295U) {
    uint32_t r, rmax;
    if ((n & (n-1)) == 0) return (irand32() % n);
    rmax = (4294967295U / n) * n;
    do {
      r = irand32();
    } while (r >= rmax);
    return r % n;
  } else {
    UV r, rmax;
    if ((n & (n-1)) == 0) return (irand64() % n);
    rmax = (UV_MAX / n) * n;
    do {
      r = irand64();
    } while (r >= rmax);
    return r % n;
  }
}

UV irandb(int nbits)
{
  if (nbits == 0) {
    return 0;
  } else if (nbits <= 32) {
    return irand32() >> (32-nbits);
#if BITS_PER_WORD == 64
  } else if (nbits <= 64) {
    return irand64() >> (64-nbits);
#endif
  }
  croak("irand64 too many bits for UV");
}
