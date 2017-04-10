/*
 * The ISAAC CSPRNG plus interface.
 * Slightly modified readable.c from Bob Jenkins 1996.
 */

#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include "ptypes.h"
#include "isaac.h"


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

void isaac_seed(uint32_t bytes, const unsigned char* data)
{
  memset(mm, 0, 4*256);
  memset(randrsl, 0, 4*256);
  if (bytes > 0 && data != 0) {
    unsigned char* rdata = (unsigned char*) randrsl;
    if (bytes > 1024) bytes = 1024;
    memcpy(rdata, data, bytes);
  }
  randinit();
}

void isaac_rand_bytes(uint32_t bytes, unsigned char* data)
{
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
}

uint32_t isaac_irand32(void)
{
  uint32_t a;
  if (randcnt > 255) isaac();
  return randrsl[randcnt++];
}
UV isaac_irand64(void)
{
  uint32_t a, b;
  if (randcnt > 255) isaac();
  a = randrsl[randcnt++];
  if (randcnt > 255) isaac();
  b = randrsl[randcnt++];
  return (((UV)b) << 32) | a;
}
