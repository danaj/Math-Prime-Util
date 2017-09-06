#include <stdio.h>
#include "entropy.h"

/* A fallback timer entropy method that will probably never be used. */
#if defined(_WIN32_WCE)
static UV timer_entropy(UV bytes, unsigned char* buf) { return 0; }
#else
#include <time.h>
static uint32_t mix32(uint32_t r0) {   /* Similar to PCG 32 */
  uint32_t word = ((r0 >> ((r0 >> 28u) + 4u)) ^ r0) * 277803737u;
  return (word >> 22u) ^ word;
}
static uint32_t timer_mix8(uint32_t acc) {
  clock_t t1;
  uint32_t bit, a;
  for (bit = a = 0; bit < 8; bit++) {
    t1 = clock(); while (t1 == clock()) a ^= 1;
    acc = (acc << 1) | a;
  }
  return mix32(acc);
}
static UV timer_entropy(UV bytes, unsigned char* buf) {
  UV byte;
  uint32_t acc = 0;

  for (byte = 0; byte < 4; byte++)
    acc = timer_mix8(acc);
  for (byte = 0; byte < bytes; byte++) {
    acc = timer_mix8( timer_mix8( acc ) );
    buf[byte] = (acc >> 24) & 0xFF;
  }
  return bytes;
}
#endif

UV get_entropy_bytes(UV bytes, unsigned char* buf)
{
  UV len = 0;

#if defined(_WIN32) || defined(_WIN32_WCE)
#ifndef _WIN32_WINNT
  #define _WIN32_WINNT 0x0400
#endif
#ifdef _WIN32_WCE
   #define UNDER_CE
   #define ARM
#endif

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <wincrypt.h>

   /* TODO: Calling RtlGenRandom is faster */

   HCRYPTPROV hProv = 0;
   if (!CryptAcquireContext(&hProv, NULL, MS_DEF_PROV, PROV_RSA_FULL,
           CRYPT_SILENT | CRYPT_VERIFYCONTEXT | CRYPT_MACHINE_KEYSET) &&
       !CryptAcquireContext (&hProv, NULL, MS_DEF_PROV, PROV_RSA_FULL,
           CRYPT_SILENT | CRYPT_VERIFYCONTEXT | CRYPT_MACHINE_KEYSET | CRYPT_NEWKEYSET))
      return 0;

   if (CryptGenRandom(hProv, bytes, buf) == TRUE)
     len = bytes;
   CryptReleaseContext(hProv, 0);

#else  /*  ^^^^^^ Windows ^^^^^^          vvvvvv /dev/urandom vvvvvvv  */

  FILE *f = fopen("/dev/urandom", "rb");

  if (f == NULL)
    f = fopen("/dev/random", "rb");

  if (f != NULL) {
    if (setvbuf(f, NULL, _IONBF, 0) == 0) { /* disable buffering */
      len = (UV)fread(buf, 1, (size_t)bytes, f);
    }
    fclose(f);
  }

#endif

  /* Do a fallback method if something didn't work right. */
  if (len != bytes)
    len = timer_entropy(bytes, buf);

  return len;
}
