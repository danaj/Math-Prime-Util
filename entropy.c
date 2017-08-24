#include <stdio.h>
#include "entropy.h"

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

  if (len == bytes)  return len;

  /* Something didn't work.  Do a fallback method. */

  /* TODO: Something better here.
   *   1. Get a decent seed, maybe some Perl_seed values
   *   2. Get a private CSPRNG context, seed with above
   *   3. Full buf from csprng
   *   4. destroy csprng context
   */
  while (len < bytes) {
    uint32_t i, s = Perl_seed();   /* TODO Perl 5.6 */
    for (i = 0; i < 4 && len < bytes; i++)
      buf[len++] = (s >> (8*i)) & 0xFF;
  }
  return len;
}
