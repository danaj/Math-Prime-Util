#ifndef MPU_THREADLOCK_H
#define MPU_THREADLOCK_H

#include "ptypes.h"

#ifdef STANDALONE
  #undef USE_ITHREADS
  #define MUTEX_INIT(x)
  #define MUTEX_LOCK(x)
  #define MUTEX_UNLOCK(x)
  #define MUTEX_DESTROY(x)
  #define COND_INIT(x)
  #define COND_DESTROY(x)
#endif

#ifndef USE_ITHREADS

  #define MUTEX_DECL(x)
  #define READ_WRITE_LOCK_DECL(x)
  #define WRITE_LOCK_START(x)
  #define WRITE_LOCK_END(x)
  #define READ_LOCK_START(x)
  #define READ_LOCK_END(x)

#else

  #define MUTEX_DECL(x) \
    static perl_mutex x ## _mutex;

  #define READ_WRITE_LOCK_DECL(x) \
    static perl_mutex x ## _mutex;   \
    static perl_cond  x ## _turn;    \
    static int        x ## _reading; \
    static int        x ## _writing; \
    static int        x ## _writers;

  #define WRITE_LOCK_START(x) \
   do { \
     MUTEX_LOCK(& x ## _mutex); \
     x ## _writers++; \
     while (x ## _reading || x ## _writing) \
       COND_WAIT(& x ## _turn, & x ## _mutex); \
     x ## _writing++; \
     MUTEX_UNLOCK(& x ## _mutex); \
   } while (0)

 #define WRITE_LOCK_END(x) \
   do { \
     MUTEX_LOCK(& x ## _mutex); \
     x ## _writing--; \
     x ## _writers--; \
     COND_BROADCAST(& x ## _turn); \
     MUTEX_UNLOCK(& x ## _mutex); \
   } while (0)

 #define READ_LOCK_START(x) \
   do { \
     MUTEX_LOCK(& x ## _mutex); \
     if (x ## _writers) \
       COND_WAIT(& x ## _turn, & x ## _mutex); \
     while (x ## _writing) \
       COND_WAIT(& x ## _turn, & x ## _mutex); \
     x ## _reading++; \
     MUTEX_UNLOCK(& x ## _mutex); \
   } while (0)

 #define READ_LOCK_END(x) \
   do { \
     MUTEX_LOCK(& x ## _mutex); \
     x ## _reading--; \
     COND_BROADCAST(& x ## _turn); \
     MUTEX_UNLOCK(& x ## _mutex); \
   } while (0)
#endif

#endif
