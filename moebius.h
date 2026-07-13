#ifndef MPU_MOEBIUS_H
#define MPU_MOEBIUS_H

#include "ptypes.h"

extern int          moebius(UV n);
extern signed char* range_moebius(UV lo, UV hi);
extern void*        start_segment_moebius(UV lo, UV hi, signed char** segmentmem);
extern bool         next_segment_moebius(void* vctx, UV* lo, UV* hi);
extern void         end_segment_moebius(void* vctx);
extern IV           mertens(UV n);
extern IV           mertens_range(UV lo, UV hi);
extern int          liouville(UV n);
extern IV           sumliouville(UV n);
extern signed char* range_liouville(UV lo, UV hi);

#endif
