#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define FUNC_log2floor 1
#include "util.h"
#include "sort.h"

#define USE_QUADSORT 0

/*
 * Sorting arrays of integers.
 *
 * We really have two quite different use cases.  The first is for internal
 * use, where it's very common to have a small number of inputs.
 * E.g. sorting roots, factors, divisors, totients, prime powers, etc.  Most
 * of these will have small arrays and the overall time is dominated by the
 * real work that created the data.
 * There are some degenerate cases that generate many inputs, but these are
 * exceptional.  Most sorts from our test suite are 32 or fewer items, with
 * the largest being 576 items.
 *
 * The second use is from vecsort, where the user or our PP code has given
 * us a possibly large array to sort.  Here we have the additional challenge
 * of making sure the overhead of Perl->C->Perl is as small as possible.
 *
 *  We have a number of possible choices.
 *
 * 1) Perl's sort.  A cache-aware merge sort, which makes a lot of sense for
 *    its use with arbitrary and complicated data structures, possibly
 *    expensive comparisons, and where a stable sort is highly desirable.
 *    Most of this is irrelevant for sorting simple integers.
 *    Problem 1: We can sort SV's but there isn't a simple UV interface.
 *    Problem 2: It's slow for shuffled inputs, like most stable merge sorts.
 *
 * 2) qsort.  Easy and works, but system dependent.
 *    Can be quite fast -- MacOS/clang is 3x faster than merge sort for
 *    shuffled inputs, and has fast behavior with sorted/reversed data.
 *
 * 3) fluxsort/quadsort/timsort/powersort/glidesort/etc.
 *    fluxsort is extremely fast and has excellent behavior with ordered data.
 *    The main reason it isn't being used here is the code size.
 *
 * 4) insertion / Shell.  Fastest on tiny arrays and very compact.  We use
 *    insertion sort for small chunks.
 *
 * 5) heapsort.  lobby99's implementation here is surprisingly fast and very
 *    consistent across a variety of inputs.  It is used as a fallback if
 *    quicksort is choosing bad partitions.
 *
 * 6) quicksort.  Yes, yet another quicksort implementation.  Fast for small
 *    inputs, competitive for larger.  This uses pseudo-median of 9
 *    partitioning, insertion sort for small partitions, and will switch to
 *    heapsort after enough bad partitions, so there is no O(n^2) disaster.
 *
 * 5) radix sort.  With enough integers, radix sort beats everything else on
 *    shuffled data.  Performance on ordered data is decent though not like
 *    fluxsort.  Uses auxiliary data equal to the input size.
 *
 * We use our quicksort for small arrays, radixsort for larger.
 *
 */


/******************************************************************************/

static void insertionsort_uv(UV *array, size_t len) {
  size_t i,j;
  for (i = 1; i < len; i++) {
    UV t = array[i];
    for (j = i; j > 0 && array[j-1] > t; j--)
      array[j] = array[j-1];
    array[j] = t;
  }
}
static void insertionsort_iv(IV *array, size_t len) {
  size_t i,j;
  for (i = 1; i < len; i++) {
    IV t = array[i];
    for (j = i; j > 0 && array[j-1] > t; j--)
      array[j] = array[j-1];
    array[j] = t;
  }
}

#if 0
static void shellsort_uv(UV *array, size_t len) {
  static unsigned short sgaps[] = {209,109,41,19,5,1}; /* Sedgewick 1986 */
  size_t i, j, gap, gi = 0;
  do {
    gap = sgaps[gi++];
    for (i = gap; i < len; i++) {
      UV t = array[i];
      for (j = i; j >= gap && array[j-gap] > t; j -= gap)
        array[j] = array[j-gap];
      array[j] = t;
    }
  } while (gap > 1);
}
static void shellsort_iv(IV *array, size_t len) {
  static unsigned short sgaps[] = {209,109,41,19,5,1}; /* Sedgewick 1986 */
  size_t i, j, gap, gi = 0;
  do {
    gap = sgaps[gi++];
    for (i = gap; i < len; i++) {
      IV t = array[i];
      for (j = i; j >= gap && array[j-gap] > t; j -= gap)
        array[j] = array[j-gap];
      array[j] = t;
    }
  } while (gap > 1);
}
#endif

/******************************************************************************/
/*                             RADIX SORT                                     */
/******************************************************************************/
#define RADIX_BIT 8
#define RADIX     (1u<<RADIX_BIT)
static bool _radixsort(UV *array, size_t n, bool is_iv)
{
  size_t i, count[RADIX];
  unsigned r;
  UV *a, *b, *ptr[RADIX];
  UV passmask = 0;
  int sh;

  memset(count, 0, sizeof count);
  for (i = 0; i < n; i++) {
    UV d = array[i];
    passmask |= d ^ (d >> RADIX_BIT);
    count[d % RADIX]++;
  }
  if (passmask < RADIX) { /* If all values < RADIX, Use *fast* counting sort */
    if (passmask) {
      size_t j = 0, lim = 0;
      for (r = 0; r < RADIX; r++)
        for (lim += count[r]; j < lim; j++)
          array[j] = r;
    }
    return 1;
  }
  /* Allocate second ping-pong buffer */
  a = array;
  b = malloc(n * sizeof(UV));
  if (b == 0) return 0;
  /* Each pass radix-sorts and counts for next pass */
  for (sh = 0; UV_MAX >> sh >= RADIX; sh += RADIX_BIT) {
    UV *p = b;
    if ((passmask >> sh) % RADIX == 0)
      continue;
    for (r = 0; r < RADIX; r++) {
      ptr[r] = p;
      p += count[r];
    }
    assert(p == b + n);
    memset(count, 0, sizeof count);
    for (i = 0; i < n; i++) {
      UV d = a[i];
      *(ptr[(d>>sh) % RADIX]++) = d;
      count[(d >> (sh + RADIX_BIT)) % RADIX]++;
    }
    p = b; b = a; a = p;
  }
  /* Last pass does no more counting */
  if (passmask >> sh) {
    UV *p = b;
    unsigned signbit = is_iv ? 1 << (BITS_PER_WORD-1)%RADIX_BIT : 0;
    for (r = 0; r < RADIX; r++) {
      ptr[r^signbit] = p;
      p += count[r^signbit];
    }
    assert(p == b + n);
    for (i = 0; i < n; i++) {
      UV d = a[i];
      *(ptr[(d>>sh) % RADIX]++) = d;
    }
    p = b; b = a; a = p;
  }
  /* Move back to input array if necessary */
  if (a != array) {
    memcpy(array, a, n * sizeof *array);
    b = a;
  }
  free(b);
  return 1;
}
#undef RADIX_BIT
#undef RADIX

/******************************************************************************/
/*                              HEAP SORT                                     */
/******************************************************************************/
static void _heapsort(UV *array, size_t len, bool is_iv)
{
  size_t a = len/2;

  if (!a)               /* Trivial cases: len < 2 */
    return;

  for (len--;;) {
    UV r;               /* Value from array[a] being sifted down */
    size_t b, c;        /* Current descendent and its child */

    /*
     * Elements [0,a) are unsorted.
     * Elements [a,n] are in the heap.
     * Elements (n,...) are sorted.
     */
    if (a > 0)          /* Building heap: sift down array[--a] */
      r = array[--a];
    else if (len > 0) { /* Extracting: Swap root<->array[n--] */
      r = array[len]; array[len--] = array[0];
    } else              /* Extraction complete */
      return;

    /* Sift element r (at "a") down into heap. */
    if (!is_iv) {
      for (b = a; (c = 2*b + 1) < len; b = c) {
        UV s = array[c];
        if (array[c+1] >= s)
          s = array[++c];
        if (r >= s)
          goto sift_done;
        array[b] = s;
      }
    } else {
      for (b = a; (c = 2*b + 1) < len; b = c) {
        IV s = array[c];
        if ((IV)array[c+1] >= s)
          s = array[++c];
        if ((IV)r >= s)
          goto sift_done;
        array[b] = s;
      }
    }
    if (c == len) { /* Corner case: last leaf with no sibling */
      if ( (!is_iv && r < array[c]) || (is_iv && (IV)r < (IV)array[c]) ) {
        array[b] = array[c];
        b = c;
      }
    }
sift_done:
    array[b] = r;
  }
}

#define radixsort_uv(L,len)  _radixsort(L, len, 0)
#define radixsort_iv(L,len)  _radixsort((UV*)L, len, 1)
#define heapsort_uv(L,len)   _heapsort(L, len, 0)
#define heapsort_iv(L,len)   _heapsort((UV*)L, len, 1)

/******************************************************************************/
/*                              QUICK SORT                                    */
/******************************************************************************/

static size_t _mid3_uv_index(UV* L, size_t a, size_t b, size_t c) {
  const UV s[3] = {a,b,c};  /* Scandum's branchless method */
  int x = L[a] > L[b];
  int y = L[a] > L[c];
  int z = L[b] > L[c];
  return s[(x == y) + (y ^ z)];
}
static size_t _mid3_iv_index(IV* L, size_t a, size_t b, size_t c) {
  const IV s[3] = {a,b,c};  /* Scandum's branchless method */
  int x = L[a] > L[b];
  int y = L[a] > L[c];
  int z = L[b] > L[c];
  return s[(x == y) + (y ^ z)];
}

static size_t _partition_uv(UV* L, size_t lo, size_t hi) {
  size_t i = lo-1, j = hi+1, len = hi-lo+1;
  UV pivot;
  if (len <= 7) {
    pivot = L[len/2];
  } else if (len <= 40) {
    pivot = L[_mid3_uv_index(L, lo, lo+(hi-lo)/2, hi)];
  } else {
    size_t s = len/8;
    size_t i1 = _mid3_uv_index(L, lo, lo+s, lo+2*s);
    size_t i2 = _mid3_uv_index(L, lo+3*s, lo+(hi-lo)/2, hi-3*s);
    size_t i3 = _mid3_uv_index(L, hi-2*s, hi-s, hi);
    pivot = L[_mid3_uv_index(L, i1, i2, i3)];
  }
  while (1) {
    do { i++; } while (L[i] < pivot);
    do { j--; } while (L[j] > pivot);
    if (i >= j) return j;
    { UV t = L[i];  L[i] = L[j];  L[j] = t; }
  }
}
static size_t _partition_iv(IV* L, size_t lo, size_t hi) {
  size_t i = lo-1, j = hi+1, len = hi-lo+1;
  IV pivot;
  if (len <= 7) {
    pivot = L[len/2];
  } else if (len <= 40) {
    pivot = L[_mid3_iv_index(L, lo, lo+(hi-lo)/2, hi)];
  } else {
    size_t s = len/8;
    size_t i1 = _mid3_iv_index(L, lo, lo+s, lo+2*s);
    size_t i2 = _mid3_iv_index(L, lo+3*s, lo+(hi-lo)/2, hi-3*s);
    size_t i3 = _mid3_iv_index(L, hi-2*s, hi-s, hi);
    pivot = L[_mid3_iv_index(L, i1, i2, i3)];
  }
  while (1) {
    do { i++; } while (L[i] < pivot);
    do { j--; } while (L[j] > pivot);
    if (i >= j) return j;
    { IV t = L[i];  L[i] = L[j];  L[j] = t; }
  }
}

static void _qs_uv(UV* L, size_t lo, size_t hi, int badpartsleft) {
  size_t p, size = hi-lo+1;

  if (size <= 16)
    return insertionsort_uv(L+lo, size);

  p = _partition_uv(L, lo, hi);

  { /* check for unbalanced partitions, same as pdqsort */
    size_t l_size = p - lo + 1;
    size_t r_size = hi - (p+1) + 1;
    bool highly_unbalanced = l_size < size / 8 || r_size < size / 8;
    if (highly_unbalanced && --badpartsleft <= 0)
      return heapsort_uv(L+lo, size);
  }

  _qs_uv(L, lo, p, badpartsleft);
  _qs_uv(L, p+1, hi, badpartsleft);
}
static void _qs_iv(IV* L, size_t lo, size_t hi, int badpartsleft) {
  size_t p, size = hi-lo+1;

  if (size <= 16)
    return insertionsort_iv(L+lo, size);

  p = _partition_iv(L, lo, hi);

  { /* check for unbalanced partitions, same as pdqsort */
    size_t l_size = p - lo + 1;
    size_t r_size = hi - (p+1) + 1;
    bool highly_unbalanced = l_size < size / 8 || r_size < size / 8;
    if (highly_unbalanced && --badpartsleft <= 0)
      return heapsort_iv(L+lo, size);
  }

  _qs_iv(L, lo, p, badpartsleft);
  _qs_iv(L, p+1, hi, badpartsleft);
}

static void quicksort_uv(UV *L, size_t len) {
  if (len > 1) _qs_uv(L, 0, len-1, log2floor(len));
}
static void quicksort_iv(IV *L, size_t len) {
  if (len > 1) _qs_iv(L, 0, len-1, log2floor(len));
}


#if USE_QUADSORT

#include "quadsortuv.h"
void sort_uv_array(UV* L, size_t len) { quadsort_uv(L, len, 0); }
void sort_iv_array(IV* L, size_t len) { quadsort_iv(L, len, 0); }

#else

void sort_uv_array(UV* L, size_t len)
{
  UV *aux;
  if (len < 800) {
    quicksort_uv(L, len);
  } else {
    /* We could use an in-place radix sort like Ska Sort.  Our radix sort
     * is traditional and uses O(n) extra memory.  If we cannot get the
     * extra memory, we fall back to an in-place sort. */
    if (!radixsort_uv(L, len))
      quicksort_uv(L, len);
  }
}

void sort_iv_array(IV* L, size_t len)
{
  IV *aux;
  if (len < 800) {
    quicksort_iv(L, len);
  } else {
    if (!radixsort_iv(L, len))   /* radixsort could fail aux allocation */
      quicksort_iv(L, len);
  }
}

#endif

/******************************************************************************/

void sort_dedup_uv_array(UV* L, bool data_is_signed, size_t *len)
{
  if (*len > 1) {
    size_t i, j;
    if (data_is_signed)  sort_iv_array((IV *)L, *len);
    else                 sort_uv_array(L, *len);
    for (i=0, j=1; j < *len; j++) {
      i += (L[i] != L[j]);
      L[i] = L[j];
    }
    *len = i+1;
  }
}

