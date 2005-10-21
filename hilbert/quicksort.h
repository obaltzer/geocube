#ifndef __QUICKSORT_H
#define __QUICKSORT_H

void quicksort(void *const context, void *const pbase,
               size_t total_elems, size_t size,
               int (*cmp)(const void*, const void*, const void*));

#endif
