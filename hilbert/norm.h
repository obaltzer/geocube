#ifndef __NORM_H
#define __NORM_H

#include "hilbert.h"
#include "sort.h"

/**
 * Context for processing in a normalized context.
 */
struct fp_norm_context
{
    /** reference to the non-normalized context */
    struct fp_context* mc;
    /** reference to the actual computation context */
    struct fp_context* zc;

    /** array of lengths of mapping arrays (cardinality of dimensions) */
    size_t* map_size;
    /** array of pointers to mapping arrays */
    fpf_t** map;

    size_t dimf;
    size_t dimz;
};

struct fp_norm_context* fp_create_norm_context(int dimz, int dimf, 
                                               int base_order);
void fp_destroy_norm_context(struct fp_norm_context* context);
void fp_normalize(struct fp_norm_context* context, void* records, 
                  size_t n, void* result);
void fp_denormalize(struct fp_norm_context* context, void* norm_records, 
                    size_t n, void* denorm_records);

#endif
