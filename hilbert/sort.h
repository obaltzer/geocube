/**
 * @file sort.h
 *
 * Defines a set of functions to sort a set of records in their Hilbert
 * order while dynamically adjusting the order of the underlying Hilbert
 * curve.
 */
#ifndef __SORT_H
#define __SORT_H

#include "hilbert.h"

struct fp_context
{
    struct fpm_env* env;
    size_t sizeof_record;
    size_t offset_coordsz;
    size_t offset_coordsf;
    size_t bbox_size;
    size_t minz;
    size_t minf;
    size_t maxz;
    size_t maxf;
    size_t fanout;
    int start_order;
    int max_order;
    int order_limit;

    int (*find_order)(const struct fp_context*,
                      const void*, const void*, int);

    /* profiling code */
    unsigned long long int build_tree_calls;
    unsigned int n_tree_nodes;
};

struct fp_context* fp_create_context(int dimz, int dimf, int base_order);
void fp_destroy_context(struct fp_context* context);
void fp_im_sort(struct fp_context* context, void* input, size_t n, 
                void** output);
void print_record(struct fp_context* context, const void* r, int k);
int fp_find_order_iterative(const struct fp_context* context, 
                            const void* r1, const void* r2,
                            int order);
int fp_find_order_constant(const struct fp_context* context, 
                           const void* r1, const void* r2,
                           int order);

#endif
