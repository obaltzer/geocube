#include "sort.h"
#include "quicksort.h"
#include "hilbert.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

void print_record(struct fp_context* context, const void* r, int k)
{
    int i;
    fpz_t coordz;
    fpf_t coordf;
    fpz_t bits = 1;
    bits = bits << k;
    
    printf("(%d, %llu: ", k, bits);
    for(i = 0; i < context->env->dimz; i++)
    {
        coordz = ((fpz_t*)(r + context->offset_coordsz))[i];
        if(i)
            printf(", %llu", coordz);
        else
            printf("%llu", coordz);
    }
    for(i = 0; i < context->env->dimf; i++)
    {
        coordf = ((fpf_t*)(r + context->offset_coordsf))[i];
        if((context->env->dimz + i))
            printf(", %f:%llu", coordf, (fpz_t)(coordf * (fpf_t)bits));
        else
            printf("%f:%llu", coordf, (fpz_t)(coordf * (fpf_t)bits));
    }
    printf(")");
}
int fp_compute_order(const struct fp_context* context, 
                     const void* r1, const void* r2)
{
    fpf_t sum = 0.0;
    int i;
    fpz_t intervals;
    fpz_t tmp;
    int o;
    
    fpf_t* r1f = (fpf_t*)(r1 + context->offset_coordsf);
    fpf_t* r2f = (fpf_t*)(r2 + context->offset_coordsf);
    fpf_t diff;
    fpf_t dist;
    
    /* compute the eucledian distance between the continiuous parts of the
     * points */
    for(i = 0; i < context->env->dimf; i++)
    {
        diff = r1f[i] - r2f[i];
        sum += diff * diff;
    }
    dist = sqrt(sum);
#ifdef WITH_PRINT
    printf("Distance: %f\n" , dist);
#endif
    /* compute number of intervals each having a diagonal of at most dist 
     * length */
    tmp = intervals = 
        (fpz_t)ceil(1.0 / sqrt(dist * dist / context->env->dimf)) ;
    
    /* compute the order that generates the smaller power of 2
     * intervals */
    o = 0;
    while(intervals >> o)
        o++;
    
    /* if intervals is still large we have to increment the power by on */
    if(intervals > ((fpz_t)1 << o))
        o++;
    /* the order should now generate a grid that is smaller then the
     * distance bewteen the two points */
    /* assert(1.0 / (fpf_t)((fpz_t)1 << o) < dist); */
#ifdef WITH_PRINT
    printf("Order: %d\n", o);
#endif
    return o;
}

int fp_find_order_constant(const struct fp_context* context, 
                           const void* r1, const void* r2,
                           int order)
{
    fpz_t* index1 = (fpz_t*)(r1 + sizeof(int));
    fpz_t* index2 = (fpz_t*)(r2 + sizeof(int));
    int new_order = order;
    
    /* If the order is uninitialized compute the Hilbert indices of the
     * records at the minimum order. */
    if(order == -1)
    {
        new_order = fp_compute_order(context, r1, r2);
        /* if the new order is smaller than the minimum order, use the
         * minimum order, which should also guarantee that the two points
         * are not sharing the same cell */
        if(new_order < context->start_order)
            new_order = context->start_order;
        fpm_c2i(context->env, new_order, 
                (fpz_t*)(r1 + context->offset_coordsz),
                (fpf_t*)(r1 + context->offset_coordsf),
                index1);
        fpm_c2i(context->env, new_order, 
                (fpz_t*)(r2 + context->offset_coordsz),
                (fpf_t*)(r2 + context->offset_coordsf),
                index2);
    }
    /* If we share the same order as well as the same indices, try to
     * compute the order via Eucledian distance. This assumes that the
     * discrete space coordinates are equal, which they have to be
     * otherwise the records would not share the same index. */
    else if(*index1 == *index2)
    {
        new_order = fp_compute_order(context, r1, r2);
        if(new_order > order)
        {
            /* make sure we do not exceed the maximum order */
            assert(new_order <= context->order_limit);
            fpm_c2i(context->env, new_order, 
                    (fpz_t*)(r1 + context->offset_coordsz),
                    (fpf_t*)(r1 + context->offset_coordsf),
                    index1);
            fpm_c2i(context->env, new_order, 
                    (fpz_t*)(r2 + context->offset_coordsz),
                    (fpf_t*)(r2 + context->offset_coordsf),
                    index2);
        }
        else
            new_order = order;
    }
    return new_order;
}

int fp_find_order_iterative(const struct fp_context* context, 
                            const void* r1, const void* r2,
                            int order)
{
    fpz_t* index1 = (fpz_t*)(r1 + sizeof(int));
    fpz_t* index2 = (fpz_t*)(r2 + sizeof(int));
   
    /* now that we know the both orders are the same we can compare both
     * records indices, if the order is -1 the indices are uninitialized
     * and we have to recompute in any case */
    while(order == -1 
        || (*index1 == *index2 && order <= context->order_limit))
    {
        if(order == -1)
            order = context->start_order;
        else
            order++;
        fpm_c2i(context->env, order, 
                (fpz_t*)(r1 + context->offset_coordsz),
                (fpf_t*)(r1 + context->offset_coordsf),
                index1);
        fpm_c2i(context->env, order, 
                (fpz_t*)(r2 + context->offset_coordsz),
                (fpf_t*)(r2 + context->offset_coordsf),
                index2);
    }
    return order;
}

int fp_compare(const void* c, const void* r1, const void* r2)
{
    struct fp_context* context = (struct fp_context*)c;
    int i;
    fpz_t* index1;
    fpz_t* index2;
    int* order1;
    int* order2;
    int order;
    int new_order;
    int done = 0;

    /* return equality of both pointers point to the same record */
    if(r1 == r2)
        return 0;

    /* check if the two records have equal values */
    for(i = 0; !done && i < context->env->dimz; i++)
        if(((fpz_t*)(r1 + context->offset_coordsz))[i] 
                != ((fpz_t*)(r2 + context->offset_coordsz))[i])
            done = 1;
    for(i = 0; !done && i < context->env->dimf; i++)
        if(((fpf_t*)(r1 + context->offset_coordsf))[i] 
                != ((fpf_t*)(r2 + context->offset_coordsf))[i])
            done = 1;
    
    if(!done)
        /* Report equality if they are in fact equal. */
        return 0;
    
    /* get the pointers to the current indices of the records */
    index1 = (fpz_t*)(r1 + sizeof(int));
    index2 = (fpz_t*)(r2 + sizeof(int));
    order1 = (int*)(r1);
    order2 = (int*)(r2);
    
    /* Determine which of the two record's has been computed with the
     * higher order. If both records have been computed at the same order,
     * all we need to do is to compare their Hilbert indices and return if
     * they differ. If they do not differ, we have to increase the order
     * and recompute. If the two records order's differ, it is necessary to
     * recompute the Hilbert index of the record with the smaller order and
     * reduce the problem to such that both records have been computed at
     * the same order. */
    /* fist make sure that both records are computed at the same order */
    if(*order1 < *order2)
    {
        /* If r1's order is smaller than r2's order, recompute r1's index
         * and set its order to the one of r2. We do not have to test for
         * order1 == -1, because we are recomputing index1 anyway and
         * index2 cannot be uninitialized, otherwise order2 would not be
         * larger then order1. */
        fpm_c2i(context->env, *order2, 
                (fpz_t*)(r1 + context->offset_coordsz), 
                (fpf_t*)(r1 + context->offset_coordsf), 
                index1);
        order = *order1 = *order2;
    }
    else if(*order1 > *order2)
    {
        /* If r2's order is smaller than r1's order, recompute r2's index
         * and set its order to the one of r1. We do not have to test for
         * order2 == -1, because we are recomputing index2 anyway and
         * index1 cannot be uninitialized, otherwise order1 would not be
         * larger then order2. */
        fpm_c2i(context->env, *order1, 
                (fpz_t*)(r2 + context->offset_coordsz),  
                (fpf_t*)(r2 + context->offset_coordsf), 
                index2);
        order = *order2 = *order1;
    }
    else
        order = *order1 = *order2;

    if(order == -1 || *index1 == *index2)
        order = fp_find_order_iterative(context, r1, r2, order);
        
    /* the order limit is the bit length of index datatype divided by the
     * number of dimensions */
    if(order > context->order_limit)
    {
        printf("Indices: %llu : %llu\n", *index1, *index2);
        printf("Order limit exceeded for: ");
        print_record(context, r1, order - 1);
        printf(" and ");
        print_record(context, r2, order - 1);
        printf("\n");
        abort();
    }
    /* set the new maximum order of we have one */
    if(context->max_order < order)
        context->max_order = order;
        
    *order1 = *order2 = order;

    assert(*index1 != *index2);
    if(*index1 < *index2)
        return -1;
    else if(*index1 > *index2)
        return 1;
    return 0;
}

void mergesort(const void* context, void* pbase,
                size_t total_elems, size_t size,
                int (*cmp)(const void*, const void*, const void*))
{
    if(total_elems >= 2)
    {
        void* newBuf1 = malloc(size * (total_elems / 2));
        void* newBuf2 = malloc(size * ((total_elems + 1) / 2));
        int i = 0;
        int j = 0;
        
        memcpy(newBuf1, pbase, size * (total_elems / 2));
        memcpy(newBuf2, pbase + (size * (total_elems / 2)),
                size * ((total_elems + 1) / 2));
        
        mergesort(context, newBuf1, total_elems / 2, size, cmp);
        mergesort(context, newBuf2, (total_elems + 1) / 2, size, cmp);

        while(i < total_elems / 2 && j < (total_elems + 1) / 2)
            if(cmp(context, newBuf1 + (size * i), newBuf2 + (size * j)) == -1)
            {
                memcpy((void*)(pbase + ((i + j) * size)), 
                        (void*)(newBuf1 + (size * i)), size);
                i++;
            }
            else
            {
                memcpy((void*)(pbase + ((i + j) * size)), 
                        (void*)(newBuf2 + (size * j)), size);
                j++;
            }

        if(i == total_elems / 2)
            while(j < (total_elems + 1) / 2)
            {
                memcpy((void*)(pbase + ((i + j) * size)), 
                        (void*)(newBuf2 + size * j), size);
                j++;
            }
        else if(j == (total_elems + 1) / 2)
            while(i < total_elems / 2)
            {
                memcpy((void*)(pbase + ((i + j) * size)), 
                        (void*)(newBuf1 + size * i), size);
                i++;
            }
        
        free(newBuf1);
        free(newBuf2);
    }
    return;
}
        
    
void fp_im_sort(struct fp_context* context, void* input,  size_t n,
                void** output)
{
    /* output array has to be NULL, because we do everything in place */
    assert(*output == NULL);
    quicksort(context, input, n, context->sizeof_record, fp_compare);
    *output = input;
}

struct fp_context* fp_create_context(int dimz, int dimf, int base_order)
{
    struct fp_context* c = 
        (struct fp_context*)malloc(sizeof(struct fp_context));

    if(c)
    {
        /* create a new mixed dimensions type */
        c->env = fpm_create_env(dimz, dimf, base_order);
        
        /* The format of a record is as follows:
         *
         * int              : order at which the current index was computed
         * fpz_t            : the current Hilbert index of the record
         * dimz * fpz_t     : discrete space dimensions
         * dimf * fpf_t     : continuous space dimensions
         */
        c->sizeof_record = sizeof(int) + sizeof(fpz_t) 
                            + dimz * sizeof(fpz_t) 
                            + dimf * sizeof(fpf_t);
        c->offset_coordsz = sizeof(int) + sizeof(fpz_t);
        c->offset_coordsf = sizeof(int) + sizeof(fpz_t)
                            + dimz * sizeof(fpz_t);
        
        /* compute the memory size of a bounding box 
         * consisting of 2 points */
        c->bbox_size = 2 * (dimz * sizeof(fpz_t) + dimf * sizeof(fpf_t));
        /* compute the offsets for each of the points */
        c->minz = 0;
        c->minf = dimz * sizeof(fpz_t);
        c->maxz = dimz * sizeof(fpz_t) + dimf * sizeof(fpf_t);
        c->maxf = dimz * sizeof(fpz_t) + dimf * sizeof(fpf_t)
                            + dimz * sizeof(fpz_t);
        /* TODO: hardcoded fanout of the tree */
        c->fanout = 16;
        /* the order limit is the bit length of index datatype divided by
         * the number of dimensions */
        c->order_limit = (sizeof(fpz_t) * 8) / (dimf + dimz);
        /* start order, which is used for unitialized records, is equal to
         * the base order of the dicrete space */
        c->start_order = base_order;
        /* the maximum order used is at least that of the base order */
        c->max_order = base_order;

        /* Profiler initialize */
        c->build_tree_calls = 0;
        c->n_tree_nodes = 0;
    }
    return c;
}

void fp_destroy_context(struct fp_context* context)
{
    if(context)
    {
        fpm_destroy_env(context->env);
        free(context);
    }
}
