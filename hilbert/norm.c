#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "norm.h"
#include "hilbert.h"
#include "sort.h"
#include "quicksort.h"

#define BLOCK_SIZE 4096

struct fp_norm_context* fp_create_norm_context(int dimz, int dimf,
                                               int base_order)
{
    size_t i;
    struct fp_norm_context* c = 
        (struct fp_norm_context*)malloc(sizeof(struct fp_norm_context));
    assert(c != NULL);
    c->mc = fp_create_context(dimz, dimf, base_order);
    c->zc = fp_create_context(dimz + dimf, 0, base_order);
    c->map_size = (size_t*)malloc(sizeof(size_t) * dimf);
    c->map = (fpf_t**)malloc(sizeof(fpf_t*) * dimf);
    for(i = 0; i < dimf; i++)
        c->map[i] = NULL;
    c->dimf = dimf;
    c->dimz = dimz;
    return c;
}

void fp_destroy_norm_context(struct fp_norm_context* context)
{
    if(context)
    {
        size_t i;
        fp_destroy_context(context->mc);
        fp_destroy_context(context->zc);
        free(context->map_size);
        for(i = 0; i < context->dimf; i++)
            free(context->map[i]);
        free(context->map);
        
        free(context);
    }
}

int fp_compare_fpf(const void* context, const void* r1, const void* r2)
{
    return *((fpf_t*)r1) < *((fpf_t*)r2) ? -1 
            : *((fpf_t*)r1) == *((fpf_t*)r2) ? 0 : 1;
}

size_t fp_find_mapping(fpf_t* map, size_t n, fpf_t v)
{
    size_t hi = n - 1;
    size_t lo = 0;
    size_t pos = (hi - lo) / 2;
    fpf_t c = map[pos];
    
    while(c != v)
    {
        if(c > v)
            hi = pos;
        else
            lo = pos + 1;
        pos = lo + (hi - lo) / 2;
        c = map[pos];
    }
    return pos;
}
    
void fp_normalize(struct fp_norm_context* context, void* records,
                  size_t n, void* result)
{
    size_t i;
    size_t j;
    size_t counter;
    size_t n_blocks = 0;
    void* record_base;
    fpf_t block[BLOCK_SIZE];
    fpf_t* last;
    
    fpf_t* dims = (fpf_t*)malloc(sizeof(fpf_t) * context->dimf * n);
    
    /* copy continuous space attributes from raw records into dims array */
    for(i = 0; i < n; i++)
    {
        record_base = records + (i * context->mc->record_size);
        for(j = 0; j < context->dimf; j++)
            dims[j * n + i] = 
                ((fpf_t*)(record_base + context->mc->coordsf_off))[j];
    }

    /* sort the attributes of each dimensions and generate a map */
    for(i = 0; i < context->dimf; i++)
    {
        counter = 0;
        last = NULL;
        quicksort(NULL, &dims[i * n], n, sizeof(fpf_t), fp_compare_fpf);
    #if 0
        printf("DimF: %d\n", i);
        for(j = 0; j < n; j++)
        {
            printf("%f\n", dims[i * n + j]);
        }
    #endif
        for(j = 0; j < n; j++)
        {
            if(last == NULL || *last != dims[i * n + j])
            {
                block[counter++] = dims[i * n + j];
                last = &dims[i * n + j];
            }
            if(counter == BLOCK_SIZE)
            {
                counter = 0;
                context->map[i] = realloc(
                    context->map[i], 
                    (n_blocks + 1) * BLOCK_SIZE * sizeof(fpf_t)
                );
                memcpy(&context->map[i][n_blocks * BLOCK_SIZE],
                       block, sizeof(fpf_t) * BLOCK_SIZE);
                n_blocks++;
            }
        }
        if(counter != 0 && counter != BLOCK_SIZE)
        {
            context->map[i] = realloc(
                context->map[i], 
                (n_blocks * BLOCK_SIZE + counter) * sizeof(fpf_t)
            );
            memcpy(&context->map[i][n_blocks * BLOCK_SIZE],
                   block, sizeof(fpf_t) * counter);
        }
        context->map_size[i] = n_blocks * BLOCK_SIZE + counter;
    #if 0
        printf("DimF map: %d\n", i);
        for(j = 0; j < context->map_size[i]; j++)
        {
            printf("%d: %f\n", j, context->map[i][j]);
        }
    #endif
    }
    /* release temporary memory to compute the mapping */
    free(dims);
    
    /* for each record map the record into discrete space */
    for(i = 0; i < n; i++)
    {
        void* src_base = records + (i * context->mc->record_size);
        void* dst_base = result + (i * context->zc->record_size);
        fpz_t* dst_dimz = (fpz_t*)(dst_base + context->zc->coordsz_off);
        *((int*)(dst_base + context->zc->order_off)) = 
            *((int*)(src_base + context->mc->order_off));
        memcpy(dst_base + context->zc->coordsz_off, 
               src_base + context->mc->coordsz_off,
               sizeof(fpz_t) * context->dimz);
        for(j = 0; j < context->dimf; j++)
        {
            dst_dimz[context->dimz + j] = fp_find_mapping(
                context->map[j],
                context->map_size[j],
                ((fpf_t*)(src_base + context->mc->coordsf_off))[j]
            );
        }
    }

    /* set the new base_order for the normalized space */
    /* set the cardinalities of the new Z dimensions */
    for(i = 0; i < context->dimf; i++)
    {
        int order = 1;
        while((context->map_size[i] - 1) >> order)
            order++;
        if(context->zc->start_order < order)
            context->zc->start_order = order;
        context->zc->cards[context->dimz + i] = context->map_size[i];
    }
    context->zc->max_order = context->zc->start_order;
}

void fp_denormalize(struct fp_norm_context* norm, void* norm_records,
                    size_t n, void* denorm_records)
{
    size_t i;
    size_t j;
    void* in_record;
    void* out_record;
    
    for(i = 0; i < n; i++)
    {
        in_record = norm_records + (i * norm->zc->record_size);
        out_record = denorm_records + (i * norm->mc->record_size);

        /* copy discrete space dimensions */
        for(j = 0; j < norm->dimz; j++)
        {
            ((fpz_t*)(out_record + norm->mc->coordsz_off))[j] =
                ((fpz_t*)(in_record + norm->zc->coordsz_off))[j];
        }

        /* remap continuous space dimensions */
        for(j = 0; j < norm->dimf; j++)
        {
            ((fpf_t*)(out_record + norm->mc->coordsf_off))[j] =
                (fpf_t)norm->map[j][
                    (int)((fpz_t*)(in_record + norm->zc->coordsz_off))
                        [norm->dimz + j]
                ];
        }
    }
}
