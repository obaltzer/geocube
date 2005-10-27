#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <sys/time.h>
#include <time.h>

#include "hilbert.h"
#include "sort.h"

int main(int argc, char** argv)
{
    FILE* f = NULL;
    FILE* of = NULL;
    int dimz;
    int dimf;
    int order;
    int i, j;
    int n;
    int coordz;
    size_t* cards;
    size_t max_card = 0;
    float coordf;
    void* record_base;
    void* results = NULL;
    struct timeval start_time;
    struct timeval end_time;
    double runtime;
#ifdef WITH_FP
    void* records = NULL;
    fpz_t index;
    fpz_t last_index = 0;
    int failures = 0;
    struct fp_context* context;
#else
    void* records = NULL;
    mpz_t index;
    struct mp_context* context;
#endif
    
    if(argc < 2)
    {
        printf("Usage: sort_test datafile [outfile]\n");
        exit(2);
    }
    
    if((f = fopen(argv[1], "r")) == NULL)
    {
        perror("Unable to open datafile:");
        exit(3);
    }

    if(argc == 3)
    {
        if((of = fopen(argv[2], "w")) == NULL)
        {
            perror("Unable to open output file:");
            exit(4);
        }
    }
    
    /* read number of discrete space dimensions */
    fscanf(f, "%d", &dimz);
    /* read the cardinalities of the dimensions and determine 
     * the largest */
    cards = (size_t*)malloc(sizeof(size_t) * dimz);
    for(i = 0; i < dimz; i++)
    {
        fscanf(f, "%d", &cards[i]);
        if(cards[i] > max_card)
            max_card = cards[i];
    }
    order = 1;
    while(max_card >> 1)
    {
        max_card >>= 1;
        order++;
    }
    printf("Using start order: %d\n", order);
    /* read number of continuous space dimensions */
    fscanf(f, "%d\n", &dimf);
    /* number of records */
    fscanf(f, "%d\n", &n);
    
#ifdef WITH_FP    
    /* create a new context */
    context = fp_create_context(dimz, dimf, order);
    /* allocate memory for the records */
    records = malloc(n * context->sizeof_record);
#else
    /* create the new context */
    context = mp_create_context(dimz, dimf);
    /* allocate memory for the records */
    records = malloc(n * (sizeof(int) + dimz * sizeof(fpz_t) 
                        + dimf * sizeof(fpf_t)));
#endif
    
    for(i = 0; i < n; i++)
    {
        record_base = records + (i * context->sizeof_record);
#ifdef WITH_FP
        /* set the starting order */
        *(int*)record_base = -1;
        /* set an invalid Hilbert index */
        /* *((fpz_t*)(record_base + sizeof(int))) = -1; */
#else
#endif
#ifdef WITH_PRINT
        printf("(");
#endif
        for(j = 0; j < dimz; j++)
        {
            fscanf(f, "%d", &coordz);
#ifdef WITH_FP
            ((fpz_t*)(record_base + context->offset_coordsz))[j]
                = coordz;
#else
            mpz_init(((mpz_t*)(record_base + context->offset_coordsz))[j]);
            mpz_set_ui(
                ((mpz_t*)(record_base + context->offset_coordsz))[j],
                coordz
            );
#endif
#ifdef WITH_PRINT
            if(j)
                printf(", %d", coordz);
            else
                printf("%d", coordz);
#endif
        }
        
        for(j = 0; j < dimf; j++)
        {
            fscanf(f, "%f", &coordf);
#ifdef WITH_FP
            ((fpf_t*)(record_base + context->offset_coordsf))[j] 
                = coordf;
#else
            mpf_init(((mpf_t*)(record_base + context->offset_coordsf))[j]);
            mpf_set_d(
                ((mpf_t*)(record_base + context->offset_coordsf))[j],
                coordz
            );
#endif
#ifdef WITH_PRINT
            if((dimz + j))
                printf(", %f", coordf);
            else
                printf("%f", coordf);
#endif
        }
#ifdef WITH_PRINT
        printf(")\n");
#endif
    }
#ifdef WITH_FP
    printf("Start sorting...\n");
    gettimeofday(&start_time, NULL);
    context->find_order = fp_find_order_constant;
    fp_im_sort(context, records, n, &results);
    gettimeofday(&end_time, NULL);
    runtime = ((double)end_time.tv_sec 
                    + ((double)end_time.tv_usec / 1000000.0)) 
                - ((double)start_time.tv_sec 
                    + ((double)start_time.tv_usec / 1000000.0));
    printf("Runtime: %f seconds\n", runtime);
    if(records == results)
        printf("Done!\n");
    printf("Max order is: %u\n", context->max_order);
    printf("Calls to index mapping: %llu\n", context->env->calls);
   
    /* write output file header */
    if(of)
    {
        fprintf(of, "%d\n", dimz);
        for(i = 0; i < dimz; i++)
            fprintf(of, "%d ", cards[i]);
        fprintf(of, "\n");
        fprintf(of, "%d\n", dimf);
        fprintf(of, "%d\n", n);
    }
    for(i = 0; i < n; i++)
    {
        fpz_t index;
        
        record_base = records + (i * context->sizeof_record);
#ifdef WITH_PRINT    
        printf("(");
#endif
        for(j = 0; j < dimz; j++)
        {
            coordz = ((fpz_t*)(record_base + context->offset_coordsz))[j];
#ifdef WITH_PRINT    
            if(j)
                printf(", %d", coordz);
            else
                printf("%d", coordz);
#endif

            if(of)
                fprintf(of, "%d ", coordz);
        }
        
        for(j = 0; j < dimf; j++)
        {
            coordf = ((fpf_t*)(record_base + context->offset_coordsf))[j];
#ifdef WITH_PRINT    
            if((dimz + j))
                printf(", %f", coordf);
            else
                printf("%f", coordf);
#endif
            
            if(of)
                fprintf(of, "%f ", coordf);
        }
        fpm_c2i(context->env, context->max_order, 
                (fpz_t*)(record_base + context->offset_coordsz),  
                (fpf_t*)(record_base + context->offset_coordsf), &index);
#ifdef WITH_PRINT    
        printf("): %llu\n", index);
#endif
        if(last_index > index)
        {
            printf("Failure!\n");
        }
        last_index = index;
        
        if(of)
            fprintf(of, "\n");
    }
#endif

#ifdef WITH_FP
    fp_destroy_context(context);
#else
    for(i = 0; i < n; i++)
    {
        record_base = records + (i * context->sizeof_record);
        for(j = 0; j < dimz; j++)
            mpz_clear(((mpz_t*)(records + context->offset_coordsz))[j]);
        for(j = 0; j < dimf; j++)
            mpf_clear(((mpf_t*)(records + context->offset_coordsf))[j]);
    }
    mp_destroy_context(context);
    mpz_clear(index);
#endif
    if(records)
        free(records);
    fclose(f);
    return 0;
} 
