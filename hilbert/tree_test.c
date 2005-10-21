#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

#include "hilbert.h"
#include "sort.h"
#include "tree.h"

void report(struct fp_context* c, size_t n, void* r)
{
    size_t i = 0;
    for(; i < n; i++)
    {
        print_record(c, (r + (i * c->sizeof_record)), 0);
        printf("\n");
    }
}
    
int main(int argc, char** argv)
{
    FILE* f = NULL;
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
#ifdef WITH_FP
    void* records = NULL;
    /* fpz_t index; */
    struct fp_context* context;
    struct fp_im_tnode* root;
    void* qr_min;
    void* qr_max;
#else
    void* records = NULL;
    mpz_t index;
    struct mp_context* context;
#endif
    
    if(argc < 2)
    {
        printf("Usage: sort_test datafile\n");
        exit(2);
    }
    
    if((f = fopen(argv[1], "r")) == NULL)
    {
        perror("Unable to open datafile:");
        exit(3);
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
    records = malloc(n * context->sizeof_record);
#endif

    for(i = 0; i < n; i++)
    {
        record_base = records + (i * context->sizeof_record);
#ifdef WITH_FP
        /* set the records order to be uninitialized */
        *(int*)record_base = -1;
#else
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
        }
    }
#ifdef WITH_FP
    printf("N is: %d\n", n);
    fp_im_sort(context, records, n, &results);
    if(records == results)
    {
        fpz_t maxi = 0;
        printf("Done!\n");
#if 0
        /* print out sorted records */
        for(i = 0; i < n; i++)
        {
            fpz_t index;
            
            record_base = records + (i * context->sizeof_record);
            printf("(");
            for(j = 0; j < dimz; j++)
            {
                coordz = ((fpz_t*)(record_base + context->offset_coordsz))[j];
                if(j)
                    printf(", %d", coordz);
                else
                    printf("%d", coordz);
            }
            
            for(j = 0; j < dimf; j++)
            {
                coordf = ((fpf_t*)(record_base + context->offset_coordsf))[j];
                if((dimz + j))
                    printf(", %f", coordf);
                else
                    printf("%f", coordf);
            }
            fpm_c2i(context->env, context->max_order, 
                    (fpz_t*)(record_base + context->offset_coordsz),  
                    (fpf_t*)(record_base + context->offset_coordsf), &index);
            /* check if the order of indices is correct */
            if(maxi > index)
                abort();
            else
                maxi = index;
            /* printf("): %llu\n", index); */
        }
#endif

        /* build a tree based on the sorted records */
        root = fp_im_build_tree(context, n, records, LEAF);
        printf("Root bbox: ");
        fp_print_bbox(context, root);
        printf("\n");
        printf("Maximum order: %d\n", context->max_order);
        printf("Number of children: %d\n", root->n_children);
        printf("Number of leaves: %d\n", root->n_leaves);
        printf("Number of levels: %d\n", root->level);
        printf("Calls to index mapping: %llu\n", context->env->calls);
        printf("Calls to build_tree: %llu\n", context->build_tree_calls);
        printf("Number of tree nodes: %d\n", context->n_tree_nodes);
        qr_min = malloc(context->sizeof_record); 
        qr_max = malloc(context->sizeof_record);
        srand(time(NULL));
        for(i = 0; i < context->env->dimz; i++)
        {
#if 1
            fpz_t a = (fpz_t)((double)rand() / (double)RAND_MAX * (double)cards[i]);
            fpz_t b = (fpz_t)((double)rand() / (double)RAND_MAX * (double)cards[i]);
            while(a == b)
                a = (fpz_t)((double)rand() / (double)RAND_MAX * (double)cards[i]);
#endif
#if 0
            fpz_t a = 0;
            fpz_t b = cards[i];
#endif
            
            if(a < b)
            {
                ((fpz_t*)(qr_min + context->offset_coordsz))[i] = a;
                ((fpz_t*)(qr_max + context->offset_coordsz))[i] = b;
            }
            else
            {
                ((fpz_t*)(qr_min + context->offset_coordsz))[i] = b;
                ((fpz_t*)(qr_max + context->offset_coordsz))[i] = a;
            }
        }
        for(i = 0; i < context->env->dimf; i++)
        {
#if 1
            fpf_t a = (fpf_t)((double)rand() / (double)RAND_MAX);
            fpf_t b = (fpf_t)((double)rand() / (double)RAND_MAX);
            while(a == b)
                a = (fpf_t)((double)rand() / (double)RAND_MAX);
#endif
#if 0
            fpf_t a = 0.0;
            fpf_t b = 1.0;
#endif
            
            if(a < b)
            {
                ((fpf_t*)(qr_min + context->offset_coordsf))[i] = a;
                ((fpf_t*)(qr_max + context->offset_coordsf))[i] = b;
            }
            else
            {
                ((fpf_t*)(qr_min + context->offset_coordsf))[i] = b;
                ((fpf_t*)(qr_max + context->offset_coordsf))[i] = a;
            }
        }
        printf("Querying tree...\n");
        printf("Query region min: ");
        print_record(context, qr_min, 0);
        printf("\nQuery region max: ");
        print_record(context, qr_max, 0);
        printf("\n");
        i = fp_im_query_tree(context, root, qr_min, qr_max, report);
        printf("Reported %d records\n", i);
        fp_im_destroy_tree(root);
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
