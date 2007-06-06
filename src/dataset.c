#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "hilbert.h"
#include "dataset.h"

metadata_t* metadata_read(FILE* f)
{
    metadata_t* meta;
    size_t max_card;
    int i;
    
    meta = malloc(sizeof(metadata_t));
    if(meta == NULL)
    {
        printf("Cannot allocate memory for meta data.\n");
        return NULL;
    }

    meta->max_card = 0;
    /* read number of discrete space dimensions */
    fscanf(f, "%d", &meta->dimz);
    /* read the cardinalities of the dimensions and determine 
     * the largest */
    meta->cards = (size_t*)malloc(sizeof(size_t) * meta->dimz);
    for(i = 0; i < meta->dimz; i++)
    {
        fscanf(f, "%d", &meta->cards[i]);
        if(meta->cards[i] > meta->max_card)
            meta->max_card = meta->cards[i];
    }
    /* read number of continuous space dimensions */
    fscanf(f, "%d\n", &meta->dimf);
    meta->start_order = 1;
    max_card = meta->max_card;
    while(max_card >> 1)
    {
        max_card >>= 1;
        meta->start_order++;
    }
    return meta;
}

int metadata_equal(metadata_t* m1, metadata_t* m2)
{
    if(m1->dimz == m2->dimz && m1->dimf == m2->dimf && m1->max_card == m2->max_card)
    {
        int i;
        for(i = 0; i < m1->dimz; i++)
        {
            if(m1->cards[i] != m2->cards[i])
                return FALSE;
        }
        return TRUE;
    }
    return FALSE;
}

void metadata_destroy(metadata_t* meta)
{
    free(meta->cards);
    free(meta);
}

dataset_t* dataset_create(metadata_t* meta, fp_context_t* context, size_t n)
{
    dataset_t* data;
    
    data = malloc(sizeof(dataset_t));
    if(data == NULL)
    {
        printf("Cannot allocate memory for dataset structure.\n");
        return NULL;
    }
    data->meta = meta;
    data->context = context;
    data->n_records = n;
    data->records = malloc(data->n_records * context->record_size);
    assert(data->records != NULL); 
    return data;
}
    

dataset_t* dataset_read(FILE* f, metadata_t* meta, fp_context_t* context)
{
    dataset_t* data;
    int i;
    int j;

    data = malloc(sizeof(dataset_t));
    if(data == NULL)
    {
        printf("Cannot allocate memory for dataset structure.\n");
        return NULL;
    }
    data->meta = meta;
    data->context = context;
    /* number of records */
    fscanf(f, "%d\n", &data->n_records);
    
    /* allocate memory for the records */
    data->records = malloc(data->n_records * context->record_size);
    printf("Allocate Memory: %d = %dM\n", 
        data->n_records * context->record_size,
        data->n_records * context->record_size / 1024 / 1024);
    assert(data->records != NULL);

    for(i = 0; i < N_RECORDS(data); i++)
    {
        /* set the records order to be uninitialized */
        INDEX(data, i) = -1;
        for(j = 0; j < DIMZ(data); j++)
        {
            fscanf(f, "%Ld", &COORDZ(data, i, j));
        }
        
        for(j = 0; j < DIMF(data); j++)
        {
            fscanf(f, "%lf", &COORDF(data, i, j));
        }
    }
    
    return data;
}

void dataset_destroy(dataset_t* data)
{
    free(data->records);
    free(data);
}

void dataset_print(dataset_t* data, int check)
{
    int i;
    fpz_t maxi = 0;
    int j;

    /* print out sorted records */
    for(i = 0; i < data->n_records; i++)
    {
        fpz_t index;
        
        printf("(");
        for(j = 0; j < DIMZ(data); j++)
        {
            if(j)
                printf(", %Ld", COORDZ(data, i, j));
            else
                printf("%Ld", COORDZ(data, i, j));
        }
        
        for(j = 0; j < DIMF(data); j++)
        {
            if((DIMZ(data) + j))
                printf(", %lf", COORDF(data, i, j));
            else
                printf("%lf", COORDF(data, i, j));
        }
        fpm_c2i(data->context->env, data->context->max_order, 
                &COORDZ(data, i, 0),  
                &COORDF(data, i, 0), &index);
         printf("): %llu\n", index);
         fflush(stdout);
        /* check if the order of indices is correct */
        if(check)
        {
            if(maxi > index)
                abort();
            else
                maxi = index;
        }
    }
}

dataset_t* dataset_merge(dataset_t* orig, dataset_t* update, int silent)
{
    dataset_t* data;
    size_t o = 0;
    size_t u = 0;
    size_t i = 0;
    int res;

    data = malloc(sizeof(dataset_t));
    if(data == NULL)
    {
        printf("Cannot allocate memory for updated datset.\n");
        return NULL;
    }
    data->meta = orig->meta;
    data->context = orig->context;
    data->records = malloc((N_RECORDS(orig) + N_RECORDS(update)) * RECORD_S(data));
    while(o < N_RECORDS(orig) && u < N_RECORDS(update))
    {
        res = orig->context->compare_func(orig->context, RECORD(orig, o), RECORD(update, u));
        if(res < 0)
        {
            if(!silent) putchar('.');
            memcpy(RECORD(data, i++), RECORD(orig, o++), RECORD_S(data));
        }
        else if(res > 0)
        {
            if(!silent) putchar('I');
            memcpy(RECORD(data, i++), RECORD(update, u++), RECORD_S(data));
        }
        else
        {
            if(!silent) putchar('U');
            memcpy(RECORD(data, i++), RECORD(orig, o++), RECORD_S(data));
            u++;
        }
        if(!silent) fflush(stdout);
    }
    if(u < N_RECORDS(update))
    {
        if(!silent) printf("[I %d]", N_RECORDS(update) - u);
        memcpy(RECORD(data, i), RECORD(update, u), RECORD_S(data) * (N_RECORDS(update) - u));
        i += N_RECORDS(update) - u;
    }
    if(o < N_RECORDS(orig))
    {
        if(!silent) printf("[. %d]", N_RECORDS(orig) - o);
        memcpy(RECORD(data, i), RECORD(orig, o), RECORD_S(data) * (N_RECORDS(orig) - o));
        i += N_RECORDS(orig) - o;
    }
    if(!silent) fflush(stdout);
    data->n_records = i;
    return data;
}

dataset_t* dataset_copy(dataset_t* orig)
{
    dataset_t* copy;

    copy = malloc(sizeof(dataset_t));
    if(copy == NULL)
    {
        printf("Cannot allocate memory for dataset structure.\n");
        return NULL;
    }
    copy->context = orig->context;
    copy->meta = orig->meta;
    copy->n_records = orig->n_records;
    copy->records = malloc(N_RECORDS(orig) * RECORD_S(orig));
    if(copy->records == NULL)
    {
        printf("Cannot allocate memory for records.\n");
        free(copy);
        return NULL;
    }
    memcpy(copy->records, orig->records, N_RECORDS(orig) * RECORD_S(orig));
    return copy;
}

void dataset_verify(dataset_t* data)
{
    int i;
    fpz_t maxi = 0;

    /* print out sorted records */
    for(i = 0; i < data->n_records; i++)
    {
        fpz_t index;
        
        fpm_c2i(data->context->env, data->context->max_order, 
                &COORDZ(data, i, 0),  
                &COORDF(data, i, 0), &index);
        /* check if the order of indices is correct */
        assert(index >= maxi);
    }
}

