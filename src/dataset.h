#ifndef __DATASET_H
#define __DATASET_H

#include "sort.h"

#define TRUE 1
#define FALSE 0

#define N_RECORDS(d) (d->n_records)
#define RECORDS(d) (d->records)
#define RECORD(d, i) (d->records + (i * d->context->record_size))
#define RECORD_S(d) (d->context->record_size)
#define DIMZ(d) (d->meta->dimz)
#define DIMF(d) (d->meta->dimf)
#define INDEX(d, i) (*((int*)(d->records + (i * d->context->record_size))))
#define COORDZ(d, i, j) (*(fpz_t*)(d->records + (i * d->context->record_size) + d->context->coordsz_off + (sizeof(fpz_t) * j)))
#define COORDF(d, i, j) (*(fpf_t*)(d->records + (i * d->context->record_size) + d->context->coordsf_off + (sizeof(fpf_t) * j)))
#define RCOORDZ(r, c, j) (*(fpz_t*)(r + c->coordsz_off + (sizeof(fpz_t) * j)))
#define RCOORDF(r, c, j) (*(fpf_t*)(r + c->coordsf_off + (sizeof(fpf_t) * j)))
#define MIN(a, b) (( a <= b ? a : b))
#define MAX(a, b) (( a > b ? a : b))

typedef struct fp_context fp_context_t;

struct metadata_s
{
    size_t dimz;
    size_t dimf;
    size_t max_card;
    size_t* cards;
    int start_order;
};
typedef struct metadata_s metadata_t;

struct dataset_s
{
    metadata_t* meta;
    fp_context_t* context;
    int n_records;
    void* records;
};
typedef struct dataset_s dataset_t;

metadata_t* metadata_read(FILE* f);
int metadata_equal(metadata_t* m1, metadata_t* m2);
void metadata_destroy(metadata_t* meta);
dataset_t* dataset_read(FILE* f, metadata_t* meta, fp_context_t* context);
void dataset_destroy(dataset_t* data);
void dataset_print(dataset_t* data, int check);
void dataset_verify(dataset_t* data);
dataset_t* dataset_merge(dataset_t* orig, dataset_t* update, int silent);
dataset_t* dataset_copy(dataset_t* orig);
dataset_t* dataset_create(metadata_t* meta, fp_context_t* context, size_t n);

#endif
