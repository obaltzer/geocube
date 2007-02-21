#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <getopt.h>
#include <string.h>

#include "hilbert.h"
#include "sort.h"
#include "tree.h"
#include "util.h"

#define TRUE 1
#define FALSE 0

typedef struct fp_context fp_context_t;

struct tree_test_config_s
{
    char datfile[256];
    int n_queries;
    int n_updates;
    char** updates;
};
typedef struct tree_test_config_s tree_test_config_t;

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

void usage()
{
    printf("tree_test [-q QUERIES] DATASET.dat\n"
           "          [UPDATE1.dat] ... [UPDATEn.dat]\n");
    exit(-1);
}

int configure(tree_test_config_t* config, int argc, char** argv)
{
    int ch;

    static struct option longopts[] = {
        {"queries",  required_argument,  NULL,   'q'},
        {NULL,          0,                  NULL,   0}
    };

    config->n_queries = 0;
    
    while((ch = getopt_long(argc, argv, "q:", longopts, NULL)) != -1)
    {
        switch(ch)
        {
            case 'q':
                config->n_queries = atoi(optarg);
                break;
            default:
                usage();
                exit(-1);
        }
    }
    argc -= optind;
    argv += optind;

    if(argc < 1)
    {
        usage();
        exit(-1);
        return FALSE;
    }   
    else
    {   
        strncpy(config->datfile, argv[0], sizeof(config->datfile));
        config->datfile[sizeof(config->datfile) - 1] = '\0';
        config->n_updates = --argc;
        config->updates = ++argv;
    }
    return TRUE;
}

void report(struct fp_context* c, size_t n, void* r)
{
#if 1
    if(n == 1)
    {
        printf("."); 
        fflush(stdout);
    }
    else
    {
        printf("[%d]", n);
        fflush(stdout);
    }
#endif
/*
        print_record(stdout, c, (r + (i * c->record_size)));
        printf("\n");
*/
}

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

#define N_RECORDS(d) (d->n_records)
#define RECORDS(d) (d->records)
#define RECORD(d, i) (d->records + (i * d->context->record_size))
#define RECORD_SIZE(d) (d->context->record_size)
#define DIMZ(d) (d->meta->dimz)
#define DIMF(d) (d->meta->dimf)
#define INDEX(d, i) (*((int*)(d->records + (i * d->context->record_size))))
#define COORDZ(d, i, j) (*(fpz_t*)(d->records + (i * d->context->record_size) + d->context->coordsz_off + (sizeof(fpz_t) * j)))
#define COORDF(d, i, j) (*(fpf_t*)(d->records + (i * d->context->record_size) + d->context->coordsf_off + (sizeof(fpf_t) * j)))
#define RCOORDZ(r, c, j) (*(fpz_t*)(r + c->coordsz_off + (sizeof(fpz_t) * j)))
#define RCOORDF(r, c, j) (*(fpf_t*)(r + c->coordsf_off + (sizeof(fpf_t) * j)))
#define MIN(a, b) (( a <= b ? a : b))
#define MAX(a, b) (( a > b ? a : b))
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

void dataset_print(dataset_t* data)
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
        if(maxi > index)
            abort();
        else
            maxi = index;
    }
}

dataset_t* dataset_merge(dataset_t* orig, dataset_t* update)
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
    data->records = malloc((N_RECORDS(orig) + N_RECORDS(update)) * RECORD_SIZE(data));
    while(o < N_RECORDS(orig) && u < N_RECORDS(update))
    {
        res = orig->context->compare_func(orig->context, RECORD(orig, o), RECORD(update, u));
        if(res < 0)
        {
            putchar('.');
            memcpy(RECORD(data, i++), RECORD(orig, o++), RECORD_SIZE(data));
        }
        else if(res > 0)
        {
            putchar('I');
            memcpy(RECORD(data, i++), RECORD(update, u++), RECORD_SIZE(data));
        }
        else
        {
            putchar('U');
            memcpy(RECORD(data, i++), RECORD(orig, o++), RECORD_SIZE(data));
            u++;
        }
        fflush(stdout);
    }
    if(u < N_RECORDS(update))
    {
        printf("[I %d]", N_RECORDS(update) - u);
        memcpy(RECORD(data, i), RECORD(update, u), RECORD_SIZE(data) * (N_RECORDS(update) - u));
        i += N_RECORDS(update) - u;
    }
    if(o < N_RECORDS(orig))
    {
        printf("[. %d]", N_RECORDS(orig) - o);
        memcpy(RECORD(data, i), RECORD(orig, o), RECORD_SIZE(data) * (N_RECORDS(orig) - o));
        i += N_RECORDS(orig) - o;
    }
    fflush(stdout);
    data->n_records = i;
    return data;
}

int main(int argc, char** argv)
{
    FILE* f = NULL;
    int i;
    void* results = NULL;
    struct sort_config config;
    /* fpz_t index; */
    struct fp_context* context;
    struct fp_im_tnode* root;
    struct runtime rt;
    double rtime;
    double qtime;
    double avgqtime;
    void* qr_min;
    void* qr_max;
    tree_test_config_t tc;
    int query;
    metadata_t* meta;
    dataset_t* data;
    int c_updates = 0;

    configure(&tc, argc, argv);
    if((f = fopen(tc.datfile, "r")) == NULL)
    {
        perror("Unable to open datafile:");
        exit(3);
    } 
   
    meta = metadata_read(f);
    
    /* create a new context */
    config.verbose = 1;
    config.normalize = 0;
    config.denormalize = 0;
    config.benchmark = 1;
    config.find_order = ITERATIVE;
    config.index = KEEP;
    config.cmp = HILBERT;
    config.print = stdout;
    context = fp_create_context(&config, meta->dimz, meta->dimf, meta->start_order);
    data = dataset_read(f, meta, context);
    fclose(f);

    printf("N is: %d\n", N_RECORDS(data));
    
    rtime = 0.0;
    start_timing(&rt);
    fp_im_sort(context, RECORDS(data), N_RECORDS(data), &results);
    stop_timing(&rt);
    rtime = get_runtime(rt);
    printf("Sort time: %f\n", rtime);
    srand(0);
    while(RECORDS(data) == results)
    {
        printf("Done!\n");
        /* dataset_print(data); */
        /* build a tree based on the sorted records */
        start_timing(&rt);
        root = fp_im_build_tree(context, N_RECORDS(data), RECORDS(data), LEAF);
        stop_timing(&rt);
        rtime = get_runtime(rt);
        printf("Build Tree Time: %f\n", rtime);
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
        qr_min = malloc(context->record_size); 
        qr_max = malloc(context->record_size);
        qtime = 0.0;
        avgqtime = 0.0;
        for(query = 0; query < tc.n_queries; query++)
        {
            for(i = 0; i < context->env->dimz; i++)
            {
                fpz_t a = (fpz_t)((double)rand() / (double)RAND_MAX * (double)meta->cards[i]);
                fpz_t b = (fpz_t)((double)rand() / (double)RAND_MAX * (double)meta->cards[i]);
                while(a == b)
                    a = (fpz_t)((double)rand() / (double)RAND_MAX * (double)meta->cards[i]);
                
                RCOORDZ(qr_min, context, i) = MIN(a, b);
                RCOORDZ(qr_max, context, i) = MAX(a, b);
            }
            for(i = 0; i < context->env->dimf; i++)
            {
                fpf_t a = (fpf_t)((double)rand() / (double)RAND_MAX);
                fpf_t b = (fpf_t)((double)rand() / (double)RAND_MAX);
                while(a == b)
                    a = (fpf_t)((double)rand() / (double)RAND_MAX);
                
                RCOORDF(qr_min, context, i) = MIN(a, b);
                RCOORDF(qr_max, context, i) = MAX(a, b);
            }
            printf("Querying tree...\n");
            printf("Query region min: ");
            print_record(stdout, context, qr_min);
            printf("\nQuery region max: ");
            print_record(stdout, context, qr_max);
            printf("\n");
            start_timing(&rt);
            i = fp_im_query_tree(context, root, qr_min, qr_max, report);
            stop_timing(&rt);
            printf("Reported %d records\n", i);
            qtime = get_runtime(rt);
            printf("Query time: %f\n", qtime);
            avgqtime += qtime;
        }
        printf("Average query time: %f\n", avgqtime / (double)tc.n_queries);
        free(qr_min);
        free(qr_max);
        fp_im_destroy_tree(root);

        if(c_updates < tc.n_updates)
        {
            FILE* u;
            metadata_t* umeta;
            dataset_t* udata;
            dataset_t* new_data;

            printf("*** Perform Update %d ***\n", c_updates + 1); 
            if((u = fopen(tc.updates[c_updates], "r")) == NULL)
            {
                printf("Cannot open update file %s.\n", tc.updates[c_updates]);
                return -1;
            }
            umeta = metadata_read(u);
            if(!metadata_equal(meta, umeta))
            {
                printf("Metadata of update dataset %s does not match.\n", tc.updates[c_updates]);
                return -1;
            }
            metadata_destroy(umeta);
            udata = dataset_read(u, meta, context);
            fclose(u);
            printf("N is: %d\n", N_RECORDS(udata));
            rtime = 0.0;
            results = NULL;
            start_timing(&rt);
            fp_im_sort(context, RECORDS(udata), N_RECORDS(udata), &results);
            stop_timing(&rt);
            rtime = get_runtime(rt);
            printf("Sort time: %f\n", rtime);
            new_data = dataset_merge(data, udata);
            if(new_data == NULL)
            {
                printf("Error merging datasets.\n");
                return -1;
            }
            dataset_destroy(udata);
            dataset_destroy(data);
            data = new_data;
            results = RECORDS(data);
            c_updates++;
        }
        else
            results = NULL;
    }
    dataset_destroy(data);
    fp_destroy_context(context);
    metadata_destroy(meta);
    return 0;
} 
