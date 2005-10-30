#include "util.h"
#include "sort.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>

struct metadata
{
    size_t dimz;
    size_t dimf;
    size_t max_card;
    size_t* cards;
    size_t n;
};
    
struct sort_config
{
    int verbose;
    enum { CONSTANT, ITERATIVE } find_order;
    enum { KEEP, FORGET } index;
    FILE* infile;
    FILE* outfile;
    FILE* print;
};

void configure(int argc, char** argv, struct sort_config* config)
{
    /*
    o -- output file
    v -- enable verbose output
    f -- find order algorithm (constant | iterative)
    c -- how to deal with Hilbert codes (keep | forget)
    i -- input file name or "-" for stdin
    */
    char* options = "o:vf:c:i:";
    static struct option long_options[] = {
        {"output-file", 1, 0, 'o'},
        {"verbose",     0, 0, 'v'},
        {"find-order",  1, 0, 'f'},
        {"code",        1, 0, 'c'},
        {"input-file",  1, 0, 'i'},
        {0,             0, 0, 0}
    };
    int option_index = 0;
    char c;
    while((c = getopt_long(argc, argv, options, 
                           long_options, &option_index)) != -1)
    {
        switch(c)
        {
            case 'i':
                if(!strcmp(optarg, "-"))
                    config->infile = stdin;
                else if((config->infile = fopen(optarg, "r")) == NULL)
                {
                    perror("Unable to open input file: ");
                    exit(3);
                }
                break;
                
            case 'o':
                if(!strcmp(optarg, "-"))
                {
                    config->outfile = stdout;
                    config->print = stderr;
                }
                else if((config->outfile = fopen(optarg, "w")) == NULL)
                {
                    perror("Unable to open output file: ");
                    exit(4);
                }
                break;
            
            case 'v':
                config->verbose = 1;
                break;
            
            case 'f':
                if(!strcmp(optarg, "constant"))
                    config->find_order = CONSTANT;
                else if(!strcmp(optarg, "iterative"))
                    config->find_order = ITERATIVE;
                else
                {
                    fprintf(stderr, "Unknown find-order algorithm: %s\n",
                            optarg);
                    exit(100);
                }
                break;

            case 'c':
                if(!strcmp(optarg, "keep"))
                    config->index = KEEP;
                else if(!strcmp(optarg, "forget"))
                    config->index = FORGET;
                else
                {
                    fprintf(stderr, "Unknown code handling method: %s\n",
                            optarg);
                    exit(101);
                }
                break;
            
            default:
                fprintf(stderr, "Unkown command line option: %s\n",
                        argv[optind]);
                exit(5);
        }
    }
}

void close_files(struct sort_config* config)
{
    if(config->infile != NULL)
        fclose(config->infile);
    if(config->outfile != NULL)
        fclose(config->outfile);
    if(config->print != NULL)
        fclose(config->print);
}

struct metadata* get_metadata(struct sort_config* config)
{
    int i = 0;
    struct metadata* meta = 
        (struct metadata*)malloc(sizeof(struct metadata));
        
    /* read the number of discrete space dimensions */
    fscanf(config->infile, "%d", &meta->dimz);
    /* read the cardinalities of the dimensions and determine 
     * the largest */
    meta->cards = (size_t*)malloc(sizeof(size_t) * meta->dimz);
    meta->max_card = 0;
    for(; i < meta->dimz; i++)
    {
        fscanf(config->infile, "%d", &meta->cards[i]);
        if(meta->cards[i] > meta->max_card)
            meta->max_card = meta->cards[i];
    }
    /* read the number of continuous space dimensions */
    fscanf(config->infile, "%d", &meta->dimf);
    /* number of records */
    fscanf(config->infile, "%d", &meta->n);
    return meta;
}

void release_metadata(struct metadata* meta)
{
    if(meta != NULL)
    {
        if(meta->cards != NULL)
            free(meta->cards);
        free(meta);
    }
}

int compute_start_order(struct sort_config* config, 
                        struct metadata* meta)
{
    int order = 1;
    while(meta->max_card >> order)
        order++;
    return order;
}

int main(int argc, char** args)
{
    struct metadata* meta;
    struct sort_config config
        = { 0, ITERATIVE, KEEP, stdin, NULL, stdout };
    int start_order;
    void* records;
    void* result;
    void* record_base;
    struct runtime loading_time;
    struct runtime sorting_time;
#ifdef WITH_FP
    struct fp_context* context;
#elif WITH_MP
#endif
    size_t i;
    size_t j;
    
    configure(argc, args, &config);
    
    /* get the metadata of the dataset */
    meta = get_metadata(&config);
    
    /* compute the start order for computing the hilbert codes */
    start_order = compute_start_order(&config, meta);
    fprintf(config.print, "Using start order: %d\n", start_order);

#ifdef WITH_FP
    /* create a new computation context */
    context = fp_create_context(meta->dimz, meta->dimf, start_order);
#elif WITH_MP
#endif
    /* allocate memory for the records */
    records = malloc(meta->n * context->record_size);
    fprintf(config.print, "Reading records...\n");
    start_timing(&loading_time);
    for(i = 0; i < meta->n; i++)
    {
        int coordz;
        float coordf;
        
        /* compute address base for the current record */
        record_base = records + (i * context->record_size);
        /* set initial order of the record to undefined */
        *(int*)(record_base + context->order_off) = -1;
        
        /* read discrete space dimensions */
        for(j = 0; j < meta->dimz; j++)
        {
            fscanf(config.infile, "%d", &coordz);
#ifdef WITH_FP
            ((fpz_t*)(record_base + context->coordsz_off))[j] = coordz;
#elif WITH_MP
#endif
        }

        /* read continuous space dimensions */
        for(j = 0; j < meta->dimf; j++)
        {
            fscanf(config.infile, "%f", &coordf);
#ifdef WITH_FP
            ((fpf_t*)(record_base + context->coordsf_off))[j] = coordf;
#elif
#endif
        }
        
        if(config.verbose)
        {
            print_record(config.print, context, record_base);
            printf("\n");
        }
    }
    stop_timing(&loading_time);
    fprintf(config.print, "Loading took %f seconds\n",
            (float)get_runtime(loading_time));

    fprintf(config.print, "Start sorting...\n");
    
    if(config.find_order == ITERATIVE)
        context->find_order = fp_find_order_iterative;
    else if(config.find_order == CONSTANT)
        context->find_order = fp_find_order_constant;
    
    start_timing(&sorting_time);
    fp_im_sort(context, records, meta->n, &result);
    assert(records == result);
    stop_timing(&sorting_time);
    fprintf(config.print, "Sorting runtime was %f seconds\n",
            (float)get_runtime(sorting_time));

    /* verify the results */
    {
#ifdef WITH_FP
        fpz_t index;
        fpz_t prev_index = 0;
#elif WITH_MP
#endif
        for(i = 0; i < meta->n; i++)
        {
            record_base = records + (i * context->record_size);
#ifdef WITH_FP
            fpm_c2i(context->env, context->max_order, 
                    (fpz_t*)(record_base + context->coordsz_off),  
                    (fpf_t*)(record_base + context->coordsf_off), &index);
            if(prev_index > index)
            {
                printf("Failure: %llu %llu\n", index, prev_index);
                abort();
            }
            prev_index = index;
#elif WITH_MP
#endif
        }
    }
    
    /* print out the sorted dataset */
    if(config.outfile)
    {
        int coordz;
        float coordf;
        
        fprintf(config.outfile, "%d\n", meta->dimz);
        for(i = 0; i < meta->dimz; i++)
            fprintf(config.outfile, "%d ", meta->cards[i]);

        fprintf(config.outfile, "\n");
        fprintf(config.outfile, "%d\n", meta->dimf);
        fprintf(config.outfile, "%d\n", meta->n);
        for(i = 0; i < meta->n; i++)
        {
            record_base = records + (i * context->record_size);
            for(j = 0; j < meta->dimz; j++)
            {
#ifdef WITH_FP
                coordz = ((fpz_t*)(record_base + context->coordsz_off))[j];
                fprintf(config.outfile, "%d ", coordz);
#elif WITH_MP
#endif
            }
            for(j = 0; j < meta->dimf; j++)
            {
#ifdef WITH_FP
                coordf = ((fpf_t*)(record_base + context->coordsf_off))[j];
                fprintf(config.outfile, "%f ", coordf);
#elif WITH_MP
#endif
            }
            fprintf(config.outfile, "\n");
        }
    }
    
#ifdef WITH_FP
    fp_destroy_context(context);
#elif WITH_MP
#endif
    
    if(records)
        free(records);

    release_metadata(meta);
    close_files(&config);
    return 0;
}
