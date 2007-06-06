#include "util.h"
#include "sort.h"
#include "norm.h"

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
    
void configure(int argc, char** argv, struct sort_config* config)
{
    /*
    o -- output file
    v -- enable verbose output
    f -- find order algorithm (constant | iterative)
    s -- Hilbert codes with record (yes/1 | no/0)
    i -- input file name or "-" for stdin
    */
    char* options = "o:vnf:s:i:dbc:";
    static struct option long_options[] = {
        {"output-file", 1, 0, 'o'},
        {"verbose",     0, 0, 'v'},
        {"normalize",   0, 0, 'n'},
        {"find-order",  1, 0, 'f'},
        {"store-code",  1, 0, 's'},
        {"input-file",  1, 0, 'i'},
        {"denormalize", 0, 0, 'd'},
        {"benchmark",   0, 0, 'b'},
        {"compare-function", 1, 0, 'c'},
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
                
            case 'b':
                config->benchmark = 1;
                if(config->outfile == NULL)
                    config->outfile = stdout;
                config->print = stderr;
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
                config->filename = optarg;
                break;
            
            case 'v':
                config->verbose = 1;
                break;
            
            case 'f':
                if(!strcmp(optarg, "compare"))
                    config->find_order = COMPARE;
                else if(!strcmp(optarg, "constant"))
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
                if(!strcmp(optarg, "hilbert"))
                    config->cmp = HILBERT;
                else if(!strcmp(optarg, "xyz"))
                    config->cmp = XYZ;
                else
                {
                    fprintf(stderr, "Unknown comparison algorithm: %s\n",
                            optarg);
                    exit(100);
                }
                break;

            case 's':
                if(!strcmp(optarg, "yes") || !strcmp(optarg, "1"))
                    config->index = KEEP;
                else if(!strcmp(optarg, "no") || !strcmp(optarg, "0"))
                    config->index = FORGET;
                else
                {
                    fprintf(stderr, "Unknown code handling method: %s\n",
                            optarg);
                    exit(101);
                }
                break;
            
            case 'n':
                config->normalize = 1;
                break; 
            
            case 'd':
                config->denormalize = 1;
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

void read_records(struct sort_config* config, 
                  struct metadata* meta,
                  struct fp_context* context,
                  void* records)
{
    void* record_base;
    struct runtime loading_time;
    size_t i;
    size_t j;

    fprintf(config->print, "Reading records...\n");
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
            fscanf(config->infile, "%d", &coordz);
#ifdef WITH_FP
            ((fpz_t*)(record_base + context->coordsz_off))[j] = coordz;
#elif WITH_MP
#endif
        }

        /* read continuous space dimensions */
        for(j = 0; j < meta->dimf; j++)
        {
            fscanf(config->infile, "%f", &coordf);
#ifdef WITH_FP
            ((fpf_t*)(record_base + context->coordsf_off))[j] = coordf;
#elif
#endif
        }
        
        if(config->verbose)
        {
            print_record(config->print, context, record_base);
            fprintf(config->print, "\n");
        }
    }
    stop_timing(&loading_time);
    fprintf(config->print, "Loading took %f seconds\n",
            (float)get_runtime(loading_time));
}
    
int main(int argc, char** args)
{
    struct metadata* meta;
    struct sort_config config
        = { 0, 0, 0, 0, ITERATIVE, KEEP, HILBERT, stdin, NULL, stdout, NULL };
    int start_order;
    void* records;
    void* result = NULL;
    void* record_base;
    struct runtime sorting_time;
    struct runtime normalization;
#ifdef WITH_FP
    struct fp_context* context;
    struct fp_norm_context* norm;
#elif WITH_MP
#endif
    size_t i;
    size_t j;
    
    configure(argc, args, &config);
    
    /* get the metadata of the dataset */
    meta = get_metadata(&config);
    
    /* compute the start order for computing the hilbert codes */
    start_order = compute_start_order(&config, meta);
    
    if(config.normalize)
    {
        void* orig_records;
        norm = fp_create_norm_context(&config, meta->dimz, meta->dimf, 
                                      start_order);
        /* allocate memory for original records and read them */
        orig_records = malloc(meta->n * norm->mc->record_size);
        read_records(&config, meta, norm->mc, orig_records);
        
        /* the working context is the one for the normalized records */
        context = norm->zc;
        records = malloc(meta->n * context->record_size);

#ifdef WITH_FP
        start_timing(&normalization);
        fp_normalize(norm, orig_records, meta->n, records);
        stop_timing(&normalization);
#elif WITH_MP
#endif
        free(orig_records);
    }
    else
    {
        /* create a new computation context */
        context = fp_create_context(&config, meta->dimz, meta->dimf, 
                                    start_order);
        /* allocate memory for the records */
        records = malloc(meta->n * context->record_size);
        if(!records)
        {
            fprintf(config.print, "Not enough memory!\n");
            fclose(config.outfile);
            if(config.filename)
            {
                remove(config.filename);
            }
            exit(123);
        }
        read_records(&config, meta, context, records);
    }
    
    /* set the input cardinalities */
    for(i = 0; i < meta->dimz; i++)
        context->cards[i] = meta->cards[i];
    
    fprintf(config.print, "Using start order: %d\n", context->start_order);
    /*************************************************/
     
    fprintf(config.print, "Start sorting...\n");
    
    start_timing(&sorting_time);
    fp_im_sort(context, records, meta->n, &result);
    stop_timing(&sorting_time);
    assert(records == result);
    fprintf(config.print, "Sorting runtime was %f seconds\n",
            (float)get_runtime(sorting_time));
    fprintf(config.print, "Maximum order is: %d\n", context->max_order);
    fprintf(config.print, "Calls to index mapping: %llu\n",
            context->env->calls);
    fprintf(config.print, "Calls to compare function: %llu\n",
            context->compare_calls);
    if(config.benchmark == 1)
    {
        if(config.normalize)
            fprintf(config.outfile, "%f\t%f\t%llu\n",
                (float)get_runtime(sorting_time), 
                (float)get_runtime(normalization),
                context->env->calls);
        else
            fprintf(config.outfile, "%f\t%f\t%llu\n",
                (float)get_runtime(sorting_time), 0.0f,
                context->env->calls);
    }
        
    /* verify the results */
#if 0
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
                fprintf(config.print, "Failure: %llu %llu\n", 
                        index, prev_index);
                abort();
            }
            prev_index = index;
#elif WITH_MP
#endif
        }
    }
#endif

    if(config.normalize && config.denormalize)
    {
        /* reverse map records */
        void* denorm_records = malloc(meta->n * norm->mc->record_size);
#ifdef WITH_FP
        fp_denormalize(norm, records, meta->n, denorm_records);
#elif WITH_MP
#endif
        context = norm->mc;
        free(records);
        records = denorm_records;
    }
    
    /* print out the sorted dataset */
    if(config.outfile && !config.benchmark)
    {
        int coordz;
        float coordf;
        
        fprintf(config.outfile, "%d\n", context->env->dimz);
        for(i = 0; i < context->env->dimz; i++)
            fprintf(config.outfile, "%d ", context->cards[i]);

        fprintf(config.outfile, "\n");
        fprintf(config.outfile, "%d\n", context->env->dimf);
        fprintf(config.outfile, "%d\n", meta->n);
        for(i = 0; i < meta->n; i++)
        {
            record_base = records + (i * context->record_size);
            for(j = 0; j < context->env->dimz; j++)
            {
#ifdef WITH_FP
                coordz = ((fpz_t*)(record_base + context->coordsz_off))[j];
                fprintf(config.outfile, "%d ", coordz);
#elif WITH_MP
#endif
            }
            for(j = 0; j < context->env->dimf; j++)
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
    
    if(config.normalize)
    {
#ifdef WITH_FP
        fp_destroy_norm_context(norm);
#elif WITH_MP
#endif
    }
    else
    {
#ifdef WITH_FP
        fp_destroy_context(context);
#elif WITH_MP
#endif
    }
    
    if(records)
        free(records);

    release_metadata(meta);
    close_files(&config);
    return 0;
}
