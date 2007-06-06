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
#include "dataset.h"
#include "norm.h"

struct tree_test_config_s
{
    char datfile[256];
    char logfile[256];
    int n_queries;
    int n_sorts;
    int n_tree_builds;
    int n_merges;
    int n_updates;
    int forget;
    int constant;
    int xyz;
    int normalize;
    char** updates;
};
typedef struct tree_test_config_s tree_test_config_t;

int reported = 0;

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
        {"queries",     required_argument,  NULL,   'q'},
        {"sorts",       required_argument,  NULL,   's'},
        {"tree-builds", required_argument,  NULL,   't'},
        {"merges",      required_argument,  NULL,   'm'},
        {"logfile",     required_argument,  NULL,   'l'},
        {"forget",      no_argument,        NULL,   'f'},
        {"constant",    no_argument,        NULL,   'c'},
        {"xyz",         no_argument,        NULL,   'x'},
        {"normalize",   no_argument,        NULL,   'n'},
        {NULL,          0,                  NULL,   0}
    };

    config->n_queries = 0;
    config->n_sorts = 1;
    config->n_tree_builds = 1;
    config->n_merges = 1;
    strcpy(config->logfile, "timelog.csv");
    config->forget = FALSE;
    config->constant = FALSE;
    config->xyz = FALSE;
    config->normalize = FALSE;
    while((ch = getopt_long(argc, argv, "q:s:t:m:l:fcxn", longopts, NULL)) != -1)
    {
        switch(ch)
        {
            case 'q':
                config->n_queries = atoi(optarg);
                break;
            case 't':
                config->n_tree_builds = atoi(optarg);
                break;
            case 's':
                config->n_sorts = atoi(optarg);
                break;
            case 'm':
                config->n_merges = atoi(optarg);
                break;
            case 'l':
                strncpy(config->logfile, optarg, sizeof(config->logfile));
                break;
            case 'f':
                config->forget = TRUE;
                break;
            case 'c':
                config->constant = TRUE;
                break;
            case 'x':
                config->xyz = TRUE;
                break;
            case 'n':
                config->normalize = TRUE;
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
void report_quiet(struct fp_context* c, size_t n, void* r)
{
    reported++;
}

int main(int argc, char** argv)
{
    FILE* f = NULL;
    int i;
    void* results = NULL;
    struct sort_config config;
    /* fpz_t index; */
    struct fp_context* context;
    struct fp_norm_context* norm;
    struct fp_im_tnode* root;
    struct runtime rt;
    struct runtime normalization;
    double rtime;
    double qtime;
    void* qr_min;
    void* qr_max;
    tree_test_config_t tc;
    int query;
    metadata_t* meta;
    dataset_t* data;
    dataset_t* orig_data;
    dataset_t* copy;
    int c_updates = 0;
    FILE* log;
    char logline[1024];
    char* loglinep = logline;
    void* queries;
    int total = 0;

    configure(&tc, argc, argv);
    if((f = fopen(tc.datfile, "r")) == NULL)
    {
        perror("Unable to open datafile:");
        exit(3);
    } 
   
    meta = metadata_read(f);
    
    /* create a new context */
    config.verbose = 1;
    config.normalize = tc.normalize ? 1 : 0;
    config.denormalize = tc.normalize ? 1 : 0;
    config.benchmark = 1;
    config.find_order = tc.constant ? CONSTANT : ITERATIVE;
    config.index = tc.forget ? FORGET : KEEP;
    config.cmp = tc.xyz ? XYZ : HILBERT;
    config.print = stdout;
    
    if(config.normalize)
    {
        norm = fp_create_norm_context(&config, meta->dimz, meta->dimf, 
                                      meta->start_order);
        orig_data = dataset_read(f, meta, norm->mc);
        /* the working context is the one for the normalized records */
        context = norm->zc;
        data = dataset_create(meta, context, N_RECORDS(orig_data));
        start_timing(&normalization);
        fp_normalize(norm, RECORDS(orig_data), N_RECORDS(orig_data), RECORDS(data));
        stop_timing(&normalization);
    }
    else
    {
        context = fp_create_context(&config, meta->dimz, meta->dimf, meta->start_order);
        data = dataset_read(f, meta, context);
        start_timing(&normalization);
        stop_timing(&normalization);
    }
    printf("Normalization time: %f\n", get_runtime(normalization));
    fclose(f);
    loglinep += sprintf(loglinep, "%d", N_RECORDS(data));
    loglinep += sprintf(loglinep, ",%d", DIMZ(data));
    loglinep += sprintf(loglinep, ",%d", DIMF(data));
    loglinep += sprintf(loglinep, ",%f", get_runtime(normalization));
    
    #if 0
    printf("Running sort benchmark");
    rtime = 0.0;
    for(i = 0; i < tc.n_sorts; i++)
    {
        copy = dataset_copy(data);
        results = NULL;
        start_timing(&rt);
        fp_im_sort(context, RECORDS(copy), N_RECORDS(copy), &results);
        stop_timing(&rt);
        rtime += get_runtime(rt);
        if(RECORDS(copy) != results)
        {
            printf("Sort failed!\n");
            exit(-1);
        }
        dataset_destroy(copy);
        putchar('.');
        fflush(stdout);
    }
    printf("\nAverage Sort time: %f\n", rtime / (double)tc.n_sorts);
    loglinep += sprintf(loglinep, ",%f", rtime / (double)tc.n_sorts);
    #endif

    results = NULL;
    printf("N is: %d\n", N_RECORDS(data));
    start_timing(&rt);
    fp_im_sort(context, RECORDS(data), N_RECORDS(data), &results);
    stop_timing(&rt);
    loglinep += sprintf(loglinep, ",%f", get_runtime(rt));
    srand(1000);
    printf("Done!\n");
    /* dataset_print(data); */
    
    if(tc.normalize)
    {
        /* reverse map records */
        fp_denormalize(norm, RECORDS(data), N_RECORDS(data), RECORDS(orig_data));
        context = norm->mc;
        dataset_destroy(data);
        data = orig_data;
    }
    dataset_verify(data);
    #if 0
    printf("Running tree build benchmark");
    rtime = 0.0;
    for(i = 0; i < tc.n_tree_builds; i++)
    {
        start_timing(&rt);
        root = fp_im_build_tree(context, N_RECORDS(data), RECORDS(data), LEAF);
        stop_timing(&rt);
        rtime += get_runtime(rt);
        fp_im_destroy_tree(root);
        putchar('.');
        fflush(stdout);
    }
    printf("\nAverage Build Tree Time: %f\n", rtime / (double)tc.n_tree_builds);
    loglinep += sprintf(loglinep, ",%f", rtime / (double)tc.n_tree_builds);
    #endif
    /* build a tree based on the sorted records */
    start_timing(&rt);
    root = fp_im_build_tree(context, N_RECORDS(data), RECORDS(data), LEAF);
    stop_timing(&rt);
    loglinep += sprintf(loglinep, ",%f", get_runtime(rt));
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
    qtime = 0.0;
    printf("Performing query benchmarks");
    rtime = 0.0;
    queries = malloc(2 * RECORD_S(data) * tc.n_queries);
    for(query = 0; query < tc.n_queries; query++)
    {
        int r1 = (int)((double)N_RECORDS(data) * ((double)rand() / (double)RAND_MAX));
        int r2 = (int)((double)N_RECORDS(data) * ((double)rand() / (double)RAND_MAX));
        qr_min = queries + (2 * RECORD_S(data) * query);
        qr_max = queries + (2 * RECORD_S(data) * query + RECORD_S(data));

        for(i = 0; i < context->env->dimz; i++)
        {
            RCOORDZ(qr_min, context, i) = MIN(COORDZ(data, r1, i), COORDZ(data, r2, i));
            RCOORDZ(qr_max, context, i) = MAX(COORDZ(data, r1, i), COORDZ(data, r2, i));
        }
        for(i = 0; i < context->env->dimf; i++)
        {
            RCOORDF(qr_min, context, i) = MIN(COORDF(data, r1, i), COORDF(data, r2, i));
            RCOORDF(qr_max, context, i) = MAX(COORDF(data, r1, i), COORDF(data, r2, i));
        }
        printf("Query region: ");
        print_record(stdout, context, qr_min);
        print_record(stdout, context, qr_max);
        printf("\n");
    }
    for(query = 0; query < tc.n_queries; query++)
    {
        qr_min = queries + (2 * RECORD_S(data) * query);
        qr_max = queries + (2 * RECORD_S(data) * query + RECORD_S(data));
        reported = 0;
        start_timing(&rt);
        i = fp_im_query_tree(context, root, qr_min, qr_max, report_quiet);
        stop_timing(&rt);
        total += reported;
        rtime += get_runtime(rt);
        if(query % 100 == 0)
        {
            putchar('.');
            fflush(stdout);
        }
    }
    printf("\nAverage query time: %f\n", rtime / (double)tc.n_queries);
    loglinep += sprintf(loglinep, ",%f", rtime / (double)tc.n_queries);
    printf("Total Number of records reported: %d\n", total);
    /*
    free(qr_min);
    free(qr_max);
    */
    fp_im_destroy_tree(root);

    while(c_updates < tc.n_updates)
    {
        FILE* u;
        metadata_t* umeta;
        dataset_t* udata;
        dataset_t* ucopy;
        dataset_t* new_data;
        double mtime;

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
        #if 0
        printf("Running merging benchmark");
        rtime = 0.0;
        mtime = 0.0;
        for(i = 0; i < tc.n_merges; i++)
        {
            ucopy = dataset_copy(udata);
            results = NULL;
            start_timing(&rt);
            fp_im_sort(context, RECORDS(ucopy), N_RECORDS(ucopy), &results);
            stop_timing(&rt);
            rtime += get_runtime(rt);
            copy = dataset_copy(data);

            start_timing(&rt);
            new_data = dataset_merge(copy, ucopy, TRUE);
            stop_timing(&rt);
            mtime += get_runtime(rt);
            
            dataset_destroy(ucopy);
            dataset_destroy(copy);
            dataset_destroy(new_data);
            putchar('.');
            fflush(stdout);
        }
        printf("\nAverage update merge time: %f\n", mtime / (double)tc.n_merges);
        loglinep += sprintf(loglinep, ",%f", mtime / (double)tc.n_merges);
        printf("Average update sort time: %f\n", rtime / (double)tc.n_merges);
        loglinep += sprintf(loglinep, ",%f", rtime / (double)tc.n_merges);
        #endif
        printf("N is: %d\n", N_RECORDS(udata));
        rtime = 0.0;
        results = NULL;
        start_timing(&rt);
        fp_im_sort(context, RECORDS(udata), N_RECORDS(udata), &results);
        stop_timing(&rt);
        rtime = get_runtime(rt);
        loglinep += sprintf(loglinep, ",%f", rtime);
        printf("Sort time: %f\n", rtime);
        start_timing(&rt);
        new_data = dataset_merge(data, udata, TRUE);
        stop_timing(&rt);
        loglinep += sprintf(loglinep, ",%f", get_runtime(rt));
        if(new_data == NULL)
        {
            printf("Error merging datasets.\n");
            return -1;
        }
        dataset_destroy(udata);
        dataset_destroy(data);
        data = new_data;
        c_updates++;
    }
    loglinep += sprintf(loglinep, ",%d", N_RECORDS(data) * RECORD_S(data));
    log = fopen(tc.logfile, "a");
    fprintf(log, "%s\n", logline);
    fclose(log);
    dataset_destroy(data);
    fp_destroy_context(context);
    metadata_destroy(meta);
    return 0;
} 
