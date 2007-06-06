#include <stdio.h>
#include <sys/time.h>
#include <time.h>

#include "util.h"

void print_record(FILE* f, struct fp_context* context, const void* r)
{
    int i;
    fpz_t coordz;
    fpf_t coordf;
    
    fprintf(f, "(");
    for(i = 0; i < context->env->dimz; i++)
    {
        coordz = ((fpz_t*)(r + context->coordsz_off))[i];
        if(i)
            fprintf(f, ", %llu", coordz);
        else
            fprintf(f, "%llu", coordz);
    }
    for(i = 0; i < context->env->dimf; i++)
    {
        coordf = ((fpf_t*)(r + context->coordsf_off))[i];
        if((context->env->dimz + i))
            fprintf(f, ", %f", coordf);
        else
            fprintf(f, "%f", coordf);
    }
    fprintf(f, ")");
}

void print_record_mapped(FILE* f, struct fp_context* context, 
                         const void* r, int k)
{
    int i;
    fpz_t coordz;
    fpf_t coordf;
    fpz_t bits = 1;
    bits = bits << k;
    
    fprintf(f, "(%d, %llu: ", k, bits);
    for(i = 0; i < context->env->dimz; i++)
    {
        coordz = ((fpz_t*)(r + context->coordsz_off))[i];
        if(i)
            fprintf(f, ", %llu", coordz);
        else
            fprintf(f, "%llu", coordz);
    }
    for(i = 0; i < context->env->dimf; i++)
    {
        coordf = ((fpf_t*)(r + context->coordsf_off))[i];
        if((context->env->dimz + i))
            fprintf(f, ", %f:%llu", coordf, (fpz_t)(coordf * (fpf_t)bits));
        else
            fprintf(f, "%f:%llu", coordf, (fpz_t)(coordf * (fpf_t)bits));
    }
    fprintf(f, ")");
}

void start_timing(struct runtime* timing)
{
    gettimeofday(&timing->start_time, NULL);
}

void stop_timing(struct runtime* timing)
{
    gettimeofday(&timing->end_time, NULL);
}

double get_runtime(struct runtime timing)
{
    return ((double)timing.end_time.tv_sec 
                    + ((double)timing.end_time.tv_usec / 1000000.0)) 
                - ((double)timing.start_time.tv_sec 
                    + ((double)timing.start_time.tv_usec / 1000000.0));
}
