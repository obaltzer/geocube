#ifndef __UTIL_H
#define __UTIL_H

#include <stdio.h>
#include <sys/time.h>
#include <time.h>

#include "sort.h"

struct runtime
{
    struct timeval start_time;
    struct timeval end_time;
};

void print_record(FILE* f, struct fp_context* context, const void* r);
void print_record_mapped(FILE* f, struct fp_context* context, 
                         const void* r, int k);
void start_timing(struct runtime* timing);
void stop_timing(struct runtime* timing);
double get_runtime(struct runtime timing);

#endif
