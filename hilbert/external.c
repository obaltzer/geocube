#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#define __USE_GNU
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include "dataset.h"
#include "util.h"

struct stats_s
{
    size_t block_read;
    size_t block_write;
    size_t seek_forward;
    size_t seek_backward;
    size_t bytes_read;
    size_t bytes_written;
};
typedef struct stats_s stats_t;

stats_t gstats;

struct config_s
{
    size_t record_size;
    size_t block_size;
    size_t memory_size;
};
typedef struct config_s config_t;

struct stream_s
{
    char filename[256];
    int fd;
    off_t pos;
    config_t* config;
};
typedef struct stream_s stream_t;

struct external_dataset_s
{
    size_t n_records;
    metadata_t* meta;
    fp_context_t* context;
    stream_t* stream;
    dataset_t* mem;
};
typedef struct external_dataset_s external_dataset_t;

void stats_print()
{
    printf("Block read: %d\n", gstats.block_read);
    printf("Block write: %d\n", gstats.block_write);
    printf("Forward seek: %d\n", gstats.seek_forward);
    printf("Backward seek: %d\n", gstats.seek_backward);
}

#define TRUE 1
#define FALSE 0

#define BLOCK_SIZE(s) (s->config->block_size)
#define RECORD_SIZE(s)  (s->config->record_size)
#define BLOCK_RECORDS(s) (s->config->block_size / s->config->record_size)
#define BLOCK_PADDING(s) (s->config->block_size % s->config->record_size)
#define MEMORY_SIZE(s) (s->config->memory_size)
#define MEMORY_BLOCKS(s) (s->config->memory_size / (BLOCK_RECORDS(s) * RECORD_SIZE(s)))
#define MEMORY_RECORDS(s) (MEMORY_BLOCKS(s) * BLOCK_RECORDS(s))
#define BLOCKS(s, n) ((n / BLOCK_RECORDS(s)) + (n % BLOCK_RECORDS(s) ? 1 : 0))
#define MEMORY_BLOCK(s, m, i) (m + (BLOCK_RECORDS(s) * RECORD_SIZE(s) * i))
#define ACTUAL_MEMORY_SIZE(s) (MEMORY_RECORDS(s) * RECORD_SIZE(s) + BLOCK_PADDING(s)) 

int block_write(stream_t* stream, void* b)
{
    ssize_t n;
    gstats.block_write++;
    printf("Block: 0x%X\n", (unsigned int)b);
    n = write(stream->fd, b, BLOCK_SIZE(stream));
    stream->pos++;
    if(n == -1)
        perror("Error writing block: ");
    gstats.bytes_written += n;
    return 0;
}

int block_read(stream_t* stream, void* b)
{
    ssize_t n;
    gstats.block_read++;
    n = read(stream->fd, b, BLOCK_SIZE(stream));
    if(n == -1)
        perror("Error reading block: ");
    gstats.bytes_read += n;
    stream->pos++;
    return 0;
}

int memory_write(stream_t* stream, void* m, size_t n)
{
    int i = 0;
    printf("Memory: 0x%X\n", (unsigned int)m);
    for(; i < BLOCKS(stream, n); i++)
        block_write(stream, MEMORY_BLOCK(stream, m, i));
    return 0;
}

int memory_read(stream_t* stream, void* m, size_t n)
{
    int i = 0;
    for(; i < BLOCKS(stream, n); i++)
        block_read(stream, MEMORY_BLOCK(stream, m, i));
    return 0;
}

off_t stream_seek(stream_t* stream, off_t offset)
{
    off_t curr;

    offset *= BLOCK_SIZE(stream);
    curr = lseek(stream->fd, 0, SEEK_CUR);
    if(offset > curr)
        gstats.seek_forward++;
    else if(offset < curr)
        gstats.seek_backward++;
    curr = lseek(stream->fd, offset, SEEK_SET);
    stream->pos = curr / BLOCK_SIZE(stream);
    return stream->pos;
}

off_t stream_tell(stream_t* stream)
{
    return lseek(stream->fd, 0, SEEK_CUR) / BLOCK_SIZE(stream);
}

stream_t* stream_create(config_t* config, const char* filename)
{
    stream_t* stream;

    stream = malloc(sizeof(stream_t));
    if(stream == NULL)
    {
        printf("Cannot allocate memory for stream %s.\n", filename);
        return NULL;
    }
    strncpy(stream->filename, filename, sizeof(stream->filename));
    stream->fd = -1;
    stream->config = config;
    return stream;
}

int stream_open(stream_t* stream, int flags)
{
    if(stream->fd != -1)
    {
        stream_seek(stream, 0);
    }
    else if((stream->fd = open(stream->filename, flags, 0666)) == -1)
    {
        printf("Cannot open Direct I/O stream %s.\n", stream->filename);
        return FALSE;
    }
    stream->pos = 0;
    return TRUE;
}

void stream_close(stream_t* stream)
{
    if(stream->fd != -1)
        close(stream->fd);
    stream->fd = -1;
}

void stream_destroy(stream_t* stream)
{
    stream_close(stream);
    free(stream);
}

void dataset_convert(dataset_t* data, stream_t* out)
{
    void* m;
    int i;

    printf("Allocating %d bytes.\n", ACTUAL_MEMORY_SIZE(out));
    m = malloc(ACTUAL_MEMORY_SIZE(out));
    if(m == NULL)
    {
        printf("Cannot allocate memory.\n");
        return;
    }
    memset(m, 0, ACTUAL_MEMORY_SIZE(out));

    for(i = 0; i < data->n_records / MEMORY_RECORDS(out); i++)
    {
        memcpy(m, data->records + (i * MEMORY_RECORDS(out) * RECORD_SIZE(out)), MEMORY_RECORDS(out) * RECORD_SIZE(out));
        memory_write(out, m, MEMORY_RECORDS(out));
        stats_print();
    }
    if(data->n_records % MEMORY_RECORDS(out))
    {
        memcpy(m, data->records + (i * MEMORY_RECORDS(out) * RECORD_SIZE(out)), (data->n_records % MEMORY_RECORDS(out)) * RECORD_SIZE(out));
        memory_write(out, m, data->n_records % MEMORY_RECORDS(out));
        stats_print();
    }
    free(m);
}

void external_dataset_destroy(external_dataset_t* dataset)
{
    if(dataset)
    {
        if(dataset->mem)
        {
            if(dataset->mem->records)
                free(dataset->mem->records);
            free(dataset->mem);
        }
        if(dataset->context)
            fp_destroy_context(dataset->context);
        free(dataset);
    }
}

external_dataset_t* external_dataset_create(stream_t* stream, metadata_t* meta, size_t n)
{
    external_dataset_t* data;
    struct sort_config sc;

    data = malloc(sizeof(external_dataset_t));
    if(data == NULL)
    {
        printf("Cannot allocate memory for external dataset.\n");
        return NULL;
    }
    memset(data, 0, sizeof(external_dataset_t));
    data->meta = meta;
    sc.verbose = 1;
    sc.normalize = 0;
    sc.denormalize = 0;
    sc.benchmark = 1;
    sc.find_order = ITERATIVE;
    sc.index = KEEP;
    sc.cmp = HILBERT;
    sc.print = stdout;
    data->n_records = n;
    data->stream = stream;
    data->context = fp_create_context(&sc, meta->dimz, meta->dimf, meta->start_order);
    data->mem = malloc(sizeof(dataset_t));
    if(data->mem == NULL)
    {
        printf("Cannot allocate memory for in-memory dataset.\n");
        external_dataset_destroy(data);
        return NULL;
    }
    memset(data->mem, 0, sizeof(dataset_t));
    data->mem->meta = meta;
    data->mem->context = data->context;
    data->mem->n_records = 0;
    data->mem->records = malloc(ACTUAL_MEMORY_SIZE(stream));
    if(data->mem->records == NULL)
    {
        printf("Cannot allocate memory for records.\n");
        external_dataset_destroy(data);
        return NULL;
    }
    return data;
}

void external_dataset_sort(external_dataset_t* data)
{
    struct runtime rt;
    void* results = NULL;
    size_t i;
    
    stream_t* tmp;
    stream_t* stream = data->stream;
    stream_seek(stream, 0);
    tmp = stream_create(stream->config, "tmp.stream");
    stream_open(tmp, O_CREAT | O_TRUNC | O_RDWR | O_SYNC);
    /* sort each individual chunk of memory size */
    for(i = 0; i < data->n_records / MEMORY_RECORDS(stream); i++)
    {
        memory_read(stream, RECORDS(data->mem), MEMORY_RECORDS(stream));
        N_RECORDS(data->mem) = MEMORY_RECORDS(stream);
        results = NULL;
        start_timing(&rt);
        fp_im_sort(data->context, RECORDS(data->mem), N_RECORDS(data->mem), &results);
        stop_timing(&rt);
        printf("Sort time: %f\n", get_runtime(rt));
        memory_write(tmp, RECORDS(data->mem), N_RECORDS(data->mem));
        dataset_print(data->mem, TRUE);

    }
    if(data->n_records % MEMORY_RECORDS(stream))
    {
        memory_read(stream, RECORDS(data->mem), data->n_records % MEMORY_RECORDS(stream));
        N_RECORDS(data->mem) = data->n_records % MEMORY_RECORDS(stream);
        start_timing(&rt);
        results = NULL;
        fp_im_sort(data->context, RECORDS(data->mem), N_RECORDS(data->mem), &results);
        stop_timing(&rt);
        memory_write(tmp, RECORDS(data->mem), N_RECORDS(data->mem));
        dataset_print(data->mem, TRUE);
    }
    /* merge memory chunks block by block */
}

int main(int argc, char** argv)
{
    FILE* f;
    metadata_t* meta;
    fp_context_t* context;
    struct sort_config sc;
    dataset_t* ds;
    char target[256];
    stream_t* stream;
    config_t conf;
    external_dataset_t* ext;
    int i; 

    if(argc < 2)
    {
        printf("Barf!\n");
        exit(-1);
    }
    
    if((f = fopen(argv[1], "r")) == NULL)
    {
        printf("Cannot open dataset file %s.\n", argv[1]);
        return -1;
    }
    meta = metadata_read(f);
    sc.verbose = 1;
    sc.normalize = 0;
    sc.denormalize = 0;
    sc.benchmark = 1;
    sc.find_order = ITERATIVE;
    sc.index = KEEP;
    sc.cmp = HILBERT;
    sc.print = stdout;
    context = fp_create_context(&sc, meta->dimz, meta->dimf, meta->start_order);
    ds = dataset_read(f, meta, context);
    fclose(f);
    dataset_print(ds, FALSE);
    
    strcpy(target, argv[1]);
    strcat(target, ".stream");
    printf("Number of records: %d\n", ds->n_records);
    printf("Record size: %d\n", context->record_size);

    conf.memory_size = 1000; 
    conf.block_size = 0x100;
    conf.record_size = context->record_size;
    stream = stream_create(&conf, target);
    stream_open(stream, O_CREAT | O_TRUNC | O_SYNC | O_WRONLY);
    dataset_convert(ds, stream);
    stream_close(stream);
    stream_open(stream, O_SYNC | O_RDONLY);
    ext = external_dataset_create(stream, meta, ds->n_records);

    for(i = 0; i < ds->n_records / MEMORY_RECORDS(stream); i++)
    {
        memory_read(stream, ext->mem->records, MEMORY_RECORDS(stream));
        ext->mem->n_records = MEMORY_RECORDS(stream);
        dataset_print(ext->mem, FALSE);
    }
    if(ds->n_records % MEMORY_RECORDS(stream))
    {
        memory_read(stream, ext->mem->records, ds->n_records % MEMORY_RECORDS(stream));
        ext->mem->n_records = ds->n_records % MEMORY_RECORDS(stream);
        dataset_print(ext->mem, FALSE);
    }
    printf("Stream position: %ld vs. %ld\n", stream->pos, stream_tell(stream));
    external_dataset_sort(ext);
    
    
    external_dataset_destroy(ext);
    stream_destroy(stream);
    stats_print();
    dataset_destroy(ds);
    fp_destroy_context(context);
    metadata_destroy(meta);
    return 0;
}
