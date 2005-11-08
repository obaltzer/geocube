#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>

struct configuration
{
    size_t dimz;
    size_t* cards;
    size_t dimf;
    size_t n;
    size_t prec;
    FILE* output;
};

int main(int argc, char** argv)
{
    struct configuration config = {0, NULL, 0, 0, 4, stdout};

    char* options = "f:z:n:p:o:";
    static struct option long_options[] = {
        {"output-file",     1, 0, 'o'},
        {"precision",       1, 0, 'p'},
        {"number",          1, 0, 'n'},
        {"dimz",            1, 0, 'z'},
        {"dimf",            1, 0, 'f'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    char c;
    char* t;
    int i;
    int d;
    
    while((c = getopt_long(argc, argv, options, 
                           long_options, &option_index)) != -1)
    {
        switch(c)
        {
            case 'o':
                if(!strcmp(optarg, "-"))
                    config.output = stdout;
                else if((config.output = fopen(optarg, "wt")) == NULL)
                {
                    perror("Unable to open output file: ");
                    exit(3);
                }
                break;

            case 'n':
                config.n = atoi(optarg);
                break;

            case 'f':
                config.dimf = atoi(optarg);
                break;

            case 'z':
                if(strcmp(optarg, "0"))
                {
                    i = 0;
                    t = optarg;
                    /* determine the number of discrete space dimensions */
                    while(*t != '\0')
                    {
                        t++;
                        if(*t == ',' || *t == '\0')
                            i++;
                    }
                    printf("dimz: %d\n", i);
                    config.dimz = i;
                    config.cards = 
                        (size_t*)malloc(sizeof(size_t) * config.dimz);
                    t = optarg;
                    i = 0;
                    while(*optarg != '\0')
                    {
                        if(*optarg == ',')
                        {
                            *optarg = '\0';
                            config.cards[i++] = atoi(t);
                            t = ++optarg;
                        }
                        else
                            optarg++;
                    }
                    config.cards[i] = atoi(t);
                }
                break;

            case 'p':
                config.prec = atoi(optarg);
                break;
                
            default:
                fprintf(stderr, "Unknow command line option: %s\n",
                        argv[optind]);
                exit(5);

        }
    }

    if(config.dimz + config.dimf < 2)
    {
        fprintf(stderr, "There are at least 2 dimensions required.\n");
        exit(6);
    }

    fprintf(config.output, "%d\n", config.dimz);
    for(i = 0; i < config.dimz; i++)
        fprintf(config.output, "%d ", config.cards[i]);
    fprintf(config.output, "\n");
    fprintf(config.output, "%d\n", config.dimf);
    fprintf(config.output, "%d\n", config.n);

    srand(time(NULL));
    for(i = 0; i < config.n; i++)
    {
        for(d = 0; d < config.dimz; d++)
            fprintf(config.output, "%d ", rand() % config.cards[d]);
        for(d = 0; d < config.dimf; d++)
            fprintf(config.output, "%.3f ", 
                    (double)rand() / (double)RAND_MAX);
        fprintf(config.output, "\n");
    }
    if(config.cards)
        free(config.cards);
    fclose(config.output);

    return -0;
}
