#include <stdio.h>
#include <stdlib.h>
#ifndef WITH_FP
#include <gmp.h>
#endif
#include "hilbert.h"


int main(int argc, char** argv)
{
    FILE* f = NULL;
    int dimz;
    int dimf;
    int order;
    int i, j;
    int n;
    int coordz;
    float coordf;
#ifdef WITH_FP
    fpz_t* coordsz = NULL;
    fpf_t* coordsf = NULL;
    fpz_t index;
    struct fpm_env* e;
#else
    mpz_t* coordsz = NULL;
    mpf_t* coordsf = NULL;
    mpz_t index;
    struct mpm_env* e;
#endif
    
    if(argc < 2)
    {
        printf("Usage: hilbert_test datafile\n");
        exit(2);
    }
    
    if((f = fopen(argv[1], "r")) == NULL)
    {
        perror("Unable to open datafile:");
        exit(3);
    }
    
    fscanf(f, "%d", &dimz);
    fscanf(f, "%d\n", &dimf);
    fscanf(f, "%d\n", &order);
    fscanf(f, "%d\n", &n);
#ifdef WITH_FP    
    if(dimz)
        coordsz = (fpz_t*)malloc(dimz * sizeof(fpz_t));
    if(dimf)
        coordsf = (fpf_t*)malloc(dimf * sizeof(fpf_t));
#else
    if(dimz)
        coordsz = (mpz_t*)malloc(dimz * sizeof(mpz_t));
    if(dimf)
        coordsf = (mpf_t*)malloc(dimf * sizeof(mpf_t));
#endif

#ifdef WITH_FP
    e = fpm_create_env(dimz, dimf, 0);
#else
    for(i = 0; i < dimz; i++)
        mpz_init(coordsz[i]);
    for(i = 0; i < dimf; i++)
        mpf_init(coordsf[i]);
    mpz_init(index);
    e = mpm_create_env(dimz, dimf);
#endif
    
    for(i = 0; i < n; i++)
    {
        printf("(");
        for(j = 0; j < dimz; j++)
        {
            fscanf(f, "%d", &coordz);
#ifdef WITH_FP
            coordsz[j] = coordz;
#else
            mpz_set_ui(coordsz[j], coordz);
#endif
            
            if(j)
                printf(", %d", coordz);
            else
                printf("%d", coordz);
        }
        for(j = 0; j < dimf; j++)
        {
            fscanf(f, "%f", &coordf);
#ifdef WITH_FP
            coordsf[j] = (double)coordf;
#else
            mpf_set_d(coordsf[j], coordf);
#endif
            if((dimz + j))
                printf(", %f", coordf);
            else
                printf("%f", coordf);
        }
#ifdef WITH_FP
        fpm_c2i(e, order, coordsz, coordsf, &index);
        printf("): %d\n", (int)index);
#else
        mpm_c2i(e, order, coordsz, coordsf, &index);
        printf("): %d\n", (int)mpz_get_ui(index));
#endif
    }
#ifdef WITH_FP
    fpm_destroy_env(e);
#else
    mpm_destroy_env(e);
    for(i = 0; i < dimz; i++)
        mpz_clear(coordsz[i]);
    for(i = 0; i < dimf; i++)
        mpf_clear(coordsf[i]);
    mpz_clear(index);
#endif
    if(dimz)
        free(coordsz);
    if(dimf)
        free(coordsf);
    fclose(f);
    return 0;
}
