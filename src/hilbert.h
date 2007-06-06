/**
 * @file hilbert.h
 *
 * This header file defines a set of functions that are used to compute
 * the k-th order Hilbert code of a point in n-dimensional space. The space
 * can be categorical and/or continuous space, where the continuous space
 * is subdivided into 2^k equally sized intervals in each dimension.
 *
 * The file contains two implementations of the functions. One for fixed
 * precision integer and floating point arithmetics and for
 * multiple/arbitrary precision integer and floating point arithmetics.
 */
#ifndef __HILBERT_H
#define __HILBERT_H

#ifndef WITHOUT_MP
#include <gmp.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifndef WITHOUT_FP
/* Those are a bunch of type definitions for fixed precision arithmetics */

/** The integer data type */
typedef unsigned long long fpz_t;
/** Floating point data type */
typedef double fpf_t;

/**
 * Computation environment for fixed precision calculation.
 */
struct fpz_env
{
    /** number of dimensions in the space */
    int dims;
};

struct fpz_env* fpz_create_env(int dims);
void fpz_destroy_env(struct fpz_env* env);
void fpz_c2i(struct fpz_env* env, int k, fpz_t coords[], fpz_t* index);
int fpz_hcmp(struct fpz_env* env, int k, char* c1, char* c2);

/**
 * Computation environment for fixed point computation but with mixed
 * dimension types, supporting discrete space and continuous space
 * dimensions.
 */
struct fpm_env
{
    /** reference to the discrete space environment */
    struct fpz_env* envz;
    
    /** number of discrete space dimensions */
    int dimz;

    /** number of continuous space dimensions */
    int dimf;

    /** the array which will store the discretized coordinates */
    fpz_t* coords1;
    fpz_t* coords2;

    /** the base order of the discrete space for skew compensation */
    int base_order;

    /** counts the number of calls */
    unsigned long long int calls;
};

struct fpm_env* fpm_create_env(int dimz, int dimf, int base_order);
void fpm_destroy_env(struct fpm_env* env);
void fpm_c2i(struct fpm_env* env, int k, fpz_t coordsz[], fpf_t coordsf[],
             fpz_t* index);
int fpm_hcmp(struct fpm_env* env, int k, fpz_t cz1[], fpf_t cf1[], 
              fpz_t cz2[], fpf_t cf2[]);
#endif

#ifndef WITHOUT_MP
/* We do not need type definitions for multiple precision since GMP already
 * supports our naming scheme */

/**
 * Computation environment for multiple precision calculation.
 */
struct mpz_env
{
    /** number of dimensions in the space */
    int dims;

    /* temporary variables used in the calculation */
    mpz_t tmp1;
    mpz_t tmp2;
    mpz_t one;
    mpz_t ndOnes;
    mpz_t nthbits;
    mpz_t reflection;
    mpz_t bits;
};

struct mpz_env* mpz_create_env(int dims);
void mpz_destroy_env(struct mpz_env* env);
void mpz_c2i(struct mpz_env* env, int k, mpz_t coords[], mpz_t* index);

/**
 * Computation environment for multiple precision calculation of mixed
 * dimensions.
 */
struct mpm_env
{
    /** Reference to the discrete space computation environment. */
    struct mpz_env* envz;

    /** Number of discrete space dimensions. */
    int dimz;

    /** Number of continuous space dimensions. */
    int dimf;

    /** Working array to hold the combined discretized coordinates. */
    mpz_t* coords;

    /** Pre-allocated temporary variable to store the number of bits in
     * each dimensions. */
    mpz_t bits;
    mpf_t bits_f;

    /** Other pre-allocated temporary variables. */
    mpf_t tmp1;
    mpf_t tmp2;
};

struct mpm_env* mpm_create_env(int dimz, int dimf);
void mpm_destroy_env(struct mpm_env* env);
void mpm_c2i(struct mpm_env* env, int k, mpz_t coordsz[], mpf_t coordsf[],
             mpz_t* index);
#endif

#endif /* hilbert.h */
