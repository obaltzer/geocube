#include "hilbert.h"

#include <stdlib.h>
#include <string.h>

#ifndef WITHOUT_FP

/**
 * Creates a new enviroment for fixed precision operation.
 * 
 * @param[in] dims number of dimensions in the view
 * @return A pointer to the environment or NULL if the memory could not be
 *         allocated.
 */
struct fpz_env* fpz_create_env(int dims)
{
    struct fpz_env* ret = malloc(sizeof(struct fpz_env));
    if(ret)
        ret->dims = dims;
    return ret;
}

/**
 * Destroys the previously created operation environment for fixed
 * precision and releases all allocated memory.
 *
 * @param[in] env a pointer to the operation environment
 */
void fpz_destroy_env(struct fpz_env* env)
{
    if(env)
        free(env);
}

/**
 * Computes the Hilber index for the point in discrete space given by the
 * coordinates coords and at the order k.
 *
 * @param[in] env a pointer to the fixed precision operation environment
 *                that should be used
 * @param[in] k the order at which the Hilbert code should be computed
 * @param[in] coords the coordinates of the point at the given order
 * @param[out] index the Hilbert index of the point that has been computed
 */
void fpz_c2i(struct fpz_env* env, int k, fpz_t coords[], fpz_t* index)
{
    /* this code has been taken from Doug Moore's implementation */
    fpz_t const one = 1;
    fpz_t const ndOnes = (one << env->dims) - 1;
    fpz_t const nthbits = (((one << env->dims * k) - one) 
                                / ndOnes) >> 1;
    int b, d;
    int rotation = 0; 
    fpz_t reflection = 0;
    *index = 0;
    
    for(b = k; b--;)
    {
        fpz_t bits = reflection;
        reflection = 0;
        for(d = 0; d < env->dims; d++)
            reflection |= ((coords[d] >> b) & 1 ) << d;
        bits ^= reflection;
        
        /* expanded rotationRight macro */ 
        bits = ((bits >> rotation) | (bits << (env->dims - rotation))) 
                    & ((1 << env->dims) - 1);
        
        *index |= bits << (env->dims * b);
        reflection ^= one << rotation;
        
        /* expanded adjust_rotation macro */
        bits &= -bits & ((1 << (env->dims - 1)) - 1);
        while(bits)
            bits >>= 1, ++rotation;
        if(++rotation >= env->dims)
            rotation -= env->dims;
    }
    *index ^= nthbits;
  
    /* take care of 32bit overflow */
    for(d = 1; ; d *= 2) 
    {
        fpz_t t;
        if(d <= 32)
        {
            t = *index >> d;
            if(!t)
                break;
        }
        else 
        {
            t = *index >> 32;
            t = t >> (d - 32);
            if(!t)
                break;
        }
        *index ^= t;
    }
}

/**
 * Compares two points in discrete space with respect to their location on
 * the Hilbert curve.
 *
 * @param[in] env a pointer to the fixed precision operation environment
 *                that should be used
 * @param[in] k the order at which the Hilbert codes should be compared
 * @param[in] c1 the coordinates of the first point at the given order
 * @param[in] c2 the coordinates of the first point at the given order
 * @return -1 if the first point is before the second in Hilbert order, 0
 *         if both points are at the same index, 1 if the first points
 *         comes after the second in Hilbert order
 */
int fpz_hcmp(struct fpz_env* env, int k, char* c1, char* c2)
{
    /* this code has been taken from Doug Moore's implementation */
    fpz_t const one = 1;
    int y, b, d;
    int rotation = 0;
    fpz_t reflection = 0;
    fpz_t index = 0;

    for(y = 0; y < sizeof(fpz_t); ++y)
    {
        for(b = 8; b >= 0; b--)
        {
            fpz_t bits = reflection;
            fpz_t bts2 = bits;
            fpz_t reflection2 = 0;
            reflection = 0;
            for(d = 0; d < env->dims; d++)
            {
                reflection |= 
                    ((((char*)c1)[y + d * sizeof(fpz_t)] >> b) & 1) << d;
                reflection2 |=
                    ((((char*)c2)[y + d * sizeof(fpz_t)] >> b) & 1) << d;
            }
            bits ^= reflection;
            /* expanded rotationRight macro */ 
            bits = ((bits >> rotation) | (bits << (env->dims - rotation))) 
                    & ((1 << env->dims) - 1);
            if(reflection != reflection2)
            {
                bts2 ^= reflection2;
                bts2 = ((bts2 >> rotation) | (bts2 << (env->dims - rotation))) 
                        & ((1 << env->dims) - 1);
                for(d = 1; d < env->dims; d <<= 1)
                {
                    index ^= index >> d;
                    bits ^= bits >> d;
                    bts2 ^= bts2 >> d;
                }
                return ((index ^ b) & 1) == (bits < bts2) ? -1 : 1;
            }
            index ^= bits;
            reflection ^= one << rotation;
            /* expansion of adjust rotation macro */
            bits &= -bits & ((1 << (env->dims - 1)) - 1);
            while(bits)
                bits >>= 1, ++rotation;
            if(++rotation >= env->dims)
                rotation -= env->dims;
        }
    }
    return 0;
}

/**
 * Creates a new computation environment to perform fixed precision
 * calculation on records with discrete and continuous space dimensions.
 * 
 * @param[in] dimz the number of discrete space dimensions
 * @param[in] dimf the number of continuous space dimensions
 * @param[in] base_order base order of discrete space for skew compensation
 * @return a reference to an appropriately allocated computation
 *         environment
 */
struct fpm_env* fpm_create_env(int dimz, int dimf, int base_order)
{
    struct fpm_env* env = (struct fpm_env*)malloc(sizeof(struct fpm_env));
    if(env)
    {
        env->dimz = dimz;
        env->dimf = dimf;
        env->envz = fpz_create_env(dimz + dimf);
        env->base_order = base_order;
        env->calls = 0;
        /* Only if there are continuous space dimensions defined we need to
         * have a working array. Otherwise we can just use the one for the
         * discrete space dimensions. */
        if(dimf)
        {
            env->coords1 = (fpz_t*)malloc(sizeof(fpz_t) * (dimz + dimf));
            env->coords2 = (fpz_t*)malloc(sizeof(fpz_t) * (dimz + dimf));
        }
        else
            env->coords1 = env->coords2 = NULL;
    }
    return env;
}

/**
 * Destroys a formally created computation environment for mixed dimensions
 * and frees all allocated memory.
 * 
 * @param[in] env a reference to the environment that should be destroyed.
 */
void fpm_destroy_env(struct fpm_env* env)
{
    if(env)
    {
        fpz_destroy_env(env->envz);
        if(env->coords1)
            free(env->coords1);
        if(env->coords2)
            free(env->coords2);
        free(env);
    }
}

/**
 * Computes the Hilbert index of a coordinate in multi-dimensional space
 * with mixed discrete and continuous dimensions.
 *
 * @param[in] env the computation environment to use
 * @param[in] k the order of the Hilbert curve to be used to compute the
 *              index
 * @param[in] coordsz the array of coordinates in the discrete space
 * @param[in] coordsf array of normalized coordinates in continuous space
 * @param[in,out] index reference to the variable which will hold the
 *                      result of the computation
 */
void fpm_c2i(struct fpm_env* env, int k, fpz_t coordsz[], fpf_t coordsf[],
             fpz_t* index)
{
    env->calls++;
    /* Only if there are continuous space dimensions defined we need to
     * perform extra computation to convert them into discrete space. */
    /* if(env->dimf) */
    if(env->dimf)
    {
        int i = 0;
        fpz_t bits = 1;

        bits = bits << k;
        
        for(i = 0; i < env->dimz; i++)
            /* compensate for discrete dimension skew */
            env->coords1[i] = coordsz[i] << (k - env->base_order);
            
        for(i = 0; i < env->dimf; i++)
        {
            /* convert continuous into discrete space, while handling
             * closed intervals gracefully */
            env->coords1[i + env->dimz] = 
                (fpz_t)(coordsf[i] * (fpf_t)bits) == bits 
                    ? bits - 1
                    : (fpz_t)(coordsf[i] * (fpf_t)bits);
        }
        
        fpz_c2i(env->envz, k, env->coords1, index);
    }
    else
        fpz_c2i(env->envz, k, coordsz, index);
}
    
/**
 * Compares two points with respect to their position on the Hilbert curve
 * in multi-dimensional space with mixed discrete and continuous space
 * dimensions. 
 *
 * @param[in] env the computation environment to use
 * @param[in] k the order of the Hilbert curve to be used to compare the
 *              indices
 * @param[in] cz1 the array of coordinates of point 1 in the discrete space
 * @param[in] cf1 array of normalized coordinates of point 1 in continuous
 *                space 
 * @param[in] cz2 the array of coordinates of point 2 in the discrete space
 * @param[in] cf2 array of normalized coordinates of point 2 in continuous
 *                space
 * @return -1 if the first point lies before the second point, 0 if both
 *         points share the same index, 1 if the second point lies before
 *         the first point
 */
int fpm_hcmp(struct fpm_env* env, int k, fpz_t cz1[], fpf_t cf1[],
              fpz_t cz2[], fpf_t cf2[])
{
    env->calls++;
    /* Only if there are continuous space dimensions defined we need to
     * perform extra computation to convert them into discrete space. */
    /* if(env->dimf) */
    if(env->dimf)
    {
        int i = 0;
        fpz_t bits = 1;

        bits = bits << k;
        
        for(i = 0; i < env->dimz; i++)
        {
            /* compensate for discrete dimension skew */
            env->coords1[i] = cz1[i] << (k - env->base_order);
            env->coords2[i] = cz2[i] << (k - env->base_order);
        }
        
        for(i = 0; i < env->dimf; i++)
        {
            /* convert continuous into discrete space, while handling
             * closed intervals gracefully */
            env->coords1[i + env->dimz] = 
                (fpz_t)(cf1[i] * (fpf_t)bits) == bits 
                    ? bits - 1
                    : (fpz_t)(cf1[i] * (fpf_t)bits);
            env->coords2[i + env->dimz] = 
                (fpz_t)(cf2[i] * (fpf_t)bits) == bits 
                    ? bits - 1
                    : (fpz_t)(cf2[i] * (fpf_t)bits);
        }
        
        return fpz_hcmp(env->envz, k, (char*)env->coords1, 
                        (char*)env->coords2);
    }
    else
        return fpz_hcmp(env->envz, k, (char*)cz1, (char*)cz2);
}
#endif

#ifndef WITHOUT_MP

/**
 * Creates a new environment for computation with multiple precision. It will
 * make all necessary memory allocations for all variables to speed up
 * computation in subsequent calls to mapping functions.
 *
 * @param[in] dims The number of dimensions used for this computation.
 * @return Pointer to the environment or NULL on failue.
 */
struct mpz_env* mpz_create_env(int dims)
{
    struct mpz_env* e = NULL;

    e = (struct mpz_env*)malloc(sizeof(struct mpz_env));
    if(e)
    {
        e->dims = dims;
        mpz_init(e->tmp1);
        mpz_init(e->tmp2);
        
        /* one = 1 */
        mpz_init_set_ui(e->one, 1);

        /* ndOnes = (1 << dims) - 1 */
        mpz_init(e->ndOnes);
        mpz_mul_2exp(e->ndOnes, e->one, dims);
        mpz_sub_ui(e->ndOnes, e->ndOnes, 1);

        mpz_init(e->nthbits);
        mpz_init(e->reflection);
        mpz_init(e->bits);
    }
    return e;
}

/**
 * Destroyes an existing environment and frees all allocated memory.
 * 
 * @param[in] env Pointer to the environment that is to be destroyed.
 */
void mpz_destroy_env(struct mpz_env* env)
{
    if(env)
    {
        mpz_clear(env->bits);
        mpz_clear(env->reflection);
        mpz_clear(env->nthbits);
        mpz_clear(env->ndOnes);
        mpz_clear(env->one);
        mpz_clear(env->tmp2);
        mpz_clear(env->tmp1);
        free(env);
    }
    return;
}

/**
 * Maps the given coordinates in the multi-dimensional space of order k to
 * a Hilbert code.
 *
 * @param[in] env Pointer to the environment that is used for computation.
 * @param[in] k The order of the grid on which the Hilbert code is
 *              computed.
 * @param[in] coords Pointer to an array of coordinates on the
 *                   multi-dimensional grid.
 * @param[out] index Pointer to a variable that holds the result index.  
 */
void mpz_c2i(struct mpz_env* env, int k, mpz_t coords[], mpz_t* index)
{
    int b, d;
    int rotation = 0; /* or (k * (dims -1)) % dims; */

    /* nthbits = (((1 << (nDims * nBits)) - 1) / ndOnes) - 1 */
    mpz_set(env->nthbits, env->one);
    mpz_mul_2exp(env->nthbits, env->nthbits, env->dims * k);
    mpz_sub(env->nthbits, env->nthbits, env->one);
    mpz_fdiv_q(env->nthbits, env->nthbits, env->ndOnes);
    mpz_tdiv_q_2exp(env->nthbits, env->nthbits, 1);
   
    /* reflection = 0 */
    mpz_set_ui(env->reflection, 0);
    
    /* index = 0 */
    mpz_set_ui(*index, 0);
    
    for(b = k; b--;)
    {
        /* bits = reflection */
        mpz_set(env->bits, env->reflection);
        /* reflection = 0 */
        mpz_set_ui(env->reflection, 0);
        
        for(d = 0; d < env->dims; d++)
        {
            /* reflection |= ((coord[d] >> b) & 1) << d */
            mpz_tdiv_q_2exp(env->tmp1, coords[d], b);
            mpz_and(env->tmp1, env->tmp1, env->one);
            mpz_mul_2exp(env->tmp1, env->tmp1, d);
            mpz_ior(env->reflection, env->reflection, env->tmp1);
        }
        /* bits ^= reflection */
        mpz_xor(env->bits, env->bits, env->reflection);
        
        /* expanded rotateRight(bits, bits, rotation, nDims): */
        mpz_tdiv_q_2exp(env->tmp1, env->bits, rotation);
        mpz_mul_2exp(env->tmp2, env->bits, env->dims - rotation);
        mpz_ior(env->tmp1, env->tmp1, env->tmp2);
        mpz_set_ui(env->tmp2, (1 << env->dims) - 1);
        mpz_and(env->bits, env->tmp1, env->tmp2);
          
        /* index |= bits << nDims * b */
        mpz_mul_2exp(env->tmp1, env->bits, env->dims * b);
        mpz_ior(*index, *index, env->tmp1);
        
        /* reflection ^= one << rotation */
        mpz_mul_2exp(env->tmp1, env->one, rotation);
        mpz_xor(env->reflection, env->reflection, env->tmp1);
            
        /* expand adjust_rotation(rotation, nDims, bits) macro: */
        mpz_mul_2exp(env->tmp1, env->one, env->dims - 1);
        mpz_sub_ui(env->tmp1, env->tmp1, 1);
        mpz_neg(env->tmp2, env->bits);
        mpz_and(env->tmp2, env->tmp2, env->tmp1);
        mpz_and(env->bits, env->bits, env->tmp2);
        /* while bits != 0 */
        while(mpz_sgn(env->bits) != 0)
        {
            mpz_tdiv_q_2exp(env->bits, env->bits, 1);
            ++rotation;
        }
        if(++rotation >= env->dims)
            rotation -= env->dims;
    }
    /* index ^= nthbits */
    mpz_xor(*index, *index, env->nthbits);
    
    for(d = 1; ; d *= 2) 
    {
        /* We do not have to worry about the 32 bits shift registers
         * since we are using GMP */
        /* t = index >> d */
        mpz_tdiv_q_2exp(env->tmp1, *index, d);
        if(!mpz_sgn(env->tmp1))
            break;
        /* index ^= t */ 
        mpz_xor(*index, *index, env->tmp1);
    }
}

/**
 * Creates a new computation environment to perform multiple precision
 * calculation on records with discrete and continuous space dimensions.
 * 
 * @param[in] dimz the number of discrete space dimensions
 * @param[in[ dimf the number of continuous space dimensions
 * @return a reference to an appropriately allocated computation
 *         environment
 */
struct mpm_env* mpm_create_env(int dimz, int dimf)
{
    struct mpm_env* env = (struct mpm_env*)malloc(sizeof(struct mpm_env));
    if(env)
    {
        env->dimz = dimz;
        env->dimf = dimf;
        env->envz = mpz_create_env(dimz + dimf);
        /* Only if there are continuous space dimensions defined we need to
         * initialize all the other variables associated with this
         * computation environment. */
        if(dimf)
        {
            int i = 0;
            env->coords = (mpz_t*)malloc(sizeof(mpz_t) * (dimf + dimz));
            /* We only need to initialize the last dimf coordinates,
             * because all coordinates from dimz will be copied. */
            for(; i < dimf; i++)
                mpz_init(env->coords[dimz + i]);
            mpz_init(env->bits);
            mpf_init(env->bits_f);
            mpf_init(env->tmp1);
            mpf_init(env->tmp2);
        }
        else
            env->coords = NULL;
    }
    return env;
}

/**
 * Destroys a formally created computation environment for mixed dimensions
 * and frees all allocated memory.
 * 
 * @param[in] env a reference to the environment that should be destroyed.
 */
void mpm_destroy_env(struct mpm_env* env)
{
    if(env)
    {
        mpz_destroy_env(env->envz);
        if(env->dimf)
        {
            int i = 0;
            /* Only clear those coordinates which actually had been
             * allocated. */
            for(; i < env->dimf; i++)
                mpz_clear(env->coords[env->dimz + i]); 
            
            free(env->coords);
            mpz_clear(env->bits);
            mpf_clear(env->bits_f);
            mpf_clear(env->tmp1);
            mpf_clear(env->tmp2);
        }
        free(env);
    }
}

/**
 * Computes the Hilbert index of a coordinate in multi-dimensional space
 * with mixed discrete and continuous dimensions.
 *
 * @param[in] env the computation environment to use
 * @param[in] k the order of the Hilbert curve to be used to compute the
 *              index
 * @param[in] coordsz the array of coordinates in the discrete space
 * @param[in] coordsf array of normalized coordinates in continuous space
 * @param[in,out] index reference to the variable which will hold the
 *                      result of the computation
 */
void mpm_c2i(struct mpm_env* env, int k, mpz_t coordsz[], mpf_t coordsf[],
             mpz_t* index)
{
    /* Only if there are continuous space dimensions defined we need to
     * perform extra computation to convert them into discrete space. */
    if(env->dimf)
    {
        int i = 0;
        /* Copy already discrete dimension coordinates into the working
         * array. */
        memcpy(env->coords, coordsz, sizeof(mpz_t) * env->dimz);
        
        /* Compute the number of bits per dimensions. */
        mpz_set_ui(env->bits, 1);
        mpz_mul_2exp(env->bits, env->bits, k);
        mpf_set_z(env->bits_f, env->bits);
        
        for(; i < env->dimf; i++)
        {
            mpf_mul(env->tmp1, coordsf[i], env->bits_f);
            mpz_set_f(env->coords[i + env->dimz], env->tmp1);
        }
        mpz_c2i(env->envz, k, env->coords, index);
    }
    else
        mpz_c2i(env->envz, k, coordsz, index);
}
 #endif
