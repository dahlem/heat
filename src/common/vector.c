/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file vector.c
 * Implementation of the method declarations in vector.h. The dot product supports
 * the parallel calculation of the dot product in an MPI environment.
 *
 * @author Dominik Dahlem
 */
#if HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef NDEBUG
# include <assert.h>
#endif /* NDEBUG */

#ifdef HAVE_MPI
# include <mpi.h>
#endif /* HAVE_MPI */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "error.h"
#include "vector.h"


void vector_alloc(vector *vec, int len)
{
    vec->len = len;
    vec->data = malloc(sizeof(double) * len);
}

void vector_calloc(vector *vec, int len)
{
    vec->len = len;
    vec->data = calloc(len, sizeof(double));
}

void vector_free(vector *vec)
{
    free(vec->data);
}

void zero(const vector *vec)
{
    int i;

    for (i = 0; i < vec->len; ++i) {
        vec->data[i] = 0.0;
    }
}

double dotProduct(const vector *const a, const vector *const b)
{
    int i;
    double result;
#ifdef HAVE_MPI
    double global_result;
#endif /* HAVE_MPI */


#ifdef NDEBUG
    /* check the assumption that both vectors have the same size */
    assert(a->len == b->len);
#endif /* NDEBUG */

    result = 0.0;

    for (i = 0; i < a->len; ++i) {
        result += a->data[i] * b->data[i];
    }

#ifdef HAVE_MPI
    /* combine all local dotproducts */
    MPI_Allreduce(&result, &global_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return global_result;
#else
    return result;
#endif /* HAVE_MPI */
}

void daxpy(double alpha, const vector *const x, const vector *y)
{
    int i;

    for (i = 0; i < x->len; ++i) {
        y->data[i] = alpha * x->data[i] + y->data[i];
    }
}

double dnrm2(const vector *const x)
{
    int i;
    double result;

    result = 0.0;

    for (i = 0; i < x->len; ++i) {
        result += (x->data[i] * x->data[i]);
    }

    return sqrt(result);
}

int add(const vector *x, const vector *const y)
{
    int i;

    if (x->len != y->len) {
        return VECTOR_DIMENSION_MISMATCH;
    }

    for (i = 0; i < x->len; ++i) {
        x->data[i] += y->data[i];
    }

    return EXIT_SUCCESS;
}

void scale(double alpha, const vector *x)
{
    int i;

    for (i = 0; i < x->len; ++i) {
        x->data[i] *= alpha;
    }
}

void vector_copy(const vector *x, const vector *const y)
{
#ifdef NDEBUG
    assert(x->len >= y->len);
#endif /* NDEBUG */

    memcpy(x->data, y->data, y->len * sizeof(double));
}

void vector_print(const vector *const x)
{
    int i;

    for (i = 0; i < x->len; ++i) {
        if (i == (x->len - 1)) {
            fprintf(stdout, "%2.2f", x->data[i]);
        } else {
            fprintf(stdout, "%2.2f, ", x->data[i]);
        }
    }

    fprintf(stdout, "\n");
    fflush(stdout);
}
