/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */
#include <assert.h>
#include <math.h>
#include <stdlib.h>

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

    /* check the assumption that both vectors have the same size */
    assert(a->len == b->len);

    result = 0.0;
    
    for (i = 0; i < a->len; ++i) {
        result += a->data[i] * b->data[i];
    }

    return result;
}

void daxpy(double alpha, const vector *const x, const vector *y)
{
    int i;

    for (i = 0; i < x->len; ++i) {
        y->data[i] = alpha * x->data[i];
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
