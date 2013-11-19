/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@gmail.com>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */
#include <stdio.h>
#include <stdlib.h>

#include <CUnit/CUnit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

#include "matrix.h"
#include "vector.h"

#include "vector_test.h"


void registerVectorTests()
{
    if (CU_register_suites(vector_suites) != CUE_SUCCESS) {
        fprintf(stderr, "Suite registration failed - %s\n", CU_get_error_msg());
        exit(CU_get_error());
    }
}

void testDotProduct()
{
    vector a, b;
    gsl_vector *gsl_a, *gsl_b;
    double dotprod;
    int i;

    vector_calloc(&a, 3);
    vector_calloc(&b, 3);

    for (i = 0; i < 3; ++i) {
        a.data[i] = i;
        b.data[i] = 2 - i;
    }

    gsl_a = gsl_vector_alloc(3);
    gsl_b = gsl_vector_alloc(3);

    for (i = 0; i < 3; ++i) {
        gsl_vector_set(gsl_a, i, i);
        gsl_vector_set(gsl_b, i, 2 - i);
    }

    gsl_blas_ddot(gsl_a, gsl_b, &dotprod);
    
    CU_ASSERT_DOUBLE_EQUAL(
        dotProduct(&a, &b),
        dotprod,
        0.01);

    vector_free(&a);
    vector_free(&b);
    gsl_vector_free(gsl_a);
    gsl_vector_free(gsl_b);
}

void testDaxpy()
{
    vector a, b;
    gsl_vector *gsl_a, *gsl_b;
    int i;

    vector_calloc(&a, 3);
    vector_calloc(&b, 3);

    for (i = 0; i < 3; ++i) {
        a.data[i] = i;
        b.data[i] = 2 - i;
    }

    gsl_a = gsl_vector_alloc(3);
    gsl_b = gsl_vector_alloc(3);

    for (i = 0; i < 3; ++i) {
        gsl_vector_set(gsl_a, i, i);
        gsl_vector_set(gsl_b, i, 2 - i);
    }

    gsl_blas_daxpy(2.5, gsl_a, gsl_b);
    daxpy(2.5, &a, &b);

    for (i = 0; i < 3; ++i) {
        CU_ASSERT_DOUBLE_EQUAL(
            a.data[i],
            gsl_vector_get(gsl_a, i),
            0.01);
    }

    vector_free(&a);
    vector_free(&b);
    gsl_vector_free(gsl_a);
    gsl_vector_free(gsl_b);
}

void testDnrm2()
{
    vector a;
    gsl_vector *gsl_a;
    int i;

    vector_calloc(&a, 3);
    gsl_a = gsl_vector_alloc(3);

    for (i = 0; i < 3; ++i) {
        a.data[i] = i;
        gsl_vector_set(gsl_a, i, i);
    }

    CU_ASSERT_DOUBLE_EQUAL(
        gsl_blas_dnrm2(gsl_a),
        dnrm2(&a),
        0.01);

    vector_free(&a);
    gsl_vector_free(gsl_a);
}

void testAdd()
{
    vector a, b;
    gsl_vector *gsl_a, *gsl_b;
    int i;

    vector_calloc(&a, 3);
    vector_calloc(&b, 3);
    gsl_a = gsl_vector_alloc(3);
    gsl_b = gsl_vector_alloc(3);

    for (i = 0; i < 3; ++i) {
        a.data[i] = i;
        b.data[i] = 2 - i;
        gsl_vector_set(gsl_a, i, i);
        gsl_vector_set(gsl_b, i, 2 - i);
    }

    gsl_vector_add(gsl_a, gsl_b);
    add(&a, &b);

    for (i = 0; i < 3; ++i) {
        CU_ASSERT_DOUBLE_EQUAL(
            a.data[i],
            gsl_vector_get(gsl_a, i),
            0.01);
    }

    vector_free(&a);
    vector_free(&b);
    gsl_vector_free(gsl_a);
    gsl_vector_free(gsl_b);
}

void testScale()
{
    vector a;
    int i;

    vector_calloc(&a, 3);

    for (i = 0; i < 3; ++i) {
        a.data[i] = i;
    }

    scale(-3.0, &a);

    for (i = 0; i < 3; ++i) {
        CU_ASSERT_DOUBLE_EQUAL(
            a.data[i],
            -3.0 * i,
            0.01);
    }

    vector_free(&a);
}
