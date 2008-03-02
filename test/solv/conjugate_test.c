/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
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

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include <CUnit/CUnit.h>

#include "conjugate.h"
#include "vector.h"
#include "matrix.h"

#include "conjugate_test.h"


void registerConjugateTests()
{
    if (CU_register_suites(conjugate_suites) != CUE_SUCCESS) {
        fprintf(stderr, "Suite registration failed - %s\n", CU_get_error_msg());
        exit(CU_get_error());
    }
}

void testConjugateSB()
{
    matrix A;
    vector u, v, x_bar;
    size_t i;
    int elems[5];
    gsl_vector *x_bar_gsl, *temp;
    
    double a_array[] = {
        1.86279, 0.47863, -0.54877,
        0.47863, 1.61609, 0.10628,
        -0.54877, 0.10628, 2.76115
    };

    gsl_matrix_view A_gsl = gsl_matrix_view_array(a_array, 3, 3);

    x_bar_gsl = gsl_vector_alloc(3);
    temp = gsl_vector_alloc(3);

    vector_calloc(&x_bar, 3);
    vector_calloc(&u, 3);
    vector_calloc(&v, 3);

    elems[0] = 3;
    elems[1] = 2;
    elems[2] = 1;

    A.mtype = SB;
    A.storage = CDS;
    A.len = 3;
    matrix_alloc(&A, elems, NULL);
    
    A.diags[0].data[0] = 1.86279;
    A.diags[0].data[1] = 1.61609;
    A.diags[0].data[2] = 2.76115;

    A.diags[1].data[0] = 0.47863;
    A.diags[1].data[1] = 0.10628;

    A.diags[2].data[0] = -0.54877;

    for (i = 0; i < 3; ++i) {
        v.data[i] = i + 1;
        u.data[i] = 1;
    }

    CU_ASSERT_EQUAL(conjugate(&A, &v, &u, &x_bar), 0);

    temp = gsl_vector_calloc(v.len);

    for (i = 0; i < x_bar.len; ++i) {
        gsl_vector_set(x_bar_gsl, i, x_bar.data[i]);
    }
    
    gsl_blas_dgemv(CblasNoTrans, 1.0, &A_gsl.matrix, x_bar_gsl, 0.0, temp);

    for (i = 0; i < temp->size; ++i) {
        CU_ASSERT_DOUBLE_EQUAL(
            v.data[i],
            gsl_vector_get(temp, i),
            0.01);
    }

    vector_free(&x_bar);
    vector_free(&u);
    vector_free(&v);
    matrix_free(&A);
    gsl_vector_free(temp);
}
