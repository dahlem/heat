/* Copyright (C) 2007 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "cds_matrix.h"
#include "mult.h"

#include "mult_test.h"


void registerMultTests()
{
    if (CU_register_suites(mult_suites) != CUE_SUCCESS) {
        fprintf(stderr, "Suite registration failed - %s\n", CU_get_error_msg());
        exit(CU_get_error());
    }
}

void testMultiplication()
{
    int i, j;
    gsl_matrix *A_gsl;
    gsl_vector *u_gsl, *temp;
    cds_matrix A;
    vector u, v;
    int elems[3];

    A_gsl = gsl_matrix_calloc(9, 9);
    u_gsl = gsl_vector_alloc(9);
    temp = gsl_vector_calloc(9);

    /* init the u vector */
    gsl_vector_set_all(u_gsl, 2);

    /* assign values */
    for (i = 0; i < A_gsl->size1; ++i) {
        for (j = 0; j < A_gsl->size2; ++j) {
            /* assign values to the diagonal */
            if (i == j) {
                gsl_matrix_set(A_gsl, i, j, 4);
            }

            /* assign upper/lower diagonal with offset 1 */
            if ((j == (i + 1)) || (i == (j + 1))) {
                gsl_matrix_set(A_gsl, i, j, -1);
            }

            /* assign the upper/lower diagonal with offset 3  */
            if ((j == (i + 3)) || (i == (j + 3))) {
                gsl_matrix_set(A_gsl, i, j, -1);
            }
        }
    }

    gsl_blas_dgemv(CblasNoTrans, 1.0, A_gsl, u_gsl, 0.0, temp);

    /* allocate memory for the bands */
    elems[0] = 9;
    elems[1] = 8;
    elems[2] = 6;
    
    cds_matrix_alloc(&A, 3, elems);
    
    /* assign values */
    for (i = 0; i < A.diags[0].len; ++i) {
        A.diags[0].data[i] = 4;
    }
    for (i = 0; i < A.diags[1].len; ++i) {
        A.diags[1].data[i] = -1;
    }
    for (i = 0; i < A.diags[2].len; ++i) {
        A.diags[2].data[i] = -1;
    }

    vector_alloc(&u, 9);
    vector_alloc(&v, 9);

    for (i = 0; i < 9; ++i) {
        u.data[i] = 2;
    }

    dcdssbmv(&A, &u, &v);

    for (i = 0; i < temp->size; ++i) {
        CU_ASSERT_DOUBLE_EQUAL(
            gsl_vector_get(temp, i),
            v.data[i],
            0.01);
    }
    
    vector_free(&u);
    vector_free(&v);
    cds_matrix_free(&A);
    
    gsl_vector_free(temp);
    gsl_vector_free(u_gsl);
    gsl_matrix_free(A_gsl);
}
