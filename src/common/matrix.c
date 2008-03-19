/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file matrix.c
 * Implementation of the method declarations in matrix.h.
 *
 * @author Dominik Dahlem
 */
#include <stdlib.h>

#include "vector.h"
#include "matrix.h"



void cds_matrix_alloc(matrix *mat, int len, int *elem)
{
    int i;

    mat->len = len;
    mat->diags = malloc(sizeof(vector) * len);

    /* allocate memory for the bands */
    for (i = 0; i < len; ++i) {
        vector_calloc(&(mat->diags[i]), elem[i]);
    }
}

void cds_matrix_free(matrix *mat)
{
    int i;

    for (i = 0; i < mat->len; ++i) {
        vector_free(&(mat->diags[i]));
    }

    free(mat->diags);
}

void cdsgb_matrix_alloc(matrix *mat, int len, int *elem, int *index)
{
    int i;

    mat->len = len;
    mat->diags = malloc(sizeof(vector) * len);
    mat->index = malloc(sizeof(int) * len);

    /* allocate memory for the bands */
    for (i = 0; i < len; ++i) {
        vector_alloc(&(mat->diags[i]), elem[i]);
        mat->index[i] = index[i];
    }
}

void cdsgb_matrix_free(matrix *mat)
{
    int i;

    for (i = 0; i < mat->len; ++i) {
        vector_free(&(mat->diags[i]));
    }

    free(mat->diags);
    free(mat->index);
}

void matrix_alloc(matrix *mat, int *elem, int *index)
{
    if (mat->mtype == SB) {
        cds_matrix_alloc(mat, mat->len, elem);
    } else if (mat->mtype == GB) {
        cdsgb_matrix_alloc(mat, mat->len, elem, index);
    }
}

void matrix_free(matrix *mat)
{
    if (mat->mtype == SB) {
        cds_matrix_free(mat);
    } else if (mat->mtype == GB) {
        cdsgb_matrix_free(mat);
    }
}

void matrix_print(matrix *mat)
{
    int i;

    for (i = 0; i < mat->len; ++i) {
        vector_print(&(mat->diags[i]));
    }
}
