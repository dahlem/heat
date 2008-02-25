/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */
#include <stdlib.h>

#include "cds_matrix.h"



void vector_alloc(vector *vec, int len)
{
    vec->len = len;
    vec->data = malloc(sizeof(double) * len);
}

void vector_free(vector *vec)
{
    free(vec->data);
}

void cds_matrix_alloc(cds_matrix *mat, int len, int *elem)
{
    int i;
    
    mat->len = len;
    mat->diags = malloc(sizeof(vector) * len);

    /* allocate memory for the bands */
    for (i = 0; i < len; ++i) {
        vector_alloc(&(mat->diags[i]), elem[i]);
    }
}

void cds_matrix_free(cds_matrix *mat)
{
    int i;

    for (i = 0; i < mat->len; ++i) {
        vector_free(&(mat->diags[i]));
    }
    free(mat->diags);
}
