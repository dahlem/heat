/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */
#include <string.h>

#include "cds_matrix.h"
#include "conjugate.h"
#include "error.h"
#include "mult.h"
#include "vector.h"



int conjugate(cdsgb_matrix *A, vector *b, vector *x, vector **x_bar)
{
    vector *r, *p, *temp;
    size_t k;
    double dp, dotprod, alpha, beta, norm;
    

    if (b->len != A->diags[0].len) {
        return MATRIX_VECTOR_UNEQUAL_ROW_DIM;
    }

    /* allocate memory for the vectors */
    vector_calloc(*x_bar, (*x_bar)->len);
    vector_calloc(p, (*x_bar)->len);
    vector_calloc(r, (*x_bar)->len);
    vector_calloc(temp, (*x_bar)->len);

    memcpy((*x_bar)->data, x->data, x->len * sizeof(double));
    memcpy(p->data, b->data, b->len * sizeof(double));
    dcdsgbmv(A, x, b);
    memcpy(r->data, p->data, p->len * sizeof(double));
    
    alpha = beta = 0.0;
    k = 1;
    norm = dnrm2(r);
    dp = norm * norm;

    while ((k <= A->len) || (dp > 1e-12)) {
        dcdsgbmv(A, p, temp);

        dotprod = dotProduct(p, temp);
        alpha = dp / dotprod;
        daxpy(alpha, p, *x_bar);
        daxpy(-alpha, temp, r);

        norm = dnrm2(r);
        beta = (norm * norm) / dp;
        memcpy(temp->data, p->data, p->len * sizeof(double));
        daxpy(beta, temp, temp);
        add(temp, r);
        memcpy(p->data, temp->data, temp->len * sizeof(double));
        norm = dnrm2(r);
        dp = norm * norm;
        k++;
    }

    vector_free(temp);
    vector_free(p);
    vector_free(r);
    
    return 0;
}
