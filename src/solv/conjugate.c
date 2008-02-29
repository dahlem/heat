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
#include <string.h>

#include "matrix.h"
#include "conjugate.h"
#include "error.h"
#include "mult.h"
#include "vector.h"



int conjugate(matrix *A, vector *b, vector *x, vector *x_bar)
{
    vector r, p, temp;
    size_t k;
    double error, dotprod, alpha, beta, norm;
    
    if (b->len != A->diags[0].len) {
        return MATRIX_VECTOR_UNEQUAL_ROW_DIM;
    }

    /* allocate memory for the vectors */
    vector_calloc(&p, b->len);
    vector_calloc(&r, b->len);
    vector_calloc(&temp, b->len);

    /* calculate the residual for an initially chosen vector
     * \f$ r_0 = b - A * x_0 \f$.
     */
    memcpy(r.data, b->data, b->len * sizeof(double));
    memcpy(x_bar->data, x->data, x->len * sizeof(double));
    dgbmv(A, x_bar, &temp);
    daxpy(-1.0, &temp, &r);

    /* set \f$ p_0 = r_0 \f$*/
    memcpy(p.data, r.data, r.len * sizeof(double));

    /* initial values \f$ \alpha = \beta = 0 \f$ */
    alpha = beta = 0.0;

    /* initial value for the while counter \f$ k = 1 \f$ */
    k = 1;

    /* calculate the error */
    /* \f$ dp = (r, r) \f$ */
    error = dotProduct(&r, &r);

    while ((k <= A->len) || (error > 1e-12)) {
        /* \f$ v = A * p \f$ */
        dgbmv(A, &p, &temp);

        /* \f$ \alpha = (r' * r)/(p' * v)  \f$ */
        dotprod = dotProduct(&p, &temp);
        alpha = error / dotprod;
        
        /* \f$ x = x + \alpha * p   \f$ */
        daxpy(alpha, &p, x_bar);

        /* \f$ new_r = r - \alpha * v  \f$ */
        daxpy(-alpha, &temp, &r);

        /* \f$ \beta = (new_r' * new_r)/(r' * r) \f$ */
        norm = dnrm2(&r);
        beta = (norm * norm) / error;
        
        /* \f$ p = new_r + \beta * p  \f$ */
        memcpy(temp.data, p.data, p.len * sizeof(double));
        scale(beta, &temp);
        add(&temp, &r);
        memcpy(p.data, temp.data, temp.len * sizeof(double));

        /* \f$ dp = (r, r) \f$ */
        error = dotProduct(&r, &r);
        k++;
    }

    /* clean-up allocated memory */
    vector_free(&temp);
    vector_free(&p);
    vector_free(&r);
    
    return EXIT_SUCCESS;
}