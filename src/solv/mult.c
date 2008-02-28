/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */
#include <gsl/gsl_math.h>

#include "matrix.h"
#include "mult.h"
#include "vector.h"



void dcdssbmv(const matrix *const mat, const vector *const u, const vector *v)
{
    int i, j;
    int p;

    for (i = 0; i < mat->diags[0].len; ++i) {
        /* 1. main diagonal part */
        v->data[i] = mat->diags[0].data[i] * u->data[i];

        for (j = 1; j < mat->len; ++j) {
            /* (0-th row) column and (last column) row index of the j-th diagonal */
            p = mat->diags[0].len - mat->diags[j].len;

            /* 2. check upper diagonal elements */
            if (i < (mat->diags[0].len - p)) {
                v->data[i] += mat->diags[j].data[i] * u->data[i + p];
            }
        
            /* 3. check lower diagonal elements */
            if ((i - p) >= 0) {
                v->data[i] += mat->diags[j].data[i - p] * u->data[i - p];
            }
        }
    }
}

void dcdsgbmv(const matrix *const mat, const vector *const u, const vector *v)
{
    int i, j;
    
    zero(v);

    for (i = 0; i < mat->len; ++i) {
        for (j = GSL_MAX(0, 0 - mat->index[i]);
             j < GSL_MIN(mat->diags[i].len, u->len - mat->index[i]);
             ++j) {
            v->data[j] +=
                mat->diags[i].data[j]
                * u->data[j + mat->index[i]];
        }
    }
}

void dgbmv(const matrix *const mat, const vector *const u, const vector *v)
{
    if (mat->mtype == SB) {
        dcdssbmv(mat, u, v);
    } else if (mat->mtype == GB) {
        dcdsgbmv(mat, u, v);
    }
}
