/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */
#if HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef HAVE_MPI
# include <mpi.h>
# include "mpi-common.h"
#endif /* HAVE_MPI */

#include <stdlib.h>

#include "matrix.h"
#include "poiss_2d.h"
#include "vector.h"


void setup_poiss_2d(matrix *A, vector *u, vector *v, vector *x_bar, int dim)
{
#ifdef HAVE_MPI
    int col_offset;
    int elems[POISS_2D_GB_MATRIX_DIAGS];
    int index[POISS_2D_GB_MATRIX_DIAGS];
#else
    int elems[POISS_2D_SB_MATRIX_DIAGS];
#endif /* HAVE_MPI */
    int num_diags;
    int i;
    int block_matrix_dim;
    int system_dim;
    int rows;

    block_matrix_dim = dim - 2;
    system_dim = block_matrix_dim * block_matrix_dim;

#ifdef HAVE_MPI
    /* configure the diagonal offsets */
    index[0] = -block_matrix_dim;
    index[1] = -1;
    index[2] = 0;
    index[3] = 1;
    index[4] = block_matrix_dim;

    num_diags = POISS_2D_GB_MATRIX_DIAGS;
    rows = system_dim/mpiArgs.num_tasks;
    col_offset = mpiArgs.rank * mpiArgs.num_tasks;

    for (i = 0; i < num_diags; ++i) {
        elems[i] = rows;
        index[i] += col_offset;
    }
#else
    num_diags = POISS_2D_SB_MATRIX_DIAGS;
    rows = system_dim;

    elems[0] = system_dim;
    elems[1] = system_dim - 1;
    elems[2] = system_dim - block_matrix_dim;
#endif /* HAVE_MPI */

    /* initialise the matrix type */
#ifdef HAVE_MPI
    A->mtype = GB;
#else
    A->mtype = SB;
#endif /* HAVE_MPI */

    A->storage = CDS;
    A->len = num_diags;

#ifdef HAVE_MPI
    matrix_alloc(A, elems, index);

    /* set up partial matrix for each processor */
    for (i = 0; i < A->diags[0].len; ++i) {
        if ((i + index[0]) < 0 ) {
            A->diags[0].data[i] = 0;
        } else {
            A->diags[0].data[i] = I_DIAG;
        }
    }
    
    for (i = 0; i < A->diags[1].len; ++i) {
        if ((i + index[1]) < 0 ) {
            A->diags[1].data[i] = 0;
        } else {
            if ((i % block_matrix_dim) == 0) {
                A->diags[1].data[i] = 0;
            } else {
                A->diags[1].data[i] = D_DIAG_1;
            }
        }
    }
    
    for (i = 0; i < A->diags[2].len; ++i) {
        A->diags[2].data[i] = D_MAIN_DIAG;
    }
    
    for (i = 0; i < A->diags[3].len; ++i) {
        if ((i + index[3]) < (mpiArgs.num_tasks * A->diags[3].len)) {
            if (((i + 1 + mpiArgs.rank * mpiArgs.num_tasks) % block_matrix_dim) == 0) {
                A->diags[3].data[i] = 0;
            } else {
                A->diags[3].data[i] = D_DIAG_1;
            }
        } else {
            A->diags[3].data[i] = 0;
        }
    }
    
    for (i = 0; i < A->diags[4].len; ++i) {
        if ((i + index[4]) < (mpiArgs.num_tasks * A->diags[4].len)) {
            A->diags[4].data[i] = I_DIAG;
        } else {
            A->diags[4].data[i] = 0;
        }
    }
#else
    matrix_alloc(A, elems, NULL);

    /* set up matrix for single processor */
    for (i = 0; i < A->diags[2].len; ++i) {
        A->diags[2].data[i] = I_DIAG;
    }
    
    for (i = 0; i < A->diags[1].len; ++i) {
        if (((i + 1) % block_matrix_dim) == 0) {
            A->diags[1].data[i] = 0;
        } else {
            A->diags[1].data[i] = D_DIAG_1;
        }
    }
    
    for (i = 0; i < A->diags[0].len; ++i) {
        A->diags[0].data[i] = D_MAIN_DIAG;
    }
#endif /* HAVE_MPI */

    /* set up the vectors for each processor */
    vector_alloc(u, rows);
    vector_alloc(v, rows);
    vector_alloc(x_bar, rows);

    for (i = 0; i < rows; ++i) {
        u->data[i] = 1;
        v->data[i] = i;
    }
}
