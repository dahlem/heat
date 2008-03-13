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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "matrix.h"
#include "poiss_2d.h"
#include "vector.h"


void init_matrix(matrix *A, int block_matrix_dim)
{
    int i;

#ifdef HAVE_MPI
    int len;
    int row_adjust, system_dim;
    int rows;

    system_dim = block_matrix_dim * block_matrix_dim;
    row_adjust = system_dim % mpiArgs.num_tasks;
    rows = (system_dim + row_adjust) / mpiArgs.num_tasks;

    /* the diagonals have equal lengths */
    len = A->diags[0].len;

    /* the last process potentially needs to adjust for unequal decomposition */
    if (mpiArgs.rank == (mpiArgs.num_tasks - 1)) {
        /* we discount row_adjust from the length of the diagonals */
        len -= row_adjust;
    }

    /* set up partial matrix for each processor */
    for (i = 0; i < len; ++i) {
        if ((i + A->index[0]) < 0 ) {
            A->diags[0].data[i] = 0;
        } else {
            A->diags[0].data[i] = I_DIAG;
        }
    }

    for (i = 0; i < len; ++i) {
        if ((i + A->index[1]) < 0 ) {
            A->diags[1].data[i] = 0;
        } else {
            if (((i + mpiArgs.rank * rows) % block_matrix_dim) == 0) {
                A->diags[1].data[i] = 0;
            } else {
                A->diags[1].data[i] = D_DIAG_1;
            }
        }
    }

    for (i = 0; i < len; ++i) {
        A->diags[2].data[i] = D_MAIN_DIAG;
    }

    for (i = 0; i < len; ++i) {
        if ((i + A->index[3]) < (((mpiArgs.num_tasks - 1) * A->diags[3].len) + len)) {
            if (((i + 1 + mpiArgs.rank * rows) % block_matrix_dim) == 0) {
                A->diags[3].data[i] = 0;
            } else {
                A->diags[3].data[i] = D_DIAG_1;
            }
        } else {
            A->diags[3].data[i] = 0;
        }
    }

    for (i = 0; i < len; ++i) {
        if ((i + A->index[4]) < (((mpiArgs.num_tasks - 1) * A->diags[4].len) + len)) {
            A->diags[4].data[i] = I_DIAG;
        } else {
            A->diags[4].data[i] = 0;
        }
    }
#else
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
}

void init_v(vector *v, int dim, double delta, double coord[2][2],
            double (*src_dens_funcPtr)(double, double),
            double (*bound_cond_funcPtr)(double, double))
{
    double d_squared;           /* delta^2 */
    double rows, cols;          /* max rows and cols; adjusted in the MPI case */
    int i, j;                   /* v-indeces; adjusted in the MPI case */
    int row_offset, col_offset; /* starting point for the v-indeces */
    int sys_dim;                /* system dimension */
    int row_adjust;             /* adjust for unequal row-decompositioning */
    int row_index;              /* row index for the v-vector */
    int row_block;              /* number of rows responsible */

    d_squared = delta * delta;
    sys_dim = dim * dim;
    row_index = 0;
    cols = dim;

#ifdef HAVE_MPI
    row_adjust = sys_dim % mpiArgs.num_tasks;
    row_offset = (double) dim / (double) mpiArgs.num_tasks * (double) mpiArgs.rank;
    rows = MIN(ceil((double) (mpiArgs.rank + 1) * (double) dim
                    / (double) mpiArgs.num_tasks),
               dim);
    row_block = ((sys_dim + row_adjust) / mpiArgs.num_tasks);
    col_offset = (mpiArgs.rank * row_block) % dim;
#else
    row_offset = 0;
    rows = sys_dim;
    row_adjust = 0;
    row_block = sys_dim;
    col_offset = 0;
#endif /* HAVE_MPI */

    /* allocate memory for the v-vector */
    vector_calloc(v, row_block);

    /* assign column-wise from left to right */
    for (i = (int) row_offset; i < rows; ++i) {
        for (j = (int) col_offset; j < cols; ++j) {
            /* we are done assigning values to the v vector */
            if (row_index == row_block) {
                return;
            }

            fprintf(stdout, "f_{%d,%d} = %f\n", i+1, j+1, src_dens_funcPtr(coord[0][0] + (delta * (i + 1)),
                                                                       coord[0][1] + (delta * (j + 1))));

            v->data[row_index] =
                d_squared
                * src_dens_funcPtr(coord[0][0] + (delta * (i + 1)),
                                   coord[0][1] + (delta * (j + 1)));

            if (i == 0) {
                /* left column */

                if (j == 0) {
                    /* bottom row */
                    fprintf(stdout, "v_{0,1} = %f\n", bound_cond_funcPtr(coord[0][0], coord[1][0] + delta));
                    fprintf(stdout, "v_{1,0} = %f\n", bound_cond_funcPtr(coord[0][0] + delta, coord[1][0]));
                    v->data[row_index] +=
                        bound_cond_funcPtr(coord[0][0], coord[1][0] + delta)
                        + bound_cond_funcPtr(coord[0][0] + delta, coord[1][0]);
                } else if (j == (dim - 1)) {
                    /* top row */
                    fprintf(stdout, "v_{0,%d} = %f\n", j+1, bound_cond_funcPtr(coord[0][0], coord[1][0] + (delta * (j + 1))));
                    fprintf(stdout, "v_{1,%d} = %f\n", j+2, bound_cond_funcPtr(coord[0][0] + delta, coord[1][0] + (delta * (j + 2))));
                    v->data[row_index] +=
                        bound_cond_funcPtr(coord[0][0], coord[1][0] + (delta * (j + 1)))
                        + bound_cond_funcPtr(coord[0][0] + delta, coord[1][0] + (delta * (j + 2)));
                } else {
                    /* middle rows */
                    fprintf(stdout, "v_{0,%d} = %f\n", j+1, bound_cond_funcPtr(coord[0][0], coord[1][0] + (delta * (j + 1))));
                    v->data[row_index] +=
                        bound_cond_funcPtr(coord[0][0], coord[1][0] + (delta * (j + 1)));
                }
            } else if (i == (dim - 1)) {
                /* right column */

                if (j == 0) {
                    /* bottom row */
                    fprintf(stdout, "v_{%d,0} = %f\n", i+1, bound_cond_funcPtr(coord[0][0] + (delta * (i + 1)), coord[1][0]));
                    fprintf(stdout, "v_{%d,1} = %f\n", i+2, bound_cond_funcPtr(coord[0][0] + (delta * (i + 2)), coord[1][0] + delta));
                    v->data[row_index] +=
                        bound_cond_funcPtr(coord[0][0] + (delta * (i + 1)), coord[1][0])
                        + bound_cond_funcPtr(coord[0][0] + (delta * (i + 2)), coord[1][0] + delta);
                } else if (j == (dim - 1)) {
                    /* top row */
                    fprintf(stdout, "v_{%d,%d} = %f\n", i+2, j+1, bound_cond_funcPtr(coord[0][0] + (delta * (i + 2)), coord[1][0] + (delta * (j + 1))));
                    fprintf(stdout, "v_{%d,%d} = %f\n", i+1, j+2, bound_cond_funcPtr(coord[0][0] + (delta * (i + 1)), coord[1][0] + (delta * (j + 2))));
                    v->data[row_index] +=
                        bound_cond_funcPtr(coord[0][0] + (delta * (i + 2)), coord[1][0] + (delta * (j + 1)))
                        + bound_cond_funcPtr(coord[0][0] + (delta * (i + 1)), coord[1][0] + (delta * (j + 2)));
                } else {
                    /* middle rows */
                    fprintf(stdout, "v_{%d,%d} = %f\n", i+2, j+1, bound_cond_funcPtr(coord[0][0] + (delta * (i + 2)), coord[1][0] + (delta * (j + 1))));
                    v->data[row_index] +=
                        bound_cond_funcPtr(coord[0][0] + (delta * (i + 2)), coord[1][0] + (delta * (j + 1)));
                }
            } else {
                /* middle column */

                if (j == 0) {
                    /* bottom row */
                    fprintf(stdout, "v_{%d,0} = %f\n", i+1, bound_cond_funcPtr(coord[0][0] + (delta * (i + 1)), coord[1][0]));
                    v->data[row_index] +=
                        bound_cond_funcPtr(coord[0][0] + (delta * (i + 1)), coord[1][0]);
                } else if (j == (dim - 1)) {
                    /* top row */
                    fprintf(stdout, "v_{%d,%d} = %f\n", i+1, j+2, bound_cond_funcPtr(coord[0][0] + (delta * (i + 1)), coord[1][0]));
                    v->data[row_index] +=
                        bound_cond_funcPtr(coord[0][0] + (delta * (i + 1)), coord[1][0] + (delta * (j + 2)));
                }
            }

            row_index++;
        }

        /* reset the column offset, because the next row starts with 0 (as usual) */
        col_offset = 0;
    }
}

void setup_poiss_2d(matrix *A, vector *u, vector *v, vector *x_bar, int dim,
                    double delta, double coord[2][2],
                    double (*src_dens_funcPtr)(double, double),
                    double (*bound_cond_funcPtr)(double, double))
{
#ifdef HAVE_MPI
    int row_adjust;         /* adjust for unequal row-decompositioning */
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

    row_adjust = system_dim % mpiArgs.num_tasks;
    num_diags = POISS_2D_GB_MATRIX_DIAGS;
    rows = (system_dim + row_adjust) / mpiArgs.num_tasks;
    col_offset = mpiArgs.rank * rows;

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
#else
    matrix_alloc(A, elems, NULL);
#endif /* HAVE_MPI */

    /* initialise the matrix */
    init_matrix(A, block_matrix_dim);

    /* initialise the x_bar vector */
    vector_alloc(x_bar, rows);

    /* initialise the u vector */
    vector_alloc(u, rows);
    for (i = 0; i < rows; ++i) {
        u->data[i] = 1;
    }

    /* initialise the v vector */
    init_v(v, block_matrix_dim, delta, coord, src_dens_funcPtr, bound_cond_funcPtr);
    vector_print(v);
}
