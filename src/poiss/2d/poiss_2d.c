/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file poiss_2d.c
 * Implementation of the method declarations in poiss_2d.h. In order to initialise
 * the matrix and vectors in an MPI environment the global row/column offsets have
 * to be taken into account, because each processor is responsible to set up its
 * partial matrix and vector structures. There is no initial communication between
 * the processes to set up any vector to avoid a huge communication overhead or
 * requiring a single processor to hold the global structures as this might stretch
 * the local resources available.
 *
 * @see poiss_2d.h
 * @author Dominik Dahlem
 */
#if HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef HAVE_MPI
# include <mpi.h>
# include "mpi-common.h"
# include "mpi-utils.h"
#endif /* HAVE_MPI */

#ifdef HAVE_OPENMP
# include <omp.h>
#endif /* HAVE_OPENMP */

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
#ifdef HAVE_OPENMP
# ifdef NDEBUG
    int tid, nthreads;
# endif /* NDEBUG */
#endif /* HAVE_OPENMP */

#ifdef HAVE_MPI
    int len;
    int row_adjust, system_dim;
    int rows;

    system_dim = block_matrix_dim * block_matrix_dim;
    row_adjust = adjustment(system_dim, mpiArgs.num_tasks);
    rows = block(system_dim, mpiArgs.num_tasks);

    /* the diagonals have equal lengths */
    len = A->diags[0].len;

    /* the last process potentially needs to adjust for unequal decomposition */
    if (mpiArgs.rank == (mpiArgs.num_tasks - 1)) {
        /* we discount row_adjust from the length of the diagonals */
        len -= row_adjust;
    }

    /* sections multi-threading:
       Each for-loop is executed in a separate thread. */
# ifdef HAVE_OPENMP
#  ifdef NDEBUG
#   pragma omp parallel shared(A, nthreads) private(i, tid)
#  else
#   pragma omp parallel shared(A) private(i)
#  endif /* NDEBUG */
    {
#  ifdef NDEBUG
        nthreads = omp_get_num_threads();
        tid = omp_get_thread_num();
        if (tid == 0)
        {
            nthreads = omp_get_num_threads();
            printf("Initializing matrices with %d threads...\n", nthreads);
        }
#  endif /* NDEBUG */

#  pragma omp sections nowait
        {

        /* set up partial matrix for each processor */
#  pragma omp section
            {
#  ifdef NDEBUG
                printf("Thread %d executes the first section\n", tid);
#  endif /* NDEBUG */
# endif /* HAVE_OPENMP */
                for (i = 0; i < len; ++i) {
                    if ((i + A->index[0]) < 0 ) {
                        A->diags[0].data[i] = 0;
                    } else {
                        A->diags[0].data[i] = I_DIAG;
                    }
                }

# ifdef HAVE_OPENMP
            }
#  pragma omp section
            {
#  ifdef NDEBUG
                printf("Thread %d executes the second section\n", tid);
#  endif /* NDEBUG */
# endif /* HAVE_OPENMP */
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

# ifdef HAVE_OPENMP
            }
#  pragma omp section
            {
#  ifdef NDEBUG
                printf("Thread %d executes the third section\n", tid);
#  endif /* NDEBUG */
# endif /* HAVE_OPENMP */
                for (i = 0; i < len; ++i) {
                    A->diags[2].data[i] = D_MAIN_DIAG;
                }

# ifdef HAVE_OPENMP
            }
#  pragma omp section
            {
#  ifdef NDEBUG
                printf("Thread %d executes the fourth section\n", tid);
#  endif /* NDEBUG */
# endif /* HAVE_OPENMP */
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

# ifdef HAVE_OPENMP
            }
#  pragma omp section
            {
#  ifdef NDEBUG
                printf("Thread %d executes the fifth section\n", tid);
#  endif /* NDEBUG */
# endif /* HAVE_OPENMP */
                for (i = 0; i < len; ++i) {
                    if ((i + A->index[4]) < (((mpiArgs.num_tasks - 1) * A->diags[4].len) + len)) {
                        A->diags[4].data[i] = I_DIAG;
                    } else {
                        A->diags[4].data[i] = 0;
                    }
                }
# ifdef HAVE_OPENMP
            } /* end last section */
        } /* end of sections */
    } /* end of parallel section */
# endif /* HAVE_OPENMP */
#else

    /* sections multi-threading:
       Each for-loop is executed in a separate thread. */
# ifdef HAVE_OPENMP
#  ifdef NDEBUG
#   pragma omp parallel shared(A,nthreads) private(i, tid)
#  else
#   pragma omp parallel shared(A) private(i)
#  endif /* NDEBUG */
    {

#  ifdef NDEBUG
        nthreads = omp_get_num_threads();
        tid = omp_get_thread_num();
        if (tid == 0)
        {
            nthreads = omp_get_num_threads();
            printf("Initializing matrices with %d threads...\n", nthreads);
        }
#  endif /* NDEBUG */

#  pragma omp sections nowait
        {
            /* set up matrix for single processor:
               only main and upper-main diagonals needed */
#  pragma omp section
            {
#  ifdef NDEBUG
                printf("Thread %d executes the first section\n", tid);
#  endif /* NDEBUG */
# endif /* HAVE_OPENMP */
                for (i = 0; i < A->diags[2].len; ++i) {
                    A->diags[2].data[i] = I_DIAG;
                }
# ifdef HAVE_OPENMP
            }

#  pragma omp section
            {
#  ifdef NDEBUG
                printf("Thread %d executes the second section\n", tid);
#  endif /* NDEBUG */
# endif /* HAVE_OPENMP */
                for (i = 0; i < A->diags[1].len; ++i) {
                    if (((i + 1) % block_matrix_dim) == 0) {
                        A->diags[1].data[i] = 0;
                    } else {
                        A->diags[1].data[i] = D_DIAG_1;
                    }
                }
# ifdef HAVE_OPENMP
            }
#  pragma omp section
            {
#  ifdef NDEBUG
                printf("Thread %d executes the third section\n", tid);
#  endif /* NDEBUG */
# endif /* HAVE_OPENMP */
                for (i = 0; i < A->diags[0].len; ++i) {
                    A->diags[0].data[i] = D_MAIN_DIAG;
                }
# ifdef HAVE_OPENMP
            } /* end last section */
        } /* end of sections */
    } /* end of parallel section */
# endif /* HAVE_OPENMP */
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
    int row_index;              /* row index for the v-vector */
    int row_block;              /* number of rows responsible */

#ifdef HAVE_MPI
    int row_adjust;             /* adjust for unequal row-decompositioning */
#endif

    d_squared = delta * delta;
    sys_dim = dim * dim;
    row_index = 0;
    cols = dim;

#ifdef HAVE_MPI
    row_adjust = adjustment(sys_dim, mpiArgs.num_tasks);
    row_offset = (double) dim / (double) mpiArgs.num_tasks * (double) mpiArgs.rank;
    rows = MIN(ceil((double) (mpiArgs.rank + 1) * (double) dim
                    / (double) mpiArgs.num_tasks),
               dim);
    row_block = block(sys_dim, mpiArgs.num_tasks);
    col_offset = (mpiArgs.rank * row_block) % dim;
#else
    row_offset = 0;
    rows = sys_dim;
    row_block = sys_dim;
    col_offset = 0;
#endif /* HAVE_MPI */

    /* allocate memory for the v-vector */
    vector_calloc(v, row_block);

    /* assign column-wise from left to right */
    for (i = row_offset; i < rows; ++i) {
        for (j = col_offset; j < cols; ++j) {
            /* we are done assigning values to the v vector */
            if (row_index == row_block) {
                return;
            }

            /* v_{i,j} += \Delta^2 * f_{i,j}  */
            v->data[row_index] = d_squared;
            v->data[row_index] *= src_dens_funcPtr(
                coord[0][0] + (delta * ((double) i + 1.0)),
                coord[1][0] + (delta * ((double) j + 1.0)));

            if (i == 0) {
                /* left column */

                if (j == 0) {
                    /* bottom row */
                    /* v_{1,1} += v_{0,1} + v_{1,0} */
                    v->data[row_index] +=
                        bound_cond_funcPtr(
                            coord[0][0],
                            coord[1][0] + delta)
                        + bound_cond_funcPtr(
                            coord[0][0] + delta,
                            coord[1][0]);
                } else if (j == (dim - 1)) {
                    /* top row */
                    /* v_{1,j} += v_{0,j} + v_{1,j+1} */
                    v->data[row_index] +=
                        bound_cond_funcPtr(
                            coord[0][0],
                            coord[1][0] + (delta * ((double) j + 1.0)))
                        + bound_cond_funcPtr(
                            coord[0][0] + delta,
                            coord[1][0] + (delta * ((double) j + 2.0)));
                } else {
                    /* middle rows */
                    /* v_{1,N} += v_{0,N} */
                    v->data[row_index] +=
                        bound_cond_funcPtr(
                            coord[0][0],
                            coord[1][0] + (delta * ((double) j + 1.0)));
                }
            } else if (i == (dim - 1)) {
                /* right column */

                if (j == 0) {
                    /* bottom row */
                    /* v_{N,1} += v_{N+1,1} + v_{N,0} */
                    v->data[row_index] +=
                        bound_cond_funcPtr(
                            coord[0][0] + (delta * ((double) i + 1.0)),
                            coord[1][0])
                        + bound_cond_funcPtr(
                            coord[0][0] + (delta * ((double) i + 2.0)),
                            coord[1][0] + delta);
                } else if (j == (dim - 1)) {
                    /* top row */
                    /* v_{N,N} += v_{N+1,1} + v_{N,N+1} */
                    v->data[row_index] +=
                        bound_cond_funcPtr(
                            coord[0][0] + (delta * ((double) i + 2.0)),
                            coord[1][0] + (delta * ((double) j + 1.0)))
                        + bound_cond_funcPtr(
                            coord[0][0] + (delta * ((double) i + 1.0)),
                            coord[1][0] + (delta * ((double) j + 2.0)));
                } else {
                    /* middle rows */
                    /* v_{i,j} += v_{i+1,j} + v_{i,j+1} */
                    v->data[row_index] +=
                        bound_cond_funcPtr(
                            coord[0][0] + (delta * ((double) i + 2.0)),
                            coord[1][0] + (delta * ((double) j + 1.0)));
                }
            } else {
                /* middle column */

                if (j == 0) {
                    /* bottom row */
                    /* v_{i,1} += v_{i,0} */
                    v->data[row_index] +=
                        bound_cond_funcPtr(
                            coord[0][0] + (delta * ((double) i + 1.0)),
                            coord[1][0]);
                } else if (j == (dim - 1)) {
                    /* top row */
                    /* v_{i,N} += v_{i,N} */
                    v->data[row_index] +=
                        bound_cond_funcPtr(
                            coord[0][0] + (delta * ((double) i + 1.0)),
                            coord[1][0] + (delta * ((double) j + 2.0)));
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
#ifdef HAVE_OPENMP
# ifdef NDEBUG
    int tid, nthreads;
# endif /* NDEBUG */
#endif /* HAVE_OPENMP */

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

    row_adjust = adjustment(system_dim, mpiArgs.num_tasks);
    rows = block(system_dim, mpiArgs.num_tasks);
    num_diags = POISS_2D_GB_MATRIX_DIAGS;
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

    /* sections multi-threading:
       Each for-loop is executed in a separate thread. */
#ifdef HAVE_OPENMP
# ifdef NDEBUG
#  pragma omp parallel shared(A,nthreads) private(i, tid)
# else
#  pragma omp parallel shared(A) private(i)
# endif /* NDEBUG */
    {

# ifdef NDEBUG
        nthreads = omp_get_num_threads();
        tid = omp_get_thread_num();
        if (tid == 0)
        {
            nthreads = omp_get_num_threads();
            printf("Initializing the linear system with %d threads...\n", nthreads);
        }
# endif /* NDEBUG */

# pragma omp sections nowait
        {
# pragma omp section
            {
# ifdef NDEBUG
                printf("Thread %d executes the first section (init matrix)\n", tid);
# endif /* NDEBUG */
#endif /* HAVE_OPENMP */
                /* initialise the matrix */
                init_matrix(A, block_matrix_dim);
#ifdef HAVE_OPENMP
            }
# pragma omp section
            {
# ifdef NDEBUG
                printf("Thread %d executes the second section (init x_bar)\n", tid);
# endif /* NDEBUG */
#endif /* HAVE_OPENMP */
                /* initialise the x_bar vector */
                vector_alloc(x_bar, rows);

                /* initialise the u vector */
                vector_alloc(u, rows);
                for (i = 0; i < rows; ++i) {
                    u->data[i] = 1;
                }
#ifdef HAVE_OPENMP
            }
# pragma omp section
            {
# ifdef NDEBUG
                printf("Thread %d executes the third section (init v)\n", tid);
# endif /* NDEBUG */
#endif /* HAVE_OPENMP */

                /* initialise the v vector */
                init_v(v, block_matrix_dim, delta, coord, src_dens_funcPtr, bound_cond_funcPtr);
#ifdef HAVE_OPENMP
            } /* end last section */
        } /* end of sections */
    } /* end of parallel section */
#endif /* HAVE_OPENMP */
}
