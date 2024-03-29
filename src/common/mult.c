/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@gmail.com>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file mult.c
 * Implementation of the method declarations in mult.h. The dgbmv method supports
 * the parallel calculation of the matrix-vector product. This method expects the
 * matrix to be in the general band format, which can be set in the matrix structure
 * declared in matrix.h.
 *
 * @author Dominik Dahlem
 */
#if HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef HAVE_MPI
# include <stdio.h>
# include <mpi.h>
# include "mpi-common.h"
#endif /* HAVE_MPI */

#ifdef HAVE_OPENMP
# include <omp.h>
#endif /* HAVE_OPENMP */

#include "macros.h"
#include "matrix.h"
#include "mult.h"
#include "vector.h"


void dcdssbmv(const matrix *const mat, const vector *const u, const vector *v)
{
    int i, j;
    int p;

#ifdef HAVE_OPENMP
# pragma omp for private(i, j, p)
#endif /* HAVE_OPENMP */
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
#ifdef HAVE_OPENMP
# pragma omp for private(j)
#endif /* HAVE_OPENMP */
        for (j = MAX(0, 0 - mat->index[i]);
             j < MIN(mat->diags[i].len, u->len - mat->index[i]);
             ++j) {
            v->data[j] +=
                mat->diags[i].data[j]
                * u->data[j + mat->index[i]];
        }
    }
}

void dgbmv(const matrix *const mat, const vector *const u, const vector *v)
{
#ifdef HAVE_MPI
    int mpi_ret;
    vector u_global;
#endif /* HAVE_MPI */

    /* allocate memory for the vectors */
#ifdef HAVE_MPI
    vector_calloc(&u_global, mpiArgs.num_tasks * u->len);
#endif /* HAVE_MPI */

    /* obtain full u-vector, if in MPI environment */
#ifdef HAVE_MPI
    mpi_ret = MPI_Allgather(u->data, u->len, MPI_DOUBLE, u_global.data,
                            u->len, MPI_DOUBLE, MPI_COMM_WORLD);

#ifdef NDEBUG
    fprintf(stdout, "Process %d got u-vector:\n", mpiArgs.rank);
    vector_print(&u_global);
#endif /* NDEBUG */

    if (mpi_ret != MPI_SUCCESS) {
        fprintf(stderr, "The MPI_Allgather primitive returned error: %d.\n", mpi_ret);
        fflush(stderr);
    }
#endif /* HAVE_MPI*/

    if (mat->mtype == SB) {
        dcdssbmv(mat, u, v);
    } else if (mat->mtype == GB) {
#ifdef HAVE_MPI
        dcdsgbmv(mat, &u_global, v);
#else
        dcdsgbmv(mat, u, v);
#endif /* HAVE_MPI*/
    }

    /* clean-up allocated memory */
#ifdef HAVE_MPI
    vector_free(&u_global);
#endif /* HAVE_MPI */
}
