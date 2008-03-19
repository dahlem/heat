/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file gnuplot.c
 * This file implements the declarations in gnuplot.h. The print_surface method
 * supports an MPI environment and serial environment, where the MPI primitive used
 * is MPI_Gather to collect all the vector blocks. Possibly a better way of doing it
 * would be to use MPI-IO instead.
 *
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cl.h"
#include "error.h"
#include "gnuplot.h"


int print_surface(vector *vec, int dim, double (*bound_cond_funcPtr)(double, double))
{
    int i, j;
    FILE *input;
    double x, y, z, error;
    double *global_buf;

#ifdef HAVE_MPI
    double sys_dim;
    int buf_len;
    int status;

    sys_dim = dim * dim;
    buf_len = block(sys_dim, mpiArgs.num_tasks);

    if ((global_buf = malloc(sys_dim * sizeof(double))) == NULL) {
        return MALLOC_ERROR;
    }

    /* the root process gathers the solution vector */
    if ((status = MPI_Gather(vec->data, buf_len, MPI_DOUBLE,
                             global_buf, buf_len, MPI_DOUBLE,
                             0, MPI_COMM_WORLD)) != MPI_SUCCESS) {
        fprintf(stderr, "The MPI_Gather primitive returned error %d.\n", status);
        fflush(stderr);
        return EXIT_FAILURE;
    }

    /* only the root process writes into the file */
    if (mpiArgs.rank == 0) {

#else
        global_buf = vec->data;
#endif /* HAVE_MPI */

        if ((input = fopen(globalArgs.f, "w")) == NULL) {
            return FILE_OPEN_FOR_WRITE_ERROR;
        }

        /* boundary conditions */
        fprintf(input, "%.6f,%.6f,%.10f,0.0\n", globalArgs.x0, globalArgs.y0,
                bound_cond_funcPtr(globalArgs.x0, globalArgs.y0));
        for (j = 0; j < dim; ++j) {
            y = globalArgs.y0 + (j + 1) * globalArgs.d;
            z = bound_cond_funcPtr(globalArgs.x0, y);
            fprintf(input, "%.6f,%.6f,%.10f,0.0\n", globalArgs.x0, y, z);
        }
        fprintf(input, "%.6f,%.6f,%.10f,0.0\n\n", globalArgs.x0, globalArgs.y1,
                bound_cond_funcPtr(globalArgs.x0, globalArgs.y1));

        for (i = 0; i < dim; ++i) {
            x = globalArgs.x0 + (i + 1) * globalArgs.d;
            z = bound_cond_funcPtr(x, globalArgs.y0);

            /* boundary conditions */
            fprintf(input, "%.6f,%.6f,%.10f,0.0\n", x, globalArgs.y0, z);

            for (j = 0; j < dim; ++j) {
                y = globalArgs.y0 + (j + 1) * globalArgs.d;
                z = global_buf[i * dim + j];
                fprintf(input, "%.6f,%.6f,%.10f,%.10f\n", x, y, z,
                        fabs(z - bound_cond_funcPtr(x, y)));
            }

            /* boundary conditions */
            fprintf(input, "%.6f,%.6f,%.10f,0.0\n\n", x, globalArgs.y1,
                    bound_cond_funcPtr(x, globalArgs.y1));
        }

        /* boundary conditions */
        fprintf(input, "%.6f,%.6f,%.10f,0.0\n", globalArgs.x1, globalArgs.y0,
                bound_cond_funcPtr(globalArgs.x1, globalArgs.y0));
        for (j = 0; j < dim; ++j) {
            y = globalArgs.y0 + (j + 1) * globalArgs.d;
            z = bound_cond_funcPtr(globalArgs.x1, y);
            fprintf(input, "%.6f,%.6f,%.10f,0.0\n", globalArgs.x1, y, z);
        }
        fprintf(input, "%.6f,%.6f,%.10f,0.0\n\n", globalArgs.x1, globalArgs.y1,
                bound_cond_funcPtr(globalArgs.x1, globalArgs.y1));

        /* close the file */
        fclose(input);

#ifdef HAVE_MPI
    }
#endif /* HAVE_MPI */

    return EXIT_SUCCESS;
}
