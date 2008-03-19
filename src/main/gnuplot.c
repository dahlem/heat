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
#include <stdio.h>
#include <stdlib.h>

#include "cl.h"
#include "error.h"
#include "gnuplot.h"


int print_surface(vector *vec, int dim, double (*bound_cond_funcPtr)(double, double))
{
    int i, j;
    FILE *input;
    double x, y, z;

    if ((input = fopen(globalArgs.f, "w")) == NULL) {
        return FILE_OPEN_FOR_WRITE_ERROR;
    }

    /* boundary conditions */
    fprintf(input, "%f,%f,%f\n", globalArgs.x0, globalArgs.y0,
            bound_cond_funcPtr(globalArgs.x0, globalArgs.y0));
    for (j = 0; j < dim; ++j) {
        y = globalArgs.y0 + (j + 1) * globalArgs.d;
        z = bound_cond_funcPtr(globalArgs.x0, y);
        fprintf(input, "%f,%f,%f\n", globalArgs.x0, y, z);
    }
    fprintf(input, "%f,%f,%f\n\n", globalArgs.x0, globalArgs.y1,
            bound_cond_funcPtr(globalArgs.x0, globalArgs.y1));

    for (i = 0; i < dim; ++i) {
        x = globalArgs.x0 + (i + 1) * globalArgs.d;

        /* boundary conditions */
        fprintf(input, "%f,%f,%f\n", x, globalArgs.y0,
                bound_cond_funcPtr(x, globalArgs.y0));

        for (j = 0; j < dim; ++j) {
            y = globalArgs.y0 + (j + 1) * globalArgs.d;
            z = vec->data[i * dim + j];
            fprintf(input, "%f,%f,%f\n", x, y, z);
        }

        /* boundary conditions */
        fprintf(input, "%f,%f,%f\n\n", x, globalArgs.y1,
                bound_cond_funcPtr(x, globalArgs.y1));
    }

    /* boundary conditions */
    fprintf(input, "%f,%f,%f\n", globalArgs.x1, globalArgs.y0,
            bound_cond_funcPtr(globalArgs.x1, globalArgs.y0));
    for (j = 0; j < dim; ++j) {
        y = globalArgs.y0 + (j + 1) * globalArgs.d;
        z = bound_cond_funcPtr(globalArgs.x1, y);
        fprintf(input, "%f,%f,%f\n", globalArgs.x1, y, z);
    }
    fprintf(input, "%f,%f,%f\n\n", globalArgs.x1, globalArgs.y1,
            bound_cond_funcPtr(globalArgs.x1, globalArgs.y1));

    return EXIT_SUCCESS;
}
