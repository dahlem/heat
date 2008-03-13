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

#include "cl.h"
#include "conjugate.h"
#include "error.h"
#include "matrix.h"
#include "mult.h"
#include "poiss_2d.h"
#include "vector.h"


void print_settings();
double src_dens(double x, double y);
double bound_cond(double x, double y);


int main(int argc, char *argv[])
{
    matrix A;
    vector u, v, x_bar;
    int i;
    int status;
    double square_pnts[2][2];

#ifdef HAVE_MPI
    setup(&argc, &argv);
#endif /* HAVE_MPI */

    /* parse the command-line arguments */
    if ((status = process_cl(argc, argv)) != 0) {
#ifdef HAVE_MPI
        finalise();
#endif /* HAVE_MPI */
        return status;
    }

    square_pnts[0][0] = globalArgs.x0;
    square_pnts[0][1] = globalArgs.x1;
    square_pnts[1][0] = globalArgs.y0;
    square_pnts[1][1] = globalArgs.y1;
    
    print_settings();

    /* set up the poisson matrix and vectors */
    setup_poiss_2d(&A, &u, &v, &x_bar, globalArgs.s, globalArgs.d,
                   square_pnts, *src_dens, *bound_cond);

#ifdef NDEBUG
    /* print the matrix and vectors */
#ifdef HAVE_MPI
    fprintf(stdout, "Partial matrix %d:\n", mpiArgs.rank);
#else
    printf("Matrix:\n");
#endif /* HAVE_MPI */

    for (i = 0; i < A.len; ++i) {
        vector_print(&(A.diags[i]));
    }
#endif /* NDEBUG */
    
    /* solve with conjugate gradient */
    conjugate(&A, &v, &u, &x_bar, globalArgs.e);

#ifdef NDEBUG
#ifdef HAVE_MPI
    fprintf(stdout, "Result from %d:\n", mpiArgs.rank);
#else
    printf("Result:\n");
#endif /* HAVE_MPI */
    vector_print(&x_bar);
#endif /* NDEBUG */

    matrix_free(&A);
    vector_free(&u);
    vector_free(&v);
    vector_free(&x_bar);

#ifdef HAVE_MPI
    finalise();
#endif /* HAVE_MPI */

    return EXIT_SUCCESS;
}

void print_settings()
{
#ifdef HAVE_MPI
    if (mpiArgs.rank == 0) {
#endif /* HAVE_MPI */
        fprintf(stdout, "(1) Application settings\n");
        fprintf(stdout, "Space Dimension  : %d\n", globalArgs.s);
        fprintf(stdout, "Time Dimension   : %d\n", globalArgs.t);
        fprintf(stdout, "Delta            : %1.4f\n\n", globalArgs.d);
        fprintf(stdout, "Input Range      : %2.2f <= x <= %2.2f; %2.2f <= y <= %2.2f\n",
                globalArgs.x0, globalArgs.x1, globalArgs.y0, globalArgs.y1);
        fprintf(stdout, "Error Threshold  : %e\n\n", globalArgs.e);

#ifdef HAVE_MPI
        fprintf(stdout, "(2) MPI settings\n");
        fprintf(stdout, "Number Processors : %d\n", mpiArgs.num_tasks);
    }
#endif /* HAVE_MPI */

    fprintf(stdout, "\n\n");
    fflush(stdout);
}

double src_dens(double x, double y)
{
    return - (4 * cos(x + y) * sin(x - y));
}

double bound_cond(double x, double y)
{
    return cos(x + y) * sin(x - y);
}

