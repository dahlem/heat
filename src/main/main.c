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
#include <stdio.h>

#include "cl.h"
#include "conjugate.h"
#include "matrix.h"
#include "mult.h"
#include "poiss_2d.h"
#include "vector.h"


void print_settings();


int main(int argc, char *argv[])
{
    matrix A;
    vector u, v, x_bar;
    int dim;
    int i;

#ifdef HAVE_MPI
    MPI_Status status;

    setup(&argc, &argv);

#endif /* HAVE_MPI */
    process_cl(argc, argv);
    
    print_settings();

    dim = 9;

    /* set up the poisson matrix and vectors */
    setup_poiss_2d(&A, &u, &v, &x_bar, dim);

    /* print the matrix and vectors */
#ifdef NDEBUG
#ifdef HAVE_MPI
    fprintf(stdout, "Partial matrix %d:\n", mpiArgs.rank);
#else
    printf("Matrix:\n");
#endif /* HAVE_MPI */

    for (i = 0; i < A.len; ++i) {
        vector_print(&(A.diags[i]));
    }
#endif /* NDEBUG */
    
    /* solv with conjugate gradient */
    conjugate(&A, &v, &u, &x_bar);

#ifdef HAVE_MPI
    fprintf(stdout, "Result from %d:\n", mpiArgs.rank);
#else
    printf("Result:\n");
#endif

    vector_print(&x_bar);

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
        fprintf(stdout, "Error Threshold  : %f\n\n", globalArgs.e);

#ifdef HAVE_MPI
        fprintf(stdout, "(2) MPI settings\n");
        fprintf(stdout, "Number Processors : %d\n", mpiArgs.num_tasks);
    }
#endif /* HAVE_MPI */

    fprintf(stdout, "\n\n");
    fflush(stdout);
}
