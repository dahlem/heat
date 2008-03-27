/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @mainpage Parallel Conjugate Gradient to solve Poisson's equation
 * The main file of the pde solver for the poisson equation
 * \f$ \bigtriangledown^2u = -f \f$ with \f$ u \equiv u(x,y) \f$ on the square region
 * \f$ ABCD, A=(-0.5, -2), B=(2, -2), C=(2, 0.5), D=(-0.5, 0.5) \f$, where the
 * source density is given by
 * \f$ f(x,y) = 4\cos{x+y}\sin{x-y} \f$.
 * The exact solution and boundary conditions are given by the formula
 * \f$ u = \cos{x+y}\sin{x-y} \f$.
 *
 * @author Dominik Dahlem
 */
#if HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef HAVE_MPI
# include <mpi.h>
# include "mpi-common.h"
#endif /* HAVE_MPI */

#ifdef HAVE_OPENMP
# include <omp.h>
#endif /* HAVE_OPENMP */

#ifdef HAVE_LIBGSL
# include <gsl/gsl_math.h>
# include <gsl/gsl_ieee_utils.h>
#endif /* HAVE_LIBGSL */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "cl.h"
#include "conjugate.h"
#include "error.h"
#include "gnuplot.h"
#include "matrix.h"
#include "poiss_2d.h"
#include "vector.h"


/** @fn void print_settings()
 * Print the settings of the current application run.
 */
void print_settings();

/** @fn double src_dens(double x, double y)
 * This function declares the source density
 * \f$ f(x,y) = 4\cos{x+y}\sin{x-y} \f$.
 *
 * @param double the x-value
 * @param double the y-value
 * @return the value of the function at the given points.
 */
double src_dens(double x, double y);

/** @fn double bound_cond(double x, double y)
 * This function declares the boundary condition given by
 * \f$ u = \cos{x+y}\sin{x-y} \f$.
 *
 * @param double the x-value
 * @param double the y-value
 * @return the value of the function at the given points.
 */
double bound_cond(double x, double y);


int main(int argc, char *argv[])
{
    matrix A;
    vector u, v, x_bar;
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

#ifdef HAVE_LIBGSL
    /* read GSL_IEEE_MODE */
    gsl_ieee_env_setup();
#endif /* HAVE_LIBGSL */

    /* specify the domain */
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
# ifdef HAVE_MPI
    fprintf(stdout, "Partial matrix %d:\n", mpiArgs.rank);
# else
    printf("Matrix:\n");
# endif /* HAVE_MPI */
    matrix_print(&A);
#endif /* NDEBUG */

    /* solve with conjugate gradient */
    if (status == 0) {
        status = conjugate(&A, &v, &u, &x_bar, globalArgs.e);
    }

    /* print the vector into a gnuplot format */
    if (status == 0) {
        status = print_surface(&x_bar, globalArgs.s - 2, *bound_cond);
    }

    matrix_free(&A);
    vector_free(&u);
    vector_free(&v);
    vector_free(&x_bar);

#ifdef HAVE_MPI
    finalise();
#endif /* HAVE_MPI */

    return status;
}

void print_settings()
{
#ifdef HAVE_MPI
    if (mpiArgs.rank == 0) {
#endif /* HAVE_MPI */
        fprintf(stdout, "(1) Application settings\n");

#ifdef HAVE_LIBGSL
        fprintf(stdout, "GSL configured        : true\n");
#else
        fprintf(stdout, "GSL configured        : false\n");
#endif /* HAVE_LIBGSL */

#ifdef HAVE_OPENMP
        fprintf(stdout, "OpenMP                : true\n");
        fprintf(stdout, "Max number of Threads : %d\n", omp_get_max_threads());
        fprintf(stdout, "Support Nesting (0/1) : %d\n\n", omp_get_nested());
#else
        fprintf(stdout, "OpenMP                : false\n");
#endif /* HAVE_OPENMP */

#ifdef NDEBUG
        fprintf(stdout, "Debug                 : true\n\n");
#else
        fprintf(stdout, "Debug                 : false\n\n");
#endif /* NDEBUG */

        fprintf(stdout, "(2) Mesh settings\n");
        fprintf(stdout, "Space Dimension       : %d\n", globalArgs.s);
        fprintf(stdout, "Time Dimension        : %d\n", globalArgs.t);
        fprintf(stdout, "Delta                 : %1.8f\n", globalArgs.d);
        fprintf(stdout, "Input Range           : %2.2f <= x <= %2.2f; %2.2f <= y <= %2.2f\n\n",
                globalArgs.x0, globalArgs.x1, globalArgs.y0, globalArgs.y1);

        fprintf(stdout, "(3) Conjugate Gradient settings\n");
        fprintf(stdout, "Error Threshold       : %e\n\n", globalArgs.e);

#ifdef HAVE_MPI
        fprintf(stdout, "(4) MPI settings\n");
        fprintf(stdout, "Number Processors     : %d\n", mpiArgs.num_tasks);
    }
#endif /* HAVE_MPI */

    fprintf(stdout, "\n\n");
    fflush(stdout);
}

double src_dens(double x, double y)
{
    return 4.0 * cos(x + y) * sin(x - y);
}

double bound_cond(double x, double y)
{
    return cos(x + y) * sin(x - y);
}
