/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file poiss_2d.h
 * Declaration of the methods to setup up the 2D poisson equations in the form
 * \f$ Au = v \f$ given the boundary conditions. Using the finite difference scheme
 * to solve the 2D poisson equation with an \f$ m x n \f$ grid, where \f$ m = n \f$, i.e.
 * there are equal grid points in the spatial and time direction.
 *
 * @author Dominik Dahlem
 */
#ifndef __POISS2D_H__
#define __POISS2D_H__

#include "matrix.h"
#include "vector.h"


/** @def POISS_2D_GB_MATRIX_DIAGS
 * The number of diagonals to be stored for the general banded form of the poisson
 * equation.
 */
#define POISS_2D_GB_MATRIX_DIAGS 5

/** @def POISS_2D_SB_MATRIX_DIAGS
 * The number of diagonals to be stored for the symmetric banded form of the poisson
 * equation.
 */
#define POISS_2D_SB_MATRIX_DIAGS 3

/** @def D_MAIN_DIAG
 * Main diagonal of the \f$m x m\f$ diagonal matrix D
 */
#define D_MAIN_DIAG 4.0

/** @def D_DIAG_1
 * Diagonal with offset 1 of the \f$m x m\f$ diagonal matrix D
 */
#define D_DIAG_1    -1.0

/** @def I_DIAG
 * Main diagonal of the \f$m x m\f$ identity matrix I
 */
#define I_DIAG      -1.0



void setup_poiss_2d(matrix *A, vector *u, vector *v, vector *x_bar, int dim,
                    double delta, double coord[4][2],
                    double (*src_dens_funcPtr)(double, double),
                    double (*bound_cond_funcPtr)(double, double));


#endif
