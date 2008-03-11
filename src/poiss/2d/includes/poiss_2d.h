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
 * \f$ Au = v \f$ given the boundary conditions.
 *
 * @author Dominik Dahlem
 */
#ifndef __POISS2D_H__
#define __POISS2D_H__

#include "matrix.h"
#include "vector.h"


void setup_poiss_2d(matrix *A, vector *u, vector *v, vector *x_bar, int dim);


#endif
