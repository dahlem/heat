/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file conjugate.h
 * Declaration of the methods for the Conjugate gradient algorithm.
 *
 * @author Dominik Dahlem
 */
#ifndef __CONJUGATE_H__
#define __CONJUGATE_H__


#include "cds_matrix.h"
#include "error.h"


/** @fn conjugate(cdsgb_matrix *A, vector *b, vector *x, vector **x_bar)
 *
 * This method solves the linear system with the Conjugate Gradient method.
 * The vector x_bar are allocated within this method and therefore
 * need to be de-allocated by the client application.
 *
 * @param cdsgb_matrix* the matrix A
 * @param vector* the vector b
 * @param vector* the vector x
 * @param vector** the solution vector x_bar
 * @return 0, if success. Otherwise the error code.
 */
int conjugate(cdsgb_matrix *A, vector *b, vector *x, vector **x_bar);



#endif
