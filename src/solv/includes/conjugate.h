/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@gmail.com>                */
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


#include "matrix.h"
#include "vector.h"


/** @fn int conjugate(matrix *A, vector *b, vector *x, vector *x_bar, double err_thres)
 *
 * This method solves the linear system with the Conjugate Gradient method.
 *
 * @param matrix* the matrix A
 * @param vector* the vector b
 * @param vector* the vector x
 * @param vector* the solution vector x_bar
 * @param double the error threshold
 * @return 0, if success. Otherwise the error code.
 */
int conjugate(matrix *A, vector *b, vector *x, vector *x_bar, double err_thres);



#endif
