/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@gmail.com>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file error.h
 * Declaration of the error codes for this application.
 *
 * @author Dominik Dahlem
 */
#ifndef __ERROR_H__
#define __ERROR_H__


/** @defgroup Errors
 * @{
 */

/** @def MATRIX_VECTOR_UNEQUAL_ROW_DIM
 * Error code if the matrix-vector dimensions do not match
 */
#define MATRIX_VECTOR_UNEQUAL_ROW_DIM   11

/** @def VECTOR_DIMENSION_MISMATCH
 * Error code if the dimensions of two or more vectors do not match
 */
#define VECTOR_DIMENSION_MISMATCH       12

/** @def GRID_DIM_MISMATCH
 * Error code if the dimensions of grid points for the poisson equations
 * does not match
 */
#define GRID_DIM_MISMATCH               13

/** @def FILE_OPEN_FOR_WRITE_ERROR
 * Error code if a file cannot be opened with write permission
 */
#define FILE_OPEN_FOR_WRITE_ERROR       14

/** @def MALLOC_ERROR
 * Error code if a malloc error occurred (i.e., the memory could not be
 * allocated).
 */
#define MALLOC_ERROR                    15

/** @def CONVERGENCE_ERROR
 * Error code if the conjugate gradient method did not converge in the
 * maximum number of iterations.
 */
#define CONVERGENCE_ERROR               16

/** @}*/



#endif
