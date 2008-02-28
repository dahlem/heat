/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file error.h
 * Declaration of the error codes for the linear algebra solvers.
 *
 * @author Dominik Dahlem
 */
#ifndef __ERROR_H__
#define __ERROR_H__


/** @defgroup Errors
 * @{
 */

/** @def
 * Error code if the matrix-vector dimensions do not match
 */
#define MATRIX_VECTOR_UNEQUAL_ROW_DIM   11


/** @def
 * Error code if the dimensions of two or more vectors do not match
 */
#define VECTOR_DIMENSION_MISMATCH       12

/** @}*/



#endif
