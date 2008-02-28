/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file mult.h
 * This header file declares the matrix-vector multiplication procedures.
 * The procedure names are modelled after the BLAS interfaces:
 * - d: elements stored as doubles
 * - cds: matrix storage format (Compressed Diagonal Storage)
 * - sb: symmetric band matrix
 * - gb: generic band matrix
 * - mv: matrix-vector operation
 *
 * @author Dominik Dahlem
 */
#ifndef __MULT_H
#define __MULT_H

#include "matrix.h"


/** @fn void dcdssbmv(const matrix *const mat, const vector *const u, const vector *v)
 *
 * This operation performs the matrix-vector multiply using a matrix in compressed
 * diagonal storage format. The multiplication algorithm assumes a
 * symmetric banded matrix with n diagonals. The algorithm does not assume any number
 * of diagonals. The offsets of each diagonal are calculated as the difference of the
 * current diagonal and the main diagonal.
 *
 * @param const matrix *const the matrix in CDS format
 * @param const vector *const the u vector
 * @param const vector* the result vector
 */
void dcdssbmv(const matrix *const mat, const vector *const u, const vector *v);

/** @fn void dcdsgbmv(const matrix *const mat, const vector *const u, const vector *v)
 *
 * This operation performs the matrix-vector multiply using a matrix in compressed
 * diagonal storage format. The multiplication algorithm assumes a
 * generic banded matrix with n diagonals. The off-main diagonals have leading
 * and trailing zeros for the lower and upper diagonals respectively. The algorithm
 * accepts any number of diagonals. The main difference to the @code dcdssbmv algorithm
 * is that an index has to be specified to declare the offsets of the left and the
 * right diagonals to the main diagonals. Left indeces are negative vs. right indeces
 * are positive. The algorithm is described on NetLib
 * http://www.netlib.org/linalg/html_templates/node99.html.
 *
 * @param const matrix *const the matrix in CDS format
 * @param const vector *const the u vector
 * @param const vector* the result vector
 */
void dcdsgbmv(const matrix *const mat, const vector *const u, const vector *v);

/** @fn void dgbmv(const matrix *const mat, const vector *const u, const vector *v)
 * This operation performs the matrix-vector multiply. It delegates the function call
 * to the respective specialised operation depending on the given matrix type.
 *
 * @param const matrix *const the matrix in CDS format
 * @param const vector *const the u vector
 * @param const vector* the result vector
 */
void dgbmv(const matrix *const mat, const vector *const u, const vector *v);



#endif
