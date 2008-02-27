/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */
#ifndef __MULT_H
#define __MULT_H

#include "cds_matrix.h"


/**
 * This operation performs the matrix-vector multiply using a matrix in compressed
 * diagonal storage format. The multiplication algorithm assumes a
 * symmetric banded matrix with n diagonals.
 *
 * @param const cds_matrix *const the matrix in CDS format
 * @param const vector *const the u vector
 * @param const vector* the result vector
 */
void dcdssbmv(const cds_matrix *const mat, const vector *const u, const vector *v);

/**
 * This operation performs the matrix-vector multiply using a matrix in compressed
 * diagonal storage format. The multiplication algorithm assumes a
 * generic banded matrix with n diagonals. The off-main diagonals have leading
 * and trailing zeros for the lower and upper diagonals respectively.
 *
 * @param const cdsgb_matrix *const the matrix in CDS format
 * @param const vector *const the u vector
 * @param const vector* the result vector
 */
void dcdsgbmv(const cdsgb_matrix *const mat, const vector *const u, const vector *v);



#endif
