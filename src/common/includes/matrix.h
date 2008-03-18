/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file matrix.h
 * This header file declares the matrix structure including its
 * allocation and clean-up procedures.
 *
 * @author Dominik Dahlem
 */

#ifndef __MATRIX_H
#define __MATRIX_H


#include "vector.h"

/** @enum matrix_type
 * This enumeration identifies the supported matrix types.
 */
enum matrix_type {
    SB, /** Symmetric Banded matrix type */
    GB  /** Generic Banded matrix type */
};

/** @enum matrix_storage
 * This enumeration identifies the supported matrix storage formats.
 */
enum matrix_storage {
    CDS  /** Compressed Diagonal Storage format */
};

/** @typedef matrix
 * This structure represents a matrix which can be represented as a generic banded
 * and symmetric banded matrix. The matrix storage format supported is the
 * compressed diagonal storage. In contrast to the symmetric case, the generic storage
 * format requires leading zeros for lower diagonals and trailing zeros for upper
 * diagonals. As a consequence all diagonals have the same size. Additionally,
 * the CDS storage format requires the specification of the relative offsets of each
 * diagonal to the main diagonal, where the main diagonal is 0.
 */
typedef struct
{
    int len; /** the number of diagonals*/
    int *index; /** the data array of indeces relative to the main diagonal */
    vector *diags; /** the data array of vectors for each diagonal */
    enum matrix_type mtype; /** the matrix type */
    enum matrix_storage storage; /** the storage format */
} matrix;


/** @fn void cds_matrix_alloc(matrix *mat, int len, int *elem)
 * Allocate memory for the symmetric banded matrix in CDS format given
 * the number of diagonals and the number of elements within each
 * diagonal.
 *
 * @param matrix* the matrix to be allocated
 * @param int the number of diagonals
 * @param int* the number of elements for each diagonal
 */
void cds_matrix_alloc(matrix *mat, int len, int *elem);

/** @fn void cds_matrix_free(matrix *mat)
 * Free the matrix structure in CDS format.
 *
 * @param matrix* the matrix structure to be free'd
 */
void cds_matrix_free(matrix *mat);

/** @fn void cdsgb_matrix_alloc(matrix *mat, int len, int *elem, int *index)
 * Allocate memory for a generic banded matrix structure in CDS format given
 * the number of diagonals, their respective elements and their relative offset
 * to the main diagonal.
 *
 * @param matrix* the matrix structure to be allocated
 * @param int the number of diagonals
 * @param int* the number of elements for each diagonal
 * @param int* the relative offsets to teh main diagonal
 */
void cdsgb_matrix_alloc(matrix *mat, int len, int *elem, int *index);

/** @fn void cdsgb_matrix_free(matrix *mat)
 * Free the generic banded matrix structure.
 *
 * @param matrix* the matrix structure to be free'd
 */
void cdsgb_matrix_free(matrix *mat);

/** @fn void matrix_alloc(matrix *mat)
 * Function to allocate the matrix. The matrix characteristics of the matrix are required
 * and have to be set within the structure before passing the matrix into this method.
 *
 * @param matrix* the matrix to be allocated
 */
void matrix_alloc(matrix *mat, int *elem, int *index);

/** @fn void matrix_free(matrix *mat)
 * Free the allocated memory for the matrix.
 */
void matrix_free(matrix *mat);

/** @fn void matrix_print(matrix *mat)
 * Print a given matrix to stdout.
 *
 * @param matrix* the matrix to be printed
 */
void matrix_print(matrix *mat);


#endif
