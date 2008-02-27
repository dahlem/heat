/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file cds_matrix.h
 * This header file declares the vector and matrix structure including their
 * allocation and clean-up procedures.
 *
 * @author Dominik Dahlem
 */

#ifndef __CDS_MATRIX_H
#define __CDS_MATRIX_H


/** @typedef vector
 * A vector structure to represent double vectors.
 */
typedef struct
{
    int len; /** the length of the vector */
    double *data; /** the data array of doubles */
} vector;

/** @typedef cds_matrix
 * This structure represents a symmetric banded matrix structure in
 * compressed diagonal storage. Only the main diagonal and the upper diagonals
 * are stored in exactly the order they appear in the matrix.
 */
typedef struct 
{
    int len; /** the number of diagonals */
    vector *diags; /** the data array of vectors for each diagonal */
} cds_matrix;

/** @typedef cdsgb_matrix
 * This structure represents a generic banded matrix structure in
 * compressed diagonal storage. In contrast to the symmetric case, this storage
 * format requires leading zeros for lower diagonals and trailing zeros for upper
 * diagonals. As a consequence all diagonals have the same size. Additionally,
 * this storage format requires the specification of the relative offsets of each
 * diagonal to the main diagonal, where the main diagonal is 0.
 */
typedef struct 
{
    int len; /** the number of diagonals*/
    int *index; /** the data array of indeces relative to the main diagonal */
    vector *diags; /** the data array of vectors for each diagonal */
} cdsgb_matrix;



/** @fn void vector_alloc(vector *vec, int len)
 * Allocate memory for the vector structure
 *
 * @param vector* the vector to be allocated
 * @param int the length of the vector
 */
void vector_alloc(vector *vec, int len);

/** @fn void vector_free(vector *vec)
 * Free the vector structure
 *
 * @param vector* the vector to be free'd
 */
void vector_free(vector *vec);

/** @fn void cds_matrix_alloc(cds_matrix *mat, int len, int *elem)
 * Allocate memory for the symmetric banded matrix in CDS format given
 * the number of diagonals and the number of elements within each
 * diagonal.
 *
 * @param cds_matrix* the matrix to be allocated
 * @param int the number of diagonals
 * @param int* the number of elements for each diagonal
 */
void cds_matrix_alloc(cds_matrix *mat, int len, int *elem);

/** @fn void cds_matrix_free(cds_matrix *mat)
 * Free the matrix structure in CDS format.
 *
 * @param cds_matrix* the matrix structure to be free'd
 */
void cds_matrix_free(cds_matrix *mat);

/** @fn void cdsgb_matrix_alloc(cdsgb_matrix *mat, int len, int *elem, int *index)
 * Allocate memory for a generic banded matrix structure in CDS format given
 * the number of diagonals, their respective elements and their relative offset
 * to the main diagonal.
 *
 * @param cdsgb_matrix* the matrix structure to be allocated
 * @param int the number of diagonals
 * @param int* the number of elements for each diagonal
 * @param int* the relative offsets to teh main diagonal
 */
void cdsgb_matrix_alloc(cdsgb_matrix *mat, int len, int *elem, int *index);

/** @fn void cdsgb_matrix_free(cdsgb_matrix *mat)
 * Free the generic banded matrix structure.
 *
 * @param cdsgb_matrix* the matrix structure to be free'd
 */
void cdsgb_matrix_free(cdsgb_matrix *mat);


#endif
