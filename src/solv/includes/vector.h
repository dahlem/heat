/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file vector.h
 * This header file declares methods related to the vector structure.
 *
 * @author Dominik Dahlem
 */
#ifndef __VECTOR_H
#define __VECTOR_H


/** @typedef vector
 * A vector structure to represent double vectors.
 */
typedef struct
{
    int len; /** the length of the vector */
    double *data; /** the data array of doubles */
} vector;


/** @fn void vector_alloc(vector *vec, int len)
 * Allocate memory for the vector structure
 *
 * @param vector* the vector to be allocated
 * @param int the length of the vector
 */
void vector_alloc(vector *vec, int len);

/** @fn void vector_calloc(vector *vec, int len)
 * Allocate memory for the vector structure and initialise the vector elements to zero.
 *
 * @param vector* the vector to be allocated
 * @param int the length of the vector
 */
void vector_calloc(vector *vec, int len);

/** @fn void vector_free(vector *vec)
 * Free the vector structure
 *
 * @param vector* the vector to be free'd
 */
void vector_free(vector *vec);

/** @fn void zero(const vector *vec)
 * Zero each element of the given vector.
 *
 * @param vector* the vector to be set to (0, 0, ..., 0).
 */
void zero(const vector *vec);

/** @fn double dotProduct(const vector *const a, const vector *const b)
 * This function calculates the dot product of two given vectors.
 *
 * @param const vector *const the left vector
 * @param const vector *const the right vector
 * @return the dot product result
 */
double dotProduct(const vector *const a, const vector *const b);

/** @fn void daxpy(double alpha, const vector *const x, const vector *y)
 * Calculate \f$ y = \alpha * x \f$.
 *
 * @param double the scale parameter
 * @param const vector *const the vector to be scaled
 * @param const vector* the vector which holds the result
 */
void daxpy(double alpha, const vector *const x, const vector *y);

/** @fn double dnrm2(const vector *const x)
 * This function calculates the euclidean norm \f$ ||x||_2 = \sqrt {\sum x_i^2} \f$.
 *
 * @param const vector *const the vector to calculate the norm of
 * @return the euclidean norm of the vector x
 */
double dnrm2(const vector *const x);

/** @fn void add(const vector *x, const vector *const y)
 * This function adds the elements of vector x to the elements of vector y,
 * \f$ x'_i = x_i + y_i \f$. The two vectors must have the same length.
 *
 * @param const vector* the left vector containing the result
 * @param const vector *const the right vector to be added to the left one.
 * @return 0, if success. Otherwise the error code.
 */
int add(const vector *x, const vector *const y);



#endif
