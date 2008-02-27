/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */
#ifndef __VECTOR_H
#define __VECTOR_H

#include "cds_matrix.h"


/**
 * Zero each element of the given vector.
 *
 * @param vector* the vector to be set to (0, 0, ..., 0).
 */
void zero(const vector *vec);


#endif
