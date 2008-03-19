/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file error.c
 * Implementation of the method declarations in error.h.
 *
 * @author Dominik Dahlem
 */
#include "error.h"
#include "vector.h"



void error(vector *x_bar, vector *x, vector *error)
{
    vector_alloc(error, x->len);
    vector_copy(error, x);
    scale(-1.0, error);
    add(error, x_bar);
    vector_abs(error);
}
