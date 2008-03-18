/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

#include "util.h"
#include "vector.h"



void gatherErrors(vector *x_bar, vector *x, vector **x_error)
{
    double elem;
    size_t i;

    vector_alloc(*x_error, x_bar->len);
    vector_copy(*x_error, x);
    gsl_vector_sub(*x_error, x_bar);
}