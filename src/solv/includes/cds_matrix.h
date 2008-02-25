/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */
#ifndef __CDS_MATRIX_H
#define __CDS_MATRIX_H


typedef struct
{
    int len;
    double *data;
} vector;

typedef struct 
{
    int len;
    vector *diags;
} cds_matrix;


void vector_alloc(vector *vec, int len);
void vector_free(vector *vec);
void cds_matrix_alloc(cds_matrix *mat, int len, int *elem);
void cds_matrix_free(cds_matrix *mat);


#endif
