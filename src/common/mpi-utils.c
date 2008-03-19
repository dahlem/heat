/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file mpi-utils.c
 * Implementation of the methods declared in mpi-utils.h
 *
 * @author Dominik Dahlem
 */
#include "mpi-utils.h"


int adjustment(int sys_dim, int num_tasks)
{
    return sys_dim % num_tasks;
}

int block(int sys_dim, int num_tasks)
{
    int adjust;

    adjust = adjustment(sys_dim, num_tasks);
    
    return ((sys_dim + adjust) / num_tasks);
}
