/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@gmail.com>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file mpi-utils.h
 * Declaration of some utility functions relevant for the parallel version of
 * the solver using MPI.
 *
 * @author Dominik Dahlem
 */
#ifndef __MPI_UTILS_H__
#define __MPI_UTILS_H__


int adjustment(int sys_dim, int num_tasks);
int block(int sys_dim, int num_tasks);


#endif
