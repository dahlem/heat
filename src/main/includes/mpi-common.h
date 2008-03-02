/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file mpi-common.h
 * Declaration of some structure and methods relevant for the parallel version of
 * the solver using MPI.
 *
 * @author Dominik Dahlem
 */
#ifndef __MPI_COMMON_H__
#define __MPI_COMMON_H__


/** @struct globalMpiArgs_t
 * A structure to capture the global MPI arguments.
 */
struct globalMpiArgs_t {
    int rank; /* MPI rank of the current process */
    int num_tasks; /* total number of MPI tasks */
} mpiArgs;


/** @fn void setup(int argc, char *argv[])
 * Display the help message for this application.
 */
void setup(int argc, char *argv[]);

/** @fn void finalise()
 * Initialise the global parameters.
 */
void finalise();


#endif
