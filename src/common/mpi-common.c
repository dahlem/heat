/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */
#if HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef HAVE_MPI
# include <mpi.h>
#endif /* HAVE_MPI */

#include "mpi-common.h"


void setup(int *argc, char **argv[])
{
    MPI_Init(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiArgs.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiArgs.num_tasks);
}

void finalise()
{
    MPI_Finalize();
}
