/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */
#include "config.h"

#ifdef HAVE_MPI
# include <mpi.h>
# include "mpi-common.h"
#endif /* HAVE_MPI */

#include <stdlib.h>
#include <stdio.h>

#include "cl.h"




int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
    MPI_Status status;
    int mpi_ret;

    setup(argc, argv);
#endif /* HAVE_MPI */

    /* initialise and configure the command line options */
    init();
    process_cl(argc, argv);


#ifdef HAVE_MPI
    finalise();
#endif /* HAVE_MPI */

    return EXIT_SUCCESS;
}
