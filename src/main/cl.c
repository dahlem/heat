/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "cl.h"


/** @var cl_arguments
 * getopt configuration of the command-line parameters.
 * All command-line arguments are optional.
 * Arguments p, e, and m are passed in by the MPI environment, we ignore
 * those, but the have to be specified.
 */
static const char *cl_arguments = "uh?s:t:r:p:e:m:";



void displayHelp()
{
    printf("pdepoiss_solv - solve Poisson's equation.\n");
    printf(" -s : Space dimension.\n");
    printf(" -t : Time dimension.\n");
    printf(" -r : Error threshold.\n");
    printf(" -? : This help message.\n");
    printf(" -h : This help message.\n");

    exit(EXIT_SUCCESS);
}


void init()
{
    globalArgs.s = DEFAULT_SPACE_DIMENSION;
    globalArgs.t = DEFAULT_TIME_DIMENSION;
    globalArgs.e = DEFAULT_ERROR;
}


/** @fn void verify_cl()
 * Verify the command line arguments and check that they are consistent with each
 * other.
 */
void verify_cl()
{
    if (globalArgs.s < 3) {
        globalArgs.s = DEFAULT_SPACE_DIMENSION;
    }
    if (globalArgs.t < 3) {
        globalArgs.t = DEFAULT_TIME_DIMENSION;
    }

    /* make sure we have equal mesh-points in space and time */
    if (globalArgs.t < globalArgs.s) {
        globalArgs.t = globalArgs.s;
    } else if (globalArgs.t > globalArgs.s) {
        globalArgs.s = globalArgs.t;
    }
    if (globalArgs.e < 0) {
        globalArgs.e = DEFAULT_ERROR;
    }
}


void process_cl(int argc, char **argv)
{
    int opt = 0;

    opt = getopt(argc, argv, cl_arguments);

    while (opt != -1) {
        switch (opt) {
            case 't':
                globalArgs.t = atoi(optarg);
                break;
            case 's':
                globalArgs.s = atoi(optarg);
                break;
            case 'r':
                globalArgs.e = atof(optarg);
                break;
            case 'p':
                break;
            case 'e':
                break;
            case 'm':
                break;
            case 'h':
            case '?':
                 displayHelp();
                 break;
            default:
                 break;
        }
        opt = getopt(argc, argv, cl_arguments);
    }

    verify_cl();
}
