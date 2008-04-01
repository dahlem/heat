/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file cl.c
 * This file implements the command-line parsing as specified in the cl.h header file.
 *
 * @author Dominik Dahlem
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "cl.h"
#include "error.h"


/** @var cl_arguments
 * getopt configuration of the command-line parameters.
 * All command-line arguments are optional.
 */
static const char *cl_arguments = "uh?s:t:d:r:1:2:3:4:f:";



void displayHelp()
{
    printf("pdepoiss_solv - solve Poisson's equation.\n");
    printf(" -s : Space dimension (Note: specify grid dimensions or the delta value).\n");
    printf(" -t : Time dimension (Note: specify grid dimensions or the delta value).\n");
    printf(" -d : Delta value (Note: specify grid dimensions or the delta value).\n");
    printf(" -r : Error threshold.\n");
    printf(" -f : File name for the result surface.\n");
    printf(" -1 : Lower bound of the domain in the x dimension.\n");
    printf(" -2 : Upper bound of the domain in the x dimension.\n");
    printf(" -3 : Lower bound of the domain in the y dimension.\n");
    printf(" -4 : Upper bound of the domain in the y dimension.\n");
    printf(" -? : This help message.\n");
    printf(" -h : This help message.\n");

    exit(EXIT_SUCCESS);
}


void init()
{
    globalArgs.s = DEFAULT_SPACE_DIMENSION;
    globalArgs.t = DEFAULT_TIME_DIMENSION;
    globalArgs.d = DEFAULT_DELTA;
    globalArgs.e = DEFAULT_ERROR;
    globalArgs.f = DEFAULT_FILENAME;
    globalArgs.x0 = DEFAULT_X0;
    globalArgs.x1 = DEFAULT_X1;
    globalArgs.y0 = DEFAULT_Y0;
    globalArgs.y1 = DEFAULT_Y1;
}


/** @fn int verify_cl()
 * Verify the command line arguments and check that they are consistent with each
 * other.
 *
 * @return 0, if sucess. Otherwise the error code.
 */
int verify_cl()
{
    double temp;

    if (globalArgs.s < 3) {
        globalArgs.s = DEFAULT_SPACE_DIMENSION;
    }
    if (globalArgs.t < 3) {
        globalArgs.t = DEFAULT_TIME_DIMENSION;
    }

    /* make sure we have equal mesh-points in space and time */
    if (globalArgs.t < globalArgs.s) {
        fprintf(stderr, "Warning: We assume an equal number of mesh points for both dimensions.\n");
        fflush(stderr);
        globalArgs.t = globalArgs.s;
    } else if (globalArgs.t > globalArgs.s) {
        fprintf(stderr, "Warning: We assume an equal number of mesh points for both dimensions.\n");
        fflush(stderr);
        globalArgs.s = globalArgs.t;
    }
    if (globalArgs.e < 0) {
        fprintf(stderr, "Warning: The error threshold has to be possitive.\n\
Use the default threshold instead.\n");

        fflush(stderr);
        globalArgs.e = DEFAULT_ERROR;
    }

    if ((globalArgs.x1 - globalArgs.x0)
        != (globalArgs.y1 - globalArgs.y0)) {
        fprintf(stderr, "Error: We assume equal input ranges for both dimensions.\n");
        fflush(stderr);
        return GRID_DIM_MISMATCH;
    }

    /* grid dimensions take precedence over delta value configuration */
    if ((globalArgs.s == DEFAULT_SPACE_DIMENSION)
        && (globalArgs.t == DEFAULT_TIME_DIMENSION)
        && (globalArgs.d != DEFAULT_DELTA)) {
        /* add a small amount to nudge the double value just over a the next integer value */
        temp = (fabs(globalArgs.x1 - globalArgs.x0) / globalArgs.d) + 0.001;
        globalArgs.s = globalArgs.t = (temp + 1.0);
    }

    /* derive the delta value which we assume is equal in both grid dimensions. */
    if (((globalArgs.s != DEFAULT_SPACE_DIMENSION)
        && (globalArgs.t != DEFAULT_TIME_DIMENSION))
        || (globalArgs.d != DEFAULT_DELTA)) {
        globalArgs.d = fabs(globalArgs.x1 - globalArgs.x0) / (globalArgs.s - 1.0);
    }

    return EXIT_SUCCESS;
}


int process_cl(int argc, char **argv)
{
    int opt = 0;

    init();

    opt = getopt(argc, argv, cl_arguments);

    while (opt != -1) {
        fprintf(stdout, "%c\n", opt);
        fflush(stdout);
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
            case 'd':
                globalArgs.d = atof(optarg);
                break;
            case 'f':
                globalArgs.f = optarg;
                break;
            case '1':
                globalArgs.x0 = atof(optarg);
                break;
            case '2':
                globalArgs.x1 = atof(optarg);
                break;
            case '3':
                globalArgs.y0 = atof(optarg);
                break;
            case '4':
                globalArgs.y1 = atof(optarg);
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

    return verify_cl();
}
