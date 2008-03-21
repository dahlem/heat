/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file cl.h
 * Declaration of some command-line arguments and a @code struct for
 * the command-line option values.
 *
 * @author Dominik Dahlem
 */
#ifndef __CL_H__
#define __CL_H__


/** @defgroup Application Definition of default application settings
 * @{
 */
/** @define DEFAULT_SPACE_DIMENSION
 * Default space dimension
 */
#define DEFAULT_SPACE_DIMENSION 26

/** @define DEFAULT_TIME_DIMENSION
 * Default time dimension
 */
#define DEFAULT_TIME_DIMENSION 26

/** @define DEFAULT_ERROR
 * Default error threshold
 */
#define DEFAULT_ERROR 1e-6

/** @def DEFAULT_X0
 * Lower bound of the range in the x dimension.
 */
#define DEFAULT_X0 -0.5

/** @def DEFAULT_X1
 * Upper bound of the range in the x dimension.
 */
#define DEFAULT_X1 2.0

/** @def DEFAULT_Y0
 * Lower bound of the range in the y dimension.
 */
#define DEFAULT_Y0 -2.0

/** @def DEFAULT_Y1
 * Upper bound of the range in the y dimension.
 */
#define DEFAULT_Y1 0.5

/** @def DEFAULT_FILENAME
 * Default file name for the result surface
 */
#define DEFAULT_FILENAME "./result.dat"

/** @}*/


/** @struct globalArgs_t
 * A structure to capture the global arguments passed into application on the
 * command-line.
 */
struct globalArgs_t {
    unsigned int s;         /* space dimension */
    unsigned int t;         /* time dimension */
    double e;               /* error threshold */
    double d;               /* the step size in the 2D grid (derived from the dimension) */
    double x0, x1, y0, y1;  /* input range of the 2D poisson equation */
    char *f;                /* filename for the result surface */
} globalArgs;


/** @fn void displayHelp()
 * Display the help message for this application.
 */
void displayHelp();

/** @fn void init()
 * Initialise the global parameters.
 */
void init();

/** @fn void process_cl(int argc, char **argv)
 * Process the command-line arguments passed into the application.
 *
 * @param int number of arguments
 * @param char** pointer to the character array representing the command-line
 *               parameters
 * @return 0, if success. Otherwise the error code.
 */
int process_cl(int argc, char **argv);



#endif
