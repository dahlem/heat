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
/**
 * Default space dimension
 */
#define DEFAULT_SPACE_DIMENSION 20

/**
 * Default time dimension
 */
#define DEFAULT_TIME_DIMENSION 20

/**
 * Default error threshold
 */
#define DEFAULT_ERROR 1e-12

/** @}*/


/** @struct globalArgs_t
 * A structure to capture the global arguments passed into application on the
 * command-line.
 */
struct globalArgs_t {
    unsigned int s; /* space dimension */
    unsigned int t; /* time dimension */
    double e; /* error threshold */
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
 */
void process_cl(int argc, char **argv);



#endif
