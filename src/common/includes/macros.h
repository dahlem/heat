/* Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>                */
/*                                                                             */
/* This file is free software; as a special exception the author gives         */
/* unlimited permission to copy and/or distribute it, with or without          */
/* modifications, as long as this notice is preserved.                         */
/*                                                                             */
/* This program is distributed in the hope that it will be useful, but         */
/* WITHOUT ANY WARRANTY, to the extent permitted by law; without even the      */
/* implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    */

/** @file macros.h
 * Declaration of some macros used in this application
 *
 * @author Dominik Dahlem
 */
#ifndef __MACROS_H__
#define __MACROS_H__


/** @defgroup Macros
 * @{
 */

#ifndef MAX
/** @def MAX
 * Given two numbers, return the bigger one.
 */
# define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
/** @def MIN
 * Given two numbers, return the smaller one.
 */
# define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

/** @}*/



#endif
