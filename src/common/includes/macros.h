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

#if HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef HAVE_LIBGSL
# include <gsl/gsl_math.h>
#endif /* HAVE_LIBGSL */


/** @defgroup Macros
 * @{
 */

#ifndef MAX
/** @def MAX
 * Given two numbers, return the bigger one.
 */
# ifdef HAVE_LIBGSL
#  define MAX(a, b) GSL_MAX(a, b)
# else
#  define MAX(a, b) ((a) > (b) ? (a) : (b))
# endif /* HAVE_LIBGSL */
#endif

#ifndef MIN
/** @def MIN
 * Given two numbers, return the smaller one.
 */
# ifdef HAVE_LIBGSL
#  define MIN(a, b) GSL_MIN(a, b)
# else
#  define MIN(a, b) ((a) < (b) ? (a) : (b))
# endif /* HAVE_LIBGSL */
#endif

/** @}*/



#endif /* __MACROS_H__ */
