# Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

AC_INIT([563_pde_poisson],
        [0.5],
        [Dominik Dahlem Dominik.Dahlem@cs.tcd.ie])


AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADER([config.h])

# extra autoconf macros
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-W -Wall foreign])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

LF_CONFIGURE_CC
LF_HOST_TYPE
LF_SET_WARNINGS

AC_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([sqrt])

# GSL library checks
# enable gsl?
AC_ARG_ENABLE(gsl, AS_HELP_STRING(
		     [--enable-gsl], 
		     [enable the GSL library @<:@default=no@:>@]), 
		     [GSLLIB=$enableval], [GSLLIB=no])

if test x$GSLLIB = xyes; then
   AC_CHECK_LIB([m],[cos])
   AC_CHECK_LIB([gslcblas],[cblas_dgemm])
   AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
fi
AM_CONDITIONAL(HAVE_GSL, test "x$GSLLIB" = xyes)


# enable cunit?
AC_ARG_ENABLE(test, AS_HELP_STRING(
		     [--enable-test], 
		     [enable test @<:@default=no@:>@]), 
		     [TESTING=$enableval], [TESTING=no])

# Check for cunit
if test x$TESTING = xyes; then
   AC_CHECK_LIB(cunit, CU_initialize_registry, HAVE_CUNIT="yes")
   
   if test x$HAVE_CUNIT != xyes; then
	AC_MSG_WARN([CUnit not found. Tests cannot be run.])
   fi
   if test "x$GSLLIB" = xno; then
	AC_MSG_ERROR([The CUnit tests require GSL to be configured.])
   fi
fi
AM_CONDITIONAL(HAVE_CUNIT, test "x$HAVE_CUNIT" = xyes)

# Check for MPI
AC_ARG_ENABLE(mpi, AS_HELP_STRING(
		   [--enable-mpi], 
		   [enable MPI @<:@default=no@:>@]), 
		   [USE_MPI=$enableval], [USE_MPI=no])
if test "x$USE_MPI" = xyes; then
   ACX_MPI([], [AC_MSG_WARN([cannot configure MPI (compile in serial)])])

   if test "x$HAVE_CUNIT" = xyes; then
	AC_MSG_WARN([CUnit is going to be disabled, if MPI is enabled.])
        AM_CONDITIONAL(HAVE_CUNIT, false)
   fi
else
   MPICC="${CC}"
   AC_SUBST(MPICC)
fi
AM_CONDITIONAL(COMPILE_PARALLEL, test "x$USE_MPI" = xyes)

# Check for OpenMP
AC_ARG_ENABLE(openmp, AS_HELP_STRING(
		   [--enable-openmp], 
		   [enable OpenMP @<:@default=no@:>@]), 
		   [USE_OPENMP=$enableval], [USE_OPENMP=no])
if test "x$USE_OPENMP" = xyes; then
   AX_OPENMP([], [AC_MSG_WARN([cannot configure OpenMP])])
   CFLAGS="$CFLAGS $OPENMP_CFLAGS"
fi
AM_CONDITIONAL(COMPILE_MULTITHREADED, test "x$USE_OPENMP" = xyes)

# Check for debug
AC_ARG_ENABLE(debug, AS_HELP_STRING(
		     [--enable-debug], 
		     [enable debug @<:@default=no@:>@]), 
		     [DEBUGGING=$enableval], [DEBUGGING=no])
if test "x$DEBUGGING" = xyes; then
  AC_DEFINE(NDEBUG, 1, [define debugging, if enabled])
fi

# enable report generation
AC_ARG_ENABLE(report, AS_HELP_STRING(
                     [--enable-report],
                     [enable report generation @<:@default=no@:>@]), 
                     [REPORT=$enableval], [REPORT=no])

# Check for latex stuff
if test x$REPORT = xyes; then
  AC_PROG_LATEX
  AC_PROG_DVIPS
  AC_PROG_PS2PDF
fi
AM_CONDITIONAL(GEN_REPORT, test "x$REPORT" = xyes)

# enable plot generation
AC_ARG_ENABLE(plot, AS_HELP_STRING(
                     [--enable-plot],
                     [enable plot generation @<:@default=no@:>@]), 
                     [PLOT=$enableval], [REPORT=no])

# Check for latex stuff
if test x$PLOT = xyes; then
  ACLTX_PROG_GNUPLOT
fi
AM_CONDITIONAL(GEN_PLOT, test "x$PLOT" = xyes)


AC_GCOV


DX_DOXYGEN_FEATURE([ON])
DX_HTML_FEATURE([ON])
DX_CHM_FEATURE([OFF])
DX_CHI_FEATURE([OFF])
DX_MAN_FEATURE([OFF])
DX_RTF_FEATURE([OFF])
DX_XML_FEATURE([OFF])
DX_PDF_FEATURE([OFF])
DX_PS_FEATURE([OFF])
DX_INIT_DOXYGEN([Poisson-Solver], [Doxyfile], [doc/doxygen])


AC_CONFIG_FILES([
   Makefile
   doc/Makefile
   doc/report/Makefile
   doc/eval/Makefile
   src/Makefile
   src/common/Makefile
   src/poiss/Makefile
   src/poiss/2d/Makefile
   src/solv/Makefile
   src/main/Makefile
   test/Makefile
   test/solv/Makefile
   test/common/Makefile
])

AC_OUTPUT


for top_builddir in . .. ../.. $ac_auxdir $ac_auxdir/..; do
  test -f $top_builddir/configure && break
done

AC_MSG_RESULT([

$PACKAGE_NAME $VERSION: Automatic configuration OK.

Type 'make' to compile this project.

])
