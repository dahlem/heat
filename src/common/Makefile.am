# Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

noinst_LTLIBRARIES = libcommon.la

if COMPILE_PARALLEL
libcommon_la_SOURCES =                                                \
	vector.c                                                      \
	matrix.c                                                      \
	mult.c                                                        \
	mpi-common.c
else
libcommon_la_SOURCES =                                                \
	vector.c                                                      \
	matrix.c                                                      \
	mult.c
endif

libcommon_la_CPPFLAGS =                                               \
	-I./includes

if HAVE_GSL
libcommon_la_LIBADD =                                                \
	$(GSL_LIBS)
else
libcommon_la_LIBADD =                                                \
	-lm
endif


include_HEADERS = includes
MAINTAINERCLEANFILES = Makefile.in

# compile with mpi
CC = $(MPICC)
LIBS += $(MPILIBS)
