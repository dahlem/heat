# Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@gmail.com>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

noinst_LTLIBRARIES = libcommon.la

libcommon_la_SOURCES =                                                \
	vector.c                                                      \
	matrix.c                                                      \
	mult.c

if COMPILE_PARALLEL
libcommon_la_SOURCES +=                                               \
	mpi-utils.c                                                   \
	mpi-common.c
endif

libcommon_la_CPPFLAGS =                                               \
	-I./includes

libcommon_la_LIBADD =                                                 \
	-lm

if HAVE_GSL
libcommon_la_LIBADD +=                                                \
	$(GSL_LIBS)
endif


include_HEADERS = includes
MAINTAINERCLEANFILES = Makefile.in

# compile with mpi
CC = $(MPICC)
LIBS += $(MPILIBS)
