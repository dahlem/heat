# Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

bin_PROGRAMS = pdepoiss_solv

if COMPILE_PARALLEL
pdepoiss_solv_SOURCES =                                               \
	main.c                                                        \
	mpi-common.c                                                  \
	cl.c
else
pdepoiss_solv_SOURCES =                                               \
	main.c                                                        \
	cl.c
endif

pdepoiss_solv_CPPFLAGS =                                              \
	-I./includes                                                  \
	-I$(top_srcdir)/src/solv/includes


LDADD =                                                               \
	$(top_srcdir)/src/solv/libsolvers.la


include_HEADERS = includes
MAINTAINERCLEANFILES = Makefile.in

# compile with mpi
CC = $(MPICC)
LIBS += $(MPILIBS)
