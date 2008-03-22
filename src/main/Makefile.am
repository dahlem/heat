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

pdepoiss_solv_SOURCES =                                               \
	gnuplot.c                                                     \
	main.c                                                        \
	stats.c                                                       \
	cl.c

pdepoiss_solv_CPPFLAGS =                                              \
	-I./includes                                                  \
	-I$(top_srcdir)/src/solv/includes                             \
	-I$(top_srcdir)/src/poiss/2d/includes                         \
	-I$(top_srcdir)/src/common/includes


LDADD =                                                               \
	-lm                                                           \
	$(top_srcdir)/src/solv/libsolvers.la                          \
	$(top_srcdir)/src/poiss/2d/libpoiss2d.la                      \
	$(top_srcdir)/src/common/libcommon.la


include_HEADERS = includes
MAINTAINERCLEANFILES = Makefile.in

# compile with mpi
CC = $(MPICC)
LIBS += $(MPILIBS)
