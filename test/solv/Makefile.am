# Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

# Rules for the test code (use `make check` to execute)
include $(top_srcdir)/mk/gcov.am

TESTS = solv_test
check_PROGRAMS = $(TESTS)

solv_test_SOURCES =                                                  \
	solv_suite.c                                                 \
	mult_test.c


solv_test_CPPFLAGS =                                                 \
	-I./includes                                                 \
	-I$(top_srcdir)/src/solv/includes


LDADD =                                                              \
	-lcunit                                                      \
	$(GSL_LIBS)                                                  \
	$(top_srcdir)/src/solv/libsolvers.la


include_HEADERS = includes


MAINTAINERCLEANFILES = Makefile.in
CLEANFILES = $(TESTS) *.log
