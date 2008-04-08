# Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>
#
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without
# modifications, as long as this notice is preserved.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
ACLOCAL_AMFLAGS = -I m4 -I /usr/share/aclocal

include mk/doxygen.am


SUBDIRS = src doc

# go into test
if HAVE_CUNIT
SUBDIRS += test
endif


MAINTAINERCLEANFILES = "configure Makefile.in semantic.cache"
MOSTLYCLEANFILES = ${DX_CLEANFILES}

EXTRA_DIST = \
	$(srcdir)/doc/coverage/ \
	$(srcdir)/m4/ \
	$(srcdir)/mk/ \
	$(srcdir)/max_error.awk \
	$(srcdir)/poiss.pbs

dist-hook:
	rm -rf `find $(distdir)/doc -name .gitignore`


if GCOV

COVOUT=./doc/coverage
LCOVOUT=$(COVOUT)/lcov-reports
GENHTMLOUT=$(LCOVOUT)/html

lcov-prepare:
	rm -rf $(COVOUT)
	mkdir -p $(GENHTMLOUT)
	lcov --zerocounters --directory . --quiet

lcov: clean lcov-prepare all
	cd test && $(MAKE) $(AM_MAKEFLAGS) lcov

clean-local:
	rm -rf doc/doxygen/html
	rm -rf $(LCOVOUT)

else
lcov: ;

clean-local:
	rm -rf doc/doxygen/html

endif



.PHONY: lcov lcov-prepare doxygen-run doxygen-doc $(DX_PS_GOAL) $(DX_PDF_GOAL)
