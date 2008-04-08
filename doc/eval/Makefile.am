# Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

eval_src =                                                           \
	common.gnu                                                   \
	poiss.gnu                                                    \
	poiss-smooth.gnu                                             \
	poiss_1e-12.gnu                                              \
	exact.gnu                                                    \
	serial.dat                                                   \
	serial_1e-12.dat                                             \
	parallel.dat

EXTRA_DIST = $(eval_src)

all-local: exact.eps poiss.eps poiss_1e-12.eps

clean-local:
	rm -f *.eps

%.eps: $(eval_src)
	( $(gnuplot) $*.gnu || { rm -f $@ ; false; }; )
