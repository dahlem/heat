# Copyright (C) 2008 Dominik Dahlem <Dominik.Dahlem@cs.tcd.ie>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

report_src = report.tex

report_images =                                                      \
	images/exact.eps                                             \
	images/poiss_1e-12.eps                                       \
	images/poiss.eps

EXTRA_DIST = $(report_src) report.pdf

images/%.eps: ../eval/%.eps
	( cp $? $@ && rm -f report.dvi; )

all-local: $(report_images) report.pdf

clean-local:
	rm -f *.aux *.log *.pdf *.dvi *.ps

dvi-local: report.dvi

report.dvi: $(report_src)
	$(latex) $<
	$(latex) $<

%.ps: %.dvi
	( unset LANG; $(dvips) $< -o $@ || { rm -f $@ ; false; }; )

%.pdf: %.ps
	( unset LANG; $(ps2pdf) $< $@ || { rm -f $@ ; false; } )
