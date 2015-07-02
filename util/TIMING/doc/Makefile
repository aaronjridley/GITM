#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
SHELL = /bin/sh

help:
	@echo 'Targets for the documentation of the TIMING utility:'
	@echo ' '
	@echo 'help              - this help (default target)'
	@echo ' '
	@echo 'MAN               - MAN_TIMING.ps.gz manual'
	@echo ' '
	@echo 'clean             - /bin/rm -f *.o *.exe'
	@echo 'distclean         - /bin/rm -f *.o *.exe *.a *~'

MAN:	MAN_TIMING.ps.gz

MAN_TIMING.dvi : MAN_TIMING.tex TIMING.tex
	latex MAN_TIMING.tex
	latex MAN_TIMING.tex

MAN_TIMING.ps : MAN_TIMING.dvi
	dvips -o MAN_TIMING.ps MAN_TIMING.dvi

MAN_TIMING.ps.gz : MAN_TIMING.ps
	gzip -f MAN_TIMING.ps

clean:
	/bin/rm -f *.log *.aux *.dvi *~

distclean: clean

