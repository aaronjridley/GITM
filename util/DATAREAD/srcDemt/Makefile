#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

include ../../../Makefile.def
include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libDEMT.a

MODULES = \
	ModLdem.o

LIB:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MODULES} 
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES} 

distclean: clean

