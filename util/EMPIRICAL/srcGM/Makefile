#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

include ../../../Makefile.def
include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libEMPIRICALGM.a

#include Makefile.RULES

OBJECTS = \
        ModTsyganenko.o

.f.o:
	${COMPILE.f77} ${CFLAGS} ${PRECISION} $<
	${MOVE_MOD}

LIB:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

distclean: clean
