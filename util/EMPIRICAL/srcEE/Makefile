#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

include ../../../Makefile.def

SEARCHDIR = -I${SHAREDIR} -I${MAGNETOGRAMDIR}

include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libEMPIRICALEE.a

include Makefile.DEPEND
#include Makefile.RULES

OBJECTS = \
	EEE_ModCommonVariables.o \
	EEE_ModMain.o \
	EEE_ModTD99.o \
	EEE_ModGL98.o \
	EEE_ModArch.o \
	EEE_ModGetB0.o \
	EEE_ModShearFlow.o \
	EEE_ModCms.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${OBJECTS}

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

distclean: clean
	rm -f Makefile.DEPEND
