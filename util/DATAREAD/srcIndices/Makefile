#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

include ../../../Makefile.def
include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libINDICES.a

include Makefile.DEPEND
#include Makefile.RULES

OBJECTS = \
	ModIndicesInterfaces.o\
	ModIndices.o\
	read_NGDC_indices.o\
	read_SWPC_indices.o\
	read_NOAAHPI_indices.o\
	read_MHDIMF_indices.o\
	insert_into_indices_array.o\
	indices_set_inputs.o\
	indices_library.o\
	IO_library.o\
	IO_set_inputs.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl -I${SHAREDIR} ${OBJECTS}

LIB:    DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

distclean: clean
	rm -f Makefile.DEPEND
