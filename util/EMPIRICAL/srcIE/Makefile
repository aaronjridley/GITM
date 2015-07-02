#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

include ../../../Makefile.def
include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libEMPIRICALIE.a

include Makefile.DEPEND
#include Makefile.RULES

OBJECTS = \
	ModEIEConductance.o\
	ModAMIE_Interface.o\
	ModErrors.o\
	ModEIE_Interface.o\
	ModMHD_Interface.o\
	ModExtras.o\
	ModEIEFiles.o\
	ED_ModInterface.o\
	ED_ModIons.o\
	ModWeimer05.o\
	ModWeimer.o\
	AMIE_Library.o\
	EIE_End.o\
	EIE_Initialize.o\
	EIE_Library.o\
	EIE_IoLibrary.o\
	EIE_UaLibrary.o\
	EIE_set_inputs.o\
	readAMIEoutput.o\
	readMHDoutput.o\
	mh86.o\
	hmr89.o\
	iz94.o\
	ihp.o\
	merge_str.o\
	MHD_Library.o\
	ED_Interface.o \
	ED_ReadIonHeat.o \
	ED_r_elec_ed_lup_subs.o

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
