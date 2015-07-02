#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
SHELL = /bin/sh

include ../../../Makefile.def
include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libTIMING.a

help:
	@echo 'Targets to make for the TIMING module:'
	@echo ' '
	@echo 'help              - this help (default target)'
	@echo ' '
	@echo 'test              - serial testing of TIMING module'
	@echo 'test_mpi          - parallel testing of TIMING module'         
	@echo 'test_empty        - serial testing of empty TIMING module'
	@echo 'test_empty_mpi    - parallel testing of empty TIMING module'
	@echo 'tests             - all four tests together'
	@echo ' '
	@echo 'LIB               - libTIMING.a the complete TIMING library'
	@echo ' '
	@echo 'clean             - /bin/rm -f *.o *.exe'
	@echo 'distclean         - /bin/rm -f *.o *.exe *.a *~'

OBJ_LIB = ModTiming.o timing_cpu.o timing.o

LIB: ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}:  $(OBJ_LIB) ${MAKEFILE_COMP_SELECT}
	@rm -f ${MY_LIB}
	${AR} ${MY_LIB} $(OBJ_LIB)

tests: test test_mpi test_empty test_empty_mpi

test: timing_test.exe
	timing_test.exe

test_mpi: timing_test_mpi.exe
	mpirun -np 2 timing_test_mpi.exe

test_empty: timing_empty.exe
	timing_empty.exe

test_empty_mpi: timing_empty_mpi.exe
	mpirun -np 2 timing_empty_mpi.exe

EXE: timing_test.exe

timing.o: ModTiming.o

LIBNOMPI = ${LIBDIR}/libNOMPI.a

${LIBNOMPI}:
	cd ../../NOMPI/src; make LIB

OBJ_TEST = timing_test.o

timing_test.exe: $(OBJ_TEST) ${MY_LIB} ${LIBNOMPI}
	$(LINK.f90) $(Lflag2) -o timing_test.exe $(OBJ_TEST) \
	-L${LIBDIR} -lTIMING -lNOMPI

timing_test_mpi.exe: $(OBJ_TEST) ${MY_LIB}
	$(LINK.f90) $(Lflag1) -o timing_test_mpi.exe $(OBJ_TEST) \
	-L${LIBDIR} -lTIMING $(MPILIB)

OBJ_EMPTY= ../srcEmpty/timing_empty.o timing_test.o

../srcEmpty/timing_empty.o: ../srcEmpty/timing_empty.f90
	cd ../srcEmpty; make timing_empty.o

timing_empty.exe: $(OBJ_EMPTY) ${LIBNOMPI}
	$(LINK.f90) ${Lflag1} -o timing_empty.exe $(OBJ_EMPTY) \
	 -L${LIBDIR} -lNOMPI

timing_empty_mpi.exe: $(OBJ_EMPTY)
	$(LINK.f90) ${Lflag2} -o timing_empty_mpi.exe $(OBJ_EMPTY) $(MPILIB)

distclean: clean

