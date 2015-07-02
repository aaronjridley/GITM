#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

include ../../../Makefile.def

include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libMAGNETOGRAM.a

OBJECTS = \
	ModMagnetogram.o

OBJECTS_HARMONICS= \
	harmonics.o \
	ModMagHarmonics.o

harmonics.o: ModMagHarmonics.o

POTENTIAL: potential_field.o ${LIBDIR}/libSHARE.a
	${LINK.f90} -o POTENTIAL.exe potential_field.o \
		-L${LIBDIR} -lSHARE ${Lflag1}

OBJECTS_FDIPS= \
	FDIPS_module.o \
	FDIPS_nohypre.o \
	FDIPS.o

FDIPS.o: FDIPS_hypre.o FDIPS_module.o

FDIPS: ${OBJECTS_FDIPS} ${LIBDIR}/libSHARE.a
	${LINK.f90} -o FDIPS.exe ${OBJECTS_FDIPS} \
		-L${LIBDIR} -lSHARE ${Lflag1}

OBJECTS_FDIPS_HYPRE= \
	FDIPS_module.o \
	FDIPS_hypre.o \
	FDIPS.o

FDIPS_HYPRE: ${OBJECTS_FDIPS_HYPRE} ${LIBDIR}/libSHARE.a
	${LINK.f90} -o FDIPS.exe ${OBJECTS_FDIPS_HYPRE} \
		-L${LIBDIR} -lSHARE ${Lflag1} ${HYPRELIB}

HYPRE: hypre.o ${LIBDIR}/libSHARE.a ${HYPRELIB}/libHYPRE.a
	${LINK.f90} -o HYPRE.exe hypre.o \
		-L${LIBDIR} -lSHARE ${Lflag1} ${HYPRELIB}

HARMONICS:
	make ${HARMONICSEXE}
	@echo ${HARMONICSEXE} has been compiled.
	ln -sf ${HARMONICSEXE} .

HARMONICSEXE = ${BINDIR}/HARMONICS.exe

${HARMONICSEXE}: ${OBJECTS_HARMONICS}
	${LINK.f90} -o ${HARMONICSEXE} ${OBJECTS_HARMONICS} \
		-L${LIBDIR} -lSHARE ${Lflag1}

DIPOLE:
	make DIPOLE.exe

DIPOLE.exe: dipole.o
	${LINK.f90} -o DIPOLE.exe dipole.o \
		-L${LIBDIR} -lSHARE ${Lflag1}

DIPOLE11:
	make DIPOLE11.exe

DIPOLE11.exe: dipole11.o
	${LINK.f90} -o DIPOLE11.exe dipole11.o \
	-L${LIBDIR} -lSHARE ${Lflag1}

LIB:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

distclean: clean
	rm -f *.dat *.out *.log *.diff FDIPS.in HARMONICS.in

### TESTS FOR HARMONICS ###

test: 
	@rm -f *.diff
	-@(make test_harmonics)
	-@(make test_potential)
	-@(make test_fdips)
	-@(make test_fdips_hypre)
	@ls -l *.diff

test_harmonics:
	@echo "test_harmonics_compile..." > test_harmonics.diff
	make HARMONICS
	rm -f HARMONICS.exe
	ln -s ${BINDIR}/HARMONICS.exe HARMONICS.exe
	make DIPOLE11
	./DIPOLE11.exe
	@echo "test_harmonics_run..." >> test_harmonics.diff
	cp HARMONICS.in.orig HARMONICS.in
	./HARMONICS.exe > harmonics.log
	@echo "test_harmonics_check..." >> test_harmonics.diff
	make test_harmonics_check

test_harmonics_check:
	${SCRIPTDIR}/DiffNum.pl -t -r=1e-5 -a=1e-12 \
		harmonics11uniform.dat \
		harmonics.ref > test_harmonics.diff
	@ls -l test_harmonics.diff

###############################################################################

test_potential:
	@echo "test_potential_compile..." > test_potential.diff
	make POTENTIAL
	make DIPOLE11
	./DIPOLE11.exe
	@echo "test_potential_run..." >> test_potential.diff
	./POTENTIAL.exe > potential.log
	@echo "test_potential_check..." >> test_potential.diff
	make test_potential_check

test_potential_check:
	${SCRIPTDIR}/DiffNum.pl -t -r=1e-5 -a=1e-12 \
		potential.log \
		potential.ref > test_potential.diff
	@ls -l test_potential.diff

###############################################################################

test_fdips:
	@echo "test_fdips_compile..." > test_fdips.diff
	make FDIPS
	make DIPOLE11
	./DIPOLE11.exe
	cp -f FDIPS.in.orig FDIPS.in
	perl -i -pe 's/real8/ascii/' FDIPS.in
	@echo "test_fdips_run..." >> test_fdips.diff
	mpirun -np 4 ./FDIPS.exe > fdips.log
	./redistribute.pl fdips_field_np010202.out fdips.out
	@echo "test_fdips_check..." >> test_fdips.diff
	make test_fdips_check

test_fdips_check:
	-${SCRIPTDIR}/DiffNum.pl -t -r=1e-5 -a=1e-12 \
		fdips.log \
		fdips.ref > test_fdips_log.diff
	gunzip -c fdips_ref.out.gz > fdips_ref.out
	-${SCRIPTDIR}/DiffNum.pl -t -r=1e-9 -a=1e-9 \
		fdips.out \
		fdips_ref.out > test_fdips.diff
	@ls -l test_fdips_log.diff test_fdips.diff

###############################################################################

test_fdips_hypre:
	@echo "test_fdips_hypre_compile..." > test_fdips_hypre.diff
	make FDIPS_HYPRE
	make DIPOLE11
	./DIPOLE11.exe
	cp -f FDIPS.in.orig FDIPS.in
	perl -i -pe 's/real8/ascii/' FDIPS.in
	perl -i -pe 's/^BiCGSTAB/GMRES/i;s/^ILU/AMG/' FDIPS.in
	@echo "test_fdips_hypre_run..." >> test_fdips_hypre.diff
	mpirun -np 4 ./FDIPS.exe > fdips_hypre.log
	./redistribute.pl fdips_field_np010202.out fdips.out
	@echo "test_fdips_hypre_check..." >> test_fdips_hypre.diff
	make test_fdips_hypre_check

test_fdips_hypre_check:
	-${SCRIPTDIR}/DiffNum.pl -t -r=1e-5 -a=1e-12 \
		fdips_hypre.log \
		fdips_hypre.ref > test_fdips_hypre_log.diff
	gunzip -c fdips_ref.out.gz > fdips_ref.out
	-${SCRIPTDIR}/DiffNum.pl -t -r=1e-9 -a=1e-9 \
		fdips.out \
		fdips_ref.out > test_fdips_hypre.diff
	@ls -l test_fdips_hypre_log.diff test_fdips_hypre.diff
