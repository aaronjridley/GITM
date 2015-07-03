include ../Makefile.def
include ../Makefile.conf
#include ./Makefile.DEPEND
#include ./Makefile.RULES

MODULES = Mod_GLOW.o\
	ModCIMPIT.o\
	ModCXPARS.o

OBJECTS = \
	egrid.o\
	ephoto.o\
	etrans.o\
	exsect.o\
	fieldm.o\
	gchem.o\
	geomag.o\
	glow.o\
	maxt.o\
	qback.o\
	rcolum.o\
	rout.o\
	solzen.o\
	ssflux.o\
	vquart.o\
	GL_library.o

EXE = GLOW.exe
SEARCH_EXTRA = -I../src 

MYLIB = libGLOW.a

LIB:	libGLOW.a
	@echo ${MYLIB} has been brought up to date.

libGLOW.a: ${MF} ${MODULES} ${OBJECTS}
	@echo 'Creating GLOW library' ; \
	${AR} ${MYLIB} ${MODULES} ${OBJECTS}

MY_LIB=libUPTOGL.a

LIBADD:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MF} ${MODULES} ${OBJECTS}
	cp -f ${LIBPREV} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS} ${MODULES}

$OBJECTS: $MODULES


#DEPEND:
#	../share/Scripts/depend.pl ./ ../src/ ../share/Library/src/ ${OBJECTS} ${MODULES}

#GLOW: DEPEND
	 @make ${EXE}
#
#
glow: 	${OBJECTS} ${MODULES} 
	${COMPILE.f77} ${Cflag3} -o ${EXE} main.o ${OBJECTS} ${MODULES}

#
#ModGITM.o:	
#	cd ../src  ;make ModGITM.f90



