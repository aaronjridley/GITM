SHELL = /bin/sh

include ../srcMake/Makefile.conf

#
#  Makefile for AB2D code:
#

MODULES = \
	AB_module.o\
	AB_XFER_module.o\
	AB_ERROR_module.o\
	AB_SPH_module.o\
	AB_COMM_module.o\
	AB_ARRAY_module.o

OBJECTS = \
	AB_XFER_array_util.o\
	AB_XFER_1blk_util.o

LIB:    libSphere.a
	@echo libSphere.a has been brought up to date.

libSphere.a:      ${MF} ${MODULES} ${OBJECTS}
	@echo 'Creating libSphere library'
	${AR} libSphere.a ${MODULES} ${OBJECTS}


exposed_objects: AB_module.o AB_XFER_module.o AB_ERROR_module.o \
	         AB_XFER_array_util.o AB_XFER_1blk_util.o AB_SPH_module.o

AB_module.o: AB_module.f90 AB_COMM_module.o AB_ARRAY_module.o \
	AB_ERROR_module.o

AB_XFER_array_util.o:  AB_XFER_array_util.f90 AB_module.f90

AB_XFER_1blk_util.o:  AB_XFER_1blk_util.f90 AB_module.f90

AB_XFER_module.o: AB_XFER_module.f90 AB_COMM_module.o AB_ARRAY_module.o \
	AB_ERROR_module.o AB_module.o

AB_SPH_module.o: AB_SPH_module.f90 AB_module.o AB_COMM_module.o \
	AB_ERROR_module.o

AB_COMM_module.o:  AB_COMM_module.f90 AB_ARRAY_module.o AB_ERROR_module.o

AB_ARRAY_module.o:  AB_ARRAY_module.f90

AB_ERROR_module.o:  AB_ERROR_module.f90

