ABDIR   = ${UADIR}/srcSphereAB
IODIR   = ${DATAREADINDICESDIR}
GLDIR   = ${UADIR}/srcGlow
SAMIDIR = ${UADIR}/srcSAMI

SEARCH_EXTRA =  -I${ABDIR} \
	 	-I${IODIR} \
		-I${EMPIRICALUADIR} \
		-I${EMPIRICALIEDIR}

include ../Makefile.conf
include ../Makefile.def

include Makefile.DEPEND
include Makefile.RULES

MODULES = \
	apexsh.o\
	ModAeAuroralModel.o\
	ModApex.o\
	ModChemistry.o\
	ModConstants.o\
	ModEUV.o\
	ModElectrodynamics.o\
	ModFtaModel.o\
	ModGITM.o\
	ModGITMImplicit.o\
	ModHwm14.o\
	ModInputs.o\
	ModLimiter.o\
	ModMagTrace.o\
	ModNewell.o\
	ModOrbital.o\
	ModOvationSME.o\
	ModPlanet.o\
	ModRates.o\
	ModRCMR.o \
	ModReadGitm3d.o\
	ModSami.o\
	ModSatellites.o\
	ModSize.o\
	ModSources.o\
	ModSphereInterface.o\
	ModTides.o\
	ModTime.o\
	ModTimingGITM.o\
	ModUser.o \
	ModVertical.o\
	ModHmeModel.o


OBJECTS = \
	add_sources.o\
	advance.o\
	advance_horizontal.o\
	advance_vertical.o\
	apex_more.o\
	apex_routines.o\
	aurora.o\
	calc_avesza.o\
	calc_chemistry.o\
	calc_conduction.o\
	calc_efield.o\
	calc_electrodynamics.o\
	calc_electron_temperature.o\
	calc_euv.o\
	calc_ion_drag.o\
	calc_ion_v.o\
	calc_ir_heating.o\
	calc_neutral_friction.o\
	calc_physics.o\
	calc_pressure.o\
	calc_rates.o\
	calc_sources.o\
	calc_tec.o \
	calc_timestep.o\
	calc_viscosity.o\
	chapman_new.o\
	check_for_nans.o\
	divergence.o\
	exchange_messages_sphere.o\
	finalize.o\
	get_glow.o\
	get_location.o \
	get_potential.o\
	gradient.o\
	initialize.o\
	init_altitude.o\
	init_b0.o\
	init_energy_deposition.o\
	init_grid.o\
	init_iri.o\
	init_mpi.o\
	init_msis.o\
	library.o\
	logfile.o\
	lubksb.o\
	ludcmp.o\
	output.o\
	output_common.o\
	overwrite.o\
	planet.o\
	RCMR_routines.o \
	read_inputs.o\
	read_ACE_data.o\
	read_f107.o\
	read_MHDIMF_Indices_new.o\
	read_NGDC_Indices_new.o\
	read_NOAAHPI_Indices_new.o\
	read_OMNIWEB_Ap_Indices_new.o\
	read_sme.o\
	restart.o\
	satellites.o\
	set_inputs.o\
	set_horizontal_bcs.o\
	set_vertical_bcs.o\
	stop_file.o\
	stretch_grid.o\
	tides.o\
	time_routines.o\
	timing.o\
	user.o \
	vertical_solver.o\
	vertical_solver_ausm.o\
	vertical_solver_rusanov.o\
	write_output.o

SAMI_EXTRAS = \
	ModSamiInterp.o \
	ModCoupSAMI3.o \
	init_mpi_w_sami.o \
	main_w_sami.o 

OBJECTS_EXE  = main.o ${OBJECTS} ${MODULES}
OBJECTS_SAMI = ${OBJECTS} ${MODULES} ${SAMI_EXTRAS}

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS_EXE}

DEPENDSAMI:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS_EXE}

EXE = ${BINDIR}/GITM.exe

SAMIEXE = ${BINDIR}/GITMSAMI.exe

GITM:	DEPEND
	@make ${EXE}
	@make POST
	@echo "${EXE} has been created"

SAMI:	DEPENDSAMI
	@make ${SAMIEXE}
	@make POST
	@echo "${SAMIEXE} has been created"

${SAMIEXE}:	${OBJECTS_SAMI}
		${LINK.f90} -o ${SAMIEXE} ${OBJECTS_SAMI} \
		-L${ABDIR} -lSphere \
		-L${LIBDIR} -lEMPIRICALIE -lEMPIRICALUA -lINDICES -lSHARE \
		-L${GLDIR} -lGLOW ${LBLAS} \
		-L${SAMIDIR}/src -lSAMI3 \
		${Lflag1}

${EXE}:	${OBJECTS_EXE}
	@(if [ ${OS} = "Darwin" ]; then make GITM_for_Darwin; \
	else                            make GITM_general;     fi)

GITM_general:
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE} \
		-L${ABDIR} -lSphere \
		-L${LIBDIR} -lEMPIRICALIE -lEMPIRICALUA -lINDICES -lSHARE \
		-L${GLDIR} -lGLOW ${LBLAS} ${Lflag1}

GITM_for_Darwin:
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ${ABDIR}/libSphere.a; \
		ar -x ${LIBDIR}/libEMPIRICALIE.a; \
		ar -x ${LIBDIR}/libEMPIRICALUA.a; \
		ar -x ${LIBDIR}/libINDICES.a; \
		ar -x ${LIBDIR}/libSHARE.a; \
		ar -x ${GLDIR}/libGLOW.a
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE} Tmp_/*.o ${LBLAS} ${Lflag1}
	rm -rf Tmp_

POST:	PostProcess.o
	${LINK.f90} ${Lflag2} -o PostProcess.exe PostProcess.o

MY_LIB = libGITM.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${LIBPREV} ${MF} ${MODULES} ${OBJECTS}
	cp -f ${LIBPREV} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS} ${MODULES}


clean:
	rm -f *~ core *.o *.mod fort.* a.out *.exe *.a *.so *.protex

distclean: clean
	rm -f 	ModSize.f90 \
		ModPlanet.f90 \
		ModChemistry.f90 \
		ModRates.f90 \
		aurora.f90 \
		calc_chemistry.f90 \
		calc_electron_temperature.f90 \
		calc_ion_density.f90 \
		calc_neutral_density.f90 \
		calc_rates.f90 \
		init_iri.f90 \
		init_msis.f90 \
		output.f90 \
		planet.f90 \
		set_vertical_bcs.f90
