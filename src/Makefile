ABDIR   = ${UADIR}/srcSphereAB
IODIR   = ${DATAREADINDICESDIR}
GLDIR   = ${UADIR}/srcGlow

SEARCH_EXTRA = -I${ABDIR} -I${IODIR} -I${EMPIRICALUADIR} -I${EMPIRICALIEDIR}

include ../srcMake/Makefile.conf

include Makefile.DEPEND
include Makefile.RULES

MODULES = \
	ModSize.o\
	ModTime.o\
	ModConstants.o\
	ModPlanet.o\
	ModGITM.o\
	ModSources.o\
	ModVertical.o\
	ModInputs.o\
	ModSphereInterface.o\
	ModTimingGITM.o\
	ModRates.o\
	ModEUV.o\
	ModElectrodynamics.o\
	ModSatellites.o\
	ModLimiter.o\
	ModTides.o\
	ModNewell.o\
	ModOvationSME.o\
	ModMagTrace.o\
	ModChemistry.o\
	ModGITMImplicit.o\
	ModHwm.o\
	ModUser.o \
	ModRCMR.o \
	apexsh.o

OBJECTS = \
	planet.o\
	satellites.o\
	calc_electron_temperature.o\
	aurora.o\
	get_potential.o\
	calc_efield.o\
	set_inputs.o\
	init_energy_deposition.o\
	get_glow.o\
	init_grid.o\
	time_routines.o\
	vertical_solver.o\
	library.o\
	stretch_grid.o\
	init_mpi.o\
	restart.o\
	initialize.o\
	init_msis.o\
	init_iri.o\
	init_b0.o\
	init_altitude.o\
	finalize.o\
	calc_neutral_friction.o\
	lubksb.o\
	ludcmp.o\
	calc_electrodynamics.o\
	divergence.o\
	calc_pressure.o\
	calc_rates.o\
	calc_timestep.o\
	exchange_messages_sphere.o\
	timing.o\
	advance.o\
	advance_vertical.o\
	advance_horizontal.o\
	set_vertical_bcs.o\
	set_horizontal_bcs.o\
	output.o\
	output_common.o\
	logfile.o\
	calc_conduction.o\
	calc_physics.o\
	calc_avesza.o\
	calc_euv.o\
	calc_chemistry.o\
	calc_ion_v.o\
	chapman_new.o\
	stop_file.o\
	read_inputs.o\
	add_sources.o\
	calc_sources.o\
	apex_routines.o\
	write_output.o\
	gradient.o\
	tides.o\
	read_MHDIMF_Indices_new.o\
	read_NGDC_Indices_new.o\
	read_NOAAHPI_Indices_new.o\
	read_ACE_data.o\
	read_sme.o\
	user.o \
	calc_tec.o \
	get_location.o \
	RCMR_routines.o

OBJECTS_EXE = main.o ${OBJECTS} ${MODULES}

DEPEND:
	@${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS_EXE} 

EXE = ${BINDIR}/GITM.exe

GITM:	DEPEND
	@make ${EXE}
	@make POST
	@echo "${EXE} has been created"

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

STATIC:
	cp -f ModGITM_static.f90    ModGITM.f90
	cp -f ModEUV_static.f90     ModEUV.f90
	cp -f ModSources_static.f90 ModSources.f90

S2D = ${SCRIPTDIR}/StaticToDynamic.pl

DYNAMIC:	
	${S2D} -l=nBlocksMax -d=4 ModGITM_static.f90    > ModGITM.f90
	${S2D} -l=nBlocksMax -d=4 ModEUV_static.f90     > ModEUV.f90
	${S2D} -l=nBlocksMax -d=4 ModSources_static.f90 > ModSources.f90
# Replaces nBlocksMax with nBlocks to make allocation flexible
	perl -pi -e 's/nBlocksMax/nBlocks/ if /allocate/' \
		ModGITM.f90 ModEUV.f90 ModSources.f90

MY_LIB = libGITM.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${LIBPREV} ${MF} ${MODULES} ${OBJECTS}
	cp -f ${LIBPREV} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS} ${MODULES}

distclean:
	rm -f 	ModSize.f90 \
		ModGITM.f90 \
		ModEUV.f90 \
		ModSources.f90 \
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

