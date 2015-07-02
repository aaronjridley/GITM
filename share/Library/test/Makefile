#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

include ../../../Makefile.def

include ../../../Makefile.conf

MPIRUN = mpirun -np 2

help:
	@echo "You can test the following modules in ../src:"
	@echo
	@echo "   test_axes             - test CON_axes"
	@echo "   test_coord_transform  - test ModTimeConvert"
	@echo "   test_freq             - test ModFreq"
	@echo "   test_geopack          - test CON_geopack"
	@echo "   test_interpolate      - test ModInterpolate"
	@echo "   test_interpolate_amr  - test ModInterpolateAMR"
	@echo "   test_io_unit          - test ModIoUnit"
	@echo "   test_initial_state    - test ModInitialState"
	@echo "   test_line_extract     - test CON_line_extract"
	@echo "   test_linear_advection - test ModLinearAdvection"
	@echo "   test_linear_solver    - test ModLinearSolver"
	@echo "   test_lookup_table     - test ModLookupTable"
	@echo "   test_planet_field     - test CON_planet_field"
	@echo "   test_plot_file        - test ModPlotFile"
	@echo "   test_process_var_name - test ModProcessVarName"
	@echo "   test_random_number    - test ModRandomNumber"
	@echo "   test_ray_trace        - test CON_ray_trace"
	@echo "   test_sort             - test ModSort"
	@echo "   test_time_convert     - test ModTimeConvert"
	@echo "   test_triangulate      - test ModTriangulate"
	@echo "   test_utilities        - test ModUtilities"
	@echo "   test                  - test all the above"
	@echo
	@echo "To test everything with easy to read output try"
	@echo
	@echo "   make -s test"
	@echo
	@echo "The parallel tests can be executed like this"
	@echo
	@echo "   make MPIRUN='mpirun -np 4' test_line_extract"
	@echo
	@echo "The default value for MPIRUN is mpirun -np 2."

test:
	-@make test_axes
	-@make test_coord_transform
	-@make test_freq
	-@make test_geopack
	-@make test_io_unit
	-@make test_interpolate
	-@make test_interpolate_amr
	-@make test_initial_state
	-@make test_line_extract
	-@make test_linear_advection
	-@make test_linear_solver
	-@make test_lookup_table
	-@make test_planet_field
	-@make test_plot_file
	-@make test_random_number
	-@make test_ray_trace
	-@make test_sort
	-@make test_time_convert
	-@make test_triangulate
	-@make test_process_var_name
	-@make test_utilities
	@echo ---------------------------------------------------------
	@ls -l *.diff

distclean: clean
	rm -rf *.tmp *.diff test*.out test_initial_state.ref \
		test_process_var_name.log

#######################################################################
OBJ_TIME_CONVERT = \
	test_time_convert.o \
	../src/ModConst.o \
	../src/ModNumConst.o \
	../src/ModTimeConvert.o

test_time_convert_exe:
	@(cd ../src; make ModTimeConvert.o)
	@make test_time_convert.o
	${LINK.f90} -o test_time_convert.exe ${OBJ_TIME_CONVERT} ${Lflag2}

test_time_convert:
	@echo test_time_convert starting        > test_time_convert.diff
	@make test_time_convert_exe            >> test_time_convert.diff
	@echo ./test_time_convert.exe          >> test_time_convert.diff
	-@(./test_time_convert.exe > test_time_convert.tmp)
	-@(diff test_time_convert.tmp test_time_convert.log \
						> test_time_convert.diff)
	@ls -l test_time_convert.diff

#######################################################################
OBJ_COORD_TRANSFORM = \
	test_coord_transform.o \
	../src/ModCoordTransform.o \
	../src/ModNumConst.o

test_coord_transform_exe:
	@(cd ../src; make ModCoordTransform.o)
	@make test_coord_transform.o
	${LINK.f90} -o test_coord_transform.exe ${OBJ_COORD_TRANSFORM} \
		${Lflag2}

test_coord_transform: 
	@echo test_coord_transform starting     > test_coord_transform.diff
	@make test_coord_transform_exe         >> test_coord_transform.diff
	@echo ./test_coord_transform.exe       >> test_coord_transform.diff
	-@(./test_coord_transform.exe > test_coord_transform.tmp)
	-@(${SCRIPTDIR}/DiffNum.pl -b -a=1e-12 \
		test_coord_transform.tmp test_coord_transform.log \
						> test_coord_transform.diff)
	@ls -l test_coord_transform.diff

#######################################################################
OBJ_FREQ = \
	test_freq.o \
	../src/ModFreq.o

test_freq_exe:
	@(cd ../src; make ModFreq.o)
	@make test_freq.o
	${LINK.f90} -o test_freq.exe ${OBJ_FREQ} ${Lflag2}

test_freq:
	@echo test_freq starting		>  test_freq.diff
	@make test_freq_exe			>> test_freq.diff
	@echo './test_freq.exe'			>> test_freq.diff
	-@(./test_freq.exe > test_freq.tmp)
	-@(diff test_freq.tmp test_freq.log)	>  test_freq.diff
	@ls -l test_freq.diff

#######################################################################
OBJ_LINEAR_SOLVER = \
	test_linear_solver.o \
	../src/blas.o \
	../src/lapack.o \
	../src/ModMpi.o \
	../src/ModMpiInterfaces.o \
	../src/ModBlasLapack.o \
	../src/ModLinearSolver.o

test_linear_solver_exe:
	@(cd ../src; make blas.o lapack.o ModLinearSolver.o)
	@make test_linear_solver.o
	${LINK.f90} -o test_linear_solver.exe ${OBJ_LINEAR_SOLVER} ${Lflag1}

test_linear_solver:
	@echo test_linear_solver starting 	>  test_linear_solver.diff
	@make test_linear_solver_exe		>> test_linear_solver.diff
	@echo ./test_linear_solver.exe		>> test_linear_solver.diff
	-@(./test_linear_solver.exe > test_linear_solver.tmp)
	-@(${SCRIPTDIR}/DiffNum.pl -b -a=1e-9 \
		test_linear_solver.tmp test_linear_solver.log \
						>  test_linear_solver.diff)
	@ls -l test_linear_solver.diff

#######################################################################
OBJ_LINEAR_ADVECTION = \
	test_linear_advection.o \
	../src/ModIoUnit.o \
	../src/ModLinearAdvection.o 

test_linear_advection_exe:
	@(cd ../src; make ModLinearAdvection.o)
	@make test_linear_advection.o
	${LINK.f90} -o test_linear_advection.exe ${OBJ_LINEAR_ADVECTION} \
		${Lflag2}

test_linear_advection:
	@echo test_linear_advection starting	>  test_linear_advection.diff
	@make test_linear_advection_exe		>> test_linear_advection.diff
	@echo ./test_linear_advection.exe	>> test_linear_advection.diff
	-@(./test_linear_advection.exe)
	-@(${SCRIPTDIR}/DiffNum.pl \
		test_linear_advection.tmp \
		test_linear_advection.log	>  test_linear_advection.diff)
	@ls -l test_linear_advection.diff

#######################################################################
OBJ_SORT = \
	test_sort.o \
	../src/ModSort.o

test_sort_exe:
	@(cd ../src; make ModSort.o)
	@make test_sort.o
	${LINK.f90} -o test_sort.exe ${OBJ_SORT} ${Lflag2}

test_sort:
	@echo test_sort starting		>  test_sort.diff
	@make test_sort_exe			>> test_sort.diff
	@echo ./test_sort.exe			>> test_sort.diff
	-@(./test_sort.exe > test_sort.tmp)
	-@(diff test_sort.tmp test_sort.log 	>  test_sort.diff)
	@ls -l test_sort.diff

#######################################################################
OBJ_IO_UNIT = \
	test_io_unit.o \
	../src/ModIoUnit.o

test_io_unit_exe:
	@(cd ../src; make ModIoUnit.o)
	@make test_io_unit.o
	${LINK.f90} -o test_io_unit.exe ${OBJ_IO_UNIT} ${Lflag2}

test_io_unit:
	@echo test_io_unit starting			>  test_io_unit.diff
	@make test_io_unit_exe				>> test_io_unit.diff
	@echo ./test_io_unit.exe			>> test_io_unit.diff
	-@(./test_io_unit.exe > test_io_unit.tmp)
	-@(diff test_io_unit.tmp test_io_unit.log	>  test_io_unit.diff)
	@ls -l test_io_unit.diff

#######################################################################
OBJ_INTERPOLATE = \
	test_interpolate.o \
	../src/ModInterpolate.o

test_interpolate_exe:
	@(cd ../src; make ModInterpolate.o)
	@make test_interpolate.o
	${LINK.f90} -o test_interpolate.exe ${OBJ_INTERPOLATE} ${Lflag2}

test_interpolate:
	@echo test_interpolate starting		>  test_interpolate.diff
	@make test_interpolate_exe		>> test_interpolate.diff
	@echo ./test_interpolate_exe		>> test_interpolate.diff
	-@(./test_interpolate.exe > test_interpolate.tmp)
	-@(diff test_interpolate.tmp test_interpolate.log \
						>  test_interpolate.diff)
	@ls -l test_interpolate.diff
#######################################################################
test_interpolate_amr.o: test_interpolate_amr.f90
	${COMPILE.f90} ${Cflag0} test_interpolate_amr.f90

OBJ_INTERPOLATE_AMR = \
	../src/ModKind.o \
	../src/ModRandomNumber.o \
	../src/ModInterpolateAMR.o \
	test_interpolate_amr.o

test_interpolate_amr_exe:
	@(cd ../src; make ModKind.o; make ModRandomNumber.o; make ModInterpolateAMR.o)
	@make test_interpolate_amr.o
	${LINK.f90} -o test_interpolate_amr.exe ${OBJ_INTERPOLATE_AMR} ${Lflag2}

test_interpolate_amr:
	@echo test_interpolate_amr starting     >  test_interpolate_amr.diff
	@make test_interpolate_amr_exe		>> test_interpolate_amr.diff
	@echo ./test_interpolate_amr.exe	>> test_interpolate_amr.diff
	-@(./test_interpolate_amr.exe > test_interpolate_amr.diff)
	@ls -l test_interpolate_amr.diff
#######################################################################
OBJ_INITIAL_STATE = \
	test_initial_state.o \
	../src/ModIoUnit.o \
	../src/ModMpiInterfaces.o \
	../src/ModUtility.o \
	../src/ModUtility_c.o \
	../src/ModReadParam.o \
	../src/ModHdf5Utils.o \
	../src/ModPlotFile.o \
	../src/ModInitialState.o

test_initial_state_exe:
	@(cd ../src; make ModInitialState.o)
	make test_initial_state.o
	${LINK.f90} -o test_initial_state.exe ${OBJ_INITIAL_STATE} ${Lflag1}

test_initial_state:
	@echo test_initial_state starting	 >  test_initial_state.diff
	@make test_initial_state_exe		 >> test_initial_state.diff
	@echo ${MPIRUN} ./test_initial_state_exe >> test_initial_state.diff
	-@(${MPIRUN} ./test_initial_state.exe > test_initial_state.tmp)
	-@gunzip -c test_initial_state.ref.gz > test_initial_state.ref
	-@(${SCRIPTDIR}/DiffNum.pl \
	    test_initial_state.out test_initial_state.ref \
						 >  test_initial_state.diff)
	@ls -l test_initial_state.diff

#######################################################################
OBJ_RANDOM_NUMBER = \
	../src/ModRandomNumber.o

test_random_number_exe:
	@(cd ../src; make ModRandomNumber.o)
	@make test_random_number.o
	${LINK.f90} -o test_random_number.exe \
		test_random_number.o ${OBJ_RANDOM_NUMBER} ${Lflag2}

test_random_number:
	@echo test_random_number starting	>  test_random_number.diff
	@make test_random_number_exe		>> test_random_number.diff
	@echo ./test_random_number.exe		>> test_random_number.diff
	-@(./test_random_number.exe > test_random_number.tmp)
	-@(diff test_random_number.tmp test_random_number.log \
						>  test_random_number.diff)
	@ls -l test_random_number.diff

#######################################################################
OBJ_RAY_TRACE = \
	../src/ModMpiInterfaces.o \
	../src/ModMpi.o \
	../src/CON_ray_trace.o

test_ray_trace_exe:
	@(cd ../src; make CON_ray_trace.o)
	@make test_ray_trace.o
	${LINK.f90} -o test_ray_trace.exe \
		test_ray_trace.o ${OBJ_RAY_TRACE} ${Lflag1}

test_ray_trace:
	@echo test_ray_trace starting		>  test_ray_trace.diff
	@make test_ray_trace_exe		>> test_ray_trace.diff
	@echo ./test_ray_trace.exe		>> test_ray_trace.diff
	-@(./test_ray_trace.exe > test_ray_trace.tmp)
	-@(${SCRIPTDIR}/DiffNum.pl -t test_ray_trace.tmp test_ray_trace.log \
						>  test_ray_trace.diff)
	@ls -l test_ray_trace.diff

#######################################################################
OBJ_LINE_EXTRACT = \
	test_line_extract.o \
	../src/ModMpiInterfaces.o \
	../src/ModMpi.o \
	../src/ModSort.o \
	../src/CON_line_extract.o

test_line_extract_exe: 
	@(cd ../src; make CON_line_extract.o)
	@make test_line_extract.o
	${LINK.f90} -o test_line_extract.exe ${OBJ_LINE_EXTRACT} ${Lflag1}

test_line_extract:
	@echo test_line_extract starting	>  test_line_extract.diff
	@make test_line_extract_exe		>> test_line_extract.diff
	@echo ${MPIRUN} test_line_extract.exe	>> test_line_extract.diff
	-@(${MPIRUN} test_line_extract.exe > test_line_extract.tmp)
	-@(diff test_line_extract.tmp test_line_extract.log \
						>  test_line_extract.diff)
	@ls -l test_line_extract.diff

#######################################################################
OBJ_GEOPACK = \
	test_geopack.o \
	../src/CON_geopack.o \
	../src/ModCoordTransform.o \
	../src/ModNumConst.o

test_geopack_exe:
	@(cd ../src; make ModCoordTransform.o CON_geopack.o)
	@make test_geopack.o
	${LINK.f90} -o test_geopack.exe ${OBJ_GEOPACK} ${Lflag2}

test_geopack:
	@echo test_geopack starting			>  test_geopack.diff
	@make test_geopack_exe				>> test_geopack.diff
	@echo ./test_geopack.exe			>> test_geopack.diff
	-@(./test_geopack.exe > test_geopack.tmp)
	-@(${SCRIPTDIR}/DiffNum.pl -b -r=1e-8 \
		test_geopack.tmp test_geopack.log	>  test_geopack.diff)
	@ls -l test_geopack.diff

#######################################################################
OBJ_AXES = \
	test_axes.o \
	../src/CON_axes.o \
	../src/CON_geopack.o \
	../src/CON_planet.o \
	../src/ModMpi.o \
	../src/ModMpiInterfaces.o \
	../src/ModCoordTransform.o \
	../src/ModIoUnit.o \
	../src/ModReadParam.o \
	../src/ModTimeConvert.o \
	../src/ModUtility.o \
	../src/ModUtility_c.o \
	../src/ModKind.o \
	../src/ModConst.o \
	../src/ModNumConst.o \
	../src/ModPlanetConst.o 

test_axes_exe:
	@(cd ../src; make CON_axes.o)
	@make test_axes.o
	${LINK.f90} -o test_axes.exe ${OBJ_AXES} ${Lflag1}

test_axes: 
	@echo test_axes starting			>  test_axes.diff
	@make test_axes_exe				>> test_axes.diff
	@echo ./test_axes.exe				>> test_axes.diff
	-@(./test_axes.exe > test_axes.tmp)
	-@(${SCRIPTDIR}/DiffNum.pl -b -a=1e-12 \
		test_axes.tmp test_axes.log 		>  test_axes.diff)
	@ls -l test_axes.diff

#######################################################################
OBJ_PLANET_FIELD = \
	test_planet_field.o \
	../src/CON_planet_field.o \
	../src/CON_planet.o \
	../src/CON_axes.o \
	../src/CON_geopack.o \
	../src/ModMpi.o \
	../src/ModMpiInterfaces.o \
	../src/ModPlanetConst.o \
	../src/ModCoordTransform.o \
	../src/ModNumConst.o \
	../src/ModConst.o \
	../src/ModIoUnit.o \
	../src/ModReadParam.o \
	../src/ModTimeConvert.o \
	../src/ModUtility.o \
	../src/ModUtility_c.o

test_planet_field_exe:
	@(cd ../src; make CON_planet_field.o)
	@make test_planet_field.o
	${LINK.f90} -o test_planet_field.exe ${OBJ_PLANET_FIELD} ${Lflag1}

test_planet_field: 
	@echo test_planet_field starting	>  test_planet_field.diff
	@make test_planet_field_exe		>> test_planet_field.diff
	@echo ./test_planet_field.exe 		>> test_planet_field.diff
	-@(./test_planet_field.exe > test_planet_field.tmp)
	-@(${SCRIPTDIR}/DiffNum.pl -t -r=1e-6 \
		test_planet_field.tmp test_planet_field.log \
						>  test_planet_field.diff)
	@ls -l test_planet_field.diff

#######################################################################

OBJ_PLOT_FILE = \
	test_plot_file.o \
	../src/ModIoUnit.o \
	../src/ModMpi.o \
	../src/ModMpiInterfaces.o \
	../src/ModUtility.o \
	../src/ModUtility_c.o \
	../src/ModHdf5Utils.o \
	../src/ModPlotFile.o

test_plot_file_exe:
	@(cd ../src; make ModPlotFile.o)
	@make test_plot_file.o
	${LINK.f90} -o test_plot_file.exe ${OBJ_PLOT_FILE} ${Lflag1}

test_plot_file:
	@echo test_plot_file starting	>  test_plot_file.diff
	@make test_plot_file_exe	>> test_plot_file.diff
	@echo ./test_plot_file.exe	>> test_plot_file.diff
	-@(./test_plot_file.exe > test_plot_file.tmp)
	-@(diff test_plot_file.tmp test_plot_file.log \
					>  test_plot_file.diff)
	@ls -l test_plot_file.diff

#######################################################################

OBJ_LOOKUP_TABLE = \
	test_lookup_table.o \
	../src/ModIoUnit.o \
	../src/ModMpi.o \
	../src/ModMpiInterfaces.o \
	../src/ModReadParam.o \
	../src/ModInterpolate.o \
	../src/ModUtility.o \
	../src/ModUtility_c.o \
	../src/ModHdf5Utils.o \
	../src/ModPlotFile.o \
	../src/ModLookupTable.o

test_lookup_table_exe:
	@(cd ../src; make ModLookupTable.o)
	@make test_lookup_table.o
	${LINK.f90} -o test_lookup_table.exe ${OBJ_LOOKUP_TABLE} ${Lflag1}

test_lookup_table:
	@echo test_lookup_table starting	>  test_lookup_table.diff
	@make test_lookup_table_exe		>> test_lookup_table.diff
	@echo ${MPIRUN} test_lookup_table.exe   >> test_lookup_table.diff
	-@(${MPIRUN} test_lookup_table.exe > test_lookup_table.tmp)
	-@(diff test_lookup_table.tmp test_lookup_table.log \
						>  test_lookup_table.diff)
	@ls -l test_lookup_table.diff

#######################################################################

OBJ_UTILITIES = \
	test_utilities.o \
	../src/ModIoUnit.o \
	../src/ModMpi.o \
	../src/ModMpiInterfaces.o \
	../src/ModUtility.o \
	../src/ModUtility_c.o

test_utilities_exe:
	@(cd ../src; make ModUtility.o ModUtility_c.o)
	@make test_utilities.o
	${LINK.f90} -o test_utilities.exe ${OBJ_UTILITIES} ${Lflag1}

test_utilities:
	@echo test_utilities starting		>  test_utilities.diff
	@make test_utilities_exe		>> test_utilities.diff
	@echo ./test_utilities.exe		>> test_utilities.diff
	-@(./test_utilities.exe > test_utilities.tmp)
	-@(diff -b test_utilities.tmp test_utilities.log \
						>  test_utilities.diff)
	@ls -l test_utilities.diff
	@rmdir xxx

#######################################################################

OBJ_TRIANGULATE = \
	test_triangulate.o \
	../src/ModTriangulate.o

test_triangulate_exe:
	@(cd ../src; make ModTriangulate.o)
	@make test_triangulate.o
	${LINK.f90} -o test_triangulate.exe ${OBJ_TRIANGULATE} ${Lflag2}

test_triangulate:
	@echo test_triangulate starting		>  test_triangulate.diff
	@make test_triangulate_exe		>> test_triangulate.diff
	@echo ./test_triangulate_exe		>> test_triangulate.diff
	-@(./test_triangulate.exe > test_triangulate.tmp)
	-@(diff test_triangulate.tmp test_triangulate.log \
						>  test_triangulate.diff)
	@ls -l test_triangulate.diff

#######################################################################

OBJ_PROCESS_VAR_NAME = \
	test_process_var_name.o \
	ModSingleFluid.o\
	ModExtraVariables.o\
	ModEquation.o \
	../src/ModIoUnit.o \
	../src/ModMpi.o \
	../src/ModMpiInterfaces.o \
	../src/ModUtility.o \
	../src/ModUtility_c.o \
	../src/ModProcessVarName.o

test_process_var_name:
	@if [ -d ${DIR}/GM/BATSRUS ]; then	\
		make test_process_var_name1 BATSDIR=${DIR}/GM/BATSRUS;	\
	else                                    \
		make test_process_var_name1 BATSDIR=${GMDIR};	\
	fi; 

test_process_var_name1:
	@rm -f test_process_var_name.tmp test_process_var_name.log
	@echo test_process_var_name starting > test_process_var_name.tmp
	@(cd ../src; make LIB                > test_process_var_name.log 2>&1)
	@cp ${BATSDIR}/src/ModExtraVariables.f90 .
	@cp ${BATSDIR}/src/ModSingleFluid.f90 .
	@make ModSingleFluid.o ModExtraVariables.o \
		>> test_process_var_name.log 2>&1
	@for i in ${BATSDIR}/srcEquation/*.f90 ; do\
		cp -f $${i} ModEquation.f90;\
		perl -i -pe 's/call init_mhd_variables/ /' ModEquation.f90;\
		echo working on $${i} >> test_process_var_name.tmp;\
		echo working on $${i} >> test_process_var_name.log;\
		sleep 1;\
		make test_process_var_name_exe;\
		make test_process_var_name_run;\
	done
	@make test_process_var_name_check
	@rm -f ModEquation* ModExtraVariables* ModSingleFluid* 

test_process_var_name_exe:
	@rm -f test_process_var_name.o ModEquation.o
	@make ModEquation.o test_process_var_name.o \
		>> test_process_var_name.log 2>&1
	@${LINK.f90} -o test_process_var_name.exe ${OBJ_PROCESS_VAR_NAME} \
		${Lflag1} >> test_process_var_name.log 2>&1

test_process_var_name_run:
	-@(./test_process_var_name.exe >> test_process_var_name.tmp)
	@echo test complete >> test_process_var_name.tmp

test_process_var_name_check:
	@rm -f test_process_var_name.diff
	@(perl -ne 'print if /error/i' \
		test_process_var_name.log test_process_var_name.tmp \
		> test_process_var_name.diff)
	@ls -l test_process_var_name.diff

#######################################################################
