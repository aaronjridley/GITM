#!/bin/sh

# This should be run from the GITM directory.
#   - need 8 processors to complete this test suite.

export RUNDIR=run.short_test
export MPIRUN=mpirun

# Make the run directory to do the tests in:
echo "-----------------------------------------------"
echo "  Making run directory"
echo "-----------------------------------------------"
rm -rf run ${RUNDIR}
make rundir
mv run ${RUNDIR}

# -----------------------------------------------------------------
# Test base-line UAM.in file:

cd ${RUNDIR}

echo "-----------------------------------------------"
echo "  Baseline test"
echo "-----------------------------------------------"
${MPIRUN} -np 4 ./GITM.exe
cd UA ; pGITM
cd data
../../../srcPython/run_plot_model_results.py -var=3 -alt=350 3DALL_t021221_000500.bin
cd ..
mv data data.baseline ; mkdir data
cd ..

# -----------------------------------------------------------------
# Save the baseline UAM.in file:
mv UAM.in UAM.base

# -----------------------------------------------------------------
# Test basic functionality of the gitm_makerun.py script:

export TEST=test_base
echo "-----------------------------------------------"
echo "  Test : ${TEST}"
echo "-----------------------------------------------"
../srcPython/gitm_makerun.py -input UAM.base -output UAM.${TEST} 20130315.01 20130315.0105
rm -f UAM.in ; ln -s UAM.${TEST} UAM.in

${MPIRUN} -np 4 ./GITM.exe
cd UA ; pGITM
cd data
../../../srcPython/run_plot_model_results.py -var=3 -alt=350 3DALL_t130315_010500.bin
cd ..
mv data data.${TEST} ; mkdir data
cd ..

# -----------------------------------------------------------------
# Test:
#   + IMF download
#   + plotting 2DGEL files (electrodynamics in 2d):
#   - note: If you add a plot type, you need to make sure to include 3dall!

export TEST=test_imf
echo "-----------------------------------------------"
echo "  Test : ${TEST}"
echo "-----------------------------------------------"
../srcPython/gitm_makerun.py -input UAM.base -output UAM.${TEST} 20130315.01 20130315.0105 -imf find -2dgel 300 -3dall 300
rm -f UAM.in ; ln -s UAM.${TEST} UAM.in

${MPIRUN} -np 4 ./GITM.exe
cd UA ; pGITM
cd data
../../../srcPython/run_plot_model_results.py -var=3 -alt=0 2DGEL_t130315_010500.bin
cd ..
mv data data.${TEST} ; mkdir data
cd ..


# -----------------------------------------------------------------
# Test IMF download and plotting 2DGEL files (electrodynamics in 2d):
# Test:
#   + IMF download
#   + plotting 2DGEL files (electrodynamics in 2d):
#   + limit latitude range to 45-90
#   - note: If you add a plot type, you need to make sure to include 3dall!

export TEST=test_grid
echo "-----------------------------------------------"
echo "  Test : ${TEST}"
echo "-----------------------------------------------"
../srcPython/gitm_makerun.py -input UAM.base -output UAM.${TEST} 20130315.01 20130315.0105 -imf find -2dgel 300 -3dall 300 -latrange 45 90 -nlats 4
rm -f UAM.in ; ln -s UAM.${TEST} UAM.in

${MPIRUN} -np 8 ./GITM.exe
cd UA ; pGITM
cd data
../../../srcPython/run_plot_model_results.py -var=7 -alt=0 2DGEL_t130315_010500.bin
cd ..
mv data data.${TEST} ; mkdir data
cd ..

# -----------------------------------------------------------------
# Test IMF download and plotting 2DGEL files (electrodynamics in 2d):
# Test:
#   + IMF download
#   + plotting 2DGEL files (electrodynamics in 2d):
#   + limit latitude range to 45-90
#   + include hemispheric power file
#   - note: If you add a plot type, you need to make sure to include 3dall!

export TEST=test_hpi
echo "-----------------------------------------------"
echo "  Test : ${TEST}"
echo "-----------------------------------------------"
../srcPython/gitm_makerun.py -input UAM.base -output UAM.${TEST} 20130315.01 20130315.0105 -imf find -2dgel 300 -3dall 300 -latrange 45 90 -nlats 4 -hpi UA/DataIn/Aurora/power_2013.txt
rm -f UAM.in ; ln -s UAM.${TEST} UAM.in

${MPIRUN} -np 8 ./GITM.exe
cd UA ; pGITM
cd data
../../../srcPython/run_plot_model_results.py -var=7 -alt=0 2DGEL_t130315_010500.bin
cd ..
mv data data.${TEST} ; mkdir data
cd ..

# -----------------------------------------------------------------
# Test IMF download and plotting 2DGEL files (electrodynamics in 2d):
# Test:
#   + IMF download
#   + plotting 2DGEL files (electrodynamics in 2d):
#   + limit latitude range to 45-90
#   + run ovation auroral model
#   - note: If you add a plot type, you need to make sure to include 3dall!

export TEST=test_newell
echo "-----------------------------------------------"
echo "  Test : ${TEST}"
echo "-----------------------------------------------"
../srcPython/gitm_makerun.py -input UAM.base -output UAM.${TEST} 20130315.01 20130315.0105 -imf find -2dgel 300 -3dall 300 -latrange 45 90 -nlats 4 -newell
rm -f UAM.in ; ln -s UAM.${TEST} UAM.in

${MPIRUN} -np 8 ./GITM.exe
cd UA ; pGITM
cd data
../../../srcPython/run_plot_model_results.py -var=7 -alt=0 2DGEL_t130315_010500.bin
cd ..
mv data data.${TEST} ; mkdir data
cd ..

# -----------------------------------------------------------------
# Test:
#   + IMF download
#   + plotting 2DGEL files (electrodynamics in 2d):
#   + limit latitude range to 45-90
#   + downloading SME data from APL
#   + running FTA model
#   - note: If you add a plot type, you need to make sure to include 3dall!

export TEST=test_fta
echo "-----------------------------------------------"
echo "  Test : ${TEST}"
echo "-----------------------------------------------"
../srcPython/gitm_makerun.py -input UAM.base -output UAM.${TEST} 20130315.01 20130315.0105 -imf find -2dgel 300 -3dall 300 -latrange 45 90 -nlats 4 -sme find -fta
rm -f UAM.in ; ln -s UAM.${TEST} UAM.in

${MPIRUN} -np 8 ./GITM.exe
cd UA ; pGITM
cd data
../../../srcPython/run_plot_model_results.py -var=7 -alt=0 2DGEL_t130315_010500.bin
cd ..
mv data data.${TEST} ; mkdir data
cd ..

# -----------------------------------------------------------------
# Test:
#   + IMF download
#   + plotting 2DGEL files (electrodynamics in 2d)
#   + low-latitude dynamo
#   + setting nlats to 4
#   + 1 hour run
#   + fta model
#   + SME download

export TEST=test_dynamo
echo "-----------------------------------------------"
echo "  Test : ${TEST}"
echo "-----------------------------------------------"
../srcPython/gitm_makerun.py -input UAM.base -output UAM.${TEST} 20130315.01 20130315.02 -imf find -2dgel 1800 -3dall 1800 -dynamo -nlats 4 -fta -sme find
rm -f UAM.in ; ln -s UAM.${TEST} UAM.in 

${MPIRUN} -np 8 ./GITM.exe
cd UA ; pGITM
cd data
../../../srcPython/run_plot_model_results.py -var=3 -alt=0 2DGEL_t130315_020000.bin
../../../srcPython/run_plot_model_results.py -var=39 -alt=250 3DALL_t130315_020000.bin
cd ..
mv data data.${TEST} ; mkdir data
cd ..




