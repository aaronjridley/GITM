#!/bin/sh

MPI=/usr/local/bin/mpirun

# ----------------------------------------------------
# 1D tests
# ----------------------------------------------------

./Config.pl -g=1,1,50,4
make
rm -rf run test_*.diff

# Plain 1D

rm -rf run_1d
make rundir
mv run run_1d
cp srcData/UAM.in.1d run_1d/UAM.in
cd run_1d
$MPI -np 1 ./GITM.exe
../share/Scripts/DiffNum.pl -b -r=1e-5 UA/data/log0000000?.dat UA/DataIn/log0000000?.1d.dat >& ../test_1d.diff
diff UA/data/run_information.txt UA/DataIn/run_information.txt >> ../test_1d.diff
cd UA ; $MPI -np 1 ./pGITM ; cd ..
cd UA/data ; idl < ../DataIn/idl_input.1d ; cd ../..
cd ..
ls -l test_1d.diff

# Eclipse

rm -rf run_eclipse
make rundir
mv run run_eclipse
cp srcData/UAM.in.eclipse run_eclipse/UAM.in
cd run_eclipse
$MPI -np 1 ./GITM.exe
../share/Scripts/DiffNum.pl -b -r=1e-5 UA/data/log0000000?.dat UA/DataIn/log0000000?.eclipse.dat >& ../test_eclipse.diff
diff UA/data/run_information.txt UA/DataIn/run_information.txt >> ../test_eclipse.diff
cd UA ; $MPI -np 1 pGITM ; cd ..
cd UA/data ; idl < ../DataIn/idl_input.1d ; cd ../..
cd ..
ls -l test_eclipse.diff

# ----------------------------------------------------
# 3D tests
# ----------------------------------------------------

./Config.pl -g=9,9,50,4
make

rm -rf run_3d
make rundir
mv run run_3d
cp srcData/UAM.in.3d run_3d/UAM.in
cd run_3d
$MPI -np 4 ./GITM.exe
../share/Scripts/DiffNum.pl -b -r=1e-5 UA/data/log0000000?.dat UA/DataIn/log0000000?.3d.dat >& ../test_3d.diff
diff UA/data/run_information.txt UA/DataIn/run_information.3d.txt >> ../test_3d.diff
cd UA ; $MPI pGITM ; cd ..
cd UA/data ; idl < ../DataIn/idl_input.3d ; cd ../..
cd ..
ls -l test_3d.diff

