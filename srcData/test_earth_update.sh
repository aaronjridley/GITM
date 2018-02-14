#!/bin/sh

# ----------------------------------------------------
# 1D tests
# ----------------------------------------------------

# Plain

cp run_1d/UA/data/log00000004.dat srcData/log00000004.1d.dat
cp run_1d/UA/data/run_information.txt srcData/run_information.1d.txt
cp run_1d/UA/data/v15.1d.ps srcData
cp run_1d/UA/data/v25.1d.ps srcData

# Eclipse

cp run_eclipse/UA/data/log00000002.dat srcData/log00000002.eclipse.dat
cp run_eclipse/UA/data/run_information.txt srcData/run_information.eclipse.txt
cp run_eclipse/UA/data/v15.1d.ps srcData/v15.eclipse.ps
cp run_eclipse/UA/data/v25.1d.ps srcData/v25.eclipse.ps

# ----------------------------------------------------
# 3D tests
# ----------------------------------------------------

cp run_3d/UA/data/log00000002.dat srcData/log00000002.3d.dat
cp run_3d/UA/data/run_information.txt srcData/run_information.3d.txt
cp run_3d/UA/data/3d.ps srcData

rm srcData/*~
#cvs commit srcData/*.1d.* srcData/*.eclipse.* srcData/*3d.*
