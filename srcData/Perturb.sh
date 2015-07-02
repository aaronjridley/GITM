#!/bin/sh

rm -f UAM.in
cp ../srcData/UAM.in.Perturb.Start ./UAM.in
./GITM.exe
cd UA
pGITM
mv restartOUT restartOUT.12
mkdir restartOUT
rm -f restartIN
ln -s restartOUT.12 restartIN
cd ..
rm -f UAM.in
cp ../srcData/UAM.in.Perturb.Restart ./UAM.in
./GITM.exe
cd UA
pGITM
cd data
../../../srcPython/gitm_diff_images.py
