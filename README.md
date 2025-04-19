# GITM
The Global Ionosphere/Thermosphere Model

Here are some V-GITM-isms to get the model up and running from a fresh
clone. I'm sure you'll notice that these could be done in better ways,
but I need to make some time to do that. Here you go:

Clone the repository:
git clone https://github.com/aaronjridley/GITM.git VGITM

Verify that you are on the Venus branch. We've never put in the work
to merge this into master:
git branch 
If necessary, git checkout Venus

Run configure script:
./Config.pl -install -compiler=gfortran10 -venus

You can use difference compilers (ifortmpif90, gfortran,
gfortran10). You can check your gfortran version with gfortran
--version. I use 10.2 which is probably out dated. 

You can try to compile now:
make 

I encountered a rank mismatch error that can be fixed by adding
-fallow-argument-mismatch to your Makefile.conf like below:
COMPILE.f90     = ${CUSTOMPATH_F}gfortran -fallow-argument-mismatch

I ran into another issue with ModEUV.f90 trying to compile and some of
the species indices (iCO_, iAr_) not being recognized. This probably
won't happen to you if you've followed the instructions, but just in
case, this is because you probably prescribed the wrong planet in your
./Config.pl step. This is because Earth doesn't have those species and
consequently aren't being defined in ModPlanet.f90.

Hopefully, it compiled. Now you can begin to set up a run directory:
make rundir

By default, it will use FISM. You must have the correct fism file used
if you want this to get past reading the inputs in your UAM.in
file. There's not a good error message if this isn't correct. Out of
the box, it will use the srcData/UAM.in.Venus as your input file and
correspondingly cp the 2009 fism file
(srcData/FISM/fismdaily_2009.dat) in
run/UA/DataIn/fismdaily.dat. I mention this, because if you want a
different date, you need to grab the correct year yourself. 

You should be able to run:
cd run
mpirun -np 4 GITM.exe
