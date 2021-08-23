# GITM
This is the home of the Global Ionosphere/Thermosphere Model (GITM).

GITM has been developed in fortran-90. It has been tested with gfortran
on linux and mac osx as well as ifort on NASA's Pleiades computer.

## Dependencies:

1. GITM needs MPI to work.

## Quick Start:

1. git clone https://github.com/aaronjridley/GITM

2. cd GITM

3. ./Config.pl -install -earth -compiler=gfortran

The biggest issue with the above command is that it assumes that you
have the gfortran compiler and things like mpif90 work ok.  If you
don't have gfortran and mpif90, then you need to get these things
for your computer.

In theory, Mars, Venus, Titan, and LV-426 should work.  These are in
various states of completion, so I wouldn't count on them being
perfect. Clearly, Earth is perfect. Ha ha ha.

4. make

5. make rundir

This creates a run directory that has all of the input files.

6. cd run

7. mpirun -np 4 ./GITM.exe

The default UAM.in file has 2 lat blocks and 2 lon blocks with 9 x 9
cells each, so the default resolution is 180 (deg lat) / (2 * 9) = 10 deg
lat, by 360 (deg lon) / (2 * 9) = 20 deg lon.

If you wanted a grid that is (for example) 1 deg (lat) by 5 deg (lon),
you would need 20 blocks (180 cells) in latitude and 8 blocks (72
cells) in longitude.  You need a total of 160 processors for this
simulation. You can change the UAM.in file for these number of cells.
(To get a simple 5 deg x 5 deg resolution, you need 4 blocks in lat
and 8 blocks in longitude, or 32 processors.  If you have fewer
processors, you can change the src/ModSize.f90 code and adjust the
number of cells in each block to compensate. For example, you have 8
processors, so you can adjust ModSize.f90 to have 18 cells in lat and
lon, then ask for 2 (lat) x 4 (lon) blocks.)

8. cd UA

9. pGITM

10. cd data

11. Make some plots:

../../../srcPython/plot_model_results.py -var=3 -alt=120 3DALL_t021221_000500.bin

Then look at the png file that is created.

## Contributing:

1. Please feel free to e-mail the development team to suggest ideas.

2. Please feel free to open an issue.  The development team gets these
issues and will review them.  We can then reach out to you to figure
out how to incorporate them.

3. Please feel free to fork this repository, make changes as you see
fit, and do a pull request.  Your suggested changes will be reviewed
and incorporated if they fit.

## External Codes:

There are a number of external codes that are not developed by the GITM
team.  For example:

1. APEX - the magnetic coordinate system in the code.  Developed at
NCAR. The IGRF code that comes with it is also not developed at UM.

2. Many electrodynamics models (Weimer, Newell's Ovation Prime,
Mitchell's Ovation SME, others in the util/EMPIRICAL/srcIE directory).

3. MSIS and IRI, which are in util/EMPIRICAL/srcUA. MSIS is used as a
lower BC at Earth. IRI is used to initialize the code.

4. The horizonal wind model (HWM) is used as a lower BC at Earth.

