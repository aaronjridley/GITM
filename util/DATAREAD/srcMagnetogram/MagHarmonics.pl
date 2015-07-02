!/usr/bin/perl 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
use strict;

my $filin; my $CR; 
my $datasource; my $n;

# Asking for filename, CR, data source, max harmonics, and number of processors

print "enter file name: \n";
chomp($filin=<STDIN>);

print "enter CR: \n";
chomp($CR=<STDIN>);

print "enter data source (MDI, WSO, GONG, etc.): \n";
chomp($datasource=<STDIN>);

print "enter mumber of processors: \n";
chomp($n=<STDIN>);

# Definning names of output files. Each output file
# is created with a default name and then overwrriten with the 
# appropriate name. 

my $filedat="CR".$CR."_".$datasource.".dat";
my $filetec="CR".$CR."_".$datasource."_tec.dat";
my $fileheader="CR".$CR."_".$datasource.".H";

# Calling idl to read the fits file. This idl routine
# creates a file contains the fits file's header (filename.H) and a 
# tecplot file of the magnetogram (filename_tec.dat).

print "\n";
print "Converting fits file to ASCII \n";

`cp $filin fitsfile.fits `;

`idl run_fits_to_ascii.pro`;
#my @commands = (".r fits_to_ascii\n","fits_to_ascii\n",
#"exit\n");
#open(IDL, "| idl");
#print IDL @commands;
#close IDL;

`mv fitsfile.H $fileheader`;
`mv fitsfile_tec.dat $filetec`;
`rm -f fitsfile.fits`;

# Run the parallel fortran routine to calculate the harmonic coefficients

print "\n";
print "Calculating harmonic coefficients \n";
my $mpirun="mpirun -np ".$n." MagHarmonics.exe";
#`$mpirun`;

`mv harmonics.dat $filedat`;

`rm -f fitsfile.dat`;
