#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
# Preplot multiple files

my $Help = $h;
my $Keep = $k;
my $Gzip = $g;
use strict;

if($Help or not @ARGV){
    print "
Purpose:

   Run preplot (and compress) for multiple Tecplot data files.
   If the data file is compressed, use gunzip before preplot.

Usage: 

   Preplot.pl [-h] [-k] [-g] FILE1 [FILE2 FILE3 ...] 

   -h    Print this help message
   -k    Keep original files
   -g    Compress .plt files

Examples:

   Preplot .dat files and keep the originals:

Preplot.pl -k *.dat

   Preplot compressed .dat files and compress results too:

Preplot.pl -g *.dat.gz

";
    exit;
}

# Loop over all files
my $file;
foreach $file (@ARGV){
    my $datfile = $file;

    # Uncompress file if necessary
    `gunzip -c $file > $datfile` if $datfile =~ s/\.dat\.gz$/.dat/;

    # Check extension
    if($datfile !~ /\.dat$/){
	warn "WARNING in Preplot.pl: extension should be .dat or .dat.gz: ".
	    "$file\n";
	next;
    }

    # Run preplot
    `preplot $datfile`;

    # Check if plt file is produced
    my $pltfile = $datfile;
    $pltfile =~ s/.dat/.plt/;
    die "ERROR in Preplot.pl: no $pltfile was produced\n" unless -s $pltfile;

    # Compress .plt file if required
    `gzip $pltfile` if $Gzip;

    # Remove uncompressed data file if any
    unlink $datfile if $file =~ /\.dat\.gz$/;

    # Remove original file
    unlink $file unless $Keep;
}
