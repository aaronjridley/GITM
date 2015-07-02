#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

#!QUOTE: \clearpage
#BOP
#!QUOTE: \subsection{Change Endianness of BATSRUS Restart Files with ConvertRestart.pl}
#!ROUTINE: ConvertRestart.pl - change the endianness of BATSRUS restart files
#!DESCRIPTION:
# When moving restart files from one machine to another, it may be
# necessary to change the endianness of the binary files.
# This script is specifically written to convert all the binary restart
# files produced by BATSRUS, as well as to copy the ASCII header file.
#
#!REVISION HISTORY:
# 07/03/2001 G. Toth - initial revision
# 03/11/2003 G. Toth - generalized to work on a machine with 
#                      arbitrary endianness
# 03/08/2006 G. Toth - converts new restart file format and added use strict.
# 06/15/2007 G. Toth - add -e option to explicitly specify target endianness.
# 06/18/2007 G. Toth - add in-place conversion
# 10/13/2010 G. Toth - handle index.rst and series of restart files.
#
#EOP

my $Endian  = ($e or $endian);
my $Quiet   = ($q or $quiat);
my $Help    = ($h or $help);

use strict;

&print_help if $Help or $#ARGV == -1;

my $indir =$ARGV[0];
$indir =~ s/\/$//;
my $outdir=($ARGV[1] or $indir."__tmp");

my $headerfile = "restart.H";
my $datafile   = "data.rst";
my $octreefile = "octree.rst";
my $indexfile  = "index.rst";

my $CODE    = "ConvertRestart.pl";
my $ERROR   = "ERROR in $CODE";
my $WARNING = "WARNING in $CODE";

$indir ne $outdir or 
    die "$ERROR: input and output directories must be different!\n";

# Check things
opendir(INDIR, $indir) or die "$ERROR: cannot open input directory $indir!\n";
my @allfiles = readdir INDIR;
close INDIR;

# No end of line in files
undef $/;

# Figure out the endianness of the machine
my $machine = (79 == unpack('V', pack('L',79))) ? 'little' : 'big';
print "$CODE: This is a $machine endian machine.\n" unless $Quiet;

# Figure out the endianness of the restart files 
# by reading the first record length from octree.rst
my @file = grep /$octreefile/, @allfiles;
my $file = "$indir/$file[0]";
open(INFILE, $file) or die "$ERROR: cannot open $file\n"; 
my $rec; read(INFILE, $rec, 4);
close INFILE;
my $from = (12 == (unpack 'V', $rec)) ? 'little' : 'big';

# Set the endianness of the converted restart files
my $to;
if($Endian =~ /^b/){
    $to = 'big';
}elsif($Endian =~ /^l/){
    $to = 'little';
}elsif($Endian =~ /^m/){
    $to = $machine;
}elsif($Endian =~ /^n/){
    $to = ($machine eq 'big') ? 'little' : 'big';
}else{
    $to = ($from eq 'big') ? 'little' : 'big';
}

if($from eq $to){
    print "$CODE: There is no need to transform\n".
	"    the $from endian restart files in $indir.\n" unless $Quiet;
    exit;
}

# Create output directory if necessary
if(not -d $outdir){
    `mkdir -p $outdir`; 
    die "$ERROR: cannot create output directory $outdir!\n" if $?;
}

# Obtain real precision from restart.H 
my @file = grep /$headerfile/, @allfiles;
print "files=@file\n";
my $file = "$indir/$file[0]";
open(INFILE, $file) or die "$ERROR: cannot open $file\n"; 
$_=<INFILE>; 
close INFILE;

# Obtain the real precision
/([48])\s*nByteReal/i or die "$ERROR: could not find nByteReal in $file\n";
my $nByteReal = $1;

# Collect all *.rst* and *restart.H files
my @rstfiles = grep /\.rst|$headerfile/, @allfiles;
print "$CODE: Number of files in $indir is ",1+$#rstfiles,
    " with $nByteReal byte reals.\n" unless $Quiet;

# Integer format for reading
my $intform = ($from eq 'little')? 'V' : 'N';

print "$CODE: Converting $from endian files in $indir\n".
    "                   to $to endian files in $outdir...\n" unless $Quiet;

my $rstfile;
foreach $rstfile (@rstfiles){

    my $file="$indir/$rstfile";
    open(INFILE, $file) or die "$ERROR: cannot open $file!\n";
    my $file = "$outdir/$rstfile";
    open OUTFILE, ">$file" or die "$ERROR: cannot open $file\n";

    while(read(INFILE, $_, 1000000)){

	if(-t $rstfile){
	    # these are ASCII files, nothing to do
	}elsif($nByteReal == 4 or $rstfile =~ /$octreefile/){
	    # single precision files are easy to convert
	    &convert4;
	}elsif($rstfile =~ /$datafile/){
	    # $datafile is a binary direct access file with fixed record length.
	    # Since there are no record markers it is trivial to swap the byte 
	    # order.
	    &convert8;
	}elsif($rstfile =~ /blk/){
	    # double precision blk files need conversion of 4 byte data markers
	    # and 8 byte data
	    &convert_dbl;
	}
	print OUTFILE $_;
    }
    close INFILE;
    close OUTFILE;
}

if($outdir eq $indir."__tmp"){
    # Remove $indir and replace it with $outdir
    print "$CODE: rm -rf $indir; mv $outdir $indir\n" unless $Quiet;
    `rm -rf $indir; mv $outdir $indir`;
}

exit 0;
###############################################################################
sub convert4{
    my $i;
    for($i=0; $i<length(); $i+=4){
	substr($_,$i,4)=reverse substr($_,$i,4);
    }
}

###############################################################################
sub convert8{
    my $i;
    for($i=0; $i<length(); $i+=8){
	substr($_,$i,8)=reverse substr($_,$i,8);
    }
}

###############################################################################
sub convert_dbl{

    # initialize pointer for the string
    my $i=0;

    while ( $i < length() ){
	# Get length of record
	my $len = unpack($intform,substr($_,$i,4));
	
	# Check if length is reasonable
	die "At position $i record length $len is too large?!\n" 
	    if $i+$len > length();

	# Reverse leading 4 byte length marker
	my $lenfixed = reverse(substr($_,$i,4));
	substr($_,$i,4)=$lenfixed;
	
	# Reverse 8 byte reals/integers
	my $j;
	for($j=$i+4; $j<$i+$len; $j+=8){
	    substr($_,$j,8)=reverse(substr($_,$j,8));
	}

	# Reverse trailing 4 byte length marker
	substr($_,$j,4)=$lenfixed;
	
	# Step to the next record
	$i = $j + 4;
    }
}

##############################################################################

sub print_help{

    print 
#BOC
"Purpose: 

  Change/fix the endianness of binary restart files produced by BATSRUS.
  Also copy the restart.H text file for sake of completeness.

Usage:

  ConvertRestart.pl [-h] [-q] [-e=ENDIAN] INPUTDIR [OUTPUTDIR]

  INPUTDIR     Directory containing the restart files to be converted.

  OUTPUTDIR    Directory for the converted restart files. The OUTPUTDIR 
               directory will be created if the code produces output and 
               the directory does not exist yet. The OUTPUT directory
               must be different from the INPUT directory. If the OUTPUT
               directory is not defined, the original directory is overwritten.

  -h -help     Print help message and exit.

  -q -quiet    Suppress all verbose information.

  -e=ENDIAN    Convert to the endianness defined by ENDIAN. 
               Possible values for ENDIAN are:
               'big', 'little', 'machine', 'not-machine' and 'swap'
               Only the first character is significant. The value 'machine'
               means that the file is converted to the endianness of the 
               machine the script is running on, while 'not-machine' converts
               to the opposite endianness of the machine. The value 'swap'
               changes the endianness of the file. If the requested endianness
               is the same as the endianness of the input files the program
               exits with a warning.
               Default is to 'swap' the endianness.

Examples:

  Swap the endianness of files in RESTART/GM to RESTART_conv/GM:

ConvertRestart.pl RESTART/GM RESTART_conv/GM

  Convert the restart files in-place to the endianness of the machine quietly:

ConvertRestart.pl -q -e=m RESTART/GM"
#EOC
    ,"\n\n";
    exit;
}
