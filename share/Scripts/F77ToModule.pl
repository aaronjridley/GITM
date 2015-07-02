#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

my $Prefix = $p;
my $Help   = $h;

use strict;

&print_help if $Help or not $Prefix or $#ARGV != 1;

my $ERROR = "Error in F77ToModule.pl";

my $indent = (" " x 6);


my $FileIn = $ARGV[0];
open(FILE, $FileIn) or die "$ERROR: could not open input file $FileIn\n";
my @text;
@text = <FILE>;
close FILE;

my $FileOut = $ARGV[1];
open(FILE, ">$FileOut") or die "$ERROR: could not open output file $FileOut\n";

# Name of the module is based on the output file name
my $Module = $Prefix.$FileOut; $Module =~ s/\.f//;

my $unit;           # name of last program unit: subroutine/function/block data
my $iEnd;           # index of last END 
my $iBlockData;     # start index of block data
my $nBlockData;     # number of lines in the block data
my $common;         # true when inside a common block declaration

for(my $i=0; $i<=$#text; $i++){
    $_ = $text[$i];
    if(/^\s[\s\w\*\(\)]*(subroutine|function|block +data) *(\w*)/i){

	my $type    = $1;
	my $name    = $2;

	if($type =~ /block/i){
	    # store location of block data unit
	    $iBlockData = $i;

	    # Make a unique name
	    if($name){
		$name = $Prefix.$name;
	    }else{
		$name = $Module."Data";
	    }
	    # put name back into current line
	    $text[$i] = "      BLOCK DATA $name\n";
	}
	$nBlockData = $iEnd - $iBlockData + 1 if $unit =~ /block/i;

	$text[$iEnd] =~ s/\s*\n/ $unit\n/ if $iEnd;
	$unit = "$type $name";

    }elsif(/^\s+end\s*$/i){
	$iEnd = $i;
    }elsif(/^\s+common\b/i){
	$text[$i] =~ s/\/(\w+)\//\/$Prefix$1\//gi;
	$common = 1;
    }elsif(/^\s\s\s\s\s\S/ and $common){
	$text[$i] =~ s/\/(\w+)\//\/$Prefix$1\//gi;
    }else{
	$common = 0;
    }
}

# Fix the last END statement
$text[$iEnd] =~ s/\s*\n/ $unit\n/ if $iEnd;
$nBlockData = $iEnd - $iBlockData + 1 if $unit =~ /block/i;


# Move block data to the beginning
print FILE splice(@text,$iBlockData,$nBlockData), "\n" if $nBlockData;

print FILE "
      module $Module

      contains

",@text,"
      end module $Module
";

close FILE;

exit 0;

###############################################################################

sub print_help{

    print "
Purpose: wrap an F77 file containing external subroutines and functions into 
         with a fixed format F90 module with contained subroutines/functions. 
         Add a prefix to the name of the module, block data, and common blocks
         to avoid name conflicts.

Usage:   F77ToModule.pl -p=PREFIX F77FILE MODFILE

    -p=PREFIX     - The string PREFIX is added to the names of the module, 
                    block data units and common blocks.

    F77FILE       - input F77 file to be processed

    MODFILE       - output F90 file (in fixed format) using a module.
                    The name of this file determines the name of the module.

Example:

    Wrap the 96 Tsyganenko model with a module ModTsyganenko96:

F77ToModule.pl -p=EGM_ t96_01.f ModTsyganenko96.f

";

    exit 0;
}

