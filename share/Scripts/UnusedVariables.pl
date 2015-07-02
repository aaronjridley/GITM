#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

my $Help = $h or $help;
my $Remove = $r;
my $Debug = $D;

use strict;

if($Help or not @ARGV){
    print '
This script can be used to remove unused variables from Fortran code.
It uses the output of the NAG Fortran compiler to find the unused variables. 
The -w compilation flag has to be removed from the CFLAG definition in 
Makefile.conf so that the compiler warnings are not suppressed.
It is a good ideat to set zero level for optimization for sake of speed.

Usage:

      UnusedVariables.pl [-h] [-r] FILE1 [FILE2 ...]

-h    This help message

-r    Remove variables. Recompile code to check for errors.

FILE1 FILE2 ... 
      Source files to be checked. Note that the script should be run in 
      the directory where the source files are.

Examples:

See the unused variables for all Fortran files:

     UnusedVariables.pl *.f90 *.f

Remove unused variables from a certain file:

     UnusedVariables.pl -r amr.f90
';
    exit;
}

my $source;
foreach $source (@ARGV){

    next if $source eq "MpiTemplate.f90";

    `touch $source`;
    
    my $object = $source;
    $object =~ s/\.\w+$/.o/;

    my @log = `make $object 2>&1 1>/dev/null`;

    # no error messages
    next unless @log;

    # read source file into an array of lines
    open(FILE, $source);
    my @lines = <FILE>;
    close(FILE);

    foreach (@log){
	my $var;
	my $msg;

	print if $Debug;

	die "error compiling original $source: $_" if /^\s*(Fatal )?Error/;

	next unless s/^Warning: $source, line (\d+): //;
	my $nline = $1;

    	if (/(\w+) (explicitly imported) into/){
	    $var = $1;
	    $msg = $2;
	}elsif(/(Unused) (symbol|local variable) (\w+)/){
	    $msg = $1;
	    $var = $3;
	}elsif(/Local variable (\w+) is (initialised but never used)/){
	    $var = $1;
	    $msg = $2;
	}else{
	    next;
	}

	$msg  = lc($msg);
	my $line = $lines[$nline-1];
	$line =~ /end (module|function|subroutine) (\w+)/ or
	    die "line=$line did not match end module/function/subroutine\n";
	my $method = "$1 $2";

	print "$source at line $nline: $var is $msg in method=$method\n";

	my $i;
	my $ilast;
	for($i = $nline-1; $i>0; $i--){
	    $line = $lines[$i-1];
	    last if $line =~ /$method/i; # stop at beginning of method
	    next if $line =~ /^\s*\!/;   # do not remove commented out code
	    next unless $line =~ /\b$var\b/i; # check if variable is present

	    $ilast = $i;
	}

	$line = $lines[$ilast-1];

	next if $line =~ /IMPLEMENTED/; # do not remove these in ModUser

	print "original line $ilast:$line";

        # remove variable, variable = initial, variable(dimensions)
	$line =~ s/\b$var\b(\s*=[^,\n]+|\([^\)]+\))?//i; 

	$line =~ s/,\s*,/,/;             # remove double comma
	$line =~ s/:\s*,/:/;             # remove comma following :
	$line =~ s/,\s*\n/\n/;           # remove trailing comma
	$line =~ s/^(\s*),\s*/$1/;       # remove leading comma
	$line =~ s/^\s*\&\s*(\!.*)?\n//; # remove line "& !comment"

	# delete empty lines, empty declarations and use statements
	if($line =~ /(^\s*|::\s*|only\s*:\s*)(\!.*)?$/i){
	    $line = "";
	    # remove continuation from previous line
	    $lines[$ilast-2] =~ s/(\s*,\s*)?\&\s*\n/\n/;
	}

	print "modified line $ilast:$line";
	print "\n" if not $line;

	$lines[$ilast-1] = $line if $Remove;
    }

    if($Remove){
	my $orig = $source."_orig_";
	`mv $source $orig` unless -f $orig;
	open(FILE, ">$source");
	print FILE @lines;
	close FILE;

	unlink($object);
	my $log = `make $object 2>&1 1>/dev/null`;
	die "error compiling modified $source\n" if $log =~ /^\s*error/i;

    }
}
