#!/usr/bin/perl -i.bak -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

my $Help =($h or $H or $Help or $help);
my $From =($f or "State_VGB");
my $To   =($t or "State_GBV");
my $Order=($o or "23451");
my $Debug=$D;

use strict;
#BOP
#!ROUTINE: ChangeIndexOrder.pl - change the index order of F90 variables
#!DESCRIPTION:
# Efficiency of Fortran code can strongly depend on data locality.
# By changing the index order of variables, data which is used together may
# get closer together in memory. Changing the index order of a variable
# which occurs thousands of times in dozens of files can be a daunting task.
# This script helps automating this process.
#!REVISION HISTORY:
# 12/05/2003 G. Toth - initial version
#EOP
if($Help){
    print 
#BOC
'Purpose: 
   Change index order and rename variable.

Usage:

   ChangeIndexOrder.pl [-h] [-D] [-f=FROM] [-t=TO] [-o=ORDER] [file1 file2..]

-h       help message

-D       debug info

-f=FROM  FROM is the name of the variable name to be modified. 
         Default value is State_VGB.

-t=TO    TO is the new name of the variable.
         Default value is State_GBV.

-o=ORDER ORDER defines the index ordering. Indexes are numbered from 1.
         Default value is 23451.

file1..  The files to be changed. 
         Default is to read from STDIN and write to STDOUT.

Examples:

   ChangeIndexOrder.pl -f=X -t=Y -o=213
X(1,2,3)=0.0
Y(2,1,3)=0.0
Ctrl-D

Limitations:

    This is a tool to minimize the mechanical typing needed 
    to change the index order. The index order in "dimension(...)"
    type declarations are not modified. Nested indirect indexing
    like "State_VGB(a(1+b(2)),2,3,4,5)" is not parsed. One level
    indirect indexing is parsed correctly. Continuation lines too.

    Even if the index order is changed correctly, there are several cases
    where the change results in an incorrect code. These need
    to be fixed by hand. It is best to check all the replacements.'
#EOC
    ,"\n\n";
    exit;
}

die "From and To variable names must be different\n" unless $From ne $To;


my @Order;
@Order = split(//,$Order);
my $nIndex = $#Order+1;

print "nIndex=$nIndex order=",join(',',@Order),"\n" if $Debug;

my $head;
my $var;
my $indexes;
my $tail;

while(<>){
    $_ .= <> while /\&\s*$/; # read continuation lines
    while(/\b$From\b/im){
	if(/\b$From[\s\&]*\(/im){
	    if(($var,$indexes) =
	       /(\b$From[\s\&]*\()((([^\(\)]*\([^)]+\))?[^\(\)]*)+)\)/im){

		print "var=$var indexes=$indexes\n" if $Debug;

		$head = $`; $tail = $'; $indexes .= ',';

		my @index=(' '); # one dummy element for index 0
		                 
		while($indexes =~ /(([^,\(\)]*(\([^\)]+\))?[^,\(\)]*)+),/gm){
		    print "index $1\n" if $Debug;
		    push(@index, $1);
		}
		print "indexes = ",join("   ",@index),"\n" if $Debug;

		if($#index == $nIndex){

		    $indexes = join(',',@index[@Order]);

		    $_ = "$head$var$indexes)$tail";
		}else{
		    print STDERR "Found $#index instead of $nIndex indexes ".
			"in $ARGV in line $_";
		}

	    }else{
		print STDERR "Could not match indexes in $ARGV in line $_";
	    }
	}
	s/\b$From\b/$To/i;
    }
    print;
}

