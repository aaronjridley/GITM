#!/usr/bin/perl
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

open(TEST, "ps xw |");

while (<TEST>) {
    if (/SWMF.exe/ or /mpirun/) {
	print "Killing processes : \n";
	print $_;
	split;
	system "kill -9 ".$_[0];
    }
}

close(TEST);


