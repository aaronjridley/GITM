#!/usr/bin/perl
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

# This script can keep connection open on machines that would 
# otherwise close it. Simply start it on the X window as
#
# share/Scripts/KeepConnection.pl

$t = 0;
{
    print "time from start=$t\n";
    sleep 10;
    $t += 10;
    redo;
}
