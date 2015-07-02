#!/usr/bin/perl
use strict;

my $pattern = $ARGV[0];

if(not $pattern){
    print "
Usage: watch.pl PATTERN >& watch.log &

Watch jobs matching pattern. If any of them start to run, delete the others.
The following jobs are available now:

",`qstat -u \$USER`,"

Use PATTERN to select a subset of these, e.g.

watch.pl SWMF >& watch.log &
";
    exit;
}

my @results;
my $running;
LOOP:{
    @results = `qstat -u \$USER | grep $pattern`;
    my $ids;
    foreach (@results){
	/^(\d+)[^:]+:\d\d ([A-Z]) (\d\d:\d\d)/;
	print "id=$1 status=$2 wait=$3\n";
	$ids .= " $1";
	$running = $1 if $2 eq "R";
    }
    print "-------------------------\n";
    if($running){
	$ids =~ s/ $running//;
	print "qdel $ids\n";
	`qdel $ids`;
	last LOOP;
    }

    sleep 5;
    redo;
}

print `qstat -u \$USER`;

print "Finished watch, job $running is running\n";
