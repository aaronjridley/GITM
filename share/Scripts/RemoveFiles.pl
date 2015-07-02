#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

my $Verbose = $v;
my $Forced  = $f;

use strict;

die "
Purpose: remove so many files that cannot be matched by the shell * or ?

Usage:   RemoveFiles.pl [-v] [-f] PATTERN
-v      - verbose information: all file names to be removed are printed
-f      - forced removal (no question is asked)
PATTERN - a Perl regular expression to match the file name.

Examples: 
   Remove all files starting with blk:
RemoveFiles.pl 'blk.*' 
   Forced removal of files ending with .rst:
RemoveFiles.pl -f '.*\\.rst'

" unless $ARGV[0];
opendir(DIR,'.');
my @files = grep {/^$ARGV[0]$/ and not -d} readdir(DIR);
closedir(DIR);
die "No file matches '$ARGV[0]' !\n" unless @files;
print "Matching files: @files\n" if $Verbose;
if($Forced){
    print "Removing all ",$#files+1," files matching '$ARGV[0]'\n";
}else{
    print "Remove all ",$#files+1," files matching '$ARGV[0]' [y/n] ? ";
    my $reply = <STDIN>;
    die "Cancelled removal !\n" if $reply !~ /^y|yes$/;
}
unlink @files;
