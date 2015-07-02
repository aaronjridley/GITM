#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
my $Help        = ($h or $help);
my $Verbose     = ($v or $verbose);
my $AbsTol      = ($a or $abs or 1e-30);
my $RelTol      = ($r or $rel or 1e-12);
my $MaxLine     = ($l or $lines or 100);
my $TextIgnore  = ($t or $text);
my $BlankIgnore = ($b or $blank);

my $TextDiff = not $TextIgnore;

my $diff = 'diff';
$diff .= ' -b' if $BlankIgnore;

use strict;

&print_help if $Help or not @ARGV;

my $WARNING = "WARNING in DiffNum.pl";
my $ERROR   = "ERROR in DiffNum.pl";

die "$ERROR: there should be two file arguments!\n" unless $#ARGV == 1;

my $File1 = $ARGV[0];
my $File2 = $ARGV[1];

if(not -e $File1){
    print "$ERROR: $File1 does not exist\n"; die "$ERROR\n";
}
if(not -e $File2){
    print "$ERROR: $File2 does not exist\n"; die "$ERROR\n";
}

if(not -T $File1){
    print "$ERROR: $File1 is not an ASCII file\n"; die "$ERROR\n";
}
if(not -T $File2){
    print "$ERROR: $File2 is not an ASCII file\n"; die "$ERROR\n";
}

if(not open(FILE1, $File1)){
    print "$ERROR: could not open $File1\n"; die "$ERROR\n";
}
if(not open(FILE2, $File2)){
    print "$ERROR: could not open $File2\n"; die "$ERROR\n";
}

# Files for text comparison. Use local directory.
my $Text1 = "_diffnum_file1_";
my $Text2 = "_diffnum_file2_";

if($TextDiff){
    open(TEXT1, ">$Text1") or 
	die "$ERROR: could not open $Text1 for writing\n";
    open(TEXT2, ">$Text2") or 
	die "$ERROR: could not open $Text2 for writing\n";
}

# Pattern for numbers

my $Number = '[\+\-]?\d+(\.\d*)?([deDE][\+\-]?\d+)?';

my $iLine1;
my $iLine2;
my $eof1;
my $eof2;
my $OrigLine1;
my $OrigLine2;
my $Line1;
my $Line2;
my $Found1;
my $Found2;
my $Num1;
my $Num2;
my $nDiff=0;

 COMPARE:{

   SEARCH: {
       $Found1 = 0;      
       if( $Line1 =~ s/($Number)/\#/){
	   $Num1   = $1;
	   $Found1 = 1;
	   last SEARCH;
       }
       print TEXT1 $Line1 if $TextDiff;
       $Line1 = <FILE1> or last SEARCH;
       $Line1 = "" if $Line1 =~ /Mellanox|RDMA devices/;
       $OrigLine1 = $Line1;
       $iLine1++;
       redo SEARCH;
   }

   SEARCH: {
       $Found2 = 0;
       if( $Line2 =~ s/($Number)/\#/){
	   $Num2   = $1;
	   $Found2 = 1;
	   last SEARCH;
       }
       print TEXT2 $Line2 if $TextDiff;
       $Line2 = <FILE2> or last SEARCH;
       $Line2 = "" if $Line2 =~ /Mellanox|RDMA devices/;
       $OrigLine2 = $Line2;
       $iLine2++;
       redo SEARCH;
   }

     print "Num1 = $Num1\n" if $Verbose and $Found1;
     print "Num2 = $Num2\n" if $Verbose and $Found2;

     last unless $Found1 and $Found2;

     # replace Fortran style D exponent with Perl's 'E'.
     $Num1 =~ s/[dD]/E/;;
     $Num2 =~ s/[dD]/E/;;

     my $Diff = abs($Num1 - $Num2);
     my $Avrg = 0.5*(abs($Num1)+abs($Num2));

     if( $Diff > $AbsTol and $Diff > $RelTol*$Avrg ){
	 print "DiffNum.pl -a=$AbsTol -r=$RelTol $File1 $File2\n" 
	     unless $nDiff;
	 $nDiff++;
	 if($MaxLine > 5){
	     print "${iLine1}n${iLine2}\n";
	     print "< $OrigLine1";
	     print "--- abs($Num1 - $Num2) = $Diff\n";
	     print "> $OrigLine2";
	     $MaxLine -= 4;
	 }
     }

     redo COMPARE;
 }

close(FILE1);
close(FILE2);

my $Message;
$Message = "$ERROR: number of significant numerical differences is $nDiff\n"
    if $nDiff;

if($TextDiff){
    close(TEXT1);
    close(TEXT2);
    my @diff = `$diff $Text1 $Text2`;
    unlink $Text1, $Text2;
    if(@diff){
	print @diff[0..$MaxLine-1] if $MaxLine > 1;
	$nDiff = grep(/^\d/, @diff);
	$Message .= "$ERROR: number of differences in the text is $nDiff\n";
    }
}

print "$WARNING: there are extra numbers in $File2\n" 
    if $Found2 and not $Found1;
print "$WARNING: there are extra numbers in $File1\n" 
    if $Found1 and not $Found2;

if($Message){
    print $Message;
    die "$ERROR\n";
}

exit 0;

##############################################################################
sub print_help{

    print "
Purpose:
    Compare two files and show significant numerical differences only. 
    Absolute and relative differences below given tolerances are ignored.
    The text between numbers can be compared with diff after the
    numbers are replaced with \#. Differences in number formats
    are ignored. 

Usage: DiffNum.pl [-a=ABSTOL] [-r=RELTOL] [-l=LINES] [-b|-t] FILE1 FILE2

   -a=ABSTOL     - ignore differences smaller than ABSTOL.
                   Default value is 1e-30

   -r=RELTOL     - ignore relative differences smaller than RELTOL.
                   The relative difference is the 2*|a-b|/(|a|+|b|).
                   Default value is 1e-12

   -l=LINES      - Show at most LINES lines of differences.
                   If -l=1 (or -l without a value) is set 
                   then only the total number of differences will be reported.
                   Default value is 100 lines. 

   -b -blank     - Ignore blanks in the text comparison (use diff -b).
                   Default is that the differences in blanks also count.

   -t -text      - Ignore the text between numbers.
                   Default is that the text is also compared.

Examples:

   Compare File1 and File2 and show numerical differences only with 
   a value exceeding 1e-10 and relative value exceeding 1e-6:

DiffNum.pl -t -a=1e-10 -r=1e-6 File1 File2

   Compare File1 and File2 with the default tolerances, ignore blanks
   in the text differences, and show at most 20 lines of comparison:

DiffNum.pl -b -l=20 File1 File2

";
    exit 0;

}
