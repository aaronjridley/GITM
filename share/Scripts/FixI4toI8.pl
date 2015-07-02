#!/usr/bin/perl
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
#BOP
#!ROUTINE: FixI4toI8.pl - replace 4 byte integers with 8 byte integers.
#!DESCRIPTION:
# Most machines use 4 byte integers as record length markers in a 
# binary Fortran file. Some Cray machines, however, use 8 byte integers 
# even for the Fortran record length markers. This script 
# can convert binary files consisting of 4 byte integers only 
# to a file consisting of 8 byte integers.
#
#!REVISION HISTORY:
# 07/03/2001 G. Toth - initial version
#EOP
if($#ARGV != 0){
    print
#BOC
"Purpose: add extra 0-s to 4 byte integers to make them 8 byte integers.

Typical usage: 

   FixI4toI8.pl inputfile > crayfile

This should be done on a machine with the same endianness as the input file, 
but the script should not be run on the Cray itself,
because the Perl interpreter does not interpret long integers correctly 
on the Cray.

The endianness of the resulting output file can be changed
with the Fix8Endian.pl script if necessary."
#EOC
   ,"\n\n";
    exit;
}

# No end of line record
undef $/;

# Define 4 byte zero for 8 byte integers in Cray
$zero=pack('L',0);

# Read the whole file into $_
$_=<>;

# Convert the file into a 4 byte integer array
@in=unpack 'L*', $_;
print STDERR "Number of 4 byte integers in file is ",$#in+1,"\n";

# Check if the length of the first record is reasonable
$first = $in[0];
die("Error: length of first record=$first is not 12\n") if $first!=12;

# Initialize array index
$i=0;

# Loop over the Fortran records
while($i <= $#in){

    # length of record in bytes
    $len=$in[$i];

    # number of 4 byte integers in record
    $n  = $len/4;

    # put 4 byte length at the beginning of record (twice the original!!!)
    push(@out,2*$len);

    # write 8 byte integers (0000xxxx) into @out
    for($j=$i+1;$j<=$i+$n;$j++){
	push(@out,$zero,$in[$j]);
    }

    # put 4 byte length at the end of record
    push(@out,2*$len);

    # next element
    $i=$i+$n+2;
}

$_=pack('L*',@out);

print;
