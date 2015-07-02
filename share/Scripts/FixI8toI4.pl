#!/usr/bin/perl
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
#BOP
#!ROUTINE: FixI8toI4.pl - replace 8 byte integers with 4 byte integers.
#!DESCRIPTION:
# Most machines use 4 byte integers as record length markers in a 
# binary Fortran file. Some Cray machines, however, use 8 byte integers 
# even for the Fortran record length markers. This script 
# can convert binary files consisting of 8 byte integers only 
# to a file consisting of 4 byte integers.
#
#!REVISION HISTORY:
# 07/03/2001 G. Toth - initial version
#EOP
if($#ARGV != 0){
    print 
#BOC
"Purpose: throw away the leading 0-s from 8 byte integers.

Typical usage: move the input file from the Cray to another machine with
the same endianness (e.g. IRIX, but not PC or DEC). Then run

   FixI8toI4.pl crayfile > outputfile

The script cannot be used on the Cray, because the Perl interpreter
does not interpret long integers correctly on the Cray.

The endianness of the resulting output file can be changed
with the Fix4Endian.pl script if necessary."
#EOC
    ,"\n\n";
    exit;
}

# No end of line record
undef $/;

# Read the whole file into $_
$_=<>;

# Convert the file into a 4 byte integer array
@in=unpack 'L*', $_;

# Check if the length of the first record is reasonable
$first = $in[0];
die("Error: length of first record=$first is not 24\n") if $first!=24;

# Initialize array index
$i=0;

# Loop over the Fortran records
while($i <= $#in){

    # length of record in bytes
    $len=$in[$i];

    # number of 4 byte integers and the extra zeros in the record 
    $n  = $len/4;

    # put 4 byte length at the beginning of record (half the original!!!)
    push(@out,$len/2);

    # write 4 byte integers into @out but skip the extra 0000-s
    for($j=$i+2;$j<=$i+$n;$j+=2){
	push(@out,$in[$j]);
    }

    # put 4 byte length at the end of record
    push(@out,$len/2);

    # next element
    $i=$i+$n+2;
}

$_=pack('L*',@out);

print;
