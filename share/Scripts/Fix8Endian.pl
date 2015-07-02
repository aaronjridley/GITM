#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
#BOP
#!ROUTINE: Fix8Endian.pl - change the byte order of a Fortran file of 8 byte reals
#!DESCRIPTION:
# Change the endianness (byte order) in a Fortran file,
# which contains 8 byte reals and integers only.
# This script only works on a machine which has the same endianness as the
# original file.
#\begin{verbatim}
# Usage:  Fix8Endian.pl InFile > Outfile
#         Fix8Endian.pl < InFile > Outfile
#         Fix8Endian.pl -i File1 File2 Fil3 ...
#\end{verbatim}
#!REVISION HISTORY:
# 07/13/2001 G. Toth - initial version
#EOP
# No end of line record

undef $/;

my $InPlace = $i;
if($InPlace){
    foreach $file (@ARGV){
	open (FILE, $file);
	$_ = <FILE>;
	close(FILE);
	print "processing $file, length=",length($_),"\n";
	&process_file;
	open(FILE, ">$file");
	print FILE $_;
	close(FILE);
    }
}else{
    # Read the whole file into $_
    $_ = <>;
    &process_file;
    print;
}
exit;

#########################################################################
sub process_file{
    # initialize pointer for the string
    $i=0;
    while ( $i < length() ){
	# Get length of record
	$len = unpack('L',substr($_,$i,4));

	# Reverse leading 4 byte length marker
	$lenfixed = reverse(substr($_,$i,4));
	substr($_,$i,4)=$lenfixed;

	# Read length from reversed string if it is unreasonably large
	$len = unpack('L',$lenfixed) if $i+$len > length();

	# Check if length is reasonable
	die "At position $i record length $len is too large?!\n" 
	    if $i+$len > length();

	# Reverse 8 byte reals/integers
	for($j=$i+4; $j<$i+$len; $j+=8){
	    substr($_,$j,8)=reverse(substr($_,$j,8));
	}

	# Reverse trailing 4 byte length marker
	substr($_,$j,4)=$lenfixed;

	# Step to the next record
	$i = $j + 4;
    }
}
