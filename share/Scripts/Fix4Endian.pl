#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
#!QUOTE: \clearpage
#BOP
#!QUOTE: \subsection{Binary Data Manipulation}
#
#!ROUTINE: Fix4Endian.pl - change the byte order of every 4 bytes
#!DESCRIPTION:
# Change the byte order in a binary file which normally consists of 
# single precision reals and 4 byte integers and logicals only. 
# The -only=A:B flag allows to limit the replacement to the bytes between
# the A-th and B-th positions (positions are starting from 1).
# The -except=M:N flag allows to exclude the bytes between 
# the B-th and C-th positions (usually a character string).
# The -i flag means that the changes are made "inplace" in possibly
# multiple files.
#\begin{verbatim}
# Usage:  
#         Fix4Endian.pl [-only=A:B] [-except=C:D] < InFile > Outfile
#         Fix4Endian.pl -i [-only=A:B] [-except=C:D] File1 [File2 ...]
#
# Example: fix endianness of the RCM restart files
#
#    share/Scripts/Fix4Endian.pl -i -except=161:240 run/IM/RCM/restartIN/*
#\end{verbatim}
#!REVISION HISTORY:
# 07/13/2001 G. Toth - initial version
# 08/18/2004 G. Toth - added -i -only=A:B -except=C:D options
#EOP
# No end of line record
undef $/;

$Inplace = $i;

# Extract range
if($only){
    $_ = $only;
    ($Start,$Finish) = /(\d+).(\d+)/;
    die "Incorrect format for -only=$only\n"
        unless 0 < $Start and $Start < $Finish;
}else{
    $Start  =  1; 
    $Finish = -1;
}

# Extract excpetion range if present
if($except){
    $_ = $except;
    ($Min,$Max) = /(\d+).(\d+)/;
    die "Incorrect format for -except=$except\n" 
	unless 0 < $Min and $Min <= $Max;
}else{
    $Min = -1;
    $Max = -1;
}

@ARGV = ("_STDIN_") if not @ARGV;

 FILE: foreach $file (@ARGV){

     # Read the whole file into $_
     if($file eq "_STDIN_"){
	 $_=<STDIN>;
     }else{
	 open(FILE,$file) 
	     or (warn "Could not open file=$file\n" and next FILE);
	 $_=<FILE>;
	 close(FILE);
     }
     $length = length();

     # Initialize variables for this file
     $start  = $Start;
     $finish = $Finish;
     $min    = $Min;
     $max    = $Max;

     # Check variables
     if($start > $length){
	 warn "Starting byte number in -only=$only exceeds length of $file\n";
	 next FILE;
     }

     if($finish < 0){
	$finish = $length;
    }elsif($finish > $length){
	$finish = $length;
	warn "Final byte number in -only=$only exceeds length of $file\n";
    }
    
    if($min > $length){
	warn "Minimum byte number in -except=$except ".
	    "exceeds length of $file\n";
    }
    
    if($max > $length){
	$max = $length;
	warn "Maximum byte number in -except=$except ".
	    "exceeds length of $file\n";
    }
    
    # Reverse every 4 bytes
    for($i = $start-1; $i < $finish; $i += 4){
	if($i+1 >= $min and $i < $max){
	    $i = $max;
	    last if $i >= $finish;
	}
	substr($_,$i,4)=reverse substr($_,$i,4);
    }

    if($Inplace){
	open(FILE,">$file") 
	    or (warn "Could not open $file for writing\n" and next);
        print FILE $_;
    }else{
        print
    }
}

