#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

$Help = $h;
$All = $a;

$Outfile = $o; $Outfile = "RenameList.pl" unless $Outfile;

#BOP
#!ROUTINE: Methods.pl - Collect subroutine, entry, function and module names
#!DESCRIPTION:
# This script can be used to avoid name conflicts between components.
# It collects the subroutine, entry, function and module names from Fortran 
# files and creates a renaming list which has the correct component prefix.
# For example 
# \begin{verbatim}
# subroutine read_param
# \end{verbatim}
# will be renamed to
# \begin{verbatim}
# subroutine IH_read_param
# \end{verbatim}
# for the IH component. This code only generates the renaming rules.
# The rules can be used by Rename.pl to do the actual renaming.
#
# There are some special rules for EE, SC, IH, OH, and GM because the
# BATSRUS code implements all of these components.
#
#!REVISION HISTORY:
# 08/03/2003 G.Toth gtoth@umich.edu - initial version
# 09/22/2005        allow a list of files as an argument
#EOP

if($Help or $help){
    print '
This script can be used to avoid name conflicts between components.
It collects the subroutine, entry, function and module names from Fortran 
files and creates a renaming list which has the correct component prefix.
The output list can be used by Rename.pl to do the actual renaming.

For components EE, SC, IH, OH, IH, and GM method names starting with
MH_ EE_ SC_ IH_ OH_ and GM_ are handled intelligently.

',
#BOC
'Usage:

       Methods.pl [-h] [-a] [-o=OUTFILE] ID [FILE1 FILE2 ...]

-h    This help message

-a    Rename ALL functions and subroutines including contained ones.
      This option must be on if the code is written in a sloppy F90 style,
      ie. modules, subroutines, and functions do not end with 
      "end module NAME", "end subroutine NAME", "end function NAME".

-o=OUTFILE
      Write the renaming list into OUTFILE. The default is RenameList.pl,
      which is the default input file for Rename.pl.

ID
      The 2 character component ID for the prefix (e.g. GM or UA).
      This argument has to be present.

FILE1 Source files to check. If no file is listed, the script checks
      all files with extensions .f, .f90, .f95, .F, .F90, and .F95

Example:

Eliminate name conflicts by renaming subroutines and functions in a 
physics module with ID "PM" and version "VERSION":

     cd PM/VERSION/src
     ../../../share/Scripts/Methods.pl PM
     ../../../share/Scripts/Rename.pl -r *.f90'
#EOC
,"\n\n";
    exit $Error;
}

if($#ARGV >= 0){
    $Comp = uc(shift @ARGV);
    if($Comp !~ /^[A-Z][A-Z]$/){
	print "Methods.pl ERROR: Component ID must be 2 capital letters!\n";
	die "Type Methods.pl -h for information on usage...\n";
    }
}else{
    print "Methods.pl ERROR: Argument component ID is missing!\n";
    die "Type Methods.pl -h for information on usage...\n";
}

if(@ARGV){
    @source = @ARGV;
}else{
    @source = glob("*.f90 *.f95 *.F90 *.F95 *.f *.F");
}

foreach $source (@source){

    $F77 = ($source =~ /\.f$/i);

    open(SRC,$source) or die "Could not open file $source\n";

    print "processing file $source F77=$F77\n" if $Debug;

    my $module;
    my $subroutine;
    my $function;

    while(<SRC>){

	### print "processing line $_\n" if $Debug;

        # end module (modules cannot be nested, name is not required at end)
	$module=0 if $module and /^\s*end\s+module/i; 

        # end subroutine XXX (subroutines may be nested, name is required)
	$subroutine='' if /^\s*end\s+subroutine\s+$subroutine/i; 

	# end function XXX (functions may be nested, name is required)
	$function='' if /^\s*end\s+function\s+$function\b/i;    

	next if ($module or $subroutine or $function) and not ($All or $F77);

        # module XXX
	if(/^\s*module\s+(\w+)/i){
	    $module = $+;
	    $method{lc($module)}=$module;           # Use lower case for key
	    push(@module,$module) if $Debug;
	}

        # [recursive] subroutine XXX
	if(/^(recursive\s+)?\s*subroutine\s+(\w+)/i){              
	    $subroutine = $+;
	    $method{lc($subroutine)}=$subroutine;   # Use lower case for key
	    push(@subroutine,$subroutine) if $Debug;
	}

        # enrtry XXX
	if(/^\s*entry\s+(\w+)/i){              
	    $entry = $+;
	    $method{lc($entry)}=$entry;             # Use lower case for key
	    push(@entry,$entry) if $Debug;
	}

	# [recursive] [type] function XXX
	if(/^\s*
	   (recursive\s+)?             # recursive
	   (
	    (double\s+precision|       # double   precision
	     real(\*\d)?|              # real*8
	     integer|                  # integer
	     logical|                  # logical
	     character(\*\d+)?         # character*12
	     )
	    (\s+|\s*\([^\)]+\)\s*     # (len=ada), (kind=ada), (2)
	     )
	    )?
	   function\s+                      # function
	   (\w+)                            # XXX
	   /ix){
	    $function = $+;
	    # print "Function=$function in line:$_";
	    $method{lc($function)}=$function;
	    push(@function,$function) if $Debug;
	}
    }
}


# Make new names for the found methods
foreach $method (sort values %method){
    $_ = $method;                    # use $_ for easy match and replace

    next if /^$Comp\_/;              # name already starts with prefix

    s/^(IH|MH)_// if $Comp eq "GM";  # IH_ and MH_ prefix removed for GM

    s/^(GM|MH)_// if $Comp eq "IH";  # GM_ and MH_ prefix removed for IH

    s/^(GM|IH|MH)_// if $Comp =~ /EE|SC|OH/; #IH_ GM_ MH_ prefix removed for EE,SC,OH

    $replace = $Comp."_$_";          # Add component prefix

    $replace = "NO$replace" 
	if $method{lc($replace)};    # Add NO prefix if replacement exists

    $replace = substr($replace,0,31) # Limit length of name to 31 characters
	if length($replace) > 31;

    $replace{$method}=$replace;      # Save replacement rule
}

open (OUT,">$Outfile") or die "Could not open $Outfile\n";

print OUT '%newname=(',"\n";

foreach $method  (sort keys %replace){
    print OUT "\t$method\t\t=>    $replace{$method},\n";
}

print OUT ");";

close OUT;

print "Modules:     @module \n"     if $Debug;
print "Subroutines: @subroutine \n" if $Debug;
print "Entries:     @entry \n"      if $Debug;
print "Functions:   @function \n"   if $Debug;

exit 0;
