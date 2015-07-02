#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
#BOP
#!ROUTINE: CheckDataName.pl
#!DESCRIPTION:
# The SWMF and some of it components are developed at the University of
# Michigan. The developers agreed in a data naming standard. 
# This scipt checks if the file, program, module, subroutine, function
# and variable names conform with the standard.
#
#!REVISION HISTORY:
# 09/06/2005 G. Toth - initial version
#EOP


# Read source code and check against data naming standard.

# Read command line options
my $Debug           = $D; undef $D;
my $Help            = $h; undef $h;
my $NoFileCheck     = $F; undef $F;
my $NoMethodCheck   = $M; undef $F;
my $NoVariableCheck = $V; undef $F;
my $Verbose         = $v; undef $v;

use strict;

# Check command line options
&print_help if $Help or not @ARGV;

if($NoFileCheck and $NoMethodCheck and $NoVariableCheck){
    warn "Specifying switches -F -M -V together means ".
	"that there is nothing to check...\n";
    exit 0;
}

######################
# Global definitions #
######################	
# Valid components:
my $ValidComp = '[A-Z][A-Z]+_';

# Valid name part with small case: calc b0
my $part = "[a-z][a-z0-9]*";

# Valid name part with upper case: Calc B0
my $Part = "[A-Z][a-z0-9]*";

# Valid first name part: i, i1, I1, Var, Var2
my $FirstPart = '([a-z][0-9]*|[A-Z][a-z0-9]+)';

# Valid method name: calc_b0
my $method = "$part(_$part)*";

# Valid subroutine name: read_parameters CON_stop IE_set_parameters
my $ValidMethodName = "($ValidComp)?$method";

# Valid module name: CON_main IE_ModMain ModSize
my $ValidModuleName = "(CON_$method|BATL_$method|($ValidComp)?Mod($Part)+)";

# Valid non-module file names: CON_main.f90 IE_set_param.F90 set_b0.f90
my $ValidMethodFileName = "$ValidMethodName\.[fF]90";

# Valid module file names: 
#     CON_main.f90 IE_ModMain.f90 ModUtil.F90 ModAdvance_static.f90
my $ValidModuleFileName = "$ValidModuleName(_$part)*\.[fF]90";

# Valid scalar variable names: VariableName IE_iGridSize
my $ValidScalarName = "($ValidComp)?$FirstPart($Part)*";

# Valid array variable names: VariableName IE_GridSize_C State_VGB
my $ValidArrayName = "${ValidScalarName}_[A-Z]+";

# Valid named index name: Rho_ x_ AnyName_ GM_
my $ValidNamedIndex = "$FirstPart($Part)*_|[A-Z][A-Z]_";

# Valid fraction name: c1over3
my $ValidFraction = 'c\d+over\d+';

# Valid first name parts depending on variable/function type:
my %ValidPart1 = ('integer'   => '(D?[i-nI-N]|Max|Min|Ijk)\d*',
		  'logical'   => '(Is|Use|Used|Unused|Do|Done)\d*',
		  'character' => '(Name|Type|String)\d*'
		  );

# Valid array index characters representing 1 index
my $ValidArrayIndex1 = '[ABCDEIPQSV]';

# Valid array index characters representing 3 indexes
my $ValidArrayIndex3 = '[CFGNXYZ]';

#########################################
# Definitions to parse the Fortran code #
#########################################

# Simple Fortran types with possible (len=..) and (kind=..) attributes:
my $SimpleType = '(real|integer|logical|character)(\s*\([^\)]+\))?\b';

# Obsolete Fortran types with *NUMBER, e.g. real*8 character*10
my $ObsoleteType = '(real|integer|logical|character)\s*\*\s*\d+';

# Derived Fortran type
my $DerivedType = 'type\s*\(\s*\w+\s*\)';

# Any Fortran Type
my $AnyType = "($SimpleType|$ObsoleteType|$DerivedType)";

# Other global variables
my $File;       # file name
my $TypeFile;   # file type based on name
my $nLine;      # line number
my $Line;       # line
my $Here;       # current location
my $Type;       # variable type
my $Vars;       # list of variables
my $Function;   # function name while parsing a function

foreach (@ARGV){
    $File = $_;
    my $Name = $File;
    # Delete directory path
    $Name =~ s/.*\///;

    unless($NoFileCheck){
	if($Name =~ /^$ValidModuleFileName$/){
	    $TypeFile = "module";
	}elsif($Name =~ /^$ValidMethodFileName$/){
	    $TypeFile = "not module";
	}else{
	    $TypeFile = "wrong";
	    warn "Invalid file name $File\n";
	    print "File names should match one of the following patterns:\n".
		"$ValidMethodFileName\n$ValidModuleFileName\n" if $Verbose;
	}
    }
    if(not (open FILE, "$File")){
	warn "Could not open file $File\n";
	next;
    }

    &check_file;
}
exit 0;

###############################################################################
sub check_file{

    $nLine = 0;
    while(<FILE>){
	$nLine++;                 # count lines
	next if /^$/;             # skip empty lines

	$Here = "at line $nLine in file $File"; # current location

	my $Length = length($_);
	warn "Line is $Length characters long $Here\n" if $Length > 80;

	next if /^\s*\!/;         # skip pure comment lines
	s/\s*[\n\r]+//;           # cut off trailing space and \n\r
	s/^\s+//;                 # remove leading spaces
	s/\s+/ /g;                # remove multiple spaces

	# Hide quotes !!!
	
	s/\s*\!.*//;              # strip comments
	
	$Line .= $_;              # collect continuation lines
	next if 
	    $Line =~ s/\s*\&$/ /; # remove & and read continuation line

	# store function name inside functions
	if($Line =~ /^function\s+(\w+)/i){
	    $Function = $1;
	}elsif($Line =~ /^end\s+function/i){
	    $Function = "";
	}

	if($Line =~ /^(program|subroutine|module|function|$AnyType\s+function)\s/i){
	    &check_methods unless $NoMethodCheck;
	}elsif($Line =~ /^$AnyType/){
	    &check_variables unless $NoVariableCheck;
	}

	$Line = '';               # delete line
    }
    
}
###############################################################################
sub check_methods{

    $_ = $Line;  # Use $_ for easier pattern matching

    # module: 'module ModSize' 'module CON_main'
    if(/^module\s+(\w+)/i){
	my $Name = $1;
	return if $Name =~ /^procedure$/i;        # 'module procedure' is OK
	if($Name !~ /^$ValidModuleName$/){
	    warn "Invalid name for 'module $Name' $Here\n";
	    print "   Valid module names match the following pattern:\n".
		"   $ValidModuleName\n" if $Verbose;
	}
	return;
    }

    # method: 'subroutine calc_flux' 'function heat_source' 'program swmf'
    if(/^(program|subroutine|function)\s+(\w+)/i){
	my $Method = $1;
	my $Name   = $2;
	if($Name !~ /^$ValidMethodName$/){
	    warn "Invalid name for '$Method $Name' $Here\n";
	    print "   Valid $Method names match the following pattern:\n".
		"   $ValidMethodName\n" if $Verbose;
	}
	return;
    }

    # typed function: 'real heat_source' 'character (len=*) sub_string'
    if(/^($AnyType)\s+function\s+(\w+)/i){
	my $TypeOrig = $1;
	my $Name = $+;

	# Simplify and unify $TypeOrig
	my $Type = lc($TypeOrig);

	# Get rid of length or kind specifier
	$Type =~ s/\s*\(.*//;

	# Get rid of obsolete *NUMBER specifier
	if($Type =~ s/\*\d+$//){
	    warn "Obsolete type declaration in '$TypeOrig function $Name' ".
		"$Here\n";
	    print "   Use kind or length attribute\n" if $Verbose;
	}

	if($Name !~ /^$ValidMethodName$/){
	    warn "Invalid name for '$TypeOrig function $Name' $Here\n";
	    print "   Valid function names match the following pattern:\n".
		"   $ValidMethodName\n" if $Verbose;
	    return;
	}

	# Ignore derived type
	return if $Type eq 'type';

	# Check if the first name part matches the type
	my $part1 = $Name;
	$part1 =~ s/^[A-Z]+_//; # remove leading IE_ or CON_
	$part1 =~ s/_.*//;      # remove trailing name parts

	# real or derived type should not match anything else
	if($Type eq 'real'){
	    my $OtherType;
	    foreach $OtherType ('integer', 'logical', 'character'){
		if($part1 =~ /^($ValidPart1{$OtherType})$/i ){
		    warn "Invalid first name part in ".
			"'$TypeOrig function $Name' $Here\n";
		    print "   The first name part '$part1' ".
			"should not match type $OtherType:\n".
			"   $ValidPart1{$OtherType}\n" if $Verbose;
		}
	    }
	}elsif($Type =~ /^character|integer|logical$/){
	    if($part1 !~ /^($ValidPart1{$Type})$/i){
		warn "Invalid first name part in ".
		    "'$TypeOrig function $Name' $Here\n";
		print "   The first name part '$part1' should match:\n".
		    "   $ValidPart1{$Type}\n" if $Verbose;
	    }
	}else{
	    die "ERROR: invalid type = $Type in '$Line' $Here\n";
	}
	return;
    }

}
###############################################################################
sub check_variables{

    $_ = $Line;  # Use $_ for easier pattern matching

    # Variable
    if(not /\s*::\s*/){
	warn "Variable declaration without '::' in '$_' $Here\n";
	print "   Fortran 90 declarations should contain ::\n" if $Verbose;
	return;
    }
    my $TypeOrig = $`;
    my $Vars     = $'; # This is here to make emacs indentation work: ';

    my $Type = lc($TypeOrig);

    if($Type =~ s/\*\d+//){
	warn "Obsolete variable type '$TypeOrig' $Here\n";
	print "   Use kind or length attribute\n" if $Verbose;
    }

    # Check if this is a parameter declaration
    my $Parameter;
    $Parameter = 1 if $Type =~ /\bparameter\b/;

    # Check for dimension(..., ..., ...) in $Type
    my $nDimAll = 0;
    if($Type =~ s/dimension\s*\(([^\(]+)\)//i){
	my $Dimension = $1; 
	$nDimAll = split(',',$Dimension);
	# Get rid of all the (..) in $Vars
	1 while $Vars =~ s/\s*\([^\(]+\)//;
    }else{
	# Replace dimension of the array variables  Var(..., ..., ...)
	# with #X where X is the number of dimension for arrays
	while($Vars =~ s/\s*\(([^\(]+)\)/\#\#/){
	    my $Dimension = $1;
	    my $nDim = split(',',$Dimension);
	    $Vars =~ s/\#\#/\#$nDim/;
	}
    }

    # Get the basic type
    $Type =~ s/^(character|integer|real|logical|type).*/$1/;

    # Remove '= reshape(#N,#M'
    $Vars =~ s/\s*=\s*reshape\([\#\d,]+//i;

    # Split up $Vars
    my @Vars = split(/\s*,\s*/, $Vars);
    my $Var;
    foreach $Var (@Vars){

        $Var =~ s/\s*=.*$//;               # get rid of initialization part

	my $nDim = $nDimAll;               # Set global dimension if any
	$nDim = $1 if $Var =~ s/\#(\d+)//; # Set individual dimension if any

	print "Type = '$Type' Var = '$Var' nDim = '$nDim'\n" if $Debug;

	# Check for a possible named index
	my $NamedIndex = $Parameter and $Type eq 'integer' and $nDim == 0;
	next if $NamedIndex and $Var =~ /^$ValidNamedIndex$/;

	# Check for special real constants like c1over3
	my $RealConstant = $Parameter and $Type eq 'real' and $nDim == 0;
	next if $RealConstant and $Var =~ /^$ValidFraction$/;

        # Check if this is a declaration for the function return value
        next if lc($Var) eq lc($Function); # ignore 

	# Check scalar or array variable name
	if(not $nDim){
	    if($Var !~ /^$ValidScalarName$/){
		warn "Invalid scalar variable name in '$TypeOrig :: $Var' ".
		    "$Here\n";
		if($Verbose){
		    print "   The variable name should match:\n".
			"   $ValidScalarName\n";
		    print "   or a named index should match:\n".
			"   $ValidNamedIndex\n" if $NamedIndex;
		    print "   or a real parameter fraction should match:\n".
			"   $ValidFraction\n" if $RealConstant;
		}
		next;
	    }
	}elsif($Var !~ /^$ValidArrayName$/){
	    warn "Invalid array name '$Var' in '$Line' $Here\n";
	    print "   The array name should match:\n".
		    "   $ValidArrayName\n" if $Verbose;
	    next;
	}

	# Check if the first name part matches the type
	my $part1 = $Var;
	$part1 =~ s/^[A-Z]+_//;     # remove leading IE_ or CON_
	$part1 =~ s/(.)[A-Z].*/$1/; # remove trailing name parts
	$part1 =~ s/_.*//;          # remove array index part

	# real type should not match anything else
	if($Type eq 'type'){
	    # No restriction on derived type
	}elsif($Type eq 'real'){
	    my $OtherType;
	    foreach $OtherType ('integer', 'logical', 'character'){
		if($part1 =~ /^($ValidPart1{$OtherType})$/ ){
		    warn "Invalid first name part in ".
			"variable '$TypeOrig :: $Var' $Here\n";
		    print "   The first name part '$part1' ".
			"should not match type $OtherType:\n".
			"   $ValidPart1{$OtherType}\n" if $Verbose;
		}
	    }
	}elsif($Type =~ /^character|integer|logical$/){
	    if($part1 !~ /^($ValidPart1{$Type})$/){
		warn "Invalid first name part in ".
		    "variable '$TypeOrig :: $Var' $Here\n";
		print "   The first name part '$part1' should match:\n".
		    "   $ValidPart1{$Type}\n" if $Verbose;
	    }
	}else{
	    die "ERROR invalid type $Type in '$Line' $Here\n";
	}

	# Check if the number of dimensions agrees with the index names
	if($nDim){
	    my $Index = $Var;
	    $Index =~ s/.*_//; # Cut everything before index names _XYZ

	    # Try both C = 1 and  C = 3 dimensions
	    $_ = $Index;
	    my $n1 = s/$ValidArrayIndex1//g + 3* s/$ValidArrayIndex3//g;

	    $_ = $Index;
	    my $n3 = 3* s/$ValidArrayIndex3//g + s/$ValidArrayIndex1//g;

	    if(length($_)){
		warn "Invalid array index character '$_' in variable '$Var' ".
		    "$Here\n";
		print "   Valid array index characters are:\n".
		    "   $ValidArrayIndex1 and $ValidArrayIndex3\n" if $Verbose;
	    }elsif($n1 != $nDim and $n3 != $nDim){
		my $Warning = "Array dimension $nDim differs from $n1";
		$Warning .= " or $n3" if $n1 != $n3;
		warn "$Warning ".
		    "implied by array index in variable '$Var' $Here\n";
	    }

	}
    }

}
###############################################################################

sub print_help{

    print
#BOC
"Purpose:
    Check the file names, subroutine, function, module and variable names
    against the data naming standard. The script can check multiple files.

Usage:

    CheckDataName.pl [-h] [-v] [-F] [-M] [-V] FILE1 [FILE2 ...]

  -h       print help message and stop

  -v       print verbose information

  -F       file names are NOT checked

  -M       method and module names are NOT checked

  -V       variable names are NOT checked

FILE1 ...  source code to be checked

Examples:

   Check one file:

share/Scripts/CheckDataName.pl CON/Control/src/swmf.f90

   Check all .f90 and .F90 files and print verbose information:

share/Scripts/CheckDataName.pl -v CON/Control/src/*.f90

   Check file names only:

share/Scripts/CheckDataName.pl -M -V GM/BATSRUS/src*/*.[fF]90"

#EOC
    ,"\n\n";
    exit 0;
}
