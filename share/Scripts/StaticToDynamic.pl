#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

my $Help  = ($h or $H or $help);
my $Debug = $D;
my $types = ($t or 'real|logical|integer|double precision');
my $large = ($l or 'nBLK|MaxBlock|MaxImplBLK|MaxImplVar'); 
my $mindim= ($d or 3);

use strict;

&print_help if $Help;

my $module;          # the name of the module
my $cont;            # true if previous line is continued
my $line;            # current line with continuation lines
my %dimension;       # hash contains dimension for variables
my $type;            # the type of the variables in the current declaration
my @allocated;       # list of allocated variables
my $allocate;        # true if allocation was written
my $deallocate;      # true if deallocation was written
my $contains;        # true if 'contains' statement has been found

my $WARNING = "StaticToDynamic.pl WARNING:";
my $ERROR = "StaticToDynamic.pl ERROR:";

while(<>){

    if(/^\s*end module/i){

	my $mod_name = $module;
	$mod_name =~ s/([A-Z])/"_".lc($1)/eg;

	print "
contains
" if not $contains;
	print "
  subroutine init$mod_name
". &print_allocate ."
  end subroutine init$mod_name
" if not $allocate;

	print "
  subroutine clean$mod_name
". &print_deallocate . "
  end subroutine clean$mod_name
" if not $deallocate;

	# Print "end module..." line itself
	print;
	last;
    }

    if($contains){
	# Print the line out
	print; 

	# Add allocate statements into subroutine allocate_* or init_*
	if(/^\s*subroutine\s+(allocate|init)_/i){
	    while(<>){
		# Skip empty lines, comments, use statements and declarations
		last unless /^\s*$|\s*!|^\s*use\s|::/;
		print;
	    }
            print &print_allocate.$_;
	}

	# Add deallocate statements into subroutine deallocate_* or clean_*
	print &print_deallocate.$_
           if $contains and /^\s*subroutine\s+(deallocate|clean)_/i;

	# Nothing else to do
	next;
    }

    # Get the name of the module
    $module = $1 if ?^\s*module\s+(\w+)?i;

    # Edit the logical IsDynamic* variable/parameter to .true.:
    s/\.false\./\.true\./i if /^\s*logical\b/i and /\bIsDynamic\w*\s*=/i;

    # Set $contains variable
    $contains = 1 if /^\s*contains\b/;

    # process lines which are before the contains statement

    # Check line for variable declaration
    if(not $cont){
	# Ignore lines that do not look like a declaration
	if(not /^(\s*($types|common))\b/i){
	    print; next;
	}
        # Store the type of the variables
	$type = $1;

	# Ignore parameter statements
	if(/^\s*$type\s*,\s*parameter/i){
	    print; next;
	}
    }

    # Collect continuation lines together
    if($cont){
	# Append continuation line
	$line .= $_;
    }else{
	# Just set $line
	$line = $_;
    }

    # Check if the current line is a continuation line
    if($line =~ /\&\s*(\!.*)?$/){
	$cont = 1; next;
    }else{
	$cont = 0;
    }

    $_=$line;

    # Remove all common statements
    next if $type =~ /common/i;

    # Do nothing if no large arrays are declared
    if(not /$large/i){
	print; next;
    }

    my @variables;

    # Check for "dimension(size1,size2,size3) :: "
    if(/\bdimension\s*(\([^\(]+\))\s*(::)?\s*/i){
	my $dimension = $1;
	my $variables = $';

	# Check if number of dimensions is larger than $mindim
	my $ndim = ($dimension =~ s/,/,/g) + 1;
	if($ndim < $mindim){
	    print; next;
	}

	# Remove comments and junk from variable list
	$variables =~ s/\!.*//ig;
	$variables =~ s/[\s\&\n]//g;
	@variables = split(/,/,$variables);

	# Set the dimensionality of the variables
	my $variable;
	foreach $variable (@variables){
	    $dimension{$variable} = $dimension;
	}
    }else{
	# Match the beginning of the variable declaration
	/^\s*$type\s*(::)?\s*/;

	# Add a comma to the end of the variable list for easier parsing
	my $variables = $'.',';

	# Remove comments and junk from variable list
	$variables =~ s/\!.*//ig;
	$variables =~ s/[\s\&\n]//g;

	# Extract variables one by one from the list
	while($variables =~ s/^\s*(\w+)\s*(\([^\)]+\))?\s*,//){
	    my $variable = $1;
	    my $dimension= $2;
	    push(@variables,$variable);
	    $dimension{$variable} = $dimension;
	}
    }

    if($Debug){
	print "VARIABLES=";
	my $variable;
	foreach $variable (@variables){
	    print $variable,'[',$dimension{$variable},']';
	}
	print "\n";
    }

    # Write out allocatable or normal declaration
    my $variable;
    foreach $variable (@variables){
	my $dimension = $dimension{$variable};
	if($dimension =~ /$large/i){
	    my $ndim = ($dimension =~ s/[^\(\),]+/:/g);
	    if($ndim >= $mindim){
		print "$type, allocatable :: $variable$dimension\n";
		# store variable for (de)allocation statements
		push(@allocated,$variable);
		next;
	    }
	}
	print "$type :: $variable","$dimension{$variable}\n";
    }
}

warn "$ERROR subroutine allocate_* or init_* was not found in $ARGV !!!\n" 
    unless $allocate;
warn "$WARNING subroutine deallocate_* or clean_* was not found in $ARGV !\n" 
    if $Debug and not $deallocate;

exit 0;

##############################################################################

sub print_allocate{

    $allocate=1;
    my $variable;
    my $result = "

    if(allocated($allocated[0])) RETURN
";
    foreach $variable (@allocated){
	$result .= "    allocate($variable".
	    "$dimension{$variable})\n";
    }

    return $result;
}

##############################################################################

sub print_deallocate{

    $deallocate=1;
    my $variable;
    my $result = "

    if(.not.allocated($allocated[0])) RETURN
";
    foreach $variable (@allocated){
	$result .= "    deallocate($variable)\n";
    }

    return $result;
}

##############################################################################

sub print_help{

#BOP
#!ROUTINE: StaticToDynamic.pl - change static F90 declarations to dynamic
#
#!DESCRIPTION:
# When more physics models are combined into a framework, the total memory
# usage my exceed the memory found on one processor. This problem can
# be avoided if the components dynamicallty allocate their variables
# only on the processors used by that component.
# 
#!REVISION HISTORY:
# 07/17/04 G.Toth - initial version
#                   a few improvements added later
#EOP
    print
#BOC
"Purpose: Transform a module with static variable declarations into a 
         module with dynamic (allocatable) declarations.
         If a logical variable/parameter named IsDynamic* = .false. is found,
         its value is modified to .true.
         The original module must contain a subroutine named 
         'allocate_*' or 'init_*' (for example 'init_mod_advance')
         and it may contain a subroutine named 'deallocate_*' or 'clean_*'.
         These subroutines can be empty in the static version of the module,
         or may contain initialization statements or a conditional 
         write statement about the (de)allocation based the logical IsDynamic*.

Usage:   StaticToDynamic.pl [-h] [-D] [-t=TYPES] [-d=MINDIM] [-l=LARGE] 
               [INFILE] [> OUTFILE]

   -h         - this help message

   -D         - print debugging information

   -t=TYPES   - TYPES contains a | separated list of the types of variables 
                to be transformed. 
                Default value is -t='real|logical|integer|double precision'.

   -d=MINDIM  - MINDIM is the minimum dimensionality of variables to be 
                transformed. Default values is -d=3.

   -l=LARGE   - LARGE is a pattern matching the large variable.
                Default values is -l='nBLK|MaxBlock|MaxImplBLK|MaxImplVar'.

   INFILE     - the file to be transformed. Default is reading from STDIN.

   OUTFILE    - the output file. By default output goes to STDOUT.

Examples:

  Convert ModAdvance_static.f90 to a dynamic ModAdvance.f90:

StaticToDynamic.pl ModAdvance_static.f90 > ModAdvance.f90

  Convert only the real and double precision arrays of at least 4 dimensions 
  matching 'nBLK' or 'MaxBlock' in their declaration:

StaticToDynamic.pl -t='real|double precision' -d=4 -l='nBLK|MaxBlock'"
#EOC
   ,"\n\n";
    exit;
}
