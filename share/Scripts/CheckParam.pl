#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

# Read a parameter file and verify its correctness based on an XML description.
# Type CheckParam.pl -h for help.

# This subroutine is outside the scope of strict and the "my" variables
# This is a safety feature, so the global variables can not be changed
sub eval_comp{eval("package COMP; $_[0]")}

# Read command line options
my $Debug       = $D; undef $D;
my $Help        = $h; undef $h;
my $HelpXmlParam= $H; undef $H;
my $HelpXml     = $X; undef $X;
my $Interactive = $i; undef $i;
my $Verbose     = $v; undef $v; 
my $XmlFile     = $x; undef $x;
my $NameComp    = $c; undef $c;
my $Components  = $C; undef $C;
my $Precision   = $p; undef $p;
my $GridSize    = $g; undef $g;
my $nProc       = $n; undef $n;
my $StandAlone  = $S; undef $S;

use strict;

# Pattern to match component ID-s
my $ValidComp = 'EE|GM|IE|IH|IM|OH|PC|PS|PT|PW|RB|SC|SP|UA';

# Error string
my $ERROR = 'CheckParam_ERROR:';

# Error status
my $IsError;

# Set default values (needed in the help message too)
my $XmlFileDefault   = 'Param/PARAM.XML';
my $InputFileDefault = 'run/PARAM.in';

# Print help message and exit if -h switch was used
&print_help if $Help;

# Print help for XML description of parameters and exit if -H switch was used
&print_help_xml_param if $HelpXmlParam;

# Print help about XML and exit if -X switch was used
&print_help_xml if $HelpXml;

# Reset filenames to defaults if needed
$XmlFile  = $XmlFileDefault unless $XmlFile;

# Check the correctness of the command line arguments
&check_arguments;

# Initialize for checking input parameter file
my $InputFile;
$InputFile       = $ARGV[0] or 
    $InputFile   = $InputFileDefault;
my $nLine        = 0; # Line number in the current file
my $IncludeLevel = 0; # Each include file increases the include level by 1
my @nLine;            # Array storing line numbers for all include levels
my @InputFile;        # Array storing file names for all include levels
my $FileHandle;       # File handle for the current file

# these global variables are needed for checking parts of strings without
# complaints and for showing possible options if all string parts fail
my $DontShowParamError = 0; 
my $optionvalues;

# Read the parameter definitions
my $commandList;
$commandList = &parse_xml($XmlFile)
    or die "$ERROR could not parse the XML file\n";

# Store values in the COMP package to be used and changed by the XML rules
&init_comp;

# Set the default values defined at the top level of the XML file
foreach my $node (@{$commandList->{content}}){
    &set_value($node) if $node->{type} eq 'e' and $node->{name} eq 'set';
}

# Find commands
my $commandName;        # Current command name (can be an alias)
my $realName;           # The real name of the command (not the alias)
my %realName;           # Hash for real names  $realName{$alias}=$realName
my %command;            # Hash for tree nodes  $command{$realName}=$node
my %commandText;        # Hash for (help) text $commandText{$realName}=$text
my @required;           # Array of required commands @required=($realName,...)

&find_commands($commandList);

# Read and check the parameter file from STDIN or $ARGV
my %defined;            # stores the line number for each command read
my %definedSessionLast; # stores the commands read in last session
my $nSession=1;         # current session number
my $paramName;          # name of the current parameter
my $paramType;          # type of the current parameter
my $paramValue;         # value of the current parameter
my $UserInput;          # set to true between #USERINPUTBEGIN and #USERINPUTEND
my $InsideComp;         # the component for which parameters are read

$InsideComp = $NameComp if $StandAlone; # For stand alone mode start reading
                                        # the component commands immediately.

# Set the file handle for the top level file (or STDIN for interactive mode)
if($Interactive){
    $FileHandle= \*STDIN;
    $InputFile = 'STDIN';
}else{
    $FileHandle="F$IncludeLevel";
    # Change into local directory if necessary
    $InputFile =~ s/(.*)\///;
    if($1){
	chdir $1 or
	    die "$ERROR Could not change directory to $1\n";
	print "chdir $1\n" if $Debug;
    }

    no strict;
    open($FileHandle,$InputFile) or
	die "$ERROR Could not open parameter input file $InputFile!\n";
}

while($_=&read_line){

    # Read command of form #COMMANDNAME
    if(/^\#(\w+)/){
	$commandName=$1;
    }else{
	$commandName='';
	next;
    }

    # Check and store the command
    next unless &check_command($commandName);

    # Get the real name for aliases (e.g. #BODY instead of #MAGNETOSPHERE)
    $realName = $realName{$commandName};

    print "Read command name=$commandName\n" if $Debug;

    # Read the parameters for the command
    &process_elements($command{$realName});

}
# Check if the final session has the required commands defined and
# if the parameters are correct
&check_session;

exit $IsError;

##############################################################################
sub check_arguments{

    # Check command line arguments

    if($XmlFile and not -f $XmlFile){
	die "$ERROR Could not find XML file $XmlFile\n" 
	    if $XmlFile and not -f $XmlFile;
    }

    if($NameComp and $NameComp !~ /^$ValidComp$/){
	die "$ERROR -c=$NameComp is not among valid component ID-s"
	    ." $ValidComp\n";
    }

    if($Components){
	foreach (split ',',$Components){
	    if(not /^$ValidComp$/){
		die "$ERROR -C=$_ is not among valid component ID-s"
		    ." $ValidComp\n";
	    }
	}
    }

    if($Precision and not $Precision =~ /^single|double$/){
	die "$ERROR -p=$Precision is not 'single' or 'double'\n";
    }

    if($GridSize and not $GridSize =~ /^\d+(,\d+)*$/){
	die "$ERROR -g=$GridSize is not"
	    ." a comma separated list of integers\n";
    }

    if(length($nProc) and not $nProc =~ /^[1-9]\d*$/){
	die "$ERROR -n=$nProc is not a positive integer\n";
    }

    if($Interactive and $ARGV[0]){
	die "$ERROR -i should not be combined"
	    ." with the file argument $ARGV[0]\n";
    }
}

##############################################################################
sub init_comp{

    # Set variables in the COMP package for use in the XML files

    # Set COMP::_IsStandAlone according to the -S switch
    $COMP::_IsStandAlone = $StandAlone;

    # Set COMP::_IsFirstSession to true
    $COMP::_IsFirstSession = 1;

    # Set the size of the grid 
    if($GridSize){
	@COMP::_GridSize = split(',',$GridSize);
    }

    # Set the number of processors
    if($nProc){
	$COMP::_nProc = $nProc;
    }

    # Set COMP::_nByteReal for easy check on precision
    $COMP::_nByteReal = $Precision eq 'single' ? 4 : 8;

    # Set the name of the component for checking
    $COMP::_NameComp = $NameComp;

    if($Components){
	# Set COMP::_Components for the XML rules to check components
	$COMP::_Components = $Components;

	# Create a COMP::_Registered hash for registered components
	foreach (split ',',$Components){
	    $COMP::_Registered{$_} = 1;
	}

	# Initialize COMP::_UsedComp hash with the registered components
	%COMP::_UsedComp = %COMP::_Registered;
    }
}
##############################################################################
sub parse_xml{
    # Parse the XML file and return a pointer to the parsed tree
    my $XmlFile=$_[0];

    push @INC, 'share/Scripts/', '../../share/Scripts/';
    require 'XmlRead.pl';

    open(XMLFILE, $XmlFile) or die "$ERROR could not open $XmlFile\n";
    my $tree = &XmlRead( join('',<XMLFILE>) );
    close XMLFILE;

    my $commandList=$tree->[0];
    $commandList=$tree->[1] if $commandList->{type} ne 'e';

    if($commandList->{type} ne 'e' or $commandList->{name} ne 'commandList'){
	die "$ERROR first node should be an element named commandList\n".
	    "but it has type='$commandList->{type}' ".
	    "and content='$commandList->{content}'\n";
    }

    return $commandList;
}

##############################################################################
sub set_value{
    # set values defined by <set name=".." type=".." value="..."> XML tag
    my $node  = $_[0];
    my $name  = &extract($node->{attrib}->{name}, "string");
    my $type  = &extract($node->{attrib}->{type}, "string");
    my $value = &extract($node->{attrib}->{value}, $type);

    &store($name,$value);
}

##############################################################################
sub find_commands{
    # a recursive traverse of the XML tree finds all the commands
    # and stores them in the global %command hash and  @required array
    my $node = $_[0];
    if($node->{type} eq 'e'){
	if($node->{name} eq 'command'){
	    $realName = $node->{attrib}->{name};
	    $realName{$realName}=$realName;
	    $command{$realName} =$node;

	    print "found command $realName\n" if $Debug;
	    
	    if($Verbose){
		# Store the text of the command
		foreach my $element (@{$node->{content}}){
		    $commandText{$realName}=$element->{content} if
			$element->{type} eq 't';
		}
	    }

	    # Store the aliases into %realName
	    foreach $commandName (split(',',$node->{attrib}->{alias})){
		$realName{$commandName}=$realName;
	    }

	    if($node->{attrib}->{required}){
		my $required=&extract($node->{attrib}->{required},"logical");
		push(@required,$realName) if $required;
	    }

	}
	foreach my $element (@{$node->{content}}){
	    &find_commands($element);
	}
    }
}

##############################################################################
sub read_line{

    # Read and return next line from the parameter file(s)

    # Loop until we run out of EOF (end of file), #END and #INCLUDE commands

    # Pre 5.6 versions of Perl do not allow <$FileHandle> with strict
    no strict;
    while((not $_=<$FileHandle>) or /^\#(END|INCLUDE)\b/){

	if($_){
	    # We found an #END or an #INCLUDE command
	    $nLine++;
	    print "At line $nLine in file $InputFile script command $_"
		if $Debug;
	}else{
	    print "EOF at line $nLine in file $InputFile\n" if $Debug;
	}

	if(/INCLUDE/){
	    # Read the parameter of the INCLUDE command
	    $commandName="INCLUDE";
	    $realName=$realName{$commandName};
	    no strict;
	    if($paramValue = <$FileHandle>){
		use strict;
		$nLine++; 
		chop $paramValue;
		$paramValue =~ s/^ +//;   # Remove leading spaces
		$paramValue =~ s/   .*//; # Remove everything after 3 spaces
		$paramValue =~ s/\t.*//;  # Remove everything after a TAB
		$paramValue =~ s/\s+$//;  # Remove trailing spaces

		print "read script parameter = $paramValue\n" if $Debug;
	    }else{
		use strict;
		&print_error(" for command \#$commandName\n".
			     "\tend of file after command");
		return 0 unless &previous_file;
	    }
	    use strict;
	    my $file;
	    # The parameter of the INCLUDE command is the file name itself
	    $file = $paramValue;
	    if(not $file){
		$paramName = "NameIncludeFile";
		$paramType = "string";
		&param_error("is missing!");
		return 1;
	    }
	    # Try to open include file with new file handle
	    my $FileHandleOld = $FileHandle;
	    $FileHandle       = "F".($IncludeLevel+1);
	    no strict;
	    if(-f $file and open($FileHandle, $file)){
		use strict;

		print "Opened include file '$file'\n" if $Debug;

		# Save current file name and line number
		$InputFile[$IncludeLevel] = $InputFile;
		$nLine[$IncludeLevel]     = $nLine;
		
		# Move to next include level
		$IncludeLevel++;
		$InputFile                = $file;
		$nLine                    = 0;
	    }else{
		use strict;
		&print_error(" for command $_".
			     "\tCould not open include file $file");
		$FileHandle=$FileHandleOld;
		return 1;
	    }
	}else{
	    # File either ended or an #END command was found
	    # Close file and return to previous file or return 0
	    return 0 unless &previous_file;
	}
    }
    use strict;

    if($UserInput and /^\#(BEGIN_COMP|END_COMP|RUN|USERINPUTBEGIN)\b/){
	&print_error( " for command $_".
		      "\tthis command cannot occur after ".
		      "#USERINPUTBEGIN at line $UserInput");
	$UserInput = 0;
    }

    # Check for BEGIN_COMP and END_COMP commands
    if(/^\#BEGIN_COMP\b/ and not $StandAlone){
	if($InsideComp){
	    &print_error(" for command $_".
			 "\talready had BEGIN_COMP $InsideComp");
	}else{
	    # Figure out which component is beginning here
	    ($InsideComp) = /BEGIN_COMP ([A-Z][A-Z])/ or
		&print_error(" for command $_".
			     "\tincorrectly formatted BEGIN_COMP command");

	    # Check if the component is registered and ON
	    if($Components){
		if(not $COMP::_Registered{$InsideComp}){
		    &print_error(" for command $_".
				 "\tcomponent $InsideComp is not registered");
		}
		if(not $COMP::_UsedComp{$InsideComp}){
		    &print_error(" for command $_".
				 "\tregistered component $InsideComp is OFF.");
		}
	    }

	    # Return an empty line
	    $_="\n";
	}
    }elsif(/^\#END_COMP\b/ and not $StandAlone){
	# Extract name of the component from #END_COMP ID
	my $Comp;
	($Comp) = /END_COMP ([A-Z][A-Z])/ or
	    &print_error(" for command $_".
			 "\tincorrectly formatted END_COMP command");

	# Check if the component name matches
	if($Comp ne $InsideComp){
	    &print_error(" for command $_".
			 "\tcomponent does not match BEGIN_COMP $InsideComp");
	}else{
	    # END_COMP matches BEGIN_COMP, so we are back to CON params
	    $InsideComp = '';
	    # Return an empty line
	    $_="\n";
	}
    }elsif(/^\#RUN\b/){ # Check if a new session has started
	# Check if the required commands are defined and
	# if the parameters are correct for the session

	print "Session $nSession is complete\n" if $Debug;

	&check_session;
	undef %definedSessionLast;
	$nSession++;
	$COMP::_IsFirstSession=0;
    }elsif(/^\#USERINPUTBEGIN\b/){
	$UserInput = $nLine+1;
    }elsif(/^\#USERINPUTEND\b/){
	if(not $UserInput){
	    &print_error(" for command $_".
			 "\tthere is no matching #USERINPUTBEGIN command");
	}
	$UserInput = 0;
    }
    
    $nLine++;

    # Return the line only for the selected component outside user input
    if($InsideComp eq $NameComp and not $UserInput){
	return $_;
    }else{
	return "\n";
    }
}
##############################################################################
sub previous_file{
    # Select previous file if it exists

    if($IncludeLevel > 0){
	no strict;
	close $FileHandle;
	use strict;

	# Restore previous include level
	$IncludeLevel--;
	$InputFile = $InputFile[$IncludeLevel];
	if($InputFile eq 'STDIN'){
	    $FileHandle = \*STDIN;
	}else{
	    $FileHandle="F$IncludeLevel";
	}
	$nLine    = $nLine[$IncludeLevel];
	return 1;
    }else{
	return 0;
    }
}

##############################################################################
sub check_if{
    my $node = $_[0];
    my $if;

    if($node->{name} =~ /^(if|rule)$/ )
    {
	# The condition is in the compulsary attribute "expr" 
	# for tags <if...> and <rule ...>

	if(not $if=$node->{attrib}->{expr})
	{
	    print "ERROR: attribute expr= is missing from tag $node->{name}".
		"in the XML description of $commandName\n";
	    $IsError = 1;
	    return 1;
	}
    }
    else
    {
	# The condition is in the optional attribute "if" for other tags
	$if=$node->{attrib}->{if};
    }

    # Evaluate the conditional expression
    if($if){
	my $check=&extract($if,"logical");

	# Rules need to be processed if condition is false
	if($node->{name} eq 'rule'){$check = 1-$check};

	print "if=$if check=$check\n" if $Debug;

	return $check;
    }else{
	return 1;
    }
}

##############################################################################
sub process_elements{
    my $node = $_[0];
    foreach my $element (@{$node->{content}}){
	next if $element->{type} eq 't';
	next unless &check_if($element);
	my $name = $element->{"name"};
	if($name eq 'parameter')
	{
	    if(&read_parameter($node,$element)){
		&check_value_format and
		    &check_param_value($element,$paramValue) and
			&store($paramName,$paramValue);
	    }else{
		&param_error("could not be read from file");
	    }
	}
	elsif($name eq 'set')
	{
	    &set_value($element);
	}
	elsif($name eq 'defined')
	{
	    &check_command($element->{attrib}->{name}); # check and store
	}
	elsif($name eq 'if')
	{
	    &process_elements($element);
	}
	elsif($name eq 'rule')
	{
	    # Substitute $\w+ with the value of the variable in package COMP
	    my $content = $element->{content}->[0]->{content};
	    $content =~ s/^\s+//; $content =~ s/\s+$//;
	    $content =~ s/\$(\w+)/'$COMP::'.$1/gee;
	    print "Rule:\t$element->{attrib}->{expr}\n" if $Verbose;
	    &print_error(" for command \#$commandName:\n\t$content");
	    print "Command description\n$commandText{$realName}\n"
		if $Verbose;
	}
	elsif($name eq 'for')
	{
	    my $index=&extract($element->{attrib}->{name},"string");
	    my $from =&extract($element->{attrib}->{from},"integer");
	    my $to   =&extract($element->{attrib}->{to},  "integer");
	    print "FOR index $index= $from to $to\n" if $Debug;
	    foreach my $value ($from..$to){
		# Set the index name in package
                &store($index,$value) if $index;
		&process_elements($element);
	    }
	}
	elsif($name eq 'foreach')
	{
	    my $index =&extract($element->{attrib}->{name},"string");
	    my $values=&extract($element->{attrib}->{values},"string");
	
	    print "FOREACH index=$index values=$values\n" if $Debug;
	    foreach my $value (split(',',$values)){
		# Set the index name in package
		&store($index,$value);
		&process_elements($element);
	    }
	}
    }
}

##############################################################################
sub read_parameter{
    my $command=$_[0];
    my $param=$_[1];

    my $attrib=$param->{attrib};

    $paramName=&extract($attrib->{name},"string");
    $paramType=lc($attrib->{type});
    my $case  =lc($attrib->{case});

    # read the value of the parameter
    print "$paramName=" if $Debug;
    $paramValue = &read_line or return 0;
    chop $paramValue;

    $paramValue =~ s/^\s+//;  # Remove leading spaces
    $paramValue =~ s/   .*//; # Remove everything after 3 spaces
    $paramValue =~ s/\t.*//;  # Remove everything after a TAB
    $paramValue =~ s/\s+$//;  # Remove trailing spaces

    $paramValue = uc($paramValue) if $case eq "upper";
    $paramValue = lc($paramValue) if $case eq "lower";

    print "reading $paramName = $paramValue\n" if $Debug;

    return 1;
}

##############################################################################
sub check_value_format{
    my $bad=0;

    if($paramType =~ /logical/)
    {
	# T F .true. .false.
	$bad = &bad_value unless 
	    $paramValue =~ s/^(f|\.false\.)$/0/i or 
	    $paramValue =~ s/^(t|\.true\.)$/1/i;
    }
    elsif($paramType =~ /integer/)
    {
	# 3 +30 -124
	$bad = &bad_value unless $paramValue =~ /^[\+\-]?\d+$/;
    }
    elsif($paramType =~ /real/)
    {

	# 3 -3 -3. .21 -3.21 -3.21e2 -3.21D+21, 3/6, 3.6/4.2, 5.0D0/3.0D0
	$bad = &bad_value
	    unless $paramValue =~ 
		/^[\+\-]?(\d+(\.\d*)?|\.\d+)([ed][\+\-]?\d+)?
		(\s*\/\s*
		  [\+\-]?(\d+(\.\d*)?|\.\d+)([ed][\+\-]?\d+)?)?$/xi;
    }
    elsif($paramType =~ /string/)
    {
	# strings cannot have wrong format
    }
    else
    {
	print "For command \#$commandName parameter $paramName ".
	    "has unknown type $paramType !?";
	$bad=1;
    }

    return(not $bad);
}

##############################################################################
sub check_param_value{
    # Check value
    my $node =$_[0];
    my $value=$_[1];

    my $content = $node->{content};
    my $attrib  = $node->{attrib};
    my $type    = &extract($attrib->{type},"string");

    # Set the values for the parts
    my $good=1;
    if($type eq 'strings')
    {
	# check if the number of parts is correct
	my @parts  =split(' ',$paramValue);
	my $nPart  =scalar @parts;
	my $min    =&extract($attrib->{min});
	if($nPart < $min){
	    &param_error("has $nPart parts ",
			 "which is less than minimum $min");
	    return 0;
	}
	my $max=&extract($attrib->{max});
	if($max and $max < $nPart){
	    &param_error("has $nPart parts ",
			 "which is greater than maximum $max");
	    return 0;
	}
	my $duplicate=&extract($attrib->{duplicate},"logical");
	if($duplicate and $nPart < $max){	    
	    foreach my $iPart ($nPart..$max-1){
		$parts[$iPart]=$parts[$nPart-1];
	    }
	    $nPart=$max;
	    print "nPart=$nPart parts=",@parts,"\n";
	}
	my $ordered=&extract($attrib->{ordered},"logical");
	my $iPart=0;
	foreach my $element (@$content)
	{
	    next unless $element->{type} eq 'e';
	    if($element->{name} eq 'part'){
		my $name     = &extract($element->{attrib}->{name},"string");
		my $required = &extract($element->{attrib}->{required},
					"logical");
		my $multiple = &extract($element->{attrib}->{multiple},
					"multiple");
		if($ordered){
		    if(not &check_param_value($element,$parts[$iPart])){
			print "\tin part $iPart named $name\n";
			$good=0;
		    }else{
			&store($name,$value);
		    }
		}else{
		    # Unordered name parts. Check one by one.
		    my $any=0;
		    $DontShowParamError=1;
		    my $store;
		    foreach my $jPart (0..$nPart-1){
			my $part=$parts[$jPart];
			if(length($part)>0 and 
			   &check_param_value($element,$part)){
			    $store.="$part ";
			    $parts[$jPart] = "";
			    $any = 1;
			    last unless $multiple;
			}
		    }
		    chop $store;
		    &store($name,$store) if $any;
		    $DontShowParamError=0;
		    if($required and not $any){
			&param_error("has value '$paramValue'",
				     " with no part matching\n",
				     "\t$optionvalues for $name");
			$good=0;
		    }
		}
		$iPart++;
	    }	    
	}
	if(join('',@parts) and not $ordered){
	    &param_error("has value '$paramValue' containing >> ",
			 join(' ',@parts),
			 " <<\n\twhich does not match anything");
	}
	return $good;
    }

    my $good=0;
    if($attrib->{input} eq 'select')
    {
	# Read options
	$optionvalues="";
	my $optionvalue;
	my $shown=0;

      OPTION: foreach my $element (@$content){
	    next OPTION unless $element->{type} eq 'e';
	    next unless &check_if($element);
	    if($element->{name} eq 'option')
	    {
		if(length($element->{attrib}->{value})>0)
		{
		    $optionvalue = 
			&extract($element->{attrib}->{value},$type);
		}
		else
		{
		    $optionvalue = 
			&extract($element->{attrib}->{name},$type);
		}
		# Collect the possible option values
		$optionvalues .= "$optionvalue,";

		if( length($value)>0 and 
		    (($type eq 'string' and $optionvalue =~/\b$value\b/) or
		     ($type eq 'integer' and $optionvalue == $value) or
		     ($type eq 'real'    and abs($optionvalue-$value)<1e-7)
		     ))
		{
		    $good=1;
		    last OPTION;
		}
	    }
	    elsif($element->{name} eq 'optioninput')
	    {
		$good=&check_range($element,$value,$type);
		if(not $good){
		    chop $optionvalues;
		    print "\tand is not found among options $optionvalues\n";
		    $shown = 1;
		}
	    }
	}
	if(not $good and not $shown){
	    chop $optionvalues;
	    if(not length($value)){
		 &param_error("is missing! Possible options are ",
			      $optionvalues);
	     }else{
		 &param_error("has value $value not listed among options ",
			      $optionvalues);
	     }
	}
    }
    else
    {
	$good=&check_range($node,$value,$type);
    }
    return $good;
}

##############################################################################
sub check_range{
    my $node =$_[0];
    my $value=$_[1];
    my $type =$_[2];

    my $min=&extract($node->{attrib}->{min},$type);
    if(length($min)>0){
	if($value < $min){
	    &param_error("has value $value which is less than minimum $min");
	    return 0;
	}
    }
    my $max=&extract($node->{attrib}->{max},$type);
    if(length($max)>0){
	if($value > $max){
	    &param_error("has value $value which is greater ",
			 "than maximum $max");
	    return 0;
	}
    }
    my $length=&extract($node->{attrib}->{length},"integer");
    if(length($length)>0){
	if(length($value) > $length){
	    &param_error("has length ",length($value),
		" which is longer than maximum length $length");
	    return 0;
	}
    }
    return 1;
}

##############################################################################
sub bad_value{

    &param_error("has an incorrectly formatted value '$paramValue'");

    return 1;
}

##############################################################################
sub check_command{
    my $commandName=$_[0];

    my $realName = $realName{$commandName};
    my $command  = $command{$realName};

    # Check if the command is descibed in the XML file or not
    if( not $command ){
	&print_error(" for command \#$commandName\n\tcommand is unknown");
	return 0;
    }

    # Check if the command is currently available
    if( not &check_if($command) ){
	&print_error(" for command \#$commandName\n".
		     "\tcommand is not available because\n".
		     "\tcondition $command->{attrib}->{if} is false");
	return 0;
    }


    # Check if command was defined before in this session
    if($definedSessionLast{$realName} and not $command->{attrib}->{multiple}){
	print "Warning at line $nLine in file $InputFile:\n".
	    "\tcommand \#$commandName has already been defined ",
	    "in this session:\n",
	    "\t$definedSessionLast{$realName}\n";
    }

    # Store the command
    $defined{$realName}="\#$commandName at line $nLine in file $InputFile";
    $definedSessionLast{$realName}=$defined{$realName};
    &store("_command",$commandName);
    
    print "Command \#$commandName defined at line $nLine in file $InputFile\n"
	if $Debug;

    return 1;
}

##############################################################################
sub check_session{

    print "InsideComp: $InsideComp\n" if $Debug;

    if($InsideComp and not $StandAlone){
	&print_error(" for session $nSession:\n".
                "\tBEGIN_COMP $InsideComp does not have matching END_COMP\n");
	$IsError = 1;
    }


    print "required commands: ",join(',',@required),"\n" if $Debug;

    foreach $realName (@required){
	if(not $defined{$realName}){
	    &print_error(" for session $nSession:\n".
		"\trequired command \#$realName was not defined\n");
	    $IsError = 1;
	}
    }

    # Check the rules
    foreach my $node (@{$commandList->{content}}){
	next unless $node->{type} eq 'e';
	if($node->{name} eq 'rule'){

	    print "checking global rule $node->{attrib}->{expr}\n"
		if $Debug;

	    next unless check_if($node);

	    # Substitute $\w+ with the value of the variable in package COMP
	    my $content = $node->{content}->[0]->{content};
	    $content =~ s/^\s+//; $content =~ s/\s+$//;
	    $content =~ s/\$(\w+)/'$COMP::'.$1/gee;

	    print "Rule:\t$node->{attrib}->{expr}\n" if $Verbose;
	    &print_error(" for session $nSession:\n\t$content");
	}
    }
}
##############################################################################
sub store{
    # Usage:
    #      store($name,$value);
    # Store value $value into variable $name in package COMP

    my $name=$_[0];
    my $value=$_[1];
    &eval_comp("\$$name='$value'");
    print "eval_comp(\$$name='$value')\n" if $Debug;
    &print_error(" for command $commandName\n".
		 "\tcould not evaluate \$$name='$value':\n\t$@")
	if $@;
}
##############################################################################
sub extract{
    # Usage:
    #        $result = &extract($expression, $type);
    #
    # $result will get the value extracted from $expression.
    # The variables containing a $ sign are substituted from package
    # COMP and the expressions are evaluated
    # The type given by the second argument can be 
    # "logical", "integer", "real" or "string"
    #
    # Example:
    #       $result = &extract('$nI+2','integer')
    #
    # the function will return the value of
    #
    #     $COMP::nI + 2
    #
    # into $result.

    my $value  = $_[0]; return if length($value) == 0;
    my $type   = $_[1];

    my $value_ = $value;

    if($value_ =~ /\$/){
	if($type eq "string"){
	    $value='"'.$value.'"';
	}else{
	    $value='('.$value.')';
	}
	&eval_comp("\$_value_=$value;");
	$value_ = $COMP::_value_;

        &print_error(" for command \#$commandName\n".
		     "\tcould not evaluate\n".
		     "\tsetting \$value_=$value:\n\t$@")
	    if $@;
    }

    if($type eq "logical"){
	# replace T or .true. with 1, F or .false. with 0
	$value_ =~ s/^(f|.false.)$/0/i;
	$value_ =~ s/^(t|.true.)$/1/i;
	# Convert possible empty string (=false) into a more readable 0
	$value_ = $value_+0;
	if($value_ !~ /^[01]$/){
	    &print_error(" for command \#$commandName\n".
			 "\tcontains incorrect logical\n".
			 "\t$value = $value_");
	    $value_=0;
	}
    }

    return $value_;

}
##############################################################################
sub print_error{
    $IsError = 1;
    print "Error at line $nLine in file $InputFile",@_,"\n";
}

##############################################################################
sub param_error{
    $IsError = 1;
    return if $DontShowParamError;
    &print_error(" for command \#$commandName:\n".
		 "\tParameter $paramName of type $paramType\n".
		 "\t@_");
    print "Command description:\n$commandText{$realName}" if $Verbose;
}
##############################################################################
sub COMP::count_split{
    # This subroutine can be used in the XML file to count the "words" 
    # separater by white space in the argument String.
    # Allow for "array syntax" using "word(30)" instead of 30 words.
    my $String = $_[0];

    # replace "word(\d\d)" with a simple "word" and count extra elements
    my $nExtra = 0;
    $nExtra += $1-1 while $String =~ s/\((\d+)\)//;

    # return total number of words
    return $nExtra + split(' ',$String);
}
##############################################################################
#!QUOTE: \clearpage
#BOP
#!QUOTE: \section{share/Scripts: for SWMF and Physics Modules}
#!QUOTE: \subsection{Checking, Converting and Editing Parameter Files}
#!ROUTINE: CheckParam.pl - check a parameter file based on an XML description
#!DESCRIPTION:
# Before running an executable, it is a good idea to check the parameter file.
# This can save time and frustration.
# This rather complex script is normally called from a customized
# TestParam.pl script. For certain applications, however, it may be
# necessary to call directly this script.
# 
#!REVISION HISTORY:
# 03/22/2004 G.Toth - initial version based on the TestParam.pl script
#                     which was developped for BATSRUS.
#EOP
sub print_help{

    print 
#BOC
"Purpose:

     Read a parameter file and verify its correctness based on an 
     XML description. In most applications this scipt is called by 
     TestParam.pl.

Usage:

  CheckParam.pl [-h] [-H] [-X] [-v] [-D] [-x=XMLFILE]
                [-S] [-c=ID] [-C=IDLIST] [-p=PRECISION] [-i] [PARAMFILE]

  -h            print help message and stop

  -H            print help about the XML tags used in PARAM.XML files and stop

  -X            print a short introduction to the XML language and stop

  -D            print debug information

  -v            verbosity adds command description to every error message.

  -x=XMLFILE    Use XMLFILE for the XML description. 
                If the -x switch is omitted the $XmlFileDefault file is used.

  -S            Check parameters assuming a stand alone mode. In stand alone 
                mode the #BEGIN_COMP and #END_COMP commands are ignored.
                Also sets \$_IsStandAlone to true to be used in the XML rules.

  -c=ID         Check the parameters for the component defined by ID,
                which is the two-character name of the component.
                The default is to check the CON parameters.

  -C=IDLIST     List of registered component ID-s separated with commas.
                The default is not checking the registration.

  -g=GRIDSIZE   List of grid dimensions separated with commas.
                The default is not checking the grid size.

  -p=PRECISION  The default precision for real numbers. Possible values are
                'single' and 'double'. Default value is 'double'.

  -i            Interactive mode. The parameters are read from STDIN.

  PARAMFILE     The file containing the parameters to test. 
                The default file name is run/PARAM.in.
                In interactive mode (-i) or when the XML tree is saved (-s)
                the PARAMFILE is ignored.

Examples:

    Check CON parameters in run/PARAM.in for correctness with verbose info:

CheckParam.pl -C='GM,IH,IE' -v

    Check GM parameters in run/PARAM_new.in for correctness:

CheckParam.pl -x=GM/BATSRUS/PARAM.XML -c=GM run/PARAM_new.in

    Check lines typed through standard input with debug info:

CheckParam.pl -D -i
#DESCRIPTION
Just a test.
...
Ctrl-D"
#EOC
    ,"\n\n";
    exit 0;
}
############################################################################
sub print_help_xml{

    print
#!QUOTE: \clearpage
#BOC
'XML - eXtended Markup Language

XML is a very simple notation used for marking up arbitrary text. 
It is a well defined standard, which is suitable for machine processing,
yet it is relatively easy to read and produce with a normal text editor.

An XML file consists of tags with attributes and text.
The names of the tags and attributes can be any alphanumeric name.
There are 3 special characters which cannot be used in the text
or in the value of the attributes:

   < > &

The < and > are used to start and finish a tag, the & sign is 
used to describe special characters, such as

  &gt;  &lt;  &amp;

which can be used in the text instead of "<", ">" and "&".

A tag without a body is enclosed between the strings "<" and "/>":

<TAGNAME ... />

A tag with an (optional) body has a starting part enclosed between "<" and ">"
and a matching finishing part with the name enclosed between "</" and ">":

<TAGNAME ...>
...
</TAGNAME>

Tags can have attributes, which are names followed by an "=" sign and
the value is quoted with double quotation marks, for example

<TAGNAME ATTRIBUTE1="VALUE1" ATTRIBUTE2="VALUE2" .../>

Tags can be nested and they can contain a text body:

<TAGNAME1>
  <TAGNAME2>
     <TAGNAME3/>
     text2
  </TAGNAME2>
  text1
</TAGNAME1>'
#EOC
    ,"\n\n";
    exit;
}
############################################################################
sub print_help_xml_param{

    print
#!QUOTE: \clearpage
#BOC
'XML Description of Input Parameters

The following shows the structure of the XML tags which can be interpreted by
the CheckParam.pl script, partially used by the XmlToLatex.pl script 
(and it will also be used by the parameter editor of the GUI in the future).

<commandList name=...>
<set name=... type=... value=.../>
<set name=... type=... value=.../>
...
<commandgroup name=...>

Description for command group.

<command name=... [alias=...] [if=...] [required="T"] [multiple="T"]>
	<parameter name=... type="logical" default=.../>
	<parameter name=... type="integer" [min=...] [max=...] [default=...]/>
	<parameter name=... type="real" [min=...] [max=...] [default=...]/>
	<parameter name=... type="string" length=... 
                 [case="upper" | case="lower"] [default=...]/>
	<parameter name=... type="integer|real|string" input="select"
							[multiple="T"]>
		<option name=... [value=...] [default="T"] 
					[exclusive="T"] [if=...]/>
		<optioninput name=... type=... [min=...] [max=...] [length=...]
						[default=...] [if=...]/>
	<parameter/>
	<parameter name=... type="strings" min="..." max="..."
					[ordered="T"] [duplicate="T"]>
		<part name=... type="string" input="select"
					[required="T"] [multiple="T"]>
			<option name=... [value=...] [default="T"] 
					[exclusive="T"] [if=...]/>
		</part>
	</parameter>
	<if expr=...>
		<parameter name=... type="real" [min=...] [max=...] 
							default=... />
	</if>
	<foreach name=... values=...>
		<parameter name=... type="logical" default=.../>
		...
	</foreach>
	<for [index=] from=... to=...>
		<parameter name=... type="integer" [min=...] [max=...]
							default=.../>
		...
	</for>
	<set name=... type=... value=.../>
	<rule expr=...>
	... error message for command rule ...
	</rule>

verbal description of command
...
</command>
...
...
...
<command name=... [required="T"] [if=...]>
...
</command>

</commandgroup>
...
...
...
...
...
...
<commandgroup name=...>
...
...
...
</commandgroup>

<rule expr=...>
... error message for global rule ...
</rule>
<rule expr=...>
... error message for global rule ...
</rule>
<commandList>'
#EOC
    ,"\n\n";
    exit 0;
}
