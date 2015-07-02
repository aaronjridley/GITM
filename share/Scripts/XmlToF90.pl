#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

my $Help     = ($h or $help); undef $h; 
my $Commands = $c;            undef $c;
my $Force    = $f;            undef $f;

use strict;

&print_help if $Help or not @ARGV;

my $ERROR   = "ERROR in XmlToF90.pl:";
my $WARNING = "WARNING in XmlToF90.pl:";

die "$ERROR cannot use -f and -c switches together!\n" if $Force and $Commands;
die "$ERROR input file argument is missing!\n" unless @ARGV;

require 'share/Scripts/XmlRead.pl';

my $InputFile  = $ARGV[0];
my $OutputFile = $ARGV[1];

open(IN, $InputFile) or die "$ERROR could not open input file $InputFile\n";
my $Input = join('',<IN>);
close(IN);

my $Update = ($OutputFile and -f $OutputFile and not $Force);

my $NewFile = ($OutputFile and not ($Update or $Commands));

my $Indent     = (' ' x 7);  # initial indentation for case statements
my $Indent1    = (' ' x 3);  # incremental indentation
my $IndentCont = (' ' x 5);  # indentation for continuation lines

my $SrcCode;            # F90 code produced by the script
my $SrcBeforeDecl;      # Source code before the parameter declarations
my $SrcDecl;            # Declaration of command parameters
my $SrcDecl2;           # Declarations after the command being processed
my $SrcAfterDecl;       # Source code after the declarations but before indexes
my $SrcIndex;           # Declaration of indexes used in loops
my $SrcBeforeCase;      # Source code after indexes and before case statements
my $SrcCase;            # Case statements reading in the parameters
my $SrcCase2;           # Case statements after the command being processed
my $SrcAfterCase;       # Source code after the case statements
my $iPart;              # part index for multi-part parameters
my $ForeachName;        # name of the index variable in a foreach loop
my $ForeachValue;       # value of the index variable in a foreach loop
my %VariableType;       # Hash for variable types
my %DefaultValue;       # Hash for default values

# put in UseStrict
$VariableType{"UseStrict"} = "logical";
$DefaultValue{"UseStrict"} = "T";

# Read in current file if this is an update. Also make a copy
&read_original if $Update;

my $Xml = &XmlRead($Input); # XML tree data structure created from XML text

&process_xml($Xml);

&check_declarations;

if($NewFile){
    my $Module = $OutputFile;
    $Module = $1 if $Module =~ /\/([^\/]+)$/;  # remove path
    $Module =~ s/\..*$//;                          # remove extension

    $SrcCode = 
"module $Module

  use ModReadParam, ONLY: i_session_read, read_line, read_command, read_var
  use ModUtilities, ONLY: split_string

  implicit none

  character(len=*), parameter:: NameMod = '$Module'
  ! GIPHT BEGIN DECLARATIONS
$SrcDecl
  ! GIPHT END DECLARATIONS

contains

  subroutine set_parameters

    character (len=100) :: NameCommand, StringPart_I(100)
    integer :: iSession, nStringPart
    logical :: UseStrict
    ! GIPHT BEGIN INDEXES
$SrcIndex
    ! GIPHT END INDEXES
    character(len=*), parameter:: NameSub = NameMod//'::set_parameters'
    !-------------------------------------------------------------------------
    iSession = i_session_read()
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       ! GIPHT BEGIN COMMANDS
       select case(NameCommand)
$SrcCase
       ! GIPHT END COMMANDS
       case default
          !if(iProc==0) then
          write(*,*) NameSub // ' WARNING: unknown command ' // &
               trim(NameCommand),' !'
          if(UseStrict)call CON_stop('Correct PARAM.in!')
          !end if
       end select
    end do

  contains

    logical function is_first_session()
      is_first_session = iSession == 1
      if(iSession > 1)then
         ! if(iProc==0) then
         write(*,*) NameSub // ' WARNING: command ',trim(NameCommand), &
              ' can be used in first session only!'
         if(UseStrict)call CON_stop('Correct PARAM.in!')
         ! end if
      end if
    end function is_first_session

  end subroutine set_parameters

end module $Module
";
}elsif($Update){
    $SrcCode  = $SrcBeforeDecl . $SrcDecl. $SrcDecl2 . $SrcAfterDecl;
    $SrcCode .= $SrcIndex if $SrcIndex;
    $SrcCode .= $SrcBeforeCase . $SrcCase . $SrcCase2 . $SrcAfterCase;
}else{
    $SrcCode .= "$SrcDecl"."  !".("-" x 75)."\n" if $SrcDecl;
    $SrcCode .= $SrcCase;
}

if($OutputFile){
    open(OUT, ">$OutputFile") or 
	die "$ERROR could not open output file $OutputFile\n";
    print OUT $SrcCode;
    close(OUT);
}else{
    print $SrcCode;
}

exit 0;

##############################################################################
sub read_original{
    my $Part = "Before decl";
    open(IN, $OutputFile) or die "$ERROR could not open file $OutputFile\n";
    `cp $OutputFile $OutputFile.orig`;
    while(<IN>){
	$Part = "Declaration" if /^ *\! *GIPHT BEGIN DECLARATIONS/;
	$Part = "After decl"  if /^ *\! *GIPHT END DECLARATIONS/;
	$Part = "Index"       if /^ *\! *GIPHT BEGIN INDEXES/;
	$Part = "Before case" if /^ *\! *GIPHT END INDEXES/;
	$Part = "Case"        if /^ *\! *GIPHT BEGIN COMMANDS/;
	$Part = "After case"  if /^ *\! *GIPHT END COMMANDS/;

	$SrcBeforeDecl .= $_ if $Part eq "Before decl";
	$SrcDecl       .= $_ if $Part eq "Declaration";
	$SrcAfterDecl  .= $_ if $Part eq "After decl";
	$SrcIndex      .= $_ if $Part eq "Index";
	$SrcBeforeCase .= $_ if $Part eq "Before case";
	$SrcCase       .= $_ if $Part eq "Case";
	$SrcAfterCase  .= $_ if $Part eq "After case";
    }
}


##############################################################################
sub process_xml{

    # recursive subroutine that processes the XML file

    my $content = shift;

    foreach my $element (@$content){

        next if $element->{type} eq 't'; # Ignore elements of type text

	my $name = lc( $element->{"name"} );

	if($name eq 'commandlist'){
	    &process_xml($element->{content});
	}elsif($name eq 'commandgroup'){
	    my $Name = $element->{attrib}->{name};

	    if($Update){
		# Find commandgroup comment and split there
		my $Src = $SrcCase . $SrcCase2;
		if($Src =~ /( *! >>> $Name <<<\n+)/){
		    $SrcCase  = $`.$1;
		    $SrcCase2 = $';
		}
		my $Src = $SrcDecl . $SrcDecl2;
		if($Src =~ /( *! >>> $Name <<<\n)/ ){
		    $SrcDecl  = $`.$1;
		    $SrcDecl2 = $'; 
		}
	    }else{
		# Remove previous command comment if there were no parameters
		$SrcDecl =~ s/  ! \"\#\w+\"\n$//;
		# Remove previous commandgroup comment if there are no commands
		$SrcDecl =~ s/  ! >>> .*\n$//;
		$SrcDecl .= "\n  ! >>> $Name <<<\n" unless $Commands;
		$SrcCase .= "\n".$Indent."! >>> $Name <<<\n\n" unless $Commands;
	    }
	    &process_xml($element->{content});
	}elsif($name eq 'command'){
	    # Remove previous command comment if there were no parameters
	    $SrcDecl =~ s/  ! \"\#\w+\"\n$//;

	    my $Attrib = $element->{attrib};
	    my $NameCommand = $Attrib->{name};
            my $Name   = "\"\#$NameCommand\"";
	    my $SrcCase1 = $SrcCase;
	    my $SrcDecl1 = $SrcDecl;
	    if($Update){
		# Find command and split SrcCase and SrcDecl
		my $Src = $SrcCase . $SrcCase2;
		if($Src =~ s/(( *)case\($Name.*\n(\2 .*\n)*)//i){
		    $SrcCase1 = $`;
		    $SrcCase  = $`.$1;
		    $SrcCase2 = $';
		}
		my $Src = $SrcDecl . $SrcDecl2;
		if($Src =~ /( *! *$Name.*\n(.*::.*\n)*)/i){
		    $SrcDecl1 = $`;
		    $SrcDecl  = $`.$1;
		    $SrcDecl2 = $'; 
		}
	    }

	    next if $Commands and $Commands !~ /\b$NameCommand\b/;

	    if($Update){
		# Remove previous version of the command if any
		$SrcCase = $SrcCase1;
		$SrcDecl = $SrcDecl1;
	    }

	    my $Alias  = $Attrib->{alias};
	    foreach my $Alias (split ",", $Attrib->{alias}){
		$Name .= ", \"\#$Alias\"";
	    }
	    # Add new comment
	    $SrcDecl .= "  ! $Name\n";
	    $SrcCase .= $Indent."case($Name)\n";
	    $Indent .= $Indent1;
	    my $If = $Attrib->{if};
	    if($If =~ /IsFirstSession/){
	    	$SrcCase .= $Indent."if(.not.is_first_session())CYCLE\n";
	    }
	    &process_xml($element->{content});
	    $Indent =~ s/$Indent1//;
	}elsif($name eq 'parameter'){
	    my $Attrib = $element->{attrib};
	    my $Name   = perl_to_f90($Attrib->{name});
	    my $Type   = $Attrib->{type};
	    my $Case   = $Attrib->{case};
	    my $If     = perl_to_f90($Attrib->{if});

            if($Type eq "string"){
		my $Length = $Attrib->{length};
		$Type = "character(len=$Length)" if $Length;
	    }

	    # Store variable
	    &add_var($Name, $Type, $Attrib->{default});

	    # Create line
	    $SrcCase .= $Indent."if($If) &\n".$IndentCont if $If;
	    $SrcCase .= $Indent."call read_var('$Name', $Name)\n";

	    $SrcCase =~ s/\)\n$/, IsUpperCase=.true.\)\n/ if $Case eq "upper";
	    $SrcCase =~ s/\)\n$/, IsLowerCase=.true.\)\n/ if $Case eq "lower";

	    if($Type eq "strings"){
		my $MaxPart = $Attrib -> {max};
		$SrcCase .= $Indent . "call split_string" . 
		    "($Name, $MaxPart, StringPart_I, nStringPart)\n";
		$iPart = 0;
		&process_xml($element->{content});
	    }
	}elsif($name eq 'part'){
	    my $Name = $element->{attrib}->{name};
	    &add_var($Name, "string");
	    $iPart++;
	    $SrcCase .= $Indent."$Name = StringPart_I($iPart)\n";
	}elsif($name eq 'if'){
	    my $Expr = perl_to_f90($element->{attrib}->{expr});
	    $SrcCase .= $Indent."if($Expr)then\n";
	    $Indent .= $Indent1;
	    &process_xml($element->{content});
	    $Indent =~ s/$Indent1//;
	    $SrcCase .= $Indent."end if\n";
	}elsif($name eq 'for'){
	    my $Attrib = $element->{attrib};
	    my $From  =  perl_to_f90($Attrib->{from});
	    my $To    =  perl_to_f90($Attrib->{to});
	    my $Index = (perl_to_f90($Attrib->{name}) or "i");

	    # Declare index variable
	    $SrcIndex .= "    integer :: $Index\n" 
		unless $SrcIndex =~ /integer :: $Index\b/i;

	    $SrcCase .= $Indent."do $Index = $From, $To\n";
	    $Indent .= $Indent1;
	    &process_xml($element->{content});
	    $Indent =~ s/$Indent1//;
	    $SrcCase .= $Indent."end do\n";
	}elsif($name eq 'foreach'){
	    my $Attrib = $element->{attrib};
	    $ForeachName = $Attrib->{name};
	    foreach (split(/,/, $Attrib->{values})){
		$ForeachValue = $_;
		&process_xml($element->{content});
	    }
	    $ForeachName = ''; $ForeachValue = '';
	}
    }
}

###############################################################################

sub perl_to_f90{

    $_ = shift;

    # replace special variables provided by the CheckParam.pl script
    s/\$_command/NameCommand/ig;
    s/\$_namecomp/NameComp/ig;

    # replace foreach variable with actual value
    s/\$$ForeachName/$ForeachValue/g if $ForeachName;

    # remove all dollar signs from variable names
    s/\$//g;

    # convert relation operator
    s/ eq / == /g;
    s/ ne / \/= /g;
    s/ and / .and. /g;
    s/ or / .or. /g;
    s/ != / \/= /g;
    s/\bnot /.not. /g;

    # replace string matching (this is not quite right!)
    s/(\w+)\s*=~\s*\/([^\/]+)\/i?/index($1, "$2") > 0/g;

    s/(\w+)\s*\!~\s*\/([^\/]+)\/i?/index($1, "$2") < 1/g;

    # Remove \b from patterns (this is not quite right!)
    s/\\b//g;

    return $_;
}

###############################################################################

sub add_var{

    my $Name  = shift;
    my $Type  = shift;
    my $Value = shift;

    $Type =~ s/strings?/character(len=100)/;

    # Fix value
    $Value =~ s/T/.true./ if $Type eq "logical";
    $Value =~ s/F/.false./ if $Type eq "logical";
    $Value =  "'$Value'" if $Type =~ /character/ and length($Value);
    $Value .= ".0" if $Type eq "real" and $Value =~ /^[\+\-\d]+$/;

    # Add declaration
    if(length($Value)){
	$SrcDecl .= "  $Type :: $Name = $Value\n";
    }else{
	$SrcDecl .= "  $Type :: $Name\n";
    }
}

###############################################################################

sub check_declarations{


    # Check and eliminate duplicate declarations
    my @Decl = split("\n",$SrcDecl . $SrcDecl2);
    foreach (@Decl){
	next unless /^ *(integer|real|logical|character)/i;
	my $Type = lc($1);

	my $Name;
	my $Value;
	if(/(\w+) *= *([^\(\)]*)$/){
	    $Name  = $1;
	    $Value = $2;
	}else{
	    /(\w+)$/;
	    $Name = $1;
	}

	my $name = lc($Name);
	my $Type2 = $VariableType{$name};
	if($Type2){
	    # Check if types agree
	    if($Type ne $Type2){
		warn "$WARNING: variable $Name has types $Type and $Type2\n";
		$_ = "!!! $Type :: $Name was declared above with $Type2\n";
		next;
	    }
	    # Check if default values agree
	    my $Value2 = $DefaultValue{$name};

	    if(length($Value) and length($Value2) and ($Value ne $Value2)){
		warn "$WARNING: variable $Name has default values ".
		    "$Value and $Value2\n";
		$_ = "!!! $Type :: $Name = $Value ".
		    "was declared above with value $Value2\n";
		next;
	    }
	    $_ = "  ! $Type :: $Name has been declared above\n";
	    next;
	}
	$DefaultValue{$name} = $Value if length($Value);
	$VariableType{$name} = $Type;
    }

    # Replace $VarXyz patterns with their default values
    foreach (@Decl){
	next unless /^ *(integer|real|logical|character)/i;

	while(/\$(\w+)/){
	    my $Var2   = $1;
	    my $Value2 = $DefaultValue{lc($Var2)};
	    if(length($Value2)){
		s/\$(\w+)/$Value2/;
	    }else{
		warn "$WARNING: default value of $Var2 could not be found\n";
		s/\$(\w+)/?$Var2?/;
	    }
	}
    }

    # Put together the declaration string
    $SrcDecl = join("\n", @Decl);
    $SrcDecl2 = '';

}

###############################################################################
#BOP
#!ROUTINE: XmlToF90.pl - generate F90 source from XML definitions of input parameters
#!DESCRIPTION:
# Generate F90 source code based on the XML definitions of the
# input parameters typically found in the PARAM.XML files.
# This script allows to store the parameter descriptions in a single XML file
# which is suitable for automated parameter checking, manual and GUI 
# generation, as well as generating the F90 code that reads in the parameters.
# The specific format of the PARAM.XML files is described by the
# share/Scripts/CheckParam.pl script and the manual.
#
#!REVISION HISTORY:
# 05/27/2008 G.Toth - initial version
#EOP
sub print_help{

    print 
#BOC
"Purpose:

   Convert XML description of input commands into F90 source code.

Usage:

   XmlToF90 [-h] [-f | -c=COMMANDS] infile [outfile]

-h          This help message

-c=COMMANDS COMMANDS is a comma separated list of commands to be transformed
            from XML to F90. Default is transforming all the commands.

-f          Force creating a new output file even if it already exists. 
            This flag cannot be combined with the -c=COMMAND switch.
            Default is to update the (selected) commands only in an existing
            file.

infile      Input XML file.

outfile     Output F90 file. Default is writing to STDOUT.

Examples:

            Generate F90 source code for a few commands and write to STDOUT:

share/Scripts/XmlToF90.pl -c=MAGNETICAXIS,ROTATIONAXIS Param/PARAM.XML

            Update all commands in an existing F90 source file

share/Scripts/XmlToF90.pl PARAM.XML src/set_parameters.f90"
#EOC
    ,"\n\n";
    exit 0;
}
