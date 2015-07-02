#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

# Read optional switches -type, -data, -sub
my $DefaultTypes = ($type or "r0");
my $DataFile     = ($file or "SWMF_MPI_routines.dat");
my $Sub          = $sub; # list of MPI subroutines used for debugging mainly
my $Verbose      = $v;

use strict;

my $InFile  = "MpiTemplate.f90";
my $OutFile = "ModMpiInterfaces.f90";

my %TypeName = (i => "integer", 
		r => "real", 
		l => "logical", 
		s => "character(len=*)"
		);

# Set string array to declare dimensions
my @Dims;
my $nDim;
foreach $nDim (1..7){
    $Dims[$nDim] = "(:" . ",:" x ($nDim-1) . ")";
}

print "dimensions=@Dims\n" if $Verbose;

# Set list of subroutines -sub=MPI_BSEND,MPI_RSEND
my @Sub;
my %Types;
if($Sub){
    @Sub = split(/,/,lc($Sub));
    print "Subs=@Sub\n" if $Verbose;
}else{
    open(IN,$DataFile) or die "Could not open data file $DataFile\n";
    while(<IN>){
	next if /^#/;
	next if /^$/;
	chop;
	my $Sub;
	my $Types;
	($Sub, $Types) = split(' ',$_,2);
	$Sub   = lc($Sub);
	$Types = lc($Types);
	print "Sub=$Sub\n" if $Verbose;
	push(@Sub, $Sub);
	$Types = $DefaultTypes unless $Types;
	$Types{$Sub}=$Types;
    }
    close IN;
}

# Create hash 
my %Sub;
foreach (@Sub){$Sub{$_} = 1};

# Read MPI templates from InFile
my %Head; # The first lines of the subroutine
my $IsHead;
my $DoRead;
my $Routine;

open(IN, $InFile) or die "Could not open input file $InFile\n";
while(<IN>) {
    if ( /^\s*SUBROUTINE\s*(\w+)\s*\(/i ) { 
	$Routine=$1;
	next unless $Sub{$Routine};
	$Sub{$Routine}="";
	$DoRead=1;
	$IsHead=1;
    }
    if ( $DoRead ) { 
	$Sub{$Routine} .= $_;
	$Head{$Routine} .= $_;
	$IsHead = 0 unless /\&/;
	$DoRead = 0 if /^\s*END\s*SUBROUTINE/i;
    }
}
close IN;

my $Interface;
my $Procedure;
my $Public;
foreach $Routine ( sort keys %Sub ) {

    $Public .= "  public:: $Routine\n";
    my $Template = $Sub{$Routine}; 

    print "Routine=$Routine\n" if $Verbose;

    print "Template=\n$Template\n\n" if $Verbose;

    # Check if subroutine contains multiple types
    if($Template !~ /<type>/){
	$Procedure .= "  interface\n$Template  end interface\n\n";
	next;
    }

    # Construct call of the original MPI subroutine
    $Template =~ /\)\s*\n/m;
    my $Call = "$`\)";
    $Call =~ s/subroutine/external $Routine\n\n       call/;
	    
    $Procedure .= "  interface $Routine\n    module procedure \&\n";

    my $Types = ($Types{$Routine} or $DefaultTypes);
    my $TypeDims;
    foreach $TypeDims (split(/,/,$Types)) {

	# Construct name of the type specific routine

	my $Type = $TypeDims; 
	my $Dims;
	$Dims = $1 if $Type =~ s/(\d)//;
	my $RoutineType = $Routine."_".$Type;

	# Create template for this variable type
	my $TemplateType = $Template;
	$TemplateType =~ s/<type>/$TypeName{$Type}/g;
	$TemplateType =~ s/$Routine/$RoutineType/g;

	if($Verbose){
	    print "Type=$Type TypeName=$TypeName{$Type}\n";
	    print "Dims=$Dims\n";
	    print "RoutineType=$RoutineType\n";
	    print "TemplateType=\n$TemplateType";
	}

	my $nDim;
	foreach $nDim (0..$Dims) {
	    my $Dim1=$Dims[$nDim];
	    my $Dim2=$Dims[$nDim+1];

	    my $RoutineTypeDim = $RoutineType.$nDim;

	    print "RoutineTypeDim=$RoutineTypeDim\n" if $Verbose;

	    # if( length($RoutineTypeDim) > 31 ) { $RoutineTypeDim =~ s/_//g;}
	    my $TemplateTypeDim = $TemplateType;
	    $TemplateTypeDim =~ s/\(dim1\)/$Dim1/ig;
	    $TemplateTypeDim =~ s/\(dim2\)/$Dim2/ig;
	    $TemplateTypeDim =~ s/$RoutineType/$RoutineTypeDim/g;

	    $TemplateTypeDim =~ s/(end subroutine)/$Call\n     $1/;
	    
	    $Procedure .= "    $RoutineTypeDim, \&\n";
	    $Interface .= "$TemplateTypeDim\n\n";
	}
    }
    $Procedure =~ s/, \&\n$/\n  end interface\n\n/;
}

open(OUT, ">$OutFile") or die "Could not open $OutFile\n";
print OUT
"module ModMPiInterfaces
  implicit none
  private

$Public

$Procedure

contains

$Interface

end module ModMpiInterfaces
";
close OUT;
