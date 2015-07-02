#!/usr/bin/perl -s

my $Help = ($h or $help);
my $CompVersion = $v;
my $OutFile = ($o or "OUTPUT.XML");


use strict;

&print_help if $Help or not @ARGV;

# Default versions
my %ParamXml = ("CON" => "Param/PARAM.XML",
		"EE"  => "GM/BATSRUS/PARAM.XML",
		"GM"  => "GM/BATSRUS/PARAM.XML",
		"IH"  => "GM/BATSRUS/PARAM.XML",
		"OH"  => "GM/BATSRUS/PARAM.XML",
		"SC"  => "GM/BATSRUS/PARAM.XML",
		"IE"  => "IE/Ridley_serial/PARAM.XML",
		"IM"  => "IM/RCM2/PARAM.XML",
		"PS"  => "PS/DGCPM/PARAM.XML",
		"PW"  => "PW/PWOM/PARAM.XML",
		"RB"  => "RB/RBE/PARAM.XML",
		"SP"  => "SP/Kota/PARAM.XML",
		"UA"  => "UA/GITM2/PARAM.XML",
    );

my $comp;
my $IsSwmf;
if(-f "Param/PARAM.XML"){
    $comp = "CON";
    $IsSwmf = 1;
}elsif(-f "PARAM.XML"){
    $comp = "MODEL";
    $IsSwmf = 0;
}else{
    die("Could not find Param/PARAM.XML or PARAM.XML file\n");
}

# Set the XML file names if -v=... argument is used
if($IsSwmf and $CompVersion){ 
    my $comp;
    my $version;
    my $compversion;
    foreach $compversion (split (/,/, $CompVersion)){
	($comp, $version) = split(/\//, $compversion);
	$ParamXml{$comp} = "$comp/$version/PARAM.XML";
    }
}

my %command;

my $InFile;
foreach $InFile (@ARGV){
    open(INFILE, $InFile);
    while(<INFILE>){
	$comp = 'CON'          if /^\#END_COMP/;
	$command{$comp}{$1}++ if /^\#(\w+)/;
	$comp = $1             if /^\#BEGIN_COMP (..)/;
    }
    close(INFILE);
}

open(OUTFILE, ">$OutFile") or die "Could not open OutFile=$OutFile\n";

my $ParamXml = "PARAM.XML"; # default for standalone code
my $command;
foreach $comp (sort keys %command){
    if($IsSwmf){
	$ParamXml = $ParamXml{$comp} or 
	    print "There is no PARAM.XML file for compomnent $comp\n" && next;
    }
    open(INFILE, $ParamXml) or die "Could not open $ParamXml\n";
    my $Xml = join("", <INFILE>);
    close(INFILE);

    print OUTFILE $1 if $Xml =~ /\n(<commandList.*\n)/;

    foreach $command (sort keys %{ $command{$comp} }){
	next if $command =~ /BEGIN_COMP|END_COMP/;
	if($Xml =~ /\n(<command name=["\w =,]*?"$command".*?<\/command>)/s){
	    print OUTFILE "$1\n";
	}else{
	    print "Command $command is not described in $ParamXml\n";
	}
    }
    print OUTFILE "</commandList>\n";
}

close(OUTFILE);

exit 0;
################################################################
sub print_help{

    print "
ParamToXml.pl can collect commands from one or more 
PARAM.in files and extract the corresponding XML
description from the PARAM.XML files. The script
should be run from the main SWMF directory
or from the main directory of one of the models, e.g. GM/BATSRUS/.
The resulting XML output can be read directly or further
processed into a LaTex file with the XmlToTex.pl script.

Usage: ParamToXml.pl [-h] [-o=OUTFILE] [-v=COMPVERSION] PARAM1 [PARAM2 ...]

-h -help        show help message.
-o=OUTFILE      set name of output XML file. Default is OUTPUT.XML
-v=COMPVERSION  set the component versions used in the PARAM.in files.
                This is only needed for components with multiple models,
                for example IM. COMPVERSION should be a comma separated list.
                Default is IM/RCM2 for the IM component.
PARAM1 ...      Names of the PARAM.in file containing the commands of interest.

Exmaples of use:

Create XML description of the commands used by SWPC and convert it to LaTex:

    share/Scripts/ParamToXml.pl -o=SWPC.XML Param/SWPC/PARAM.in*
    share/Scripts/XmlToTex.pl SWPC.XML > SWPC.tex

Create XML description of the commands in run/PARAM.in when IM/CRCM is used:

    share/Scripts/ParamToXml.pl -v=IM/CRCM run/PARAM.in

";
    exit 0;
}
