#!/usr/bin/perl -i
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
use strict;

our $Component       = 'UA';
our $Code            = 'GITM2';
our $MakefileDefOrig = 'srcMake/Makefile.def';
our @Arguments       = @ARGV;

my $config     = "share/Scripts/Config.pl";
if(-f $config){
    require $config;
}else{
    require "../../$config";
}

# These are inherited from $config
our %Remaining;   # Arguments not handled by share/Scripts/Config.pl
our $Show;
our $Help;
our $ERROR;
our $WARNING;
our $NewGridSize;
our $ShowGridSize;
our $Install;
our $Compiler;
our $OS;

&print_help if $Help;

# src/Makefile (needs fixing for Linux/ifort compiler)
my $SrcMakefile     = "src/Makefile";

# Grid size variables
my $NameGridFile="src/ModSize.f90";          # File with actual grid size
my $GridSize;
my ($nLon, $nLat, $nAlt, $MaxBlock);

# Planet variables
my $ModPlanet = "ModPlanet.f90";
my $Planet;
my $PlanetOrig;

my $NoFlush;

&get_settings;

$Planet = "Earth" if $Install; # Default planet
$NoFlush = 0;

foreach (@Arguments){
    if(/^-LV-426$/i)          {$Planet="LV-426";               next};
    if(/^-Titan$/i)           {$Planet="Titan";                next};
    if(/^-mars$/i)            {$Planet="Mars";                 next};
    if(/^-earth$/i)           {$Planet="Earth";                next};
    if(/^-noflush$/i)         {$NoFlush=1     ;                next};
    if(/^-s$/)                {$Show=1;                        next};

    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}

&modify_utilities if $NoFlush;

&install_code if $Install;

&set_grid_size if $NewGridSize and $NewGridSize ne $GridSize;

&set_planet if $Planet and $Planet ne $PlanetOrig;

&show_settings if $Show; 

print "Config.pl -g=$GridSize\n" if $ShowGridSize and not $Show;

exit 0;

############################################################################
sub modify_utilities{

    my $IsDone;
    my $DoWrite;
    my $line;

    print "Altering ModUtilities.f90!\n";

    system "mv share/Library/src/ModUtilities.F90 share/Library/src/ModUtilities.orig.F90";

    open(OUTFILE,">share/Library/src/ModUtilities.f90");
    open(INFILE, "<share/Library/src/ModUtilities.orig.F90");

    while(<INFILE>) {

	if (/subroutine flush_unit/) {
	    print OUTFILE $_;
	    $IsDone = 0;
	    $DoWrite = 1;
	    while (!$IsDone) {
		$line = <INFILE>;
		$DoWrite = 0 if ($line =~ /ifdef/);
		print OUTFILE $line if ($DoWrite);
		$DoWrite = 1 if ($line =~ /endif/);
		$IsDone = 1 if ($line =~ /end subroutine/);
	    }
	} else {
	    print OUTFILE $_;
	}
    }

    close(OUTFILE);
    close(INFILE);

}

############################################################################
sub get_settings{

    # Read size of the grid from $NameGridFile
    open(FILE,$NameGridFile) or die "$ERROR could not open $NameGridFile\n";
    while(<FILE>){
	next if /^\s*!/; # skip commented out lines
        $nLon=$1         if /\bnLons\s*=\s*(\d+)/i;
	$nLat=$1         if /\bnLats\s*=\s*(\d+)/i;
	$nAlt=$1         if /\bnAlts\s*=\s*(\d+)/i;
	$MaxBlock=$1     if /\bnBlocksMax\s*=\s*(\d+)/i;
	last if $nLon and $nLat and $nAlt and $MaxBlock;
    }
    close FILE;

    die "$ERROR could not read nLon from $NameGridFile\n" unless length($nLon);
    die "$ERROR could not read nLat from $NameGridFile\n" unless length($nLat);
    die "$ERROR could not read nAlt from $NameGridFile\n" unless length($nAlt);
    die "$ERROR could not read MaxBlock from $NameGridFile\n" 
	unless length($MaxBlock);

    $GridSize = "$nLon,$nLat,$nAlt,$MaxBlock";

    # Get the current planet
    if(-l "src/$ModPlanet"){
        my $link = `ls -l src/$ModPlanet`;
	$link =~ /Mod(\w+)\.f90$/ or 
	    warn "GITM2/config: Could not find planet in $link";
        $PlanetOrig = $1;
    }
}

#############################################################################

sub install_code{

    return unless $Compiler =~ /ifort/ and $OS =~ /Linux/;
    # Unfix object list for Linux/ifort compiler (this is not kosher)
    @ARGV = ($SrcMakefile);
    while (<>){
	s/-lSHARE.*ModIons.o/-lSHARE/;
	print;
    }
}

#############################################################################

sub set_grid_size{

    print "Writing new grid size $NewGridSize into $NameGridFile...\n";

    $GridSize = $NewGridSize;

    if($GridSize=~/^\d+,\d+,\d+,\d+$/){
	($nLon, $nLat, $nAlt, $MaxBlock)= split(',', $GridSize);
    }elsif($GridSize){
	die "$ERROR -g=$GridSize should be 4 integers separated with commas\n";
    }

    if($MaxBlock and $MaxBlock !~ /^\d+$/){
	die "$ERROR -b=$MaxBlock should be a positive integer\n";
    }

    @ARGV = ($NameGridFile);
    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(nLons\s*=[^\d]*)(\d+)/$1$nLon/i;
	s/\b(nLats\s*=[^\d]*)(\d+)/$1$nLat/i;
	s/\b(nAlts\s*=[^\d]*)(\d+)/$1$nAlt/i;
	s/\b(nBlocksMax\s*=[^\d]*)(\d+)/$1$MaxBlock/i;
	print;
    }

}

############################################################################

sub set_planet{
    $PlanetOrig = $Planet;

    chdir "src" or die "Could not change directory to src\n";

    print "Configuring GITM for $Planet!!\n"; 

    &shell_command("rm -f ModPlanet.f90 ModPlanet.o planet.f90");
    &shell_command("ln -s Mod$Planet.f90 ModPlanet.f90");
    &shell_command("ln -s $Planet.f90 planet.f90");

    my $file;
    foreach $file (glob("*.$Planet.f90")) {
	$file =~ /(.*)\.$Planet\.f90/;
	my $name = $1;
	unlink "$name.f90";
	&shell_command("ln -s $name.$Planet.f90 $name.f90");
    }

    chdir "..";

    &shell_command("cd srcData ; cp UAM.in.$Planet UAM.in");

    if($Planet eq 'Earth'){
          $nLon = 9;
          $nLat = 9;
          $nAlt = 50;
          $MaxBlock = 4;
    }elsif($Planet eq 'LV-426'){
          $nLon = 9;
          $nLat = 9;
          $nAlt = 100;
          $MaxBlock = 4;
    }elsif($Planet eq 'Mars'){
          $nLon = 9;
          $nLat = 9;
          $nAlt = 120;
          $MaxBlock = 4;
    } elsif($Planet eq 'Titan'){
          $nLon = 2;
          $nLat = 2;
          $nAlt = 100;
          $MaxBlock = 4;
    }

    @ARGV = ($NameGridFile);
    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(nLons\s*=[^\d]*)(\d+)/$1$nLon/i;
	s/\b(nLats\s*=[^\d]*)(\d+)/$1$nLat/i;
	s/\b(nAlts\s*=[^\d]*)(\d+)/$1$nAlt/i;
	s/\b(nBlocksMax\s*=[^\d]*)(\d+)/$1$MaxBlock/i;
	print;
    }
}

############################################################################

sub show_settings{

    print "Number of cells in a block: nLon=$nLon, nLat=$nLat, nAlt=$nAlt\n";
    print "Max. number of blocks     : MaxBlock=$MaxBlock\n";
    print "Planet=$PlanetOrig\n";
}

############################################################################

sub print_help{
    print "Additional options for GITM2/Config.pl:

-Titan      Configure GITM2 for Titan. This flag is case insensitive.

-Mars       Configure GITM2 for Mars. This flag is case insensitive.

-LV-426     Configure GITM2 for Testing. This flag is case insensitive.

-Earth      Configure GITM2 for Earth. This flag is case insensitive.

-s          Show current planet.

Additional examples for GITM2/Config.pl:

Install for Titan:

    Config.pl -install -Titan

Install for Mars:

    Config.pl -install -Mars

Reconfigure GITM for Earth:

    Config.pl -Earth

Set grid to nLon=36, nLat=36, nAlt=50 and the number of blocks to 16:

    Config.pl -g=36,36,50,16

Show settings for UA/GITM2:

    Config.pl -s
 ";
    exit 0;
}
