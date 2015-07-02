#!/usr/bin/perl -i
#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
use strict;

# Default compiler per machine or OS
my %Compiler = (
		"Linux"               => "nagfor",
		"Darwin"              => "nagfor",
		"OSF1"                => "f90",
		"IRIX64"              => "f90",
		"AIX"                 => "xlf90",
		"palm"                => "ifort",
		"cfe"                 => "ifort",
		"pfe"                 => "ifort,icc",
		"sysx"                => "xlf90",
		"nyx-login-intel"     => "ifortmpif90",
		"nyx-login-amd"       => "ifortmpif90",
		"flux-login"          => "ifortmpif90",
		"hera"                => "mpiifort",
		"ubgl"                => "mpxlf90,mpxlc",
		"jaguarpf-ext"        => "ifortftn",
                "kraken-gsi"          => "ifortftn",
                "yslogin"             => "ifortmpif90,icc",
                "h2ologin"            => "crayftn,cc",
		);

my $WARNING_='share/Scripts/Config.pl WARNING:';
my $ERROR_  ='share/Scripts/Config.pl ERROR:';

my $ErrorCode;   # Return value of system(...)
my $IsStrict=1;  # If true, shell_command will stop on error

# Obtain $OS, $DIR, and the machine name and provide it to caller script
our $OS  = `uname`    or die "$ERROR_ could not obtain OS\n"; chop $OS;
our $DIR = `/bin/pwd` or die "$ERROR_ could not obtain DIR\n"; chop $DIR;
our $Machine = `hostname`; chop($Machine); $Machine =~ s/\..*//;

# remove numbers from the machine name
$Machine =~ s/\d+$//; 

# These are either obtained from the calling script or set here
our $Component;             # The SWMF component the code is representing
our $Code;                  # The name of the code
($Component, $Code) = ($DIR =~ /([A-Z][A-Z])\/([^\/]+)$/)
    unless $Code;

# Strings for the error and warning messages for the caller script
our $ERROR   = "$Code/Config.pl ERROR:";
our $WARNING = "$Code/Config.pl WARNING:";

# Obtain the default Fortran compiler part of Makefile.conf
our $Compiler;
$Compiler = $Compiler{$Machine} or $Compiler = $Compiler{$OS} or
    die "$ERROR_ default compiler is not known for OS=$OS\n";

# Default C compiler part of Makefile.conf
our $CompilerC = "gcc_mpicc";
$CompilerC = $1 if $Compiler =~ s/,(.+)//;

# These are always obtained from the calling script
our $MakefileDefOrig;       # Original Makefile.def 
our @Arguments;             # Arguments obtained from the caller script

# The arguments not handled by this script are provided to the caller
our %Remaining;

# These file names are provided to the calling script
our $MakefileDef      = 'Makefile.def';
our $MakefileConf     = 'Makefile.conf';
our $MakefileConfOrig = 'share/build/Makefile';
our $MakefileRules    = 'Makefile.RULES';
our $MakefileDepend   = 'Makefile.DEPEND';

# These options are set here and provided to the calling script
our $Help;                  # Print help message
our $Verbose;               # Verbose information is printed if true
our $Show;                  # Show information
our $DryRun;                # True if no change is actually made
our $Precision='unknown';   # Precision set in $MakefileConf
our $Installed;             # true if code is installed ($MakefileConf exists)
our $Install;               # True if code is (re)installed
our $Uninstall;             # True if code is uninstalled
our $ShowGridSize;          # Show grid size for caller code
our $NewGridSize;           # New grid size to be set in caller code
our $Hdf5;                  # True if HDF5  lib is enabled
our $Hypre;                 # True if HYPRE lib is enabled
our $Spice;                 # True if SPICE lib is enabled
our $Fishpak;               # True if Fishpak lib is enabled

# The name of the parallel HDF5 Fortran and C compilers
my $H5pfc = "h5pfc";
my $H5pcc = "h5pcc";

# This string should be added into Makefile.conf when HYPRE is enabled
my $HypreDefinition = "# HYPRE library definitions
HYPRELIB     = -L\${UTILDIR}/HYPRE/lib -lHYPRE
HYPRESEARCH  = -I\${UTILDIR}/HYPRE/include
";             	    

# This string should be added into Makefile.conf when Fishpak is enabled
my $FishpakDefinition = "# Fishpak library definitions                               
FISHPAKSEARCH  = -I\${UTILDIR}/FISHPAK/lib
FISHPAKLIB     = -L\${UTILDIR}/FISHPAK/lib -lFISHPAK           
";



# Default precision for installation
my $DefaultPrecision = 'double';

# Global variables for the settings
my $IsComponent=0;         # True if code is installed as a component of SWMF
my $NewPrecision;
my $NewOptimize;
my $NewDebug;
my $NewMpi;
my $NewHdf5;
my $NewHypre;
my $NewFishpak;
my $NewSpice;
my $IsCompilerSet;
my $Debug;
my $Mpi;
my $Fcompiler;
my $Ccompiler;
my $MpiCompiler;
my $MpiHeaderFile = "share/Library/src/mpif.h";
my $Optimize;
my $ShowCompiler;

# Obtain current settings
&get_settings_;

# Show current settings if no arguments are given.
$Show = 1 if not @Arguments;

# Set actions based on the switches
foreach (@Arguments){
    if(/^-dryrun$/)           {$DryRun=1;                       next};
    if(/^-verbose$/i)         {$Verbose=1;                      next};
    if(/^-h(elp)?$/i)         {$Help=1;                         next};
    if(/^-show$/i)            {$Show=1;                         next};
    if(/^-(single|double)$/i) {$NewPrecision=lc($1);            next};
    if(/^-install(=.*)?$/)    {my $value=$1;
			       $IsComponent=1 if $value =~ /^=c/i;
			       $IsComponent=0 if $value =~ /^=s/i;
			       $Install=1;                      next};
    if(/^-uninstall$/i)       {$Uninstall=1;                    next};
    if(/^-compiler=(.*)$/i)   {$Compiler=$1; 
			       $CompilerC=$1 if $Compiler =~ s/,(.+)//;
			       $IsCompilerSet=1;  next};
    if(/^-compiler$/i)        {$ShowCompiler=1;                 next};
    if(/^-mpi$/i)             {$NewMpi="yes";                   next};
    if(/^-nompi$/i)           {$NewMpi="no";                    next};
    if(/^-debug$/i)           {$NewDebug="yes";                 next};
    if(/^-nodebug$/i)         {$NewDebug="no";                  next};
    if(/^-hdf5$/i)            {$NewHdf5="yes";                  next};
    if(/^-nohdf5$/i)          {$NewHdf5="no";                   next};
    if(/^-hypre$/i)           {$NewHypre="yes";                 next};
    if(/^-fishpak$/i)         {$NewFishpak="yes";               next};
    if(/^-nohypre$/i)         {$NewHypre="no";                  next};
    if(/^-nofishpak$/i)       {$NewFishpak="no";                next};
    if(/^-spice=(.*)$/i)      {$NewSpice=$1;                    next};
    if(/^-nospice$/i)         {$NewSpice="no";                  next};
    if(/^-O[0-5]$/i)          {$NewOptimize=$_;                 next};  
    if(/^-g(rid)?$/)          {$ShowGridSize=1;                 next};
    if(/^-g(rid)?=([\d,]+)$/) {$NewGridSize=$+;                 next};

    if(/^-g(rid)?=(.*)$/ and $IsComponent){
	die "$ERROR: incorrect grid size -g=$+\n"};

    $Remaining{$_}=1;
}

if(not $MakefileDefOrig and not $IsComponent){
    warn "$WARNING $Code cannot be used in stand alone mode!\n ".
	"Switching to component mode...\n";
    $IsComponent = 1;
}

&print_help_ if $Help;

if($Uninstall){
    if(not $Installed){
	warn "$ERROR_ $Code is not installed.\n";
	exit 1;
    }else{
	&shell_command("cd share; make distclean")
	    if -f "share/Makefile" and not $IsComponent;
	&shell_command("cd util; make distclean")
	    if -d "util" and not $IsComponent;
	&shell_command("make allclean");
	&shell_command("rm -f Makefile.def Makefile.conf dataCRASH ".
		       "src*/$MakefileDepend src*/$MakefileRules");
	&shell_command("rm -f data") if -l "data";
	exit 0;
    }
}

if($ShowCompiler){
    my @File= glob($MakefileConfOrig.'*');
    print "List of Fortran compilers for OS=$OS:\n";
    foreach (@File){
	next unless s/$MakefileConfOrig\.$OS\.//;
	print "  $_\n";
	$_ = ""; # remove the compiler from the list
    }
    print "List of C compilers:\n";
    foreach (@File){
	next unless $_ and not /$MakefileConfOrig\.(\w+\.|conf|doc)/;
	s/$MakefileConfOrig\.//;
        print "  $_\n";
    }
    exit 0;
}

# Execute the actions in the appropriate order
&install_code_ if $Install;

# Change to the main directory of the SWMF if called from a component
chdir "../.." if $IsComponent;

# Check if Makefile.def is up to date
if(-f $MakefileDef){
    my @Stat = stat $MakefileDef;
    my $Time = $Stat[9];
    @Stat = stat $MakefileDefOrig;
    my $TimeOrig = $Stat[9];
    die "$ERROR $MakefileDefOrig is newer than $MakefileDef !\n".
	"   Reinstall or merge changes into $MakefileDef !\n"
	if $Time < $TimeOrig;
}

# Check if Makefile.conf is up to date
if(-f $MakefileConf){
    my @Stat = stat $MakefileConf;
    my $Time = $Stat[9];
    foreach ("$OS.$Compiler", $CompilerC){
	my $Makefile = $MakefileConfOrig.".".$_;
	next unless -f $Makefile;
	my @Stat = stat $Makefile;
	my $TimeOrig = $Stat[9];
	die "$ERROR $Makefile is newer than $MakefileConf !\n".
	    "   Reinstall or merge changes into $MakefileConf !\n"
	    if $Time < $TimeOrig;
    }
}

# Change precision of reals if required
if($NewPrecision and $NewPrecision ne $Precision){
    &shell_command("make clean");
    &set_precision_;
}

# Change debugging flags if required
&set_debug_ if $NewDebug and $NewDebug ne $Debug;

# Change optimization level if required
&set_optimization_ if $NewOptimize and $NewOptimize ne $Optimize;

# Link with MPI vs. NOMPI library if required
&set_mpi_ if $NewMpi and $NewMpi ne $Mpi;

# Link with HDF5 library is required
&set_hdf5_ 
    if ($Install and not $IsComponent) or ($NewHdf5 and $NewHdf5 ne $Hdf5);

# Link with HYPRE library is required
&set_hypre_ if $NewHypre and $NewHypre ne $Hypre;

# Link with FISHPAK library is required 
&set_fishpak_ if $NewFishpak and $NewFishpak ne $Fishpak;

# Link with SPICE library is required
&set_spice_ 
    if ($Install and not $IsComponent) or ($NewSpice and $NewSpice ne $Spice);

# Get new settings
&get_settings_;

# Show settings if required
&show_settings_ if $Show;

# Recreate Makefile.RULES with the current settings
&create_makefile_rules;

# Return into the component directory
chdir $DIR if $IsComponent;

# DO NOT USE exit HERE as this code is called from other perl scripts!

##############################################################################
sub get_settings_{

    $Installed = (-e $MakefileConf and -e $MakefileDef);

    return if not $Installed;

    # Set defaults/initial values
    $Precision   = "unknown";

  TRY:{
      # Read information from $MakefileDef
      open(MAKEFILE, $MakefileDef)
	  or die "$ERROR_ could not open $MakefileDef\n";

      while(<MAKEFILE>){
	  if(/^\s*include\s+(.*$MakefileDef)\s*$/){
	      $MakefileDef = $1;
	      $IsComponent = 1;
	      close MAKEFILE;
	      redo TRY;
	  }
	  $OS         = $1 if /^\s*OS\s*=\s*(\w+)/;
      }
      close(MAKEFILE);
  }

    $Debug     = "no";
    $Mpi       = "yes";
    $Hdf5      = "no";
    $Hypre     = "no";
    $Fishpak   = "no";
    $Spice     = "no";
  TRY:{
      # Read information from $MakefileConf
      open(MAKEFILE, $MakefileConf)
	  or die "$ERROR_ could not open $MakefileConf\n";

      while(<MAKEFILE>){
	  if(/^\s*include\s+(.*$MakefileConf)\s*$/){
	      $MakefileConf = $1;
	      close MAKEFILE;
	      redo TRY;
	  }
	  $Compiler = $1 if 
	      /^\#\s*Fortran language.*Makefile\.$OS\.(\S+)/i;
	  $CompilerC= $1 if /^\#\s*C language.*:\s*Makefile\.(\S+)/;

	  $Fcompiler = $+ if 
	      /^\s*COMPILE\.f90\s*=\s*(\$\{CUSTOMPATH_F\})?(\S+)/;
	  $Ccompiler   = $1 if /^\s*COMPILE\.c\s*=\s*(\S+)/;
	  $MpiCompiler = $1 if /^\s*LINK\.f90\s*=\s*(.*)/;

	  $Precision = lc($1) if /^\s*PRECISION\s*=.*(SINGLE|DOUBLE)PREC/;
          $Debug = "yes" if /^\s*DEBUG\s*=\s*\$\{DEBUGFLAG\}/;
	  $Mpi   = "no"  if /^\s*MPILIB\s*=.*\-lNOMPI/;
	  $Hdf5  = "yes" if /^\s*LINK\.f90\s*=.*$H5pfc/;
	  $Hypre = "yes" if /^\s*HYPRELIB/;
	  $Fishpak = "yes" if /^\s*FISHPAKLIB/;
	  $Spice = "$1"  if /^\s*SPICELIB\s*=\s*(\S*)/;
          $Optimize = $1 if /^\s*OPT[0-5]\s*=\s*(-O[0-5])/;
      }
    }
    close(MAKEFILE);

    # Fix these if the Fortran language and C language lines were missing
    $Compiler  = $Fcompiler if not $Compiler;
    $CompilerC = $Ccompiler if not $CompilerC;

    # Fix MpiCompiler definition if needed
    $MpiCompiler =~ s/\{COMPILE.f90\}\#\s*/\{CUSTOMPATH_MPI\}/;

    # Remove the commented out name of the original linker when h5pfc is used
    $MpiCompiler =~ s/$H5pfc \#.*$/$H5pfc/;

}

##############################################################################

sub show_settings_{

    if(not $Installed){
	print "$Code is not installed\n";
	exit 0;
    }

    print "\n";
    print "$Code is installed in directory $DIR\n";

    if($IsComponent){
	print "    as the $Component component.\n";
    }else{
	if(-e $MpiHeaderFile){
	    print "    as a stand-alone code for serial execution.\n";
	}else{
	    print "    as a stand-alone code.\n";
	}
    }
    print "The installation is for the $OS operating system.
Makefile.conf was created from $MakefileConfOrig.$OS.$Compiler 
                           and $MakefileConfOrig.$CompilerC
The selected F90 compiler is $Fcompiler.
The selected C   compiler is $Ccompiler.
The default precision for reals is $Precision precision.
The maximum optimization level is $Optimize
Debugging flags:   $Debug
Linked with MPI:   $Mpi
Linked with HDF5:  $Hdf5
Linked with HYPRE: $Hypre
Linked with FISHPAK: $Fishpak
Linked with SPICE: $Spice
";

}

##############################################################################
sub install_code_{

    my $Text = $Installed ? "Reinstalling $Code" : "Installing $Code";
    $Text .= " as $Component component" if $IsComponent;  
    print "$Text\n";

    if($IsComponent){
	my $dir = $DIR; $dir =~ s|/[^/]*/[^/]*$||;  # go two directories up
	my $makefile = "$dir/Makefile.def";          # makefile to be included
	die "$ERROR_ could not find file $makefile\n" unless -f $makefile;
	&shell_command("echo include $makefile > Makefile.def");

	$makefile = "$dir/Makefile.conf"; # makefile to be included
	die "$ERROR_ could not find file $makefile\n" unless -f $makefile;
	&shell_command("echo include $makefile > Makefile.conf");
    }else{
	die "$ERROR_ original $MakefileDef is not given\n" unless
	    $MakefileDefOrig;
	die "$ERROR_ $MakefileDefOrig is missing\n" unless
	    -f $MakefileDefOrig;
	&shell_command("echo OS=$OS > $MakefileDef");
	&shell_command("echo MYDIR=$DIR >> $MakefileDef");
	&shell_command("echo ${Component}DIR=$DIR >> $MakefileDef");
	&shell_command("echo COMPILER=$Compiler >> $MakefileDef");
	&shell_command("cat $MakefileDefOrig >> $MakefileDef");

	my $Makefile = "$MakefileConfOrig.$OS.$Compiler";
	if(-f $Makefile){
	    &shell_command("cat $Makefile > $MakefileConf");
	}else{
	    # Try to use generic Makefile with provided compiler
	    warn "$WARNING_: $Makefile was not found,".
		" using generic $MakefileConfOrig.conf\n";
	    sleep 10;
	    $Makefile = "$MakefileConfOrig.conf";
	    open(IN, $Makefile) or die "$ERROR_ $Makefile is missing\n";
	    open(OUT, ">$MakefileConf") 
		or die "$ERROR_ could not open $MakefileConf\n";
	    while(<IN>){
		s/_COMPILER_/$Compiler/;
		s/_OS_/$OS/;
		print OUT $_;
	    }
	    close IN; close OUT;
	}

	# Append the C compiler
	$Makefile = "$MakefileConfOrig.$CompilerC";
	if(-f $Makefile){
            &shell_command("cat $Makefile >> $MakefileConf");
	}else{
	    die "$ERROR_ could not find $Makefile\n";
	}

	# Remove -lmpicxx from CPPLIB definition in Makefile.conf if not needed
	my $remove_mpicxx = (`mpicxx -show` !~ /\-lmpi_cxx/);
	if($remove_mpicxx){
	    @ARGV = ($MakefileConf);
	    while(<>){
		s/ -lmpi_cxx// if /^CPPLIB/;
		print;
	    }
	}
    }

    # Read info from main Makefile.def
    &get_settings_;

    # Set initial precision for reals
    $NewPrecision = $DefaultPrecision unless $NewPrecision;
    &set_precision_ if $NewPrecision ne $Precision;

    # Create Makefile.RULES as needed
    &create_makefile_rules;

    # Link data and dataCRASH directories if possible
    &link_swmf_data;

    # Install the code
    &shell_command("cd share; make install") 
	if -f "share/Makefile" and not $IsComponent;
    &shell_command("cd util; make install") 
	if -d "util" and not $IsComponent;
    &shell_command("make install");

    # Now code is installed
    $Installed = 1 unless $DryRun;
}

##############################################################################

sub set_precision_{

    # Set the precision for reals in $MakefileConf

    # Precision will be NewPrecision after changes
    $Precision = $NewPrecision;

    my $PREC = '${'.uc($Precision).'PREC}';
    print "Setting PRECISION variable to $PREC in $MakefileConf\n";
    if(not $DryRun){
	@ARGV = ($MakefileConf);
	while(<>){
	    s/^(\s*PRECISION\s*=\s*).*/$1$PREC/;
	    print;
	}
    }
}

##############################################################################

sub set_debug_{

    # Set the debug compilation flags in $MakefileConf

    # Debug will be NewDebug after changes
    $Debug = $NewDebug;

    my $DEBUG; $DEBUG = '${DEBUGFLAG}' if $Debug eq "yes";
    print "Setting debugging flags to '$Debug' in $MakefileConf\n";
    if(not $DryRun){
	@ARGV = ($MakefileConf);
	while(<>){
	    s/^(\s*DEBUG\s*=).*/$1 $DEBUG/;
	    print;
	}
    }
}

##############################################################################

sub set_mpi_{

    if(-e $MpiHeaderFile){
	warn "$WARNING code was installed with -nompi, cannot switch on MPI\n";
	return;
    }

    # $Mpi will be $NewMpi after changes
    $Mpi = $NewMpi;

    if($Mpi eq "no" and $Install){
	&shell_command("cp share/include/nompif.h $MpiHeaderFile");
	$MpiCompiler = '${COMPILE.f90}';
    }

    # Select the MPI or NOMPI library in $MakefileConf

    print "Selecting MPI library in $MakefileConf\n" if $Mpi eq "yes";
    print "Selecting NOMPI library in $MakefileConf\n" if $Mpi eq "no";
    if(not $DryRun){
	@ARGV = ($MakefileConf);
	while(<>){
	    # Modify LINK.f90 definition
	    if(/^\s*LINK.f90\s*=/){
		s/\{CUSTOMPATH_MPI\}/\{COMPILE.f90\}\# \t/ if $Mpi eq "no";
		s/\{COMPILE.f90\}\#\s*/\{CUSTOMPATH_MPI\}/ if $Mpi eq "yes";
		$_ = 'LINK.f90 = ${COMPILE.f90}'."\n" if $Install;
	    }

	    # Comment/uncomment MPILIB definitions
	    if(/MPILIB\s*=/){
		s/^\s*M/\#M/ if /lNOMPI/ eq ($Mpi eq "yes");
		s/^\#\s*M/M/ if /lNOMPI/ eq ($Mpi eq "no");
	    }

	    # Comment/uncomment mpi_cxx library
	    if(/mpi_cxx/){
		s/ \-lmpi_cxx/ \#\-lmpi_cxx/ if $Mpi eq "no";
		s/ \#\-lmpi_cxx/ \-lmpi_cxx/ if $Mpi eq "yes";
	    }
	    print;
	}
    }

    if($Mpi eq "no"){
	# Make sure Makefile.RULES are up to date, then compile NOMPI
	&create_makefile_rules;
	&shell_command("make NOMPI") if $Mpi eq "no";
    }

    print "Remove executable and make it to link with the (NO)MPI library!\n";
}

##############################################################################

sub set_hdf5_{

    $NewHdf5=$Hdf5 if $Install and not $NewHdf5;

    # Check if HDF5 module is loaded
    if($NewHdf5 eq "yes" and not `which $H5pfc`){
        print "Warning: $H5pfc is not in path. Load parallel hdf5 module!/\n";
        return;
    }

    if($NewHdf5 eq "yes" and not `which $H5pcc`){
        print "Warning: $H5pcc is not in path. Load parallel hdf5 module!/\n";
        return;
    }

    # $Hdf5 will be $NewHdf5 after changes
    $Hdf5 = $NewHdf5;

    print "Enabling HDF5 library in $MakefileConf\n" if $Hdf5 eq "yes";
    print "Disabling Hdf5 library in $MakefileConf\n" if $Hdf5 eq "no";
    if(not $DryRun){

	# For the NAG compiler find the HDF5 include directory from h5pfc -show
	my $H5include;
	$H5include = $1 if ($Compiler eq "f95" or $Compiler eq "nagfor") 
	    and `$H5pfc -show` =~ /( \-I\S+)/;

	@ARGV = ($MakefileConf);
	while(<>){
	    if($Hdf5 eq "yes"){
		# Modify linker definition to use h5pfc
		s/^(LINK\.f90\s*=\s*\$\{CUSTOMPATH_\w+\})(.*)/$1$H5pfc \#$2/
		    unless /\#/;

		# For pgf90 the F90 compiler has to be changed too
		s/^(COMPILE\.f90\s*=.*)(pgf90)/$1$H5pfc \#$2/
		    unless /\#/;

		# Change the parallel C++ compiler too
		s/^(COMPILE\.mpicxx\s*=\s*)(.*)/$1$H5pcc \#$2/
		    unless /\#/;

		# Add the h5pfc include directory to the search path for hdf5.mod
		s/\s+$/$H5include\n/ if /^SEARCH\b/ and $H5include
		    and not /$H5include/;
	    }else{
		# Undo the modifications
		s/($H5pfc|$H5pcc) \#//;
		s/$H5include// if $H5include;
	    }
	    print;
	}
    }

    # PGF90 modules include HDF5 info if compiled with h5pfc
    &shell_command("make clean") if not $Install and $Compiler eq 'pgf90';

    my @files = glob("src/*Hdf5_orig.f90 ".
		     "??/*/src/*Hdf5_orig.f90 ".
		     "share/*/src/ModHdf5Utils_orig.f90 ");
    foreach my $file (@files){
	my $outfile = $file;
	$outfile =~ s/_orig//;
	my $infile  = $file;
	$infile =~ s/_orig/_empty/ if $Hdf5 eq "no";
	print "set_hdf5_: cp $infile $outfile\n";
	&shell_command("cp $infile $outfile");
    }

}

##############################################################################

sub set_hypre_{

    # Check if library is present
    if($NewHypre eq "yes" and not -d "util/HYPRE"){
	print "Warning: util/HYPRE is missing. Use cd util; cvs co HYPRE\n";
	return;
    }

    if($NewHypre eq "yes" and not -e "util/HYPRE/lib/libHYPRE.a"){
	$IsStrict = 0;
	&shell_command("cd util/HYPRE; make install");
	$IsStrict = 1;
	if($ErrorCode){
	    print "$ERROR cd util/HYPRE; make install failed with ".
		"error $ErrorCode\n";
	    print "!!! renaming util/HYPRE to util/HYPRE_FAILED !!!\n";
	    shell_command("rm -rf util/HYPRE_FAILED; ",
			  "mv util/HYPRE util/HYPRE_FAILED");
	    return;
	}
    }

    # $Hypre will be $NewHypre after changes
    $Hypre = $NewHypre;

    print "Enabling HYPRE library in $MakefileConf\n" if $Hypre eq "yes";
    print "Disabling HYPRE library in $MakefileConf\n" if $Hypre eq "no";
    if(not $DryRun){
	@ARGV = ($MakefileConf);
	while(<>){
	    # Add/remove HYPRE related definitions after MPILIB
	    $_ .= $HypreDefinition if $Hypre eq "yes" and /-lNOMPI/;
	    $_ = "" if $Hypre eq "no" and /HYPRE/i;
	    print;
	}
    }

    my @files = glob("src/*Hypre_orig.f90 ??/*/src/*Hypre_orig.f90");
    foreach my $file (@files){
	my $outfile = $file;
	$outfile =~ s/_orig//;
	my $infile  = $file;
	$infile =~ s/_orig/_empty/ if $Hypre eq "no";
	print "set_hypre_: cp $infile $outfile\n";
	&shell_command("cp $infile $outfile");
    }

}

##############################################################################

sub set_fishpak_{

    # Check if library is present 
    if($NewFishpak eq "yes" and not -d "util/FISHPAK"){
        print "Warning: util/FISHPAK is missing. Use cd util; cvs co FISHPAK\n";
        return;
    }

    if($NewFishpak eq "yes" and not -e "util/FISHPAK/lib/libFISHPAK.a"){
        $IsStrict = 0;
        &shell_command("cd util/FISHPAK; make install");
        $IsStrict = 1;
        if($ErrorCode){
            print "$ERROR cd util/FISHPAK; make install failed with ".
                "error $ErrorCode\n";
            print "!!! renaming util/FISHPAK to util/FISHPAK_FAILED !!!\n";
            shell_command("rm -rf util/FISHPAK_FAILED; ",
                          "mv util/FISHPAK util/FISHPAK_FAILED");
            return;
        }
    }

    # $Fishpak will be $NewFishpak after changes                                     
    $Fishpak = $NewFishpak;

    print "Enabling FISHPAK library in $MakefileConf\n" if $Fishpak eq "yes";
    print "Disabling FISHPAK library in $MakefileConf\n" if $Fishpak eq "no";
    if(not $DryRun){
        @ARGV = ($MakefileConf);
        while(<>){
            # Add/remove Fishpak related definitions after MPILIB                                          
            $_ .= $FishpakDefinition if $Fishpak eq "yes" and /-lNOMPI/;
            $_ = "" if $Fishpak eq "no" and /FISHPAK/i;
            print;
        }
    }

}

############################################################################## 

sub set_spice_{

    # By default switch off SPICE for installation
    $NewSpice="no" if $Install and not $NewSpice;

    # Check if library is present, but how?

    my $EnableSpice = ($Spice eq "no" and $NewSpice ne "no");
    my $RemoveSpice = ($Spice ne "no" and $NewSpice eq "no");
    my $ChangeSpice = ($Spice ne "no" and $NewSpice ne "no");

    print "Enabling SPICE library in $MakefileConf\n"  if $EnableSpice;
    print "Disabling SPICE library in $MakefileConf\n" if $RemoveSpice;
    print "Changing SPICE library in $MakefileConf\n"  if $ChangeSpice;

    # $Spice will be $NewSpice after changes
    $Spice = $NewSpice;

    if(not $DryRun){
	@ARGV = ($MakefileConf);
	while(<>){
	    # Remove SPICE related definitions
	    $_ = "" if $RemoveSpice and /SPICE/i;

	    # Add SPICE related definitions after MPILIB
	    $_ .= "# SPICE library\nSPICELIB = $Spice\n" 
		if $EnableSpice and /-lNOMPI/;
	    
	    # Change SPICE related definition
	    s/(SPICELIB\s*=\s*)\S*/$1$Spice/ if $ChangeSpice;

	    print;
	}
    }

    if(not $ChangeSpice){
	# Copy the _orig or _empty versions of the Spice files
	my @files = glob("src/*Spice_orig.f90 */*/src/*Spice_orig.f90");
	foreach my $file (@files){
	    my $outfile = $file;
	    $outfile =~ s/_orig//;
	    my $infile  = $file;
	    $infile =~ s/_orig/_empty/ if $Spice eq "no";
	    print "set_spice_: cp $infile $outfile\n";
	    &shell_command("cp $infile $outfile");
	}
    }
}

##############################################################################

sub set_optimization_{

    # Set the optimization flags in $MakefileConf
    $Optimize = $NewOptimize;

    my $Level=$Optimize; $Level =~ s/-O//;
    print "Setting maximum optimization flag to $Optimize in $MakefileConf\n";
    if(not $DryRun){
	@ARGV = ($MakefileConf);
	while(<>){
	    if (/^\s*OPT([0-5])\s*=\s*/){
		if($1 > $Level){
		    $_ = "OPT$1 = -O$Level\n";
		}else{
		    $_ = "OPT$1 = -O$1\n";
		}
	    }
	    print;
	}
    }
}

##############################################################################

sub create_makefile_rules{

    my @InFile = glob("src*/$MakefileRules.all */*/src*/$MakefileRules.all");

    return unless @InFile;

    # Hash for general configuration settings
    my %Settings = (OS         => $OS, 
		    Compiler   => $Compiler,
		    Mpi        => $Mpi,
		    Debug      => $Debug,
		    Hdf5       => $Hdf5,
		    Machine    => $Machine,
		    Precision  => $Precision);

    # Add settings from the caller Config.pl script
    my $Settings = shift;
    $Settings{$1}=$2 while $Settings =~ s/(\w+)\s*=\s*([^,\n]+)//;

    # Create Makefile.RULES from Makefile.RULES.all in all src* directories
    my $InFile;
    foreach $InFile (@InFile){

	$InFile =~ /([\w\.\/]*)\//;
	my $SrcDir = $1;

	# Open Makefile.RULES.all for reading
	open(INFILE, $InFile) or die "$ERROR_: could not open $InFile\n";

	# Open Makefile.RULES for writing
	my $OutFile = "$SrcDir/$MakefileRules";
	open(OUTFILE, ">$OutFile") or die "$ERROR_: could not open $OutFile\n";

	# Evaluate conditional rules in Makefile.RULES.all
	my $Condition;
	while(<INFILE>){
	    next if /^\#/ or /^\s/; 
	    $Condition = $_;

	    # Replace $xxx variables with their actually set values
	    my $key;
	    foreach $key (keys %Settings){
		$Condition =~ s/\$$key/"$Settings{$key}"/g;
	    }

	    # Skip rules unless condition is true
	    next unless eval($Condition);

	    # Create compilation rules
	    my $Rule;
	    while($Rule = <INFILE>){
		last unless $Rule =~ s/^\t//;

		# Find source file name in the rule
		$Rule =~ /([\w\.]+)\s*$/;

		my $SrcFile = $1;
		my $ObjectFile = $SrcFile; $ObjectFile =~ s/\.\w+$/.o/;

		# Replace ${LINK.f90} with the $MpiCompiler
		$Rule =~ s/\$\{LINK.f90\}/$MpiCompiler/;

		print OUTFILE "$ObjectFile: $SrcFile\n\t$Rule\n";
	    }
	}
	close INFILE;
	close OUTFILE;
    }
}

##############################################################################

sub link_swmf_data{

    if($Code eq "BATSRUS" and not -d "dataCRASH"){
	my $CrashData = "";
	foreach ("$ENV{HOME}/CRASH_data", "/csem1/CRASH_data"){
	    next unless -d $_;
	    $CrashData = $_;
	    last;
	}
	&shell_command("ln -s $CrashData dataCRASH") if $CrashData;
    }

    return if -d "data";
    my $SwmfData = "";
    foreach ("$DIR/SWMF_data", "$DIR/../../SWMF_data", 
	     "$ENV{HOME}/SWMF_data", "/csem1/SWMF_data"){
	next unless -d $_;
	$SwmfData = $_;
	last;
    }
    my $DataDir;
    if($Code eq "SWMF"){
	$DataDir = "$SwmfData/data";
    }else{
	$DataDir = "$SwmfData/$Component/$Code/data";
    }
    &shell_command("ln -s $DataDir data") if -d $DataDir;

}

##############################################################################

sub shell_command{

    my $command = join(' ',@_);
    print "$command\n" if $Verbose;

    return if $DryRun;

    $ErrorCode = system($command) / 256; 

    die "$ERROR Could not execute command=$command: code = $ErrorCode\n"
	if $ErrorCode and $IsStrict;
}

##############################################################################
#BOP
#!QUOTE: \subsection{Installation and Configuration with Config.pl}
#!ROUTINE: Config.pl - (un)installation and configuration of SWMF/components
#!DESCRIPTION:
# The Config.pl provides a single uniform interface towards 
# installation, configuration and uninstallation for the SWMF and its
# components.
#
#!REVISION HISTORY:
# 12/16/2006 G. Toth - initial version based on SetSWMF.pl
#EOP
sub print_help_{

    print 
#BOC
"Config.pl can be used for installing and setting various options for SWMF
or its components. The core of the script is in share/Scripts/Config.pl,
and this is used by the Config.pl scripts in the main SWMF and component 
directories. This help message starts with the options/features/examples 
of the core script, and continues with the additional features (if any)
of the SWMF/component script (starting with the text 'Additional ...').

This script edits the appropriate Makefile-s, copies files and executes 
shell commands. The script can also show the current settings.

Usage: Config.pl [-help] [-verbose] [-dryrun] [-show] [-compiler] 
                 [-install[=s|=c] [-compiler=FC[,CC] [-nompi]]
                 [-uninstall]
                 [-single|-double] [-debug|-nodebug] [-mpi|-nompi]
                 [-hdf5|-nohdf5] [-hypre|-nohypre] [-spice=SPICELIB|-nospice]
                 [-O0|-O1|-O2|-O3|-O4|-O5]

If called without arguments, the current settings are shown.

Information:

-h  -help       show help message.
-dryrun         dry run (do not modify anything, just show actions).
-show           show current settings.
-verbose        show verbose information.

(Un/Re)installation:

-uninstall      uninstall code (make distclean)

-install=c      (re)install code as an SWMF component (c)
-install=s      (re)install code as a stand-alone (s) code
-install        install code as a stand-alone if it is not yet installed,
                or reinstall the same way as it was installed originally:
                (re)creates Makefile.conf, Makefile.def, make install

-compiler       show available compiler choices for this operating system (OS)
-compiler=FC    create Makefile.conf from a non-default F90 compiler FC
                and the default C compiler.
                This setting only works together with the -install flag!
-compiler=FC,CC create Makefile.conf with a non-default F90 compiler FC
                and non-default C compiler CC.
                This setting only works together with the -install flag!
-nompi          install the code on a machine with no MPI library
                This setting only works this way with the -install flag!

Compilation:

-single         set precision to single in Makefile.conf and make clean
-double         set precision to double in Makefile.conf and make clean
-debug          select debug options for the compiler in Makefile.conf
-nodebug        do not use debug options for the compiler in Makefile.conf
-mpi            compile and link with the MPI library for parallel execution
-nompi          compile and link with the NOMPI library for serial execution
-hdf5           compile and link with HDF5 library for HDF5 plot output
-nohdf5         do not compile with HDF5 library
-hypre          link with HYPRE library for linear solver
-nohypre        do not link with HYPRE library
-spice=SPICELIB link with SPICE library SPICELIB for coordinate transforms
-nospice        do not link with SPICE library
-O0             set all optimization levels to -O0
-O1             set optimization levels to at most -O1
-O2             set optimization levels to at most -O2
-O3             set optimization levels to at most -O3
-O4             set optimization levels to at most -O4
-O5             set optimization levels to at most -O5

Examples of use:

Show current settings: 

    Config.pl

Show current settings with more detail: 

    Config.pl -show

Show available compiler choices:

    Config.pl -compiler

Install code with the mpxlf90 Fortran compiler and mpxlc C compiler

    Config.pl -install -compiler=mpxlf90,mpxlc

Install code with the gfortran compiler and no MPI library on the machine

    Config.pl -install -compiler=gfortran -nompi

Use the HDF5 plotting library, the HYPRE linear solver library, and SPICE:

    Config.pl -hdf5 -hypre -spice=/usr/local/lib/spicelib.a

Do not link with the HDF5, HYPRE and SPICE libraries:

    Config.pl -nohdf5 -nohypre -nospice

Set optimization level to -O0, switch on debugging flags and link with NOMPI:

    Config.pl -debug -O0 -nompi

Set optimization level to -03, switch off debugging flags and link with MPI:

    Config.pl -nodebug -O3 -mpi

Uninstall code (if this fails, run Config.pl -install first):

    Config.pl -uninstall"
#EOC
    ,"\n\n";
}

##############################################################################

1;
