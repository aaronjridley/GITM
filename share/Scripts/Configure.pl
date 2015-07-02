#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
# This code is a copyright protected software.
# (c) 2002- University of Michigan
#^CMP FILE CONFIGURE

# List of option files in the CVS distributions. First one is the default.
my @Optionfiles=('Configure.options','Configure.public','Spherical.options');

# Pass switches to reasonably named variables
my $Cfg        = ($c or "CMP");
my $Debug      = $D;
my $Dir        = $d; 
my $Help       = $h;
my $Interactive= $i;
my $Optionfile = $o; 
my $On         = $on;
my $Off        = $off;
my $Show       = $s;
my $Test       = $t;
my $Verbose    = $v;
my $Update     = $u;
my $Exe        = $exe;
my $KeepAll    = $keepall;

use strict;

# Set default destination directory
$Dir = "Build" unless length($Dir)>0;

# Define unique string for directives
my $cfg='^'.$Cfg;

if($Help){&print_help};

$Optionfile=$Optionfiles[0] unless $Optionfile;

# Define a pattern for matching the option files including the one used now
push @Optionfiles, $Optionfile;
my $Optionfiles=join('|',@Optionfiles);

# Read in options into a %switch hash with values ON or OFF
# Read the dependencies into the %needs hash
my $option;
my %switch;
my %needs;

if(-f $Optionfile){
    &read_optionfile;
}else{
    die "ERROR: Could not find option file $Optionfile !\n";
}

# Read the -on switch
foreach $option (split(/,/,$On)){
    $switch{uc($option)}="ON";
}

# Read the -off switch
foreach $option (split(/,/,$Off)){
    $switch{uc($option)}="OFF";
}

# If -keepall and/or -exe is used set CONFIGURE to ON
$switch{"CONFIGURE"} = "ON" if $KeepAll or $Exe;

if (not %switch){
    print "ERROR: No options were provided at all!\n\n";
    print "Help: Configure.pl -h\n";
    exit 1;
}

#  Check relationships between options
foreach $option (keys %needs){
    &check_option($option);
    my $neededoption;
    foreach $neededoption (split(',',$needs{$option})){
	&check_option($neededoption);
	if($Debug){print "checking $option needing $neededoption\n";}
	if($switch{$option} eq "ON" and $switch{$neededoption} eq "OFF"){
	    die "ERROR: $option can be ON only if $neededoption is ON !\n";
	}
    }
}

# Show options if required
if($Show){
    foreach $option (sort keys %switch){
	if($needs{$option}){
	    printf "%-24s %-3s NEEDS %s\n",$option,$switch{$option},
	           $needs{$option};
	}else{
	    printf "%-24s %s\n",$option,$switch{$option};
	}
    }
    exit 0;
}

# Define some global variables for the configuration stage

my $infile;   # the name of the original file
my $outfile;  # the name of the configured file 
my $nline;    # the line number in the original file
my $line;     # the actual line read
my $before;   # the part of the line before the directive
my $linekept; # the part of the line kept after configuration
my %begin;    # $begin{"OPTION"} gives the line number of $Cfg ... OPTION BEGIN
my $option_skip; # the name of the option for which a $Cfg END OPTION 
                 # will close the skipped section
my $output;   # the content of the configured file 

# If CONFIGURE option is ON, make a pattern $KeptOptions matching
# the options which remain in the code.

my $KeepConfig = $switch{CONFIGURE} eq "ON";
my $KeptOptions='nothing_at_all';

# Make sure that $Exe contains an unmatchable pattern if not given

$Exe = ( uc($Exe) or 'nothing_at_all' );

if($KeepConfig){
    my %KeepConfig; # hash for collecting kept options

    if($KeepAll){
	# Set all options to be saved
	foreach $option (keys %switch){$KeepConfig{$option}=1};
    }else{
	# Check which options remain in the configured option file
	$infile = $Optionfile;
	$outfile="$Optionfile.tmp";
	&process_file($infile,$outfile);
	open(OPTION,$outfile) or 
	    die "Cannot open temporary option file $outfile";
	while(<OPTION>){
	    # Ignore empty lines and lines starting with # 
	    next if /^\s*$/ or /^\s*\#/;
	    # Read options of the form "OPTIONNAME    (ON|OFF)"
	    if(/^\s*(\w+)\s+(on|off)\b/i){
		# Save the options in the configured option file into the 
		# KeepConfig hash since these can be further processed
		$option=uc($1);
		$KeepConfig{$option}=1;
	    }else{
		die "ERROR: wrong option definition ".
		    "in configured option file $outfile at line $.:\n$_";
	    }
	}
	close(OPTION);
	print "Configured temporary option file $outfile\n" if $Debug;
	unlink $outfile unless $Debug;
    }
    # Remove options listed in the -exe switch from %KeepConfig
    foreach $option (split /,/,$Exe){
	delete $KeepConfig{uc($option)};
    }
    # Make a pattern for the options to be kept intact
    $KeptOptions = (join('|',sort keys %KeepConfig) or $KeptOptions);

    print "Recursively configurable build:\n";
    print "Options kept: $KeptOptions\n";

    # Make a pattern from $Exe for option line removal
    $Exe =~ s/,/|/g;
}

# Add _FALSE_ to the %switch hash with value OFF
$switch{"_FALSE_"} = "OFF";

# If Copyright is on, we need to insert copyright notices
# Create a hash with the short version of copyright messages 
# indexed by the copyright holder
my $InsertCopyright = $switch{COPYRIGHT} eq "ON";
my $Holder;
my %Copyright;
if($InsertCopyright){
    # Read copyright files from Copyrights/*.short 
    # into an associative array %Copyright indexed by the * part
    my $CopyrightFile;
    foreach $CopyrightFile (glob("Copyrights/*.short")){
	open(FILE,$CopyrightFile) or die "Could not open $CopyrightFile!\n";
	# Extract the holder of the copyright from the file name
	if($CopyrightFile=~/Copyrights\/(\w+)\.short/){
	    $Holder=uc($1);
	}else{
	    die "Could not match filename $CopyrightFile!\n";
	}
	# read text from Copyright/$Holder.short into $Copyright{$Holder}
	while(<FILE>){
	    next if /^#/; # skip lines starting with "#"
	    $Copyright{$Holder}.="# $_";
	}
	close(FILE);
    }
}

# Interactive usage
if($Interactive){
    &process_file("__STDIN__","__STDOUT__");
    exit;
}

# A list of files is given in the argument
if(@ARGV){
    my $infile;
    foreach $infile (@ARGV){

	# Only text files can be processed
	next unless -T $infile;

        # Omit *~ files
        next if $infile =~ /\~$/;

	&process_file($infile,"$Dir/$infile");
    }
    exit;
}

# Process the whole directory tree
if(not -d $Dir){
    print "mkdir $Dir\n" if $Verbose;
    mkdir $Dir,0777 or die "ERROR cannot make $Dir";
}

&process_dir(".");

if($switch{MAKEPDF}      eq "ON" or
   $switch{MAKEHTML}     eq "ON" or
   $switch{REMOVEDOCTEX} eq "ON"){

    chdir $Dir or die "Could not cd $Dir\n";

    # Try installation if SetSWMF.pl is in the configured directory
    if(-f "Config.pl"){
	my $result=system('./Config.pl -install');
	warn "./Config.pl -install: $result\n" if $result;
    }

    # Process documentation if required
    my $texdir;
    TEXDIR: foreach $texdir ("Doc/Tex","doc/Tex","missing"){

	die "Could not find LaTex directory $Dir/[Dd]oc/Tex\n"
	    if $texdir eq "missing";

	# Try to enter the LaTex directory
	chdir $texdir or next TEXDIR;

	# Make HTML manual if required and possible
	if($switch{MAKEHTML} eq "ON"){
	    my $result=system("make HTML");
	    die "Could not make HTML manuals\n" if $result;
	}

        # Make PDF manual if required
        if($switch{MAKEPDF} eq "ON"){
	    my $result=system("make PDF");
	    die "Could not make PDF manuals\n" if $result;
	}

        # Remove $texdir directory if required
	chdir "../.." or die "Could not cd ../..\n";
	if($switch{REMOVEDOCTEX} eq "ON"){
	    print "removing $texdir\n";### if $Verbose;
	    my $result=system("rm -rf $texdir");
	    warn "Could not rm -rf Doc/Tex\n" if $result;
	}

	last TEXDIR;
    }
}

exit 0;

############################################################################
sub process_dir{

    my $dir = $_[0];

    -d $dir or die "ERROR: $dir does not seem to be a directory\n";

    # Check if directory contains a CVS directory or CVS_Entries file
    return unless -d "$dir/CVS" or -f "$dir/CVS_Entries";

    # Check if there is a $Cfg file. 
    # If the $Cfg file contains OPTION     then skipvalue is OFF
    # if the $Cfg file contains NOT OPTION then skipvalue is ON.
    my $skipvalue;
    $infile = "$dir/$Cfg";
    if (-f $infile) {
	# Read OPTION or NOT OPTION from the file
	open(CFGDIRFILE,"$dir/$Cfg") or die "ERROR: cannot open $infile\n";
	$option=uc(<CFGDIRFILE>); chomp $option;
	close(CFGDIRFILE);
	my $skipvalue = ($option =~ s/NOT +//) ? "ON" : "OFF";
	$nline=1;                 # set this for the error message
	&check_option($option);   # in check_options
	# Check if the directory should be skipped.
	# The directory is not skipped if the option is kept configurable
	if($switch{$option} eq $skipvalue and $option !~ /$KeptOptions/i){
            printf "%-45s %-25s %s\n","Skipping directory $dir",
	           " because $option"," is $skipvalue"
		       unless $Update;
	    return;
	}
    }

    if(not -d "$Dir/$dir"){
	print "mkdir $Dir/$dir\n" if $Verbose;
	mkdir "$Dir/$dir",0777 or die "ERROR cannot make $Dir/$dir";
    }

    # Read entries listed in CVS directory or in the CVS_Entries file
    my @entries;
    my $entryfile;
    foreach $entryfile ("$dir/CVS_Entries",
			"$dir/CVS/Entries",
			"$dir/CVS/Entries.Log"){
	next unless -f $entryfile;
	open(ENTRIES,$entryfile) or die "ERROR: cannot open $entryfile\n";
	my $line;
	while($line=<ENTRIES>){
	    # Extract entry names between /.../ delimiters:
	    push (@entries, $1) if $line=~ m|/([^/]+)/|	;
	}
	close(ENTRIES);
    }
    print "$dir ENTRIES:",join(',',@entries),"\n" if $Debug;

    # Process each file listed as an entry
    my @files;
    my $file;
    ENTRY: foreach $file (@entries){

	#skip $Cfg file itself
	if($file eq "$Cfg"){
	    # The $Cfg file should not be kept if the build is not configurable
	    next if not $KeepConfig;
	    # The $Cfg file should not be kept if the option is not kept
	    open(CFGDIRFILE,"$dir/$Cfg");
	    $option=uc(<CFGDIRFILE>); chomp $option;
	    close(CFGDIRFILE);
	    next if $option !~ /$KeptOptions/i;
	}

        my $infile="$dir/$file";       # recursive local variable
	my $outfile="$Dir/$dir/$file"; # recursive local variable
	if(-d $infile){
	    &process_dir($infile);
	}elsif(not -f $infile){
	    warn "WARNING: entry $infile is missing!\n" unless $Update;
	    next ENTRY;
	}else{
	    if($Update and -f $outfile){
		my $mtime_in =(stat($infile ))[9];
		my $mtime_out=(stat($outfile))[9];
		next ENTRY if $mtime_in < $mtime_out;
	    }
	    # Text files except for ps and eps files are processed
	    if(-T $infile and $file !~ /\.(ps|eps)$/){
		&process_file($infile,$outfile);
	    }else{
		print "copy $infile --> $outfile\n" if $Verbose;
		system("cp $infile $outfile");
	    }
	}
	# Save list of files which were produced
	push(@files,$file) if -e $outfile;
    }

    # Produce a CVS_Entries file in target directory
    if(@files){
	open(ENTRIES,">$Dir/$dir/CVS_Entries");
	foreach $file (@files){
	    print ENTRIES "/$file/\n";
	}
	close(ENTRIES);
    }
}

############################################################################
sub process_file{

    # Configure $infile into $outfile or STDIN into STDOUT.

    ($infile,$outfile) = @_;  # global variables !
    print "infile=$infile, outfile=$outfile\n" if $Debug;
    if($infile ne '__STDIN__'){
	open(IN, $infile) or die "ERROR: cannot open input file $infile!\n";
    }

    # Add a backslash for pattern matching   
    my $CFG="\\$cfg";

    # Reinitialize global variables for every file
    $line="";
    $nline=0;
    $output="";
    $option_skip="";
    undef %begin;

    my $commentchar; # this local variable is used several times
    my $skip;        # true while skipping lines inside a BEGIN ... END section

    # process input line by line
    while( ($infile eq '__STDIN__' and $line=<STDIN>) or $line=<IN>){

	$nline++;
	print "line=$line" if $Debug;

	# Remove options listed after the -exe switch
	# from the configured option files
	if($infile =~ /^(\.\/)?($Optionfiles)$/){
	    if($line =~ s/^($Exe)\s*[\s\w]*/#/i){
	       print "REMOVED: $& from $_" if $Debug;
	   }
	}

	# process line
	if($line !~ /$CFG/oi or $line =~ /$CFG[ A-Z]+($KeptOptions)\b/i){

	    # line without directives or with a kept directive

	    print "NO $Cfg or KEPT OPTION: skip=$skip\n" if $Debug;

	    &output($line) unless $skip;

	}elsif($line =~ /(.?)$CFG +COPYRIGHT +(\w+)/oi){

	    $commentchar=$1;
	    $Holder=uc($2);
	    print "matches $Cfg COPYRIGHT $Holder\n" if $Debug;

	    if($InsertCopyright){
		# Find the Copyright message belonging to $Holder
		my $Copyright;
		if($Copyright=$Copyright{$Holder}){
		    # replace the '#' with the comment character or nothing
		    $Copyright=~s/\#/$commentchar/g;
		    # print the copyright notice
		    &output($Copyright);
		}else{
		    #unknown Holder
		    die "ERROR in $Cfg COPYRIGHT $Holder:\n".
                        "   Wrong name $Holder or".
                        " missing file Copyrights/$Holder.short !\n";
		}
	    }

	}elsif($line =~ /$CFG +FILE +((NOT +)?\w+)/oi){

	    print "matches $Cfg FILE $2\n" if $Debug;

	    $option=uc($1);
	    my $skipvalue = ($option =~ s/NOT +//) ? "ON" : "OFF";
            &check_option($option);
            if($switch{$option} eq $skipvalue){
		# Omit this file completely
		$output="";
		close(IN);
		return;
	    }

	}elsif($line =~ /\s*.?$CFG +END +(\w+)/oi
	    or $line =~ /\s*.?$CFG +IF +(\w+) +END/oi){

	    print "matches $Cfg END skip=$skip\n" if $Debug;
	    $before=$`;
	    $option=uc($1);
	    &end_option($option);
	    &set_line_kept;

	    &output($linekept) if length($linekept)>1 and not $skip;

	    if($option eq $option_skip){
		print "skip set to 0 for $option\n" if $Debug;
		$skip=0;
		$option_skip="";
	    }
	}elsif($line=~/\s*.?$CFG +IF +((NOT +)?\w+) +BEGIN/oi){
	    print "matches $Cfg IF $2 ... BEGIN\n" if $Debug;

	    $before=$`;
	    $option=uc($1);
	    my $skipvalue = ($option =~ s/NOT +//) ? "ON" : "OFF";
	    &begin_option($option);
	    &set_line_kept;

	    if(not $skip){
		if($switch{$option} eq $skipvalue){
		    # Skip lines from this point on
		    $option_skip=$option;
		    $skip=1;
		    print "skip set to 1 for $option\n" if $Debug;
		}elsif(length($linekept)>1){
		    &output($linekept);
		}
	    }

	}elsif($line=~/\s*(.)?$CFG +IF +((NOT +)?\w+)/oi){
	    $before=$`;
	    $commentchar=$1;
	    $option=uc($2);
	    my $not=uc($3);
	    my $keepvalue = ($option =~ s/NOT +//) ? "OFF" : "ON";

	    print "matches $Cfg IF $not...\n" if $Debug;

	    &check_option($option);
	    &set_line_kept;

            # remove leading commentchar if present
	    $linekept =~ s/^(\s*)$commentchar/$1/ if $not and $commentchar;

	    &output($linekept) if not $skip and $switch{$option} eq $keepvalue;

	}elsif($line=~/\s*$CFG +UNCOMMENT +IF +((NOT +)?\w+)/oi){
	    $before=$`;
	    $option=uc($1);
	    my $not=uc($2);

	    print "matches $Cfg UNCOMMENT IF $not ...\n" if $Debug;

	    my $uncommentvalue = ($option =~ s/NOT +//) ? "OFF" : "ON";
	    &check_option($option);
	    &set_line_kept;

	    # remove leading commentchar
	    $linekept =~ s/^(\s*)(<!-- ?|\W)/$1/
		or die "ERROR: $cfg UNCOMMENT IF $not $option ".
		"with missing comment character ".
		"in file $infile at line $nline:\n$line";

	    &output($linekept) if not $skip and 
		$switch{$option} eq $uncommentvalue;

	}else{
	    die "ERROR: Unknown $cfg command ".
		"in file $infile at line $nline:\n$line";
	}
    }
    close(IN);

    foreach $option (sort keys %begin){
	warn "WARNING: $option BEGIN found without matching END ".
	    "in file $infile at line $begin{$option}\n"
	    if $begin{$option};
    }

    if($output =~ /^\s*$/){
	print "File $infile produced no output!\n" 
	    unless $Update or $Interactive;
	return;
    }

    print "processed $infile --> $outfile\n" if $Verbose;

    return if $Test;

    if($outfile eq "__STDOUT__"){
	print $output;
    }else{
	open(OUTPUT,">$outfile") or 
	    die "ERROR: cannot open output file $outfile!\n";
	print OUTPUT $output;
	close(OUTPUT);

	# Set the same permissions as for the infile
	my $stat=(stat($infile ))[2];
	chmod $stat, $outfile;
    }

}
############################################################################
sub read_optionfile{

    # Initialize %switch and %needs associative arrays
    undef %switch;
    undef %needs;

    open(OPTION,$Optionfile) or die "Cannot open option file $Optionfile";

    while(<OPTION>){
	# Ignore empty lines and lines starting with # 
	next if /^\s*$/ or /^\s*\#/;

        # Read options of the form 
        # OPTIONNAME    (ON|OFF)
	if(/^\s*(\w+)\s+(on|off)\b/i){
	    my $option=uc($1);
	    my $switch=uc($2);
	    $switch{$option} = $switch;
	    if(/needs\s+([\w\,]+)/i){
		$needs{$option}=uc($1);
	    }
	}else{
	    die "ERROR: wrong option definition ".
		"in option file $Optionfile at line $.:\n$_";
	}
    }

    close(OPTION);
}
############################################################################
sub check_option{
    my $option = $_[0];
    if( length( $switch{$option} ) == 0 ){
	die "ERROR: unknown option $option ".
	    " in file $infile at line $nline:\n$line";
    }
}
############################################################################
sub begin_option{
    my $option = $_[0];

    # Check if option name is valid
    &check_option($option);

    # Check if option had a BEGIN already
    if( $begin{$option} ){
	die "ERROR: second $cfg $option BEGIN without $cfg END $option ",
	      "in file $infile at line $nline: \n$line\n";
    }

    # Remember that option has begun
    $begin{$option}=$nline;
}
############################################################################
sub end_option{
    my $option = $_[0];

    # Check if option name is valid
    &check_option($option);

    # Check if option has begun
    if( not $begin{$option} ){
	die "ERROR: $cfg END $option found without $cfg $option BEGIN ",
	      "in file $infile at line $nline:\n$line";
    }

    # Close the begin
    $begin{$option}=0;
}
############################################################################
sub set_line_kept{
    # based on the previous matching of the CFG directive 
    # set the part of the line which should be kept after configuration

    $linekept = "$before\n";
    $linekept =~ s/\s*<!--?\s*$/\n/; # Remove XML comment <!--
}
############################################################################
sub output{
    my $data=$_[0];
    if($outfile eq '__STDOUT__'){
	print "$data";
    }else{
	$output .= "$data";
    }
}
############################################################################
#BOP
#!ROUTINE: Configure.pl - create a configured software package based on directives
#!DESCRIPTION:
# This script together with the directives inserted into the source files
# (this can be any ASCII file, not just source code) allows a flexible
# and versatile configuration of a software package. The software package
# is assumed to be under CVS control, and normally only those files are
# configured which belong to the CVS repositry. Adding a special 
# config file with a directory can configure out the whole directory.
# The configuration can be recursive, i.e. the configured software package
# can be further configured.
#
#!REVISION HISTORY:
# 01/24/2002 G. Toth - initial version
#                      many-many extensions, fixes, improvements
#EOP
sub print_help{

    print 
#BOC
"   Purpose: 

Build a new code distribution containing a subset of the code.
The subset is selected according to the switches set to ON and OFF
and the corresponding directives in the files and directories.

   Usage:

Configure.pl [-c=STR] [-D] [-d=DIR] [-exe=option1,option2...] [-h] [-i] 
             [-keepall]
             [-o=optionfile] [-on=option1,option2...] [-off=option3,option4...]
             [-s] [-t] [-v]
             [inputfile1] [inputfile2] ...

      NOTE: all option names, 'on' or 'off' strings and all \^CXX directives
      are capitalized by Configure.pl, so the syntax is case insensitive !

      -c=STR  Sets the string used as the prefix of the directives to ^STR 
              and the name of the directory config files to STR.
              The default value for STR is 'CMP',

      -D      print debug info

      -d=DIR  write output into directory DIR with the same filename as
              the input file(s). The default directory is ./Build

      -h      print this help message

      -i      read from STDIN and flush result to STDOUT

      -keepall keep (and ignore) all directives except for option _FALSE_ 
              and the options listed with the -exe switch.
              The -keepall switch also sets CONFIGURE to ON.

      -exe=option1,option2...

              Execute the directives belonging to the listed options and
              remove the directives from the code and the options from the 
              configured option files. The -exe switch sets CONFIGURE to ON.

      -o=OPTIONFILE

              If the option file OPTIONFILE is not defined with the -o flag, 
              the default Configure.options file is read if available.

              The option file contains the options, the ON/OFF info,
              and the dependencies in the following format

              # some comment
              #
              option1 ON
              option2 OFF NEEDS option1
              option3 ON  NEEDS option1,option6
              ...

              Here 'NEEDS option1' means that option2 can be ON only 
              if option1 is ON, while 'NEEDS option1,option6' means
              that option3 can be ON only if option1 and option6 are ON.
              Empty lines and comments following # are permitted.

              The options in combination with the directives and the
              directory config files can be used to remove or insert some 
              text, to omit files, or to skip whole directories during 
              configuration.

              Special options:

              If CONFIGURE is ON then the directives are kept 
              without any processing for the options which remain 
              configurable. This is based on the directives in the option file 
              (default is Configure.options): if the configured
              option file contains an option it is regarded as configurable.
              When -keepall is used, all options are taken as configurable.
              A subset of the configurable options can be removed with
              the -exe switch if required. 

              If COPYRIGHT is ON, the COPYRIGHT directives are 
              replaced with the appropriate COPYRIGHT message.

              If MAKEPDF is ON, Configure.pl will execute 'make PDF'
              in the configured distribution.

              If MAKEHTML is ON, then HTML documentation is made 
              in the configured distribution. 

              If REMOVEDOCTEX is ON, the Doc/Tex directory will be removed
              from the configured distribution.

              The _FALSE_ options is always OFF and it cannot be kept.
              This option marks unconfigurable features.

      -on=option1,option2,...

              list of options being ON separated with commas. 
              Adds new options or overwrites values read from optionfile.

      -off=option1,option2,...

              list of options being OFF separated with commas. 
              Adds new options or overwrites values read from optionfile.

      -s      show list of options and their values read from the option file 
              and modified by the -on and -off flags but do not process input.
              The list is in the same format as the option file.

      -t      test input for syntax errors but no output is printed.

      -v      verbose (for complete configuration)

      inputfile(s): the files to be preprocessed. 

                    If no inputfile is given, the content of the 
                    current directory is processed recursively into
                    Build or the directory defined by -d=DIR

                    Only files listed in CVS/Entries, CVS/Entries.Log
                    or CVS_Entries are included. The CVS_Entries file
                    is produced by Configure.pl for reference.

                    Binary, .ps, and .eps files are copied 
                    without processing.

                    If a directory config file exists in a directory and c
                    ontains an OPTIONNAME (or NOT OPTIONNAME), the directory 
                    is omitted when that option is OFF (or ON). 
                    Otherwise the directory is processed, but 
                    the config file is only kept if the CONFIGURE option is ON.

   Examples:

      Show list of options in Configure.options:

Configure.pl -s

      Modify option list read from Configure.options and put it into a 
      new option file:

Configure.pl -on=roeflux -off=implicit,lindeflux -s > new.options

      Convert STDIN to STDOUT using Configure.options:

Configure.pl -i

      Convert infile to outfile using new.options:

Configure.pl -i -o=new.options < infile > outfile

      Process files in indir and put processed files into outdir:

Configure.pl -d=outdir indir/*

      Configure BATSRUS and build a complete distribution:

Configure.pl

      Make a configurable build with the implicit scheme:

Configure.pl -on=configure,implicit

      Update all files that changed and be verbose about it:

Configure.pl -on=configure,implicit -u -v

      Make a configurable but purely Cartesian build into CARTESIAN directory:

Configure.pl -keepall -exe=cartesian -on=cartesian -d=CARTESIAN

      Build a complete distribution with manuals but no Doc/Tex directory:

Configure.pl -on=DOC,DOCHTML,MAKEPDF,MAKEHTML,REMOVEDOCTEX


   Syntax for directives:

      There can be only one directive per line, but 
      the IF (NOT) OPTION BEGIN ... END OPTION
      constructs can be arbitrarily nested. 
      The following example shows all the directives with the directive
      prefix set to the STR string:

!^STR COPYRIGHT UM
!^STR FILE OPTION1
!^STR FILE NOT OPTION4
text0
text1 !^STR IF OPTION1
text2 !^STR IF OPTION2
text3 !^STR IF OPTION1 BEGIN
text4
text5 !^STR END OPTION1
text6 !^STR IF OPTION2 BEGIN
text7
text8 !^STR END OPTION2
text9
text10
!text11 ^STR UNCOMMENT IF OPTION1
<!-- text12 ^STR UNCOMMENT IF NOT OPTION2 -->
text13 !^STR IF NOT OPTION3 BEGIN
text14
text15 !^STR END OPTION3
text16 !^STR IF NOT OPTION4 BEGIN
text17
text18 !^STR END OPTION4
text19 !^STR IF OPTION3

      If this is proccessed with 

Configure.pl -on=OPTION1,OPTION3,COPYRIGHT -off=OPTION2,OPTION4

       we obtain:

! Copyright of Univeristy of Michigan (c) 2002
! See the file Copyrights/UM.long for complete copyright information
!
text0
text1
text3
text4
text5
text9
text10
text11
text12
text16
text17
text18
text19

       Note that the comment character in front of ^STR COPYRIGHT is used
       to comment out the text of the COPYRIGHT.

       Also note that the leading comment character 
       (a non-space-or-alphanumeric character) 
       as well as the XML comment string (<--!) 
       in front of text11 and text12, respectively, were removed! 

       If OPTION1 was OFF or OPTION4 was ON, the configured file would be 
       empty and omitted due to the ^STR FILE directives. 

       Recursive configuration is possible if the CONFIGURE option is ON.
       In this case all options which remain in the option file after
       configuration will keep their directives for later configuration.
       Let's assume that we have the following option file 'options':

CONFIGURE                OFF            ^STR IF CONFIGURE
COPYRIGHT                OFF            ^STR IF NOT COPYRIGHT
OPTION1                  ON             ^STR IF OPTION1
OPTION2                  OFF            ^STR IF OPTION2
OPTION3                  ON             ^STR IF OPTION3
OPTION4                  OFF            ^STR IF NOT OPTION4

       If the same input file is processed with

Configure.pl -o=options -on=CONFIGURE

       we obtain:

!^STR COPYRIGHT UM
!^STR FILE OPTION1
!^STR FILE NOT OPTION4
text0
text1 !^STR IF OPTION1
text3 !^STR IF OPTION1 BEGIN
text4
text5 !^STR END OPTION1
text9
text10
!text11 ^STR UNCOMMENT IF OPTION1
text12
text13 !^STR IF NOT OPTION3 BEGIN
text14
text15 !^STR END OPTION3
text16 !^STR IF NOT OPTION4 BEGIN
text17
text18 !^STR END OPTION4
text19 !^STR IF OPTION3

       The output remains configurable in terms of OPTION1, OPTION3, OPTION4,
       and COPYRIGHT. Note that line 11 remains commented out and 
       lines 13 to 15 are kept although OPTION1 and OPTION3 were ON. 
       If we want to keep only the lines that correspond to
       OPTION3 with ON then

Configure.pl -o=options -exe=OPTION3

       and we obtain:

!^STR COPYRIGHT UM
!^STR FILE OPTION1
!^STR FILE NOT OPTION4
text0
text1 !^STR IF OPTION1
text3 !^STR IF OPTION1 BEGIN
text4
text5 !^STR END OPTION1
text9
text10
!text11 ^STR UNCOMMENT IF OPTION1
text12
text16 !^STR IF NOT OPTION4 BEGIN
text17
text18 !^STR END OPTION4
text19

       Note that in this case all the OPTION3 directives were 
       executed and removed. The OPTION3 line will also be removed 
       from the configured option files."
#EOC
    ,"\n\n";
    exit 0;
}
