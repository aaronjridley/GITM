#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
#BOP
#!ROUTINE: Rename.pl - rename many variables/strings in many source files
#
#!DESCRIPTION:
# This script makes variable renaming relatively safe and easy in 
# a large number of source files. It can rename many variables at
# the same time. Case insensitivity is taken care of, sub strings
# are not replaced, and conflicting renaming rules are checked for.
# The renaming rules are stored in an associative array, which
# is an easy to edit Perl code. 
#
#!REVISION HISTORY:
# 07/13/2001 G.Toth gtoth@umich.edu - initial version
# 07/23/2004 G.Toth allow protecting some lines
# 09/23/2005        added -common=XX option
#EOP

$Help=$h; 
$Input=$i;
$List=$l; 
$Check=$c;  
$Rename=$r; 
$Undo=$u; 
$Debug=$d;
$Quiet=$q;
$Warning=$w;
$Common=$common;
$nThread= ($n or 4);

if($Check and $Warning){
    print "Either check (-c) or no warnings (-w) \n\n";
    $Help=1;
}

if($List + $Check + $Rename + $Undo != 1){
    print "Specify exactly one of the option out of  -l, -c, -r, -u !\n\n";
    $Help=1;
}

$Input="RenameList.pl" unless $Input;

if($Help or not ($List or $Check or $Rename or $Undo) ){
    print 
#BOC
'Purpose:

   Replace multiple variable/method/module names in multiple files.

Usage: 

   Rename.pl [options...] file1 file2 ...

Options (specify them separately as  -a -b  and not as  -ab !):

-h             help (this message)

-i=INPUTFILE   input file with the renaming rules. Default is RenameList.pl
               The INPUTFILE contains a single Perl statement, which sets 
               the associative array %newname:

               %newname = (oldvar1 => newvar1,
                           oldvar2 => newvar2,
                           ...
                          );

-l             list names to be changed
-n=NTHREAD     process files in parallel. Default is 4 threads.
-c             check input file and source files, do not rename variables
-r             replace old names with new
-u             undo replacements
-common=XX     replace common/NAME/ blocks with common/XX_NAME/ blocks
               Only works if -r or -u are present.

-d             debug info is printed
-q             quiet run, errors and warnings printed only
-w             warning messages are suppressed

You have to specify exactly one of -l, -c, -r, or -u. 

For replace (-r) and undo (-u) the original file is put into filename~
if any replacements were done. Individual lines can be protected against
renaming by adding a "!do not rename" trailing comment. 

Typical usage:

Rename.pl -c                   #  check the replacement rules
Rename.pl -c -n=1 *.f90        #  check the source files on 1 core
Rename.pl -r -n=8 *.f90        #  do replacements on 8 cores
Rename.pl -r -common=SC *.f90  #  do replacements including common blocks'
#EOC
,"\n";

    exit;
}

require $Input or die "Could not read input file $Input !!!\n";

@oldname=sort keys %newname;

# Define upper case version of old names as the keys of %NEWNAME
foreach $oldname (@oldname){
    $NEWNAME{uc($oldname)}=$newname{$oldname};
}

# Check for multiple occurances of new names and produce inverse lookup
$error=0;
foreach $oldname (@oldname){
    $newname=$newname{$oldname};
    $NEWNAME=uc($newname);
    if($otheroldname=$OLDNAME{$NEWNAME}){
	print "ERROR: both $oldname \& $otheroldname --> $newname !!!\n";
	$error=1;
    }else{
	$oldname{$newname}=$oldname;
        $OLDNAME{$NEWNAME}=$oldname;
    }
}

exit if $error;

# Check for new names that would be renamed next time
foreach $oldname (@oldname){
    $newname=$newname{$oldname};
    $NEWNAME=uc($newname);

    if($newnewname=$NEWNAME{$NEWNAME}){
        $NEWNEWNAME=uc($newnewname);
        # Do not worry if only the case is different
	next if $NEWNEWNAME eq $NEWNAME;

        # Get the old name that should really be renamed into $newnewname 
        $oldnewnewname=$OLDNAME{$NEWNEWNAME};

	if($oldnewnewname eq $newname){
            # Even the cases are the same
	    print "WARNING: $oldname  --> $newname -->  $newnewname !!!\n"
		unless $Warning;
	}else{
	    print "Warning: $oldname --> $newname, ".
		"$oldnewnewname --> $newnewname !!!\n" unless $Warning;
	    $DANGER{$NEWNAME}=1;
	}
    }
}

if($List){
    # list old and new names in alphabetical order 
    foreach $oldname (@oldname){
	print "$oldname => $newname{$oldname}\n";
    }
}

if($Rename or $Check){
    &rename_files(%newname);
}elsif($Undo){
    &rename_files(%oldname);
};

exit;

##########################################################################
sub rename_files{

    %rename=@_;

    @old=sort keys   %rename;
    @new=sort values %rename;

    @Files = @ARGV;
    $nFile = scalar(@Files);

    print "Number of files to process=$nFile\n" unless $Quiet or $nFile<=1;

    $nfile = 0; $countall = 0;
    if($nThread > 1 and $nFile > 1){
	foreach my $iThread (1..$nThread){
	    # parent process does nothing
	    next if fork();
	    for (my $iFile = $iThread-1; $iFile < $nFile; $iFile+=$nThread){
		$File = $Files[$iFile];
		&rename($File);
	    }
	    print "Thread $iThread finished $countall replacement(s) ".
		"in $nfile file(s)\n" if $countall and not $Quiet;
	    exit;
	}
	foreach (1..$nThread){wait};
    }else{
	foreach $File (@Files){
	    &rename($File);
	}
	print "Finished $countall replacement(s) in $nfile file(s)\n"
	    if $countall and not $Quiet;
    }
}

##########################################################################
sub rename{

    $file = shift;

    if(not open(FILE,$file)){
	print "Error opening file $file !!!\n";
	return;
    }
    read(FILE,$text,-s $file);
    close(FILE);
    print "old text=\n$text\n" if $Debug;

    # Protect lines with containing the '!DO NOT RENAME' string and 
    # case(' or case(" statements.
    $icase=0; @case=();
    while($text =~ s/^(.*\!\s*do\ not\ rename.*|
			   \s*case\s*\(\s*[\'\"].*)/_\[CASE$icase\]_/imx){
	print " replacing case $icase\n" if $Debug;
	$case[$icase++]=$1;
    }

    foreach $oldname (@old){
	$newname=$rename{$oldname};
	next if lc($newname) eq lc($oldname); # Only capitalization changes
	if($text=~/\b$newname\b/i){
	    print "  warning: variable $newname occurs in file $file !\n"
		unless $Warning;
	}
    }

    # Finished with this file if check only
    return if $Check;

    print "Renaming variables...\n" if $Debug;

    # Replace old names with tokens of the form _[NUMBER]_
    $count=0;

    # Replace common block names if required
    if($Common){
	if($Undo){
	    # Undo: common/XX_.../ --> common/.../
	    $count+=($text =~ s[^(\s*common\s*/)$Common\_(\w+/)][$1$2]mgi);
	}else{
	    # Replace: common/.../ --> common/XX_.../
	    $count +=
		($text =~ s[^(\s*common\s*/)(\w+/)][$1$Common\_$2]mgi);
	    # Avoid double change: common/XX_XX_ --> common/XX_
	    $count -=
		($text =~ 
		 s[^(\s*common\s*/)$Common\_$Common\_][$1$Common\_]mgi);
	}
    }

    for($i=0;$i<=$#old;$i++){
	$oldname=$old[$i];
	if($Undo or $DANGER{uc($oldname)}){
	    # Replace only with strict case agreement
	    $count += ($text=~s/\b$oldname\b/_\[$i\]_/g);
	}else{
	    # Replace with relaxed case checking
	    $count += ($text=~s/\b$oldname\b/_\[$i\]_/ig);
	}
    }
    if($count==0){
	print "No variables to rename in file $file\n" unless $Quiet;
	return;
    }
    print "tok text=\n$text\n" if $Debug;
    # Replace tokens with new names
    for($i=0;$i<=$#old;$i++){
	$text=~s/_\[$i\]_/$rename{$old[$i]}/g;
    }
    print "New text=\n$text\n" if $Debug;

    # Put back the case(' and case(" lines
    for($icase=0; $icase<=$#case; $icase++){
	$text =~ s/_\[CASE$icase\]_/$case[$icase]/i;
    }

    # Replace the file with the modified text
    rename $file,"$file~";
    open(FILE,">$file");
    print FILE $text;
    close(FILE);

    print "Finished $count replacements in file $file\n" unless $Quiet;
    $nfile++;
    $countall += $count;

}
