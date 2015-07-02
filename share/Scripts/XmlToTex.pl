#!/usr/bin/perl -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
#!QUOTE: \clearpage
#BOP
#!QUOTE: \subsection{Automated Documentation}
#!ROUTINE: XmlToTex.pl - generate Latex documentation from XML definitions of input parameters
#!DESCRIPTION:
# Create Latex documentation based on the XML definitions of the 
# input parameters typically found in the PARAM.XML files.
# This script allows to store the parameter descriptions in a single XML file
# which is suitable for automated parameter checking and GUI generation,
# yet provide the same information in a well formatted and printable manual
# as well. The specific format of the PARAM.XML files is described by the
# share/Scripts/CheckParam.pl script and the manual.
#
#!REVISION HISTORY:
# 03/22/2004 G.Toth - initial version
# 08/13/2004          preserve _{.. for subscripts in Latex.
# 01/25/2005          added generation of index terms for commands.
# 01/27/2005          Better handling of verbatim environment.
#EOP
if($h|$help|$H|$Help){
    print 
#BOC
"Purpose:

   Convert XML description of input commands into LaTex description.

Usage:

   XmlToTex [-h] [infile] [> outfile]

-h       This help message

infile   Input file. Default is reading from STDIN.

outfile  Output file. Default is writing to STDOUT.

Example: 

  share/Scripts/XmlToTex.pl Param/PARAM.XML > Param/PARAM.xmltex"
#EOC
    ,"\n\n";
    exit 0;
}

use strict;
my ($verbatim, $Verbatim, $comment, $rule, $start);

while(<>){

    # Skip lines with lots of exclamation marks used for human reading
    next if /!!!!!!!!/;

    # commandList --> section
    if(/^\s*<commandList\s+name=[\'\"]([^\'\"]+)/){
	$_="\\section\{Input Commands for the $1\}\n\n";
        $start = 1;
    }

    # commandgroup --> subsection
    if(/^\s*<commandgroup\s+name=[\'\"]([^\'\"]+)/){
	$_="\\subsection\{\u\L$1\}\n\n";
    }

    # skip things before the commandList
    next unless $start; 

    # command --> subsubsection with an index term
    if(/^\s*<command/){
	/name=[\'\"]([^\'\"]+)/;
	my $command   = "\#$1";
	my @path = split('/',$ARGV);
	my $comp = 'CON'; 
	$comp="$path[$#path-2]/$path[$#path-1]"  # UA/GITM2
	    if $path[$#path-1] ne 'Param';       # Param/PARAM.XML is for CON

	$comp =~ s/GM\/BATSRUS/EE,GM,SC,IH,OH\/BATSRUS/; # GM --> EE,GM,SC,IH,OH

	my $index = "$comp\!$command";           # Form the index term
	$_="\\subsubsection\{$command command\}\\index\{$index\}\n\n";
    }

    # \begin{verbatim} starts forced verbatim mode
    $Verbatim = 1 if /\\begin\{verbatim\}/;

    # \end{verbatim} finishes forced verbatim mode
    $Verbatim = 0 if /\\end\{verbatim\}/;

    # #COMMAND or #COMMAND ID --> verbatim unless in forced Verbatim mode
    if(/^\#[A-Z0-9_]+( [A-Z][A-Z])?\s*$/ and not $Verbatim){
	$_='\begin'.'{verbatim}'."\n$_";
	$verbatim = 1;
    }

    # verbatim part ends with an empty line unless in forced Verbatim mode
    if($verbatim and /^\s*$/ and not $Verbatim){
	$_='\end'.'{verbatim}'."\n";
	$verbatim = 0;
    }

    # Remove one line XML comments
    s/<!--(.*?)-->//;

    # Beginning of multiline XML comment
    $comment = 1 if s/<!--.*//;

    # End of XML comment
    $comment = 0 if s/.*-->//;

    # Beginning of rule
    $rule=1 if /\s*<rule/;

    # End of rule
    $rule=0 if s/.*<\/rule>//;

    # Skip XML lines, XML comments and rules
    next if /[<>]/ or $comment or $rule;

    # Replace special XML characters with Tex
    if($verbatim or $Verbatim){
	s/\&lt;/</g;
	s/\&gt;/>/g;
	s/\&amp;/\&/g;
    }else{
	s/\&lt;/\$<\$/g;
	s/\&gt;/\$>\$/g;
	s/\&amp;/\\&/g;
    }

    # Replace TAB characters with spaces
    1 while s/\t+/' ' x (length($&) * 8 - length($`) % 8)/e;

    if(not ($verbatim or $Verbatim)){
	# Put \ before special Tex character #
	s/\#/\\\#/g;
	# Put \ before special Tex character _ but not before _{
	# This allows the use of _{1} subscript in math mode
	s/(\_[^\{])/\\$1/g;
	# Do not put two \ before # and _
	s/\\\\([\#_])/\\$1/g;
    }

    # Remove ! from the beginning of the lines
    s/^! ?//;

    # Print out the line
    print;
}

# Finish with closing the page
print "\\clearpage\n";
