#!/usr/bin/perl -pi
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
#BOP
#!ROUTINE: StopToConStop.pl
#!DESCRIPTION:
# This script helps eliminating stop statements from the 
# F90 and F77 source code. Avoiding STOP statements is a 
# requirement in SWMF.
# It replaces STOP statements with call CON\_stop(string) in F90 and F77 codes.
# The string will start with the 'ERROR in FILENAME' message.
# In case the script is called from the main SWMF directory,
# the filename will contain the path as well. This is the
# recommended usaga.
#
# Usage: 
#\begin{verbatim}
#     StopToConStop.pl [file1, file2, ...]
#\end{verbatim}
# Example:
#\begin{verbatim}
#     share/Scripts/StopToConStop.pl IM/RCM/src/RCM_advec.f90
#\end{verbatim}
#
# Limitations: 
#     "command; stop" and 
#     stop statements with continuation lines are not recognized.
#
#!REVISION HISTORY:
# 08/18/2004: G. Toth - initial version
# 01/05/2006;           added code to handle if(condition) stop...
#EOP
#BOC

# remove initial if statement if any
my $if;
$if=$1 if s/^(\s*if\s*\(.*\))\s*stop/ stop/i;

# Plain stop statement with an optional semi-colon
s/^(\s*)stop\s*;?\s*$/$1call CON_stop('ERROR in $ARGV')\n/i;

# Plain stop statement with trailing remark
if(/^(\s*)stop\s;?\s*!(.*)$/i){
    $remark = $2; 
    $remark =~ s/[\s\n]+$//;  # get rid of trailing spaces
    $remark =~ s/\'//g;       # get rid of quotation marks
    $_="$1call CON_stop('ERROR in $ARGV:$remark')\n";
}

# Stop statement with single quoted string
s/^(\s*)stop\s*'([^']*)'/$1call CON_stop('ERROR in $ARGV:$2')/i;

# Stop statement with double quoted string
s/^(\s*)stop\s*"([^"]*)"/$1call CON_stop("ERROR in $ARGV:$2")/i;

# Put back the initial if(...) if present
$_ = "$if$_" if $if;

# Shorten F90 line if necessary
s|^(\s*)(.*call CON_stop\()(["'])(ERROR in $ARGV:)|$1$2$3$4$3// &\n    $1$3|
    if length($_) > 79 and $ARGV =~ /\.(f90|f95)$/i;

# Shorten F77 line if necessary
s{^      (\s*)(.*call CON_stop\()(["'])(ERROR in $ARGV:)}
 {      $1$2$3$4$3//\n     &    $1$3}
    if length($_) > 72 and $ARGV =~ /\.(f|f77|for)$/i;

#EOC
