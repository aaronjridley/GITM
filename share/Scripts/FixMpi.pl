#!/usr/bin/perl -pi -s
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

#BOP
#!ROUTINE: FixMpi.pl
#!DESCRIPTION:
# It is customary to include the MPI header file into a Fortran 90 program.
# But using a module is a nicer solution. This script replaces the include
# statements with 'use ModMpi' or 'use ModMpiOrig' statements. 
# Special care is taken for placing the use statement before the
# 'implicit none' statement.
# Usage:
#
#\begin{verbatim}
#    FixMpi.pl [-include=INCLUDEFILE] [-module=MODULE] *.f90 *.F90
#\end{verbatim}
#  The default MODULE is 'ModMpi', the default INCLUDEFILE is 'mpif90.h'.
#  Here are a few examples of the effect of running this code:
#\begin{verbatim}
#   implicit none
#   ! at most 2 lines here
#   include '../Common/mpif90.h'
#\end{verbatim}
# is replaced with
#\begin{verbatim}
#   use ModMpi
#   implicit none
#   ! at most 2 lines here
#\end{verbatim}
# and any remaining 
#\begin{verbatim}
#   include '../Common/mpif90.h' 
#\end{verbatim}
# is replaced with
#\begin{verbatim}
#   use ModMpi
#\end{verbatim}
#!REVISION HISTORY:
# 09/01/2003 G. Toth - initial version
#EOP
#BOC
BEGIN{$module = "ModMpi" unless $module; $include ="mpif90.h" unless $include}

if(/^ *implicit +none/i){
    $implnone=$_;
    for $i (1..3){
	$_=<>;
	if(/^( *)include +['"].*$include["']/i){
	    $_="$1use $module\n$implnone";
	    last;
	}else{
	    $implnone.=$_;
	    $_=$implnone;
	}
    }
}

# Replace the remaining strings
s/include +['"].*$include["']/use $module/gi;
#EOC
