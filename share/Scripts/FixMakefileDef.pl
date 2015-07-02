#!/usr/bin/perl -pi
#BOP
#!ROUTINE: FixMakefileDef.pl
#!DESCRIPTION:
# Fix the definition of OS and SWMF\_ROOT in a file (typically Makefile.def).
# This script can be called from the main Makefile to edit Makefile.def
# when the directory was renamed, moved or transferred to another platform.
#
#!REVISION HISTORY:
# 03/22/2004 G. Toth - initial version
#EOP
#BOC
s/^OS *=.*\n/"OS = ".`uname`/e;
s/^DIR *=.*\n/"DIR = ".`\/bin\/pwd`/e;
#EOC
