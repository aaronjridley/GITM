#BOP
#!QUOTE: \newpage \subsection{Scripts and Makefiles} 
#!QUOTE: Here are the scripts and Makefile-s
#!ROUTINE: Makefile - an example how to produce Protex documentation

#!DESCRIPTION:
# This Makefile is a template for creating documentation from 
# Protex documented F90 source code and scripts/Makefile-s.
#
# Describe here the script or Makefile with Latex syntax
# Note that in scripts/Makefile the Protex tags require a comment string 
# and an exclamation mark except for \#BOP, \#EOP, \#BOC and \#EOC.
# 
# It is a good idea to include the help message of the script 
# or the help target of the Makefile between the \#BOC and \#EOC,
# so the manual contains the same information as the interactive help.
#
# The !QUOTE:$\backslash$newpage $\backslash$subsection\{..\} tag is 
# only required in the first script in
# the section.

#!USES:
include ../../Makefile.def
include ../../Makefile.conf
include ../build/Makefile.doc

#!REVISION HISTORY:
# 04/28/2004 G.Toth - initial version
# 03/04/2005 G.Toth - include ../build/Makefile.doc
#EOP

# SCRIPTDIR = ../Scripts

PROTEXFILES = \
	program.protex \
	subroutine.protex \
	module.protex \
	internal_subroutine.protex \
	Makefile.protex

OBJECTS = \
	module.o \
	subroutine.o \
	internal_subroutine.o \
	program.o

README.dvi: README.tex PROTEX.tex

PROTEX.tex: ${PROTEXFILES}
	cat ${PROTEXFILES} > PROTEX.tex

# TEST TARGETS

#BOC
help:
	@echo 'The following targets are available:'
	@echo
	@echo 'test      - make DVI PS PDF COMPILE RUN'
	@echo 'DVI       - create README.dvi'
	@echo 'PS        - create README.ps'
	@echo 'PDF       - create README.pdf'
	@echo 'COMPILE   - compile the Fortran files into program_example.exe'
	@echo 'RUN       - run program_example.exe'
	@echo 'clean     - remove *~ and *.protex files'
	@echo 'distclean - remove all temporary files'
	@echo
#EOC

test:
	make DVI
	make PS
	make PDF
	make COMPILE
	make RUN
	@echo
	@echo 'All tests have been completed. You can read these files now:'
	@ls -l README.ps README.pdf
	@ls -l *.f90
	@echo

DVI: README.dvi

PS:  README.ps

PDF: README.pdf

RUN: COMPILE
	./program_example.exe

COMPILE: program_example.exe


# The compilation of object files is defined Makefile.conf
# Only the linking needs to be specified here
program_example.exe: ${OBJECTS}
	${LINK.f90} -o program_example.exe ${OBJECTS}


# CLEANING

# The clean target is defined in Makefile.conf

distclean: clean
	rm -f *.dvi *.ps *.pdf *.aux *.log *.toc PROTEX.tex program_example.exe

# RULES

# Since Makefile does not have an extension, this rule is explicit,
# but otherwise the rule is the same as for scripts with a .pl extension.

Makefile.protex: Makefile
	${SCRIPTDIR}/protex.pl -b -n -f -S Makefile > Makefile.protex

# The general rule Makefile.doc puts the PDF into the ../ directory. 
# To get the REAME.pdf in the local directory, explicit rule is used
README.pdf: README.ps
	ps2pdf README.ps README.pdf
