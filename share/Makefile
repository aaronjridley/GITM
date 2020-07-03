include ../Makefile.def

INSTALL_FILES = \
	Library/src/Makefile.DEPEND \
	Library/src/Makefile.RULES

install:
	touch ${INSTALL_FILES}

cleanall:
	cd Library/src; make clean
	cd Library/test; make clean
	rm -f include/*.mod

distclean: clean
	cd Library/src; make distclean
	cd Library/test; make distclean
	rm -f Library/src/mpif*.h *~ */*~ ${INSTALL_FILES}

