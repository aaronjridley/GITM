include ../Makefile.def

INSTALL_FILES = \
	Library/src/Makefile.DEPEND \
	Library/src/Makefile.RULES

install:
	touch ${INSTALL_FILES}

clean:
	cd Library/src; make clean
	cd Library/test;make clean
	cd Prologs;     make clean
	rm -f include/*.mod

distclean: clean
	cd Library/src; make distclean
	cd Library/test;make distclean
	cd Prologs;     make distclean
	rm -f Library/src/mpif*.h *~ */*~ ${INSTALL_FILES}

