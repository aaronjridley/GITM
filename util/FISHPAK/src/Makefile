include ../../../Makefile.def
include ../../../Makefile.conf

help:
	@echo "    "
	@echo "make MODULE             - create ModFish.o"
	@echo "make LIB                - create libFISHPAK.a in ${LIBDIR}"
	@echo "make test               - run all tests"
	@echo "make test TEST=blktri   - run a specific test"
	@echo "clean                   - remove object files etc"

MODULE: ModFish.o
	
LIB: ModFish.o
	${AR} ${LIBDIR}/libFISH.a ModFish.o

TEST = blktri cblktri cmgnbn genbun hstcrt hstcsp hstcyl hstplr hstssp hw3crt hwscrt hwscsp hwscyl hwsplr hwsssp pois3d poistg sepeli sepx4

test: ModFish.o
	for i in ${TEST}; do \
		make t$$i.o; \
		${COMPILE.f90} -o test.exe t$$i.o ModFish.o; \
		./test.exe; \
	done
