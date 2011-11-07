CMP = mpic++
CXXFLAGS = -O3
FLAGS =

AR = ar

default: lib

clean:
	rm -rf *~ *.o *.a vlsv2silo

DEPS_COMMON = muxml.h vlsv_common.h
DEPS_MUXML = muxml.h muxml.cpp
DEPS_VLSVCOMMON = vlsv_common.h vlsv_common.cpp
DEPS_READER = ${DEPS_VLSCOMMON} vlsvreader.h vlsvreader2.cpp
DEPS_WRITER = ${DEPS_VLSCOMMON} vlsvwriter.h vlsvwriter2.cpp
DEPS_VLSV2SILO = vlsvreader.o muxml.o vlsv_common.o vlsv2silo.cpp

OBJS=muxml.o vlsv_common.o vlsvreader.o vlsvwriter.o

INC=-I/home/sandroos/codes/util

lib: ${OBJS}
	${AR} r libvlsv.a ${OBJS}
	make vlsv2silo

muxml.o: ${DEPS_MUXML}
	${CMP} ${CXXFLAGS} ${FLAGS} -o muxml.o -c muxml.cpp

vlsv_common.o: ${DEPS_VLSVCOMMON}
	${CMP} ${CXXFLAGS} ${FLAGS} -o vlsv_common.o -c vlsv_common.cpp

vlsvreader.o: ${DEPS_READER}
	${CMP} ${CXXFLAGS} ${FLAGS} -o vlsvreader.o -c vlsvreader2.cpp ${INC}

vlsvwriter.o: ${DEPS_WRITER}
	${CMP} ${CXXFLAGS} ${FLAGS} -o vlsvwriter.o -c vlsvwriter2.cpp ${INC}

vlsv2silo: ${DEPS_VLSV2SILO}
	${CMP} ${CXXFLAGS} ${FLAGS} -o vlsv2silo vlsv2silo.cpp -L${CURDIR} -lvlsv -lsilo
