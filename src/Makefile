# GCC options (Linux assumed, should work on all gcc platforms.)
CXX = g++
CC = gcc
DEBUGOPT = 
WARNOPT =
OPTOPT = -O3 #-fno-rtti -fno-exceptions #-mfpmath=sse
LIBS = 
LINKOPT = -s

DEBUG=1

ifdef DEBUG
DEBUGOPT = -g
WARNOPT = -Wall -pedantic
OPTOPT =
LINKOPT =
endif

# MIPSPRO (Irix) options (no debug available in this case)

ifdef MIPS
CXX = CC
CC = cc
DEBUGOPT =
WARNOPT = 
OPTOPT = -O3 -mips4
LIBS = -L /share/local/libn32 -lm
endif

# OSF1 (Tru64) options (again, only an optimized build)
ifdef OSF1
CXX = cxx
CC = cc
DEBUGOPT = -pg
WARNOPT = -w0
OPTOPT = -arch host -tune host -O4 -assume nomath_errno -om 
LIBS = -lm
endif

LIBOBJS = node.o etree.o datacache.o 

EXEFILES = ../etreetrain ../etreequery ../etreenewick
EXEOBJS = etreetrain.o etreequery.o etreenewick.o

CXXFLAGS = $(DEBUGOPT) $(WARNOPT) $(OPTOPT)

all:   	$(EXEFILES)



../etreetrain: $(LIBOBJS) etreetrain.o
	$(CXX) $(LINKOPT) -o ../etreetrain $(LIBOBJS) etreetrain.o $(LIBS)

../etreequery: $(LIBOBJS) etreequery.o
	$(CXX) $(LINKOPT) -o ../etreequery $(LIBOBJS) etreequery.o $(LIBS)

../etreenewick: $(LIBOBJS) etreenewick.o
	$(CXX) $(LINKOPT) -o ../etreenewick $(LIBOBJS) etreenewick.o $(LIBS)


docs:
	doxygen

clean :: 
	rm -f *.o
	rm -rf ii_files
	rm -rf cxx_repository

realclean:: clean ; rm -f $(EXEFILES)
