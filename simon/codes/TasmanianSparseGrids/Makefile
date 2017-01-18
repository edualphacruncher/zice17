UNAME = $(shell uname)
ifeq ($(UNAME), Linux)
COMPILE_OPTIONS = -fopenmp
else
COMPILE_OPTIONS =
endif

# You can overwrite COMPILE_OPTIONS here
# COMPILE_OPTIONS =

CC = g++
OPTC = -O3 -ffast-math -mtune=native $(COMPILE_OPTIONS) -Wall -Wextra -Wshadow -Wno-unused-parameter -pedantic -Wno-vla
OPTL = $(OPTC)
IADD =
LADD =

LHEADERS = tsgIndexSets.hpp tsgCoreOneDimensional.hpp tsgIndexManipulator.hpp tsgGridGlobal.hpp tsgCacheLagrange.hpp tsgSequenceOptimizer.hpp \
           tsgEnumerates.hpp tsgOneDimensionalWrapper.hpp tsgGridSequence.hpp tsgGridCore.hpp tsgLinearSolvers.hpp \
           tsgRuleLocalPolynomial.hpp tsgHardCodedTabulatedRules.hpp tsgGridLocalPolynomial.hpp \
           tsgRuleWavelet.hpp tsgGridWavelet.hpp \
           tasgridTestFunctions.hpp tasgridExternalTests.hpp tasgridWrapper.hpp \
           TasmanianSparseGrid.hpp

LIBOBJ = tsgIndexSets.o tsgCoreOneDimensional.o tsgIndexManipulator.o tsgGridGlobal.o tsgCacheLagrange.o tsgSequenceOptimizer.o tsgOneDimensionalWrapper.o \
         tsgGridCore.o tsgLinearSolvers.o tsgGridSequence.o tsgRuleLocalPolynomial.o tsgHardCodedTabulatedRules.o \
         tsgGridLocalPolynomial.o tsgRuleWavelet.o tsgGridWavelet.o \
         TasmanianSparseGrid.o 
		
WROBJ = tasgrid_main.o tasgridTestFunctions.o tasgridExternalTests.o tasgridWrapper.o

LIBNAME = libtasmaniansparsegrid.a

EXECNAME = tasgrid


%.o: %.cpp $(LHEADERS) $(TASGRIDHEADERS)
	$(CC) $(OPTC) -c $< -o $@ 

all: $(LIBNAME) $(EXECNAME)

$(LIBNAME): $(LIBOBJ) 
	ar rcs $(LIBNAME) $(LIBOBJ)

$(EXECNAME):  $(LIBNAME) $(WROBJ)
	$(CC) $(OPTL) $(LADD) -L. $(WROBJ) -o $(EXECNAME) $(LIBNAME) -lm


clean:
	rm -fr *.o
	rm -fr $(LIBNAME)
	rm -fr $(EXECNAME)
