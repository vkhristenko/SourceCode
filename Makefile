# MY_PROGRAM = quickextract_wBetaGamma
# MY_PROGRAM = quickextract_wBetaGamma_wE
MY_PROGRAM = analyze
# MY_PROGRAM = quickreco
# MY_PROGRAM = quickextract_0

# MY_PROGRAM = photonanalyzer
# MY_PROGRAM = photonplotter

##############
# ROOT STUFF #
##############

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTCINT      = $(ROOTSYS)/bin/rootcint

####################
# BASE DIRECTORIES #
####################

BASE_DIR = `pwd`

############
# INCLUDES #
############

ROOT_INC = -I$(ROOTSYS)/include

EVERYTHING_INCLUDE = -I$(BASE_DIR) $(ROOT_INC) 

#############
# LIBRARIES #
#############

EVERYTHING_LIBS =  $(ROOTLIBS)
# -lCDFROOT

#############
# Compilers #
#############

      CC  = cc
      CPP = `which g++`

####################
# Compiler Options #
####################

       CPP_OPTIONS = -g -O0 $(ROOTCFLAGS)

################
# Suffix rules #
################
#
.SUFFIXES: .o .C .cpp .cc
.C.o:
	@echo
	@echo  COMPILING $<
	@echo
	rm -f $@
	$(CPP) -c $(CPP_OPTIONS) $(EVERYTHING_INCLUDE) $<
.cc.o:
	@echo
	@echo  COMPILING $<
	@echo
	rm -f $@
	$(CPP) -c $(CPP_OPTIONS) $(EVERYTHING_INCLUDE) $<

.cpp.o:
	@echo
	@echo  COMPILING $<
	@echo
	rm -f $@
	$(CPP) -c $(CPP_OPTIONS) $(EVERYTHING_INCLUDE) $<

MY_OBJECTS = $(MY_PROGRAM)

$(MY_PROGRAM) : $(MY_OBJECTS:=.o)
	@echo
	@echo Linking $@
	@echo  	
	@echo  	
	$(CPP) $(CPP_OPTIONS) -o $@ $(MY_OBJECTS:=.o) $(EVERYTHING_INCLUDE) $(EVERYTHING_LIBS)

clean :
	rm -f $(MY_PROGRAM)
	rm -f  *.o
	rm -f  *.so 
	rm -f  *.l 
	rm -f  *.a
	rm -f  *.exe 
	rm -f  core


