# Generated automatically from Makefile.in by configure.
#+
# File: Makefile.in
# Description : Makefile Template for fsim_analysis
#-

# 1. User Specifications

OBJS = reco.o userinfo.o

MODULE	= Reco.so

LIBS =	-L$(BELLE_RUN_DIR)/lib/so \
	-ltuple -lmdst -lip -lkid -leid -lparticle -lkfitter -lbenergy


# -lekpfullrecon -lekpcontsuppress

#	-lpntdb \
#	-L/belle/local/depot/postgresql~6321/lib -lpq


# 2. System Specifications
#    --- Do not change without knowledge

# Compiler Setup with machine dependence

FC = f77
CC = gcc
CXX = g++

DEFS =  -DHAVE_LIBCURSES=1 -DHAVE_LIBREADLINE=1 -DHAVE_POSTGRES=1 -DHAVE_LIBCURSES=1 -DHAVE_LIBTERMCAP=1 -DHAVE_LIBHISTORY=1 -DHAVE_LIBREADLINE=1 -DHAVE_HISTORY=1 -DHAVE_LIBBSD=1 -DHAVE_LIBM=1 -DHAVE_LIBDL=1 -DHAVE_LIBNSL=1 -DHAVE_LIBCRYPT=1 -DHAVE_LIBNSL=1 -DHAVE_LIBDL=1 -DFORTRAN_PPU=1 -DHAVE_LIBCRYPT=1  -DCERNLIB_TYPE
CPPFLAGS = 
DEPCPPFLAGS = -MM

FFLAGS = -DBELLE_TARGET_H=\"belle-default.h\"  -fno-second-underscore -fno-automatic -finit-local-zero -fno-emulate-complex
CFLAGS = -O3 -DBELLE_TARGET_H=\"belle-default.h\"
CXXFLAGS = -O3 -DHEP_SHORT_NAMES -DBELLE_SHORT_NAMES -DDSTXX_NOINLINE -DBELLE_TARGET_H=\"belle-default.h\" -fPIC

SOFLAGS = -shared -Wl,-export-dynamic
LDFLAGS = 

SYSLIB = -lcrypt   -L/sw/belle/local/lib/gcc-lib/i686-pc-linux-gnu/3.0.4 -L/sw/belle/local/lib/gcc-lib/i686-pc-linux-gnu/3.0.4/../../.. -lm -lgcc -lgcc

CLHEPLIB = -lbelleCLHEP

MOTIF_LIBS = -L/usr/X11R6/lib -lXm


# Include directories

POSTGRESINC = /sw/belle/local/include

INCLUDES_C = $(BELLE_TOP_DIR)/include
INCLUDES_FORTRAN = $(BELLE_TOP_DIR)/inc

# Dependence description

include $(BELLE_RUN_DIR)/src/config/Makefile.panther

COMPILE_FCPP := $(FC) -c $(PANTHER_FMACROS) $(INCLUDES_FORTRAN:%=-I%) $(CPPFLAGS) $(FFLAGS)
COMPILE_FC := $(FC) -c  $(INCLUDES_FORTRAN:%=-I%) $(FFLAGS)
COMPILE_CC := $(CC) -c  $(PANTHER_CMACROS) $(INCLUDES_C:%=-I%) $(CPPFLAGS) $(CFLAGS)
COMPILE_CXX := $(CXX) -c  $(PANTHER_CMACROS) $(INCLUDES_C:%=-I%) $(CPPFLAGS) $(CXXFLAGS)

LINK_FCPP := $(FC)
LINK_FC := $(FC)
LINK_CC := $(CC)
LINK_CXX := $(CXX)

DEPEND_FCPP := /lib/cpp -M $(DEFS) $(PANTHER_FMACROS) $(INCLUDES_FORTRAN:%=-I%) $(CPPFLAGS) $(FFLAGS)
DEPEND_CC := $(CC) $(DEPCPPFLAGS) $(DEFS) $(PANTHER_CMACROS) $(INCLUDES_C:%=-I%) $(CPPFLAGS) $(CFLAGS)
DEPEND_CXX := $(CXX) $(DEPCPPFLAGS) $(DEFS) $(PANTHER_CMACROS) $(INCLUDES_C:%=-I%) $(CPPFLAGS) $(CXXFLAGS)

LINUX_G77_BUG = -e 's/\.F\.o:/.o:/'

%.o:%.c
	$(COMPILE_CC) $<

%.d:%.c
	$(SHELL) -ec '$(DEPEND_CC) $< | sed -e "s/$*.o[ :]*/$@ &/g" -e 's/\.[12][0-9][0-9][0-9][0-9][0-9][0-9][0-9][a-z]\.h/.tdf/g' > $@'

%.o:%.cc
	$(COMPILE_CXX) $<

%.d:%.cc
	$(SHELL) -ec '$(DEPEND_CXX) $< | sed -e "s/$*.o[ :]*/$@ &/g" -e 's/\.[12][0-9][0-9][0-9][0-9][0-9][0-9][0-9][a-z]\.h/.tdf/g'> $@'


%.o:%.F
	$(COMPILE_FCPP) $<

%.d:%.F
	$(SHELL) -ec '$(DEPEND_FCPP) $< | sed $(LINUX_G77_BUG) -e "s/$*.o[ :]*/$@ &/g" -e 's/\.[12][0-9][0-9][0-9][0-9][0-9][0-9][0-9][a-z]\.inc/.tdf/g'> $@'


# CERNLIB

ifeq "$(CERN)/$(CERN_LEVEL)" "/"
  CERNLIB_LIB_DIR = /sw/belle/cern/2006/lib
else
  CERNLIB_LIB_DIR = $(CERN)/$(CERN_LEVEL)/lib
endif


# Dependencies

all::	$(OBJS)
	$(LINK_CXX) -o $(MODULE) $(SOFLAGS) $(OBJS) $(LIBS) \
	$(CLHEPLIB) $(CERNLIB) $(SYSLIB)

check:
	./run.csh $(BELLE_LEVEL) > ./out 2>&1

clean::
	rm -f $(OBJS) $(MODULE) *~ fort.* fpda_pid.* test.hbk out

