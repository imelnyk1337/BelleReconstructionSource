# Generated automatically from Makefile.in by configure.
#+
# File: Makefile.in
# Description : Makefile Template for gsim_analysis
#-

# 1. User Specifications


OBJS = combination.o utility.o exkfitter.o

MODULE = exkfitter.so

LIBS = -L$(BELLE_RUN_DIR)/lib/so -ltuple

# 2. System Specifications
#    --- Do not change without knowledge

# Compiler Setup with machine dependence

FC = f77
CC = gcc
CXX = c++

DEFS =  -DWORDS_BIGENDIAN=1 -DHAVE_LIBCURSES=1 -DHAVE_LIBREADLINE=1 -DHAVE_LIBNSL=1 -DHAVE_LIBSOCKET=1 -DHAVE_LIBDL=1 -DFORTRAN_PPU=1 -DHAVE_LIBCRYPT=1  -DCERNLIB_TYPE
CPPFLAGS = 
DEPCPPFLAGS = -MM

FFLAGS = -O -g -DBELLE_TARGET_H=\"belle-sparc-sun-solaris2.5.1-c++.h\" -PIC
CFLAGS = -O -g -DBELLE_TARGET_H=\"belle-sparc-sun-solaris2.5.1-c++.h\" -fPIC
CXXFLAGS = -g -DHEP_SHORT_NAMES -DBELLE_SHORT_NAMES -DDSTXX_NOINLINE -DBELLE_TARGET_H=\"belle-sparc-sun-solaris2.5.1-c++.h\" -fPIC
SOFLAGS = -G
LDFLAGS = 

SYSLIB = -L/opt/kekb/programs/SUNWspro/lib -L/opt/kekb/programs/SUNWspro/SC4.0/lib -L/opt/SUNWspro/lib -L/usr/ccs/lib -L/usr/lib /opt/kekb/programs/SUNWspro/SC4.0/lib/libM77.a -lF77 -lsunmath -lm -lcx

PARTICLELIB = -lparticle -lkfitter -lhelix -lmdst -L/belle/local/lib -lCLHEP-$(CXX)
CLHEPINC = /belle/local/include
CLHEPLIB = -L/$(BELLE_TOP_DIR)/lib -lkfitter -lhelix -lqqutil -L/belle/local/lib -lCLHEP-$(CXX)

# Include directories

INCLUDES_C = $(BELLE_TOP_DIR)/include $(CLHEPINC)
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

LINUX_G77_BUG = @LINUX_G77_BUG@

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
  CERNLIB_LIB_DIR = /bn10/belle/cern/pro/lib
else
  CERNLIB_LIB_DIR = $(CERN)/$(CERN_LEVEL)/lib
endif


# Dependencies

all::	$(OBJS)
	$(LINK_CXX) -o $(MODULE) $(SOFLAGS) $(OBJS) $(PARTICLELIB) $(LIB) \
	$(CLHEPLIB) $(CERNLIB) $(SYSLIB)

check:
	./run.csh $(BELLE_LEVEL) > ./out 2>&1

clean::
	rm -f $(OBJS) $(MODULE) *~ fort.* fpda_pid.* test.hbk panther.dat \
	cdcgeo.bdat* panther.dat out cdchit.paw





