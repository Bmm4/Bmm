# ======================================================================
# NOTE: make all
# -----
#
# ======================================================================


ROOTCINT      = $(ROOTSYS)/bin/rootcint
ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

ROOT6 = $(shell $(ROOTSYS)/bin/root-config --version)
ifneq ($(findstring 6.,$(ROOT6)),)
ROOT6         = "yes"
else 
ROOT6         = "no"
endif

ifdef SCRAM_ARCH
  CXX         := $(shell scramv1 tool info gcc-cxxcompiler | grep CXX= | sed s/CXX=//)
else 
  ifdef CXXCOMPILER
    CXX         := $(CXXCOMPILER)
  else 
    CXX         := c++
  endif
endif

CXXFLAGS      = -g -Os -Wall -fPIC -pipe
LD            = $(CXX)
LDFLAGS       = -g
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS)
GLIBS         = $(filter-out -lz, $(ROOTGLIBS))

# -- Files
UTIL     := PidTable.o PidData.o \
            AnalysisCuts.o AnalysisDistribution.o selpoint.o dataset.o \
            util.o initFunc.o hpl.o JSON.o
UTILDICT := ${UTIL:.o=Dict.o}

DICTHEADERS := ${UTIL:.o=Dict.h}


# -- Default rules
$(addprefix obj/,%.o) : %.cc %.hh %.icc
	$(CXX) $(CXXFLAGS)  -c $< -o $@

$(addprefix obj/,%.o) : %.cc %.hh
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(addprefix obj/,%.o) : %.cc 
	$(CXX) $(CXXFLAGS)  -c $< -o $@

%Dict.cc : %.hh %LinkDef.h
	$(ROOTCINT) -f $@ -c  $^ 
#ifneq ($(findstring root/6,$(ROOTSYS)),)
ifeq ($(ROOT6),"yes")
	mv $(subst Dict.cc,Dict_rdict.pcm,$@) lib
	cd ../../../../lib/$(SCRAM_ARCH)/ && rm -f $(subst Dict.cc,Dict_rdict.pcm,$@) && \
	ln -s ../../src/Bmm/RootAnalysis/common/lib/$(subst Dict.cc,Dict_rdict.pcm,$@) && \
	cd -; 
endif

%Dict.cc : %.hh
	$(ROOTCINT) -f $@ -c  $< 
ifeq ($(ROOT6),"yes")
	mv $(subst Dict.cc,Dict_rdict.pcm,$@) lib
	cd ../../../../lib/$(SCRAM_ARCH)/ && rm -f $(subst Dict.cc,Dict_rdict.pcm,$@) && \
	ln -s ../../src/Bmm/RootAnalysis/common/lib/$(subst Dict.cc,Dict_rdict.pcm,$@) && \
	cd -; 
endif

.PHONY: prep all clean vars

# ----------------------------------------------------------------------
all: vars prep lib 

lib: lib/libUtil.so 

lib/libUtil.so: $(addprefix obj/,$(UTILDICT) $(UTIL))
	$(CXX) $(SOFLAGS) $(addprefix obj/,$(UTIL) $(UTILDICT)) $(GLIBS) -o lib/libUtil.so

# -- preparatory setup
prep:
	mkdir -p obj bin lib
	mkdir -p ../../../../lib/$(SCRAM_ARCH)

# -- clean up
clean:
	rm -f $(addprefix obj/,$(UTIL) $(UTILDICT))
	mkdir -p ../../../../lib/$(SCRAM_ARCH)
	cd ../../../../lib/$(SCRAM_ARCH)/ && rm -f libUtil.so *_rdict.pcm && cd - 
	rm -f lib/*

# -- ensure that cmsenv has been run
vars:
ifndef CMSSW_BASE
    $(error CMSSW_BASE is undefined, please run cmsenv somewhere)
endif
