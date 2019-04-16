# ======================================================================
# NOTE: make all links
# -----
#
# ======================================================================

ROOTBASE      = $(BASE)/RootAnalysis

ROOTCINT      = $(ROOTSYS)/bin/rootcint
ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

ifdef SCRAM_ARCH
  CXX         := $(shell scramv1 tool info gcc-cxxcompiler | grep CXX= | sed s/CXX=//)
else
  ifdef CXXCOMPILER
    CXX         := $(CXXCOMPILER)
  else
    CXX         := c++
  endif
endif

ROOT6 = $(shell $(ROOTSYS)/bin/root-config --version)
ifneq ($(findstring 6.,$(ROOT6)),)
ROOT6         = "yes"
else
ROOT6         = "no"
endif



CXXFLAGS      = -g -O0 -Wall -fPIC -pipe
LD            = $(CXX)
LDFLAGS       = -g
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS)
GLIBS         = $(filter-out -lz, $(ROOTGLIBS))
EXTHEADERS    = -I..

# -- Files
UTIL     := AnalysisCuts.o AnalysisDistribution.o  \
            PidTable.o PidData.o \
            HistCutEfficiency.o selpoint.o dataset.o \
            util.o initFunc.o hpl.o JSON.o PdTrigger.o Lumi.o ana.o fitPsYield.o
UTILDICT := ${UTIL:.o=Dict.o}

DICTHEADERS := ${UTIL:.o=Dict.h}


# -- Default rules
$(addprefix obj/,%.o) : %.cc %.hh %.icc
	$(CXX) $(CXXFLAGS) $(EXTHEADERS) -c $< -o $@

$(addprefix obj/,%.o) : %.cc %.hh
	$(CXX) $(CXXFLAGS)  $(EXTHEADERS) -c $< -o $@

$(addprefix obj/,%.o) : %.cc
	$(CXX) $(CXXFLAGS)  $(EXTHEADERS) -c $< -o $@

%Dict.cc : %.hh %LinkDef.h
	$(ROOTCINT) -f $@ -c $(EXTHEADERS) $^

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
	$(CXX) $(SOFLAGS) $(addprefix obj/,$(UTIL) $(UTILDICT)) $(GLIBS) -L$(ROOTBASE)/rootio/lib -lAna00 -o lib/libUtil.so
##	$(CXX) $(SOFLAGS) $(addprefix obj/,$(UTIL) $(UTILDICT)) $(GLIBS) -o lib/libUtil.so

# ======================================================================
links:
	mkdir -p ../../../../lib/$(SCRAM_ARCH)
	cd ../../../../lib/$(SCRAM_ARCH)/ && rm -f libUtil.so  && ln -s ../../src/Bmm/RootAnalysis/common/lib/libUtil.so && cd -

copy:
	mkdir -p ../../../../lib/$(SCRAM_ARCH)
	cd ../../../../lib/$(SCRAM_ARCH)/ && rm -f libUtil.so  && cp ../../src/Bmm/RootAnalysis/common/lib/libUtil.so  . && cd -


# -- preparatory setup
prep:
	mkdir -p obj bin lib
	mkdir -p ../../../../lib/$(SCRAM_ARCH)
	cd lib && ln -f -s ../../rootio/lib/libAna00.so && cd -

# -- clean up
clean:
	rm -f $(addprefix obj/,$(UTIL) $(UTILDICT))
	rm -f obj/* lib/* bin/*

# -- ensure that cmsenv has been run
vars:
ifndef CMSSW_BASE
    $(error CMSSW_BASE is undefined, please run cmsenv somewhere)
endif
