# ======================================================================
# NOTE: make all
# -----
#
# ======================================================================


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

CXXFLAGS      = -g -Os -Wall -fPIC -pipe
LD            = $(CXX)
LDFLAGS       = -g
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS)
GLIBS         = $(filter-out -lz, $(ROOTGLIBS))

# -- Files
ANA00 := TAna01Event.o TGenCand.o TAnaTrack.o TSimpleTrack.o TAnaMuon.o TTrgObj.o TAnaCand.o TAnaVertex.o TAnaJet.o
ANA00DICT := ${ANA00:.o=Dict.o}

UTIL := PidTable.o PidData.o functions.o util.o hpl.o JSON.o
UTILDICT := ${UTIL:.o=Dict.o}

DICTHEADERS := ${ANA00:.o=Dict.h} ${UTIL:.o=Dict.h}


# -- Default rules
$(addprefix obj/,%.o) : %.cc %.hh %.icc
	$(CXX) $(CXXFLAGS)  -c $< -o $@

$(addprefix obj/,%.o) : %.cc %.hh
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(addprefix obj/,%.o) : %.cc 
	$(CXX) $(CXXFLAGS)  -c $< -o $@

%Dict.cc : %.hh %LinkDef.h
	$(ROOTCINT) -f $@ -c  $^ 
	mv $(subst Dict.cc,Dict_rdict.pcm,$@) lib
	cd ../../../../lib/$(SCRAM_ARCH)/ && rm -f $(subst Dict.cc,Dict_rdict.pcm,$@) && ln -s ../../src/Bmm/RootAnalysis/rootio/lib/$(subst Dict.cc,Dict_rdict.pcm,$@) && cd -; 

%Dict.cc : %.hh
	$(ROOTCINT) -f $@ -c  $< 
	mv $(subst Dict.cc,Dict_rdict.pcm,$@) lib
	cd ../../../../lib/$(SCRAM_ARCH)/ && rm -f $(subst Dict.cc,Dict_rdict.pcm,$@) && ln -s ../../src/Bmm/RootAnalysis/rootio/lib/$(subst Dict.cc,Dict_rdict.pcm,$@) && cd -; 

.PHONY: prep all clean vars

# ----------------------------------------------------------------------
all: vars prep lib 

lib: lib/libAna00.so lib/libUtil.so 

lib/libAna00.so: $(addprefix obj/, $(ANA00DICT) $(ANA00))
	$(CXX) $(SOFLAGS) $(GLIBS) $(addprefix obj/,$(ANA00) $(ANA00DICT)) -o lib/libAna00.so

lib/libUtil.so: $(addprefix obj/,$(UTILDICT) $(UTIL))
	$(CXX) $(SOFLAGS) $(addprefix obj/,$(UTIL) $(UTILDICT)) $(GLIBS) -o lib/libUtil.so

# ================================================================================
links:
	mkdir -p ../../../../lib/$(SCRAM_ARCH)
	cd ../../../../lib/$(SCRAM_ARCH)/ && rm -f libAna00.so && ln -s ../../src/Bmm/RootAnalysis/rootio/lib/libAna00.so && cd -
	cd ../../../../lib/$(SCRAM_ARCH)/ && rm -f libUtil.so && ln -s ../../src/Bmm/RootAnalysis/rootio/lib/libUtil.so && cd -


copy:
	cd ../../../../lib/$(SCRAM_ARCH)/ && rm -f libAna00.so && cp ../../src/Bmm/RootAnalysis/rootio/lib/libAna00.so . && cd -
	cd ../../../../lib/$(SCRAM_ARCH)/ && rm -f libUtil.so && cp ../../src/Bmm/RootAnalysis/rootio/lib/libUtil.so . && cd -


# -- preparatory setup
prep:
	mkdir -p obj bin lib

# -- clean up
clean:
	rm -f $(addprefix obj/,$(ANA00) $(ANA00DICT) $(UTIL) $(UTILDICT)) 
	cd ../../../../lib/$(SCRAM_ARCH)/ && rm -f libAna00.so libUtil.so *_rdict.pcm && cd - 
#	rm -f $(DICTHEADERS) 
	rm -f bin/runBmm
	rm -f lib/*

# -- ensure that cmsenv has been run
vars:
ifndef CMSSW_BASE
    $(error CMSSW_BASE is undefined, please run cmsenv somewhere)
endif
