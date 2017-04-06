# ======================================================================
# NOTE: make all links
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
ANA00 := TAna01Event.o TGenCand.o TAnaTrack.o TSimpleTrack.o TAnaMuon.o TTrgObj.o TTrgObjv2.o TAnaCand.o TAnaVertex.o TAnaJet.o
ANA00DICT := ${ANA00:.o=Dict.o}

DICTHEADERS := ${ANA00:.o=Dict.h}


# -- Default rules
$(addprefix obj/,%.o) : %.cc %.hh %.icc
	$(CXX) $(CXXFLAGS)  -c $< -o $@

$(addprefix obj/,%.o) : %.cc %.hh
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(addprefix obj/,%.o) : %.cc
	$(CXX) $(CXXFLAGS)  -c $< -o $@

%Dict.cc : %.hh %LinkDef.h
	$(ROOTCINT) -f $@ -c  $^

ifeq ($(ROOT6),"yes")
	mv $(subst Dict.cc,Dict_rdict.pcm,$@) lib
	cd ../../../../lib/$(SCRAM_ARCH)/ && rm -f $(subst Dict.cc,Dict_rdict.pcm,$@) && \
	ln -s ../../src/Bmm/RootAnalysis/rootio/lib/$(subst Dict.cc,Dict_rdict.pcm,$@) && \
	cd -;
endif

%Dict.cc : %.hh
	$(ROOTCINT) -f $@ -c  $<
ifeq ($(ROOT6),"yes")
	mv $(subst Dict.cc,Dict_rdict.pcm,$@) lib
	cd ../../../../lib/$(SCRAM_ARCH)/ && rm -f $(subst Dict.cc,Dict_rdict.pcm,$@) && \
	ln -s ../../src/Bmm/RootAnalysis/rootio/lib/$(subst Dict.cc,Dict_rdict.pcm,$@) && \
	cd -;
endif

.PHONY: prep all clean vars

# ----------------------------------------------------------------------
all: vars prep lib

lib: lib/libAna00.so

lib/libAna00.so: $(addprefix obj/, $(ANA00DICT) $(ANA00))
	$(CXX) $(SOFLAGS) $(GLIBS) $(addprefix obj/,$(ANA00) $(ANA00DICT)) -o lib/libAna00.so

# ================================================================================
links:
	mkdir -p ../../../../lib/$(SCRAM_ARCH)
	cd ../../../../lib/$(SCRAM_ARCH)/ && rm -f libAna00.so && ln -s ../../src/Bmm/RootAnalysis/rootio/lib/libAna00.so && cd -
a
copy:
	mkdir -p ../../../../lib/$(SCRAM_ARCH)
	cd ../../../../lib/$(SCRAM_ARCH)/ && rm -f libAna00.so && cp ../../src/Bmm/RootAnalysis/rootio/lib/libAna00.so . && cd -


# -- preparatory setup
prep:
	mkdir -p obj bin lib
	mkdir -p ../../../../lib/$(SCRAM_ARCH)

# -- clean up
clean:
	rm -f $(addprefix obj/,$(ANA00) $(ANA00DICT))
	mkdir -p ../../../../lib/$(SCRAM_ARCH)
	cd ../../../../lib/$(SCRAM_ARCH)/ && rm -f libAna00.so  *_rdict.pcm && cd -
	rm -f lib/*

# -- ensure that cmsenv has been run
vars:
ifndef CMSSW_BASE
    $(error CMSSW_BASE is undefined, please run cmsenv somewhere)
endif
