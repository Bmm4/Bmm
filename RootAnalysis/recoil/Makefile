# ======================================================================
# NOTE: make all
# -----
#
# ======================================================================
.PHONY: prep clean vars

ROOTCINT      = $(ROOTSYS)/bin/rootcint
ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
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
GLIBS         = $(filter-out -lz, $(ROOTGLIBS)) -lTMVA -lTMVAGui -lRooFitCore -lRooFit
EXTHEADERS    = -I..


# -- Files
READER = t1Reader.o recoilReader.o genAnalysis.o \
         candAna.o candAnaRecoil.o tmva1.o

DICT = ${READER:.o=Dict.o}
DICTHEADERS = ${READER:.o=Dict.h}
DICTSOURCES = ${READER:.o=Dict.cc}

#OTHERPCM = $(shell find ../rootio/lib/*.pcm ../common/lib/*.pcm)


# -- Default rules
$(addprefix obj/,%.o) : %.cc %.hh %.icc
	$(CXX) $(CXXFLAGS) $(EXTHEADERS) -c $< -o $@

$(addprefix obj/,%.o) : %.cc %.hh
	$(CXX) $(CXXFLAGS) $(EXTHEADERS) -c $< -o $@

$(addprefix obj/,%.o) : %.cc
	$(CXX) $(CXXFLAGS) $(EXTHEADERS) -c $< -o $@

%Dict.cc : %.hh %LinkDef.h
	$(ROOTCINT) -f $@ -c $(EXTHEADERS) $^
ifeq ($(ROOT6),"yes")
	mv $(subst Dict.cc,Dict_rdict.pcm,$@) lib
endif

%Dict.cc : %.hh
	$(ROOTCINT) -f $@ -c $(EXTHEADERS) $<
ifeq ($(ROOT6),"yes")
	mv $(subst Dict.cc,Dict_rdict.pcm,$@) lib
endif

$(addprefix lib/,%.pcm) :
	cd lib && ln -s $< && cd -


# -- Targets
all: vars prep lib bin

lib: $(addprefix obj/,$(CSTUFF) $(READER)  $(DICT))
	$(CXX) $(SOFLAGS) $(addprefix obj/,$(CSTUFF) $(READER) $(DICT)) -o lib/libRecoilClasses.so $(GLIBS) -lMinuit lib/libAna00.so lib/libUtil.so

bin: lib obj/runT1Reader.o obj/runPlot.o
	$(LD) $(LDFLAGS) -o bin/runT1Reader $(GLIBS) obj/runT1Reader.o lib/libRecoilClasses.so lib/libAna00.so lib/libUtil.so
	$(LD) $(LDFLAGS) -o bin/runPlot $(GLIBS) obj/runPlot.o  lib/libRecoilClasses.so lib/libAna00.so lib/libUtil.so

# -- preparatory setup
prep:
	mkdir -p obj bin lib
	cd lib && ln -f -s ../../rootio/lib/libAna00.so && cd -
	cd lib && /bin/ls -1 ../../rootio/lib/*.pcm | perl -ne 'use File::Basename; chop; symlink($$_, basename($$_))' && cd -
	cd lib && ln -f -s ../../common/lib/libUtil.so && cd -
	cd lib && /bin/ls -1 ../../common/lib/*.pcm | perl -ne 'use File::Basename; chop; symlink($$_, basename($$_))' && cd -

# -- clean up
clean:
	rm -f $(addprefix obj/,$(READER) $(DICT) $(CSTUFF) runT1Reader.o runPlot.o)
	rm -f $(DICTHEADERS)  $(DICTSOURCES)
	rm -f bin/*
	rm -f lib/*

cleanresults:
	rm -f results/*.root results/*.pdf results/*.tex results/*.tex.old

# -- ensure that the environment variable DELPHES is set
vars:
ifndef CMSSW_BASE
    $(error CMSSW_BASE is undefined, please run cmsenv somewhere)
endif
