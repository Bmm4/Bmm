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
GLIBS         = $(filter-out -lz, $(ROOTGLIBS)) -lTMVA
EXTHEADERS    = -I.. 
#EXTHEADERS    = -I.. -I../fastjet-install/include/
#FASTJETLIBS   = `../fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`

# -- Files
CSTUFF = preselection.o 
READER = treeReader01.o inclbReader.o genAnalysis.o candAna.o  plotClass.o plotWork.o
DICT = ${READER:.o=Dict.o}
DICTHEADERS = ${READER:.o=Dict.h}
DICTSOURCES = ${READER:.o=Dict.cc}


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
	$(CXX) $(SOFLAGS) $(addprefix obj/,$(CSTUFF) $(READER) $(DICT)) -o lib/libAnaInclBClasses.so $(GLIBS) -lMinuit lib/libAna00.so lib/libUtil.so 
#	$(CXX) $(SOFLAGS) $(addprefix obj/,$(CSTUFF) $(READER) $(DICT)) -o lib/libAnaInclBClasses.so $(GLIBS) -lMinuit lib/libAna00.so lib/libUtil.so $(FASTJETLIBS)


bin: lib obj/runInclB.o
	$(LD) $(LDFLAGS) -o bin/runInclB $(GLIBS) obj/runInclB.o lib/libAnaInclBClasses.so lib/libAna00.so lib/libUtil.so 
#	$(LD) $(LDFLAGS) -o bin/runInclB $(GLIBS) obj/runInclB.o lib/libAnaInclBClasses.so lib/libAna00.so lib/libUtil.so $(FASTJETLIBS)


# -- preparatory setup
prep:
	mkdir -p obj bin lib
	cd lib && ln -f -s ../../rootio/lib/libAna00.so && cd - 
	cd lib && /bin/ls -1 ../../rootio/lib/*.pcm | perl -ne 'use File::Basename; chop; symlink($$_, basename($$_))' && cd - 
	cd lib && ln -f -s ../../common/lib/libUtil.so && cd - 
	cd lib && /bin/ls -1 ../../common/lib/*.pcm | perl -ne 'use File::Basename; chop; symlink($$_, basename($$_))' && cd - 


# -- clean up
clean:
	rm -f $(addprefix obj/,$(READER) $(DICT) $(CSTUFF) runInclB.o) 
	rm -f $(DICTHEADERS)  $(DICTSOURCES) 
	rm -f lib/* obj/* bin/*

# -- ensure that the environment variable DELPHES is set
vars:
ifndef CMSSW_BASE
    $(error CMSSW_BASE is undefined, please run cmsenv somewhere)
endif
