#include "TH1.h"
#include "rootio/TAna01Event.hh"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ function hpl(TH1 *, const char *);
#pragma link C++ function hpl(const char *, const char *);

#pragma link C++ function numberOfBPixLayers(TAnaTrack *);
#pragma link C++ function numberOfPixLayers(TAnaTrack *);
#pragma link C++ function numberOfPixelHits(TAnaTrack *);
#pragma link C++ function numberOfBPixLayer1Hits(TAnaTrack *);
#pragma link C++ function numberOfTrackerLayers(TAnaTrack *);

#pragma link C++ function tightMuon(TAnaMuon *, bool, int);
#pragma link C++ function highPurity(TAnaTrack *);

#pragma link C++ function dumpCand(TAnaCand *, TAna01Event*);

#endif
