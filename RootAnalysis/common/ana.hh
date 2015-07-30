#ifndef ANA_H
#define ANA_H

#include "rootio/TAnaTrack.hh"
#include "rootio/TAnaMuon.hh"

int  numberOfBPixLayers(TAnaTrack *t);
int  numberOfPixLayers(TAnaTrack *t);
int  numberOfPixelHits(TAnaTrack *pTrack);
int  numberOfBPixLayer1Hits(TAnaTrack *t);
int  numberOfTrackerLayers(TAnaTrack *t);

bool tightMuon(TAnaMuon *pM, bool hadronsPass = false, int year = -1); 


#endif
