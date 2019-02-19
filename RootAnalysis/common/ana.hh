#ifndef ANA_H
#define ANA_H

#include "rootio/TAna01Event.hh"

int  numberOfBPixLayers(TAnaTrack *t);
int  numberOfPixLayers(TAnaTrack *t);
int  numberOfPixelHits(TAnaTrack *pTrack);
int  numberOfBPixLayer1Hits(TAnaTrack *t);
int  numberOfTrackerLayers(TAnaTrack *t);

bool highPurity(TAnaTrack *pt);

bool tightMuon(TAnaMuon *pM, bool hadronsPass = false, int year = -1);


void dumpCand(TAnaCand *pCand, TAna01Event *pEvt);

#endif
