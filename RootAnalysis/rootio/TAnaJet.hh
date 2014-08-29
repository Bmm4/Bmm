#ifndef TANAJET
#define TANAJET

#define TANAJET_MAXTRK 200

#include "TObject.h"
#include "TArrayC.h"
#include "TArrayF.h"
#include "TVector3.h"
#include "TAnaVertex.hh"

class TAnaJet: public TObject {

public:
  TAnaJet(int index = -1);
  ~TAnaJet() { };
  void clear(); 
  void dump();

  // ----------------------------------------------------------------------

  // ----------------------------------------------------------------------
  int            fIndex;    // pointer into jet block (e.g. in signal block)

  // -- Tracking parameters
  int            fQ;                            // charge
 
  // -- Momentum
  TVector3       fPlab;
  double         fE; 
  double         fEt;
  double         fM;
  double         fMt;

  float          fEMEnergy; 
  float          fHADEnergy;
  float          finvisibleEnergy;
 
  int            fn60; 
  int            fn90;

  double         fD1,fD2,fD3,fD4,fD5,fD6,fD7; //b tagging discriminator

  int            fJetFlavorAlgo; 
  int            fJetFlavorPhys;
  int            fJetFlavorEne; //not used

  void           addTrack(int index);  
  int            getNtracks()    {return fNtracks;}
  int            getTrack(int n) {return fTracksIndex[n];}
  
  // -- Vertex
  TAnaVertex  fVtx;


private: 

  int            fNtracks;
  int            fTracksIndex[TANAJET_MAXTRK];

  ClassDef(TAnaJet,1)

};

#endif
