#ifndef TANATRACK
#define TANATRACK


#include "TObject.h"
#include "TArrayC.h"
#include "TArrayF.h"
#include "TVector3.h"

class TAnaTrack: public TObject {

public:
  TAnaTrack(int index = -1);
  ~TAnaTrack() { };
  void clear(); 
  void dump();

  // ----------------------------------------------------------------------

  // ----------------------------------------------------------------------
  int            fIndex;    // pointer into (simple or rec) track block (e.g. for tracks in signal block)
  int            fMCID;     // MC truth ID
  int            fGenIndex; // index in generator block (if matched)

  // -- Momentum
  TVector3       fPlab;
  double         fPtE, fPhiE, fEtaE; 

  // -- Tracking parameters
  int            fQ;                            // charge

  int            fDof, fValidHits, fValidHitFraction, fLayersWithHits;
  double         fChi2;
  int            fTrackQuality; 
  int            fAlgorithm;                    // see DataFormats/TrackReco/interface/TrackBase.h
  int            fHitPattern[40];               // see DataFormats/TrackReco/interface/HitPattern.h

  // -- refitted track parameters (for sig tracks)
  TVector3       fRefPlab; 
  double         fRefChi2; 
  int            fRefDof, fRefValidHits; 

  // verbatim track parameters
  double         fdz, fdzE, fd0, fd0E, fdsz, fdszE, fdxy, fdxyE;  
  // longitudinal/transverse impact parameters wrt PV/SV
  double         fLip, fLipE, fTip, fTipE;
  double         fBsLip, fBsLipE, fBsTip, fBsTipE;

  int            fMuID, fMuIndex; 

  double         fDeDx, fDeDx2; 
  int            fDeDxNmeas, fDeDxNsat; 

  int            fPvIdx; 

  // -- reserve
  int            fInt1,    fInt2,    fInt3; 
  double         fDouble1, fDouble2, fDouble3;

  int            fTag;

private:

  ClassDef(TAnaTrack,1)

};

#endif
