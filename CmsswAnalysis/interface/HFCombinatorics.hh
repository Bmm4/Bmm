/*
 *  HFCombinatorics.hh
 *  include tracking information to allow cutting on
 *  - maxdoca
 *  - total charge
 */

#ifndef HFCOMBINATORICS_H
#define HFCOMBINATORICS_H

#include <set>
#include <utility>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/View.h"

class HFCombinatorics {
public:

  class triplet{

  public:
    triplet(int a, int b, int c) : fp1(a), fp2(b), fp3(c) {
    };

    bool operator==(triplet t) {return (fp1 == t.p1() && fp2 == t.p2() && fp3 == t.p3()); };
    int p1(){ return fp1;};
    int p2(){ return fp2;};
    int p3(){ return fp3;};

  private:
    int fp1, fp2, fp3; // particles!
  };

  class doublet{

  public:
    doublet(int a, int b) : fp1(a), fp2(b) {
    };

    bool operator==(doublet t) {return (fp1 == t.p1() && fp2 == t.p2() ); };
    int p1(){ return fp1;};
    int p2(){ return fp2;};

  private:
    int fp1, fp2; // particles!
  };


  HFCombinatorics(edm::Handle<edm::View<reco::Track> > &hTracks,
		     edm::ESHandle<TransientTrackBuilder> &TTB,
		     int verbose) : fhTracks(hTracks), fTTB(TTB),
				    fVerbose(verbose) {}

  // -- Note:
  //    o qCorr = -1: no cut on charges
  //      qCorr =  0: require unlike charges
  //      qCorr =  1: require like charges
  void combine2Tracks(std::vector<doublet> &result,
		      std::vector<int> &tlist1, double mass1, std::vector<int> &tlist2, double mass2,
		      double loMass = 0.4, double hiMass = 2., double maxDoca = 0.1, int qCorr = -1);

  // -- Note:
  //    o cut on charge = +/-1!
  //    o all tracks have same mass
  void combine3Tracks(std::vector<triplet> &result, std::vector<int> &tlist1, double mass = 0.,
		      double loMass = 0.4, double hiMass = 2., double maxDoca = 0.1);

  // -- Note:
  //    o kaon is first, then the two pions!
  void combineDp2Km2Pip(std::vector<triplet> &result, std::vector<int> &tlist1,
			double loMass = 1.6, double hiMass = 2.1, double maxDoca = 0.1);

private:
  edm::Handle<edm::View<reco::Track> > &fhTracks;
  edm::ESHandle<TransientTrackBuilder> fTTB;
  int fVerbose;
};

#endif
