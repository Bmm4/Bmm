#include <iostream>

#include "Bmm/CmsswAnalysis/interface/HFCombinatorics.hh"

#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

#include "Bmm/RootAnalysis/common/HFMasses.hh"

#include <iostream>
#include <TLorentzVector.h>

using namespace std;
using namespace reco;
using namespace edm;

// ----------------------------------------------------------------------
void HFCombinatorics::combine2Tracks(vector<doublet> &result,
				     vector<int> &tlist1, double t1mass,
				     vector<int> &tlist2, double t2mass,
				     double loMass, double hiMass, double maxDoca, int qCorr) {
  int verbose(fVerbose);
  //  verbose = 1;
  TLorentzVector lv1, lv2, lv;
  double mass(0.), doca(0.);
  TwoTrackMinimumDistance md;

  if (verbose > 0) {
    cout << "----------------------------------------------------------------------" << endl;
    cout << "combine2Tracks: tlist1.size() = " << tlist1.size() << " tlist2.size() = " << tlist2.size() << endl;
  }

  bool sameMass(false);
  if (TMath::Abs(t1mass - t2mass) < 1.e-3) {
    sameMass = true;
  }
  if (verbose) cout << "sameMass = " << sameMass << endl;

  for (unsigned int it1 = 0; it1 < tlist1.size(); ++it1) {
    TrackBaseRef pi1TrackView(fhTracks, tlist1[it1]);
    TransientTrack pi1TT = fTTB->build(*pi1TrackView);
    FreeTrajectoryState fs1 = pi1TT.initialFreeState();;
    lv1.SetXYZM((*fhTracks)[tlist1[it1]].px(), (*fhTracks)[tlist1[it1]].py(), (*fhTracks)[tlist1[it1]].pz(), t1mass);

    for (unsigned int it2 = 0; it2 < tlist2.size(); ++it2) {
      if (tlist1[it1] == tlist2[it2]) {
	if (verbose) cout << Form("%4d %4d: ", tlist1[it1], tlist2[it2]) << "same indices " << endl;
	continue;
      }
      if (sameMass) {
	if (result.end() != find(result.begin(), result.end(), doublet(tlist2[it2], tlist1[it1]))) {
	  if (verbose) cout << Form("%4d %4d: ", tlist1[it1], tlist2[it2]) << "flipped indices for same mass " << endl;
	  continue;
	}
      }

      TrackBaseRef pi2TrackView(fhTracks, tlist2[it2]);
      TransientTrack pi2TT = fTTB->build(*pi2TrackView);
      if (-1 != qCorr) {
	if (0 == qCorr) {
	  if (pi1TT.charge() == pi2TT.charge()) {
	    if (verbose) cout << Form("%4d %4d: ", tlist1[it1], tlist2[it2]) << "qCorr = " << qCorr << " q1 = " << pi1TT.charge() << " q2 = " << pi2TT.charge() << endl;
	    continue;
	  }
	} else if (1 == qCorr) {
	  if (pi1TT.charge() != pi2TT.charge()) {
	    if (verbose) cout << Form("%4d %4d: ", tlist1[it1], tlist2[it2]) << "qCorr = " << qCorr << " q1 = " << pi1TT.charge() << " q2 = " << pi2TT.charge() << endl;
	    continue;
	  }
	}
      }
      FreeTrajectoryState fs2 = pi2TT.initialFreeState();
      lv2.SetXYZM((*fhTracks)[tlist2[it2]].px(), (*fhTracks)[tlist2[it2]].py(), (*fhTracks)[tlist2[it2]].pz(), t2mass);
      // -- cut on mass
      lv = lv1 + lv2;
      mass = lv.M();
      if ((mass < loMass) || (mass > hiMass)) {
	if (verbose) cout << Form("%4d %4d: ", tlist1[it1], tlist2[it2]) << "loMass = " << loMass << " hiMass = " << hiMass << ", but mass = " << mass << endl;
	continue;
      }
      md.calculate(fs1, fs2);
      // -- cut on doca
      doca = md.distance();
      if (doca > maxDoca) {
	if (verbose) cout << Form("%4d %4d: ", tlist1[it1], tlist2[it2]) << "doca = " << doca << " > " << maxDoca << endl;
	continue;
      }
      if (verbose) cout << Form("%4d %4d: ", tlist1[it1], tlist2[it2]) << " adding to doublet vector! mass = " << mass
			 << " doca = " << doca << " q1*q2 = " << pi1TT.charge() << "*" << pi2TT.charge()
			 << endl;
      result.push_back(doublet(tlist1[it1], tlist2[it2]));
    }
  }

  if (verbose) cout << " doublet vector size = " << result.size() << endl;

}


// ----------------------------------------------------------------------
void HFCombinatorics::combineDp2Km2Pip(std::vector<triplet> &result, std::vector<int> &tlist1, double loMass, double hiMass, double maxDoca) {
  int verbose(fVerbose);
  // verbose = 1;
  TLorentzVector lv1, lv2, lv3, lv;
  double mass(0.), doca(0.), doca1(0.), doca2(0.);
  TwoTrackMinimumDistance md;

  if (verbose > 0) {
    cout << "----------------------------------------------------------------------" << endl;
    cout << "combineDp2Km2Pip: tlist1.size() = " << tlist1.size() << endl;
  }

  // -- kaon loop
  for (unsigned int it1 = 0; it1 < tlist1.size(); ++it1) {
    //    if (verbose) cout << "LOOP kaon it1 = " << it1 << " -> tlist1[it1] = " << tlist1[it1] << endl;
    TrackBaseRef kaTrackView(fhTracks, tlist1[it1]);
    TransientTrack kaTT = fTTB->build(*kaTrackView);
    FreeTrajectoryState fs1 = kaTT.initialFreeState();;
    lv1.SetXYZM((*fhTracks)[tlist1[it1]].px(), (*fhTracks)[tlist1[it1]].py(), (*fhTracks)[tlist1[it1]].pz(), MKAON);

    // -- pion1 loop
    for (unsigned int it2 = 0; it2 < tlist1.size(); ++it2) {
      //      if (verbose) cout << "LOOP pion1 it2 = " << it2 << " -> tlist1[it2] = " << tlist1[it2] << endl;
      if (tlist1[it1] == tlist1[it2]) {
	if (verbose) cout << Form("%4d %4d ----: ", tlist1[it1], tlist1[it2]) << "same indices (K, pi1)" << endl;
	continue;
      }

      TrackBaseRef pi1TrackView(fhTracks, tlist1[it2]);
      TransientTrack pi1TT = fTTB->build(*pi1TrackView);
      if (pi1TT.charge() == kaTT.charge()) {
	if (verbose) cout << Form("%4d %4d ----: ", tlist1[it1], tlist1[it2]) << "q(K) = " << kaTT.charge() << " q(pi1) = " << pi1TT.charge() << endl;
	continue;
      }

      FreeTrajectoryState fs2 = pi1TT.initialFreeState();
      md.calculate(fs1, fs2);
      // -- cut on doca
      doca = md.distance();
      if (doca > maxDoca) {
	if (verbose) cout << Form("%4d %4d ----: ", tlist1[it1], tlist1[it2]) << "doca = " << doca << " > " << maxDoca << " for (K,pi1) " << endl;
	continue;
      }
      lv2.SetXYZM((*fhTracks)[tlist1[it2]].px(), (*fhTracks)[tlist1[it2]].py(), (*fhTracks)[tlist1[it2]].pz(), MPION);


      // -- pion 2 loop
      for (unsigned int it3 = it2+1; it3 < tlist1.size(); ++it3) {
	//	if (verbose) cout << "LOOP pion2 it3 = " << it3 << " -> tlist1[it3] = " << tlist1[it3] << endl;
	if (tlist1[it1] == tlist1[it3]) {
	  if (verbose) cout << Form("%4d %4d %4d: ", tlist1[it1], tlist1[it2], tlist1[it3]) << "same indices (ka, pi2) " << endl;
	  continue;
	}

	if (tlist1[it2] == tlist1[it3]) {
	  if (verbose) cout << Form("%4d %4d %4d: ", tlist1[it1], tlist1[it2], tlist1[it3]) << "same indices (pi1, pi2) " << endl;
	  continue;
	}

	TrackBaseRef pi2TrackView(fhTracks, tlist1[it3]);
	TransientTrack pi2TT = fTTB->build(*pi2TrackView);
	if (pi2TT.charge() == kaTT.charge()) {
	  if (verbose) cout << Form("%4d %4d %4d: ", tlist1[it1], tlist1[it2], tlist1[it3])
			    << "q(K) = " << kaTT.charge() << " q(pi2) = " << pi2TT.charge() << endl;
	  continue;
	}

	FreeTrajectoryState fs3 = pi2TT.initialFreeState();
	md.calculate(fs1, fs3);
	// -- cut on doca
	doca1 = md.distance();
	if (doca1 > maxDoca) {
	  if (verbose) cout << Form("%4d %4d %4d: ", tlist1[it1], tlist1[it2], tlist1[it3])
			    << "doca1 = " << doca1 << " > " << maxDoca << " for (K,pi2) " << endl;
	  continue;
	}

	md.calculate(fs2, fs3);
	doca2 = md.distance();
	if (doca2 > maxDoca) {
	  if (verbose) cout << Form("%4d %4d %4d: ", tlist1[it1], tlist1[it2], tlist1[it3])
			    << "doca2 = " << doca2 << " > " << maxDoca << " for (pi1,pi2) " << endl;
	  continue;
	}

	lv3.SetXYZM((*fhTracks)[tlist1[it3]].px(), (*fhTracks)[tlist1[it3]].py(), (*fhTracks)[tlist1[it3]].pz(), MPION);
	// -- cut on mass
	lv = lv1 + lv2 + lv3;
	mass = lv.M();
	if ((mass < loMass) || (mass > hiMass)) {
	  if (verbose) cout << Form("%4d %4d %4d: ", tlist1[it1], tlist1[it2], tlist1[it3])
			    << "loMass = " << loMass << " hiMass = " << hiMass
			    << ", but mass = " << mass
			    << endl;
	  continue;
	}


	if (verbose) cout << Form("%4d %4d %4d: ", tlist1[it1], tlist1[it2], tlist1[it3]) << " adding to triplet vector! mass = " << mass
			  << " doca = " << doca << " doca1 = " << doca1 << " doca2 = " << doca2
			  << " q1*q2*q3 = " << kaTT.charge() << "*" << pi1TT.charge() << "*" << pi2TT.charge()
			  << endl;
	result.push_back(triplet(tlist1[it1], tlist1[it2], tlist1[it3]));

      }
    }
  }

  if (verbose) cout << " triplet vector size = " << result.size() << endl;

}

// ----------------------------------------------------------------------
void HFCombinatorics::combine3Tracks(vector<triplet> &result, vector<int> &tlist1, double tmass,
				 double loMass, double hiMass, double maxDoca) {

  TLorentzVector lv1, lv2, lv3, lv;
  double mass(0.), doca(0.);
  TwoTrackMinimumDistance md;

  for (unsigned int it1 = 0; it1 < tlist1.size(); ++it1) {
    TrackBaseRef pi1TrackView(fhTracks, tlist1[it1]);
    TransientTrack pi1TT = fTTB->build(*pi1TrackView);
    FreeTrajectoryState fs1 = pi1TT.initialFreeState();;
    lv1.SetXYZM((*fhTracks)[tlist1[it1]].px(), (*fhTracks)[tlist1[it1]].py(), (*fhTracks)[tlist1[it1]].pz(), tmass);
    for (unsigned int it2 = it1+1; it2 < tlist1.size(); ++it2) {
      TrackBaseRef pi2TrackView(fhTracks, tlist1[it2]);
      TransientTrack pi2TT = fTTB->build(*pi2TrackView);
      FreeTrajectoryState fs2 = pi2TT.initialFreeState();
      lv2.SetXYZM((*fhTracks)[tlist1[it2]].px(), (*fhTracks)[tlist1[it2]].py(), (*fhTracks)[tlist1[it2]].pz(), tmass);
      md.calculate(fs1, fs2);
      doca = md.distance();
      if (doca > maxDoca) continue;
      for (unsigned int it3 = it2+1; it3 < tlist1.size(); ++it3) {
	TrackBaseRef pi3TrackView(fhTracks, tlist1[it3]);
	TransientTrack pi3TT = fTTB->build(*pi3TrackView);
	FreeTrajectoryState fs3 = pi3TT.initialFreeState();
	md.calculate(fs1, fs3);
	doca = md.distance();
	if (doca > maxDoca) continue;
	md.calculate(fs1, fs2);
	doca = md.distance();
	if (doca > maxDoca) continue;
  	lv3.SetXYZM((*fhTracks)[tlist1[it3]].px(), (*fhTracks)[tlist1[it3]].py(), (*fhTracks)[tlist1[it3]].pz(), tmass);
	int qtot = (*fhTracks)[tlist1[it1]].charge() + (*fhTracks)[tlist1[it2]].charge() + (*fhTracks)[tlist1[it3]].charge();
	if (TMath::Abs(qtot) > 1) continue;
	lv = lv1 + lv2 + lv3;
	mass = lv.M();
	if (mass < loMass) continue;
	if (mass > hiMass) continue;
	result.push_back(triplet(tlist1[it1], tlist1[it2], tlist1[it3]));
      }
    }
  }
}
