// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// HFDebug
// ------------
//
// 2016/01/27  Urs Langenegger      first shot
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDebug.h"

#include <sstream>
#include <iostream>
#include <bitset>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "DataFormats/Math/interface/deltaR.h"


#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/MuonReco/interface/MuonEnergy.h"
#include "DataFormats/MuonReco/interface/MuonTime.h"
#include "CondFormats/AlignmentRecord/interface/TrackerSurfaceDeformationRcd.h"



// Transient tracks (for extrapolations)
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"


#include "TH1.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"
#include "TRegexp.h"

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"
#include "Bmm/RootAnalysis/rootio/TAnaTrack.hh"
#include "Bmm/RootAnalysis/rootio/TAnaCand.hh"
#include "Bmm/RootAnalysis/rootio/TGenCand.hh"
#include "Bmm/RootAnalysis/rootio/TAnaVertex.hh"
#include "Bmm/RootAnalysis/rootio/TAnaMuon.hh"
#include "Bmm/RootAnalysis/rootio/TTrgObjv2.hh"



// -- Yikes!
extern TAna01Event *gHFEvent;
extern TFile       *gHFFile;

using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;


// ----------------------------------------------------------------------
HFDebug::HFDebug(const edm::ParameterSet& iConfig):
  HFVirtualDecay(iConfig),
  fTriggerEventLabel(iConfig.getUntrackedParameter<InputTag>("TriggerEventLabel")),
  fTokenTriggerEvent(consumes<TriggerEvent>(fTriggerEventLabel)),

  fHLTResultsLabel(iConfig.getUntrackedParameter<InputTag>("HLTResultsLabel")),
  fTokenTriggerResults(consumes<TriggerResults>(fHLTResultsLabel)),
  fHLTProcessName(iConfig.getUntrackedParameter<string>("HLTProcessName")),
  fTriggerNames(iConfig.getParameter<std::vector<std::string> > ("triggerNames"))
{
  dumpConfiguration();
}


// ----------------------------------------------------------------------
HFDebug::~HFDebug() {
}


// ----------------------------------------------------------------------
void HFDebug::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDebug configuration" << endl;
  cout << "--- HLT process name            : " << fHLTProcessName << endl;
  cout << "--- HLTResultsLabel             : " << fHLTResultsLabel << endl;
  cout << "--- Trigger Event Label         : " << fTriggerEventLabel << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
void HFDebug::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  if (1) cout << "--- HFDebug -------------------------------------------------------------------" << endl;
  float pig = TMath::Pi();
  try {
    HFVirtualDecay::analyze(iEvent,iSetup);
  } catch(HFSetupException e) {
    cout << "==>HFDebug> " << e.fMsg << endl;
    return;
  }

  // -- get trigger information
  Handle<TriggerResults> triggerResults;
  try {
    iEvent.getByToken(fTokenTriggerResults, triggerResults);
    triggerResults.isValid();
  } catch (cms::Exception &ex) {
    if (fVerbose > 0) cout << "==>HFDumpTrigger> Triggerresults  " << fHLTResultsLabel.encode() << " not found " << endl;
    return;
  }

  Handle<trigger::TriggerEvent> triggerEvent;
  try {
    iEvent.getByToken(fTokenTriggerEvent, triggerEvent);
    triggerEvent.isValid();
  } catch (const cms::Exception& e) {
    cout << "Error!! No TriggerEvent with label " << fTriggerEventLabel << endl;
    return;
  }
  const trigger::TriggerObjectCollection triggerObjects = triggerEvent->getObjects();


  iSetup.get<TrackingComponentsRecord>().get("SmartPropagatorAny", fPropagatorAlong);
  iSetup.get<TrackingComponentsRecord>().get("SmartPropagatorAnyOpposite", fPropagatorOpposite);


  // -- HLT L1 muon candidates
  string L1NameCollection("hltL1extraParticles");
  trigger::size_type Index(0);
  Index = triggerEvent->collectionIndex(edm::InputTag(L1NameCollection, "", fTriggerEventLabel.process()));
  vector<TVector3> l1p3;
  TVector3 x;
  if (Index < triggerEvent->sizeCollections()) {
    TString label = TString(L1NameCollection.c_str());
    const trigger::Keys& Keys(triggerEvent->collectionKeys());
    const trigger::size_type n0 (Index == 0? 0 : Keys.at(Index-1));
    const trigger::size_type n1 (Keys.at(Index));
    for (trigger::size_type i = n0; i != n1; ++i) {
      const trigger::TriggerObject& obj(triggerEvent->getObjects().at(i) );
      cout << "===> L1:   "
	   << Form(" %2d id: %+2d m = %4.1f pT/eta/phi = %6.3f/%+5.4f/%+5.4f",
		   i, obj.id(), obj.mass(), obj.pt(), obj.eta(), obj.phi())
	   << " label = " << label
	   << endl;
      x.SetPtEtaPhi(obj.pt(), obj.eta(), obj.phi());
      l1p3.push_back(x);
    }
  }
  // -- HLT L2 muon candidates
  string L2NameCollection("hltL2MuonCandidates");
  Index = triggerEvent->collectionIndex(edm::InputTag(L2NameCollection, "", fTriggerEventLabel.process()));

  if (Index < triggerEvent->sizeCollections()) {
    TString label = TString(L2NameCollection.c_str());
    const trigger::Keys& Keys(triggerEvent->collectionKeys());
    const trigger::size_type n0 (Index == 0? 0 : Keys.at(Index-1));
    const trigger::size_type n1 (Keys.at(Index));
    for (trigger::size_type i = n0; i != n1; ++i) {
      const trigger::TriggerObject& obj( triggerEvent->getObjects().at(i) );
      cout << "===> L2:   "
	   << Form(" %2d id: %+2d m = %4.1f pT/eta/phi = %6.3f/%+5.4f/%+5.4f",
		   i, obj.id(), obj.mass(), obj.pt(), obj.eta(), obj.phi())
	   << " label = " << label
	   << endl;

      x.SetPtEtaPhi(obj.pt(), obj.eta(), obj.phi());
      for (unsigned int i = 0; i < l1p3.size(); ++i) {
	((TH1D*)gHFFile->Get("df2"))->Fill(l1p3[i].DeltaPhi(x));
	((TH1D*)gHFFile->Get("de2"))->Fill(l1p3[i].Eta() - x.Eta());
	((TH1D*)gHFFile->Get("dr2"))->Fill(l1p3[i].DeltaR(x));
      }

    }
  }

  // -- HLT L3 muon candidates
  string L3NameCollection("hltL3MuonCandidates");
  Index = triggerEvent->collectionIndex(edm::InputTag(L3NameCollection, "", fTriggerEventLabel.process()));
  if (Index < triggerEvent->sizeCollections()) {
    TString label = TString(L3NameCollection.c_str());
    const trigger::Keys& Keys(triggerEvent->collectionKeys());
    const trigger::size_type n0 (Index == 0? 0 : Keys.at(Index-1));
    const trigger::size_type n1 (Keys.at(Index));
    for (trigger::size_type i = n0; i != n1; ++i) {
      const trigger::TriggerObject& obj( triggerEvent->getObjects().at(i) );
      cout << "===> L3:   "
	   << Form(" %2d id: %+2d m = %4.1f pT/eta/phi = %6.3f/%+5.4f/%+5.4f",
		   i, obj.id(), obj.mass(), obj.pt(), obj.eta(), obj.phi())
	   << " label = " << label
	   << endl;
    }
  }

  // -- Extrapolation
  double triggerMaxDeltaR(0.1);
  for (MuonCollection::const_iterator imu = fMuonCollection->begin(); imu != fMuonCollection->end(); ++imu) {
    if (imu->pt() < 5) continue;
    bool isSA = (!imu->isGlobalMuon() && imu->isStandAloneMuon() && !imu->isTrackerMuon());
    bool isTR = (!imu->isGlobalMuon() && imu->isTrackerMuon() && !imu->isStandAloneMuon());
    bool isGL = (imu->isGlobalMuon());//&&!(imu->isStandAloneMuon())&&!(imu->isTrackerMuon()));
    bool isTRSA  = (!imu->isGlobalMuon() && imu->isStandAloneMuon()&&imu->isTrackerMuon());

    double matchDeltaR = 9999.;
    int hasTriggered = 0;

    matchDeltaR = matchTrigger(fTriggerIndices, triggerObjects, triggerEvent, (*imu));
    if (matchDeltaR < triggerMaxDeltaR) hasTriggered = 1;
    if (hasTriggered) {
      cout << "muon " << imu-fMuonCollection->begin() << " with pT = " << imu->pt();
      if (isSA)   cout << " STA muon" << endl;
      if (isTR)   cout << " TRA muon" << endl;
      if (isGL)   cout << " GLB muon" << endl;
      if (isTRSA) cout << " STA and TRA muon" << endl;
      cout << " triggered, matchDeltaR = " << matchDeltaR << endl;
    } else {
      cout << " NOT trigger matched"  << endl;
    }

    if (isSA || isGL){
    } else {
      cout << " skipping muon with pt/eta/phi = " << imu->pt() << "/" << imu->eta() << "/" << imu->phi() << endl;
      continue;
    }
    TrackRef tr_mu  = imu->outerTrack();
    if (isSA)   cout << " STA muon ";
    if (isGL)   cout << " GLB muon ";
    cout << " muon pt = " << imu->pt() << endl;

    double ptX(0.), etaX(0.), phiX(0.);
    TrajectoryStateOnSurface tsos = cylExtrapTrkSam(tr_mu, 500);  // track at MB2 radius - extrapolation
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double zz = tsos.globalPosition().z();
      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      double eta0   = tsos.globalPosition().eta();
      double phi0   = tsos.globalPosition().phi();
      double phi(0.);
      if (yy>=0)
	phi = acos(cosphi);
      else
	phi = 2*pig-acos(cosphi);
      // cout << " at MB2: pt = " << imu->pt() << " eta0 = " << eta0 << "phi =  " << phi << " phi0 = " << phi0
      // 	   << " x/y/r = " << xx << " " << yy << " " << rr << " cosphi = " << cosphi
      // 	   << endl;
      if (TMath::Abs(zz) < 790) {
	cout << "===> MB2:  "
	     << Form(" %2d id: %+2d m = %4.1f pT/eta/phi = %6.3f/%+5.4f/%+5.4f, expol: %6.3f/%+5.4f/%+5.4f",
		     1, (isSA?2:1), 0.105, tr_mu->pt(), tr_mu->eta(), tr_mu->phi(), tr_mu->pt(), eta0, phi0)
	     << " phi = " << phi
	     << " z = " << tsos.globalPosition().z()
	     << " r = " << rr
	     << endl;
	etaX = eta0;
	phiX = phi0;
	ptX  = tr_mu->pt();
      }
    }

    tsos = surfExtrapTrkSam(tr_mu, 790);   // track at ME2+ plane - extrapolation
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      double eta0   = tsos.globalPosition().eta();
      double phi0   = tsos.globalPosition().phi();
      double phi(0.);
      if (yy >= 0)
	phi = acos(cosphi);
      else
	phi = 2*pig-acos(cosphi);

      if ((tr_mu->eta() > 0) && (rr < 500)) {
	cout << "===> ME2+: "
	     << Form(" %2d id: %+2d m = %4.1f pT/eta/phi = %6.3f/%+5.4f/%+5.4f, expol: %6.3f/%+5.4f/%+5.4f",
		     1, (isSA?2:1), 0.105, tr_mu->pt(), tr_mu->eta(), tr_mu->phi(), tr_mu->pt(), eta0, phi0)
	     << " phi = " << phi
	     << " z = " << tsos.globalPosition().z()
	     << " r = " << rr
	     << endl;
	etaX = eta0;
	phiX = phi0;
	ptX  = tr_mu->pt();
      }
    }

    tsos = surfExtrapTrkSam(tr_mu, -790);   // track at ME2- plane - extrapolation
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      double eta0   = tsos.globalPosition().eta();
      double phi0   = tsos.globalPosition().phi();
      double phi(0.);
      if (yy>=0)
	phi = acos(cosphi);
      else
	phi = 2*pig-acos(cosphi);
      if ((tr_mu->eta() < 0) && (rr < 500)) {
	cout << "===> ME2-: "
	     << Form(" %2d id: %+2d m = %4.1f pT/eta/phi = %6.3f/%+5.4f/%+5.4f, expol: %6.3f/%+5.4f/%+5.4f",
		     1, (isSA?2:1), 0.105, tr_mu->pt(), tr_mu->eta(), tr_mu->phi(), tr_mu->pt(), eta0, phi0)
	     << " phi = " << phi
	     << " z = " << tsos.globalPosition().z()
	     << " r = " << rr
	     << endl;
	etaX = eta0;
	phiX = phi0;
	ptX  = tr_mu->pt();
      }
    }

    TVector3 xpp3;
    xpp3.SetPtEtaPhi(ptX, etaX, phiX);
    int idx = imu->innerTrack().index();
    for (int im = 0; im < gHFEvent->nMuons(); ++im) {
      TAnaMuon *pm = gHFEvent->getMuon(im);
      TVector3 hfdm;
      if (pm->fIndex == idx) {
	hfdm.SetPtEtaPhi(pm->fOuterPlab.Perp(), pm->fPositionAtM2.Eta(), pm->fPositionAtM2.Phi());
	cout << "     HFDM:                         "
	     << Form(" outer:  %6.3f/%+5.4f/%+5.4f", pm->fOuterPlab.Perp(), pm->fOuterPlab.Eta(), pm->fOuterPlab.Phi())
	     << Form("  expol: %6.3f/%+5.4f/%+5.4f", pm->fOuterPlab.Perp(), pm->fPositionAtM2.Eta(), pm->fPositionAtM2.Phi())
	     << endl;

	for (unsigned int i = 0; i < l1p3.size(); ++i) {
	  ((TH1D*)gHFFile->Get("df0"))->Fill(l1p3[i].DeltaPhi(xpp3));
	  ((TH1D*)gHFFile->Get("df1"))->Fill(l1p3[i].DeltaPhi(hfdm));

	  ((TH1D*)gHFFile->Get("de0"))->Fill(l1p3[i].Eta() - etaX);
	  ((TH1D*)gHFFile->Get("de1"))->Fill(l1p3[i].Eta() - pm->fPositionAtM2.Eta());

	  ((TH1D*)gHFFile->Get("dr0"))->Fill(l1p3[i].DeltaR(xpp3));
	  ((TH1D*)gHFFile->Get("dr1"))->Fill(l1p3[i].DeltaR(hfdm));

	}

      }
    }

  }



}


// ----------------------------------------------------------------------
void  HFDebug::beginRun(const Run &run, const EventSetup &iSetup) {

  bool changed = true;
  if (!fHltConfig.init(run, iSetup, fHLTProcessName, changed)) {
    // if you can't initialize hlt configuration, crash!
    cout << "Error: didn't find process" << fHLTProcessName << endl;
    assert(false);
  }

  bool enableWildcard = true;
  for (size_t iTrig = 0; iTrig < fTriggerNames.size(); ++iTrig) {
    // prepare for regular expression (with wildcards) functionality:
    TString tNameTmp = TString(fTriggerNames[iTrig]);
    TRegexp tNamePattern = TRegexp(tNameTmp, enableWildcard);
    int tIndex = -1;
    // find the trigger index:
    for (unsigned ipath = 0; ipath < fHltConfig.size(); ++ipath) {
      // use TString since it provides reg exp functionality:
      TString tmpName = TString(fHltConfig.triggerName(ipath));
      if (tmpName.Contains(tNamePattern)) {
	tIndex = int(ipath);
	fTriggerIndices.push_back(tIndex);
      }
    }
    if (tIndex < 0) { // if can't find trigger path at all, give warning:
      std::cout << "Warning: Could not find trigger" << fTriggerNames[iTrig] << std::endl;
      //assert(false);
    }
  } // end for triggerNames

  TDirectory *pDir = gDirectory;
  gHFFile->cd();
  TH1D *h1 = new TH1D("dr0", "dr0", 100, 0.0, 0.5);
  h1->SetDirectory(gHFFile);
  h1 = new TH1D("dr1", "dr1", 100,  0.0, 0.5);
  h1->SetDirectory(gHFFile);
  h1 = new TH1D("dr2", "dr2", 100,  0.0, 1.5);
  h1->SetDirectory(gHFFile);
  h1 = new TH1D("de0", "de0", 100, -0.5, 0.5);
  h1->SetDirectory(gHFFile);
  h1 = new TH1D("de1", "de1", 100, -0.5, 0.5);
  h1->SetDirectory(gHFFile);
  h1 = new TH1D("de2", "de2", 100, -0.5, 0.5);
  h1->SetDirectory(gHFFile);
  h1 = new TH1D("df0", "df0", 100, -0.5, 0.5);
  h1->SetDirectory(gHFFile);
  h1 = new TH1D("df1", "df1", 100, -0.5, 0.5);
  h1->SetDirectory(gHFFile);
  h1 = new TH1D("df2", "df2", 100, -1.0, 1.0);
  h1->SetDirectory(gHFFile);
  pDir->cd();

}


// ----------------------------------------------------------------------
void HFDebug::endRun(Run const&, EventSetup const&) {
}


// ----------------------------------------------------------------------
void  HFDebug::beginJob() {
}


// ----------------------------------------------------------------------
void  HFDebug::endJob() {
}


// ----------------------------------------------------------------------
double HFDebug::matchTrigger(vector<int> &trigIndices, const TriggerObjectCollection &trigObjs,
			      Handle<TriggerEvent>  &triggerEvent, const Muon &mu) {
  double matchDeltaR = 9999;

  int iIdx(-1), iObj(-1);
  for(size_t iTrigIndex = 0; iTrigIndex < trigIndices.size(); ++iTrigIndex) {
    int triggerIndex = trigIndices[iTrigIndex];
    const std::vector<std::string> moduleLabels(fHltConfig.moduleLabels(triggerIndex));
    // find index of the last module:
    const unsigned moduleIndex = fHltConfig.size(triggerIndex)-2;
    // find index of HLT trigger name:
    const unsigned hltFilterIndex = triggerEvent->filterIndex(edm::InputTag(moduleLabels[moduleIndex], "", fHLTProcessName));

    if (hltFilterIndex < triggerEvent->sizeFilters()) {
      const trigger::Keys triggerKeys(triggerEvent->filterKeys(hltFilterIndex));
      const trigger::Vids triggerVids(triggerEvent->filterIds(hltFilterIndex));

      const unsigned nTriggers = triggerVids.size();
      for (size_t iTrig = 0; iTrig < nTriggers; ++iTrig) {
        // loop over all trigger objects:
        const trigger::TriggerObject trigObject = trigObjs[triggerKeys[iTrig]];

        double dRtmp = deltaR( mu, trigObject );

        if (dRtmp < matchDeltaR) {
          matchDeltaR = dRtmp;
	  iIdx = iTrigIndex;
	  iObj = iTrig;
        }

      } // loop over different trigger objects
    } // if trigger is in event (should apply hltFilter with used trigger...)
  } // loop over muon candidates
  if (matchDeltaR < 0.1) {
    cout << " matched to iIdx = " << iIdx << " with name = " << fTriggerNames[iIdx] << endl;
    cout << " matched to iObj = " << iObj << endl;
  }
  return matchDeltaR;
}


// ----------------------------------------------------------------------
// to get the track position info at a particular rho
TrajectoryStateOnSurface HFDebug::cylExtrapTrkSam(reco::TrackRef track, double rho) {
  Cylinder::PositionType pos(0, 0, 0);
  Cylinder::RotationType rot;
  Cylinder::CylinderPointer myCylinder = Cylinder::build(pos, rot, rho);
  FreeTrajectoryState recoStart = freeTrajStateMuon(track);
  TrajectoryStateOnSurface recoProp;
  recoProp = fPropagatorAlong->propagate(recoStart, *myCylinder);
  if (!recoProp.isValid()) {
    recoProp = fPropagatorOpposite->propagate(recoStart, *myCylinder);
  }
  return recoProp;
}

// ----------------------------------------------------------------------
// to get track position at a particular (xy) plane given its z
TrajectoryStateOnSurface HFDebug::surfExtrapTrkSam(reco::TrackRef track, double z) {
  Plane::PositionType pos(0, 0, z);
  Plane::RotationType rot;
  Plane::PlanePointer myPlane = Plane::build(pos, rot);
  FreeTrajectoryState recoStart = freeTrajStateMuon(track);
  TrajectoryStateOnSurface recoProp;
  recoProp = fPropagatorAlong->propagate(recoStart, *myPlane);
  if (!recoProp.isValid()) {
    recoProp = fPropagatorOpposite->propagate(recoStart, *myPlane);
  }
  return recoProp;
}


// ----------------------------------------------------------------------
FreeTrajectoryState HFDebug::freeTrajStateMuon(reco::TrackRef track) {
  GlobalPoint  innerPoint(track->innerPosition().x(), track->innerPosition().y(),  track->innerPosition().z());
  GlobalVector innerVec  (track->innerMomentum().x(),  track->innerMomentum().y(),  track->innerMomentum().z());
  FreeTrajectoryState recoStart(innerPoint, innerVec, track->charge(), fMagneticField);
  return recoStart;
}


//define this as a plug-in
DEFINE_FWK_MODULE(HFDebug);
