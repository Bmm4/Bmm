#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFInclB.h"

#include <iostream>
#include <utility>


#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"


#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"
#include "Bmm/RootAnalysis/rootio/TAnaTrack.hh"
#include "Bmm/RootAnalysis/rootio/TAnaCand.hh"
#include "Bmm/RootAnalysis/rootio/TGenCand.hh"

#include "Bmm/CmsswAnalysis/interface/HFDumpUtilities.hh"
#include "Bmm/CmsswAnalysis/interface/HFSequentialVertexFit.hh"
#include "Bmm/CmsswAnalysis/interface/HFTwoParticleCombinatoricsNew.hh"
#include "Bmm/RootAnalysis/common/HFMasses.hh"
#include "Bmm/CmsswAnalysis/interface/HFDecayTree.hh"
#include "Bmm/CmsswAnalysis/interface/HFTrackListBuilder.hh"

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"


// -- Yikes!
extern TAna01Event *gHFEvent;

using namespace std;
using namespace reco;
using namespace edm;
using namespace fastjet;


typedef boost::shared_ptr<fastjet::ClusterSequence>  ClusterSequencePtr;
typedef boost::shared_ptr<fastjet::JetDefinition>    JetDefPtr;

// ----------------------------------------------------------------------
HFInclB::HFInclB(const ParameterSet& iConfig) :
  HFVirtualDecay(iConfig),
  fNmuons(iConfig.getUntrackedParameter<int>("nMuons", 1)) {
  dumpConfiguration();
}

// ----------------------------------------------------------------------
void HFInclB::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFInclB configuration" << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "---  nMuons:                   " << fNmuons << endl;
  cout << "----------------------------------------------------------------------" << endl;
} 

// ----------------------------------------------------------------------
void HFInclB::analyze(const Event& iEvent, const EventSetup& iSetup) {

  int printout(1);

  cout << "----------------------------------------------------------------------" << endl;

  try {
    HFVirtualDecay::analyze(iEvent,iSetup);
  } catch(HFSetupException e) {
    cout << "==>HFInclB> " << e.fMsg << endl;
    return;
  }

  fListBuilder->setMinPt(fMuonPt); // work with muon pt
  vector<int> muonList = fListBuilder->getMuonList();
  if (muonList.size() < static_cast<unsigned int>(fNmuons)) return; // not enough muons
	
  // -- loop over muons, build close-track list for each
  int pv0(-2), pv1(-2); 
  double pv0Lip(-2.), pv1Lip(-2.);
  double pv0LipE(-2.), pv1LipE(-2.);
  vector<int> oneMuon; 
  for (unsigned int im = 0; im < muonList.size(); ++im) {
    oneMuon.clear(); 
    oneMuon.push_back(muonList[im]);
    fListBuilder->setMinPt(fTrackPt);
    fListBuilder->setMaxDocaToTracks(fMaxDoca);
    fListBuilder->setCloseTracks(&oneMuon);
    vector<int> trkList = fListBuilder->getTrackList(); 

    pv0 = pv1 = -2; 
    pv0Lip = pv1Lip = -2.;
    pv0LipE = pv1LipE = -2.;

    reco::TrackBaseRef mTrackView(fTracksHandle, muonList[im]);
    reco::Track tM(*mTrackView);
    if (printout) {
      int pvmidx = getPv(muonList[im], &fVertexCollection); 
      cout << "muon " << muonList[im] << " with pt = " << tM.pt() << " from PV " << pvmidx << endl;
      
      reco::TransientTrack ttM = fTTB->build(tM);
      TwoTrackMinimumDistance md;
      
      for (unsigned int it = 0; it < trkList.size(); ++it) {
	if (trkList[it] == muonList[im]) continue; // skip muon
	reco::TrackBaseRef tTrackView(fTracksHandle, trkList[it]);
	reco::Track tT(*tTrackView);
	
	reco::TransientTrack ttT = fTTB->build(tT);
	md.calculate(ttM.initialFreeState(),ttT.initialFreeState());
	
	
	int pvidx = getPv(trkList[it], &fVertexCollection); 
	
	cout << " track " << trkList[it] << " with pt = " << tT.pt() << " and doca(m,t) = " << md.distance() 
	     << " from PV " << pvidx
	     << endl;
      }
    }

    // -- vertex these tracks to obtain a 'best' PV
    HFDecayTree theTree(fType, true, 0, false);
    theTree.addTrack(muonList[im], 13);
    for (unsigned int it = 0; it < trkList.size(); ++it) {
      if (trkList[it] == muonList[im]) continue; // skip muon
      theTree.addTrack(trkList[it], 11);
    }
    fSequentialFitter->doFit(&theTree);
    fSequentialFitter->pvIndices(pv0, pv1); 
    fSequentialFitter->pvLip(pv0Lip, pv1Lip); 
    fSequentialFitter->pvLipE(pv0LipE, pv1LipE); 
    
    cout << "--> PV choices: pv0 = " << pv0  << " with lip/lipE = " << pv0Lip << "/" << pv0LipE 
	 << " pv1 = " << pv1  << " with lip/lipE = " << pv1Lip << "/" << pv1LipE 
	 << endl;

    // -- jet finding with tracks from best vertex plus the close tracks
    vector<PseudoJet> particles;
    int already(0); 
    for (unsigned int ix = 0; ix < fTracksHandle->size(); ix++) {
      int pvidx = getPv(ix, &fVertexCollection); 
      reco::TrackBaseRef tTrackView(fTracksHandle, ix);
      reco::Track tT(*tTrackView);
      already = 0; 
      if (pvidx == pv0) {
	if (ix == muonList[im]) continue; // muon is added below
	if (tT.pt() < 0.5) continue;
	if (tT.pt() > 500.) continue;
	already = 1; 
	PseudoJet a(tT.px(), tT.py(), tT.px(), tT.p());
	a.set_user_index(ix); 
	particles.push_back(a);
	cout << "  adding bestPV trk " << ix << endl;
      }
      if (trkList.end() !=  find(trkList.begin(), trkList.end(), ix)) {
	if (1 == already) {
	  continue;
	}
	if (tT.pt() < 0.5) continue;
	if (tT.pt() > 500.) continue;
	PseudoJet a(tT.px(), tT.py(), tT.px(), tT.p());
	a.set_user_index(ix); 
	particles.push_back(a); 
	cout << "  adding close trk  " << ix << endl;
      }
    }
    // -- add muon
    PseudoJet am(tM.px(), tM.py(), tM.px(), tM.p());
    am.set_user_index(muonList[im]); 
    particles.push_back(am); 


    // choose a jet definition
    double R(0.5);
    JetDefinition jet_def(antikt_algorithm, R);
    
    // run the clustering, extract the jets
    ClusterSequence cs(particles, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    
    // print the jets
    bool bjet(false); 
    for (unsigned i = 0; i < jets.size(); i++) {
      bjet = false; 
      vector<PseudoJet> constituents = jets[i].constituents();
      for (unsigned j = 0; j < constituents.size(); j++) {
	if (muonList[im] ==  constituents[j].user_index()) {
	  bjet = true; 
	  break;
	}
      }
      cout << (bjet?"  bjet ": "   jet ") << i << ": "<< jets[i].pt() << " "  << jets[i].rap() << " " << jets[i].phi() << endl;
      for (unsigned j = 0; j < constituents.size(); j++) {
	cout << "       constituent " << j << " with trk index " << constituents[j].user_index() << " and pt: " << constituents[j].pt() << endl;
      }
    }
    

  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFInclB);
