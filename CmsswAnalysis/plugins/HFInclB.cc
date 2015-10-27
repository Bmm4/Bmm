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

#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"


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

  int printout(0);
  if (printout) cout << "----------------------------------------------------------------------" << endl;

  try {
    HFVirtualDecay::analyze(iEvent,iSetup);
  } catch(HFSetupException e) {
    cout << "==>HFInclB> " << e.fMsg << endl;
    return;
  }

  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  TransverseImpactPointExtrapolator transverseExtrapolator(fMagneticField);


  fListBuilder->setMinPt(fMuonPt); // work with muon pt
  vector<int> muonList = fListBuilder->getMuonList();
  if (muonList.size() < static_cast<unsigned int>(fNmuons)) return; // not enough muons


  if (printout) {
    cout << "npv = " << fVertexCollection.size() 
	 << ": ";
    for (unsigned int i = 0; i < fVertexCollection.size(); ++i) {
      cout << fVertexCollection[i].position().z() << ", ";
    }
    cout << " nmu = " << muonList.size() 
	 << endl;
  }

  // -- loop over muons, build close-track list for each
  int pv0(-2), pv1(-2); 
  double pv0Lip(-2.), pv1Lip(-2.);
  double pv0LipE(-2.), pv1LipE(-2.);
  vector<int> oneMuon; 
  vector<int> trkList;
  for (unsigned int im = 0; im < muonList.size(); ++im) {
    oneMuon.clear(); 
    oneMuon.push_back(muonList[im]);
    fListBuilder->setMinPt(fTrackPt);
    fListBuilder->setMaxDocaToTracks(fMaxDoca);
    fListBuilder->setCloseTracks(&oneMuon);
    trkList.clear();
    trkList = fListBuilder->getTrackList(); 

    pv0 = pv1 = -2; 
    pv0Lip = pv1Lip = -2.;
    pv0LipE = pv1LipE = -2.;

    reco::TrackBaseRef mTrackView(fTracksHandle, muonList[im]);
    reco::Track tM(*mTrackView);
    int pvmidx = getPv(muonList[im], &fVertexCollection); 
    reco::TransientTrack ttM = fTTB->build(tM);

    if (printout) {
      cout << "muon " << muonList[im] << " with pt = " << tM.pt() << " from PV " << pvmidx << endl;
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
    if (trkList.size() > 0) {
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
      if (printout) cout << "--> PV choices: pv0 = " << pv0  << " with lip/lipE = " << pv0Lip << "/" << pv0LipE 
			 << " pv1 = " << pv1  << " with lip/lipE = " << pv1Lip << "/" << pv1LipE 
			 << endl;
    } else {
      if (fVertexCollection.size() < 2) {
	pv0 = 0; 
	if (printout) cout << "--> trkList.size() = " << trkList.size() << " and npv = " << fVertexCollection.size() << ",  simple PV choice: pv0 = " << pv0  << endl;
      } else {
	double distZ(0.), distZmin(99.); 
	int bestPV(-2);
	for (unsigned int iv = 0; iv < fVertexCollection.size(); ++iv) {
	  TransientTrack tTrk = fTTB->build(tM);
	  std::pair<bool, Measurement1D> currentIp = IPTools::signedDecayLength3D(tTrk, GlobalVector(0,0,1), fVertexCollection[iv]);
	  distZ = currentIp.second.value();
	  if (distZ < distZmin) {
	    distZmin = distZ; 
	    bestPV = iv; 
	  }
	}
	if (printout) cout << "--> trkList.size() = " << trkList.size() << ",  npv = " << fVertexCollection.size() << ", chose pv0 = " << bestPV << " with distZ = " << distZmin << endl;
      }
    }

    // -- jet finding with tracks from best vertex plus the close tracks
    vector<PseudoJet> particles;
    particles.clear();
    int already(0); 
    for (unsigned int ix = 0; ix < fTracksHandle->size(); ++ix) {
      int pvidx = getPv(ix, &fVertexCollection); 
      reco::TrackBaseRef tTrackView(fTracksHandle, ix);
      reco::Track tT(*tTrackView);
      already = 0; 
      if (pvidx == pv0) {
	if (ix == static_cast<unsigned int>(muonList[im])) continue; // muon is added below
	if (tT.pt() < 0.3) continue;
	if (tT.pt() > 500.) continue;
	already = 1; 
	PseudoJet a(tT.px(), tT.py(), tT.pz(), TMath::Sqrt(tT.p()*tT.p() + MPION*MPION));
	a.set_user_index(ix); 
	particles.push_back(a);
	if (printout) cout << "  adding to jet particles bestPV trk " << ix << " with pt/eta/phi = " << tT.pt() << "/" << tT.eta() << "/" << tT.phi() << endl;
      }
      if (trkList.end() !=  find(trkList.begin(), trkList.end(), ix)) {
	if (1 == already) {
	  continue;
	}
	if (tT.pt() < 0.3) continue;
	if (tT.pt() > 500.) continue;
	PseudoJet a(tT.px(), tT.py(), tT.pz(), TMath::Sqrt(tT.p()*tT.p() + MPION*MPION));
	a.set_user_index(ix); 
	particles.push_back(a); 
	if (printout) cout << "  adding to jet particles close trk  " << ix << " with pt/eta/phi = " << tT.pt() << "/" << tT.eta() << "/" << tT.phi() << endl;
      }
    }
    // -- add muon
    PseudoJet am(tM.px(), tM.py(), tM.pz(), TMath::Sqrt(tM.p()*tM.p() + MMUON*MMUON));
    am.set_user_index(muonList[im]); 
    if (printout) cout << "  adding to jet particles muon  " << muonList[im] << " with pt/eta/phi = " << tM.pt() << "/" << tM.eta() << "/" << tM.phi() << endl;
    particles.push_back(am); 


    // -- jet finding
    double R(0.5);
    JetDefinition jet_def(antikt_algorithm, R);
    
    ClusterSequence cs(particles, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(1.0));
    
    // -- print the jets
    int ibjet(-1), bJetSize(-1); 
    for (unsigned i = 0; i < jets.size(); i++) {
      vector<PseudoJet> constituents = jets[i].constituents();
      for (unsigned j = 0; j < constituents.size(); j++) {
	if (muonList[im] ==  constituents[j].user_index()) {
	  ibjet = i;
	  bJetSize = constituents.size();
	  break;
	}
      }
      if (ibjet > -1) break;
    }

    if (printout) {
      for (unsigned i = 0; i < jets.size(); i++) {
	vector<PseudoJet> constituents = jets[i].constituents();
	if (ibjet != static_cast<int>(i)) {
	  cout << "   jet " << i << ": "<< jets[i].pt() << " "  << jets[i].pseudorapidity() << " " << jets[i].phi() << endl;
	  for (unsigned j = 0; j < constituents.size(); j++) {
	    cout << "       constituent " << j << " with trk index " << constituents[j].user_index() << " and pt: " << constituents[j].pt() << endl;
	  }
	}
      }
    
      cout << "  bjet " << ibjet << ": "<< jets[ibjet].pt() << " "  << jets[ibjet].pseudorapidity() << " " << jets[ibjet].phi() << endl;
      vector<PseudoJet> constituents = jets[ibjet].constituents();
      for (unsigned j = 0; j < constituents.size(); j++) {
	cout << "       constituent " << j << " with trk index " << constituents[j].user_index() << " and pt: " << constituents[j].pt() << endl;
      }
    }

    if (ibjet < 0) continue;
    if (bJetSize < 2) continue;

    TAnaJet *pTrackJet(0);
    pTrackJet = gHFEvent->addTrackJet();
    pTrackJet->fIndex            = gHFEvent->nTrackJets() - 1;  
    pTrackJet->fPlab.SetPtEtaPhi(jets[ibjet].pt(),
				 jets[ibjet].pseudorapidity(),
				 jets[ibjet].phi()
				 );
    pTrackJet->fQ                = 0.;
    pTrackJet->fE                = jets[ibjet].E();
    pTrackJet->fEt               = jets[ibjet].Et();
    pTrackJet->fM                = jets[ibjet].m();
    pTrackJet->fMt               = jets[ibjet].mt();
    
    TLorentzVector vect;
    vect.SetPtEtaPhiE(jets[ibjet].pt(), jets[ibjet].pseudorapidity(), jets[ibjet].phi(), jets[ibjet].E());
    
    TAnaTrack *pTrack; 
    pTrack           = gHFEvent->addSigTrack();
    pTrack->fInt1    = 100101; // type
    pTrack->fInt2    = pTrackJet->fIndex;
    
    pTrack->fDouble1 = -99.;   // ptrel
    TVector3 a; 
    a.SetPtEtaPhi(tM.pt(), tM.eta(), tM.phi());
    pTrack->fPlab    = a;
    
    TLorentzVector muvect;
    muvect.SetPtEtaPhiM(tM.pt(), tM.eta(), tM.phi(), MMUON);
    vect = vect - muvect;
    pTrack->fDouble1 = muvect.Perp(vect.Vect()); //Ptrel in respct to the corrected jets direction
    if (printout) cout << "--> matched muon with ptrel = " << muvect.Perp(vect.Vect()) << endl;
    
    //define direction
    GlobalVector direction(vect.X(),vect.Y(),vect.Z());
    
    pTrack->fIndex = muonList[im];
    pTrack->fQ     = tM.charge(); 
    TAnaMuon *pM = gHFEvent->getSimpleTrackMuon(muonList[im]); 
    if (pM) {
      pTrack->fMuID = pM->fMuID;
    } else {
      pTrack->fMuID = 0;
    }

  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFInclB);
