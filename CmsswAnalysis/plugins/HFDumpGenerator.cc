#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDumpGenerator.h"

#include <iostream>

#include <TRandom.h>

#include "HepMC/GenVertex.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"

// -- Yikes!
extern TAna01Event *gHFEvent;

using namespace std;
using namespace edm;
using namespace reco;

// ----------------------------------------------------------------------
HFDumpGenerator::HFDumpGenerator(const ParameterSet& iConfig):
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fGenCandidatesLabel(iConfig.getUntrackedParameter<string>("generatorCandidates", string("MCCandidate"))),
  fGenEventLabel(iConfig.getUntrackedParameter<string>("generatorEvent", string("EvtGenProducer")))  {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpGenerator constructor: " << fGenCandidatesLabel << "  " << fGenEventLabel << endl;
  cout << "---  verbose:                         " << fVerbose << endl;
  cout << "---  candidates label:                " << fGenCandidatesLabel << endl;
  cout << "---  event label:                     " << fGenEventLabel << endl;
  cout << "----------------------------------------------------------------------" << endl;

  fTokenGenParticle  = consumes<vector<reco::GenParticle> >(fGenCandidatesLabel);
  fTokenHepMCProduct = consumes<HepMCProduct>(fGenEventLabel);

}


// ----------------------------------------------------------------------
HFDumpGenerator::~HFDumpGenerator() {

}


// ----------------------------------------------------------------------
void HFDumpGenerator::analyze(const Event& iEvent, const EventSetup& iSetup) {

  // -- CAREFUL!
  // keep this Clear()! If a genfilter based on gHFEvent is used in the job, this clear will
  // make sure that the old GenBlock is cleared before filling it again.
  gHFEvent->Clear();

  static int nevt(0);
  ++nevt;

  if (fVerbose > 3) {
    // FIXME: This is broken! AGAIN?!
    // cf http://cmslxr.fnal.gov/source/GeneratorInterface/GenFilters/src/PythiaFilter.cc?v=CMSSW_7_6_2
    cout << "=================HEPMC===================" << endl;
    Handle<HepMCProduct> evt;
    iEvent.getByToken(fTokenHepMCProduct, evt);
    const HepMC::GenEvent *genEvent = evt->GetEvent();
    genEvent->print();

    cout << "=================HEPMC===================" << endl;
  }

  TGenCand  *pGen;
  // -- From PhysicsTools/HepMCCandAlgos/plugins/ParticleListDrawer.cc
  int iMo1(-1), iMo2(-1), iDa1(-1), iDa2(-1);

  vector<const GenParticle *> cands;
  cands.clear();
  vector<const GenParticle *>::const_iterator found = cands.begin();

  Handle<GenParticleCollection> genParticlesH;
  genParticlesH.clear();
  try {
    iEvent.getByToken(fTokenGenParticle, genParticlesH);
  } catch(cms::Exception ce) {
    cout << "==> HFDumpGenerator caught std::exception " << ce.what() << endl;
  }

  for (GenParticleCollection::const_iterator p = genParticlesH->begin(); p != genParticlesH->end(); ++p) {
    cands.push_back( & * p );
  }

  if (fVerbose > 1) cout << Form("Number of genParticles = %i", (int)genParticlesH->size()) << endl;

  int i(-1);
  for(GenParticleCollection::const_iterator p  = genParticlesH->begin(); p != genParticlesH->end();  p++) {
    ++i;
    pGen = gHFEvent->addGenCand();
    pGen->fID     = p->pdgId();
    pGen->fStatus = p->status();
    pGen->fNumber = i;
    pGen->fQ      = p->charge();
    pGen->fMass   = p->mass();
    double vx = p->vx(), vy = p->vy(), vz = p->vz();
    pGen->fP.SetXYZT(p->px(),
		     p->py(),
		     p->pz(),
		     p->energy());
    pGen->fV.SetXYZ(vx, vy, vz);

    // Particles Mothers and Daighters
    iMo1 = -1;
    iMo2 = -1;
    iDa1 = -1;
    iDa2 = -1;
    int nMo = p->numberOfMothers();
    int nDa = p->numberOfDaughters();

    pGen->fDau1 = 99999;
    pGen->fDau2 = -1;

    found = find(cands.begin(), cands.end(), p->mother(0));
    if (found != cands.end()) {
      iMo1 = found - cands.begin();
      pGen->fMom1 = iMo1;
    } else {
      pGen->fMom1 = -1;
    }

    found = find(cands.begin(), cands.end(), p->mother(nMo-1));
    if (found != cands.end()) {
      iMo2 = found - cands.begin();
      pGen->fMom2 = iMo2;
    } else {
      pGen->fMom2 = -1;
    }

    found = find(cands.begin(), cands.end(), p->daughter(0));
    if (found != cands.end()) {
      iDa1 = found - cands.begin();
      pGen->fDau1 = iDa1;
    }

    found = find(cands.begin(), cands.end(), p->daughter(nDa-1));
    if (found != cands.end()) {
      iDa2 = found - cands.begin();
      pGen->fDau2 = iDa2;
    }

    if (fVerbose > 2) pGen->dump();
  }

  genParticlesH.clear(); // WHY?

  if (fVerbose > 0) cout << "==> HFDumpGenerator> Event " << nevt << ", dumped  " << gHFEvent->nGenCands() << " generator cands" << endl;
}


// ------------ method called once each job just before starting event loop  ------------
void  HFDumpGenerator::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDumpGenerator::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpGenerator);
