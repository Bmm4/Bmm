#ifndef HFRGENRECOILPRODUCER_H
#define HFRGENRECOILPRODUCER_H

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Bmm/CmsswAnalysis/plugins/HfrBaseProducer.h"

#include <memory>
#include <vector>
#include <set>

// ----------------------------------------------------------------------
// -- Input: gen-level signal decay
//    Algo:  (1) find specified gen-level B decay
//           (2) find OTHER gen-level B candidate (Btag) from same PV
//           (3) create brecoIdx with all simple track indices for OTHER b candidate (Btag)
//           (4) do the same as in HfrSimpleRecoilProducer
// ----------------------------------------------------------------------

class HfrGenRecoilProducer : public HfrBaseProducer {
public:

  void dumpConfiguration();
  explicit HfrGenRecoilProducer(const edm::ParameterSet&);

private:
  void beginRun(const edm::Run&,const edm::EventSetup&) override;
  void endRun(const edm::Run&,const edm::EventSetup&) override;

  void produce(edm::Event&, const edm::EventSetup&) override;
  void put(edm::Event&, std::unique_ptr<std::vector<int> >&, std::unique_ptr<int>&);

  double             fVtxProb;

  int                fMotherID;
  bool               fPartialDecayMatching;
  int                fStableDaughters;
  std::multiset<int> fDaughtersSet;
  std::vector<int>   fDaughtersID;
};
#endif
