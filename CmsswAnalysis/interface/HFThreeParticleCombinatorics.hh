#ifndef GUARD_HFTHREEPARTICLECOMBINATORICS_H
#define GUARD_HFTHREEPARTICLECOMBINATORICS_H

#include <vector>
#include <utility>

#include <TLorentzVector.h>

class triplet{

public:
  triplet(int a, int b, int c) : kaon(a) {       // first one is always kaon
     if(b<c){                                    // pions are ordered
           pion1=b; pion2=c;
     } else {
	   pion1=c; pion2=b;}
  };

  bool operator==(triplet t) {return (kaon==t.ka() && pion1==t.pi1() && pion2==t.pi2()); };
  bool different(triplet a){
    if (pion1 == a.kaon) return false;
    if (pion1 == a.pion1) return false;
    if (pion1 == a.pion2) return false;
    if (pion2 == a.kaon) return false;
    if (pion2 == a.pion1) return false;
    if (pion2 == a.pion2) return false;
    if (kaon == a.kaon) return false;
    if (kaon == a.pion1) return false;
    if (kaon == a.pion2) return false;
    return true;
  };
  int pi1(){ return pion1;};
  int pi2(){ return pion2;};
  int ka(){ return kaon;};
  int pi3(){ return kaon;};

private:
   int kaon;
   int pion1;
   int pion2;
};

// ----------------------------------------------------------------------
class HFThreeParticleCombinatorics {

public:

  HFThreeParticleCombinatorics(int verbose = 0);
  ~HFThreeParticleCombinatorics();

// combining particle from KList with two particles from piList into a combination
// with mass between loMass and hiMass
//
  void combine(std::vector<triplet > &combList,
	       std::vector<std::pair<int, TLorentzVector> > &KList,
	       std::vector<std::pair<int, TLorentzVector> > &piList,
	       double loMass = 0.2, double hiMass = 0.8);

// combine three particles from piList into a combination with mass between loMass and hiMass
// if loResMass and hiResMass >0, the first two particles form a resonance
// with mass between loResMass and hiResMass
//
  void combine(std::vector<triplet > &combList,
	       std::vector<std::pair<int, TLorentzVector> > &piList,
	       double loMass = 0.4, double hiMass = 20., double loResMass = 0., double hiResMass = 0.);

  int fVerbose;
};

#endif
