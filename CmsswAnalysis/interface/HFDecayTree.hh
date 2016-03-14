#ifndef HFDECAYTREE_H
#define HFDECAYTREE_H

// CMSSW
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "DataFormats/GeometrySurface/interface/ReferenceCounted.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "Bmm/RootAnalysis/rootio/TAnaCand.hh"
#include "Bmm/RootAnalysis/rootio/TAnaVertex.hh"

// STL
#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <utility>

// ROOT
#include <TVector3.h>

class HFDecayTree;
typedef bool (HFDecayTree::*pMemFn)(void);
#define CALL_MEMBER_FN(object,ptrToMember) ((object).*(ptrToMember))

// ----------------------------------------------------------------------
struct track_entry_t {
  track_entry_t(int ix, int pid, bool mfit) : trackIx(ix), particleID(pid), massFit(mfit) {}
  int trackIx;
  int particleID;
  bool massFit;
};
bool operator<(const track_entry_t &t1, const track_entry_t &t2);


// ----------------------------------------------------------------------
struct ImpactParameters {
  ImpactParameters() {
    lip = Measurement1D();
    tip = Measurement1D();
    ip3d = Measurement1D();
  }
  ImpactParameters(Measurement1D plip, Measurement1D ptip, Measurement1D pip3d) {
    lip = plip;
    tip = ptip;
    ip3d = pip3d;
  }
  Measurement1D lip;
  Measurement1D tip;
  Measurement1D ip3d;
};

typedef ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3> > cov33_t;
typedef ROOT::Math::SMatrix<double,6,6,ROOT::Math::MatRepSym<double,6> > cov66_t;
typedef ROOT::Math::SMatrix<double,7,7,ROOT::Math::MatRepSym<double,7> > cov77_t;
typedef ROOT::Math::SMatrix<double,9,9,ROOT::Math::MatRepSym<double,9> > cov99_t;
typedef ROOT::Math::SVector<double,9> jac9_t;


// ----------------------------------------------------------------------
// This struct combines variables
// - used for node cuts (HFSimpleCuts), together with the flag that the variables were filled and are valid (V) to use
// - that depend on two vertices (flight length and pointing angle)
// - that are used in the choice of the PV
class treeVariables {
public:
  treeVariables() {}; 
  ~treeVariables() {}; 
  bool valid; 
  // -- this is always zero to allow failing candidates (avoid saving them, e.g. when only computing the mass-constrained version)
  double zero; 
  bool zeroV; 
  // -- variables and flags for node cuts
  void setMass(double x) {
    mass = x;
    //    std::cout << "setting mass " << x << " for treeVariables = " << this << " &m = " << &mass << std::endl;
  }
  double mass; 
  double massLo, massHi;
  void setPt(double x) {pt = x;}
  double pt; 
  double ptLo, ptHi;
  void setMassE(double x) {masserr = x;}
  double masserr; 
  void setChi2(double x) {chi2 = x;}
  double chi2; 
  double chi2Lo, chi2Hi; 
  void setPvips(double x) {pvips = x;}
  double pvips; 
  double pvipsLo, pvipsHi;
  void setMaxDoca(double x) {maxDoca = x;}
  void setMinDoca(double x) {minDoca = x;}
  double maxDoca, maxDocaLo, maxDocaHi;
  double minDoca, minDocaLo, minDocaHi;
  void setFlxy(double x)  {flxy  = x;}
  void setFlsxy(double x) {flsxy = x;}
  void setFls3d(double x) {fls3d = x;}
  double flxy,  flxyLo,  flxyHi; 
  double flsxy, flsxyLo, flsxyHi; 
  double fls3d, fls3dLo, fls3dHi; 
  int    pvIx, pvIx2;
  // -- PV choice
  ImpactParameters pvImpParams, pvImpParams2nd;   
  double pvLip, pvLip2;   
  double pvLipE, pvLipE2; 
  // -- candidate vertexing
  double  vtxDistanceCosAlphaPlab;
  cov99_t vtxDistanceCov;
  jac9_t  vtxDistanceJac3d, vtxDistanceJac2d;
  double  diffChi2, tau3d, tau3dE, tauxy, tauxyE;
};


// ----------------------------------------------------------------------
class HFSimpleCut {
public:
  HFSimpleCut(double *var, bool *val, double loCut, double hiCut, std::string name = "unset"):
    fVar(var), fVal(val), fLoCut(loCut), fHiCut(hiCut), fName(name) { }
  HFSimpleCut() {clear();}
  int pass() {return ((*fVal)? ((fLoCut < (*fVar)) && ((*fVar) < fHiCut)): -1); }
  void clear() {fVar = 0; fVal = 0; fLoCut = fHiCut = 0.; fName = "unset"; }
  double *fVar;
  bool   *fVal; // valid?
  double fLoCut, fHiCut; 
  std::string fName; 
};


class HFDecayTree;
typedef std::vector<HFDecayTree>::iterator HFDecayTreeIterator;
typedef std::set<track_entry_t>::iterator HFDecayTreeTrackIterator;
// ----------------------------------------------------------------------
class HFDecayTree {
public:
  HFDecayTree(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma = -1.0, bool daughtersToPV = false);
  virtual ~HFDecayTree() { delete fpKinTree;}

  // -- NodeCuts
  void addNodeCut(bool (HFDecayTree::*)(), double lo, double hi, const char *name = ""); 
  bool passAllCuts();
  bool passMass()    {
    // std::cout << "passMass " << fTV.massLo << " < " << fTV.mass << " < " << fTV.massHi
    // 	      << " for tree " << particleID() << " at " << this
    // 	      << std::endl;
    return ((fTV.massLo < fTV.mass) && (fTV.mass < fTV.massHi));
  }
  bool passPt()      {
    // std::cout << "passPt " << fTV.ptLo << " < " << fTV.pt << " < " << fTV.ptHi
    // 	      << " for tree " << particleID() << " at " << this
    // 	      << std::endl;
    return ((fTV.ptLo < fTV.pt) && (fTV.pt < fTV.ptHi));
  }
  bool passMaxDoca() {
    // std::cout << "passMaxDoca " << fTV.maxDocaLo << " < " << fTV.maxDoca << " < " << fTV.maxDocaHi
    // 	      << " for tree " << particleID() << " at " << this
    // 	      << std::endl;
    return ((fTV.maxDocaLo < fTV.maxDoca) && (fTV.maxDoca < fTV.maxDocaHi));
  }
  bool passPvips()   {
    // std::cout << "passPvips " << fTV.pvipsLo << " < " << fTV.pvips << " < " << fTV.pvipsHi
    // 	      << " for tree " << particleID() << " at " << this
    // 	      << std::endl;
    return ((fTV.pvipsLo < fTV.pvips) && (fTV.pvips < fTV.pvipsHi));
  }
  bool passFlsxy()   {
    // std::cout << "passFlsxy " << fTV.flsxyLo << " < " << fTV.flsxy << " < " << fTV.flsxyHi
    // 	      << " for tree " << particleID() << " at " << this
    // 	      << std::endl;
    return ((fTV.flsxyLo < fTV.flsxy) && (fTV.flsxy < fTV.flsxyHi));
  }
  bool passFlxy()    {
    // std::cout << "passPt " << fTV.flxyLo << " < " << fTV.flxy << " < " << fTV.flxyHi
    // 	      << " for tree " << particleID() << " at " << this
    // 	      << std::endl;
    return ((fTV.flxyLo  < fTV.flxy)  && (fTV.flxy  < fTV.flxyHi));
  }
  bool passNever()   {
    // std::cout << "passNever"
    // 	      << " for tree " << particleID() << " at " << this
    // 	      << std::endl;
    return false;
  }

  
  // Constructing the tree structure: add a track with a given type and massFit
  void addTrack(int trackIx, int trackID, bool massFit = true); 

  // To append an already constructed decay tree
  void appendDecayTree(HFDecayTree subTree); 
  // To get a reference to the subvertex
  HFDecayTreeIterator addDecayTree(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma = -1.0, bool daughtersToPV = false);
		
  void clear();
  // Variant to clear and initialize the tree with same signature as constructor
  void clear(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma = -1.0, bool daughtersToPV = false);
		
  // Accessing the track data
  HFDecayTreeTrackIterator getTrackBeginIterator();
  HFDecayTreeTrackIterator getTrackEndIterator();
		
  HFDecayTreeIterator getVerticesBeginIterator();
  HFDecayTreeIterator getVerticesEndIterator();
		
  void getAllTracks(std::vector<track_entry_t> *out_vector, int onlyThisVertex = 0);
  std::vector<track_entry_t> getAllTracks(int onlyThisVertex = 0);
  std::set<int> getAllTracksIndices(int onlyThisVertex = 0);
		
  // Kinematic Tree associated stuff
  std::map<int,int> *getKinParticleMap();
  void setKinParticleMap(std::map<int,int> newMap);
  RefCountedKinematicTree *getKinematicTree();
  void setKinematicTree(RefCountedKinematicTree newTree);
  void resetKinematicTree(int recursive = 0);
		
  // Reconstruction
  TAnaCand *getAnaCand();
  void setAnaCand(TAnaCand *cand);
  
  // Debugging!
  void dump(unsigned indent = 0);

  // Accessors to previously public data members
  bool   vertexing() { return fVertexing; };
  double particleID() { return fParticleID; };
  bool   massConstraint() { return fMassConstraint; };
  double mass() { return fMass; };
  double mass_tracks() { return fMassTracks; }
  double massSigma() { return fMassSigma; };
  double maxDoca() { return fMaxDoca; };
  double minDoca() { return fMinDoca; };
  bool   daughtersToPV() { return fDaughtersToPV; }
		
  void set_vertexing(bool vertexing) { fVertexing = vertexing; };
  void set_particleID(double particleID) { fParticleID = particleID; };
  void set_massConstraint(bool massConstraint) { fMassConstraint = massConstraint; };
  void set_mass(double mass) { fMass = mass; };
  void set_mass_tracks(double mass_tracks) { fMassTracks = mass_tracks; }
  void set_massSigma(double massSigma) { fMassSigma = massSigma; };
  void set_maxDoca(double maxDoca) { fMaxDoca = maxDoca; };
  void set_minDoca(double minDoca) { fMinDoca = minDoca; };
  void set_daughtersToPV(bool daughtersToPV) { fDaughtersToPV = daughtersToPV; }
  
  void clearTreeVariables();
  treeVariables fTV; 
  TAnaVertex fAnaVertex;
  double readBackMass() {return fTV.mass;}
  
private:
  double fParticleID; // if == 0, then no TAnaCandidate should be created.
  bool   fVertexing; // do a vertexing at this node
  double fMass;
  double fMassTracks; // apply a mass constraint to the tracks
  bool   fMassConstraint; // false: no massconstraint at this vertex
  double fMassSigma;
  double fMaxDoca;
  double fMinDoca;
  bool   fDaughtersToPV;

  void dumpTabs(unsigned indent); // used by dump()
		
  std::set<track_entry_t> fTrackIndices; // added tracks
  std::map<int,int> fKinParticleMap; // map: trackIx -> entry in the daughter kinematic particles...
  std::vector<HFDecayTree> fSubVertices;
  RefCountedKinematicTree *fpKinTree;
  TAnaCand *fpAnaCand;
  std::vector<bool (HFDecayTree::*)()> fNodeCuts;
  std::vector<std::string>             fNodeCutNames;
  void        *fTVLocation;
};

#endif
