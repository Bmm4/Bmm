#ifndef HFDECAYTREE_H
#define HFDECAYTREE_H

// CMSSW
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "DataFormats/GeometrySurface/interface/ReferenceCounted.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "Bmm/RootAnalysis/rootio/TAnaCand.hh"
#include "Bmm/RootAnalysis/rootio/TAnaVertex.hh"

// STL
#include <set>
#include <map>
#include <vector>
#include <utility>

// ROOT
#include <TVector3.h>


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
struct treeVariables {
  bool valid; 
  // -- this is always zero to allow failing candidates (avoid saving them, e.g. when only computing the mass-constrained version)
  double zero; 
  bool zeroV; 
  // -- variables and flags for node cuts
  void setMass(double x) {mass = x; massV = true;}
  double mass; 
  bool   massV; // valid?
  void setPt(double x) {pt = x; ptV = true;}
  double pt; 
  bool   ptV;
  void setMassE(double x) {masserr = x; masserrV = true;}
  double masserr; 
  bool   masserrV; 
  void setChi2(double x) {chi2 = x; chi2V = true;}
  double chi2; 
  bool   chi2V; 
  void setPvips(double x) {pvips = x; pvipsV = true;}
  double pvips; 
  bool   pvipsV;
  void setMaxDoca(double x) {maxDoca = x; maxDocaV = true;}
  void setMinDoca(double x) {minDoca = x; minDocaV = true;}
  double maxDoca, minDoca; 
  bool   maxDocaV, minDocaV; 
  void setFlxy(double x)  {flxy  = x; flxyV  = true;}
  void setFlsxy(double x) {flsxy = x; flsxyV = true;}
  void setFls3d(double x) {fls3d = x; fls3dV = true;}
  double flxy,  flsxy,  fls3d; 
  bool   flxyV, flsxyV, fls3dV;
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
  HFSimpleCut(double *var, bool *val, double loCut, double hiCut, std::string name = "unset"): fVar(var), fVal(val), fLoCut(loCut), fHiCut(hiCut), fName(name) { }
  int pass() {return ((*fVal)? ((fLoCut < (*fVar)) && ((*fVar) < fHiCut)): -1); }
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

  // HFSimpleCuts
  void addSimpleCut(HFSimpleCut s);
  HFSimpleCut* getSimpleCut(int i);
  unsigned int nSimpleCuts() {return fSimpleCuts.size();}
  bool passSimpleCuts();
  void clearTreeVariables();
  treeVariables fTV; 
  TAnaVertex fAnaVertex;
  
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
  std::vector<HFSimpleCut> fSimpleCuts;
};

#endif
