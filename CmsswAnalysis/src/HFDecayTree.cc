// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// HFDecayTree
// -----------
//
// 2016/02/28 Urs Langenegger      replace nodeCuts with simpleCuts
// 2011/03/31 Frank Meier          added flag for massConstraint, 
//                                 sign of mass no longer determines behavior
// 2010/04/28 Christoph Naegeli    first shot
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include <iostream>

#include "Bmm/CmsswAnalysis/interface/HFDecayTree.hh"

using namespace std;

// ----------------------------------------------------------------------
// track_entry operator for including in set
bool operator<(const track_entry_t &t1, const track_entry_t &t2) {
  bool result = false;
  if (t1.massFit && !t2.massFit) // only t1 has mass fit
    result = true;
  else if (!t1.massFit && t2.massFit) // only t2 has massFit
    result = false;
  else if (t1.trackIx < t2.trackIx)
    result = true;
  else if (t1.trackIx > t2.trackIx)
    result = false;
  else if (t1.particleID < t2.particleID)
    result = true;
  return result;
}


// ----------------------------------------------------------------------
HFDecayTree::HFDecayTree(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma, bool daughtersToPV) :
  fParticleID(pID),
  fVertexing(doVertexing),
  fMass(mass),
  fMassTracks(0),
  fMassConstraint(massConstraint),
  fMassSigma(massSigma),
  fMaxDoca(0),
  fMinDoca(0),
  fDaughtersToPV(daughtersToPV),
  fpKinTree(0) {
  if (massConstraint && massSigma <= 0.0) fMassSigma = 0.0001 * mass;
  clearTreeVariables(); 
}


// ----------------------------------------------------------------------
void HFDecayTree::clearTreeVariables() {
  fTV.valid       = false; 
  fTV.zero        = 0.;
  fTV.zeroV       = true;
  
  fTV.mass      = -9999.; 
  fTV.massV     = false;
  fTV.pt        = -9999.; 
  fTV.ptV       = false;
  fTV.masserr   = -9999.; 
  fTV.masserrV  = false;
  fTV.chi2      = -9999.; 
  fTV.chi2V     = false;
  fTV.pvips     = -9999.; 
  fTV.pvipsV    = false;
  fTV.maxDoca   = -9999.; 
  fTV.maxDocaV  = false;
  fTV.minDoca   = -9999.; 
  fTV.minDocaV  = false;
  fTV.flxy   = -9999.;
  fTV.flsxy  = -9999.;
  fTV.fls3d  = -9999.; 
  fTV.flxyV  = false;
  fTV.flsxyV = false;
  fTV.fls3dV = false; 
  fTV.pvIx   = -1;
  fTV.pvIx2  = -1; 

  fTV.pvLip    = -9999.;
  fTV.pvLip2   = -9999.;
  fTV.pvLipE   = -9999.;
  fTV.pvLipE2  = -9999.;
  fTV.diffChi2 = -9999.;
  fTV.vtxDistanceCosAlphaPlab = -9999.;

  fTV.tau3d   = -9999.;
  fTV.tau3dE  = -9999.;
  fTV.tauxy   = -9999.;
  fTV.tauxyE  = -9999.;

  
  fTV.pvImpParams.lip     = Measurement1D(-9999.,-9999.);
  fTV.pvImpParams.tip     = Measurement1D(-9999.,-9999.);
  fTV.pvImpParams.ip3d    = Measurement1D(-9999.,-9999.);
  fTV.pvImpParams2nd.lip  = Measurement1D(-9999.,-9999.);   
  fTV.pvImpParams2nd.tip  = Measurement1D(-9999.,-9999.);   
  fTV.pvImpParams2nd.ip3d = Measurement1D(-9999.,-9999.);   

  // FIXME clear them!
  //  cov99_t vtxDistanceCov;
  //  jac9_t  vtxDistanceJac3d, fVtxDistanceJac2d;


}


// ----------------------------------------------------------------------
void HFDecayTree::addTrack(int trackIx, int trackID, bool massFit) {
  fTrackIndices.insert(track_entry_t(trackIx, trackID, massFit));
} 


// ----------------------------------------------------------------------
void HFDecayTree::appendDecayTree(HFDecayTree subTree) {
  fSubVertices.push_back(subTree);
} 


// ----------------------------------------------------------------------
HFDecayTreeIterator HFDecayTree::addDecayTree(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma, bool daughtersToPV) {
  return fSubVertices.insert(fSubVertices.end(), HFDecayTree(pID, doVertexing, mass, massConstraint, massSigma, daughtersToPV));
}


// ----------------------------------------------------------------------
void HFDecayTree::clear() {
  clear(-1, true, 0.0, false, -1.0);
}


// ----------------------------------------------------------------------
void HFDecayTree::clear(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma, bool daughtersToPV) {
  fParticleID = pID;
  fVertexing = doVertexing;
  fMass = mass;
  fMassConstraint = massConstraint;
  fMassSigma = massSigma;
  fMaxDoca = -1.0;
  fMinDoca = -1.0;
  fDaughtersToPV = daughtersToPV;
  fMassTracks = 0.0;
  
  // clear the containers
  fTrackIndices.clear();
  fKinParticleMap.clear();
  fSubVertices.clear();

  // clear the kinematic tree
  delete fpKinTree;
  fpKinTree = NULL;
  fpAnaCand = NULL;

  clearTreeVariables();
  fSimpleCuts.clear();
}


// ----------------------------------------------------------------------
HFDecayTreeTrackIterator HFDecayTree::getTrackBeginIterator() {
  return fTrackIndices.begin();
}


// ----------------------------------------------------------------------
HFDecayTreeTrackIterator HFDecayTree::getTrackEndIterator() {
  return fTrackIndices.end();
} 


// ----------------------------------------------------------------------
HFDecayTreeIterator HFDecayTree::getVerticesBeginIterator() {
  return fSubVertices.begin();
}


// ----------------------------------------------------------------------
HFDecayTreeIterator HFDecayTree::getVerticesEndIterator() {
  return fSubVertices.end();
}


// ----------------------------------------------------------------------
void HFDecayTree::getAllTracks(vector<track_entry_t> *out_vector, int onlyThisVertex) {
  HFDecayTreeTrackIterator trackIt;
  HFDecayTreeIterator treeIt;
  
  for (trackIt = fTrackIndices.begin(); trackIt!=fTrackIndices.end(); ++trackIt)
    out_vector->push_back(*trackIt);
  
  for (treeIt = fSubVertices.begin(); treeIt!=fSubVertices.end(); ++treeIt) {
    if (!treeIt->fVertexing || !onlyThisVertex)
      treeIt->getAllTracks(out_vector,onlyThisVertex);
  }
} 


// ----------------------------------------------------------------------
vector<track_entry_t> HFDecayTree::getAllTracks(int onlyThisVertex) {
  vector<track_entry_t> tracks;
  getAllTracks(&tracks,onlyThisVertex);
  return tracks;
} 


// ----------------------------------------------------------------------
set<int> HFDecayTree::getAllTracksIndices(int onlyThisVertex) {
  vector<track_entry_t> tracks;
  set<int> result;
  getAllTracks(&tracks,onlyThisVertex);
  
  for(vector<track_entry_t >::const_iterator it = tracks.begin(); it != tracks.end();++it)
    result.insert(it->trackIx);
  
  return result;
} 


// ----------------------------------------------------------------------
map<int,int> *HFDecayTree::getKinParticleMap() {
  return &fKinParticleMap;
}


// ----------------------------------------------------------------------
void HFDecayTree::setKinParticleMap(map<int,int> newMap) {
  fKinParticleMap = newMap;
} 


// ----------------------------------------------------------------------
RefCountedKinematicTree* HFDecayTree::getKinematicTree() {
  return fpKinTree;
}


// ----------------------------------------------------------------------
void HFDecayTree::setKinematicTree(RefCountedKinematicTree newTree) {
  if(!fpKinTree) fpKinTree = new RefCountedKinematicTree;
  *fpKinTree = newTree; // make a copy from the reference counting pointer
} 


// ----------------------------------------------------------------------
void HFDecayTree::resetKinematicTree(int recursive) {
  HFDecayTreeIterator treeIt;
  
  if(recursive) {
    for(treeIt = getVerticesBeginIterator(); treeIt!=getVerticesEndIterator(); ++treeIt)
      treeIt->resetKinematicTree(recursive);
  }
  
  fKinParticleMap.clear();
  delete fpKinTree;
  fpKinTree = NULL;
  fpAnaCand = NULL;
} 


// ----------------------------------------------------------------------
TAnaCand *HFDecayTree::getAnaCand() {
  return fpAnaCand;
}


// ----------------------------------------------------------------------
void HFDecayTree::setAnaCand(TAnaCand *cand) {
  fpAnaCand = cand;
} 


// ----------------------------------------------------------------------
void HFDecayTree::addSimpleCut(HFSimpleCut a) {
  fSimpleCuts.push_back(a); 
}


// ----------------------------------------------------------------------
HFSimpleCut* HFDecayTree::getSimpleCut(int i) {
  return &(fSimpleCuts[i]);
}


// ----------------------------------------------------------------------
bool HFDecayTree::passSimpleCuts() {
  bool result(true); 
  int cpass(0); 
  HFSimpleCut *c(0);
  for (unsigned int i = 0; i < nSimpleCuts(); ++i) {
    c = getSimpleCut(i); 
    cpass = c->pass(); 
    if (0) cout << "HFSimpleCut " << i << " at " << c << ", tree = " << this << " name: " << c->fName
		<< " " << c->fLoCut << " < " << *(c->fVar) << " < " << c->fHiCut
		<< " at fVar = " << c->fVar
		<< ", (*fVal) = " << *(c->fVal)
		<< " pass: " << cpass
		<< endl;
    if (cpass > -1) {
      if (0 == cpass) {
	result = false;
	if (0) cout << "   failed " << c->fLoCut << " < " << *(c->fVar) << " < " << c->fHiCut << endl;
	break;
	}
    }
  }
  if (0 && result) cout << "-> all passed! " << endl;
  return result; 
}



// ----------------------------------------------------------------------
void HFDecayTree::dump(unsigned indent) {
  HFDecayTreeIterator treeIt;
  HFDecayTreeTrackIterator trackIt;
  
  dumpTabs(indent);
  cout << "HFDecayTree (particleID = " << fParticleID << "), at " << this
       << ", vertexing = " << fVertexing
       << ", massConstraint = " << fMassConstraint << ", mass = " << fMass << ", massSigma = " << fMassSigma 
       << " {" << endl;
  

  for (unsigned int i = 0; i < nSimpleCuts(); ++i) {
    HFSimpleCut *c = getSimpleCut(i); 
    dumpTabs(indent+1);
    cout << "HFSimpleCut " << i << " at " << c << ", tree = " << this << " name: " << c->fName
	 << " " << c->fLoCut << " < ... < " << c->fHiCut
	 << " at fVar = " << c->fVar
	 << endl;
  }
  
  for (trackIt = fTrackIndices.begin(); trackIt!=fTrackIndices.end(); ++trackIt) {
    dumpTabs(indent+1);
    cout << "trackIx = " << trackIt->trackIx << ", trackParticleID = " << trackIt->particleID << ", massFit = " << trackIt->massFit << endl;
  }
  
  for (treeIt = fSubVertices.begin(); treeIt != fSubVertices.end(); ++treeIt)
    treeIt->dump(indent+1);
  
  dumpTabs(indent);
  cout << '}' << endl;
} 


// ----------------------------------------------------------------------
void HFDecayTree::dumpTabs(unsigned indent) {
  for (unsigned j = 0; j < indent; j++) cout << '\t';
} 
