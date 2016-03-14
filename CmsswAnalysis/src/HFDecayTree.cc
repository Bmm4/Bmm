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
  fpKinTree(0)
{
  if (massConstraint && massSigma <= 0.0) fMassSigma = 0.0001 * mass;
  clearTreeVariables();
}


// ----------------------------------------------------------------------
void HFDecayTree::clearTreeVariables() {
  // for (int i = 0; i < MAXHFSIMPLECUTS; ++i) {
  //   fSimpleCuts[i].clear();
  // }

  fNodeCuts.clear();
  
  fTV.valid       = false; 
  fTV.zero        = 0.;
  
  fTV.mass      = -9999.; 
  fTV.pt        = -9999.; 
  fTV.masserr   = -9999.; 
  fTV.chi2      = -9999.; 
  fTV.pvips     = -9999.; 
  fTV.maxDoca   = -9999.; 
  fTV.minDoca   = -9999.; 
  fTV.flxy   = -9999.;
  fTV.flsxy  = -9999.;
  fTV.fls3d  = -9999.; 
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


// // ----------------------------------------------------------------------
// void HFDecayTree::addSimpleCut(HFSimpleCut a) {
//   if (fNSimpleCuts < MAXHFSIMPLECUTS) {
//     fSimpleCuts[fNSimpleCuts] = a;
//     ++fNSimpleCuts;
//   } else {
//     cout << "ERROR: too many HFSimpleCuts for tree " << particleID() << " at " << this << endl;
//   }
// }


// // ----------------------------------------------------------------------
// HFSimpleCut* HFDecayTree::getSimpleCut(int i) {
//   return &(fSimpleCuts[i]);
// }


// // ----------------------------------------------------------------------
// bool HFDecayTree::passSimpleCuts() {
//   bool result(true); 
//   int cpass(0); 
//   HFSimpleCut *c(0);
//   for (unsigned int i = 0; i < nSimpleCuts(); ++i) {
//     c = getSimpleCut(i); 
//     cpass = c->pass(); 
//     if (1) cout << "HFSimpleCut " << i << " at " << c << ", tree = " << this << " fTV = " << &fTV
// 		<< " name: " << c->fName
// 		<< " " << c->fLoCut << " < " << *(c->fVar) << " < " << c->fHiCut
// 		<< " at fVar = " << c->fVar  
// 		<< ", (*fVal) = " << (*(c->fVal)?"true":"false")
// 		<< " pass: " << cpass
// 		<< endl;
//     if (cpass > -1) {
//       if (0 == cpass) {
// 	result = false;
// 	if (1) cout << "   failed " << c->fLoCut << " < " << *(c->fVar) << " < " << c->fHiCut << endl;
// 	break;
// 	}
//     }
//   }
//   if (0 && result) cout << "-> all passed! " << endl;
//   return result; 
// }


// ----------------------------------------------------------------------
void HFDecayTree::addNodeCut(bool (HFDecayTree::*f)(), double lo, double hi, const char *name) {
  //this does not work here? it did in my toy example??
  //if (&HFDecayTree::passMass == f) {

  // char line[200];
  // sprintf(line, "addNodeCut: f = %p, passMass = %p, passMaxDoca = %p",
  //         (void*)f, (void*)(&HFDecayTree::passMass), (void*)(&HFDecayTree::passMaxDoca)); 
  // cout << line << endl;
  
  if (!strcmp(name, "mass")) {
    //    cout << "passMass setting mass limits" << endl;
    fTV.massLo = lo;
    fTV.massHi = hi;
  } else if (!strcmp(name, "pt")) {
    //    cout << "passPt setting pt limits" << endl;
    fTV.ptLo = lo;
    fTV.ptHi = hi;
  } else if (!strcmp(name, "maxdoca")) {
    //    cout << "passMaxDoca setting maxDoca limits" << endl;
    fTV.maxDocaLo = lo;
    fTV.maxDocaHi = hi;
  } else if (!strcmp(name, "flsxy")) {
    //    cout << "passFlsxy setting limits" << endl;
    fTV.flsxyLo = lo;
    fTV.flsxyHi = hi;
  } else if (!strcmp(name, "fls3d")) {
    //    cout << "passFls3d setting limits" << endl;
    fTV.fls3dLo = lo;
    fTV.fls3dHi = hi;
  } else if (!strcmp(name, "flxy")) {
    //    cout << "passFlxy setting limits" << endl;
    fTV.flxyLo = lo;
    fTV.flxyHi = hi;
  } else if (!strcmp(name, "pvips")) {
    //    cout << "passPvips setting limits" << endl;
    fTV.pvipsLo = lo;
    fTV.pvipsHi = hi;
  } else if (!strcmp(name, "never")) {
    // do nothing
    //    cout << "passNever setting limits" << endl;
  } else {
    cout << "problems determining pointer to member function with name " << name << endl;
  }
  fNodeCutNames.push_back(name); 
  fNodeCuts.push_back(f); 
}


// ----------------------------------------------------------------------
bool HFDecayTree::passAllCuts() {
  bool result(true); 
  for (unsigned int i = 0; i < fNodeCuts.size(); ++i) {
    bool cut = CALL_MEMBER_FN(*this, fNodeCuts[i])();
    //    cout << "passAllCuts " << i << ": " << cut << endl;
    result = result && cut; 
    if (false == result) return false; 
  }
  return result; 
}


// ----------------------------------------------------------------------
void HFDecayTree::dump(unsigned indent) {
  HFDecayTreeIterator treeIt;
  HFDecayTreeTrackIterator trackIt;
  
  dumpTabs(indent);
  cout << "HFDecayTree (particleID = " << fParticleID << "), at " << this << " fTV = " << &fTV
       << ", n(nodecuts) = " << fNodeCuts.size()
       << ", vtx = " << fVertexing
       << ", mconstr = " << fMassConstraint << ", m = " << fMass << ", mSigma = " << fMassSigma 
       << " {" << endl;
  

  for (unsigned int i = 0; i < fNodeCuts.size(); ++i) {
    dumpTabs(indent+1);
    cout << "node cut " << i << ": " << (void*)(fNodeCuts[i]) << "  "; 
    if (fNodeCutNames[i] == "mass") {
      cout << "passMass: " << fTV.massLo << " < " << fTV.mass << " < " << fTV.massHi << endl;
    } else if (fNodeCutNames[i] == "pt") {
      cout << "passPt: " << fTV.ptLo << " < " << fTV.pt << " < " << fTV.ptHi << endl;
    } else if (fNodeCutNames[i] == "maxdoca") {
      cout << "passMaxDoca: " << fTV.maxDocaLo << " < " << fTV.maxDoca << " < " << fTV.maxDocaHi << endl;
    } else if (fNodeCutNames[i] == "flsxy") {
      cout << "passFlsxy: " << fTV.flsxyLo << " < " << fTV.flsxy << " < " << fTV.flsxyHi << endl;
    } else if (fNodeCutNames[i] == "flxy") {
      cout << "passFlxy: " << fTV.flxyLo << " < " << fTV.flxy << " < " << fTV.flxyHi << endl;
    } else if (fNodeCutNames[i] == "pvips") {
      cout << "passPvips: " << fTV.pvipsLo << " < " << fTV.pvips << " < " << fTV.pvipsHi << endl;
    } else if (fNodeCutNames[i] == "never") {
      cout << "passNever: always false" << endl;
    } else {
      cout << endl;
    }
  }
  
  // for (unsigned int i = 0; i < nSimpleCuts(); ++i) {
  //   HFSimpleCut *c = getSimpleCut(i); 
  //   dumpTabs(indent+1);
  //   cout << "HFSimpleCut " << i << " at " << c << ", tree = " << this << " name: " << c->fName
  // 	 << " " << c->fLoCut << " < " << ((*(c->fVal))? Form("%f", *(c->fVar)): "...") << " < " << c->fHiCut
  // 	 << " at fVar = " << c->fVar;
  //   if (*(c->fVal)) {
  //     cout <<  (c->pass() ? " -> true": " -> false");
  //   } else {
  //     cout << " (undefined)";
  //   }
  //   cout << endl;
  // }
  
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
