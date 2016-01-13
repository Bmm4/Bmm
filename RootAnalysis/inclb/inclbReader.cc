#include "inclbReader.hh"
#include <cmath>
#include <string>

#include "common/PidTable.hh"
#include "common/HFMasses.hh"

#include "candAna.hh"

using namespace std;

// ----------------------------------------------------------------------
inclbReader::inclbReader(TChain *tree, TString evtClassName): treeReader01(tree, evtClassName) {
  cout << "==> inclbReader: constructor..." << endl;
  fVerbose = 0; 
}


// ----------------------------------------------------------------------
inclbReader::~inclbReader() {
  cout << "==> inclbReader: destructor..." << endl;
}


// ----------------------------------------------------------------------
void inclbReader::startAnalysis() {
  cout << "==> inclbReader: fVerbose = " << fVerbose << endl;
  fpJSON = new JSON(JSONFILE.c_str(), 1); 
}

// ----------------------------------------------------------------------
void inclbReader::endAnalysis() {
  cout << "==> inclbReader: endAnalysis() destroying cand analysis modules" << endl;
  for (unsigned int i = 0; i < lCandAnalysis.size(); ++i) {
    lCandAnalysis[i]->endAnalysis();
  }
  
}


// ----------------------------------------------------------------------
void inclbReader::eventProcessing() {
  ((TH1D*)fpHistFile->Get("monEvents"))->Fill(0); 

  bool json = false;  
   
  if (fIsMC > 0) {
    json = 1; 
    processTypePythia8(); 
  } else {
    json = fpJSON->good(fRun, fLS); 
    if (fVerbose > 100 && !json) {
      cout << "JSON = 0 for run = " << fRun << " and LS = " << fLS << endl;
    }
    fProcessType = -98; 
  }
 
  // -- call candidate analyses
  for (unsigned int i = 0; i < lCandAnalysis.size(); ++i) {
    lCandAnalysis[i]->fIsMC        = fIsMC;
    lCandAnalysis[i]->fJSON        = json; 
    lCandAnalysis[i]->fRun         = fRun; 
    lCandAnalysis[i]->fEvt         = fEvt; 
    lCandAnalysis[i]->fLS          = fLS; 
    lCandAnalysis[i]->fEvent       = fEvent; 
    lCandAnalysis[i]->fProcessType = fProcessType; 

    lCandAnalysis[i]->evtAnalysis(fpEvt);
  }

}


// ----------------------------------------------------------------------
void inclbReader::bookHist() {
  fpHistFile->cd();
  TH1D *h(0);
  h = new TH1D("monEvents", "monEvents", 20, 0., 20.);

  (void)h; // make compiler warning go away

  for (unsigned int i = 0; i < lCandAnalysis.size(); ++i) {
    lCandAnalysis[i]->bookHist();
  }

}


// ----------------------------------------------------------------------
void inclbReader::readCuts(TString filename, int dump) {
  if (dump) cout << "==> inclbReader: Reading " << filename << " for classes setup" << endl;
  
  ifstream is(filename.Data());
  char buffer[1000]; 
  char className[200], cutFile[200]; 
  while (is.getline(buffer, 1000, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %s", className, cutFile);
    string sclassName(className); 
    // -- set up candidate analyzer classes
    if (string::npos != sclassName.find("candAna")) {
      candAna *a = new candAna(this, className, cutFile); 
      lCandAnalysis.push_back(a); 
    }

    // -- all the rest ...
    if (!strcmp(className, "JSON")) {
      char json[1000];
      sscanf(buffer, "%s %s", className, json);
      JSONFILE = string(json);
      if (dump) cout << "JSON FILE:           " << JSONFILE << endl;
    }
  }

}


// ----------------------------------------------------------------------
// http://home.thep.lu.se/~torbjorn/pythia82html/ParticleProperties.html
void inclbReader::processTypePythia8() {

  TGenCand *pG;
  
  // hard-scatter partons (entries { d, u, s, c, b, t } )
  double hsPartCnt[6];
  double hsAntiCnt[6];

  // partons
  double parPartCnt[6];
  double parAntiCnt[6];    
    
  for (int i = 0; i < 6; i++) {
    hsPartCnt[i] = 0; 
    hsAntiCnt[i] = 0; 
    parPartCnt[i] = 0; 
    parAntiCnt[i] = 0; 
  }

  int aid(0);
  for (int i = 0; i < fpEvt->nGenCands(); ++i) {
    pG = fpEvt->getGenCand(i);
    // -- in PYTHIA8 the hard scatter particles are 21-29
    if (pG->fStatus > 20 && pG->fStatus < 30) {
      aid = TMath::Abs(pG->fID); 
      
      for (int j = 0; j < 6; j++) {
	if (pG->fID == j+1) {  
	  hsPartCnt[j]++;
	}
	if (pG->fID == -(j+1)) {  
	  hsAntiCnt[j]++;
	}
      }
      //      pG->dump(); 
    }
  }

  // -- beauty
  if (hsPartCnt[4] >= 1 && hsAntiCnt[4] >= 1) {
    fProcessType = 40; // gluon fusion
    //    cout << Form("====> b: GGF (%i)", fProcessType) << endl;
    return;
  } 
  
  if ((hsPartCnt[4] >= 1 && hsAntiCnt[4] == 0) || (hsPartCnt[4] == 0 && hsAntiCnt[4] >= 1) ) {
    fProcessType = 41; // flavor excitation
    //    cout << Form("====> b: FEX (%i)", fProcessType) << endl;
    return;
  }

  if (hsPartCnt[4] == 0 && hsAntiCnt[4] == 0) {
    fProcessType = 42; // gluon splitting
    //    cout << Form("====> b: GSP (%i)", fProcessType) << endl;

    return;
  }

  fpEvt->dumpGenBlock();
  
}


// ----------------------------------------------------------------------
void inclbReader::processTypePythia6() {

  TGenCand *pG;
  
  // documentation line partons (entries { d, u, s, c, b, t } )
  double docPartCnt[6];
  double docAntiCnt[6];

  // partons
  double parPartCnt[6];
  double parAntiCnt[6];    
    
  for (int i = 0; i < 6; i++) {
    docPartCnt[i] = 0; 
    docAntiCnt[i] = 0; 
    parPartCnt[i] = 0; 
    parAntiCnt[i] = 0; 
  }

  int aid(0);
  for (int i = 0; i < fpEvt->nGenCands(); ++i) {

    pG = fpEvt->getGenCand(i);

    aid = TMath::Abs(pG->fID); 
    if ( aid == 1 || aid == 2 ||
         aid == 3 || aid == 4 || 
         aid == 5 || aid == 6 || 
         aid == 21) {
      if ( pG->fStatus == 3 ) {
        //      cout << "quark/gluon from documentation #" << i << "(ID: " << pG->fID << ")" << endl;
      }
      if ( pG->fStatus == 2 &&  TMath::Abs(pG->fID) != 21) {
        //      cout << "decayed quark/gluon #" << i << " (ID: " << pG->fID << ")" << endl;
      }
      if ( pG->fStatus == 1 ) {
        //      cout << "undecayed (?) quark/gluon #" << i << " (ID: " << pG->fID  << ")" << endl;
      }
    }

    for (int j = 0; j < 6; j++) {

      if ( pG->fStatus == 3 ) {
        if ( pG->fID == j+1 ) {  
          docPartCnt[j]++;
        }
        if ( pG->fID == -(j+1) ) {  
          docAntiCnt[j]++;
        }
      }

      if ( pG->fStatus == 2 ) {
        if ( pG->fID == j+1 ) {  
          parPartCnt[j]++;
        }
        if ( pG->fID == -(j+1) ) {  
          parAntiCnt[j]++;
        }
      }
    }
  }

  fProcessType = -99;
  // -- top 
  if (docPartCnt[5] >= 1 && docAntiCnt[5] >= 1) {
    fProcessType = 50; 
    //    printf("====> t: GGF (%i)\n", fProcessType);
    return;
  }
  
  if ((docPartCnt[5] >= 1 && docAntiCnt[5] == 0) || (docPartCnt[5] == 0 && docAntiCnt[5] >= 1) ) {
    fProcessType = 51; 
    //    printf("====> t: FEX (%i)\n", fProcessType);
    return;
  }

  if (docPartCnt[5] == 0 && docAntiCnt[5] == 0 && (parPartCnt[5] >= 1 || parAntiCnt[5] >= 1)) {
    fProcessType = 52;
    //    printf("====> t: GSP (%i)\n", fProcessType); 
    return;
  }
  
  // -- beauty
  if (docPartCnt[4] >= 1 && docAntiCnt[4] >= 1) {
    fProcessType = 40; 
   //    printf("====> b: GGF (%i)\n", fProcessType);
    return;
  } 
  
  if ((docPartCnt[4] >= 1 && docAntiCnt[4] == 0) || (docPartCnt[4] == 0 && docAntiCnt[4] >= 1) ) {
    fProcessType = 41; 
    //    printf("====> b: FEX (%i)\n", fProcessType);
    return;
  }

  if (docPartCnt[4] == 0 && docAntiCnt[4] == 0 && (parPartCnt[4] >= 1 || parAntiCnt[4] >= 1)) {
    fProcessType = 42; 
    //    printf("====> b: GSP (%i)\n", fProcessType);

    return;
  }

  if (docPartCnt[3] >= 1 && docAntiCnt[3] >= 1) {
    fProcessType = 30; 
    //    printf("====> c: GGF (%i)\n", fProcessType);
    return;
  }
  
 
  if ((docPartCnt[3] >= 1 && docAntiCnt[3] == 0) || (docPartCnt[3] == 0 && docAntiCnt[3] >= 1) ) {
    fProcessType = 31; 
    //    printf("====> c: FEX (%i)\n", fProcessType);
    return;
  }
  
  if (docPartCnt[3] == 0 && docAntiCnt[3] == 0 && (parPartCnt[3] >= 1 || parAntiCnt[3] >= 1)) {
    fProcessType = 32; 
    //    printf("====> c: GSP (%i)\n", fProcessType);
    return;
  }
  
  // light flavors
  if ((docPartCnt[5] == 0 && docAntiCnt[5] == 0) && (parPartCnt[5] == 0 && parAntiCnt[5] == 0)
      && (docPartCnt[4] == 0 && docAntiCnt[4] == 0) && (parPartCnt[4] == 0 && parAntiCnt[4] == 0)
      && (docPartCnt[3] == 0 && docAntiCnt[3] == 0) && (parPartCnt[3] == 0 && parAntiCnt[3] == 0)
      ) {
    fProcessType = 1; 
    //    printf("====> UDS: light flavors (%i)\n", fProcessType);
    return;
  }

  // if no process type was determined
  //  printf("====> Could not determine process type !!!\n");

  fpEvt->dumpGenBlock();     

}
