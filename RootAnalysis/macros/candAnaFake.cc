#include "candAnaFake.hh"
#include <cmath>
#include <string>

#include "common/HFMasses.hh"


using namespace std;


// ----------------------------------------------------------------------
candAnaFake::candAnaFake(bmmReader *pReader, std::string name, std::string cutsFile) : 
  candAna(pReader, name, cutsFile) {
  cout << "==> candAnaFake: name = " << name << ", reading cutsfile " << cutsFile << endl;
  BLIND = 0;
  readCuts(cutsFile, 1); 

}


// ----------------------------------------------------------------------
candAnaFake::~candAnaFake() {
  cout << "==> candAnaFake: destructor..." << endl;
}


// ----------------------------------------------------------------------
bool candAnaFake::anaMC(TAna01Event *evt) {
  return true;
}


//------------------------------------------------------------------------------------------
void candAnaFake::evtAnalysis(TAna01Event *evt) {
  candAna::evtAnalysis(evt);
  return;
}


// ----------------------------------------------------------------------
void candAnaFake::candAnalysis() {
  if (0 == fpCand) {
    cout << "No candidate with type = " << TYPE << " found, skipping" << endl;
    return;
  }

  candAna::candAnalysis(); 

  bool skip(false); 
  if (!fJSON) skip = true; 
  if (fCandChi2 > 3) skip = true; 
  if (fCandDoca > 0.005) skip = true; 
  if (!fGoodTracks) skip = true;
  if (!fGoodQ) skip = true; 
  if (fMu1Pt <  3.5) skip = true;
  if (fMu2Pt <  3.5) skip = true;
  // FIXME: These should be replaced with the proper signed 3D impact parameter (to be filled!)
  //  if (fpMuon1->fLip/fpMuon1->fLipE < 1.) skip = true; 
  //  if (fpMuon2->fLip/fpMuon2->fLipE < 1.) skip = true; 
  
  if (skip) return; 

  if (310 == TRUTHCAND) {
    // root [172] events->Draw("m", "0.45<m&&m<0.55&&fls3d>15&&maxdoca<0.005&&flxy<4&&pvips<5&&pvlips<1&&alpha<0.008", "samep")
    if (fCandFLxy > 4.0)      skip = true; 
    if (fCandFLSxy < 15.0)    skip = true; 
    if (fCandFLS3d < 15.0)    skip = true; 
    if (fCandA > 0.1)         skip = true; 
    if (fCandDoca > 0.005)    skip = true; 
    if (fCandPvIpS > 5.0)     skip = true; 
    if (fCandPvLipS > 1.0)    skip = true;
    // -- calculate Lambda mass
    TLorentzVector pr, pi, la;
    pr.SetPtEtaPhiM(fpMuon1->fPlab.Perp(), fpMuon1->fPlab.Eta(), fpMuon1->fPlab.Phi(), MPROTON); 
    pi.SetPtEtaPhiM(fpMuon2->fPlab.Perp(), fpMuon2->fPlab.Eta(), fpMuon2->fPlab.Phi(), MPION);
    la = pr + pi;
    double mla = la.M();
    if (TMath::Abs(mla - MLAMBDA_0) < 0.006) skip = true; 
  } else if (333 == TRUTHCAND) {
  } else if (3122 == TRUTHCAND) {
    // root [227] events->Draw("m", "1.0<m&&m<1.2&&fls3d>15&&flxy<4&&maxdoca<0.005&&pvips<5", "")
    if (fCandFLxy > 4.0)   skip = true; 
    if (fCandFLSxy < 15)   skip = true; 
    if (fCandFLS3d < 15)   skip = true;  
    if (fCandDoca > 0.005) skip = true; 
    if (fCandPvIpS > 5.0)  skip = true; 

    // -- calculate KS masses
    TLorentzVector pi1, pi2, ks;
    pi1.SetPtEtaPhiM(fpMuon1->fPlab.Perp(), fpMuon1->fPlab.Eta(), fpMuon1->fPlab.Phi(), MPION); 
    pi2.SetPtEtaPhiM(fpMuon2->fPlab.Perp(), fpMuon2->fPlab.Eta(), fpMuon2->fPlab.Phi(), MPION);
    ks = pi1 + pi2;
    double mks = ks.M();
    if (TMath::Abs(mks - MKSHORT) < 0.025) skip = true; 
    pi2.SetPtEtaPhiM(fpMuon1->fPlab.Perp(), fpMuon1->fPlab.Eta(), fpMuon1->fPlab.Phi(), MPION); 
    pi1.SetPtEtaPhiM(fpMuon2->fPlab.Perp(), fpMuon2->fPlab.Eta(), fpMuon2->fPlab.Phi(), MPION);
    ks = pi1 + pi2;
    mks = ks.M();
    if (TMath::Abs(mks - MKSHORT) < 0.025) skip = true; 

    // -- proton has higher pT than pion
    if (211 == TMath::Abs(fpMuon1->fMCID)) skip = true; 
  }
  if (skip) return; 

  ((TH1D*)fHistDir->Get("mall"))->Fill(fpCand->fMass); 
  // cout << "MUBDT = " << MUBDT << " fMu1rMvaId = " << fMu1rMvaId << " fMu1rBDT = " << fMu1rBDT
  //      << " fMu2rMvaId = " << fMu2rMvaId << " fMu2rBDT = " << fMu2rBDT
  //      << endl;

  // -- leg 1
  if (fMu1rMvaId > MUBDT) {
    ((TH1D*)fHistDir->Get("pt1muo"))->Fill(fMu1Pt, fpCand->fMass); 
    ((TH1D*)fHistDir->Get("pt1all"))->Fill(fMu1Pt, fpCand->fMass); 
  } else {
    ((TH1D*)fHistDir->Get("pt1nmu"))->Fill(fMu1Pt, fpCand->fMass); 
    ((TH1D*)fHistDir->Get("pt1all"))->Fill(fMu1Pt, fpCand->fMass); 
  }

  // -- leg 2
  if (fMu2rMvaId > MUBDT) {
    ((TH1D*)fHistDir->Get("pt2muo"))->Fill(fMu2Pt, fpCand->fMass); 
    ((TH1D*)fHistDir->Get("pt2all"))->Fill(fMu2Pt, fpCand->fMass); 
  } else {
    ((TH1D*)fHistDir->Get("pt2nmu"))->Fill(fMu2Pt, fpCand->fMass); 
    ((TH1D*)fHistDir->Get("pt2all"))->Fill(fMu2Pt, fpCand->fMass); 
  }
  
  return;
}

// ----------------------------------------------------------------------
void candAnaFake::dumpHFTruthCand(TAnaCand *pC) {
  TAnaTrack *pT(0); 
  if( pC->fSig1 == -1 && pC->fSig2==-1 ) return;
  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getRecTrack(fpEvt->getSigTrack(i)->fIndex); 
    pT->dump(); 
  }
}

// ----------------------------------------------------------------------
int candAnaFake::truthMatch(TAnaCand *pCand, int verbose) {
  return 0; 
}
  

// ----------------------------------------------------------------------
void candAnaFake::bookHist() {
  cout << "==>candAnaFake: bookHist" << endl;

  candAna::bookHist();

  const double xbins[] = {0., 4., 6., 10., 15., 25., 50.};
  const int nbins(6); 
  
  fHistDir->cd();  
  string name("K_{S}");
  double xmin(0.45), xmax(0.55);
  double nbin(50); 
  if (310 == TRUTHCAND) {
    new TH1D("mall", name.c_str(), nbin, xmin, xmax);
    new TH2D("pt1all", "pt1 pion (all)", 6, xbins, nbin, xmin, xmax);
    new TH2D("pt1muo", "pt1 pion (muon)", 6, xbins, nbin, xmin, xmax);
    new TH2D("pt1nmu", "pt1 pion (non-muon)", 6, xbins, nbin, xmin, xmax);
    new TH2D("pt2all", "pt2 pion (all)", nbins, xbins, nbin, xmin, xmax);
    new TH2D("pt2muo", "pt2 pion (muon)", nbins, xbins, nbin, xmin, xmax);
    new TH2D("pt2nmu", "pt2 pion (non-muon)", nbins, xbins, nbin, xmin, xmax);
  } else if (333 == TRUTHCAND) {
    name = "#phi";
    xmin = 0.98;
    xmax = 1.06;
    new TH1D("mall", name.c_str(), nbin, xmin, xmax);
    new TH2D("pt1all", "pt1 kaon (all)", nbins, xbins, nbin, xmin, xmax);
    new TH2D("pt1muo", "pt1 kaon (muon)", nbins, xbins, nbin, xmin, xmax);
    new TH2D("pt1nmu", "pt1 kaon (non-muon)", nbins, xbins, nbin, xmin, xmax);
    new TH2D("pt2all", "pt2 kaon (all)", nbins, xbins, nbin, xmin, xmax);
    new TH2D("pt2muo", "pt2 kaon (muon)", nbins, xbins, nbin, xmin, xmax);
    new TH2D("pt2nmu", "pt2 kaon (non-muon)", nbins, xbins, nbin, xmin, xmax);
  } else if (3122 == TRUTHCAND) {
    name = "#Lambda";
    xmin = 1.05;
    xmax = 1.20;
    new TH1D("mall", name.c_str(), nbin, xmin, xmax);
    new TH2D("pt1all", "pt1 proton (all)", nbins, xbins, nbin, xmin, xmax);
    new TH2D("pt1muo", "pt1 proton (muon)", nbins, xbins, 50,  1.0, 1.2);
    new TH2D("pt1nmu", "pt1 proton (non-muon)", nbins, xbins, nbin, xmin, xmax);
    new TH2D("pt2all", "pt2 pion (all)", nbins, xbins, nbin, xmin, xmax);
    new TH2D("pt2muo", "pt2 pion (muon)", nbins, xbins, 50,  1.0, 1.2);
    new TH2D("pt2nmu", "pt2 pion (non-muon)", nbins, xbins, nbin, xmin, xmax);
  } 

}

// ----------------------------------------------------------------------
void candAnaFake::readCuts(string filename, int dump) {
  candAna::readCuts(filename, dump); 

  fCutFile = filename;

  if (dump) cout << "==> candAnaFake: Reading " << fCutFile << " for cut settings" << endl;
  vector<string> cutLines; 
  readFile(fCutFile, cutLines);

  char CutName[100];
  float CutValue;

  char  buffer[200];
  fHistDir->cd();
  TH1D *hcuts = (TH1D*)fHistDir->Get("hcuts");
  hcuts->GetXaxis()->SetBinLabel(200, fCutFile.c_str());
  string cstring = "B cand"; 
  int ibin;

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 
    
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

  } // end for 

}


// ----------------------------------------------------------------------
void candAnaFake::genMatch() {
  cout << "genMatch()  function" << endl;
}


// ----------------------------------------------------------------------
void candAnaFake::recoMatch() {
  cout << "recoMatch() function" << endl;
}


// ----------------------------------------------------------------------
void candAnaFake::candMatch() {
  cout << "candMatch()  function" << endl;
}


// ----------------------------------------------------------------------
void candAnaFake::efficiencyCalculation() {
  cout << "efficiencyCalculation() function" << endl;
}

// ----------------------------------------------------------------------
void candAnaFake::processType() {
  cout << "processType() function" << endl;
}
