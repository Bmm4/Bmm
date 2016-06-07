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
  if (fCandDoca > 0.01) skip = true;
  if (!fGoodTracks) skip = true;
  if (!fGoodQ) skip = true;
  // FIXME: These should be replaced with the proper signed 3D impact parameter (to be filled!)
  //  if (fpMuon1->fLip/fpMuon1->fLipE < 1.) skip = true;
  //  if (fpMuon2->fLip/fpMuon2->fLipE < 1.) skip = true;

  if (skip) return;

  if (310 == TRUTHCAND) {
    // see bmm4.org::160217a
    //    if (fCandA > 0.01)        skip = true;
    if (fCandFLxy > 4.0)      skip = true;
    if (fCandFLSxy < 10.0)    skip = true;
    if (fCandPvIpS > 5.0)     skip = true;

    // -- calculate Lambda mass
    TLorentzVector pr, pi, la;
    pr.SetPtEtaPhiM(fpMuon1->fPlab.Perp(), fpMuon1->fPlab.Eta(), fpMuon1->fPlab.Phi(), MPROTON);
    pi.SetPtEtaPhiM(fpMuon2->fPlab.Perp(), fpMuon2->fPlab.Eta(), fpMuon2->fPlab.Phi(), MPION);
    la = pr + pi;
    double mla = la.M();
    if (TMath::Abs(mla - MLAMBDA_0) < 0.006) skip = true;
  } else if (333 == TRUTHCAND) {
    // see bmm4.org::160217a
    if (fCandFLxy > 1.0)      skip = true;
    if (fCandFLSxy < 1.0)    skip = true;
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

  bool muid(false);
  // -- leg 1
  if (fMu1Pt >  3.5) {
    muid = ((fpMuon1->fMuID & 2) == 2); // global muon!
    ((TH1D*)fHistDir->Get("pt1all"))->Fill(fMu1Pt, fpCand->fMass);
    ((TH1D*)fHistDir->Get("eta1all"))->Fill(TMath::Abs(fMu1Eta), fpCand->fMass);
    if (muid) {
      ((TH1D*)fHistDir->Get("pt1muo"))->Fill(fMu1Pt, fpCand->fMass);
      ((TH1D*)fHistDir->Get("eta1muo"))->Fill(TMath::Abs(fMu1Eta), fpCand->fMass);
    } else {
      ((TH1D*)fHistDir->Get("pt1nmu"))->Fill(fMu1Pt, fpCand->fMass);
      ((TH1D*)fHistDir->Get("eta1nmu"))->Fill(TMath::Abs(fMu1Eta), fpCand->fMass);
    }
  }

  // -- leg 2
  if (fMu2Pt >  3.5) {
    muid = ((fpMuon2->fMuID & 2) == 2); // global muon!
    ((TH1D*)fHistDir->Get("pt2all"))->Fill(fMu2Pt, fpCand->fMass);
    ((TH1D*)fHistDir->Get("eta2all"))->Fill(TMath::Abs(fMu2Eta), fpCand->fMass);
    if (muid) {
      ((TH1D*)fHistDir->Get("pt2muo"))->Fill(fMu2Pt, fpCand->fMass);
      ((TH1D*)fHistDir->Get("eta2muo"))->Fill(TMath::Abs(fMu2Eta), fpCand->fMass);
    } else {
      ((TH1D*)fHistDir->Get("pt2nmu"))->Fill(fMu2Pt, fpCand->fMass);
      ((TH1D*)fHistDir->Get("eta2nmu"))->Fill(TMath::Abs(fMu2Eta), fpCand->fMass);
    }
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
  double mmin(0.45), mmax(0.55);
  double mbin(50);
  TH1::SetDefaultSumw2(kTRUE);
  vector<string> mode;
  mode.push_back("all");
  mode.push_back("nmu");
  mode.push_back("muo");
  if (310 == TRUTHCAND) {
    new TH1D("mall", Form("%s", name.c_str()), mbin, mmin, mmax);
    for (unsigned i = 0; i < mode.size(); ++i) {
      new TH2D(Form("pt1%s", mode[i].c_str()), Form("pt1 pion (%s)", mode[i].c_str()), nbins, xbins, mbin, mmin, mmax);
      new TH2D(Form("pt2%s", mode[i].c_str()), Form("pt2 pion (%s)", mode[i].c_str()), nbins, xbins, mbin, mmin, mmax);
      new TH2D(Form("eta1%s", mode[i].c_str()), Form("eta1 pion (%s)", mode[i].c_str()), 6, 0., 2.4, mbin, mmin, mmax);
      new TH2D(Form("eta2%s", mode[i].c_str()), Form("eta2 pion (%s)", mode[i].c_str()), 6, 0., 2.4, mbin, mmin, mmax);
    }
  } else if (333 == TRUTHCAND) {
    name = "#phi";
    mmin = 0.98;
    mmax = 1.06;
    new TH1D("mall", Form("%s", name.c_str()), mbin, mmin, mmax);
    for (unsigned i = 0; i < mode.size(); ++i) {
      new TH2D(Form("pt1%s", mode[i].c_str()), Form("pt1 pion (%s)", mode[i].c_str()), nbins, xbins, mbin, mmin, mmax);
      new TH2D(Form("pt2%s", mode[i].c_str()), Form("pt2 pion (%s)", mode[i].c_str()), nbins, xbins, mbin, mmin, mmax);
      new TH2D(Form("eta1%s", mode[i].c_str()), Form("eta1 pion (%s)", mode[i].c_str()), 6, 0., 2.4, mbin, mmin, mmax);
      new TH2D(Form("eta2%s", mode[i].c_str()), Form("eta2 pion (%s)", mode[i].c_str()), 6, 0., 2.4, mbin, mmin, mmax);
    }
  } else if (3122 == TRUTHCAND) {
    name = "#Lambda";
    mmin = 1.05;
    mmax = 1.20;
    new TH1D("mall", Form("%s", name.c_str()), mbin, mmin, mmax);
    for (unsigned i = 0; i < mode.size(); ++i) {
      new TH2D(Form("pt1%s", mode[i].c_str()), Form("pt1 pion (%s)", mode[i].c_str()), nbins, xbins, mbin, mmin, mmax);
      new TH2D(Form("pt2%s", mode[i].c_str()), Form("pt2 pion (%s)", mode[i].c_str()), nbins, xbins, mbin, mmin, mmax);
      new TH2D(Form("eta1%s", mode[i].c_str()), Form("eta1 pion (%s)", mode[i].c_str()), 6, 0., 2.4, mbin, mmin, mmax);
      new TH2D(Form("eta2%s", mode[i].c_str()), Form("eta2 pion (%s)", mode[i].c_str()), 6, 0., 2.4, mbin, mmin, mmax);
    }
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
void candAnaFake::genMatchOld() {
  cout << "genMatchOld()  function" << endl;
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
