#include "plotWork.hh"

#include <fstream>
#include <iostream>
#include <sstream>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TFitResult.h"

#include "dataset.hh"
#include "util.hh"

ClassImp(plotWork)

using namespace std; 

// ----------------------------------------------------------------------
plotWork::plotWork(string dir,  string files, string setup): plotClass(dir, files, setup) {
  loadFiles(files);

  if (setup == "") {
    fHistFileName = Form("%s/plotWork.root", dir.c_str()); 
  } else {
    fHistFileName = Form("%s/plotWork-%s.root", dir.c_str(), setup.c_str()); 
  }

  fTexFileName = fHistFileName; 
  replaceAll(fTexFileName, ".root", ".tex"); 
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));

}


// ----------------------------------------------------------------------
plotWork::~plotWork() {

}


// ----------------------------------------------------------------------
void plotWork::makeAll(int bitmask) {

}



// ----------------------------------------------------------------------
void plotWork::triggerEff(string ds, string dir) {

  //  TFile *f = TFile::Open(fHistFileName.c_str(), "UPDATE");
  //  cout << "file: " << f << endl;
 
  fCds = ds; 
  fSetup = dir; 
  fHists.clear();

  char hist[200], thist[200], ahist[200];

  // -- book histograms
  sprintf(hist, "triggerEff_%s_%s", fCds.c_str(), fSetup.c_str()); 
  sprintf(thist, "%s", fDS[fCds]->fName.c_str()); 
  sprintf(ahist, "p_{T} [GeV]"); 
  if (fHists.count(hist) > 0) {
    fHists[hist]->Reset();
  } else {
    fHists.insert(make_pair(hist, new TH1D(hist, thist, 50, 0, 25.))); 
    fHists[hist]->Sumw2();
    setHistTitles(fHists[hist], fDS[fCds], ahist, "Entries/bin");
  }
  
  sprintf(hist, "ntriggerEff_%s_%s", fCds.c_str(), fSetup.c_str()); 
  if (fHists.count(hist) > 0) {
    fHists[hist]->Reset();
  } else {
    fHists.insert(make_pair(hist, new TH1D(hist, thist, 50, 0, 25.))); 
    fHists[hist]->Sumw2();
    setHistTitles(fHists[hist], fDS[fCds], ahist, "Efficeincy");
  }
  
  sprintf(hist, "NtriggerEff_%s_%s", fCds.c_str(), fSetup.c_str()); 
  if (fHists.count(hist) > 0) {
    fHists[hist]->Reset();
  } else {
    fHists.insert(make_pair(hist, new TH1D(hist, thist, 50, 0, 25.))); 
    fHists[hist]->Sumw2();
  }
  


  TTree *t = getTree(ds, dir); 
  cout << "tree: " << t << endl;
  setupTree(t); 
  loopOverTree(t, 1); 

  string nname = Form("ntriggerEff_%s_%s", fCds.c_str(), fSetup.c_str());
  fHists[nname]->Draw();
  string Nname = Form("NtriggerEff_%s_%s", fCds.c_str(), fSetup.c_str());
  fHists[Nname]->Draw("same");
  
  string ename = Form("triggerEff_%s_%s", fCds.c_str(), fSetup.c_str());
  fHists[ename]->Divide(fHists[nname], fHists[Nname], 1., 1., "B");
  fHists[ename]->Draw();

//   f->Write();
//   f->Close();
}


// ----------------------------------------------------------------------
void plotWork::loopFunction1() {
  if (fb.m1pt < 3.5) return;
  if (TMath::Abs(fb.m1eta) > 2.1) return;
      
  if (fb.m2pt < 3.5) return;
  if (TMath::Abs(fb.m2eta) > 2.1) return;

  if (!fb.gmuid) return;

  string hname = Form("NtriggerEff_%s_%s", fCds.c_str(), fSetup.c_str());
  fHists[hname]->Fill(fb.pt);
  hname = Form("ntriggerEff_%s_%s", fCds.c_str(), fSetup.c_str());
  if (fb.hlt) {
    fHists[hname]->Fill(fb.pt);
  }

}

// ----------------------------------------------------------------------
void plotWork::bookHist(int mode) {
  char hist[200], thist[200], ahist[200];

  // -- triggerEff and loopfunction1
  if (1 == mode) {
    sprintf(hist, "triggerEff_%s_%s", fCds.c_str(), fSetup.c_str()); 
    sprintf(thist, "triggerEff %s %s", fCds.c_str(), fSetup.c_str()); 
    sprintf(ahist, "pT [GeV]"); 
    if (fHists.count(hist) > 0) {
      fHists[hist]->Reset();
    } else {
      fHists.insert(make_pair(hist, new TH1D(hist, thist, 50, 0, 25.))); 
      setHistTitles(fHists[hist], fDS[fCds], ahist, "Entries/bin");
    }

    sprintf(hist, "ntriggerEff_%s_%s", fCds.c_str(), fSetup.c_str()); 
    sprintf(thist, "ntriggerEff %s %s", fCds.c_str(), fSetup.c_str()); 
    if (fHists.count(hist) > 0) {
      fHists[hist]->Reset();
    } else {
      fHists.insert(make_pair(hist, new TH1D(hist, thist, 50, 0, 25.))); 
      setHistTitles(fHists[hist], fDS[fCds], ahist, "Entries/bin");
    }

    sprintf(hist, "NtriggerEff_%s_%s", fCds.c_str(), fSetup.c_str()); 
    sprintf(thist, "NtriggerEff %s %s", fCds.c_str(), fSetup.c_str()); 
    if (fHists.count(hist) > 0) {
      fHists[hist]->Reset();
    } else {
      fHists.insert(make_pair(hist, new TH1D(hist, thist, 50, 0, 25.))); 
      setHistTitles(fHists[hist], fDS[fCds], ahist, "Entries/bin");
    }

  }    

  cout << "booked all hists" << endl;
}



// ----------------------------------------------------------------------
void plotWork::candAnalysis() {
  fGoodCand = true; 
}


// ----------------------------------------------------------------------
void plotWork::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
  int nentries = Int_t(t->GetEntries());
  int nbegin(0), nend(nentries); 
  if (nevts > 0 && nentries > nevts) {
    nentries = nevts;
    nbegin = 0; 
    nend = nevts;
  }
  if (nevts > 0 && nstart > 0) {
    nentries = nstart + nevts;
    nbegin = nstart; 
    if (nstart + nevts < t->GetEntries()) {
      nend = nstart + nevts; 
    } else {
      nend = t->GetEntries();
    }
  }
  
  nentries = nend - nstart; 
  
  int step(1000000); 
  if (nentries < 5000000)  step = 500000; 
  if (nentries < 1000000)  step = 100000; 
  if (nentries < 100000)   step = 10000; 
  if (nentries < 10000)    step = 1000; 
  if (nentries < 1000)     step = 100; 
  step = 500000; 
  cout << "==> plotWork::loopOverTree> loop over dataset " << fCds << " in file " 
       << t->GetDirectory()->GetName() 
       << " with " << nentries << " entries" 
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotWork::*pF)(void);
  if (ifunc == 1) pF = &plotWork::loopFunction1;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;
   
    candAnalysis();
    (this->*pF)();
  }

}


// ----------------------------------------------------------------------
void plotWork::setupTree(TTree *t) {
  if (string::npos != fCds.find("Mc")) {
    fIsMC = true; 
  } else {
    fIsMC = false; 
  }

  t->SetBranchAddress("pt", &fb.pt);
  t->SetBranchAddress("q", &fb.q);

  t->SetBranchAddress("tau", &fb.tau);
  t->SetBranchAddress("gtau", &fb.gtau);

  t->SetBranchAddress("bdt",&fb.bdt);
  t->SetBranchAddress("bdt",&fb.bdt);
  t->SetBranchAddress("lip",&fb.lip);
  t->SetBranchAddress("lipE",&fb.lipE);
  t->SetBranchAddress("tip",&fb.tip);
  t->SetBranchAddress("tipE",&fb.tipE);

  t->SetBranchAddress("closetrk",&fb.closetrk);
  t->SetBranchAddress("pvlip",   &fb.pvlip);
  t->SetBranchAddress("pvlips",  &fb.pvlips);
  t->SetBranchAddress("pvlip2",  &fb.pvlip2);
  t->SetBranchAddress("pvlips2", &fb.pvlips2);
  t->SetBranchAddress("maxdoca", &fb.maxdoca);
  t->SetBranchAddress("pvip",    &fb.pvip);
  t->SetBranchAddress("pvips",   &fb.pvips);
  t->SetBranchAddress("pvip3d",  &fb.pvip3d);
  t->SetBranchAddress("pvips3d", &fb.pvips3d);
  t->SetBranchAddress("pvw8",    &fb.pvw8);

  t->SetBranchAddress("m1pix",    &fb.m1pix);
  t->SetBranchAddress("m2pix",    &fb.m2pix);
  t->SetBranchAddress("m1bpix",   &fb.m1bpix);
  t->SetBranchAddress("m2bpix",   &fb.m2bpix);
  t->SetBranchAddress("m1bpixl1", &fb.m1bpixl1);
  t->SetBranchAddress("m2bpixl1", &fb.m2bpixl1);

  t->SetBranchAddress("rr",     &fb.rr);
  t->SetBranchAddress("pvn",    &fb.pvn);
  t->SetBranchAddress("run",    &fb.run);
  t->SetBranchAddress("evt",    &fb.evt);
  t->SetBranchAddress("hlt",    &fb.hlt);
  t->SetBranchAddress("hltm",   &fb.hltm);
  t->SetBranchAddress("hltm2",  &fb.hltm2);
  t->SetBranchAddress("ls",     &fb.ls);
  t->SetBranchAddress("cb",     &fb.cb);
  t->SetBranchAddress("json",   &fb.json);
  t->SetBranchAddress("gmuid",  &fb.gmuid);
  t->SetBranchAddress("gmutmid", &fb.gmutmid);
  t->SetBranchAddress("gmumvaid", &fb.gmumvaid);
  t->SetBranchAddress("gtqual", &fb.gtqual);
  t->SetBranchAddress("tm",     &fb.tm);
  t->SetBranchAddress("procid", &fb.procid);
  t->SetBranchAddress("m",      &fb.m);
  t->SetBranchAddress("m3",     &fb.m3);
  t->SetBranchAddress("m4",     &fb.m4);
  t->SetBranchAddress("me",     &fb.me);
  t->SetBranchAddress("cm",     &fb.cm);
  t->SetBranchAddress("pt",     &fb.pt);
  t->SetBranchAddress("phi",    &fb.phi);
  t->SetBranchAddress("eta",    &fb.eta);
  t->SetBranchAddress("cosa",   &fb.cosa);
  t->SetBranchAddress("alpha",  &fb.alpha);
  t->SetBranchAddress("iso",    &fb.iso);
  t->SetBranchAddress("chi2",   &fb.chi2);
  t->SetBranchAddress("dof",    &fb.dof);
  t->SetBranchAddress("prob",   &fb.pchi2dof);
  t->SetBranchAddress("chi2dof",&fb.chi2dof);
  t->SetBranchAddress("flsxy",  &fb.flsxy);
  t->SetBranchAddress("fls3d",  &fb.fls3d);
  t->SetBranchAddress("fl3d",   &fb.fl3d);
  t->SetBranchAddress("fl3dE",  &fb.fl3dE);
  t->SetBranchAddress("m1pt",   &fb.m1pt);
  t->SetBranchAddress("m1gt",   &fb.m1gt);
  t->SetBranchAddress("m1eta",  &fb.m1eta);
  t->SetBranchAddress("m1phi",  &fb.m1phi);
  t->SetBranchAddress("m1q",    &fb.m1q);
  t->SetBranchAddress("m2pt",   &fb.m2pt);
  t->SetBranchAddress("m2gt",   &fb.m2gt);
  t->SetBranchAddress("m2eta",  &fb.m2eta);
  t->SetBranchAddress("m2phi",  &fb.m2phi);
  t->SetBranchAddress("m2q",    &fb.m2q);
  t->SetBranchAddress("docatrk",&fb.docatrk);

  t->SetBranchAddress("m1id",     &fb.m1id);
  t->SetBranchAddress("m1rmvaid", &fb.m1rmvaid);
  t->SetBranchAddress("m1trigm",  &fb.m1trigm);
  t->SetBranchAddress("m1rmvabdt",&fb.m1rmvabdt);
  t->SetBranchAddress("m1tmid",   &fb.m1tmid);
 
  t->SetBranchAddress("m2id",     &fb.m2id);
  t->SetBranchAddress("m2rmvaid", &fb.m2rmvaid);
  t->SetBranchAddress("m2trigm",  &fb.m2trigm);
  t->SetBranchAddress("m2rmvabdt",&fb.m2rmvabdt);
  t->SetBranchAddress("m2tmid",   &fb.m2tmid);


  t->SetBranchAddress("m1iso",     &fb.m1iso);
  t->SetBranchAddress("m2iso",     &fb.m2iso);
  t->SetBranchAddress("closetrks1",&fb.closetrks1);
  t->SetBranchAddress("closetrks2",&fb.closetrks2);
  t->SetBranchAddress("closetrks3",&fb.closetrks3);
  t->SetBranchAddress("othervtx",  &fb.othervtx);
  t->SetBranchAddress("pvdchi2",   &fb.pvdchi2);

  t->SetBranchAddress("g1pt",   &fb.g1pt);
  t->SetBranchAddress("g2pt",   &fb.g2pt);
  t->SetBranchAddress("g1eta",  &fb.g1eta);
  t->SetBranchAddress("g2eta",  &fb.g2eta);
  t->SetBranchAddress("g1id",   &fb.g1id);
  t->SetBranchAddress("g2id",   &fb.g2id);
  if (string::npos != fCds.find("No")) {
    if (string::npos != fCds.find("Mc")) {
      t->SetBranchAddress("g3pt", &fb.g3pt);
      t->SetBranchAddress("g3eta",&fb.g3eta);
    }
    t->SetBranchAddress("kpt",  &fb.k1pt);
    t->SetBranchAddress("kgt",  &fb.k1gt);
    t->SetBranchAddress("keta", &fb.k1eta);
    t->SetBranchAddress("mpsi", &fb.mpsi);
    t->SetBranchAddress("psipt",&fb.psipt); //FIXME
  }

  if (string::npos != fCds.find("Cs")) {
    if (string::npos != fCds.find("Mc")) {
      t->SetBranchAddress("g3pt", &fb.g3pt);
      t->SetBranchAddress("g3eta",&fb.g3eta);
      t->SetBranchAddress("g4pt", &fb.g4pt);
      t->SetBranchAddress("g4eta",&fb.g4eta);
    }
    t->SetBranchAddress("psipt",&fb.psipt);   //FIXME
    t->SetBranchAddress("mpsi", &fb.mpsi);
    t->SetBranchAddress("mkk",  &fb.mkk);
    t->SetBranchAddress("dr",   &fb.dr);
    t->SetBranchAddress("k1pt", &fb.k1pt);
    t->SetBranchAddress("k1gt", &fb.k1gt);
    t->SetBranchAddress("k1eta",&fb.k1eta);
    t->SetBranchAddress("k2pt", &fb.k2pt);
    t->SetBranchAddress("k2gt", &fb.k2gt);
    t->SetBranchAddress("k2eta",&fb.k2eta);
  } else {
    fb.mkk = 999.;
    fb.dr = 999.;
  }

  if (string::npos != fCds.find("DstarPi")) {
    t->SetBranchAddress("md0",&fb.md0);
    t->SetBranchAddress("dm",&fb.dm);
    t->SetBranchAddress("ptd0",&fb.ptd0);
  }

}


// ----------------------------------------------------------------------
void plotWork::setCuts(string cuts) {
  cout << "==> plotWork::setCuts: " << cuts << endl;

  istringstream ss(cuts);
  string token, name, sval;

  while (getline(ss, token, ',')) {
    
    string::size_type m1 = token.find("="); 
    name = token.substr(0, m1);
    sval = token.substr(m1+1);

    if (string::npos != name.find("PTLO")) {
      float val; 
      val = atof(sval.c_str()); 
      PTLO = val;
    }

  }
}
