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
#include "TFractionFitter.h"

#include "common/dataset.hh"
#include "common/util.hh"

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
void plotWork::makeAll(int bitmask, string what) {

  if (0 == bitmask) {
    makeAll(0x2);
    makeAll(0x4);
    makeAll(0x8);
    makeAll(0x10); 
    makeAll(0x11); 
    
  } else if (bitmask == 0x2) {
    zone(2, 2); 
    validation("RECO_5", "ptrelvsmuonpt", 4, 6, what, "bSignalHIL2Mu3", "cSignalHIL2Mu3");
    c0->cd(2);
    validation("RECO_5", "ptrelvsmuonpt", 6, 8, what, "bSignalHIL2Mu3", "cSignalHIL2Mu3");
    c0->cd(3);
    validation("RECO_5", "ptrelvsmuonpt", 8, 12, what, "bSignalHIL2Mu3", "cSignalHIL2Mu3");
    c0->cd(4);
    validation("RECO_5", "ptrelvsmuonpt",12, 30, what, "bSignalHIL2Mu3", "cSignalHIL2Mu3");
    c0->SaveAs(Form("ptrelvsmuonpt-%s.pdf", what.c_str()));
  } else if (bitmask == 0x4) {
    zone(2,2); 
    validation("RECO_5", "ptrelvsmuonpt", 8, 10, "dataMu8", "bSignalMu8", "cSignalMu8");
    c0->cd(2);
    validation("RECO_5", "ptrelvsmuonpt",10, 12, "dataMu8", "bSignalMu8", "cSignalMu8");
    c0->cd(3);
    validation("RECO_5", "ptrelvsmuonpt",12, 15, "dataMu8", "bSignalMu8", "cSignalMu8");
    c0->cd(4);
    validation("RECO_5", "ptrelvsmuonpt",15, 25, "dataMu8", "bSignalMu8", "cSignalMu8");
    c0->SaveAs("ptrelvsmuonpt-Mu8.pdf");
  } else if (bitmask == 0x8) {
    zone(2,2); 
    validation("RECO_5", "ptrelvsmuonpt", 24, 30, "dataMu24", "bSignalMu24", "cSignalMu24");
    c0->cd(2);
    validation("RECO_5", "ptrelvsmuonpt", 30, 40, "dataMu24", "bSignalMu24", "cSignalMu24");
    c0->cd(3);
    validation("RECO_5", "ptrelvsmuonpt", 40, 50, "dataMu24", "bSignalMu24", "cSignalMu24");
    c0->cd(4);
    validation("RECO_5", "ptrelvsmuonpt", 50, 100, "dataMu24", "bSignalMu24", "cSignalMu24");
    c0->SaveAs("ptrelvsmuonpt-Mu24.pdf");
  } else if (bitmask == 0x10) {
    zone(2,2); 
    validation("RECO_5", "ptrelvsmuoneta",-2, -1, what, "bSignalHIL2Mu3", "cSignalHIL2Mu3");
    c0->cd(2);
    validation("RECO_5", "ptrelvsmuoneta",-1,  0, what, "bSignalHIL2Mu3", "cSignalHIL2Mu3");
    c0->cd(3);
    validation("RECO_5", "ptrelvsmuoneta", 0,  1, what, "bSignalHIL2Mu3", "cSignalHIL2Mu3");
    c0->cd(4);
    validation("RECO_5", "ptrelvsmuoneta", 1,  2, what, "bSignalHIL2Mu3", "cSignalHIL2Mu3");
    c0->SaveAs(Form("ptrelvsmuonpt-%s.pdf", what.c_str()));
  } else if (bitmask == 0x11) {
    zone(2,2); 
    TH1D *h0 = getPtRel("RECO_5_1_ptrelvsmuonpt", "candAnaMuHIL2Mu3", "bSignalHIL2Mu3", 6, 10); 
    h0->Draw();

    c0->cd(2);
    TH1D *h1 = getPtRel("RECO_5_8_ptrelvsmuonpt", "candAnaMuHIL2Mu3", "bSignalHIL2Mu3", 6, 10); 
    h1->Draw();

    c0->cd(3);
    TH1D *h2 = getPtRel("RECO_5_9_ptrelvsmuonpt", "candAnaMuHIL2Mu3", "bSignalHIL2Mu3", 6, 10); 
    h2->Draw();

    c0->cd(4);
    TH1D *h3 = getPtRel("RECO_5_10_ptrelvsmuonpt", "candAnaMuHIL2Mu3", "bSignalHIL2Mu3", 6, 10); 
    h3->Draw();

    c0->SaveAs("ptrelvsmuonpt-procid.pdf");
    
    
  }

}



// ----------------------------------------------------------------------
void plotWork::hinValidation() {

  TH1D *hd(0), *hb(0), *hc(0);

  makeCanvas(4);
  //  c3->SetCanvasSize(1000, 400); 
  c3->SetWindowSize(1000, 400); 
  zone(4, 2, c3);
  
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(0); 

  int i(1); 
  
  double minPt(4.), maxPt(20.); 
  tl->SetNDC(kTRUE);
  tl->SetTextSize(0.07);
  numbers *a(0);
  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    c3->cd(i); 
    if (string::npos == it->first.find("run")) continue;
    a = new numbers;
    
    hd = getPtRel("RECO_5_0_ptrelvsmuonpt", "candAnaMuHIL2Mu3", it->first, minPt, maxPt); 
    cout << "hd = " << hd << ": " << hd->GetName() << endl;
    hb = getPtRel("RECO_5_1_ptrelvsmuonpt", "candAnaMuHIL2Mu3", "bSignalMuHIL2Mu3", minPt, maxPt); 
    cout << "hb = " << hb << ": " << hb->GetName() << endl;
    hc = getPtRel("RECO_5_2_ptrelvsmuonpt", "candAnaMuHIL2Mu3", "cSignalMuHIL2Mu3", minPt, maxPt); 
    cout << "hc = " << hc << ": " << hc->GetName() << endl;

    fitPtRel(a, hd, hb, hc);
    a->minPt = minPt;
    a->maxPt = maxPt;
    
    a->hD->SetTitle("");
    a->hD->SetMarkerStyle(24);
    a->hD->SetMarkerSize(0.7);
    a->hD->Draw();
    a->hSum->Draw("samehist");
    a->hB->SetLineWidth(1);
    setHist(a->hB, fDS["bSignalMuHIL2Mu3"]);
    a->hB->Draw("samehist");
    a->sname = it->first;
    tl->DrawLatex(0.2, 0.92, it->first.c_str()); 
    tl->DrawLatex(0.4, 0.7, Form("N_{B} = %4.1f #pm %4.1f", a->nB, a->nBE0)); 
    ++i;
    fVectorResult.push_back(a);     
  }  
  c3->SaveAs(Form("%s/hinValidation-fits.pdf", fDirectory.c_str()));

  
  TH1D *hruns = new TH1D("hruns", "", fVectorResult.size(), 0., fVectorResult.size());
  hruns->SetMinimum(0.); 
  setTitles(hruns, "runs in HINMuon", "N_{B}/nb^{-1}", 0.06, 1.1, 1.3, 0.06); 
  TH1D *hlumi = new TH1D("hlumi", "", 50, 0., 50.);
  hlumi->SetMinimum(0.); 
  setTitles(hlumi, "integrated lumi in HINMuon [nb^{-1}]", "N_{B}/nb^{-1}", 0.06, 1.1, 1.3, 0.06); 
  double ilumi(0);
  for (unsigned int i = 0; i < fVectorResult.size(); ++i) {
    hruns->GetXaxis()->SetBinLabel(i+1, fVectorResult[i]->sname.c_str());
    hruns->SetBinContent(i+1, 0.001*fVectorResult[i]->nB/fDS[fVectorResult[i]->sname]->fLumi); 
    hruns->SetBinError(i+1, 0.001*fVectorResult[i]->nBE0/fDS[fVectorResult[i]->sname]->fLumi); 

    ilumi += fDS[fVectorResult[i]->sname]->fLumi; 
    cout << "ilumi = " << ilumi*1.e3 << endl;
    hlumi->SetBinContent(hlumi->FindBin(ilumi*1.e3), 0.001*fVectorResult[i]->nB/fDS[fVectorResult[i]->sname]->fLumi); 
    hlumi->SetBinError(hlumi->FindBin(ilumi*1.e3), 0.001*fVectorResult[i]->nBE0/fDS[fVectorResult[i]->sname]->fLumi); 
  }

  makeCanvas(1);
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(0); 

  hruns->Fit("pol1");
  tl->DrawLatex(0.23, 0.4, Form("p0 = %4.1f #pm %4.1f",
				hruns->GetFunction("pol1")->GetParameter(0),
				hruns->GetFunction("pol1")->GetParError(0)));
  tl->DrawLatex(0.23, 0.32, Form("p1 = %3.2f #pm %3.2f",
				hruns->GetFunction("pol1")->GetParameter(1),
				 hruns->GetFunction("pol1")->GetParError(1)));
  c1->SaveAs(Form("%s/hinValidation-result-run.pdf", fDirectory.c_str()));

  hlumi->Fit("pol1");
  tl->DrawLatex(0.23, 0.4, Form("p0 = %4.1f #pm %4.1f",
				hlumi->GetFunction("pol1")->GetParameter(0),
				hlumi->GetFunction("pol1")->GetParError(0)));
  tl->DrawLatex(0.23, 0.32, Form("p1 = %3.2f #pm %3.2f",
				hlumi->GetFunction("pol1")->GetParameter(1),
				 hlumi->GetFunction("pol1")->GetParError(1)));
  c1->SaveAs(Form("%s/hinValidation-result-lumi.pdf", fDirectory.c_str()));

}


// ----------------------------------------------------------------------
void plotWork::dSigmadPt() {

  TH1D *hd(0), *hb(0), *hc(0);

  makeCanvas(4);
  c3->SetWindowSize(1000, 400); 
  zone(4, 3, c3);
  
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(0); 

  int i(1); 

  fVectorResult.clear();
  numbers *a(0);
  a = new numbers();  a->minPt =  4.;  a->maxPt =  5.;  a->sname = "dataMuHIL2Mu3"; a->dname = "candAnaMuHIL2Mu3";
  fVectorResult.push_back(a);
  a = new numbers();  a->minPt =  5.;  a->maxPt =  6.;  a->sname = "dataMuHIL2Mu3"; a->dname = "candAnaMuHIL2Mu3";
  fVectorResult.push_back(a);
  a = new numbers();  a->minPt =  6.;  a->maxPt =  7.;  a->sname = "dataMuHIL2Mu3"; a->dname = "candAnaMuHIL2Mu3";
  fVectorResult.push_back(a);
  a = new numbers();  a->minPt =  6.;  a->maxPt =  8.;  a->sname = "dataMuHIL2Mu3"; a->dname = "candAnaMuHIL2Mu3";
  fVectorResult.push_back(a);
  a = new numbers();  a->minPt =  8.;  a->maxPt = 10.;  a->sname = "dataMuHIL2Mu3"; a->dname = "candAnaMuHIL2Mu3";
  fVectorResult.push_back(a);
  a = new numbers();  a->minPt = 10.;  a->maxPt = 15.;  a->sname = "dataMuHIL2Mu3"; a->dname = "candAnaMuHIL2Mu3";
  fVectorResult.push_back(a);
  a = new numbers();  a->minPt = 15.;  a->maxPt = 20.;  a->sname = "dataMu8"; a->dname = "candAnaMu8";
  fVectorResult.push_back(a);
  a = new numbers();  a->minPt = 20.;  a->maxPt = 25.;  a->sname = "dataMu8"; a->dname = "candAnaMu8";
  fVectorResult.push_back(a);
  a = new numbers();  a->minPt = 25.;  a->maxPt = 30.;  a->sname = "dataMu24"; a->dname = "candAnaMu24";
  fVectorResult.push_back(a);
  a = new numbers();  a->minPt = 30.;  a->maxPt = 40.;  a->sname = "dataMu24"; a->dname = "candAnaMu24";
  fVectorResult.push_back(a);
  a = new numbers();  a->minPt = 40.;  a->maxPt = 50.;  a->sname = "dataMu24"; a->dname = "candAnaMu24";
  fVectorResult.push_back(a);
  a = new numbers();  a->minPt = 50.;  a->maxPt =100.;  a->sname = "dataMu50"; a->dname = "candAnaMu50";
  fVectorResult.push_back(a);

  tl->SetNDC(kTRUE);
  tl->SetTextSize(0.05);
  string dString, bString, cString, tString;
  for (unsigned int i = 0; i < fVectorResult.size(); ++i) {
    c3->cd(i+1); 
    a = fVectorResult[i];
    tString = a->sname;
    replaceAll(tString, "dataMu", ""); 

    cout << "**********************************************************************" << endl;
    cout << "***  " << i << ": " << a->sname << " " << a->minPt << " .. " << a->maxPt << " ***" << endl;
    cout << "RECO_5_0_ptrelvsmuonpt" << " .. " << a->dname << " .. " <<  a->sname<< endl;
    hd = getPtRel("RECO_5_0_ptrelvsmuonpt", a->dname, a->sname, a->minPt, a->maxPt); 
    if (0 == hd) {
      return;
    }
    cout << "hd = " << hd << ": " << hd->GetName() << endl;
    bString = "bSignalMu" + tString;
    cout << "RECO_5_1_ptrelvsmuonpt" << " .. " << a->dname << " .. " <<  bString << endl;
    hb = getPtRel("RECO_5_1_ptrelvsmuonpt", a->dname, bString, a->minPt, a->maxPt); 
    efficiency(a, bString); 
    if (0 == hb) {
      return;
    }
    cout << "hb = " << hb << ": " << hb->GetName() << endl;
    cString = "cSignalMu" + tString;
    cout << "RECO_5_2_ptrelvsmuonpt" << " .. " << a->dname << " .. " <<  cString << " "
	 <<  a->minPt << " " <<  a->maxPt << endl;
    hc = getPtRel("RECO_5_2_ptrelvsmuonpt", a->dname, cString, a->minPt, a->maxPt); 
    if (0 == hc) {
      return;
    }
    cout << "hc = " << hc << ": " << hc->GetName() << endl;

    fitPtRel(a, hd, hb, hc);
    
    a->hD->SetTitle("");
    a->hD->SetMarkerStyle(24);
    a->hD->SetMarkerSize(0.7);
    a->hD->Draw();
    tl->DrawLatex(0.2, 0.92, a->sname.c_str()); 
    cout << "pT: " <<  a->minPt << " " <<  a->maxPt << endl;
    tl->DrawLatex(0.4, 0.77, Form(" %3.0f < p_{T} < %3.0f", a->minPt, a->maxPt)); 
    if (0 == a->status) {
      a->hSum->Draw("samehist");
      a->hB->SetLineWidth(1);
      setHist(a->hB, fDS[bString]);
      a->hB->Draw("samehist");

      a->hC->SetLineWidth(1);
      setHist(a->hC, fDS[cString]);
      a->hC->Draw("samehist");
      tl->DrawLatex(0.4, 0.7, Form("N_{B} = %4.1f #pm %4.1f", a->nB, a->nBE0)); 
    }

  }  
  c3->SaveAs(Form("%s/dSigmadPt-fits.pdf", fDirectory.c_str()));

  
}




// ----------------------------------------------------------------------
// E.g. validation("RECO_5", "ptrelvsmuonpt", 4, 6, "run251721", "bSignalHIL2Mu3", "cSignalHIL2Mu3");
void plotWork::validation(string hist1, string hist2, double xmin, double xmax, string dname, string bname, string cname) {

  string dir("candAnaMu8");
  if (string::npos != dname.find("Mu24")) dir = "candAnaMu24";
  if (string::npos != dname.find("Mu50")) dir = "candAnaMu50";
  if (string::npos != dname.find("run")) dir = "candAnaMuHIL2Mu3";
  
  // -- my histograms
  TH2D *h2d = fDS[dname]->getHist2(Form("%s/%s_0_%s", dir.c_str(), hist1.c_str(), hist2.c_str()));
  TH2D *h2b = fDS[bname]->getHist2(Form("%s/%s_1_%s", dir.c_str(), hist1.c_str(), hist2.c_str()));
  TH2D *h2c = fDS[cname]->getHist2(Form("%s/%s_2_%s", dir.c_str(), hist1.c_str(), hist2.c_str()));

  int bin1 = h2d->GetXaxis()->FindBin(xmin); 
  int bin2 = h2d->GetXaxis()->FindBin(xmax); 

  cout << "xmin .. xmax: " << xmin << " .. " << xmax << " bin1 .. bin2: " << bin1 << " .. " << bin2 << endl;
  
  TH1D *hd = h2d->ProjectionY(Form("hd_%s_py_%d_%d", hist2.c_str(), bin1, bin2), bin1, bin2);   setHist(hd, fDS[dname]); 
  TH1D *hb = h2b->ProjectionY(Form("hb_%s_py_%d_%d", hist2.c_str(), bin1, bin2), bin1, bin2);   setHist(hb, fDS[bname]); 
  TH1D *hc = h2c->ProjectionY(Form("hc_%s_py_%d_%d", hist2.c_str(), bin1, bin2), bin1, bin2);   setHist(hc, fDS[cname]); 

  TObjArray *mc = new TObjArray(2);        // MC histograms are put in this array
  mc->Add(hb);
  mc->Add(hc);
  TFractionFitter* fit = new TFractionFitter(hd, mc); // initialise
  //  fit->Constrain(1, 0.0, 1.0);               // constrain fraction 1 to be between 0 and 1
  fit->SetRangeX(1, 50);                    // use only the first 15 bins in the fit
  Int_t status = fit->Fit();               // perform the fit
  cout << "fit status: " << status << endl;

  tl->SetNDC(kTRUE);
  if (status == 0) {                       // check on fit status
    double fracB, fracC, err;
    hd->DrawCopy("Ep");
    fit->GetResult(0, fracB, err);
    fit->GetResult(1, fracC, err);
    hb->Scale(fracB*hd->GetSumOfWeights()/hb->GetSumOfWeights());
    hc->Scale(fracC*hd->GetSumOfWeights()/hc->GetSumOfWeights());
    hb->DrawCopy("samehist");
    hc->DrawCopy("samehist");
    TH1D *result = (TH1D*)hb->Clone(Form("sum_%s_%d_%d", hist2.c_str(), bin1, bin2)); 
    result->SetLineColor(kBlack); 
    result->SetFillStyle(0);
    result->Add(hc); 
    result->Draw("samehist");
    if (string::npos != hist2.find("muonpt")) tl->DrawLatex(0.4, 0.85, Form("%3.1f < p_{T} < %3.1f GeV", xmin, xmax)); 
    if (string::npos != hist2.find("muoneta")) tl->DrawLatex(0.4, 0.85, Form("%3.1f < #eta < %3.1f", xmin, xmax)); 
    tl->DrawLatex(0.5, 0.78, dname.c_str()); 
    tl->DrawLatex(0.6, 0.71, Form("N_{B} ~ %4.0f", fracB*hd->GetSumOfWeights())); 
  }


}


// ----------------------------------------------------------------------
void plotWork::loopFunction1() {

}

// ----------------------------------------------------------------------
void plotWork::bookHist(int mode) {

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

}


// ----------------------------------------------------------------------
TH1D* plotWork::getPtRel(string histname, string dir, string dname, double xmin, double xmax) {
  cout << "trying to get " << Form("%s/%s", dir.c_str(), histname.c_str()) << endl;
  TH2D *h2 = fDS[dname]->getHist2(Form("%s/%s", dir.c_str(), histname.c_str()));
  if (0 == h2) return 0; 
  int bin1 = h2->GetXaxis()->FindBin(xmin); 
  int bin2 = h2->GetXaxis()->FindBin(xmax); 

  cout << "projecting into " << Form("proj_%s_%d_%d_%s_%s", h2->GetName(), bin1, bin2, dir.c_str(), dname.c_str()) << endl;
  
  TH1D *h = h2->ProjectionY(Form("proj_%s_%d_%d_%s_%s", h2->GetName(), bin1, bin2, dir.c_str(), dname.c_str()), bin1, bin2);
  return h;
  
}

// ----------------------------------------------------------------------
void plotWork::efficiency(numbers *a, string bString) {
  // FIXME
  // TH1D *hRec = getPtRel("RECO_5_1_muon_pt", a->dname, bString, a->minPt, a->maxPt); 
  // TH1D *hGen = getPtRel("GEN_5_1_muon_pt", a->dname, bString, a->minPt, a->maxPt); 

  cout << "*** Efficiency: rec = " << hRec->GetSumOfWeights() << endl;
  cout << "*** Efficiency: gen = " << hGen->GetSumOfWeights() << endl;
   

}



// ----------------------------------------------------------------------
void plotWork::fitPtRel(numbers* n, TH1D* hd, TH1D* hb, TH1D* hc, TH1D* hl) {

  n->status = -1; 

  TObjArray *mc(0);
  if (0 == hl) {
    mc = new TObjArray(2);
    mc->Add(hb);
    mc->Add(hc);
  } else {
    mc = new TObjArray(3);
    mc->Add(hb);
    mc->Add(hc);
    mc->Add(hl);
  }

  TFractionFitter* fit = new TFractionFitter(hd, mc); 
  fit->SetRangeX(1, 50); 
  Int_t status = fit->Fit();
  cout << "fitted " << hd->GetName() << ", status: " << status << endl;

  n->nData = hd->GetSumOfWeights(); 
  n->hD = hd; 
  if (status == 0) {                       // check on fit status
    n->status  = 0; 
    n->nData   = hd->GetSumOfWeights();
    n->nData2  = hd->Integral();
    
    double fracB, fracC, fracL, err;
    fit->GetResult(0, fracB, err);
    n->fracB   = fracB; 
    n->fracBE0 = err; 
    n->nB      = fracB*hd->GetSumOfWeights();
    n->nBE0    = (err/fracB)*n->nB;
    n->nB2     = fracB*hd->Integral();
    
    fit->GetResult(1, fracC, err);
    n->fracC  = fracC; 
    n->fracCE0 = err; 
    n->nC      = fracC*hd->GetSumOfWeights();
    n->nCE0    = (err/fracC)*n->nC;
    n->nC2     = fracC*hd->Integral();

    if (0 != hl) {
      fit->GetResult(1, fracL, err);
      n->fracL  = fracL; 
      n->fracLE0 = err; 
      n->nL      = fracL*hd->GetSumOfWeights();
      n->nLE0    = (err/fracL)*n->nL;
      n->nL2     = fracL*hd->Integral();
    }
    
    hb->Scale(fracB*hd->GetSumOfWeights()/hb->GetSumOfWeights());
    hc->Scale(fracC*hd->GetSumOfWeights()/hc->GetSumOfWeights());
    if (hl) hl->Scale(fracL*hd->GetSumOfWeights()/hl->GetSumOfWeights());
    n->hB = hb; 
    n->hC = hc; 
    n->hL = hl; 
    
    n->hSum = (TH1D*)hb->Clone(Form("sum_%s", hd->GetName())); 
    n->hSum->SetLineColor(kBlack); 
    n->hSum->SetFillStyle(0);
    n->hSum->Add(hc); 
    if (hl) n->hSum->Add(hl); 
    n->hPlot = (TH1D*)fit->GetPlot();
  }


  //FIXME??  delete mc;
  delete fit;
  
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


// ----------------------------------------------------------------------
void plotWork::loadFiles(string afiles) {
  
  string files = fDirectory + string("/") + afiles;
  cout << "==> Loading files listed in " << files << endl;

  char buffer[1000];
  ifstream is(files.c_str());
  while (is.getline(buffer, 1000, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    
    string sbuffer = string(buffer); 
    replaceAll(sbuffer, " ", ""); 
    replaceAll(sbuffer, "\t", ""); 
    if (sbuffer.size() < 1) continue;
    
    string::size_type m1 = sbuffer.find("lumi="); 
    string stype = sbuffer.substr(5, m1-5); 
    
    string::size_type m2 = sbuffer.find("file="); 
    string slumi = sbuffer.substr(m1+5, m2-m1-6); 
    string sfile = sbuffer.substr(m2+5); 
    string sname, sdecay; 
    
    cout << "stype: ->" << stype << "<-" << endl;
    
    TFile *pF(0); 
    if (string::npos != stype.find("data")) {
      // -- DATA
      pF = loadFile(sfile); 
      
      dataset *ds = new dataset(); 
      ds->fSize = 1; 
      ds->fWidth = 2;
      ds->fLumi = atof(slumi.c_str()); 
      
      if (string::npos != stype.find("HIL2Mu3")) {
        sname = "dataMuHIL2Mu3"; 
        sdecay = "HIL2Mu3"; 
	ds->fColor = kBlack; 
	ds->fSymbol = 20; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3350; 
      }

      if (string::npos != stype.find("run251721")) {
	sname = "run251721"; 
        sdecay = "HIL2Mu3 (251721)"; 
	ds->fColor = kBlack; 
	ds->fSymbol = 20; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3350; 
      }

      if (string::npos != stype.find("run254987")) {
        sname = "run254987"; 
        sdecay = "HIL2Mu3 (254987)"; 
	ds->fColor = kBlack; 
	ds->fSymbol = 20; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3350; 
      }

      if (string::npos != stype.find("run254989")) {
        sname = "run254989"; 
        sdecay = "HIL2Mu3 (254989)"; 
	ds->fColor = kBlack; 
	ds->fSymbol = 20; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3350; 
      }

      if (string::npos != stype.find("run254993")) {
        sname = "run254993"; 
        sdecay = "HIL2Mu3 (254993)"; 
	ds->fColor = kBlack; 
	ds->fSymbol = 20; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3350; 
      }

      if (string::npos != stype.find("run255019")) {
        sname = "run255019"; 
        sdecay = "HIL2Mu3 (255019)"; 
	ds->fColor = kBlack; 
	ds->fSymbol = 20; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3350; 
      }

      if (string::npos != stype.find("run255029")) {
        sname = "run255029"; 
        sdecay = "HIL2Mu3 (255029)"; 
	ds->fColor = kBlack; 
	ds->fSymbol = 20; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3350; 
      }

      if (string::npos != stype.find("run255031")) {
        sname = "run255031"; 
        sdecay = "HIL2Mu3 (255031)"; 
	ds->fColor = kBlack; 
	ds->fSymbol = 20; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3350; 
      }

      if (string::npos != stype.find("Mu8")) {
        sname = "dataMu8"; 
        sdecay = "Mu8"; 
	ds->fColor = kBlack; 
	ds->fSymbol = 20; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3350; 
      }

      if (string::npos != stype.find("Mu24")) {
        sname = "dataMu24"; 
        sdecay = "Mu24"; 
	ds->fColor = kBlack; 
	ds->fSymbol = 20; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3350; 
      }

      if (string::npos != stype.find("Mu50")) {
        sname = "dataMu50"; 
        sdecay = "Mu50"; 
	ds->fColor = kBlack; 
	ds->fSymbol = 20; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3350; 
      }
      
      ds->fLcolor = ds->fColor; 
      ds->fFcolor = ds->fColor; 
      ds->fName   = sdecay; 
      ds->fFullName = sname; 
      fDS.insert(make_pair(sname, ds)); 
      
      
    } else {
      // -- MC
      pF = loadFile(sfile); 
      
      dataset *ds = new dataset(); 
      ds->fSize = 1; 
      ds->fWidth = 2; 
 
      if (string::npos != stype.find("bsignal,HIL2Mu3")) {
        sname = "bSignalMuHIL2Mu3"; 
        sdecay = "bSignalMuHIL2Mu3"; 
	ds->fColor = kBlue-7; 
	ds->fSymbol = 24; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365; 
      }

      if (string::npos != stype.find("csignal,HIL2Mu3")) {
        sname = "cSignalMuHIL2Mu3"; 
        sdecay = "cSignalMuHIL2Mu3"; 
	ds->fColor = kGreen+3; 
	ds->fSymbol = 24; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3356; 
      }


      if (string::npos != stype.find("bsignal,Mu8")) {
        sname = "bSignalMu8"; 
        sdecay = "bSignalMu8"; 
	ds->fColor = kBlue-7; 
	ds->fSymbol = 24; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365; 
      }

      if (string::npos != stype.find("csignal,Mu8")) {
        sname = "cSignalMu8"; 
        sdecay = "cSignalMu8"; 
	ds->fColor = kGreen+3; 
	ds->fSymbol = 24; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3356; 
      }

      if (string::npos != stype.find("bsignal,Mu24")) {
        sname = "bSignalMu24"; 
        sdecay = "bSignalMu24"; 
	ds->fColor = kBlue-7; 
	ds->fSymbol = 24; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365; 
      }

      if (string::npos != stype.find("csignal,Mu24")) {
        sname = "cSignalMu24"; 
        sdecay = "cSignalMu24"; 
	ds->fColor = kGreen+3; 
	ds->fSymbol = 24; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3356; 
      }

      if (string::npos != stype.find("bsignal,Mu50")) {
        sname = "bSignalMu50"; 
        sdecay = "bSignalMu50"; 
	ds->fColor = kBlue-7; 
	ds->fSymbol = 24; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365; 
      }

      if (string::npos != stype.find("csignal,Mu50")) {
        sname = "cSignalMu50"; 
        sdecay = "cSignalMu50"; 
	ds->fColor = kGreen+3; 
	ds->fSymbol = 24; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3356; 
      }


      ds->fLcolor = ds->fColor; 
      ds->fFcolor = ds->fColor; 
      ds->fName   = sdecay; 
      ds->fFullName = sname; 
      fDS.insert(make_pair(sname, ds)); 
      
    }
  }

  is.close();
  cout << "Summary: " << endl;
  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    cout << "===>" << it->first << "<==" << endl;
    cout << "       " << it->second->fName << endl;
    //    cout << "       " << it->second->fF->GetName() << endl;
    cout << "       " << it->first << ": " << it->second->fName << ", " << it->second->fF->GetName() << endl;
  }
}

  
