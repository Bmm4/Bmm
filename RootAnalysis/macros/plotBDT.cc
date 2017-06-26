#include "plotBDT.hh"

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>


#include "TROOT.h"
#include "TStyle.h"
#include "TKey.h"
#include "TMath.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TPad.h"
#include "TF1.h"
#include "THStack.h"
#include "TFitResult.h"

#include "common/dataset.hh"
#include "common/util.hh"
#include "common/Lumi.hh"
#include "common/fitPsYield.hh"

#include "preselection.hh"

ClassImp(plotBDT)

using namespace std;

// ----------------------------------------------------------------------
plotBDT::plotBDT(string dir, string files, string cuts, string setup) : plotClass(dir, files, cuts, setup) {
  cout << "plotBDT: plotClass::loadFiles(" << files << ")" << endl;
  plotClass::loadFiles(files);
  changeSetup(dir, "plotBDT", setup);
  init();

  fTEX.close();
  cout << Form("/bin/rm -f %s", fTexFileName.c_str()) << endl;
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  cout << Form("open for TeX output: %s", fTexFileName.c_str()) << endl;
  fTEX.open(fTexFileName.c_str(), ios::app);

  string hfname  = fDirectory + "/plotBDT." + fSuffix + ".root";

}



// ----------------------------------------------------------------------
plotBDT::~plotBDT() {
  if (fHistFile) {
    fHistFile->Close();
  }
}


// ----------------------------------------------------------------------
void plotBDT::init() {
  fTEX.close();
  cout << Form("/bin/rm -f %s", fTexFileName.c_str()) << endl;
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  cout << Form("open for TeX output: %s", fTexFileName.c_str()) << endl;
  fTEX.open(fTexFileName.c_str(), ios::app);

}

// ----------------------------------------------------------------------
void plotBDT::makeAll(string what) {

  if ("all" == what) {
    setBdtStrings(0);
    dumpParameters("Events0");

    string old("nada");
    // -- loop over channels in case they use separate BDTs
    for (int ichan = 0; ichan < fNchan; ++ichan) {
      setBdtStrings(ichan);
      if (old != fCuts[ichan]->bdtXml) {
	old = fCuts[ichan]->bdtXml;
	tmvaControlPlots(ichan);
	getTLFRanking("BDT", "Events0");
	getTLFRanking("BDT", "Events1");
	getTLFRanking("BDT", "Events2");
      } else {
	cout << "SKIPPING " << ichan << " because it is the same as the previous one!" << endl;
      }
    }
    getTLFParameters();
    getTLFEventNumbers();
  }

  if (("all" == what) || ("apply" == what)) {
    apply("fill");
  }

  if (("all" == what) || ("opspt" == what)) {
    opspt("hi");
  }

}

// ----------------------------------------------------------------------
void plotBDT::setBdtStrings(int ichan) {
  string XmlName = fCuts[ichan]->bdtXml;
  fBdtString = XmlName.substr(0, fBdtString.find("events"));
  fBdtLogFile = "weights/" + fCuts[ichan]->bdtXml + ".log";
  if (0 == fLogFileLines.size()) readLogFile(fBdtLogFile);
  cout << "fBdtString: " << fBdtString << endl;
  cout << "fBdtLogFile: " << fBdtLogFile << ", with " << fLogFileLines.size() << " lines" << endl;
}


// ----------------------------------------------------------------------
void plotBDT::readLogFile(string fname) {
  char  buffer[2000];
  cout << "getRanking: open file " << fname << endl;
  ifstream is(Form("%s", fname.c_str()));
  while (is.getline(buffer, 2000, '\n')) fLogFileLines.push_back(string(buffer));
  cout << "read " << fLogFileLines.size() << " lines" << endl;
  is.close();
}


// ----------------------------------------------------------------------
void plotBDT::tmvaControlPlots(int ichan) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  string type[3] = {"Events0", "Events1", "Events2"};
  for (int j = 0; j < 3; ++j) {
    string rootfile = "weights/" + fBdtString + "-" + type[j] + ".root";
    fRootFile = TFile::Open(rootfile.c_str());
    cout << "fRootFile: " << rootfile << endl;

    dumpParameters(type[j]);
    tmvaPlots(type[j]);
    fRootFile->Close();
  }

}


// ----------------------------------------------------------------------
void plotBDT::dumpParameters(string type) {
  fTEX << "% -- dumpParameters " << fBdtString << " " << fBdtLogFile << endl;
  TH1D *h = getPreselectionNumbers();
  for (int i = 1; i < h->GetNbinsX(); ++i) {
    if (!strcmp("", h->GetXaxis()->GetBinLabel(i))) continue;
    fTEX << formatTex(h->GetBinContent(i), Form("%s:%s:%s",  fSuffix.c_str(), fBdtString.c_str(), h->GetXaxis()->GetBinLabel(i)), 2) << endl;
  }
}


// ----------------------------------------------------------------------
void plotBDT::tmvaPlots(string type) {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  zone();
  shrinkPad(0.15, 0.15, 0.15);
  TH2 *h2 = (TH2*)fRootFile->Get("CorrelationMatrixS");
  h2->SetLabelSize(0.05, "x");
  h2->SetLabelSize(0.05, "y");
  h2->GetXaxis()->LabelsOption("v");

  h2->Draw("colztext");
  c0->SaveAs(Form("%s/%s-%s-CorrelationMatrixS.pdf", fDirectory.c_str(), fBdtString.c_str(), type.c_str()));

  h2 = (TH2*)fRootFile->Get("CorrelationMatrixB");
  h2->SetLabelSize(0.05, "x");
  h2->SetLabelSize(0.05, "y");
  h2->GetXaxis()->LabelsOption("v");
  h2->Draw("colztext");
  c0->SaveAs(Form("%s/%s-%s-CorrelationMatrixB.pdf", fDirectory.c_str(), fBdtString.c_str(), type.c_str()));

  // -- from "mvas"
  // --------------
  Int_t width = 600;   // size of canvas

  // this defines how many canvases we need
  TCanvas *c = new TCanvas( Form("canvas%d", 1), "canvas1",  200, 20, width, static_cast<int>(width*0.78) );

  // search for the right histograms in full list of keys
  TIter next(fRootFile->GetListOfKeys());
  TKey *key(0);
  while ((key = (TKey*)next())) {

    if (!TString(key->GetName()).BeginsWith("Method_")) continue;
    if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TDirectory")) continue;

    TString methodName;
    GetMethodName(methodName,key);

    TDirectory* mDir = (TDirectory*)key->ReadObj();

    TIter keyIt(mDir->GetListOfKeys());
    TKey *titkey;
    while ((titkey = (TKey*)keyIt())) {

      if (!gROOT->GetClass(titkey->GetClassName())->InheritsFrom("TDirectory")) continue;
      c->Clear();

      TDirectory *titDir = (TDirectory *)titkey->ReadObj();
      TString methodTitle;
      GetMethodTitle(methodTitle,titDir);

      cout << "--- Found directory for method: " << methodName << "::" << methodTitle << flush;
      TString hname = "MVA_" + methodTitle;
      TH1* sig = dynamic_cast<TH1*>(titDir->Get( hname + "_S" ));
      TH1* bgd = dynamic_cast<TH1*>(titDir->Get( hname + "_B" ));

      cout << " containing " << hname << "_S/_B" << endl;
      // chop off useless stuff
      sig->SetTitle( Form("TMVA overtraining check for classifier: %s", methodTitle.Data()) );

      // create new canvas
      TString ctitle = Form("TMVA comparison %s",methodTitle.Data()) ;

      // set the histogram style
      SetSignalAndBackgroundStyle( sig, bgd );

      // normalise both signal and background
      NormalizeHists( sig, bgd );

      // frame limits (choose judicuous x range)
      Float_t nrms = 10;
      cout << "--- Mean and RMS (S): " << sig->GetMean() << ", " << sig->GetRMS() << endl;
      cout << "--- Mean and RMS (B): " << bgd->GetMean() << ", " << bgd->GetRMS() << endl;
      Float_t xmin = TMath::Max( TMath::Min(sig->GetMean() - nrms*sig->GetRMS(),
                                            bgd->GetMean() - nrms*bgd->GetRMS() ),
                                 sig->GetXaxis()->GetXmin() );
      Float_t xmax = TMath::Min( TMath::Max(sig->GetMean() + nrms*sig->GetRMS(),
                                            bgd->GetMean() + nrms*bgd->GetRMS() ),
                                 sig->GetXaxis()->GetXmax() );
      Float_t ymin = 0;
      Float_t maxMult = 2.0;
      Float_t ymax = TMath::Max( sig->GetMaximum(), bgd->GetMaximum() )*maxMult;

      xmin = -1.;
      xmax = 1.;

      // build a frame
      Int_t nb = 500;
      TString hFrameName(TString("frame") + methodTitle);
      TObject *o = gROOT->FindObject(hFrameName);
      if(o) delete o;
      TH2F* frame = new TH2F( hFrameName, sig->GetTitle(),
                              nb, xmin, xmax, nb, ymin, ymax );
      frame->GetXaxis()->SetTitle( methodTitle + " response"  );
      frame->GetYaxis()->SetTitle("(1/N) dN^{ }/^{ }dx");
      SetFrameStyle( frame );

      // eventually: draw the frame
      frame->Draw();

      c->GetPad(0)->SetLeftMargin( 0.105 );
      frame->GetYaxis()->SetTitleOffset( 1.2 );

      // Draw legend
      TLegend *legend= new TLegend( c->GetLeftMargin(), 1 - c->GetTopMargin() - 0.12,
                                    c->GetLeftMargin() +  0.40, 1 - c->GetTopMargin() );
      legend->SetFillStyle( 1 );
      legend->AddEntry(sig,TString("Signal")     + " (test sample)", "F");
      legend->AddEntry(bgd,TString("Background") + " (test sample)", "F");
      legend->SetBorderSize(1);
      legend->SetMargin(0.2);
      legend->Draw("same");

      // overlay signal and background histograms
      sig->Draw("samehist");
      bgd->Draw("samehist");

      TH1* sigOv = 0;
      TH1* bgdOv = 0;

      TString ovname = hname += "_Train";
      sigOv = dynamic_cast<TH1*>(titDir->Get( ovname + "_S" ));
      bgdOv = dynamic_cast<TH1*>(titDir->Get( ovname + "_B" ));

      if (sigOv == 0 || bgdOv == 0) {
        cout << "+++ Problem in \"mvas.C\": overtraining check histograms do not exist" << endl;
      }
      else {
        cout << "--- Found comparison histograms for overtraining check" << endl;

        TLegend *legend2= new TLegend( 1 - c->GetRightMargin() - 0.42, 1 - c->GetTopMargin() - 0.12,
                                       1 - c->GetRightMargin(), 1 - c->GetTopMargin() );
        legend2->SetFillStyle( 1 );
        legend2->SetBorderSize(1);
        legend2->AddEntry(sigOv,"Signal (training sample)","P");
        legend2->AddEntry(bgdOv,"Background (training sample)","P");
        legend2->SetMargin( 0.1 );
        legend2->Draw("same");
      }
      // normalise both signal and background
      NormalizeHists( sigOv, bgdOv );

      Int_t col = sig->GetLineColor();
      sigOv->SetMarkerColor( col );
      sigOv->SetMarkerSize( 0.7 );
      sigOv->SetMarkerStyle( 20 );
      sigOv->SetLineWidth( 1 );
      sigOv->SetLineColor( col );
      sigOv->Draw("e1same");

      col = bgd->GetLineColor();
      bgdOv->SetMarkerColor( col );
      bgdOv->SetMarkerSize( 0.7 );
      bgdOv->SetMarkerStyle( 20 );
      bgdOv->SetLineWidth( 1 );
      bgdOv->SetLineColor( col );
      bgdOv->Draw("e1same");

      ymax = TMath::Max(ymax,
                        TMath::Max( static_cast<Float_t>(sigOv->GetMaximum()), static_cast<Float_t>(bgdOv->GetMaximum()) )*maxMult);
      frame->GetYaxis()->SetLimits( 0, ymax );

      // for better visibility, plot thinner lines
      sig->SetLineWidth( 1 );
      bgd->SetLineWidth( 1 );

      // perform K-S test
      cout << "--- Perform Kolmogorov-Smirnov tests" << endl;
      Double_t kolS = sig->KolmogorovTest( sigOv );
      Double_t kolB = bgd->KolmogorovTest( bgdOv );
      cout << "--- Goodness of signal (background) consistency: " << kolS << " (" << kolB << ")" << endl;
      fTEX << formatTex(kolS, Form("%s:%s:KS:Sg",  fSuffix.c_str(), fBdtString.c_str()), 3) << endl;
      fTEX << formatTex(kolB, Form("%s:%s:KS:Bg",  fSuffix.c_str(), fBdtString.c_str()), 3) << endl;

      TString probatext = Form( "Kolmogorov-Smirnov test: signal (background) probability = %5.3g (%5.3g)", kolS, kolB );
      TText* tt = new TText( 0.12, 0.96, probatext );
      tt->SetNDC(); tt->SetTextSize( 0.032 ); tt->AppendPad();

      // redraw axes
      frame->Draw("sameaxis");

      // text for overflows
      Int_t    nbin = sig->GetNbinsX();
      Double_t dxu  = sig->GetBinWidth(0);
      Double_t dxo  = sig->GetBinWidth(nbin+1);
      TString uoflow = Form( "U/O-flow (S,B): (%.1f, %.1f)%% / (%.1f, %.1f)%%",
                             sig->GetBinContent(0)*dxu*100, bgd->GetBinContent(0)*dxu*100,
                             sig->GetBinContent(nbin+1)*dxo*100, bgd->GetBinContent(nbin+1)*dxo*100 );
      TText* t = new TText( 0.975, 0.115, uoflow );
      t->SetNDC();
      t->SetTextSize( 0.030 );
      t->SetTextAngle( 90 );
      t->AppendPad();

      c->Update();

      // save canvas to file
      c->SaveAs(Form("%s/%s-%s-overtrain0.pdf", fDirectory.c_str(), fBdtString.c_str(), type.c_str()));
      frame->GetYaxis()->SetLimits(1.e-5, 5.e1);
      c->SetLogy(1);
      c->SaveAs(Form("%s/%s-%s-overtrain1.pdf", fDirectory.c_str(), fBdtString.c_str(), type.c_str()));

    }
    cout << "";
  }



  // -- from "variables"
  // -------------------

  c->SetLogy(0);

  TString title = "TMVA Input Variables";
  // obtain shorter histogram title
  TString htitle = title;
  htitle.ReplaceAll("variables ","variable");
  htitle.ReplaceAll("and target(s)","");
  htitle.ReplaceAll("(training sample)","");

  TString dirName = "Method_BDT/BDT";
  TDirectory* dir = (TDirectory*)fRootFile->Get(dirName);
  if (dir==0) {
    cout << "No information about " << title << " available in directory " << dirName << " of file " << fRootFile << endl;
    return;
  }
  dir->cd();


  // -- create empty pdfs to be overwritten for the real variables
  c->Clear();
  vector<string> allvars;
  allvars.push_back("pt");
  allvars.push_back("eta");
  allvars.push_back("pvips");
  allvars.push_back("pvip");
  allvars.push_back("alpha");
  allvars.push_back("chi2dof");
  allvars.push_back("maxdoca");
  allvars.push_back("fls3d");
  allvars.push_back("closetrk");
  allvars.push_back("docatrk");
  allvars.push_back("iso");
  allvars.push_back("m1iso");
  allvars.push_back("m2iso");
  for (unsigned int i = 0; i < allvars.size(); ++i) {
    c->Clear();
    c->SaveAs(Form("%s/%s-%s-%s.pdf", fDirectory.c_str(), fBdtString.c_str(), type.c_str(), allvars[i].c_str()));
  }



  // loop over all objects in directory
  TIter next1(dir->GetListOfKeys());
  while ((key = (TKey*)next1())) {
    // make sure, that we only look at histograms
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1 *sig = (TH1*)key->ReadObj();
    TString hname(sig->GetName());
    if (hname.Contains("__Background")) continue;
    TString htitle(sig->GetTitle());
    if (htitle.Contains("Input Variables")) continue;
    if (htitle.Contains(" ")) continue;

    if (htitle.Contains("MVA_BDT")) {
      //  cout << "--> SKIPPING " << htitle << endl;
      continue;
    }

    // find the corredponding backgrouns histo
    TString bgname = hname;
    bgname.ReplaceAll("__Signal","__Background");
    TH1 *bgd = (TH1*)dir->Get(bgname);
    if (bgd == NULL) {
      cout << "ERROR!!! couldn't find background histo for" << hname << endl;
      return;
    }

    // this is set but not stored during plot creation in MVA_Factory
    SetSignalAndBackgroundStyle(sig, bgd);

    sig->SetTitle(TString( htitle ) + ": " + sig->GetTitle() );
    SetFrameStyle(sig, 1.2);

    // normalise both signal and background
    NormalizeHists(sig, bgd);

    // finally plot and overlay
    Float_t sc = 1.4;
    sig->SetMaximum( TMath::Max( sig->GetMaximum(), bgd->GetMaximum() )*sc );
    sig->Draw( "hist" );
    gPad->SetLeftMargin( 0.17 );

    sig->GetYaxis()->SetTitleOffset( 1.70 );
    bgd->Draw("histsame");
    TString ytit = TString("(1/N) ") + sig->GetYaxis()->GetTitle();
    sig->GetYaxis()->SetTitle( ytit ); // histograms are normalised

    // Draw legend
    TLegend *legend= new TLegend( gPad->GetLeftMargin(),
                                  1-gPad->GetTopMargin()-.15,
                                  gPad->GetLeftMargin()+.4,
                                  1-gPad->GetTopMargin() );
    legend->SetFillStyle(1);
    legend->AddEntry(sig,"Signal","F");
    legend->AddEntry(bgd,"Background","F");
    legend->SetBorderSize(1);
    legend->SetMargin( 0.3 );
    legend->Draw("same");

    // redraw axes
    sig->Draw("sameaxis");

    // text for overflows
    Int_t    nbin = sig->GetNbinsX();
    Double_t dxu  = sig->GetBinWidth(0);
    Double_t dxo  = sig->GetBinWidth(nbin+1);
    TString uoflow = Form( "U/O-flow (S,B): (%.1f, %.1f)%% / (%.1f, %.1f)%%",
                           sig->GetBinContent(0)*dxu*100, bgd->GetBinContent(0)*dxu*100,
                           sig->GetBinContent(nbin+1)*dxo*100, bgd->GetBinContent(nbin+1)*dxo*100 );

    TText* t = new TText( 0.98, 0.14, uoflow );
    t->SetNDC();
    t->SetTextSize( 0.040 );
    t->SetTextAngle( 90 );
    t->AppendPad();

    c->SaveAs(Form("%s/%s-%s-%s.pdf", fDirectory.c_str(), fBdtString.c_str(), type.c_str(), htitle.Data()));

    delete legend;
  }




  // "BoostMonitor","BoostWeight","BoostWeightVsTree","ErrFractHist","NodesBeforePruning"
  dirName = "Method_BDT/BDT";
  dir = (TDirectory*)fRootFile->Get(dirName);
  if (dir==0) {
    cout << "No information about " << title << " available in directory " << dirName << " of file " << fRootFile << endl;
    return;
  }
  dir->cd();

  TCanvas *cc = new TCanvas("cc", "", 300, 200, 1000, 400);

  cc->cd();
  shrinkPad(0.20, 0.15, 0.1, 0.);
  TH1F *hf = (TH1F*)((TH1F*)dir->Get("BoostWeightVsTree"))->Clone("hf");
  setTitles(hf, hf->GetXaxis()->GetTitle(), hf->GetYaxis()->GetTitle(), 0.09, 1.02, 0.8, 0.09);
  hf->Draw();
  cc->SaveAs(Form("%s/%s-%s-BoostWeightVsTree.pdf", fDirectory.c_str(), fBdtString.c_str(), type.c_str()));

  cc->Clear();
  shrinkPad(0.20, 0.15, 0.1, 0.);
  delete hf;
  hf = (TH1F*)((TH1F*)dir->Get("ErrFractHist"))->Clone("hf");
  setTitles(hf, hf->GetXaxis()->GetTitle(), hf->GetYaxis()->GetTitle(), 0.09, 1.02, 0.8, 0.09);
  hf->Draw();
  cc->SaveAs(Form("%s/%s-%s-ErrFractHist.pdf", fDirectory.c_str(), fBdtString.c_str(), type.c_str()));

  cc->Clear();
  shrinkPad(0.20, 0.15, 0.1, 0.);
  delete hf;
  hf = (TH1F*)((TH1F*)dir->Get("NodesBeforePruning"))->Clone("hf");
  setTitles(hf, hf->GetXaxis()->GetTitle(), hf->GetYaxis()->GetTitle(), 0.09, 1.02, 0.8, 0.09);
  hf->Draw();
  cc->SaveAs(Form("%s/%s-%s-NodesBeforePruning.pdf", fDirectory.c_str(), fBdtString.c_str(), type.c_str()));

}




// ----------------------------------------------------------------------
void plotBDT::loopFunction1() {

  if (fChan < 0) return;

  //  if ((1 == fMode) && ((fb.m < 4.9) || (fb.m > 5.9) || ((5.2 < fb.m) && fb.m < 5.45))) return;
  if ((1 == fMode) && ((fb.m < 4.9) || (fb.m > 5.9))) return;

  if (preselection(fb)) {
    frd.pt = fb.pt;
    frd.eta = fb.eta;
    frd.m1eta = fb.m1eta;
    frd.m2eta = fb.m2eta;
    frd.m1pt = fb.m1pt;
    frd.m2pt = fb.m2pt;
    frd.fls3d = fb.fls3d;
    frd.alpha = fb.alpha;
    frd.maxdoca = fb.maxdoca;
    frd.pvip = fb.pvip;
    frd.pvips = fb.pvips;
    frd.iso = fb.iso;
    frd.docatrk = fb.docatrk;
    frd.chi2dof = fb.chi2dof;
    frd.closetrk = fb.closetrk;

    frd.m1iso = fb.m1iso;
    frd.m2iso = fb.m2iso;

    frd.closetrks1 = fb.closetrks1;
    frd.closetrks2 = fb.closetrks2;
    frd.closetrks3 = fb.closetrks3;

    frd.pv2lip  = fb.pv2lip;
    frd.pv2lips = fb.pv2lips;

    frd.m  = fb.m;
    fBDT0   = fReaderEvents0[fChan]->EvaluateMVA("BDT");
    fBDT1   = fReaderEvents1[fChan]->EvaluateMVA("BDT");
    fBDT2   = fReaderEvents2[fChan]->EvaluateMVA("BDT");

    int remainder = TMath::Abs(fb.evt%3);
    calcBDT();

    fHists[Form("cand%dChan%dBdt0Evt%d", fMode, fChan, remainder)]->Fill(fb.m, fBDT0);
    fHists[Form("cand%dChan%dBdt1Evt%d", fMode, fChan, remainder)]->Fill(fb.m, fBDT1);
    fHists[Form("cand%dChan%dBdt2Evt%d", fMode, fChan, remainder)]->Fill(fb.m, fBDT2);
    fHists[Form("cand%dChan%dBdtEvt%d",  fMode, fChan, remainder)]->Fill(fb.m, fBDT);
  } else {
    fBDT0 = -10.;
    fBDT1 = -10.;
    fBDT2 = -10.;
  }


}

// ----------------------------------------------------------------------
void plotBDT::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
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
  cout << "==> plotBDT::loopOverTree> loop over dataset " << (fCds?fCds->fName:"undefined") << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << " nbegin = " << nbegin << " nend = " << nend
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotBDT::*pF)(void);
  if (ifunc == 1) pF = &plotBDT::loopFunction1;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;
    candAnalysis();
    (this->*pF)();
  }

}


// ----------------------------------------------------------------------
// -- stuff from tmvaglob.C
void plotBDT::GetMethodName( TString & name, TKey * mkey ) {
  if (mkey==0) return;
  name = mkey->GetName();
  name.ReplaceAll("Method_","");
}

void plotBDT::GetMethodTitle( TString & name, TKey * ikey ) {
  if (ikey==0) return;
  name = ikey->GetName();
}

void plotBDT::GetMethodTitle( TString & name, TDirectory * idir ) {
  if (idir==0) return;
  name = idir->GetName();
}


void plotBDT::SetSignalAndBackgroundStyle( TH1* sig, TH1* bkg, TH1* all )    {
      //signal
      // const Int_t FillColor__S = 38 + 150; // change of Color Scheme in ROOT-5.16.
      // convince yourself with gROOT->GetListOfColors()->Print()
  static Int_t c_SignalLine     = TColor::GetColor( "#0000ee" );
  static Int_t c_SignalFill     = TColor::GetColor( "#7d99d1" );
  static Int_t c_BackgroundLine = TColor::GetColor( "#ff0000" );
  static Int_t c_BackgroundFill = TColor::GetColor( "#ff0000" );

      Int_t FillColor__S = c_SignalFill;
      Int_t FillStyle__S = 3356; //1001;
      Int_t LineColor__S = c_SignalLine;
      Int_t LineWidth__S = 2;

      // background
      //Int_t icolor = UsePaperStyle ? 2 + 100 : 2;
      Int_t FillColor__B = c_BackgroundFill;
      Int_t FillStyle__B = 3365; //3554;
      Int_t LineColor__B = c_BackgroundLine;
      Int_t LineWidth__B = 2;

      Int_t FillColor__A = kBlack;
      Int_t FillStyle__A = 3350; //3554;
      Int_t LineColor__A = kBlack;
      Int_t LineWidth__A = 2;

      if (sig != NULL) {
         sig->SetMarkerColor(LineColor__S );
         sig->SetLineColor( LineColor__S );
         sig->SetLineWidth( LineWidth__S );
         sig->SetFillStyle( FillStyle__S );
         sig->SetFillColor( FillColor__S );
      }

      if (bkg != NULL) {
	bkg->SetMarkerColor(LineColor__B );
	bkg->SetLineColor( LineColor__B );
	bkg->SetLineWidth( LineWidth__B );
	bkg->SetFillStyle( FillStyle__B );
	bkg->SetFillColor( FillColor__B );
      }

      if (all != NULL) {
	all->SetMarkerColor(LineColor__A );
	all->SetLineColor( LineColor__A );
	all->SetLineWidth( LineWidth__A );
	all->SetFillStyle( FillStyle__A );
	all->SetFillColor( FillColor__A );
      }
   }

void plotBDT::NormalizeHists( TH1* sig, TH1* bkg )   {
      if (sig->GetSumw2N() == 0) sig->Sumw2();
      if (bkg && bkg->GetSumw2N() == 0) bkg->Sumw2();

      if(sig->GetSumOfWeights()!=0) {
         Float_t dx = (sig->GetXaxis()->GetXmax() - sig->GetXaxis()->GetXmin())/sig->GetNbinsX();
         sig->Scale( 1.0/sig->GetSumOfWeights()/dx );
      }
      if (bkg != 0 && bkg->GetSumOfWeights()!=0) {
         Float_t dx = (bkg->GetXaxis()->GetXmax() - bkg->GetXaxis()->GetXmin())/bkg->GetNbinsX();
         bkg->Scale( 1.0/bkg->GetSumOfWeights()/dx );
      }
   }

void plotBDT::SetFrameStyle( TH1* frame, Float_t scale)   {
  frame->SetLabelOffset( 0.012, "X" );// label offset on x axis
  frame->SetLabelOffset( 0.012, "Y" );// label offset on x axis
  frame->GetXaxis()->SetTitleOffset( 1.25 );
  frame->GetYaxis()->SetTitleOffset( 1.22 );
  frame->GetXaxis()->SetTitleSize( 0.045*scale );
  frame->GetYaxis()->SetTitleSize( 0.045*scale );
  Float_t labelSize = 0.04*scale;
  frame->GetXaxis()->SetLabelSize( labelSize );
  frame->GetYaxis()->SetLabelSize( labelSize );

  // global style settings
  gPad->SetTicks();
  gPad->SetLeftMargin  ( 0.108*scale );
  gPad->SetRightMargin ( 0.050*scale );
  gPad->SetBottomMargin( 0.120*scale  );
}


int plotBDT::GetNumberOfTargets( TDirectory *dir )   {
      TIter next(dir->GetListOfKeys());
      TKey* key    = 0;
      Int_t noTrgts = 0;

      while ((key = (TKey*)next())) {
         if (key->GetCycle() != 1) continue;
         if (TString(key->GetName()).Contains("__Regression_target")) noTrgts++;
      }
      return noTrgts;
}

int plotBDT::GetNumberOfInputVariables( TDirectory *dir )   {
  TIter next(dir->GetListOfKeys());
  TKey* key    = 0;
  Int_t noVars = 0;

  while ((key = (TKey*)next())) {
    if (key->GetCycle() != 1) continue;

    // count number of variables (signal is sufficient), exclude target(s)
    if (TString(key->GetName()).Contains("__Signal") || (TString(key->GetName()).Contains("__Regression") && !(TString(key->GetName()).Contains("__Regression_target")))) noVars++;
  }

  return noVars;
}



// ----------------------------------------------------------------------
int plotBDT::bdtString2Channel(string s) {
  int channel(-1);
  if (string::npos != s.find("TMVA-0")) channel = 0;
  if (string::npos != s.find("TMVA-1")) channel = 1;
  if (string::npos != s.find("TMVA-2")) channel = 0;
  if (string::npos != s.find("TMVA-3")) channel = 1;
  if (string::npos != s.find("TMVA-4")) channel = 1;
  return channel;
}


// ----------------------------------------------------------------------
void plotBDT::opspt(string what) {
  fHistFile = TFile::Open(fHistFileName.c_str());

  TH2D *h2(0), *hCombBg(0), *hCombSg(0);

  h2 = (TH2D*)fHistFile->Get("cand1Chan0Bdt0Evt0");
  hCombSg = (TH2D*)h2->Clone("hCombSg");
  hCombSg->Reset();
  hCombBg = (TH2D*)h2->Clone("hCombBg");
  hCombBg->Reset();

  int nchan(2);
  map<string, TH1D*> hSSBs, hSSBd;
  zone(2,2);
  double sLumi(0.), bLumi(0.);
  for (int ibdt = 0; ibdt < 3; ++ibdt) {
    for (int ievt = 0; ievt < 3; ++ievt) {
      hCombSg->Reset();
      hCombBg->Reset();
      for (int ichan = 0; ichan < nchan; ++ichan) {
	h2 = (TH2D*)fHistFile->Get(Form("cand1Chan%dBdt%dEvt%d", ichan, ibdt, ievt));
	hCombBg->Add(h2);
	if (0 == ichan) {
	  string hname = h2->GetTitle();
	  string pd = hname.substr(hname.rfind("PD:") + 4);
	  bLumi = fDS[pd]->fLumi;
	  hname = h2->GetName();
	  replaceAll(hname, "cand1Chan0", "Bg");
	  hCombBg->SetName(hname.c_str());
	}


	h2 = (TH2D*)fHistFile->Get(Form("cand3Chan%dBdt%dEvt%d", ichan, ibdt, ievt));
	hCombSg->Add(h2);
	if (0 == ichan) {
	  string hname = h2->GetTitle();
	  string pd = hname.substr(hname.rfind("PD:") + 4);
	  sLumi = fDS[pd]->fLumi;
	  hname = h2->GetName();
	  replaceAll(hname, "cand1Chan0", "Sg");
	  hCombSg->SetName(hname.c_str());
	}
      }
      // -- scale signal with lumi ratio
      hCombSg->Scale(bLumi/sLumi);

      // -- scan BDT cut
      hSSBs.insert(make_pair(Form("bdt%d_evt%d", ibdt, ievt), ssb(hCombSg, hCombBg, what)));
      TH2D *hCombB0 = (TH2D*)hCombSg->Clone(Form("B0_%s", hCombSg->GetName()));
      hCombB0->Scale(0.1);
      hSSBd.insert(make_pair(Form("bdt%d_evt%d", ibdt, ievt), ssb(hCombB0, hCombBg, what)));
    }
  }


  zone(3,3);
  map<string, TH1D*>::iterator is = hSSBs.begin();
  map<string, TH1D*>::iterator id = hSSBd.begin();
  int cnt(1);
  TF1 *f1(0), *f2(0);
  TH1D *h(0);
  double xmax(0.), rms(0.);
  gStyle->SetOptFit(0);
  for (; is != hSSBs.end(); ++is, ++id) {
    c0->cd(cnt);
    shrinkPad(0.12, 0.12, 0.1, 0.12);
    h = is->second;
    setTitles(h, "BDT>", "S/#sqrt{S+B}");
    h->SetLineWidth(2);
    f1 = fIF->pol2local(h);
    xmax = h->GetBinCenter(h->GetMaximumBin());
    rms  = 0.4*h->GetRMS();
    h->Fit(f1, "", "", xmax-rms, xmax+rms);
    tl->SetTextColor(kBlack);
    tl->DrawLatexNDC(0.2, 0.7, Form("B_{max}(B_{s}) = %3.2f", f1->GetParameter(2)));

    h = id->second;
    h->SetLineColor(kBlue);
    h->SetLineWidth(2);
    f2 = fIF->pol2local(h);
    f2->SetLineColor(kCyan);
    xmax = h->GetBinCenter(h->GetMaximumBin());
    rms  = 0.4*h->GetRMS();
    h->Fit(f2, "same", "same", xmax-rms, xmax+rms);
    tl->SetTextColor(kBlue);
    tl->DrawLatexNDC(0.2, 0.6, Form("B_{max}(B^{0}) = %3.2f", f2->GetParameter(2)));

    tl->SetTextColor(kBlack);
    tl->DrawLatexNDC(0.2, 0.8, is->first.c_str());
    ++cnt;
  }


  savePad(Form("%s-bdtops-%s.pdf", fBdtString.c_str(), what.c_str()), c0);

}


// ----------------------------------------------------------------------
TH1D* plotBDT::ssb(TH2D *hS, TH2D *hB, string what) {
  TH1D *hssb = hS->ProjectionY(Form("hSSB_%s", hS->GetName()));
  hssb->Reset();

  TH1D *hs(0), *hb(0);
  cout << "pointers: " << hS << " " << hB << " nbins = " << hS->GetNbinsY() << endl;
  double nS(-1.), nB(-1.);
  double scale(0.);
  for (int ibdt = 1; ibdt < hS->GetNbinsY(); ++ibdt) {
    hs = hS->ProjectionX("hs", ibdt, hS->GetNbinsY());
    nS = hs->Integral();

    hb = hB->ProjectionX("hb", ibdt, hB->GetNbinsY());
    if (1 == ibdt) {
      scale = hb->Integral(hb->FindBin(5.45001), hb->GetNbinsX())/hb->Integral(hb->FindBin(5.20001), hb->FindBin(5.4999));
      scale *= (5.45-5.20)/(5.9-5.45);
    }
    if ("hi" == what) {
      for (int ibin = 1; ibin <= hb->FindBin(5.44999); ++ibin) {
	hb->SetBinContent(ibin, 0.);
	hb->SetBinError(ibin, 0.);
      }
      nB = hb->Integral()*(5.45-5.20)/(5.9-5.45);
    } else if ("sidebands" == what) {
      for (int ibin = hb->FindBin(5.200001); ibin < hb->FindBin(5.450001); ++ibin) {
	hb->SetBinContent(ibin, 0.);
	hb->SetBinError(ibin, 0.);
      }
    } else if ("signalregion" == what) {
      for (int ibin = 1; ibin <= hb->FindBin(5.199999); ++ibin) {
	hb->SetBinContent(ibin, 0.);
	hb->SetBinError(ibin, 0.);
      }
      double hiIntegral(0.);
      int hiCnt(0);
      for (int ibin = hb->FindBin(5.450001); ibin <= hb->GetNbinsX(); ++ibin) {
	++hiCnt;
	hiIntegral += hb->GetBinContent(ibin);
	hb->SetBinContent(ibin, 0.);
	hb->SetBinError(ibin, 0.);
      }
      cout << " hB integral: " << hb->Integral() << " hiIntegral/hiCnt = " << hiIntegral/hiCnt
	   << " (hb->FindBin(5.44999)-hb->FindBin(5.20001)) = " << (hb->FindBin(5.44999)-hb->FindBin(5.20001)) << endl;
      if (hb->Integral() > 0) {
	nB = hb->Integral() * scale;
      } else {
	nB = hb->Integral();
      }
    }
    cout << "ibdt = " << ibdt << " nS = " << nS << " nB = " << nB << endl;
    if (nS+nB > 0) {
      hssb->SetBinContent(ibdt, nS/TMath::Sqrt(nS+nB));
    }
  }

  return hssb;
}


// ----------------------------------------------------------------------
void plotBDT::apply(string what) {

  if (what == "fill") {
    fHistFile = TFile::Open(fHistFileName.c_str(), "RECREATE");

    TDirectory *hDir = gDirectory;

    string pd = "bsmmMcOff";
    setup(pd);
    TTree *t = getTree(fSetup, fTreeDir);
    if (0 == t) {
      cout << "tree for sample = " << fSetup << " not found" << endl;
      return;
    }
    setupTree(t, fSetup);
    fCds = fDS[fSetup];
    cout << "fSetup = " << fSetup << " with fMode = " << fMode << endl;
    int nchan(5);
    for (int ichan = 0; ichan < nchan; ++ichan) {
      for (int ibdt = 0; ibdt < 3; ++ibdt) {
	for (int jevt = 0; jevt < 3; ++jevt) {
	  fHists.insert(make_pair(Form("cand%dChan%dBdt%dEvt%d", fMode, ichan, ibdt, jevt),
				  new TH2D(Form("cand%dChan%dBdt%dEvt%d", fMode, ichan, ibdt, jevt),
					   Form("cand%dChan%dBdt%dEvt%d PD: %s", fMode, ichan, ibdt, jevt, pd.c_str()),
					   50, 4.9, 5.9, 120, -0.6, 0.6)));
	  if (0 == ibdt) {
	    fHists.insert(make_pair(Form("cand%dChan%dBdtEvt%d", fMode, ichan, jevt),
				    new TH2D(Form("cand%dChan%dBdtEvt%d", fMode, ichan, jevt),
					     Form("cand%dChan%dBdtEvt%d PD: %s", fMode, ichan, jevt, pd.c_str()),
					     50, 4.9, 5.9, 120, -0.6, 0.6)));
	  }
	}
      }
    }
    loopOverTree(t, 1);

    pd = "bmmData";
    setup(pd);
    t = getTree(fSetup, fTreeDir);
    setupTree(t, fSetup);
    fCds = fDS[fSetup];
    cout << "fSetup = " << fSetup << " with fMode = " << fMode << endl;
    for (int ichan = 0; ichan < nchan; ++ichan) {
      for (int ibdt = 0; ibdt < 3; ++ibdt) {
	for (int jevt = 0; jevt < 3; ++jevt) {
	  fHists.insert(make_pair(Form("cand%dChan%dBdt%dEvt%d", fMode, ichan, ibdt, jevt),
				  new TH2D(Form("cand%dChan%dBdt%dEvt%d", fMode, ichan, ibdt, jevt),
					   Form("cand%dChan%dBdt%dEvt%d PD: %s", fMode, ichan, ibdt, jevt, pd.c_str()),
					   50, 4.9, 5.9, 120, -0.6, 0.6)));
	  if (0 == ibdt) {
	    fHists.insert(make_pair(Form("cand%dChan%dBdtEvt%d", fMode, ichan, jevt),
				    new TH2D(Form("cand%dChan%dBdtEvt%d", fMode, ichan, jevt),
					     Form("cand%dChan%dBdtEvt%d PD: %s", fMode, ichan, jevt, pd.c_str()),
					     50, 4.9, 5.9, 120, -0.6, 0.6)));
	  }
	}
      }
    }

    loopOverTree(t, 1);

    for (map<string, TH1*>::iterator it = fHists.begin(); it != fHists.end(); ++it) {
      cout << "name = " << it->first << " hist " << it->second->GetName() << endl;
      it->second->SetDirectory(hDir);
      it->second->Write();
      delete it->second;
    }
    fHistFile->Close();
    return;
  } else {
    fHistFile = TFile::Open(fHistFileName.c_str());

  }

}


// ----------------------------------------------------------------------
void plotBDT::getTLFEventNumbers() {
  // -- get number of events in all trees used for training/testing/validation
  int sg0(-1), sg1(-1), sg2(-1);
  int bg0(-1), bg1(-1), bg2(-1);
  string snum("");
  for (unsigned int i = 0; i < fLogFileLines.size(); ++i) {
    if (string::npos != fLogFileLines[i].find("==============> applySg =  signalAllEvents0")) {
      snum = fLogFileLines[i].substr(9 + fLogFileLines[i].rfind("entries: "));
      sg0 = atoi(snum.c_str());
    }
    if (string::npos != fLogFileLines[i].find("==============> trainSg =  signalAllEvents1")) {
      snum = fLogFileLines[i].substr(9 + fLogFileLines[i].rfind("entries: "));
      sg1 = atoi(snum.c_str());
    }
    if (string::npos != fLogFileLines[i].find("==============> testSg  =  signalAllEvents2")) {
      snum = fLogFileLines[i].substr(9 + fLogFileLines[i].rfind("entries: "));
      sg2 = atoi(snum.c_str());
    }
    if (string::npos != fLogFileLines[i].find("==============> applyBg =  sidebandAllEvents0")) {
      snum = fLogFileLines[i].substr(9 + fLogFileLines[i].rfind("entries: "));
      bg0 = atoi(snum.c_str());
    }
    if (string::npos != fLogFileLines[i].find("==============> trainBg =  sidebandAllEvents1")) {
      snum = fLogFileLines[i].substr(9 + fLogFileLines[i].rfind("entries: "));
      bg1 = atoi(snum.c_str());
    }
    if (string::npos != fLogFileLines[i].find("==============> testBg  =  sidebandAllEvents2")) {
      snum = fLogFileLines[i].substr(9 + fLogFileLines[i].rfind("entries: "));
      bg2 = atoi(snum.c_str());
    }

    if (sg0 > 0 && sg1 > 0 && sg2 > 0 && bg0 > 0 && bg1 > 0 && bg2 > 0) break;
  }
  fTEX << formatTex(sg0, Form("%s:%s-Events0:sgcnt",  fSuffix.c_str(), fBdtString.c_str()), 0) << endl;
  fTEX << formatTex(sg1, Form("%s:%s-Events1:sgcnt",  fSuffix.c_str(), fBdtString.c_str()), 0) << endl;
  fTEX << formatTex(sg2, Form("%s:%s-Events2:sgcnt",  fSuffix.c_str(), fBdtString.c_str()), 0) << endl;

  fTEX << formatTex(bg0, Form("%s:%s-Events0:bgcnt",  fSuffix.c_str(), fBdtString.c_str()), 0) << endl;
  fTEX << formatTex(bg1, Form("%s:%s-Events1:bgcnt",  fSuffix.c_str(), fBdtString.c_str()), 0) << endl;
  fTEX << formatTex(bg2, Form("%s:%s-Events2:bgcnt",  fSuffix.c_str(), fBdtString.c_str()), 0) << endl;

}

// ----------------------------------------------------------------------
void plotBDT::getTLFRanking(string prefix, string type) {
  cout << "% -- getTLFRanking: "  << prefix << " " << type << endl;
  fTEX << "% -- getTLFRanking: "  << prefix << " " << type << endl;
  string::size_type m1, m2;
  string varn, vars, varp;
  int bail(0), istart(0);

  string after;
  if (type == "Events0") {
    after = "Apply on events0, train on events1, test on events2";
  } else if (type == "Events1") {
    after = "Apply on events1, train on events2, test on events0";
  } else if (type == "Events2") {
    after = "Apply on events2, train on events0, test on events1";
  }

  for (unsigned int i = 0; i < fLogFileLines.size(); ++i) {
    if (string::npos != fLogFileLines[i].find(after)) {
      istart = i;
      break;
    }
  }

  cout << "start after: " << after << " at line " << istart << endl;
  int ibin(1);
  for (unsigned int i = istart; i < fLogFileLines.size(); ++i) {
    // -- method unspecific classification
    if ((string::npos != fLogFileLines[i].find(Form("--- %s", prefix.c_str())))
	&& (string::npos != fLogFileLines[i].find(": Rank : Variable "))
	) {
      bail = 0;
      for (unsigned int j = i+2; j < i+100; ++j) {
	if (string::npos != fLogFileLines[j].find(": ---------------------------------")) {
	  bail = 1;
	  cout << "  -> breaking out " << endl;
	  break;
	}

	m1 = fLogFileLines[j].find(":");
	m2 = fLogFileLines[j].find(":", m1+1);
	varn = fLogFileLines[j].substr(m1+2, m2-m1-2);
	m1 = m2;
	m2 = fLogFileLines[j].find(":", m1+1);
	vars = fLogFileLines[j].substr(m1+2, m2-m1-2);
	m1 = m2;
	m2 = fLogFileLines[j].find(":", m1+1);
	varp = fLogFileLines[j].substr(m1+2, m2-m1-2);
	cout << varn << "-> " << vars << " -> " << varp << endl;
	//	cout << fLogFileLines[j] << endl;
	replaceAll(vars, " " , "");
	string texlabel = fVarToTexSymbol[vars];
	cout << "       -> " << texlabel << endl;
	fTEX << Form("\\vdef{%s:%s-%s_%s:%d:name} {%s}",
		     fSuffix.c_str(), fBdtString.c_str(), type.c_str(), prefix.c_str(), ibin, texlabel.c_str())
	     << endl;
	fTEX << formatTex(atof(varp.c_str()), Form("%s:%s-%s_%s:%d:val",
						   fSuffix.c_str(), fBdtString.c_str(), type.c_str(), prefix.c_str(), ibin), 4)
	     << endl;
	++ibin;
      }
      if (1 == bail) break;
    }
  }
}


// ----------------------------------------------------------------------
void plotBDT::getTLFParameters(string prefix, string type) {
  cout << "% -- getTLFParameters: " << prefix << " " << type << endl;
  fTEX << "% -- getTLFParameters: " << prefix << " " << type << endl;


  string::size_type m1, m2;
  string varn, vars, varp;
  int bail(0), istart(0);
  string after;

  if (type == "Events0") {
    after = "Apply on events0, train on events1, test on events2";
  } else if (type == "Events1") {
    after = "Apply on events1, train on events2, test on events0";
  } else if (type == "Events2") {
    after = "Apply on events2, train on events0, test on events1";
  }

  for (unsigned int i = 0; i < fLogFileLines.size(); ++i) {
    if (string::npos != fLogFileLines[i].find(after)) {
      istart = i;
      break;
    }
  }
  cout << "start after: " << after << " at line " << istart << endl;

  string searchFor0 = "BDT";
  string searchFor1 = "The following options are set";
  for (unsigned int i = istart; i < fLogFileLines.size(); ++i) {
    if ((string::npos != fLogFileLines[i].find(searchFor0)) && (string::npos != fLogFileLines[i].find(searchFor1))){
      istart = i+1;
      break;
    }
  }

  cout << "start after: " << searchFor0 << "+" << searchFor1 << " at line " << istart << endl;
  for (unsigned int i = istart; i < fLogFileLines.size(); ++i) {
    if (string::npos != fLogFileLines[i].find("deprecated")) continue;
    if (string::npos != fLogFileLines[i].find(": - ")) continue;
    if (string::npos != fLogFileLines[i].find(" DataSetFactory           : Parsing option string:")) break;
    m1 = fLogFileLines[i].find(":");
    m2 = fLogFileLines[i].find(":", m1+1);
    varn = fLogFileLines[i].substr(m1+1, m2-m1-1);
    replaceAll(varn, " ", "");
    m1 = fLogFileLines[i].find("\"", m1);
    m2 = fLogFileLines[i].find("\"", m1+1);
    vars = fLogFileLines[i].substr(m1+1, m2-m1-1);
    replaceAll(vars, " ", "");
    replaceAll(vars, "%", "&");
    replaceAll(vars, "&", "\\%");
    if (varn == "V") continue;
    if (varn == "H") continue;
    if (varn == "VerbosityLevel") continue;
    if (varn == "CreateMVAPdfs") continue;
    cout << "varn: ->" << varn << "<- vars: ->" << vars << "<-" << endl;
    fTEX << Form("\\vdef{%s:%s:bdtPar_%s:name} {%s}",
		 fSuffix.c_str(), fBdtString.c_str(), varn.c_str(), varn.c_str())
	 << endl;
    fTEX << Form("\\vdef{%s:%s:bdtPar_%s:val} {%s}",
		 fSuffix.c_str(), fBdtString.c_str(), varn.c_str(), vars.c_str())
	 << endl;
  }
}
