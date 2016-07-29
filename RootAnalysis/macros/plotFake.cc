#include "plotFake.hh"

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

#include "common/dataset.hh"
#include "common/util.hh"

ClassImp(plotFake)

using namespace std;

// ----------------------------------------------------------------------
plotFake::plotFake(string dir, string files, string cuts, string setup): plotClass(dir, files, cuts, setup) {
  loadFiles(files);

  if (setup == "") {
    fHistFileName = Form("%s/plotFake.root", dir.c_str());
  } else {
    fHistFileName = Form("%s/plotFake-%s.root", dir.c_str(), setup.c_str());
  }

  fTexFileName = fHistFileName;
  replaceAll(fTexFileName, ".root", ".tex");
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));

}


// ----------------------------------------------------------------------
plotFake::~plotFake() {

}


// ----------------------------------------------------------------------
void plotFake::makeAll(int bitmask) {

}

// ----------------------------------------------------------------------
void plotFake::fakeRate(string var, string dataset, string particle) {

  tl->SetNDC(kTRUE);

  c0->Clear();
  //  c0->Divide(4, 4);
  c0->Divide(2, 3);
  int ipad(0);

  vector<string> mode;
  mode.push_back("nmu");
  //  mode.push_back("muo");

  map<string, TH1D*> hmode;

  vector<string> histos;
  if (particle == "pion") {
    histos.push_back(Form("candAnaFake310/%s1", var.c_str()));
    histos.push_back(Form("candAnaFake310/%s2", var.c_str()));
  } else if (particle == "kaon") {
    histos.push_back(Form("candAnaFake333/%s1", var.c_str()));
    //    histos.push_back(Form("candAnaFake333/%s2", var.c_str()));
  }  else if (particle == "proton") {
    histos.push_back(Form("candAnaFake3122/%s1", var.c_str()));
  } else {
    cout << "particle " << particle << " not known, returning" << endl;
  }

  // -- get "default" to properly initialize results histograms
  string hname = histos[0] + mode[0];
  TH2D *h2 = fDS[dataset]->getHist2(hname, false);

  for (unsigned int imode = 0; imode < mode.size(); ++imode) {
    hname = particle + "_" + mode[imode];
    TH1D *hr = new TH1D(hname.c_str(), hname.c_str(), h2->GetNbinsX(), h2->GetXaxis()->GetXbins()->GetArray());
    hr->Sumw2();
    hmode[hname] = hr;
    for (unsigned int ihist = 0; ihist < histos.size(); ++ihist) {
      hname = histos[ihist] + mode[imode];
      cout << "====> getting " << hname << " from dataset " << dataset << endl;
      TH2D *h2 = fDS[dataset]->getHist2(hname, false);
      int nbins(h2->GetNbinsX());
      cout << "x bins: " << nbins  << endl;
      TH1D *h1(0);
      //      for (int i = 1; i <= nbins; ++i) {
      for (int i = 2; i <= 3; ++i) {
	c0->cd(++ipad);
	h1 = h2->ProjectionY(Form("hist%dpt%d", ihist, i), i, i);
	h1->SetTitle(Form("hist%dpt%d %s", ihist, i, h1->GetTitle()));
	if (string::npos != hname.find("Fake310")) {
	  fitKs(h1);
	} else if (string::npos != hname.find("Fake333")) {
	  fitPhi(h1);
	} else if (string::npos != hname.find("Fake3122")) {
	  fitLambda(h1);
	}
	tl->DrawLatex(0.2, 0.7, Form("%4.2f#pm%4.2f", fYield, fYieldE));
	hr->SetBinContent(i, hr->GetBinContent(i) + fYield);
	hr->SetBinError(i, TMath::Sqrt(hr->GetBinError(i)*hr->GetBinError(i) + fYieldE*fYieldE));
      }
      c0->cd(++ipad);
      h2->Draw("colz");
    }
  }
  // c0->cd(++ipad);
  // hmode[Form("%s_nmu", particle.c_str())]->Draw("e1");
  // for (int i = 1; i <= hmode["pion_nmu"]->GetNbinsX(); ++i) {
  //   cout << hmode["pion_nmu"]->GetBinLowEdge(i) << ": " << hmode["pion_nmu"]->GetBinContent(i)
  // 	 << " +/- " << hmode["pion_nmu"]->GetBinError(i)
  // 	 << endl;
  // }

}


// ----------------------------------------------------------------------
void plotFake::fitKs(TH1D *h) {
  fIF->fVerbose = false;
  fIF->resetLimits();
  fIF->limitPar(1, 0.495, 0.503);
  fIF->limitPar(2, 0.004, 0.008);
  TF1* f1 = fIF->pol0gauss(h, 0.498, 0.006);
  f1->FixParameter(4, 0.);
  TFitResultPtr r = h->Fit(f1, "lsq", "e");
  double bwidth = h->GetBinWidth(h->FindBin(1));
  double peak = f1->GetParameter(1);
  double sigma = f1->GetParameter(2);
  double xmin = peak - 2.*sigma;
  double xmax = peak + 2.*sigma;
  double aintegral  = f1->Integral(xmin, xmax)/bwidth;
  // -- (slight?) overestimate of signal integral error by including the background
  double fintegralE = f1->IntegralError(xmin, xmax, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())/bwidth;
  // -- now set constant of pol0 to zero and integrate over signal only
  f1->SetParameter(3, 0.);
  double fintegral  = f1->Integral(xmin, xmax)/bwidth;
  if (fintegral > 0) {
    double err1 = f1->GetParError(0)/f1->GetParameter(0);
    double err2 = f1->GetParError(2)/f1->GetParameter(2);
    double errT = TMath::Sqrt(err1*err1 + err2*err2);
    fYieldE = errT*fintegral;
    fYieldE = fintegralE;
    fYield  = fintegral;
  } else {
    double fallbackXmin(0.485);
    double fallbackXmax(0.511);
    int nsgbins = h->FindBin(fallbackXmax) - h->FindBin(fallbackXmin) + 1;
    double sg = h->Integral(h->FindBin(fallbackXmin), h->FindBin(fallbackXmax));
    int nbgbins = h->FindBin(fallbackXmax) - 1 - 1 + 1;
    nbgbins    += h->GetNbinsX() - (h->FindBin(fallbackXmax) + 1) + 1;
    double bg = h->Integral(1, h->FindBin(fallbackXmin)-1) +  h->Integral(h->FindBin(fallbackXmax)+1, h->GetNbinsX());

    fYield  = sg - nsgbins*(bg/nbgbins);
    fYieldE = (sg > 1? TMath::Sqrt(sg): 1.);
  }
}



// ----------------------------------------------------------------------
void plotFake::fitPhi(TH1D *h) {
  fIF->fVerbose = true;
  fIF->resetLimits();
  fIF->limitPar(1, 1.01, 1.03);
  fIF->limitPar(2, 0.002, 0.005);
  fIF->limitPar(3, 0.001, 0.150);
  fIF->limitPar(4, 0.0051, 0.030);

  TF1 *f1 = fIF->phiKK(h);
  h->Fit(f1, "l", "e");

  TF1 *f2 = fIF->argus(h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1));
  f2->SetLineColor(kBlue);
  f2->SetLineStyle(kDashed);
  f2->SetParameter(0, f1->GetParameter(5));
  f2->SetParameter(1, f1->GetParameter(6));
  f2->SetParameter(2, f1->GetParameter(7));
  f2->Draw("same");
}


// ----------------------------------------------------------------------
void plotFake::fitLambda(TH1D *h) {

}


// ----------------------------------------------------------------------
void plotFake::loopFunction1() {

}

// ----------------------------------------------------------------------
void plotFake::bookHist(int mode) {

}



// ----------------------------------------------------------------------
void plotFake::candAnalysis() {
  fGoodCand = true;
}


// ----------------------------------------------------------------------
void plotFake::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
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
  cout << "==> plotFake::loopOverTree> loop over dataset " << fCds << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotFake::*pF)(void);
  if (ifunc == 1) pF = &plotFake::loopFunction1;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;

    candAnalysis();
    (this->*pF)();
  }

}


// ----------------------------------------------------------------------
void plotFake::setupTree(TTree *t) {
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
  t->SetBranchAddress("pv2lip",  &fb.pv2lip);
  t->SetBranchAddress("pv2lips", &fb.pv2lips);
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
void plotFake::setCuts(string cuts) {
  cout << "==> plotFake::setCuts: " << cuts << endl;

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
void plotFake::loadFiles(string afiles) {

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

      if (string::npos != stype.find("charmonium")) {
        sname = "data_charmonium";
        sdecay = "bmm";
	ds->fColor = kBlack;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      ds->fLcolor = ds->fColor;
      ds->fFcolor = ds->fColor;
      ds->fName   = sdecay;
      ds->fFullName = sname;
      fDS.insert(make_pair(sname, ds));


    } else {
      // -- MC
      pF = loadFile(sfile);
      cout << "  " << sfile << ": " << pF << endl;

      dataset *ds = new dataset();
      ds->fSize = 1;
      ds->fWidth = 2;

      string filter = "nofilter";
      if (string::npos != stype.find("etaptfilter")) filter = "etaptfilter";

      if (string::npos != stype.find("bu2jpsik")) {
        sname = "bu2jpsik_" + filter;
        sdecay = "bu2jpsik";
	ds->fColor = kBlue-7;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      cout << "  inserting as " << sname << " and " << sdecay << endl;
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
    cout << "===> " << it->first << endl;
    cout << "       " << it->second->fName << endl;
    cout << "       " << it->second->fF->GetName() << endl;
    cout << "       " << it->first << ": " << it->second->fName << ", " << it->second->fF->GetName() << endl;
  }
}
