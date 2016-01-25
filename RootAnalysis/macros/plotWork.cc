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
void plotWork::makeAll(int bitmask) {

  if (bitmask & 0x1) {
    prodValidation("all", "default_nofilter", "inelastic_nofilter", "pythia6_nofilter"); 
    prodValidation("all", "default_etaptfilter", "inelastic_etaptfilter", "pythia6_etaptfilter"); 
    prodValidation("all", "default_etaptfilter", "inelastic_etaptfilter", "hardqcd8_etaptfilter");
    prodValidation("all", "default_nofilter", "pythia6_nofilter", "noevtgen_nofilter");
    prodValidation("all", "default_etaptfilter", "pythia6_etaptfilter", "noevtgen_etaptfilter");
  } else if (bitmask & 0x2) {
    prodValidation("all", "default_etaptfilter", "nophotos_etaptfilter", "nophotonnoevtgen_etaptfilter");
    prodValidation("tau", "default_etaptfilter", "dgs_etaptfilter",      "dgs0_etaptfilter");
  }

}



// ----------------------------------------------------------------------
void plotWork::prodValidation(string hist, string ds1, string ds2, string ds3, bool loga, bool legend, double xleg, double yleg) {
  cout << endl << "==>prodValidation: " << hist << " ds: " << ds1 << " vs. " << ds2 << " vs. " << ds3 << endl;  

  if ("all" == hist) {
    zone(2,2, c0); 
    
    vector<int> particles = defVector(1, 531); 

    for (unsigned int i = 0; i< particles.size(); ++i) {
      cout << "loop: " << particles[i] << endl;
      c0->cd(1);
      prodValidation(Form("pt%d", particles[i]), ds1, ds2, ds3, true, true, 0.4, 0.75);     
      c0->cd(2);
      prodValidation(Form("cpt%d", particles[i]), ds1, ds2, ds3, false, true, 0.4, 0.75); 
      c0->cd(3);
      prodValidation(Form("eta%d", particles[i]), ds1, ds2, ds3, false); 
      c0->cd(4);
      prodValidation(Form("iso%d", particles[i]), ds1, ds2, ds3, false); 
      c0->SaveAs(Form("%s/prodValidation-%d-%s-%s-%s.pdf", fDirectory.c_str(), particles[i], ds1.c_str(), ds2.c_str(), ds3.c_str())); 
    }
    return;
  } else if ("tau" == hist) {
    zone(2,2, c0); 
    vector<int> particles = defVector(1, 531); 
    for (unsigned int i = 0; i< particles.size(); ++i) {
      cout << "loop: " << particles[i] << endl;
      c0->cd(1);
      prodValidation(Form("t%d", particles[i]), ds1, ds2, ds3, true, true, 0.4, 0.75);     
      c0->cd(2);
      prodValidation(Form("tpos%d", particles[i]), ds1, ds2, ds3, true, true, 0.4, 0.75); 
      c0->cd(3);
      prodValidation(Form("tneg%d", particles[i]), ds1, ds2, ds3, true, true, 0.4, 0.75); 
      c0->SaveAs(Form("%s/tauValidation-%d-%s-%s-%s.pdf", fDirectory.c_str(), particles[i], ds1.c_str(), ds2.c_str(), ds3.c_str())); 
    }
    return;
  }

  overlay(Form("%s", hist.c_str()), ds1, 
	  Form("%s", hist.c_str()), ds2, 
	  (ds3 != ""? Form("%s", hist.c_str()) : ""), ds3, 
	  UNITY, loga, legend, xleg, yleg); 

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

      if (string::npos != stype.find("bmm")) {
        sname = "data_bmm"; 
        sdecay = "bmm"; 
	ds->fColor = kBlack; 
	ds->fSymbol = 24; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365; 
      }

      if (string::npos != stype.find("bu2jpsik")) {
        sname = "data_bu2jpsik"; 
        sdecay = "bu2jpsik"; 
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

      // -------------------------
      // -- genAnalyis files below
      // -------------------------
      if (string::npos != stype.find("default,")) {
        sname = "default_" + filter; 
        sdecay = "default"; 
	ds->fColor = kBlue-7; 
	ds->fSymbol = 24; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365; 
      } 

      if (string::npos != stype.find("inelastic,")) {
        sname = "inelastic_" + filter; 
        sdecay = "inelastic"; 
	ds->fColor = kMagenta-1; 
	ds->fSymbol = 25; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365; 
      } 

      if (string::npos != stype.find("hardqcd8,")) {
        sname = "hardqcd8_" + filter; 
        sdecay = "hardqcd8"; 
	ds->fColor = kRed-2; 
	ds->fSymbol = 30; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365; 
      } 

      if (string::npos != stype.find("msel1,")) {
        sname = "pythia6_" + filter; 
        sdecay = "pythia6"; 
	ds->fColor = kGreen+3; 
	ds->fSymbol = 26; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365; 
      } 

      if (string::npos != stype.find("noevtgen,")) {
        sname = "noevtgen_" + filter; 
        sdecay = "noevtgen"; 
	ds->fColor = kYellow+2; 
	ds->fSymbol = 28; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365; 
      } 

      if (string::npos != stype.find("dgs0,")) {
        sname = "dgs0_" + filter; 
        sdecay = "dgs0"; 
	ds->fColor = kGreen+3; 
	ds->fSymbol = 28; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3356; 
      } 
      
      if (string::npos != stype.find("dgs,")) {
        sname = "dgs_" + filter; 
        sdecay = "dgs"; 
	ds->fColor = kRed+2; 
	ds->fSymbol = 28; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3356; 
      } 

      if (string::npos != stype.find("nophotos,")) {
        sname = "nophotos_" + filter; 
        sdecay = "nophotos"; 
	ds->fColor = kGreen+3; 
	ds->fSymbol = 28; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3356; 
      } 

      if (string::npos != stype.find("nophotonnoevtgen,")) {
        sname = "nophotonnoevtgen_" + filter; 
        sdecay = "nophotonnoevtgen"; 
	ds->fColor = kRed+2; 
	ds->fSymbol = 28; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3356; 
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

