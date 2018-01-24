#include <iostream>
#include "preselection.hh"
#include "common/util.hh"
#include "TString.h"
#include "TH1.h"

using namespace std;


#define PTMIN 5.00
#define PTMAX 9999.0

#define M1PTMIN 4.0
#define M1PTMAX 999.0

#define M2PTMIN 4.0
#define M2PTMAX 999.0

#define FL3DMAX 2.0
#define CHI2DOFMAX 5.0

#define PVIPMAX 0.1
#define PVIPSMAX 5.0

#define PVLIPMAX  1.0
#define PVLIPSMAX 5.0

#define MAXDOCAMAX 0.08
#define CLOSETRKMAX 21

#define FLSXYMIN   4.0

#define FLS3DMIN  4.0
#define FLS3DMAX 200.0

#define DOCATRKMAX 2.5

#define ISOMIN 0.0
#define ALPHAMAX 1.0

#define MASSERRORMAX 0.2

// ----------------------------------------------------------------------
presel::presel() {
  fPTMIN        = PTMIN;
  fPTMAX        = PTMAX;
  fM1PTMIN      = M1PTMIN;
  fM1PTMAX      = M1PTMAX;
  fM2PTMIN      = M2PTMIN;
  fM2PTMAX      = M2PTMAX;

  fFL3DMAX      = FL3DMAX;
  fCHI2DOFMAX   = CHI2DOFMAX;

  fPVIPMAX      = PVIPMAX;
  fPVIPSMAX     = PVIPSMAX;
  fPVLIPMAX     = PVLIPMAX;
  fPVLIPSMAX    = PVLIPSMAX;

  fMAXDOCAMAX   = MAXDOCAMAX;
  fCLOSETRKMAX  = CLOSETRKMAX;

  fFLSXYMIN     = FLSXYMIN;
  fFLS3DMIN     = FLS3DMIN;
  fFLS3DMAX     = FLS3DMAX;
  fDOCATRKMAX   = DOCATRKMAX;
  fISOMIN       = ISOMIN;
  fALPHAMAX     = ALPHAMAX;
  fMASSERRORMAX = MASSERRORMAX;

  setupMap();
}


// ----------------------------------------------------------------------
presel::presel(TH1* h1) {
  fPTMIN        = getValueByLabel(h1, "PTMIN");
  fPTMAX        = getValueByLabel(h1, "PTMAX");
  fM1PTMIN      = getValueByLabel(h1, "M1PTMIN");
  fM1PTMAX      = getValueByLabel(h1, "M1PTMAX");
  fM2PTMIN      = getValueByLabel(h1, "M2PTMIN");
  fM2PTMAX      = getValueByLabel(h1, "M2PTMAX");

  fFL3DMAX      = getValueByLabel(h1, "FL3DMAX");
  fCHI2DOFMAX   = getValueByLabel(h1, "CHI2DOFMAX");

  fPVIPMAX      = getValueByLabel(h1, "PVIPMAX");
  fPVIPSMAX     = getValueByLabel(h1, "PVIPSMAX");
  fPVLIPMAX     = getValueByLabel(h1, "PVLIPMAX");
  fPVLIPSMAX    = getValueByLabel(h1, "PVLIPSMAX");

  fMAXDOCAMAX   = getValueByLabel(h1, "MAXDOCAMAX");
  fCLOSETRKMAX  = static_cast<int>(getValueByLabel(h1, "CLOSETRKMAX"));

  fFLSXYMIN     = getValueByLabel(h1, "FLSXYMIN");
  fFLS3DMIN     = getValueByLabel(h1, "FLS3DMIN");
  fFLS3DMAX     = getValueByLabel(h1, "FLS3DMAX");
  fDOCATRKMAX   = getValueByLabel(h1, "DOCATRKMAX");
  fISOMIN       = getValueByLabel(h1, "ISOMIN");
  fALPHAMAX     = getValueByLabel(h1, "ALPHAMAX");
  fMASSERRORMAX = getValueByLabel(h1, "MASSERRORMAX");

  setupMap();
}

// ----------------------------------------------------------------------
void presel::setupMap() {
  fCuts.clear();
  fCuts.insert(make_pair("PTMIN", PTMIN));
  fCuts.insert(make_pair("PTMAX", PTMAX));
  fCuts.insert(make_pair("M1PTMIN", M1PTMIN));
  fCuts.insert(make_pair("M1PTMAX", M1PTMAX));
  fCuts.insert(make_pair("M2PTMIN", M2PTMIN));
  fCuts.insert(make_pair("M2PTMAX", M2PTMAX));

  fCuts.insert(make_pair("FL3DMAX", FL3DMAX));
  fCuts.insert(make_pair("CHI2DOFMAX", CHI2DOFMAX));

  fCuts.insert(make_pair("PVIPMAX", PVIPMAX));
  fCuts.insert(make_pair("PVIPSMAX", PVIPSMAX));
  fCuts.insert(make_pair("PVLIPMAX", PVLIPMAX));
  fCuts.insert(make_pair("PVLIPSMAX", PVLIPSMAX));

  fCuts.insert(make_pair("MAXDOCAMAX", MAXDOCAMAX));
  fCuts.insert(make_pair("CLOSETRKMAX", CLOSETRKMAX));

  fCuts.insert(make_pair("FLSXYMIN", FLSXYMIN));
  fCuts.insert(make_pair("FLS3DMIN", FLS3DMIN));
  fCuts.insert(make_pair("FLS3DMAX", FLS3DMAX));
  fCuts.insert(make_pair("DOCATRKMAX", DOCATRKMAX));
  fCuts.insert(make_pair("ISOMIN", ISOMIN));
  fCuts.insert(make_pair("ALPHAMAX", ALPHAMAX));
  fCuts.insert(make_pair("MASSERRORMAX", MASSERRORMAX));
}


// ----------------------------------------------------------------------
string presel::preselection() {
  string cut = Form(" (%3.2f<pt)&&(pt<%3.2f) && (%3.2f<m1pt)&&(m1pt<%3.2f) && (%3.2f<m2pt)&&(m2pt<%3.2f)",
			 fPTMIN, fPTMAX, fM1PTMIN, fM1PTMAX, fM2PTMIN, fM2PTMAX);
  cut += Form(" && (flsxy>%3.2f) && (fl3d<%3.2f) && (pvip<%3.2f) && !(TMath::IsNaN(pvips)) && (pvips>0) && (pvips<%3.2f)",
	      fFLSXYMIN, fFL3DMAX, fPVIPMAX, fPVIPSMAX);
  cut += Form(" && TMath::Abs(pvlip) < %3.2f && TMath::Abs(pvlips) < %3.2f", fPVLIPMAX, fPVLIPSMAX);
  cut += Form(" && (closetrk<%d) && (fls3d>%3.2f) && (fls3d<%3.2f) && (docatrk<%3.2f) && (maxdoca<%3.2f)",
	      fCLOSETRKMAX, fFLS3DMIN, fFLS3DMAX, fDOCATRKMAX, fMAXDOCAMAX);
  cut += Form(" && (chi2dof<%3.2f)  && (alpha<%3.2f) && (me<%3.2f)", fCHI2DOFMAX, fALPHAMAX, fMASSERRORMAX);
  cut += Form(" && (iso>%3.2f) && (m1iso>%3.2f) && (m2iso>%3.2f)", fISOMIN, fISOMIN, fISOMIN);
  cut += Form(" && (m1q*m2q<0)");
  return cut;
}

// ----------------------------------------------------------------------
void presel::setCut(string name, double value) {
  int verbose(0);
  if (name == "PTMIN") {
    fPTMIN = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }
  if (name == "PTMAX") {
    fPTMAX = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

  if (name == "M1PTMIN") {
    fM1PTMIN = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

  if (name == "M1PTMAX") {
    fM1PTMAX = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

  if (name == "M2PTMIN") {
    fM2PTMIN = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

  if (name == "M2PTMAX") {
    fM2PTMAX = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }


  if (name == "FL3DMAX") {
    fFL3DMAX = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

  if (name == "CHI2DOFMAX") {
    fCHI2DOFMAX = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

  if (name == "PVIPMAX") {
    fPVIPMAX = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

  if (name == "PVIPSMAX") {
    fPVIPSMAX = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

  if (name == "PVLIPMAX") {
    fPVLIPMAX = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

  if (name == "PVLIPSMAX") {
    fPVLIPSMAX = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

  if (name == "MAXDOCAMAX") {
    fMAXDOCAMAX = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

  if (name == "CLOSETRKMAX") {
    fCLOSETRKMAX = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

  if (name == "FLSXYMIN") {
    fFLSXYMIN = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

  if (name == "FLS3DMIN") {
    fFLS3DMIN = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

  if (name == "FLS3DMAX") {
    fFLS3DMAX = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

  if (name == "DOCATRKMAX") {
    fDOCATRKMAX = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

  if (name == "ISOMIN") {
    fISOMIN = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

  if (name == "ALPHAMAX") {
    fALPHAMAX = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

  if (name == "MASSERRORMAX") {
    fMASSERRORMAX = value;
    if (verbose) cout << "presel::setCut(" << name << ": " << value << endl;
  }

}

// ----------------------------------------------------------------------
bool presel::passCut(string name, double value) {
  if (name == "PTMIN")        return (value > fPTMIN);
  if (name == "PTMAX")        return (value < fPTMAX);
  if (name == "M1PTMIN")      return (value > fM1PTMIN);
  if (name == "M1PTMAX")      return (value < fM1PTMAX);
  if (name == "M2PTMIN")      return (value > fM2PTMIN);
  if (name == "M2PTMAX")      return (value < fM2PTMAX);

  if (name == "FL3DMAX")      return (value < fFL3DMAX);
  if (name == "CHI2DOFMAX")   return (value < fCHI2DOFMAX);

  if (name == "PVIPMAX")      return (value < fPVIPMAX);
  if (name == "PVIPSMAX")     return (value < fPVIPSMAX);
  if (name == "PVLIPMAX")     return (value < fPVLIPMAX);
  if (name == "PVLIPSMAX")    return (value < fPVLIPSMAX);

  if (name == "MAXDOCAMAX")   return (value < fMAXDOCAMAX);
  if (name == "CLOSETRKMAX")  return (value < fCLOSETRKMAX);

  if (name == "FLSXYMIN")     return (value > fFLSXYMIN);
  if (name == "FLS3DMIN")     return (value > fFLS3DMIN);
  if (name == "FLS3DMAX")     return (value < fFLS3DMAX);
  if (name == "DOCATRKMAX")   return (value < fDOCATRKMAX);
  if (name == "ISOMIN")       return (value > fISOMIN);
  if (name == "ALPHAMAX")     return (value < fALPHAMAX);
  if (name == "MASSERRORMAX") return (value < fMASSERRORMAX);
  return false;
}


// ----------------------------------------------------------------------
bool presel::preselection(redTreeData &b, int verbose) {
  //NO?!?!?!?!  if (!b.gmugmid) return false;

  if (b.m1q*b.m2q > 0) return false;

  if (b.pt < fPTMIN) return false;
  if (b.pt > fPTMAX) return false;
  if (verbose > 9) cout << "passed pT" << endl;

  if (b.m1pt < fM1PTMIN) return false;
  if (b.m1pt > fM1PTMAX) return false;
  if (b.m2pt < fM2PTMIN) return false;
  if (b.m2pt > fM2PTMAX) return false;
  if (verbose > 9) cout << "passed muon pT" << endl;

  if (b.flsxy < fFLSXYMIN) return false;
  if (b.fl3d > fFL3DMAX) return false;
  if (b.pvip > fPVIPMAX) return false;
  if (TMath::IsNaN(b.pvips)) return false;
  if (b.pvips < 0) return false;
  if (b.pvips > fPVIPSMAX) return false;
  if (verbose > 8) cout << "passed vertexing cuts" << endl;

  if (verbose > 8) cout << "pvlip* = " << b.pvlip << " " << b.pvlips << endl;
  if (TMath::Abs(b.pvlip) > fPVLIPMAX) return false;
  if (TMath::Abs(b.pvlips) > fPVLIPSMAX) return false;
  if (verbose > 7) cout << "passed pvlip* cuts" << endl;

  if (b.closetrk >= fCLOSETRKMAX) return false;
  if (b.fls3d < fFLS3DMIN) return false;
  if (b.fls3d > fFLS3DMAX) return false;
  if (b.docatrk > fDOCATRKMAX) return false;
  if (b.maxdoca > fMAXDOCAMAX) return false;
  if (b.me > fMASSERRORMAX) return false;
  if (verbose > 6) cout << "passed misc cuts" << endl;

  // NONONO if we want to use the onia for any x-checks!
  //   if (b.m < 4.9) return false;
  //   if (b.m > 5.9) return false;
  //   if (verbose > 4) cout << "passed mass cuts" << endl;
  if (verbose > 5) {
    cout << "chi2dof = " << b.chi2dof
	 << " iso = " << b.iso
	 << " m1iso = " << b.m1iso
	 << " m2iso = " << b.m2iso
	 << " alpha = " << b.alpha
	 << endl;
  }

  // -- physics preselection
  if (b.chi2dof > fCHI2DOFMAX) return false;
  if (b.iso < fISOMIN) return false;
  if (b.m1iso < fISOMIN) return false;
  if (b.m2iso < fISOMIN) return false;
  if (b.alpha > fALPHAMAX) return false;
  if (verbose > 0) cout << "passed physics cuts" << endl;

  return true;
}


// ----------------------------------------------------------------------
// this version is  slower than the hard-coded version! Some numbers:
// changing m2pt:  34" vs 28"
// changing alpha: 36" vs 30"
bool presel::preselection1(redTreeData &b) {
  const int verbose(-1);

  //NO?!?!?!?!  if (!b.gmugmid) return false;

  if (b.m1q*b.m2q > 0) return false;

  if (b.pt < fCuts["PTMIN"]) return false;
  if (b.pt > fCuts["PTMAX"]) return false;
  if (verbose > 9) cout << "passed pT" << endl;

  if (b.m1pt < fCuts["M1PTMIN"]) return false;
  if (b.m1pt > fCuts["M1PTMAX"]) return false;
  if (b.m2pt < fCuts["M2PTMIN"]) return false;
  if (b.m2pt > fCuts["M2PTMAX"]) return false;
  if (verbose > 9) cout << "passed muon pT" << endl;

  if (b.flsxy < fCuts["FLSXYMIN"]) return false;
  if (b.fl3d > fCuts["FL3DMAX"]) return false;
  if (b.pvip > fCuts["PVIPMAX"]) return false;
  if (TMath::IsNaN(b.pvips)) return false;
  if (b.pvips < 0) return false;
  if (b.pvips > fCuts["PVIPSMAX"]) return false;
  if (verbose > 8) cout << "passed vertexing cuts" << endl;

  if (verbose > 8) cout << "pvlip* = " << b.pvlip << " " << b.pvlips << endl;
  if (TMath::Abs(b.pvlip) > fCuts["PVLIPMAX"]) return false;
  if (TMath::Abs(b.pvlips) > fCuts["PVLIPSMAX"]) return false;
  if (verbose > 7) cout << "passed pvlip* cuts" << endl;

  if (b.closetrk >= fCuts["CLOSETRKMAX"]) return false;
  if (b.fls3d < fCuts["FLS3DMIN"]) return false;
  if (b.fls3d > fCuts["FLS3DMAX"]) return false;
  if (b.docatrk > fCuts["DOCATRKMAX"]) return false;
  if (b.maxdoca > fCuts["MAXDOCAMAX"]) return false;
  if (b.me > fCuts["MASSERRORMAX"]) return false;
  if (verbose > 6) cout << "passed misc cuts" << endl;

  // NONONO if we want to use the onia for any x-checks!
  //   if (b.m < 4.9) return false;
  //   if (b.m > 5.9) return false;
  //   if (verbose > 4) cout << "passed mass cuts" << endl;
  if (verbose > 5) {
    cout << "chi2dof = " << b.chi2dof
	 << " iso = " << b.iso
	 << " m1iso = " << b.m1iso
	 << " m2iso = " << b.m2iso
	 << " alpha = " << b.alpha
	 << endl;
  }

  // -- physics preselection
  if (b.chi2dof > fCuts["CHI2DOFMAX"]) return false;
  if (b.iso < fCuts["ISOMIN"]) return false;
  if (b.m1iso < fCuts["ISOMIN"]) return false;
  if (b.m2iso < fCuts["ISOMIN"]) return false;
  if (b.alpha > fCuts["ALPHAMAX"]) return false;
  if (verbose > 0) cout << "passed physics cuts" << endl;

  return true;
}

// ----------------------------------------------------------------------
TH1D* presel::getPreselectionNumbers() {
  TH1D *h  = new TH1D("hpresel", "", 100, 0., 100.);
  int ibin(1);
  ibin = 1; h->SetBinContent(ibin, fPTMIN); h->GetXaxis()->SetBinLabel(ibin, "cut:PTMIN");
  ibin = 2; h->SetBinContent(ibin, fPTMAX); h->GetXaxis()->SetBinLabel(ibin, "cut:PTMAX");
  ibin = 3; h->SetBinContent(ibin, fM1PTMIN); h->GetXaxis()->SetBinLabel(ibin, "cut:M1PTMIN");
  ibin = 4; h->SetBinContent(ibin, fM1PTMAX); h->GetXaxis()->SetBinLabel(ibin, "cut:M1PTMAX");
  ibin = 5; h->SetBinContent(ibin, fM2PTMIN); h->GetXaxis()->SetBinLabel(ibin, "cut:M2PTMIN");
  ibin = 6; h->SetBinContent(ibin, fM2PTMAX); h->GetXaxis()->SetBinLabel(ibin, "cut:M2PTMAX");

  ibin = 10; h->SetBinContent(ibin, fFL3DMAX); h->GetXaxis()->SetBinLabel(ibin, "cut:FL3DMAX");
  ibin = 11; h->SetBinContent(ibin, fFLS3DMIN); h->GetXaxis()->SetBinLabel(ibin, "cut:FLS3DMIN");
  ibin = 12; h->SetBinContent(ibin, fFLS3DMAX); h->GetXaxis()->SetBinLabel(ibin, "cut:FLS3DMAX");
  ibin = 13; h->SetBinContent(ibin, fFLSXYMIN); h->GetXaxis()->SetBinLabel(ibin, "cut:FLSXYMIN");
  ibin = 14; h->SetBinContent(ibin, fCHI2DOFMAX); h->GetXaxis()->SetBinLabel(ibin, "cut:CHI2DOFMAX");

  ibin = 20; h->SetBinContent(ibin, fPVIPMAX); h->GetXaxis()->SetBinLabel(ibin, "cut:PVIPMAX");
  ibin = 21; h->SetBinContent(ibin, fPVIPSMAX); h->GetXaxis()->SetBinLabel(ibin, "cut:PVIPSMAX");
  ibin = 22; h->SetBinContent(ibin, fPVLIPMAX); h->GetXaxis()->SetBinLabel(ibin, "cut:PVLIPMAX");
  ibin = 23; h->SetBinContent(ibin, fPVLIPSMAX); h->GetXaxis()->SetBinLabel(ibin, "cut:PVLIPSMAX");

  ibin = 30; h->SetBinContent(ibin, fMAXDOCAMAX); h->GetXaxis()->SetBinLabel(ibin, "cut:MAXDOCAMAX");
  ibin = 31; h->SetBinContent(ibin, fCLOSETRKMAX); h->GetXaxis()->SetBinLabel(ibin, "cut:CLOSETRKMAX");
  ibin = 32; h->SetBinContent(ibin, fDOCATRKMAX); h->GetXaxis()->SetBinLabel(ibin, "cut:DOCATRKMAX");

  ibin = 40; h->SetBinContent(ibin, fISOMIN); h->GetXaxis()->SetBinLabel(ibin, "cut:ISOMIN");
  ibin = 41; h->SetBinContent(ibin, fALPHAMAX); h->GetXaxis()->SetBinLabel(ibin, "cut:ALPHAMAX");
  ibin = 42; h->SetBinContent(ibin, MASSERRORMAX); h->GetXaxis()->SetBinLabel(ibin, "cut:MASSERRORMAX");

  return h;
}




// ----------------------------------------------------------------------
std::string preselection() {
  std::string cut = Form(" (%3.2f<pt)&&(pt<%3.2f) && (%3.2f<m1pt)&&(m1pt<%3.2f) && (%3.2f<m2pt)&&(m2pt<%3.2f)",
			 PTMIN, PTMAX, M1PTMIN, M1PTMAX, M2PTMIN, M2PTMAX);
  // cut += Form(" && (flsxy>%3.2f) && (fl3d<%3.2f) && (pvip<%3.2f) && !(TMath::IsNaN(pvips)) && (pvips>0) && (pvips<%3.2f)",
  // 	      FLSXYMIN, FL3DMAX, PVIPMAX, PVIPSMAX);
  // cut += Form(" && TMath::Abs(pvlip) < %3.2f && TMath::Abs(pvlips) < %3.2f", PVLIPMAX, PVLIPSMAX);
  // cut += Form(" && (closetrk<%d) && (fls3d>%3.2f) && (fls3d<%3.2f) && (docatrk<%3.2f) && (maxdoca<%3.2f)",
  // 			  CLOSETRKMAX, FLS3DMIN, FLS3DMAX, DOCATRKMAX, MAXDOCAMAX);
  // cut += Form(" && (chi2dof<%3.2f)  && (alpha<%3.2f) && (me<%3.2f)", CHI2DOFMAX, ALPHAMAX, MASSERRORMAX);
  // cut += Form(" && (iso>%3.2f) && (m1iso>%3.2f) && (m2iso>%3.2f)", ISOMIN, ISOMIN, ISOMIN);
  // cut += Form(" && (m1q*m2q<0)");
  return cut;
}

// ----------------------------------------------------------------------
bool preselection(redTreeData &b) {
  // const int verbose(-1);

  // //NO?!?!?!?!  if (!b.gmugmid) return false;

  // if (b.m1q*b.m2q > 0) return false;

  // if (b.pt < PTMIN) return false;
  // if (b.pt > PTMAX) return false;
  // if (verbose > 9) cout << "passed pT" << endl;

  // if (b.m1pt < M1PTMIN) return false;
  // if (b.m1pt > M1PTMAX) return false;
  // if (b.m2pt < M2PTMIN) return false;
  // if (b.m2pt > M2PTMAX) return false;
  // if (verbose > 9) cout << "passed muon pT" << endl;

  // if (b.flsxy < FLSXYMIN) return false;
  // if (b.fl3d > FL3DMAX) return false;
  // if (b.pvip > PVIPMAX) return false;
  // if (TMath::IsNaN(b.pvips)) return false;
  // if (b.pvips < 0) return false;
  // if (b.pvips > PVIPSMAX) return false;
  // if (verbose > 8) cout << "passed vertexing cuts" << endl;

  // if (verbose > 8) cout << "pvlip* = " << b.pvlip << " " << b.pvlips << endl;
  // if (TMath::Abs(b.pvlip) > PVLIPMAX) return false;
  // if (TMath::Abs(b.pvlips) > PVLIPSMAX) return false;
  // if (verbose > 7) cout << "passed pvlip* cuts" << endl;

  // if (b.closetrk >= CLOSETRKMAX) return false;
  // if (b.fls3d < FLS3DMIN) return false;
  // if (b.fls3d > FLS3DMAX) return false;
  // if (b.docatrk > DOCATRKMAX) return false;
  // if (b.maxdoca > MAXDOCAMAX) return false;
  // if (b.me > MASSERRORMAX) return false;
  // if (verbose > 6) cout << "passed misc cuts" << endl;

  // // NONONO if we want to use the onia for any x-checks!
  // //   if (b.m < 4.9) return false;
  // //   if (b.m > 5.9) return false;
  // //   if (verbose > 4) cout << "passed mass cuts" << endl;
  // if (verbose > 5) {
  //   cout << "chi2dof = " << b.chi2dof
  // 	 << " iso = " << b.iso
  // 	 << " m1iso = " << b.m1iso
  // 	 << " m2iso = " << b.m2iso
  // 	 << " alpha = " << b.alpha
  // 	 << endl;
  // }

  // // -- physics preselection
  // if (b.chi2dof > CHI2DOFMAX) return false;
  // if (b.iso < ISOMIN) return false;
  // if (b.m1iso < ISOMIN) return false;
  // if (b.m2iso < ISOMIN) return false;
  // if (b.alpha > ALPHAMAX) return false;
  // if (verbose > 0) cout << "passed physics cuts" << endl;

  return true;
}


// ----------------------------------------------------------------------
TH1D* getPreselectionNumbers() {
  TH1D *h  = new TH1D("hpresel", "", 100, 0., 100.);
  // int ibin(1);
  // ibin = 1; h->SetBinContent(ibin, PTMIN); h->GetXaxis()->SetBinLabel(ibin, "PTMIN");
  // ibin = 2; h->SetBinContent(ibin, PTMAX); h->GetXaxis()->SetBinLabel(ibin, "PTMAX");
  // ibin = 3; h->SetBinContent(ibin, M1PTMIN); h->GetXaxis()->SetBinLabel(ibin, "M1PTMIN");
  // ibin = 4; h->SetBinContent(ibin, M1PTMAX); h->GetXaxis()->SetBinLabel(ibin, "M1PTMAX");
  // ibin = 5; h->SetBinContent(ibin, M2PTMIN); h->GetXaxis()->SetBinLabel(ibin, "M2PTMIN");
  // ibin = 6; h->SetBinContent(ibin, M2PTMAX); h->GetXaxis()->SetBinLabel(ibin, "M2PTMAX");

  // ibin = 10; h->SetBinContent(ibin, FL3DMAX); h->GetXaxis()->SetBinLabel(ibin, "FL3DMAX");
  // ibin = 11; h->SetBinContent(ibin, FLS3DMIN); h->GetXaxis()->SetBinLabel(ibin, "FLS3DMIN");
  // ibin = 12; h->SetBinContent(ibin, FLS3DMAX); h->GetXaxis()->SetBinLabel(ibin, "FLS3DMAX");
  // ibin = 13; h->SetBinContent(ibin, FLSXYMIN); h->GetXaxis()->SetBinLabel(ibin, "FLSXYMIN");
  // ibin = 14; h->SetBinContent(ibin, CHI2DOFMAX); h->GetXaxis()->SetBinLabel(ibin, "CHI2DOFMAX");

  // ibin = 20; h->SetBinContent(ibin, PVIPMAX); h->GetXaxis()->SetBinLabel(ibin, "PVIPMAX");
  // ibin = 21; h->SetBinContent(ibin, PVIPSMAX); h->GetXaxis()->SetBinLabel(ibin, "PVIPSMAX");
  // ibin = 22; h->SetBinContent(ibin, PVLIPMAX); h->GetXaxis()->SetBinLabel(ibin, "PVLIPMAX");
  // ibin = 23; h->SetBinContent(ibin, PVLIPSMAX); h->GetXaxis()->SetBinLabel(ibin, "PVLIPSMAX");

  // ibin = 30; h->SetBinContent(ibin, MAXDOCAMAX); h->GetXaxis()->SetBinLabel(ibin, "MAXDOCAMAX");
  // ibin = 31; h->SetBinContent(ibin, CLOSETRKMAX); h->GetXaxis()->SetBinLabel(ibin, "CLOSETRKMAX");
  // ibin = 32; h->SetBinContent(ibin, DOCATRKMAX); h->GetXaxis()->SetBinLabel(ibin, "DOCATRKMAX");

  // ibin = 40; h->SetBinContent(ibin, ISOMIN); h->GetXaxis()->SetBinLabel(ibin, "ISOMIN");
  // ibin = 41; h->SetBinContent(ibin, ALPHAMAX); h->GetXaxis()->SetBinLabel(ibin, "ALPHAMAX");
  // ibin = 42; h->SetBinContent(ibin, MASSERRORMAX); h->GetXaxis()->SetBinLabel(ibin, "MASSERRORMAX");


  return h;
}

// ----------------------------------------------------------------------
void printRedTreeEvt(redTreeData &b) {
  std::cout << b.run << "/" << b.evt << " mu&hlt: " << b.gmuid << "/" << b.hlt1
	    << std::endl << "  "
	    << " B: " << " " << b.pt << "/" << b.eta << " "
	    << " mu: " << " " << b.m1q << "/" << b.m2q << " "
	    << b.m1pt << "/" << b.m2pt << ", "  << b.m1eta << "/" << b.m2eta
	    << " ->c" << ((TMath::Abs(b.m1eta) < 1.4 && TMath::Abs(b.m2eta) < 1.4)?0 :1)
	    << " fl: " << b.flsxy << " " << b.fls3d << " (" << b.fl3d
	    << ") ip: " << b.pvip << " " << b.pvips << " " << b.pvlip << " " << b.pvlips << " "
	    << std::endl << "  "
	    << " iso: " << b.iso << " " << b.closetrk  << " " << b.docatrk << " "
	    << " miso: " << b.m1iso << " " << b.m2iso
	    << " vtx: dca " << b.maxdoca << " chi " << b.chi2dof << " a " << b.alpha << " me " << b.me
	    << std::endl << "  "
	    << " pvw8: " << b.pvw8 << " gt: " << b.gtqual
	    << std::endl << "   "
	    << preselection()
	    << std::endl;
}
