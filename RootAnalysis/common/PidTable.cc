#include "PidTable.hh"

#include "TF1.h"
#include "TProfile.h"
#include "TMath.h"

#include <algorithm>
#include <set>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

ClassImp(PidTable)

// ----------------------------------------------------------------------
PidTable::PidTable(int mode) {
  fPidTablesDir = TString("/u/ec/ursl/pidtables");
  fMinRun = 0;
  fMaxRun = 1000000;
  fMinEff = 99.;
  fMaxEff = -1.;
  fHistMin = 0.;
  fHistMax = 1.1;
  fMinErr = 99.;
  fMaxErr = -1.;
  fVerbose = 1;
  fDataVector = new TList();
  fMode = mode;

 Tmax = 2.4;
 Tmin = -2.4;
 Tbin = 48;
 Pmax = 20.;
 Pmin = 0.;
 Pbin = 40;


  fRandom = new TF1("f1", "gaus", -10., 10.);
  fRandom->SetParameters(1., 0., 1.);
  fFlat = new TF1("fPIDTable", "pol0", 0., 1.);
  fFlat->SetParameter(0, 1.);
}


// ----------------------------------------------------------------------
PidTable::PidTable(const char *s, int mode) {
  fPidTablesDir = TString("/u/ec/ursl/pidtables");
  fMinRun = 0;
  fMaxRun = 1000000;
  fVerbose = 1;
  fMinEff = 99.;
  fMaxEff = -1.;
  fHistMin = 0.;
  fHistMax = 1.1;
  fMinErr = 99.;
  fMaxErr = -1.;
  fDataVector = new TList();
  fMode = mode;

  fRandom = new TF1("f1", "gaus", -10., 10.);
  fRandom->SetParameters(1., 0., 1.);
  fFlat = new TF1("fPIDTable", "pol0", 0., 1.);
  fFlat->SetParameter(0, 1.);

  readFromFile(s, mode);
}


// ----------------------------------------------------------------------
void PidTable::readFromEffHist(TDirectory *f, const char *histName, double phiMin, double phiMax, int mode) {
  setHistName(histName);

  TH2D *h     = (TH2D*)f->Get(histName);
  if (!h) {
    if (fVerbose > -1) cout << "Error: Could not retrieve histogram " << histName << endl;
  }
  else {
    if (fVerbose > 0) cout << "Reading from " << histName << endl;
    for (Int_t ix = 1; ix <= h->GetNbinsX(); ++ix) {
      for (Int_t iy = 1; iy <= h->GetNbinsY(); ++iy) {
	PidData *datum = new PidData(h->GetYaxis()->GetBinLowEdge(iy), h->GetYaxis()->GetBinLowEdge(iy+1),
				     h->GetXaxis()->GetBinLowEdge(ix), h->GetXaxis()->GetBinLowEdge(ix+1),
				     phiMin, phiMax,
				     h->GetBinContent(ix, iy), h->GetBinError(ix, iy),
				     -99., -99., -99., -99., mode); // ?? should this maybe be initialized to zero?
        fDataVector->Add(datum);
      }
    }
  }
}


// ----------------------------------------------------------------------
void PidTable::readFromHist(TDirectory *f, const char *passName, const char *totName, double phiMin, double phiMax, int mode) {
  setHistName(passName);

  TH2D *hPass = (TH2D*)f->Get(passName);
  TH2D *hTot  = (TH2D*)f->Get(totName);
  if (!hPass) {
    if (fVerbose > -1) cout << "Error: Could not retrieve histogram " << passName << endl;
  }
  if (!hTot) {
    if (fVerbose > -1) cout << "Error: Could not retrieve histogram " << totName << endl;
  }
  if (hTot && hPass) {
    if (fVerbose > 0) cout << "Reading from " << passName  << " and " << totName << endl;
    double n, N, nE, NE;

    for (Int_t ix = 1; ix <= hPass->GetNbinsX(); ++ix) {
      for (Int_t iy = 1; iy <= hPass->GetNbinsY(); ++iy) {

	n = hPass->GetBinContent(ix, iy);
	N = hTot->GetBinContent(ix, iy);

	nE = hPass->GetBinError(ix, iy);
	NE = hTot->GetBinError(ix, iy);

	PidData *datum = new PidData(hPass->GetYaxis()->GetBinLowEdge(iy), hPass->GetYaxis()->GetBinLowEdge(iy+1),
				     hPass->GetXaxis()->GetBinLowEdge(ix), hPass->GetXaxis()->GetBinLowEdge(ix+1),
				     phiMin, phiMax,
				     -99., -99., n, nE, N, NE, mode);
        fDataVector->Add(datum);
	if (datum->getE() < fMinEff) fMinEff = datum->getE();
	if (datum->getE() > fMaxEff) fMaxEff = datum->getE();
	if (datum->getS() < fMinErr) fMinErr = datum->getS();
	if (datum->getS() > fMaxErr) fMaxErr = datum->getS();
      }
    }
  }
}


// ----------------------------------------------------------------------
void PidTable::readFromFile(const char *filename, int mode) {

  TString fullName;
  if (filename[0] == '/') {
    fullName = filename;
  }
  else {
    fullName = filename;
    //    fullName = fPidTablesDir + "/" + filename;
  }

  fFileName = fullName;

  float pmin, pmax, tmin, tmax, fmin, fmax, e, se;
  float pass, total, passE, totalE;
  if (fVerbose > 0) cout << "Read table contents from " << fullName << endl;

  char  buffer[200];
  ifstream is(fullName);
  int readLines(0);
  TString tcom;
  string scom;
  while (is.getline(buffer, 200, '\n')) {
    scom = buffer;
    if (string::npos != scom.find("# Comment")) {
      tcom = buffer;
      tcom.ReplaceAll("# Comment:            ", "");
      fComment = tcom;
      continue;
    }
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%f %f %f %f %f %f %f %f %f %f %f %f",
	   &pmin, &pmax, &tmin, &tmax, &fmin, &fmax, &e, &se, &pass, &passE, &total, &totalE);
    PidData *datum = new PidData(pmin, pmax, tmin, tmax, fmin, fmax, e, se, pass, passE, total, totalE, mode);
    fDataVector->Add(datum);
    if (e < fMinEff)  fMinEff = e;
    if (e > fMaxEff)  fMaxEff = e;
    if (se < fMinErr) fMinErr = se;
    if (se > fMaxErr) fMaxErr = se;
    readLines++;
  }
  if (readLines < 1) {
    if (fVerbose > -1) cout << "  Nothing read from " << fullName << endl;
  }
  is.close();
}


// ----------------------------------------------------------------------
void PidTable::printAll(int i) {
  TIter next(fDataVector); PidData *iTable;
  if (i > 0) {print();}
  while ((iTable = (PidData*)next())) {
    cout << *iTable << endl;
  }

}

// ----------------------------------------------------------------------
double PidTable::getMean() {
  int cells(0);
  double pass(0.), total(0.);
  double a(0.), wi(0.), w(0.), meanE(0.), meanS(0.), meanU(0.);
  double tmin(9999.), tmax(-9999.), pmin(9999.), pmax(-9999.);

  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    pass  += iTable->getPass();
    total += iTable->getTot();
    cells++;
    if (iTable->getPmin() < pmin) pmin = iTable->getPmin();
    if (iTable->getPmax() > pmax) pmax = iTable->getPmax();
    if (iTable->getTmin() < tmin) tmin = iTable->getTmin();
    if (iTable->getTmax() > tmax) tmax = iTable->getTmax();

    a = iTable->getS();
    if (a > 1e-8) {
      wi = 1./(a*a);
      w     += wi;
      meanE += wi*iTable->getE();
      meanS += wi*iTable->getS();
    }
  }

  if (w) fMeanEff = meanE / w;
  if (w) fMeanEffErr  = 1./TMath::Sqrt(w);  // uncertainty of mean efficiency
  if (w) fMeanErr = meanS / w;              // mean of efficiency uncertainties

  return fMeanEff;
}


// ----------------------------------------------------------------------
void PidTable::scale(double scaleFactor) {

  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    iTable->setPass(scaleFactor*(iTable->getPass()));
    iTable->setTot(iTable->getTot());
    iTable->calcEffAndErr();
  }
}


// ----------------------------------------------------------------------
void PidTable::print(ostream &OUT) {

  int cells(0);
  double pass(0.), total(0.);
  double a(0.), wi(0.), w(0.), meanE(0.), meanS(0.), meanU(0.);
  double tmin(9999.), tmax(-9999.), pmin(9999.), pmax(-9999.);

  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    pass  += iTable->getPass();
    total += iTable->getTot();
    cells++;
    if (iTable->getPmin() < pmin) pmin = iTable->getPmin();
    if (iTable->getPmax() > pmax) pmax = iTable->getPmax();
    if (iTable->getTmin() < tmin) tmin = iTable->getTmin();
    if (iTable->getTmax() > tmax) tmax = iTable->getTmax();

    a = iTable->getS();
    if (a > 1e-8) {
      wi = 1./(a*a);
      w     += wi;
      meanE += wi*iTable->getE();
      meanS += wi*iTable->getS();
    }
  }

  if (w) meanE /= w;
  if (w) meanU  = 1./TMath::Sqrt(w);  // uncertainty of mean efficiency
  if (w) meanS /= w;                  // mean of efficiency uncertainties

  char line[200];
  sprintf(line, "# ============"); for (int i = 0; i<getFileName().Length(); ++i) sprintf(line, "%s=", line);
  OUT << line << endl;
  OUT << "# Summary for " << getFileName() << endl;
  sprintf(line, "# ------------");  for (int i = 0; i<getFileName().Length(); ++i) sprintf(line, "%s-", line);
  OUT << line << endl;

  OUT << "# Comment:            " << fComment << endl;
  OUT << "# Number of cells:    " << cells << endl;
  OUT << "# Number of pass :    " << Form("%9.0f", pass)  << endl;
  OUT << "# Number of total:    " << Form("%9.0f", total) << endl;
  OUT << "# -------------------------------" << endl;

  OUT << "# min momentum:       " << pmin << endl;
  OUT << "# max momentum:       " << pmax << endl;
  OUT << "# min polar angle:    " << tmin << endl;
  OUT << "# max polar angle:    " << tmax << endl;
  OUT << "# -------------------------------" << endl;

  OUT << "# min efficiency:       " << fMinEff << endl;
  OUT << "# max efficiency:       " << fMaxEff << endl;
  OUT << "# min error:            " << fMinErr << endl;
  OUT << "# max error:            " << fMaxErr << endl;
  OUT << "# -------------------------------" << endl;

  OUT << "# <Efficiency>:       " << meanE << endl;
  OUT << "#                 +/- " << meanU << endl;
  OUT << "# <Delta Efficiency>: " << meanS << endl;
  OUT << "# -------------------------------" << endl;

}


// ----------------------------------------------------------------------
void PidTable::setEffAndErrMode(int mode) {
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    iTable->setEffAndErrMode(mode);
  }
}

// ----------------------------------------------------------------------
void PidTable::recalculate() {
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    iTable->calcEffAndErr();
  }
}

// ----------------------------------------------------------------------
void PidTable::setAxes(int pA, double *pAxis, int tA, double *tAxis, int fA, double *fAxis) {

  clear();

  for (Int_t ip = 0; ip < pA; ++ip) {
    for (Int_t it = 0; it < tA; ++it) {
      for (Int_t ifi = 0; ifi < fA; ++ifi) {
	PidData *datum = new PidData(pAxis[ip], pAxis[ip+1], tAxis[it], tAxis[it+1], fAxis[ifi],fAxis[ifi+1],
				     0., 0., 0., 0., 0., 0.);
	fDataVector->Add(datum);
      }
    }
  }
}


// ----------------------------------------------------------------------
void PidTable::setAxes(int nP, double pmin, double pmax, int nT, double tmin, double tmax,
		       int nF, double fmin, double fmax) {
  double pbin = (pmax-pmin)/nP;
  double tbin = (tmax-tmin)/nT;
  double fbin = (fmax-fmin)/nF;

  clear();

  for (Int_t ip = 0; ip < nP; ++ip) {
    for (Int_t it = 0; it < nT; ++it) {
      for (Int_t ifi = 0; ifi < nF; ++ifi) {
	PidData *datum = new PidData(pmin + ip*pbin, pmin + (ip+1)*pbin,
				     tmin + it*tbin, tmin + (it+1)*tbin,
				     fmin +ifi*fbin, fmin +(ifi+1)*fbin,
				     0., 0., 0., 0., 0., 0.);
	fDataVector->Add(datum);
      }
    }
  }
}


// ----------------------------------------------------------------------
void PidTable::fillTable(PidTable &b) {
  PidData *a;
  TIter bnext = b.next(); PidData *iTable;
  while ((iTable = (PidData*)bnext())) {
    if (fVerbose > 1)  cout << "fillFrom: getting  " << *iTable << endl;
    a = new PidData(*iTable);
    a->setPass(iTable->getPass());
    a->setTot(iTable->getTot());
    a->calcEffAndErr();
    fDataVector->Add(a);
  }
}


// ----------------------------------------------------------------------
void PidTable::fillEff(PidTable &b) {
  PidData *a;
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    if (fVerbose > 1)  cout << "fillEff for  " << *iTable << endl;
    a = b.getDataRange(*iTable);
    iTable->setPass(a->getPass());
    iTable->setTot(a->getTot());
    iTable->calcEffAndErr();
    delete a;
  }
}


// ----------------------------------------------------------------------
void PidTable::smear(int mode) {
  double bla(0.);
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    if (mode == 0) {
      // do nothing
    }
    else if (mode == 99) {
      bla = iTable->getE() + fRandom->GetRandom()*iTable->getS();
      if (bla < 0.) {
	bla = 0.;
	cout << "Reset to 0 in cell (" << iTable->getPmin() << "," << iTable->getTmin() << "," << iTable->getFmin()<<")"<<endl;
      } else if (bla > 1.) {
	bla = 1.;
	cout << "Reset to 1 in cell (" << iTable->getPmin() << "," << iTable->getTmin() << "," << iTable->getFmin()<<")"<<endl;
      }
      iTable->setE(bla);
    } else {
      bla = iTable->getE() + mode*iTable->getS();
      if (bla < 0.) bla = 0.;
      iTable->setE(bla);
    }
  }
}


// ----------------------------------------------------------------------
void PidTable::shiftRel(double shift) {
  double bla(0.);
  if (TMath::Abs(shift) > 0.0001) {
    TIter next(fDataVector); PidData *iTable;
    while ((iTable = (PidData*)next())) {
      bla = iTable->getE() + shift*iTable->getE();
      if (bla < 0.) {
	bla = 0.;
	cout << "Reset to 0 in cell (" << iTable->getPmin() << "," << iTable->getTmin() << "," << iTable->getFmin()<<")"<<endl;
      } else if (bla > 1.) {
	bla = 1.;
	cout << "Reset to 1 in cell (" << iTable->getPmin() << "," << iTable->getTmin() << "," << iTable->getFmin()<<")"<<endl;
      }
      iTable->setE(bla);
    }
  }  else {
    cout << fFileName << " no shifting, because shift = " << shift << endl;
  }
}


// ----------------------------------------------------------------------
int PidTable::idD(double p, double t, double f) {
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    if (iTable->isCell(p, t, f)) {
      double bla = fFlat->GetRandom();
      if (bla < iTable->getE()) {
	return 1;
      } else {
	return 0;
      }
    }
  }
  return 0;
}

// ----------------------------------------------------------------------
// ?? double check that phi works in CMS??
int PidTable::idR(double p, double t, double f) {
  return (idD(p, t*57.2958, f*57.2958/2));
}


// ----------------------------------------------------------------------
double PidTable::effD(double p, double t, double f, int sigma) {
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    if (iTable->isCell(p, t, f)) {
      return iTable->getE();
    }
  }
  return -99.;
}


// ----------------------------------------------------------------------
double PidTable::effR(double momentum, double theta, double phi, int sigma) {
  return effD(momentum, theta*57.2958, phi*57.2958, sigma);
}


// ----------------------------------------------------------------------
double PidTable::errD(double p, double t, double f, int run) {
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    if (iTable->isCell(p, t, f)) return iTable->getS();
  }
  return -99.;
}

// ----------------------------------------------------------------------
double PidTable::effRange(double pmin,double pmax,double tmin,double tmax,double fmin,double fmax) {
  PidData a;
  a.setEffAndErrMode(fMode);
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    if ((pmin < iTable->getPmax()) && (iTable->getPmin() < pmax)
	&& (tmin < iTable->getTmax()) && (iTable->getTmin() < tmax)
	&& (fmin < iTable->getFmax()) && (iTable->getFmin() < fmax)
	) {
      a.merge(iTable);
      if (fVerbose > 1)  cout << "effRange: merging  " << *iTable << endl;
    }
  }
  return a.getE();
}


// ----------------------------------------------------------------------
double PidTable::effPhiRange(double momentum, double theta, double fmin, double fmax) {
  PidData a;
  a.setEffAndErrMode(fMode);
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    if ((fmin < iTable->getFmax()) && (iTable->getFmin() < fmax)
	&& (iTable->getPmin() <= momentum) && (momentum < iTable->getPmax())
	&& (iTable->getTmin() <= theta)    && (theta < iTable->getTmax())
	) {
      a.merge(iTable);
      if (fVerbose > 1)  cout << "effPhiRange: merging  " << *iTable << endl << "into " << a << endl;
    }
  }
  return a.getE();
}


// ----------------------------------------------------------------------
double PidTable::passPhiRange(double momentum, double theta, double fmin, double fmax) {
  PidData a;
  a.setEffAndErrMode(fMode);
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    if ((fmin < iTable->getFmax()) && (iTable->getFmin() < fmax)
	&& (iTable->getPmin() <= momentum) && (momentum < iTable->getPmax())
	&& (iTable->getTmin() <= theta)    && (theta < iTable->getTmax())
	) {
      a.merge(iTable);
      if (fVerbose > 1)  cout << "passPhiRange: merging  " << *iTable << endl;
    }
  }
  return a.getPass();
}


// ----------------------------------------------------------------------
double PidTable::totPhiRange(double momentum, double theta, double fmin, double fmax) {
  PidData a;
  a.setEffAndErrMode(fMode);
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    if ((fmin < iTable->getFmax()) && (iTable->getFmin() < fmax)
	&& (iTable->getPmin() <= momentum) && (momentum < iTable->getPmax())
	&& (iTable->getTmin() <= theta)    && (theta < iTable->getTmax())
	) {
      a.merge(iTable);
      if (fVerbose > 1)  cout << "totPhiRange: merging  " << *iTable << endl;
    }
  }
  return a.getTot();
}


// ----------------------------------------------------------------------
double PidTable::errPhiRange(double momentum, double theta, double fmin, double fmax) {
  PidData a;
  a.setEffAndErrMode(fMode);
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    if ((fmin < iTable->getFmax()) && (iTable->getFmin() < fmax)
	&& (iTable->getPmin() <= momentum) && (momentum < iTable->getPmax())
	&& (iTable->getTmin() <= theta)    && (theta < iTable->getTmax())
	) {
      a.merge(iTable);
      if (fVerbose > 1)  cout << "errPhiRange: merging  " << *iTable << endl;
    }
  }
  return a.getS();
}


// ----------------------------------------------------------------------
double PidTable::errR(double momentum, double theta, double phi, int run) {
  return errD(momentum, theta*57.2958, phi*57.2958, run);
}


// ----------------------------------------------------------------------
void PidTable::dumpToFile(const char *psname) {
  if (fVerbose > 0) cout << "Dumping table to file " << psname << endl;
  fFileName = psname;
  ofstream EFF(psname);
  print(EFF);
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    EFF << *iTable << endl;
  }
  EFF.close();
}


// ----------------------------------------------------------------------
// ?? missing implementation ??
void PidTable::fitEmptyBins(double cut) {
  if (fVerbose > 0) cout << "continue table" << endl;
  double oldTheta(-10.);
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    if (iTable->getTctr() > oldTheta) {
      oldTheta = iTable->getTmax();
      // do what needs to be done
    }
    cout << "do something with oldTheta: " << oldTheta << endl;
  }
}

// ----------------------------------------------------------------------
void PidTable::contLowEfficiencyBins(double cut) {
  if (fVerbose > 0) cout << "continue table with contLowEfficiencyBins, cut: " << cut << endl;
  double oldEff(0.), oldErr(0.);
  TIter next(fDataVector); PidData *iTable;
  PidData *begin = (PidData*)fDataVector->First();
  double plo = begin->getPmax();
  cout << "maximum momentum of first cell: " << plo << endl;
  while ((iTable = (PidData*)next())) {
    if ((iTable->getE() < cut) && (iTable != begin)) {
      cout << "hello: " << iTable->getE() << " . " << oldEff << endl;
      iTable->print();
      cout << "hello: " << iTable->getE() << " . " << oldEff << endl;
      iTable->setE(oldEff);
      iTable->setS(oldErr);
    }
    // -- update oldEff only with good statistics values
    if (iTable->getPctr() > plo) {
      if ((iTable->getS() < 0.1)) {
	oldEff = iTable->getE();
	oldErr = iTable->getS();
      }
    } else {
      // -- when looking at first momentum cell, reset to zero.
      oldEff = 0.;
      oldErr = 0.;
    }
  }
}


// ----------------------------------------------------------------------
void PidTable::contEmptyBins(double cut) {
  if (fVerbose > 0) cout << "continue table" << endl;
  double oldEff(0.), oldErr(0.);
  TIter next(fDataVector); PidData *iTable;
  PidData *begin = (PidData*)fDataVector->First();
  while ((iTable = (PidData*)next())) {
    if ((iTable->getPass() < 1.) || (iTable == begin)) {
      iTable->setE(oldEff);
      iTable->setS(oldErr);
    }
    // -- update oldEff only with good statistics values
    if (iTable != begin) {
      PidData *previous = iTable - 1;
      if ((iTable->getPass() >= cut) || (previous->getPmin() < 0.01)) {
	oldEff = iTable->getE();
	oldErr = iTable->getS();
      }
      // -- when looking at last momentum cell, reset to zero.
      if (iTable->getPmax() > 24.9) {
	oldEff = 0.;
	oldErr = 0.;
      }
    }
  }
}

// ----------------------------------------------------------------------
void PidTable::compare(PidTable &cmp) {
  double eCmp;
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    eCmp = cmp.effD(iTable->getPmax(), iTable->getTmax(), iTable->getFmax());
    if (TMath::Abs(eCmp - iTable->getE()) > 0.01) {
      if (fVerbose > 1) cout << "Disagreement in " << *iTable << endl;
    }
  }
}


// ----------------------------------------------------------------------
void PidTable::setPidTablesDir(const char*s) {
  fPidTablesDir = s;
}


// ----------------------------------------------------------------------
PidData* PidTable::getData(double p, double t, double f) {
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    if (iTable->isCell(p, t, f)) return iTable;
  }
  if (fVerbose > -1) cout << "did not find PidData for given values of p, t, f." << endl;
  return 0;
}


// ----------------------------------------------------------------------
PidData* PidTable::getDataRange(PidData d) {
  PidData *a = new PidData();
  double eps(1e-6);
  a->setEffAndErrMode(fMode);
  TIter next(fDataVector); PidData *iTable;
  if (fVerbose > 1)  cout << "d: " << d.getPmin() << "/" << d.getPmax() << " .. " << d.getTmin() << "/" << d.getTmax() << endl;
  while ((iTable = (PidData*)next())) {
    if ((d.getPmin() < iTable->getPmax()-eps) && (iTable->getPmin()+eps < d.getPmax())
	&& (d.getTmin() < iTable->getTmax()-eps) && (iTable->getTmin()+eps < d.getTmax())
	&& (d.getFmin() < iTable->getFmax()-eps) && (iTable->getFmin()+eps < d.getFmax())
	) {
      a->merge(iTable);
      if (fVerbose > 1) cout << "getDataRange: merging  " << *iTable << endl;
    }
  }
  return a;
}


// ----------------------------------------------------------------------
PidData* PidTable::getDataRange(double pmin, double pmax, double tmin, double tmax, double fmin, double fmax) {
  double eps(1e-6);
  PidData *a = new PidData(pmin, pmax, tmin, tmax, fmin, fmax);
  a->setEffAndErrMode(fMode);
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    if ((pmin < iTable->getPmax()-eps) && (iTable->getPmin()+eps < pmax)
	&& (tmin < iTable->getTmax()-eps) && (iTable->getTmin()+eps < tmax)
	&& (fmin < iTable->getFmax()-eps) && (iTable->getFmin()+eps < fmax)
	) {
      a->merge(iTable);
      if (fVerbose > 4)  cout << "getDataRange: merging  " << *iTable << endl;
    }
  }
  if (fVerbose > 1)  cout << "getDataRange: returning  " << *a << endl;
  return a;
}


// ----------------------------------------------------------------------
// ?? Why is this called this way??
TH1D* PidTable::getMomentumHist(double tmin, double tmax, double fmin, double fmax, int absVal) {
  double pmin, pmax;
  char name[100];
  sprintf(name, "theta: %5.2f .. %5.2f, phi = %5.2f .. %5.2f", tmin, tmax, fmin, fmax);
  TH1D *h = new TH1D("hPidTable", name, Pbin, Pmin, Pmax);
  PidData *a;
  for (int ibins = 1; ibins <= Pbin; ++ibins) {
    pmin = (ibins-1)*(Pmax-Pmin)/Pbin;
    pmax = (ibins)*(Pmax-Pmin)/Pbin;
    a = getDataRange(pmin, pmax, tmin, tmax, fmin, fmax);
    cout << *a << endl;
    cout << " p = " << pmin << "  " << pmax
	 << " (" << (int)(a->getPctr()/Pmax*Pbin) << ")"
	 << " (" << a->getPctr() << ")"
	 << " t = " << tmin << "  " << tmax
	 << " f = " << fmin << "  " << fmax
	 << " e = " << a->getE() << " +/- " << a->getS() << endl;
    h->SetBinContent((int)(a->getPctr()/Pmax*Pbin), a->getE());
    h->SetBinError((int)(a->getPctr()/Pmax*Pbin), a->getS());
    delete a;
  }
  return h;
}


// ----------------------------------------------------------------------
void PidTable::combine(PidTable &b) {
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    PidData *pData = b.getData(iTable->getPctr(), iTable->getTctr(), iTable->getFctr());
    if (pData) {
      if (pData->getPass() > -1) {
	iTable->combine(pData->getPass(), pData->getTot());
      } else {
	iTable->merge(pData);
      }
      if (iTable->getE() < fMinEff) fMinEff = iTable->getE();
      if (iTable->getE() > fMaxEff) fMaxEff = iTable->getE();
      if (iTable->getS() < fMinErr) fMinErr = iTable->getS();
      if (iTable->getS() > fMaxErr) fMaxErr = iTable->getS();
    }
    else {
      if (fVerbose > 1)
	cout << "b does not contain PidData for p = " << (iTable->getPmin() + iTable->getPmax())/2.
	     << ", t = " << (iTable->getTmin() + iTable->getTmax())/2.
	     << ", f = " << (iTable->getFmin() + iTable->getFmax())/2.
	     << endl;
    }
  }
}


// ----------------------------------------------------------------------
void PidTable::divide(PidTable &b) {
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    PidData *pData = b.getData(iTable->getPctr(), iTable->getTctr(), iTable->getFctr());
    if (pData) {
      iTable->divide(pData);
      if (iTable->getE() < fMinEff) fMinEff = iTable->getE();
      if (iTable->getE() > fMaxEff) fMaxEff = iTable->getE();
      if (iTable->getS() < fMinErr) fMinErr = iTable->getS();
      if (iTable->getS() > fMaxErr) fMaxErr = iTable->getS();
    }
    else {
      if (fVerbose > 1)
	cout << "other Table does not contain PidData for p = " << (iTable->getPmin() + iTable->getPmax())/2.
	     << ", t = " << (iTable->getTmin() + iTable->getTmax())/2.
	     << ", f = " << (iTable->getFmin() + iTable->getFmax())/2.
	     << endl;
    }
  }
  fFileName += "/";
  fFileName += b.getFileName();
}

// ----------------------------------------------------------------------
void PidTable::subtract(PidTable &b) {
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    PidData *pData = b.getData(iTable->getPctr(), iTable->getTctr(), iTable->getFctr());
    if (pData) {
      iTable->subtract(pData);
      if (iTable->getE() < fMinEff) fMinEff = iTable->getE();
      if (iTable->getE() > fMaxEff) fMaxEff = iTable->getE();
      if (iTable->getS() < fMinErr) fMinErr = iTable->getS();
      if (iTable->getS() > fMaxErr) fMaxErr = iTable->getS();
    }
    else {
      if (fVerbose > 1)
	cout << "other Table does not contain PidData for p = " << (iTable->getPmin() + iTable->getPmax())/2.
	     << ", t = " << (iTable->getTmin() + iTable->getTmax())/2.
	     << ", f = " << (iTable->getFmin() + iTable->getFmax())/2.
	     << endl;
    }
  }
}


// ----------------------------------------------------------------------
void PidTable::clear() {
  if (fVerbose > 0) cout << "Data size before clearance: " << fDataVector->LastIndex() << endl;
  fDataVector->Delete();
  if (fVerbose > 0) cout << "Data size after clearance: " << fDataVector->LastIndex() << endl;
}


// ----------------------------------------------------------------------
void PidTable::flush() {
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    iTable->setPass(0.);
    iTable->setPassE(0.);
    iTable->setTot(0.);
    iTable->setTotE(0.);
    iTable->setE(0.);
    iTable->setS(0.);
  }
}


// ----------------------------------------------------------------------
void PidTable::compress() {
  TList *newVector = new TList();
  PidData *n, *nn, *t;
  TIter next(fDataVector); PidData *iTable;
  PidData *end = (PidData*)fDataVector->Last();
  while ((iTable = (PidData*)next())) {
    if (iTable->isZero()) {
      double theta = iTable->getTmax();
      double phi   = iTable->getFmax();
      n = iTable + 1;
      if (fVerbose > 1) cout << "Zero: " << *iTable << endl
			     << "Next: " << *n << endl;

      while (n->isZero()
	     && (n->getTctr() < theta)
	     && (n->getFctr() < phi)
	     && (iTable != end)
	     ) {
	if (fVerbose > 2) cout << "   next one too: " << *n << endl;
	n += 1;
      }
      nn = n - 1;
      t = new PidData(iTable->getPmin(), nn->getPmax(),
		      iTable->getTmin(), nn->getTmax(),
		      iTable->getFmin(), nn->getFmax(),
		      0., 0., 0., 0., 0., 0.);
      if (iTable != end) {
	iTable = nn;
      }
      else {
	newVector->Add(t);
	break;
      }
    }
    else {
      t = new PidData(*iTable);
    }
    newVector->Add(t);
  }
  int bla = fVerbose;
  fVerbose = 0;
  clear();
  fVerbose = bla;
  fDataVector = newVector;
  cout << "Here is probably a memory leak" << endl;
}


// ----------------------------------------------------------------------
void PidTable::effHist(TH1 *h, double pmin, double pmax, double tmin, double tmax,
		       double fmin, double fmax) {
  int drawIt(0);
  if (0 == h) {
    drawIt = 1;
    double min(0.), max(1.);
    if (fMinEff > 0.5) min = 0.5;
    if (fMaxEff < 0.2) max = 0.2;
    if (fMaxEff < 0.1) max = 0.1;
    if (fMaxEff < 0.02) max = 0.02;
    if (fMaxEff > 1.0) min = 1.1;
    if (fMaxEff > 2.0) min = 2.1;
    if (fMaxEff > 3.0) min = 3.1;
    h = new TH1D("hefftmp", "", 100, min, max);
  }
  h->Reset();
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    if ((pmin < iTable->getPmax()) && (iTable->getPmin() < pmax)
	&& (tmin < iTable->getTmax()) && (iTable->getTmin() < tmax)
	&& (fmin < iTable->getFmax()) && (iTable->getFmin() < fmax)
	) {
      if (iTable->getS() > 0.000001) h->Fill(iTable->getE());
    }
  }
  if (drawIt) h->Draw();
}


// ----------------------------------------------------------------------
void PidTable::errHist(TH1 *h, double pmin, double pmax, double tmin, double tmax,
		       double fmin, double fmax) {
  int drawIt(0);
  if (0 == h) {
    drawIt = 1;
    double min(0.), max(1.);
    if (fMinErr > 0.5)  min = 0.5;
    if (fMaxErr < 0.5)  max = 0.5;
    if (fMaxErr < 0.3)  max = 0.3;
    if (fMaxErr < 0.2)  max = 0.2;
    if (fMaxErr < 0.1)  max = 0.1;
    if (fMaxErr < 0.02) max = 0.02;
    h = new TH1D("herrtmp", "", 100, min, max);
  }
  h->Reset();
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    if ((pmin < iTable->getPmax()) && (iTable->getPmin() < pmax)
	&& (tmin < iTable->getTmax()) && (iTable->getTmin() < tmax)
	&& (fmin < iTable->getFmax()) && (iTable->getFmin() < fmax)
	) {
      if (iTable->getE() > 0.000001) h->Fill(iTable->getS());
    }
  }
  if (drawIt) h->Draw();
}


// ----------------------------------------------------------------------
void PidTable::relErrHist(TH1 *h, double pmin, double pmax, double tmin, double tmax,
			  double fmin, double fmax) {
  int drawIt(0);
  if (0 == h) {
    drawIt = 1;
    double min(0.), max(1.);
    if (fMinErr > 0.5) min = 0.5;
    if (fMaxErr < 0.3) max = 0.3;
    if (fMaxErr < 0.2) max = 0.2;
    if (fMaxErr < 0.1) max = 0.1;
    if (fMaxErr < 0.02) max = 0.02;
    h = new TH1D("hrelerrtmp", "", 100, min, max);
  }
  h->Reset();
  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    if ((pmin < iTable->getPmax()) && (iTable->getPmin() < pmax)
	&& (tmin < iTable->getTmax()) && (iTable->getTmin() < tmax)
	&& (fmin < iTable->getFmax()) && (iTable->getFmin() < fmax)
	) {
      if (iTable->getE() > 0.000001) h->Fill(iTable->getS()/iTable->getE());
    }
  }
  if (drawIt) h->Draw();
}


// ----------------------------------------------------------------------
void PidTable::projectP(TH1 *h, double tmin, double tmax, double fmin, double fmax, int absVal) {

  int drawIt(0);

  if (0 == h) {
    drawIt = 1;
    double min(0.), max(1.);
    if (fMinEff > 0.5) min = 0.5;
    if (fMaxEff < 0.2) max = 0.2;
    if (fMaxEff < 0.1) max = 0.1;
    if (fMaxEff < 0.02) max = 0.02;
    h = new TH1D("hpprojtmp", "", Pbin, Pmin, Pmax);  h->Sumw2();
    h->SetMinimum(min);  h->SetMaximum(max);
  }

  PidData *a, *b;
  for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
    a = getDataRange(h->GetBinLowEdge(ix), h->GetBinLowEdge(ix+1),
		     tmin, tmax, fmin, fmax);
    if (absVal == 1) {
      b = getDataRange(h->GetBinLowEdge(ix), h->GetBinLowEdge(ix+1),
		       -tmax, -tmin, fmin, fmax);
      a->merge(b);
      delete b;
    }
    h->SetBinContent(ix, a->getEff());
    h->SetBinError(ix, a->getErr());
    delete a;
  }

  if (drawIt) {
    h->Draw();
  }

}


// ----------------------------------------------------------------------
void PidTable::projectT(TH1 *h, double pmin, double pmax, double fmin, double fmax) {
  int drawIt(0);
  if (0 == h) {
    drawIt = 1;
    double min(0.), max(1.);
    if (fMinEff > 0.5) min = 0.5;
    if (fMaxEff < 0.2) max = 0.2;
    if (fMaxEff < 0.1) max = 0.1;
    if (fMaxEff < 0.02) max = 0.02;
    h = new TH1D("htprojtmp", "", Tbin, Tmin, Tmax);
    h->SetMinimum(min);
    h->SetMaximum(max);
  }

  PidData *a;
  for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
    a = getDataRange(pmin, pmax,
		     h->GetBinLowEdge(ix), h->GetBinLowEdge(ix+1),
		     fmin, fmax);
    h->SetBinContent(ix, a->getEff());
    h->SetBinError(ix, a->getErr());
    delete a;
  }

  if (drawIt) h->Draw();
}


// ----------------------------------------------------------------------
void PidTable::eff2d(TH2 *h, double fmin, double fmax) {
  double e;
  h->Reset();
  for (Int_t ix = 1; ix <= h->GetNbinsX(); ++ix) {
    for (Int_t iy = 1; iy <= h->GetNbinsY(); ++iy) {
      e = effPhiRange(0.5*(h->GetYaxis()->GetBinLowEdge(iy) + h->GetYaxis()->GetBinLowEdge(iy+1)),
		      0.5*(h->GetXaxis()->GetBinLowEdge(ix) + h->GetXaxis()->GetBinLowEdge(ix+1)),
		      fmin, fmax);
      h->SetBinContent(ix, iy, e);
    }
  }
}


// ----------------------------------------------------------------------
void PidTable::err2d(TH2 *h, double fmin, double fmax) {
  double e;
  for (Int_t ix = 1; ix <= h->GetNbinsX(); ++ix) {
    for (Int_t iy = 1; iy <= h->GetNbinsY(); ++iy) {
      e = errPhiRange(h->GetYaxis()->GetBinLowEdge(iy), h->GetXaxis()->GetBinLowEdge(ix), fmin, fmax);
      h->SetBinContent(ix, iy, e);
    }
  }
}


// ----------------------------------------------------------------------
void PidTable::relErr2d(TH2 *h, double fmin, double fmax) {
  double e, s;
  for (Int_t ix = 1; ix <= h->GetNbinsX(); ++ix) {
    for (Int_t iy = 1; iy <= h->GetNbinsY(); ++iy) {
      e = effPhiRange(h->GetYaxis()->GetBinLowEdge(iy), h->GetXaxis()->GetBinLowEdge(ix), fmin, fmax);
      s = errPhiRange(h->GetYaxis()->GetBinLowEdge(iy), h->GetXaxis()->GetBinLowEdge(ix), fmin, fmax);
      if (e > 1e-6) {
	h->SetBinContent(ix, iy, s/e);
      }
    }
  }
}


// ----------------------------------------------------------------------
void PidTable::pass2d(TH2 *h, double fmin, double fmax) {
  double e;
  h->Reset();
  for (Int_t ix = 1; ix <= h->GetNbinsX(); ++ix) {
    for (Int_t iy = 1; iy <= h->GetNbinsY(); ++iy) {
      e = passPhiRange(h->GetYaxis()->GetBinLowEdge(iy), h->GetXaxis()->GetBinLowEdge(ix), fmin, fmax);
      h->SetCellContent(ix, iy, e);
    }
  }
}


// ----------------------------------------------------------------------
void PidTable::tot2d(TH2 *h, double fmin, double fmax) {
  double e;
  h->Reset();
  for (Int_t ix = 1; ix <= h->GetNbinsX(); ++ix) {
    for (Int_t iy = 1; iy <= h->GetNbinsY(); ++iy) {
      e = totPhiRange(h->GetYaxis()->GetBinLowEdge(iy), h->GetXaxis()->GetBinLowEdge(ix), fmin, fmax);
      h->SetCellContent(ix, iy, e);
    }
  }
}


// ----------------------------------------------------------------------
TH2D* PidTable::get2dHist(const char *hname, const char *title, int mode) {
  if (0 == mode) {
    TH2D *h = new TH2D(hname, title, Tbin, Tmin, Tmax, Pbin, Pmin, Pmax);
    h->SetMinimum(fHistMin);
    h->SetMaximum(fHistMax);
    return h;
  } else {

    double xBins[50], yBins[50];
    set<double> etaBins, ptBins;
    TIter next(fDataVector); PidData *iTable;
    while ((iTable = (PidData*)next())) {
      etaBins.insert(iTable->getTmin());
      etaBins.insert(iTable->getTmax());

      ptBins.insert(iTable->getPmin());
      ptBins.insert(iTable->getPmax());
    }

    set<double>::iterator it;
    int i(0);
    for (it = etaBins.begin(); it != etaBins.end(); it++) {
      xBins[i] = *it;
      ++i;
    }

    i = 0;
    for (it = ptBins.begin(); it != ptBins.end(); it++) {
      yBins[i] = *it;
      ++i;
    }

    for (unsigned int j = 0; j < etaBins.size(); ++j) {
      cout << j << " " << xBins[j] << endl;
    }

    for (unsigned int j = 0; j < ptBins.size(); ++j) {
      cout << j << " " << yBins[j] << endl;
    }

    TH2D *h = new TH2D(hname, title, etaBins.size()-1, xBins, ptBins.size()-1, yBins);
    return h;
  }
}


// ----------------------------------------------------------------------
void PidTable::drawEff(const char *s, double fmin, double fmax) {
  char name[100];
  sprintf(name, "%s", fHistName.Data());
  TH2D *h = new TH2D("PidTableH2", name, Tbin, Tmin, Tmax, Pbin, Pmin, Pmax);
  eff2d(h, fmin, fmax);
  h->SetMinimum(fHistMin);
  h->SetMaximum(fHistMax);
  h->Draw(s);
}


// ----------------------------------------------------------------------
// ?? What is this?? This is wrong!!
void PidTable::drawTableEff(double fi, const char *s) {
  char name[100];
  sprintf(name, "Table Eff: phi = %5.3f", fi);
  TH2D *h = new TH2D("PidTableH2", name, Tbin, Tmin, Tmax, Pbin, Pmin, Pmax);
  double bin_size_T = (Tmax + fabs(Tmin))/Tmin;
  double bin_size_P = (Pmax + fabs(Pmin))/Pmax;
  double p, t;
  for (int ix = -Tbin/2 + 1; ix <= Tbin/2; ++ix) {
    t = bin_size_T*ix - bin_size_T/2.;
    for (int iy = 1; iy <= Pbin; ++iy) {
      p = bin_size_P*iy - bin_size_P/2.;
      if (fVerbose > 2) cout << "t = " << t << " p = " << p << " e = " << effD(p, t, fi) << endl;
      if (effD(p, t, fi) > -1.) h->SetCellContent(ix, iy, effD(p, t, fi));
    }
  }
  h->Draw(s);
}


// ----------------------------------------------------------------------
void PidTable::drawErr(const char *s, double fmin, double fmax) {
  char name[100];
  sprintf(name, "Err: phi = %5.3f ... %5.3f", fmin, fmax);
  TH2D *h = new TH2D("PidTableH2", name, Tbin, Tmin, Tmax, Pbin, Pmin, Pmax);
  err2d(h, fmin, fmax);
  h->Draw(s);
}


// ----------------------------------------------------------------------
void PidTable::drawTot(const char *s, double fmin, double fmax) {
  char name[100];
  sprintf(name, "Total: phi = %5.3f ... %5.3f", fmin, fmax);
  TH2D *h = new TH2D("PidTableH2", name, Tbin, Tmin, Tmax, Pbin, Pmin, Pmax);
  tot2d(h, fmin, fmax);
  h->Draw(s);
}


// ----------------------------------------------------------------------
void PidTable::drawPass(const char *s, double fmin, double fmax) {
  char name[100];
  sprintf(name, "Pass: phi = %5.3f ... %5.3f", fmin, fmax);
  TH2D *h = new TH2D("PidTableH2", name, Tbin, Tmin, Tmax, Pbin, Pmin, Pmax);
  pass2d(h, fmin, fmax);
  h->Draw(s);
}

// ----------------------------------------------------------------------
void PidTable::shiftPmax(double oldMax, double newMax) {

  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    if (iTable->getPmax() > oldMax) iTable->setPmax(newMax);
  }

}

// ----------------------------------------------------------------------
void PidTable::shiftTmax(double oldMax, double newMax) {

  TIter next(fDataVector); PidData *iTable;
  while ((iTable = (PidData*)next())) {
    if (iTable->getTmax() > oldMax) iTable->setTmax(newMax);
  }

}
