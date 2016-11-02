#ifndef UMLLIFETIME_h
#define UMLLIFETIME_h

#include "plotClass.hh"
#include "common/dataset.hh"
#include "common/AnalysisDistribution.hh"
#include "redTreeData.hh"

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooLognormal.h"
#include "RooPolynomial.h"
#include "RooTruthModel.h"
#include "RooDecay.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"

#define  MLO  4.8
#define  MHI  6.0
#define  TLO  0.
#define  THI  15.
#define  TAU0 1.68
#define  NCHAN 2

// ----------------------------------------------------------------------
class model1 {
public:
  model1(std::string name = "") : fName(name) {}
  ~model1() {
    std::cout << "deleting model1 " << fName << std::endl;
    delete bsN;
    delete bdN;
    delete bgN;

    delete bsMassPeak;
    delete bsMassSigma;
    delete bdMassPeak;
    delete bdMassSigma;
    delete bgMassSlope;

    delete bsTau;
    delete bdTau;
    delete bgTau;

    delete tTruth;
    delete bsPdfT;
    delete bdPdfT;
    delete bgPdfT;

    delete bsPdfM;
    delete bdPdfM;
    delete bgPdfM;

    delete bsPdf;
    delete bdPdf;
    delete bgPdf;

    delete modelPdf;
  };

  std::string fName;

  // -- fit (fixed) parameters:
  RooRealVar *bsTau, *bdTau, *bgTau;
  RooRealVar *bsN, *bdN, *bgN;
  RooRealVar *bsMassPeak, *bsMassSigma, *bdMassPeak, *bdMassSigma, *bgMassSlope;

  RooTruthModel  *tTruth;
  RooDecay       *bsPdfT, *bdPdfT, *bgPdfT;
  RooGaussian    *bsPdfM, *bdPdfM;
  RooExponential *bgPdfM;

  RooAbsPdf      *bsPdf, *bdPdf, *bgPdf;
  RooAbsPdf      *modelPdf;
};


// ----------------------------------------------------------------------
class model2 {
public:
  model2(std::string name = "") : fName(name) {}
  ~model2() {
    std::cout << "deleting model1 " << fName << std::endl;
    // delete bsN;
    // delete bdN;
    // delete bgN;

    // delete bsMassPeak;
    // delete bsMassSigma;
    // delete bdMassPeak;
    // delete bdMassSigma;
    // delete bgMassSlope;

    // delete bsTau;
    // delete bdTau;
    // delete bgTau;

    // delete tTruth;
    // delete bsPdfT;
    // delete bdPdfT;
    // delete bgPdfT;

    // delete bsPdfM;
    // delete bdPdfM;
    // delete bgPdfM;

    // delete bsPdf;
    // delete bdPdf;
    // delete bgPdf;

    // delete modelPdf;
  };

  std::string fName;

  // -- fit (fixed) parameters:
  RooRealVar *bsTau, *bdTau, *bgTau;
  RooRealVar *bsN[NCHAN], *bdN[NCHAN], *bgN[NCHAN];
  RooRealVar *bsMassPeak[NCHAN], *bsMassSigma[NCHAN], *bdMassPeak[NCHAN], *bdMassSigma[NCHAN], *bgMassSlope[NCHAN];

  RooTruthModel  *tTruth[NCHAN];
  RooDecay       *bsPdfT[NCHAN], *bdPdfT[NCHAN], *bgPdfT[NCHAN];
  RooGaussian    *bsPdfM[NCHAN], *bdPdfM[NCHAN];
  RooExponential *bgPdfM[NCHAN];

  RooAbsPdf      *bsPdf[NCHAN], *bdPdf[NCHAN], *bgPdf[NCHAN];
  RooAbsPdf      *modelPdf[NCHAN];

  RooSimultaneous *simPdf;
  RooCategory     *sample;

};



// ----------------------------------------------------------------------
class umlLifetime: public plotClass {

public :
  umlLifetime(std::string dir = "results",
	      std::string files = "plotResults.2016.files",
	      std::string cuts = "baseCuts.cuts",
	      std::string setup = "");
  virtual ~umlLifetime();

  // -- Main analysis methods
  virtual void makeAll(std::string what = "all");
  virtual void init();
  virtual void loadFiles(std::string afiles);
  virtual void bookHist(std::string dsname);

  // -- models
  model1*  createModel1(std::string name, int mode = 0);
  model2*  createModel2(std::string name, int mode = 0);

  RooDataSet *createData2(model2 *m, int nsig, int nbg, bool channelWise = false);

  void   runToys1(std::string toy , int ntoys, int nsig, int nbg);
  void   runToys2(std::string toy , int ntoys, int nsig, int nbg);

  // -- code for loops
  void   loopFunction1();

  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);

  int fRndmSeed;

protected:

  RooRealVar *fm, *ft;   // variables: mass and reconstructed decay time

  double fMeanValue, fMeanValueWidth;
  double fMeanError, fMeanErrorWidth;

  // ----------------------------------------------------------------------
  ClassDef(umlLifetime,1)

};


#endif
