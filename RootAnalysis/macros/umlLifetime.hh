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

#define  MLO  4.8
#define  MHI  6.0
#define  TLO  0.
#define  THI  15.



class model {
public:
  model() {};
  ~model() {};
  std::string name;

  RooRealVar *m, *t;   // variables: mass and reconstructed decay time

  // -- fit (fixed) parameters:
  RooRealVar *sgMassPeak, *sgMassSigma, *bgMassSlope;
  RooRealVar *sgTau, *bgTau;
  RooRealVar *sgN, *bgN;

  RooTruthModel  *tTruth;
  RooDecay       *sgPdfT, *bgPdfT;
  RooGaussian    *sgPdfM;
  RooExponential *bgPdfM;

  RooAbsPdf      *sgPdf, *bgPdf;
  RooAbsPdf      *modelPdf;
  RooArgSet      poi;
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
  model* createModel1(std::string name, int mode = 0);

  void   runToy(model *pM, int nsig, int nbg);

  // -- code for loops
  void   loopFunction1();

  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);

protected:

  // ----------------------------------------------------------------------
  ClassDef(umlLifetime,1)

};


#endif
