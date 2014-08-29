#ifndef PIDTABLE
#define PIDTABLE

#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TF1.h"
#include "TH2.h"

#include "TList.h"
#include "TSortedList.h"
#include "TObjArray.h"
#include "TCollection.h"
#include "TOrdCollection.h"

#include <fstream>

#include "PidData.hh"

class PidTable: public TObject {

public:


  // -- Constructors
  // ---------------
  PidTable(int mode = 0);
  // Instantiation based on a PidTable in an ascii file
  PidTable(const char *filename, int mode = 0);

  // -- Input 
  // --------
  // from a file containing a PidTable
  void readFromFile(const char *filename = "bla.dat", int mode = 0);
  // open a rootfile containing a 2d efficiency histogram
  void readFromEffHist(TDirectory *f, const char *histname = "e6003", double fMin = -180., double fMax = 180., int mode = 0);
  // open a rootfile containing 2d histograms of 'pass' and 'total' 
  void readFromHist(TDirectory *f, const char *pass = "h6003f1", const char *total = "h1003f1", 
		    double fMin = -180., double fMax = 180., int mode = 0);
  // define grid for PidTables with array 
  void setAxes(int nP, double *pAxis, int nT, double *tAxis, int nF, double *fAxis);
  // define regular grid for PidTables with nbins, min, and max
  void setAxes(int nP, double pmin, double pmax, int nT, double tmin, double tmax, 
	       int nF, double fmin, double fmax);
  // fill an empty table (i.e. no cells defined) with input from another table
  void fillTable(PidTable &);                       
  // keep the format, fill efficiency from another table
  void fillEff(PidTable &);                         

  // -- Management
  // -------------
  // remove all cells and reset the data container
  void     clear();    
  // clear all cells and fill them with 0-0-0-0
  void     flush();    
  // concatenate adjacent (in momentum) 0-0-0-0 bins
  void     compress();
  // output into ascii file
  void dumpToFile(const char *filename = "bla.dat"); 
  void scale(double scaleFactor); 

  // -- Utilities
  void     setVerbosity(int t) {fVerbose = t;}
  void     printAll(int i = 0);
  void     print(ostream &OUT = std::cout);
  void     setPidTablesDir(const char *s = "/u/ec/ursl/epidtables");
  void     setHistName(const char *s = "e6003") {fHistName = s;}
  void     setHistMinMax(double min, double max) {fHistMin = min; fHistMax = max;}
  TH2D*    get2dHist(const char *hname, const char *title, int mode = 0);
  void     setEffAndErrMode(int mode = 0);
  void     recalculate();
  void     shiftPmax(double oldMax, double newMax); 
  void     shiftTmax(double oldMax, double newMax); 
  TIter    next() {return TIter(fDataVector);} 
  
  double   getMean(); 
  
  TString  getHistName() {return fHistName;}
  TString  getFileName() {return fFileName;}
  PidData* getData(double p, double t, double f);
  PidData* getDataRange(PidData);
  PidData* getDataRange(double pmin, double pmax, double tmin, double tmax, double fmin, double fmax);

  TH1D*    getMomentumHist(double tmin, double tmax, double fmin, double fmax, int absVal = 0);


  // -- Operations 
  // -------------
  // modify table, cell efficiency no longer 'consistent' with cell content
  void contEmptyBins(double cut = 1);
  void contLowEfficiencyBins(double cut = 0.02);

  // modify table, cell efficiency no longer 'consistent' with cell content
  void fitEmptyBins(double cut = 1);

  void compare(PidTable &);

  // combine tables, e.g. to merge different data sets (2000 + 2001 + 2002 -> 200x)
  // eff and err are recalculated from pass/tot of new table if (tot > 0)
  // if (tot == 0), then eff and err are the weighted mean and error on it
  void combine(PidTable &);

  // divide this table by another table 
  // this is useful only for eff and err
  void divide(PidTable &);

  // subtract another table from this
  // this is useful only for eff and err
  void subtract(PidTable &);

  // replace cell efficiency by smeared (gaussian with sigma = err) value
  void smear(int mode = 0);                         

  // replace cell efficiency by shifted (oldvalue +/- shift*oldvalue) value
  void shiftRel(double shift = 0.1);                


  // -- 1-d Histograms of ...
  // ------------------------
  // ... efficiency
  void effHist(TH1 *h = 0, double pmin = 3.0, double pmax = 15., double tmin = -2.4, double tmax = 2.4, 
	       double fmin = -180., double fmax = 180.);

  // ... error 
  void errHist(TH1 *h = 0, double pmin = 3.0, double pmax = 15., double tmin = -2.4, double tmax = 2.4, 
	       double fmin = -180., double fmax = 180.);

  // ... relative error
  void relErrHist(TH1 *h = 0, double pmin = 3.0, double pmax = 15., double tmin = -2.4, double tmax = 2.4, 
		  double fmin = -180., double fmax = 180.);


  // -- Graphical displays of the table vs momentum-(th)eta
  // ------------------------------------------------------
  // efficiency 
  void eff2d(TH2 *h = 0, double fmin = -180., double fmax = 180.);

  // error 
  void err2d(TH2 *h = 0, double fmin = -180., double fmax = 180.);

  // relative error 
  void relErr2d(TH2 *h = 0, double fmin = -180., double fmax = 180.);

  // the number of 'passing' candidates
  void pass2d(TH2 *h = 0, double fmin = -180., double fmax = 180.);

  // normalization
  void tot2d(TH2 *h = 0, double fmin = -180., double fmax = 180.);

  void drawEff(const char *s = "colz", double fmin = -180., double fmax = 180.);
  void drawErr(const char *s = "colz", double fmin = -180., double fmax = 180.);
  void drawTot(const char *s = "colz", double fmin = -180., double fmax = 180.);
  void drawPass(const char *s = "colz", double fmin = -180., double fmax = 180.);

  void drawTableEff(double fi = 0., const char *s = "colz");


  // -- Projections
  // --------------
  void projectP(TH1 *h = 0, double tmin=-2.4, double tmax=2.4, double fmin=-180., double fmax=180., int absVal = 0);
  void projectT(TH1 *h = 0, double pmin=3., double pmax=15., double fmin=-180., double fmax=180.);

  // -- Access to table contents (efficiency and error) indexes with p/theta/phi
  // ---------------------------------------------------------------------------
  int idD(double momentum, double theta, double phi);  
  int idR(double momentum, double theta, double phi);  

  // default (think 'eta' for CMS and 'degrees' for BABAR)
  double effD(double momentum, double theta, double phi, int sigma = 0);  
  double errD(double momentum, double theta, double phi, int sigma = 0);  
  // in radians (not meaningful for CMS)
  double effR(double momentum, double theta, double phi, int sigma = 0);
  double errR(double momentum, double theta, double phi, int sigma = 0);

  // -- eff and err integrating over some variables
  double effRange(double pmin, double pmax, double tmin, double tmax, double fmin, double fmax);
  double effPhiRange(double momentum, double theta, double fmin= -180., double fmax = 180.);
  double errPhiRange(double momentum, double theta, double fmin= -180., double fmax = 180.);
  double totPhiRange(double momentum, double theta, double fmin= -180., double fmax = 180.);
  double passPhiRange(double momentum, double theta, double fmin= -180., double fmax = 180.);

private:

  TList *fDataVector;

  TString fHistName;
  TString fFileName;
  // -- base directory for pidtables
  TString fPidTablesDir; 

  // -- maximum and minimum efficiency
  double fMaxEff, fMinEff, fMeanEff, fMeanEffErr;
  // -- maximum and minimum efficiency error
  double fMaxErr, fMinErr, fMeanErr;
  // -- SetMinimum and SetMaximum ranges
  double fHistMin, fHistMax;
  // -- validity range
  long fMinRun, fMaxRun;  
  int fVerbose, fMode;
  
  double Tmax;  
  double Tmin;
  int Tbin;
  double Pmax;
  double Pmin;
  int Pbin;
 
  TF1 *fRandom;
  TF1 *fFlat;
  
  ClassDef(PidTable,1) //Testing PidTable

};

#endif
