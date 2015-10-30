#ifndef UTIL_H
#define UTIL_H


#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TMinuit.h"
#include "TDirectory.h"
#include "TColor.h"
#include "TLatex.h"

#include <vector>


class dataset;

// extern TROOT  *gROOT;
// extern TStyle *gStyle;
// extern TFile  *gFile;
// extern TDirectory  *gDirectory;

// ======================================================================
// General utilities
// ======================================================================
void printNonZero(TH1*); 
void stampAndSave(TCanvas *fC, const char *s = "bla.eps");
void shrinkPad(double b = 0.1, double l = 0.1, double r = 0.1, double t = 0.1);
void babar(double x, double y, double scale = 1.0, int prel = 0);
void zone(int x = 1, int y = 1, TCanvas *c = 0);
int  wait();
void colors(int choice = 0);
void replaceAll(std::string &s, const std::string &from, const std::string &to);
std::vector<int> defVector(int n, ...);

std::string formatTex(double n, std::string name, int digits, int sgn = 0);
std::string formatTex(double n, std::string name, std::string tag);
void stamp(double x1, std::string text1, double x2, std::string text2, TLatex *tl = 0);

void rmSubString(std::string &sinput, const std::string &remove);
void rmPath(std::string &sInput);


std::vector<std::string> glob(std::string basedir, std::string basename);

// ======================================================================
// Utitilities for displays and typing shortcuts
// ====================================================================== 
// -- Histograms
void setTitles(TH1 *h, const char *sx, const char *sy, 
	       float size = 0.05, float xoff = 1.1, float yoff = 1.1, float lsize = 0.05, int font = 42);
void setHist(TH1 *h, int color = kBlack, int symbol = 20, double size = 1., double width = 2.);
void setHist(TH1 *h, dataset *ds);
void setGraph(TGraph *g, int color = kBlack, int symbol = 20, double size = 1., double width = 2.);
void setFilledHist(TH1 *h, int lcol = kBlack, int fcol = kYellow, int fstyle = 1000, int width = 1);
void setMaximum(double scale = 1.2, TH1 *h1 = 0, TH1 *h2 = 0);
void showOverflow(TH1 *h);

TH1D *unmix(TH1D *rightSign, TH1D *wrongSign, double chid = 0.181); 

// ======================================================================
// Utitilities for calculations
// ====================================================================== 
// -- Chi2 probability
// ----------------------------------------------------------------------
double chi2Prob(double chi2, double ndof); 

// -- Chi^2 test for two histograms, where the errors are sqrt(n_bin)
//    If they are constrained to the same area, constrain = 0
//    If they are NOT constrained, set constrain = -1 
double chi2Test(TH1*, TH1*, double& chi2, double& ndof, int constrain = 0);

// -- Chi^2 test for two histograms, where the errors are taken from the histogram
//    If they are constrained to the same area, constrain = 0
//    If they are NOT constrained, set constrain = -1 
//    This function is equivalent to chi2Test() if the errors are sqrt(n_bin)
double chi2TestErr(TH1*, TH1*, double& chi2, double& ndof, int constrain = 0);

void average(double &av, double &error, int n, double *val, double *err); 

// TGraph* scanParameter(int parnum, int nsig, TMinuit &a, void (*func)(int &, double *, double &, double *, int)); 

// ======================================================================
// Utitilities for error calculations
// ====================================================================== 
// -- d'Agostini page 86
double dEff(int, int);
double dEff(int, double);
double dEff(double, int);
double dEff(double, double, double, double);
// -- from the ROOT source
double dBinomial(double, double);
// -- normal error propagation, assumes error = sqrt(..)
double dRatio(double, double);
// -- normal error propagation:  n, error, N, error
double dRatio(double, double, double, double);
// -- error on branching fraction n/(epsilon N)
double dBF(double, double, double, double, double, double);
// -- Add up in quadrature the errors of all bins
double getError(TH1*);
double getErrorRange(TH1*, int = -1, int = -1);

// ======================================================================
// -- gen utilities
// ======================================================================
bool isBeautyMeson(int i); 
bool isBeautyMesonWeak(int i); 
bool isBeautyBaryon(int i); 
bool isCharmMeson(int i); 
bool isCharmMesonWeak(int i); 
bool isCharmBaryon(int i); 
bool isLightMeson(int i); 
bool isStableCharged(int i); 



#endif
