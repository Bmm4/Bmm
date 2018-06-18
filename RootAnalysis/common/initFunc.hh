#ifndef INITFUNC
#define INITFUNC

#include "TString.h"
#include "TObject.h"
#include "TH1.h"
#include "TF1.h"

#include <iostream>


class initFunc: public TObject {

public:

  initFunc();
  initFunc(std::string name);
  ~initFunc();

  void resetLimits();
  void limitPar(int ipar, double lo, double hi);
  void fixPar(int ipar, double fix);
  void applyLimits(TF1 *f, std::string name);
  void dumpParameters(TF1 *f);
  void dumpLimits(std::string print = "");

  // -- background functions
  TF1* pol0(TH1 *h);
  TF1* pol0(double lo, double hi);

  TF1* err(TH1 *h);
  TF1* err(double lo, double hi);

  TF1* err2(TH1 *h);
  TF1* err2(double lo, double hi);

  TF1* pol1(double lo, double hi);
  TF1* pol1(TH1 *h);
  TF1* pol1(TH1 *h, double lo, double hi);

  TF1* pol2local(TH1 *h, double width = 0.1);

  TF1* expo(double lo, double hi);
  TF1* expo(TH1 *h);
  TF1* expo(TH1 *h, double lo, double hi);


  TF1* argus(double lo, double hi);
  TF1* argus(TH1 *h);
  TF1* argus(TH1 *h, double lo, double hi);

  TF1* landau(double lo, double hi);
  TF1* landau(TH1 *h);

  TF1* bigauss(double lo, double hi);
  TF1* bigauss(TH1 *h);

  TF1* gauss(double lo, double hi);
  TF1* gauss(TH1 *h);
  TF1* Gauss(double lo, double hi);

  TF1* gauss2(double lo, double hi);
  TF1* gauss2(TH1 *h);

  TF1* gauss2c(double lo, double hi);
  TF1* gauss2c(TH1 *h);

  TF1* gauss3(double lo, double hi);
  TF1* gauss3(TH1 *h);

  TF1* pisat(double norm);
  static double iF_int_pisat(double norm);

  TF1* kstarsat(double norm);
  static double iF_int_kstarsat(double norm);

  TF1* pol1Err(double lo, double hi);
  TF1* expoErr(double lo, double hi);
  TF1* expoErr2(double lo, double hi);

  // -- signal+background functions
  TF1* pol0gauss(TH1 *h, double peak = 0.5, double sigma = 0.006);
  TF1* pol0Gauss(TH1 *h, double peak = 0.5, double sigma = 0.006);
  TF1* pol1gauss(TH1 *h, double peak = 5.3, double sigma = 0.04);
  TF1* pol1Gauss(TH1 *h, double peak = 5.3, double sigma = 0.04);

  TF1* pol1gauss2c(TH1 *h, double peak = 5.3, double sigma = 0.04);
  TF1* pol1gauss2(TH1 *h, double peak = 5.3, double sigma = 0.04, double deltaPeak = 0.1, double deltaSigma = 0.1);

  TF1* expoGauss(TH1 *h, double peak = 5.3, double sigma = 0.04);
  TF1* expogauss2(TH1 *h, double peak = 5.3, double sigma = 0.04, double deltaPeak = 0.1, double deltaSigma = 0.1);
  TF1* expogauss2c(TH1 *h, double peak = 5.3, double sigma = 0.04, double scaleSigma = 1.6);

  TF1* expoErrGauss(TH1 *h, double peak = 5.3, double sigma = 0.04, double preco = 5.14);
  TF1* expoErrgauss2(TH1 *h, double peak = 5.3, double sigma = 0.04, double peak2 = 5.3, double sigma2 = 0.1, double preco = 5.14);
  TF1* expoErrgauss2c(TH1 *h, double peak = 5.3, double sigma1 = 0.04, double sigma2 = 0.1, double preco = 5.14);
  TF1* expoErrgauss2f(TH1 *h, double peak = 5.3, double sigma1 = 0.04, double peak2 = 5.425, double sigma2 = 0.079, double fraction = -1.,
		      double preco = -1.);

  TF1* expoErr2Gauss(double lo, double hi);
  TF1* expoErr2Gauss(TH1 *h, double peak = 5.3, double sigma = 0.04, double preco = 5.14);

  TF1* pol1ErrGauss(TH1 *h, double peak = 5.3, double sigma = 0.04, double preco = 5.14);

  TF1* pol1Err2Gauss(double lo, double hi);
  TF1* pol1Err2Gauss(TH1 *h, double peak = 5.3, double sigma = 0.04, double preco = 5.14);

  TF1* pol1Err2gauss2c(double lo, double hi);

  TF1* crystalBall(TH1 *h, double peak = 5.3, double sigma = 0.04, double alpha = 1., double tailLength = 1.);
  TF1* crystalBallGauss(TH1 *h, double peak = 5.3, double sigma = 0.04, double alpha = 1., double tailLength = 1.);
  TF1* crystalBallBiGauss(TH1 *h, double peak = 5.3, double sigma = 0.04, double alpha = 1., double tailLength = 1.);
  TF1* pol1CrystalBall(TH1 *h, double peak = 5.3, double sigma = 0.04, double alpha = 1., double tailLength = 1.);
  TF1* pol1Landau(TH1 *h, double peak = 5.3, double sigma = 0.04);

  TF1* pol0BsBlind(TH1 *h);
  TF1* pol1BsBlind(TH1 *h);
  TF1* expoBsBlind(TH1 *h);

  // -- Added to include landau functions
  TF1* land(TH1 *h, double peak=5.3, double sigma=0.04);
  TF1* landsimp(TH1 *h, double peak=5.3, double sigma=0.04);
  TF1* expoErrgauss2Landau(TH1 *h, double peak = 5.3, double sigma = 0.04, double peak2 = 5.3, double sigma2 = 0.1, double preco = 5.14);
  TF1* expoErrGaussLandau(TH1 *h, double peak = 5.3, double sigma = 0.04, double preco = 5.14);

  // -- Specific models
  TF1* phiKK(TH1 *h);
  TF1* bupsik(TH1 *h, double sigma = 0.01);
  TF1* bupsik1(TH1 *h, double sigma = 0.01);
  TF1* bspsiphi(TH1 *h, double lo = 5.2, double hi = 5.8, double sigma = 0.01);

  void initPol0(double &p0, TH1 *h);
  void initPol1(double &p0, double &p1, TH1 *h);
  void initExpo(double &p0, double &p1, TH1 *h);
  void initExpoHS(double &p0, double &p1, double &p2, TH1 *h);

  double fLo, fHi;
  double fBgFractionLo, fBgFractionHi;
  double fLimitLo[20], fLimitHi[20];
  bool   fLimit[20], fFix[20];
  bool   fVerbose;

  std::string fName;

private:

  ClassDef(initFunc,1) //Testing initFunc


};

#endif
