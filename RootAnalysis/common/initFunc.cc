#include "initFunc.hh"

#include "TROOT.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"

#include <iostream>
#include <iomanip>


#include "HFMasses.hh"

ClassImp(initFunc)

using namespace std;

namespace {

  //-------------------------------------------------------------------------
  // This is the usual Landau from root
  double iF_landau(double *x, double *par) {
    //     const                        mpv   sigma
    return par[2] * TMath::Landau(x[0],par[0],par[1]);
  }

  // ----------------------------------------------------------------------
  // This is a simplfied "landau" using an analytical formula.
  // I use it because it was easiter for me to understand the integral (area) of this formula,
  double iF_landausimp(double *x, double *par) {  // simplified landau
    double xx = ( x[0]- par[0] ) / par[1];
    return par[2] * 0.3989 * TMath::Exp(-0.5 * ( xx + TMath::Exp(-xx) ) );
  }

  // ----------------------------------------------------------------------
  double iF_expo(double *x, double *par) {
    return par[0]*TMath::Exp(x[0]*par[1]);
  }

  // ----------------------------------------------------------------------
  double iF_expo_HS(double *x, double *par) {
    if (x[0] > par[2]) {
      return par[0]*TMath::Exp((x[0]-par[2])*par[1]);
    } else {
      return 0.;
    }
  }

  // ----------------------------------------------------------------------
  double iF_err(double *x, double *par) {
    // from DK: a4*(TMath::Erf((a1-x)/a2))+a3)
    return par[3]*(TMath::Erf((par[0]-x[0])/par[1])+par[2]);
  }

  // ----------------------------------------------------------------------
  double iF_err2(double *x, double *par) {
    // par[0]: step
    // par[1]: resolution
    // par[2]: level of plateau
    return par[2]*(TMath::Erf((par[0]-x[0])/par[1]) + 1.);
  }

  // ----------------------------------------------------------------------
  double iF_cb(double *x, double *par) {
    // par[0]:  mean
    // par[1]:  sigma
    // par[2]:  alpha, crossover point
    // par[3]:  n, length of tail
    // par[4]:  N, normalization

    Double_t cb = 0.0;
    Double_t exponent = 0.0;

    if (x[0] > par[0] - par[2]*par[1]) {
      exponent = (x[0] - par[0])/par[1];
      cb = TMath::Exp(-exponent*exponent/2.);
    } else {
    double nenner  = TMath::Power(par[3]/par[2], par[3])*TMath::Exp(-par[2]*par[2]/2.);
    double zaehler = (par[0] - x[0])/par[1] + par[3]/par[2] - par[2];
    zaehler = TMath::Power(zaehler, par[3]);
    cb = nenner/zaehler;
    }

    if (par[4] > 0.) {
      cb *= par[4];
    }

    return cb;
  }

  // ----------------------------------------------------------------------
  double iF_argus(double *x, double *par) {
    // par[0] -> normalization
    // par[1] -> exponential factor
    // par[2] -> endpoint
    //           > 0: normal situation (zero above endpoint)
    //           < 0: inverted situation (zero below endpoint)

    //  double ebeam = 10.58/2;
    double ebeam = par[2];
    double ebeam2 = ebeam*ebeam;
    double background = 0.;
    double x2 = x[0]*x[0];
    double ratio = x2/ebeam2;
    if (par[2] < 0) ratio = ebeam2/x2;
    if (ratio < 1.) {
      background = par[0]*x[0] * sqrt(1. - (ratio)) * exp(par[1] * (1. - (ratio)));
    } else {
      background = 0.;
    }
    return background;
  }

  // ----------------------------------------------------------------------
  double iF_argus_gauss2(double *x, double *par) {
    // par[0] -> Gaussian const
    // par[1] -> Gaussian mean
    // par[2] -> Gaussian sigma
    // par[3] -> fraction in second gaussian
    // par[4] -> sigma of second gaussian
    // par[5] -> Argus normalization
    // par[6] -> Argus exponential factor
    // par[7] -> Argus endpoint (fixed!)
    //           > 0: normal situation (zero above endpoint)
    //           < 0: inverted situation (zero below endpoint)

    //  double ebeam = 10.58/2;
    double ebeam = par[7];
    double ebeam2 = ebeam*ebeam;
    double background = 0.;
    double x2 = x[0]*x[0];
    double ratio = x2/ebeam2;
    if (par[7] < 0) ratio = ebeam2/x2;
    if (ratio < 1.) {
      background = par[5]*x[0] * TMath::Sqrt(1. - ratio) * TMath::Exp(par[6] * (1. - ratio));
    } else {
      background = 0.;
    }


    Double_t arg1(0.), arg2(0.), sig1(0.), sig2(0.);
    if (par[2] > 0) {
      arg1 = (x[0] - par[1]) / par[2];
      sig1 = par[0]*TMath::Exp(-0.5*arg1*arg1);
    }
    if (par[5] > 0.) {
      arg2 = (x[0] - par[1]) / par[4];
      sig2 = par[0]*TMath::Exp(-0.5*arg2*arg2);
    }

    return ((1-par[3])*sig1 + par[3]*sig2 + background);
  }




  // ----------------------------------------------------------------------
  double iF_gauss(double *x, double *par) {
    // par[0] -> const
    // par[1] -> mean
    // par[2] -> sigma

    if (par[2] > 0.) {
      Double_t arg = (x[0] - par[1]) / par[2];
      Double_t fitval =  par[0]*TMath::Exp(-0.5*arg*arg);
      return fitval;
    }
    else {
      return -1.;
    }
  }

  // ----------------------------------------------------------------------
  double iF_bigauss(double *x, double *par) {
    // par[0] -> const
    // par[1] -> mean
    // par[2] -> sigmaL
    // par[3] -> sigmaH

    if (x[0] > par[1]) {
      if (par[2] > 0.) {
	Double_t arg = (x[0] - par[1]) / par[2];
	Double_t fitval =  par[0]*TMath::Exp(-0.5*arg*arg);
	return fitval;
      }
      else {
	return -1.;
      }
    } else {
      if (par[3] > 0.) {
	Double_t arg = (x[0] - par[1]) / par[3];
	Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
	return fitval;
      }
      else {
	return -1.;
      }
    }
    return -1.;
  }


  // ----------------------------------------------------------------------
  double iF_Gauss(double *x, double *par) {
    // par[0] -> area
    // par[1] -> mean
    // par[2] -> sigma

    double sqrt2pi = 2.506628275;

    if (par[2] > 0.) {
      Double_t arg = (x[0] - par[1]) / par[2];
      Double_t fitval =  (par[0]/(sqrt2pi*par[2])) * TMath::Exp(-0.5*arg*arg);
      return fitval;
    }
    else {
      return -1.;
    }
  }

  // ----------------------------------------------------------------------
  double iF_gauss2(double *x, double *par) {
    // par[0] -> const
    // par[1] -> mean
    // par[2] -> sigma
    // par[3] -> fraction in second gaussian
    // par[4] -> mean of second gaussian
    // par[5] -> sigma of second gaussian
    Double_t arg1(0.), arg2(0.), fitval1(0.), fitval2(0.);
    if (par[2] > 0) {
      arg1 = (x[0] - par[1]) / par[2];
      fitval1 =  par[0]*TMath::Exp(-0.5*arg1*arg1);
    }
    if (par[5] > 0.) {
      arg2 = (x[0] - par[4]) / par[5];
      fitval2 =  par[3]*par[0]*TMath::Exp(-0.5*arg2*arg2);
    }
    Double_t fitval = fitval1 + fitval2;
    return fitval;
  }

  // ----------------------------------------------------------------------
  double iF_gauss2f(double *x, double *par) {
    // par[0] -> const
    // par[1] -> mean
    // par[2] -> sigma
    // par[3] -> fraction,  area ratio second/first gaussian
    // par[4] -> mean of second gaussian
    // par[5] -> sigma of second gaussian
    Double_t arg1(0.), arg2(0.), fitval1(0.), fitval2(0.);
    if (par[2] > 0) {
      arg1 = (x[0] - par[1]) / par[2];
      fitval1 =  par[0]*TMath::Exp(-0.5*arg1*arg1);
    }

    if (par[5] > 0.) {
      // relate the area of the 2 gaussians
      double fraction = par[0] * par[2]/par[5];
      arg2 = (x[0] - par[4]) / par[5];
      fitval2 =  par[3]*fraction*TMath::Exp(-0.5*arg2*arg2);
    }
    Double_t fitval = fitval1 + fitval2;
    return fitval;
  }

  // ----------------------------------------------------------------------
  double iF_gauss2c(double *x, double *par) {
    // constrained to have the same mean in the second gaussian
    // par[0] -> const
    // par[1] -> mean
    // par[2] -> sigma
    // par[3] -> fraction in second gaussian
    // par[4] -> sigma of second gaussian
    Double_t arg1(0.), arg2(0.), fitval1(0.), fitval2(0.);
    if (par[2] > 0) {
      arg1 = (x[0] - par[1]) / par[2];
      fitval1 =  par[0]*TMath::Exp(-0.5*arg1*arg1);
    }
    if (par[4] > 0.) {
      arg2 = (x[0] - par[1]) / par[4];
      fitval2 =  par[3]*par[0]*TMath::Exp(-0.5*arg2*arg2);
    }
    Double_t fitval = fitval1 + fitval2;
    return fitval;
  }

  // ----------------------------------------------------------------------
  double iF_gauss3(double *x, double *par) {
    // par[0] -> const
    // par[1] -> mean
    // par[2] -> sigma
    // par[3] -> fraction in second gaussian
    // par[4] -> mean of second gaussian
    // par[5] -> sigma of second gaussian
    // par[6] -> fraction in third gaussian
    // par[7] -> mean of third gaussian
    // par[8] -> sigma of third gaussian
    Double_t arg1(0.), arg2(0.), arg3(0.), fitval1(0.), fitval2(0.), fitval3(0.);
    if (par[2] > 0) {
      arg1 = (x[0] - par[1]) / par[2];
      fitval1 =  par[0]*TMath::Exp(-0.5*arg1*arg1);
    }
    if (par[5] > 0.) {
      arg2 = (x[0] - par[4]) / par[5];
      fitval2 =  par[3]*par[0]*TMath::Exp(-0.5*arg2*arg2);
    }
    if (par[8] > 0.) {
      arg3 = (x[0] - par[7]) / par[8];
      fitval3 =  par[6]*par[0]*TMath::Exp(-0.5*arg3*arg3);
    }

    Double_t fitval = fitval1 + fitval2 + fitval3;
    return fitval;
  }

  // ----------------------------------------------------------------------
  double iF_pol0(double *x, double *par) {
    return par[0];
  }


  // ----------------------------------------------------------------------
  double iF_pol1(double *x, double *par) {
    return par[0] + par[1]*x[0];
  }


  // ----------------------------------------------------------------------
  // double iF_pol2(double *x, double *par) {
  //   return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
  // }

  // ----------------------------------------------------------------------
  double iF_pol2local(double *x, double *par) {
    return par[0] + par[1]*(x[0]-par[2])*(x[0]-par[2]);
  }


  // ----------------------------------------------------------------------
  double iF_pol0_BsBlind(double *x, double *par) {
    // FIXME fixed limits!!!
    if (x[0] >= 5.2 && x[0] <= 5.45) {
      TF1::RejectPoint();
      return 0;
    }
    return par[0];
  }


  // ----------------------------------------------------------------------
  double iF_pol1_BsBlind(double *x, double *par) {
    // FIXME fixed limits!!!
    if (x[0] >= 5.2 && x[0] <= 5.45) {
      TF1::RejectPoint();
      return 0;
    }
    return par[0] + par[1]*x[0];
  }

  // ----------------------------------------------------------------------
  double iF_expo_BsBlind(double *x, double *par) {
    // FIXME fixed limits!!!
    if (x[0] >= 5.2 && x[0] <= 5.45) {
      TF1::RejectPoint();
      return 0;
    }
    return par[0]*TMath::Exp(x[0]*par[1]);
  }

  // ----------------------------------------------------------------------
  // pol1 and crystal ball
  double iF_pol1_cb(double *x, double *par) {
    // par[0]:  mean
    // par[1]:  sigma
    // par[2]:  alpha, crossover point
    // par[3]:  n, length of tail
    // par[4]:  N, normalization
    // par[5] = par 0 of pol1
    // par[6] = par 1 of pol1
    return  (iF_pol1(x, &par[5]) + iF_cb(x, &par[0]));
  }

  // ----------------------------------------------------------------------
  // gauss and crystal ball
  double iF_gauss_cb(double *x, double *par) {
    // par[0]:  norm  Gauss, ratio to xball normalization
    // par[1]:  peak  Gauss, ratio to xball peak
    // par[2]:  sigma Gauss, ratio to xball sigma
    // par[3]:  mean
    // par[4]:  sigma
    // par[5]:  alpha, crossover point
    // par[6]:  n, length of tail
    // par[7]:  N, normalization

    double fitval(-1.);
    if (par[2] > 0.) {
      Double_t arg = (x[0] - par[1]*par[3]) / (par[2]*par[4]);
      fitval = par[0]*par[7]*TMath::Exp(-0.5*arg*arg);
    } else {
      return -1.;
    }

    return  (fitval + iF_cb(x, &par[3]));
  }


  // ----------------------------------------------------------------------
  // bigauss and crystal ball
  double iF_bigauss_cb(double *x, double *par) {
    // par[0]:  norm  Gauss, ratio to xball normalization
    // par[1]:  peak  Gauss, ratio to xball peak
    // par[2]:  sigmaL Gauss, ratio to xball sigma
    // par[3]:  sigmaH Gauss, ratio to xball sigma
    // par[4]:  mean
    // par[5]:  sigma
    // par[6]:  alpha, crossover point
    // par[7]:  n, length of tail
    // par[8]:  N, normalization

    double fitval(-1.);
    if (x[0] < par[1]) {
      if (par[2] > 0.) {
	Double_t arg = (x[0] - par[1]*par[4]) / (par[2]*par[5]);
	fitval = par[0]*par[8]*TMath::Exp(-0.5*arg*arg);
      } else {
	return -1.;
      }
    } else {
      if (par[3] > 0.) {
	Double_t arg = (x[0] - par[1]*par[4]) / (par[3]*par[5]);
	fitval = par[0]*par[8]*TMath::Exp(-0.5*arg*arg);
      } else {
	return -1.;
      }
    }
    return  (fitval + iF_cb(x, &par[4]));
  }


  // ----------------------------------------------------------------------
  // bigauss and crystal ball and satellite peak
  double iF_bigauss_cb_sat(double *x, double *par) {
    // par[0]:  norm  Gauss, ratio to xball normalization
    // par[1]:  peak  Gauss, ratio to xball peak
    // par[2]:  sigmaL Gauss, ratio to xball sigma
    // par[3]:  sigmaH Gauss, ratio to xball sigma
    // par[4]:  mean
    // par[5]:  sigma
    // par[6]:  alpha, crossover point
    // par[7]:  n, length of tail
    // par[8]:  N, normalization
    // par[9]:  fraction in satellite peak

    // FIXME: This is likely a call for trouble to use par[4] before it's been used in 'its' function call below
    double fitval(-1.);
    if (x[0] < par[4]) {
      if (par[2] > 0.) {
	Double_t arg = (x[0] - /*par[1]* */par[4]) / (par[2]*par[5]);
	fitval = par[0]*par[8]*TMath::Exp(-0.5*arg*arg);
      } else {
	return -1.;
      }
    } else {
      if (par[3] > 0.) {
	Double_t arg = (x[0] - /*par[1]* */par[4]) / (par[3]*par[5]);
	fitval = par[0]*par[8]*TMath::Exp(-0.5*arg*arg);
      } else {
	return -1.;
      }
    }

    double cbval = iF_cb(x, &par[4]);
    double cbintegral = 2.50663 * par[8] * TMath::Abs(par[5]);

    // -- apppoximate integrals of 'signal' analytically
    // TMath::Sqrt(TMath::Pi()) * TMath::Sqrt(2.) = 2.50663;
    double gintegral = 2.50663 * par[0] * TMath::Abs(0.5*(par[2]+par[3]));


    // -- calculate B+ -> J/psi pi contribution:
    double g3par[] = {par[9]*cbintegral*gintegral, 5.36416e+00, 4.83167e-02,
		      3.39117e-01,   5.46405e+00, 9.83832e-02,
		      9.36584e-02,   5.52549e+00, 2.69386e-01};
    double peakingBg = iF_gauss3(x, g3par);

    return  (fitval + cbval + peakingBg);
  }


  // ----------------------------------------------------------------------
  // pol0 and gauss
  double iF_pol0_gauss(double *x, double *par) {
    //   par[0] = normalization of gaussian
    //   par[1] = mean of gaussian
    //   par[2] = sigma of gaussian
    //   par[3] = par 0 of pol0
    return  (iF_pol0(x, &par[3]) + iF_gauss(x, &par[0]));
  }

  // ----------------------------------------------------------------------
  // pol1 and gauss
  double iF_pol1_gauss(double *x, double *par) {
    //   par[0] = normalization of gaussian
    //   par[1] = mean of gaussian
    //   par[2] = sigma of gaussian
    //   par[3] = par 0 of pol1
    //   par[4] = par 1 of pol1
    return  (iF_pol1(x, &par[3]) + iF_gauss(x, &par[0]));
  }


  // ----------------------------------------------------------------------
  // expo and gauss
  double iF_expo_Gauss(double *x, double *par) {
    //   par[0] = normalization of gaussian
    //   par[1] = mean of gaussian
    //   par[2] = sigma of gaussian
    //   par[3] = base of expo
    //   par[4] = exponent of expo
    return  (iF_expo(x, &par[3]) + iF_Gauss(x, &par[0]));
  }

  // ----------------------------------------------------------------------
  // expo and err and gauss
  double iF_expo_err_Gauss(double *x, double *par) {
    //   par[0] = normalization of gaussian
    //   par[1] = mean of gaussian
    //   par[2] = sigma of gaussian
    //   par[3] = base
    //   par[4] = exp
    //   par[5] = par[0] of err
    //   par[6] = par[1] of err
    //   par[7] = par[2] of err
    //   par[8] = par[3] of err
    return  (iF_err(x, &par[5]) + iF_expo(x, &par[3]) + iF_Gauss(x, &par[0]));
  }

  // ----------------------------------------------------------------------
  // expo and err2 and gauss
  double iF_expo_err2_Gauss(double *x, double *par) {
    //   par[0] = normalization of gaussian
    //   par[1] = mean of gaussian
    //   par[2] = sigma of gaussian
    //   par[3] = base
    //   par[4] = exp
    //   par[5] = par[0] of err2
    //   par[6] = par[1] of err2
    //   par[7] = par[2] of err2
    return  (iF_err2(x, &par[5]) + iF_expo(x, &par[3]) + iF_Gauss(x, &par[0]));
  }

  // ----------------------------------------------------------------------
  // expo and err and gauss2c
  double iF_expo_err_gauss2c(double *x, double *par) {
    // par[0] -> const
    // par[1] -> mean
    // par[2] -> sigma
    // par[3] -> fraction in second gaussian
    // par[4] -> sigma of second gaussian
    //   par[5]  = base
    //   par[6]  = exp
    //   par[7]  = par[0] of err
    //   par[8]  = par[1] of err
    //   par[9]  = par[2] of err
    //   par[10] = par[3] of err
    return  (iF_err(x, &par[7]) + iF_expo(x, &par[5]) + iF_gauss2c(x, &par[0]));
  }

  // ----------------------------------------------------------------------
  // expo and err and gauss2c
  double iF_expo_err_gauss2(double *x, double *par) {
    // par[0] -> const
    // par[1] -> mean
    // par[2] -> sigma
    // par[3] -> fraction in second gaussian
    // par[4] -> mean of second gaussian
    // par[5] -> sigma of second gaussian
    //   par[6]  = base
    //   par[7]  = exp
    //   par[8]  = par[0] of err
    //   par[9]  = par[1] of err
    //   par[10] = par[2] of err
    //   par[11] = par[3] of err


    return  (iF_err(x, &par[8]) + iF_expo(x, &par[6]) + iF_gauss2(x, &par[0]));
  }

  // ----------------------------------------------------------------------
  // expo and err and gauss2f (2 gaussians correlated through the area)
  double iF_expo_err_gauss2f(double *x, double *par) {
    // par[0] -> const
    // par[1] -> mean
    // par[2] -> sigma
    // par[3] -> fraction in second gaussian
    // par[4] -> mean of second gaussian
    // par[5] -> sigma of second gaussian
    //   par[6]  = base
    //   par[7]  = exp
    //   par[8]  = par[0] of err
    //   par[9]  = par[1] of err
    //   par[10] = par[2] of err
    //   par[11] = par[3] of err


    return  (iF_err(x, &par[8]) + iF_expo(x, &par[6]) + iF_gauss2f(x, &par[0]));
  }

  // ----------------------------------------------------------------------
  // expo and err and gauss2 and landau
  double iF_expo_err_gauss2_landau(double *x, double *par) {
    // par[0] -> const
    // par[1] -> mean
    // par[2] -> sigma
    // par[3] -> fraction in second gaussian
    // par[4] -> mean of second gaussian
    // par[5] -> sigma of second gaussian
    //   par[6]  = base
    //   par[7]  = exp
    //   par[8]  = par[0] of err
    //   par[9]  = par[1] of err
    //   par[10] = par[2] of err
    //   par[11] = par[3] of err
    //   par[12] = par[0] of landau (mpv)
    //   par[13] = par[1] of landau (sigma)
    //   par[14] = par[2] of landau (const) WILL BE FIXED BY THE AREA RATIO

    // paramter 3 of the landau, the constant, has to be fixed by the BuJpsiPi / BuJpsiK ratio
    const double branchingRatio = 0.049;
    double area = 2.507 * par[0] * (par[2] + par[3]*par[5]); // area of the 2 gaussians
    double area_landau = branchingRatio * area;  // expected area of the landau
    par[14] = area_landau/par[13];  // landau const = area/sigma, feed it to the function
    //cout<<par[14]<<" "<<par[13]<<" "<<par[12]<<endl;

    return  ( iF_err(x, &par[8]) + iF_expo(x, &par[6]) + iF_gauss2(x, &par[0]) + iF_landausimp(x, &par[12]) );
  }

  // ----------------------------------------------------------------------
  // expo and err and gauss and landau
  double iF_expo_err_Gauss_landau(double *x, double *par) {
    // par[0] -> const
    // par[1] -> mean
    // par[2] -> sigma
    //   par[3]  = base
    //   par[4]  = exp
    //   par[5]  = par[0] of err
    //   par[6]  = par[1] of err
    //   par[7] = par[2] of err
    //   par[8] = par[3] of err
    //   par[9] = par[0] of landau (mpv)
    //   par[10] = par[1] of landau (sigma)
    //   par[11] = par[2] of landau (const) WILL BE FIXED BY THE AREA RATIO

    // paramter 3 of the landau, the constant, has to be fixed by the BuJpsiPi / BuJpsiK ratio
    const double branchingRatio = 0.049;
    double area = 2.507 * par[0] * (par[2]); // area of the gaussian
    double area_landau = branchingRatio * area;  // expected area of the landau
    par[11] = area_landau/par[10];  // landau const = area/sigma, feed it to the function

    return  (iF_err(x, &par[5]) + iF_expo(x, &par[3]) + iF_Gauss(x, &par[0]) + iF_landausimp(x, &par[9]) );
  }


  // ----------------------------------------------------------------------
  // expo and err
  double iF_expo_err(double *x, double *par) {
    //   par[0] = base
    //   par[1] = exp
    //   par[2] = par[0] of err
    //   par[3] = par[1] of err
    //   par[4] = par[2] of err
    //   par[5] = par[3] of err
    return  (iF_err(x, &par[2]) + iF_expo(x, &par[0]));
  }

  // ----------------------------------------------------------------------
  // expo and err2
  double iF_expo_err2(double *x, double *par) {
    //   par[0] = base
    //   par[1] = exp
    //   par[2] = par[0] of err2
    //   par[3] = par[1] of err2
    //   par[4] = par[2] of err2
    return  (iF_err2(x, &par[2]) + iF_expo(x, &par[0]));
  }


  // ----------------------------------------------------------------------
  // expo and err and Gauss
  double iF_pol1_err_Gauss(double *x, double *par) {
    //   par[0] = normalization of gaussian
    //   par[1] = mean of gaussian
    //   par[2] = sigma of gaussian
    //   par[3] = const
    //   par[4] = slope
    //   par[5] = par[0] of err
    //   par[6] = par[1] of err
    //   par[7] = par[2] of err
    //   par[8] = par[3] of err
    return  (iF_err(x, &par[5]) + iF_pol1(x, &par[3]) + iF_Gauss(x, &par[0]));
  }

  // ----------------------------------------------------------------------
  // expo and err2 and Gauss
  double iF_pol1_err2_Gauss(double *x, double *par) {
    //   par[0] = normalization of gaussian
    //   par[1] = mean of gaussian
    //   par[2] = sigma of gaussian
    //   par[3] = const
    //   par[4] = slope
    //   par[5] = par[0] of err
    //   par[6] = par[1] of err
    //   par[7] = par[2] of err
    return  (iF_err2(x, &par[5]) + iF_pol1(x, &par[3]) + iF_Gauss(x, &par[0]));
  }

  // ----------------------------------------------------------------------
  // pol1 and err2 and gauss2
  double iF_pol1_err2_gauss2c(double *x, double *par) {
    //   par[0]  = const of gaussian
    //   par[1]  = mean of gaussian
    //   par[2]  = sigma of gaussian
    //   par[3]  = fraction in second gaussian
    //   par[4]  = sigma of second gaussian
    //   par[5]  = const
    //   par[6]  = slope
    //   par[7]  = par[0] of err2
    //   par[8]  = par[1] of err2
    //   par[9]  = par[2] of err2
    // cout << x[0] << ": " << iF_err2(x, &par[7]) << " " << iF_pol1(x, &par[5]) << " " << iF_gauss2c(x, &par[0])
    // 	 << " -> " << (iF_err2(x, &par[7]) + iF_pol1(x, &par[5]) + iF_gauss2c(x, &par[0]))
    // 	 << endl;
    return  (iF_err2(x, &par[7]) + iF_pol1(x, &par[5]) + iF_gauss2c(x, &par[0]));
  }

  // ----------------------------------------------------------------------
  // expo and err and gauss
  double iF_pol1_err(double *x, double *par) {
    //   par[0] = const
    //   par[1] = slope
    //   par[2] = par[0] of err
    //   par[3] = par[1] of err
    //   par[4] = par[2] of err
    //   par[5] = par[3] of err
    return  (iF_err(x, &par[2]) + iF_pol1(x, &par[0]));
  }


  // ----------------------------------------------------------------------
  // pol1 and Gauss
  double iF_pol1_Gauss(double *x, double *par) {
    //   par[0] = area of gaussian
    //   par[1] = mean of gaussian
    //   par[2] = sigma of gaussian
    //   par[3] = par 0 of pol1
    //   par[4] = par 1 of pol1
    return  (iF_pol1(x, &par[3]) + iF_Gauss(x, &par[0]));
  }

  // ----------------------------------------------------------------------
  // pol0 and Gauss
  double iF_pol0_Gauss(double *x, double *par) {
    //   par[0] = area of gaussian
    //   par[1] = mean of gaussian
    //   par[2] = sigma of gaussian
    //   par[3] = par 0 of pol1
    return  (iF_pol1(x, &par[3]) + iF_Gauss(x, &par[0]));
  }

  // ----------------------------------------------------------------------
  // pol1 and gauss2c
  double iF_pol1_gauss2c(double *x, double *par) {
    //   par[0] = norm of gaussian
    //   par[1] = mean of gaussian
    //   par[2] = sigma of gaussian
    //   par[3] = fraction in second gaussian
    //   par[4] = sigma of gaussian
    //   par[5] = par 0 of pol1
    //   par[6] = par 1 of pol1
    return  (iF_pol1(x, &par[5]) + iF_gauss2c(x, &par[0]));
  }

  // ----------------------------------------------------------------------
  // pol1 and gauss2
  double iF_pol1_gauss2(double *x, double *par) {
    // par[0] -> const
    // par[1] -> mean
    // par[2] -> sigma
    // par[3] -> fraction in second gaussian
    // par[4] -> mean of second gaussian
    // par[5] -> sigma of second gaussian
    // par[6] = par 0 of pol1
    // par[7] = par 1 of pol1
    return  (iF_pol1(x, &par[6]) + iF_gauss2(x, &par[0]));
  }

  // ----------------------------------------------------------------------
  // expo and gauss2
  double iF_expo_gauss2(double *x, double *par) {
    // par[0] -> const
    // par[1] -> mean
    // par[2] -> sigma
    // par[3] -> fraction in second gaussian
    // par[4] -> mean of second gaussian
    // par[5] -> sigma of second gaussian
    // par[6] = base
    // par[7] = exp
    return  (iF_expo(x, &par[6]) + iF_gauss2(x, &par[0]));
  }

  // ----------------------------------------------------------------------
  // expo and gauss2c
  double iF_expo_gauss2c(double *x, double *par) {
    //   par[0] = norm of gaussian
    //   par[1] = mean of gaussian
    //   par[2] = sigma of gaussian
    //   par[3] = fraction in second gaussian
    //   par[4] = sigma of gaussian
    //   par[5] = base
    //   par[6] = exp
    return  (iF_expo(x, &par[5]) + iF_gauss2c(x, &par[0]));
  }


  // ----------------------------------------------------------------------
  // pol1 and landau
  double iF_pol1_landau(double *x, double *par) {
    //   par[0] = mpv
    //   par[1] = sigma
    //   par[2] = norm
    //   par[3] = par 0 of pol1
    //   par[4] = par 1 of pol1
    return  (iF_pol1(x, &par[3]) + iF_landau(x, &par[0]));
  }


  // ----------------------------------------------------------------------
  // satellite peak from B+ -> J/psi pi+
  double iF_pisat(double *x, double *par) {
    double g3par[] = {par[0], 5.36416e+00, 4.83167e-02, 3.39117e-01, 5.46405e+00,
		      9.83832e-02, 9.36584e-02, 5.52549e+00, 2.69386e-01};
    return iF_gauss3(x, g3par);
  }

  // ----------------------------------------------------------------------
  // satellite peak from B0 -> J/psi Kstar0 (for Bs -> J/psi phi)
  double iF_kstarsat(double *x, double *par) {
    // double g2pol1par[] = {par[0], 5.381, 0.04847, 0.4206, 5.415,
    // 			  0.1104, 19.81, -3.269};
    double g2pol1par[] = {par[0], 5.381, 0.04847, 0.4206, 5.415,
			  0.1104, 0., 0.};
    return iF_pol1_gauss2(x, g2pol1par);
  }

  // ----------------------------------------------------------------------
  // double gaussian with same mean
  double iF_gauss2c_sat(double *x, double *par) {
    // constrained to have the same mean in the second gaussian
    // par[0] -> const
    // par[1] -> mean
    // par[2] -> sigma
    // par[3] -> fraction in second gaussian
    // par[4] -> sigma of second gaussian
    double fracSat(0.04);

    Double_t arg1(0.), arg2(0.), fitval1(0.), fitval2(0.);
    if (par[2] > 0) {
      arg1 = (x[0] - par[1]) / par[2];
      fitval1 =  par[0]*TMath::Exp(-0.5*arg1*arg1);
    }
    if (par[4] > 0.) {
      arg2 = (x[0] - par[1]) / par[4];
      fitval2 =  par[3]*par[0]*TMath::Exp(-0.5*arg2*arg2);
    }
    Double_t fitval = fitval1 + fitval2;

    // -- double Gaussian integral:
    double sqrt2pi = 2.506628275;
    double gintegral = sqrt2pi*par[0]*(par[2] + par[3]*par[4]);

    // -- B+ -> J/psi pi normalization:
    double satintegral = initFunc::iF_int_pisat(1.);
    double norm        = fracSat*gintegral/satintegral;

    double g3par[] = {norm, 5.36416e+00, 4.83167e-02, 3.39117e-01, 5.46405e+00,
		      9.83832e-02, 9.36584e-02, 5.52549e+00, 2.69386e-01};
    double peakingBg = iF_gauss3(x, g3par);
    double sintegral = initFunc::iF_int_pisat(norm);

    if (0) cout << " satintegral = " << satintegral
		<< " gintegral = " << gintegral
		<< " norm = " << norm
		<< " sintegral = " << sintegral
		<< " frac = " << sintegral/gintegral
		<< endl;

    return  (fitval + peakingBg);
  }


  // ----------------------------------------------------------------------
  // double gaussian with same mean plus kstar satellite peak
  double iF_gauss2c_kstarsat(double *x, double *par) {
    // constrained to have the same mean in the second gaussian
    // par[0] -> const
    // par[1] -> mean
    // par[2] -> sigma
    // par[3] -> fraction in second gaussian
    // par[4] -> sigma of second gaussian
    // par[5] -> normalization of K*0 satellite peak

    Double_t arg1(0.), arg2(0.), fitval1(0.), fitval2(0.);
    if (par[2] > 0) {
      arg1 = (x[0] - par[1]) / par[2];
      fitval1 =  par[0]*TMath::Exp(-0.5*arg1*arg1);
    }
    if (par[4] > 0.) {
      arg2 = (x[0] - par[1]) / par[4];
      fitval2 =  par[3]*par[0]*TMath::Exp(-0.5*arg2*arg2);
    }
    Double_t fitval = fitval1 + fitval2;

    double peakingBg = 0.;
    double parArr[] = {par[5]};
    if (4.5 < x[0] && x[0] < 6.) {
      peakingBg = iF_kstarsat(x, parArr);
    }

    return  (fitval + peakingBg);
  }


  // ----------------------------------------------------------------------
  // B+ -> J/psi K+ pdf
  double iF_bupsik(double *x, double *par) {
    return (iF_gauss2c_sat(x, &par[0]) + iF_expo(x, &par[5]) + iF_err2(x, &par[7]));
  }

  // ----------------------------------------------------------------------
  // B+ -> J/psi K+ pdf (too complicated)
  double iF_bupsik1(double *x, double *par) {
    return (iF_bigauss_cb_sat(x, &par[0]) + iF_expo(x, &par[10]) + iF_err2(x, &par[12]));
  }


  // ----------------------------------------------------------------------
  // Bs-> J/psi phi pdf
  double iF_bspsiphi(double *x, double *par) {
    return (iF_gauss2c_kstarsat(x, &par[0]) + iF_expo_HS(x, &par[6]));
  }



}


// ----------------------------------------------------------------------
initFunc::initFunc(): fLo(99.), fHi(-99.), fBgFractionLo(-99), fBgFractionHi(-99.), fVerbose(true), fName("IF") {
  resetLimits();
}

// ----------------------------------------------------------------------
initFunc::initFunc(string name): fLo(99.), fHi(-99.), fBgFractionLo(-99), fBgFractionHi(-99.), fVerbose(true), fName(name) {
  resetLimits();
}


// ----------------------------------------------------------------------
initFunc::~initFunc() {
  //  cout << "dtor initFunc" << endl;
}


// ----------------------------------------------------------------------
void initFunc::limitPar(int i, double lo, double hi) {
  fLimit[i]   = true;
  fFix[i]     = false;
  fLimitLo[i] = lo;
  fLimitHi[i] = hi;
}

// ----------------------------------------------------------------------
void initFunc::fixPar(int i, double fix) {
  fLimit[i]   = false;
  fFix[i]     = true;
  fLimitLo[i] = fix;
  fLimitHi[i] = fix;
}


// ----------------------------------------------------------------------
void initFunc::resetLimits() {
  for (int i = 0; i < 20; ++i) {
    fLimitLo[i] = 0.;
    fLimitHi[i] = -1.;
    fLimit[i]   = false;
    fFix[i]     = false;
  }
}


// ----------------------------------------------------------------------
void initFunc::dumpParameters(TF1 *f1) {
  double lo, hi;
  cout << "initFunc(" << fName << ") dumpParameters: " << f1->GetName() << endl;
  for (int i = 0; i < f1->GetNpar(); ++i) {
    f1->GetParLimits(i, lo, hi);
    cout << Form("%2d: %f (%f .. %f)", i, f1->GetParameter(i), lo, hi)
	 << endl;
  }
}



// ----------------------------------------------------------------------
void initFunc::applyLimits(TF1 *f, string name) {
  if (fVerbose) cout << "initFunc(" << fName << ") applyLimits:" << endl;
  for (int i = 0; i < f->GetNpar(); ++i) {
    if (fFix[i]) {
      if (fVerbose) cout << "initFunc::" << name << "> fixing par " << i << " to " << fLimitLo[i] << endl;
      f->FixParameter(i, fLimitLo[i]);
    } else {
      if (fLimit[i]) {
	if (fVerbose) cout << "initFunc::" << name << "> limiting par " << i << " from " << fLimitLo[i] << " .. " << fLimitHi[i]
			   << ", value = " << f->GetParameter(i)
			   << endl;
	f->SetParLimits(i, fLimitLo[i], fLimitHi[i]);
      } else {
	if (fVerbose) cout << "initFunc::" << name << "> releasing par " << i << endl;
	f->ReleaseParameter(i);
      }
    }
  }
}


// ----------------------------------------------------------------------
void initFunc::dumpLimits(string bla) {
  cout << "initFunc(" << fName << ") dumpLimits(" << bla << "):" << endl;
  for (int i = 0; i < 20; ++i) {
    cout << Form("%2d limits: %+5.3f ... %+5.3f, limited: %c, fixed: %c",
		 i, fLimitLo[i], fLimitHi[i], (fLimit[i]?'y':'n'), (fFix[i]?'y':'n'))
	 << endl;
  }

}

// ----------------------------------------------------------------------
TF1* initFunc::err(double lo, double hi) {
  TF1 *f = new TF1(Form("%s_err", fName.c_str()), iF_err, lo, hi, 4);
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::err(TH1 *h) {
  if (0 == h) {
    cout << "empty histogram pointer" << endl;
    return 0;
  }
  TF1 *f = new TF1(Form("%s_err", fName.c_str()), iF_err, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 4);

  return f;
}

// ----------------------------------------------------------------------
TF1* initFunc::err2(double lo, double hi) {
  TF1 *f = new TF1(Form("%s_err2", fName.c_str()), iF_err2, lo, hi, 3);
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::err2(TH1 *h) {
  if (0 == h) {
    cout << "empty histogram pointer" << endl;
    return 0;
  }
  TF1 *f = new TF1(Form("%s_err2", fName.c_str()), iF_err2, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 3);

  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::pol0(double lo, double hi) {
  TF1 *f = new TF1(Form("%s_pol0", fName.c_str()), iF_pol0, lo, hi, 1);
  return f;
}

// ----------------------------------------------------------------------
TF1* initFunc::pol1(double lo, double hi) {
  TF1 *f = new TF1(Form("%s_pol1", fName.c_str()), iF_pol1, lo, hi, 2);
  return f;
}

// ----------------------------------------------------------------------
TF1* initFunc::pol1Err(double lo, double hi) {
  TF1 *f = new TF1(Form("%s_pol1Err", fName.c_str()), iF_pol1_err, lo, hi, 6);
  return f;
}

// ----------------------------------------------------------------------
TF1* initFunc::expo(double lo, double hi) {
  TF1 *f = new TF1(Form("%s_expo", fName.c_str()), iF_expo, lo, hi, 2);
  f->SetParNames("base", "expo");
  return f;
}

// ----------------------------------------------------------------------
TF1* initFunc::expoErr(double lo, double hi) {
  TF1 *f = new TF1(Form("%s_expoErr", fName.c_str()), iF_expo_err, lo, hi, 6);
  return f;
}

// ----------------------------------------------------------------------
TF1* initFunc::expoErr2(double lo, double hi) {
  TF1 *f = new TF1(Form("%s_expoErr2", fName.c_str()), iF_expo_err2, lo, hi, 5);
  return f;
}

// ----------------------------------------------------------------------
TF1* initFunc::argus(double lo, double hi) {
  TF1 *f = new TF1(Form("%s_argus", fName.c_str()), iF_argus, lo, hi, 3);
  f->SetParNames("norm.", "expo.", "endpoint");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::landau(double lo, double hi) {
  TF1 *f = new TF1(Form("%s_landau", fName.c_str()), iF_landau, lo, hi, 3);
  f->SetParName(0, "mpvl");
  f->SetParName(1, "sigl");
  f->SetParName(2, "norm");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::landau(TH1 *h) {
  if (0 == h) {
    cout << "empty histogram pointer" << endl;
    return 0;
  }
  TF1 *f = new TF1(Form("%s_landau", fName.c_str()), iF_landau, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 3);
  f->SetParName(0, "mpvl");
  f->SetParName(1, "sigl");
  f->SetParName(2, "norm");
  f->SetParameter(0, h->GetBinCenter(h->GetMaximumBin()));
  f->SetParameter(1, h->GetRMS());
  f->SetParameter(2, h->GetMaximum());
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::bigauss(double lo, double hi) {
  TF1 *f = new TF1(Form("%s_bigauss", fName.c_str()), iF_bigauss, lo, hi, 4);
  f->SetParName(0, "norm");
  f->SetParName(1, "peak");
  f->SetParName(2, "sigl");
  f->SetParName(3, "sigh");
  return f;
}

// ----------------------------------------------------------------------
TF1* initFunc::bigauss(TH1 *h) {
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_bigauss", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_bigauss", fName.c_str()), iF_bigauss, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 4);
  f->SetParName(0, "norm");
  f->SetParName(1, "peak");
  f->SetParName(2, "sigl");
  f->SetParName(3, "sigh");
  f->SetParameters(h->GetMaximum(), h->GetMean(), h->GetRMS(), h->GetRMS());
  return f;

}

// ----------------------------------------------------------------------
TF1* initFunc::Gauss(double lo, double hi) {
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_Gauss", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_gauss", fName.c_str()), iF_Gauss, lo, hi, 3);
  f->SetParName(0, "area");
  f->SetParName(1, "peak");
  f->SetParName(2, "sigma");
  return f;
}

// ----------------------------------------------------------------------
TF1* initFunc::gauss(double lo, double hi) {
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_gauss", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_gauss", fName.c_str()), iF_gauss, lo, hi, 3);
  f->SetParName(0, "const");
  f->SetParName(1, "peak");
  f->SetParName(2, "sigma");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::gauss(TH1 *h) {
  if (0 == h) {
    cout << "empty histogram pointer" << endl;
    return 0;
  }
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_gauss", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_gauss", fName.c_str()), iF_gauss, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1),  3);
  f->SetParName(0, "const");
  f->SetParName(1, "peak");
  f->SetParName(2, "sigma");

  f->SetParameters(h->GetMaximum(), h->GetMean(), h->GetRMS());
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::gauss2(double lo, double hi) {
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_gauss2", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_gauss2", fName.c_str()), iF_gauss2, lo, hi, 6);
  f->SetParNames("const", "peak", "sigma", "fraction", "mean2", "sigma2");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::gauss2(TH1 *h) {
  if (0 == h) {
    cout << "empty histogram pointer" << endl;
    return 0;
  }
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_gauss2", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_gauss2", fName.c_str()), iF_gauss2, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 6);
  f->SetParNames("const", "peak", "sigma", "fraction", "mean2", "sigma2");
  f->SetParameters(h->GetMaximum(), h->GetMean(), 0.2*h->GetRMS(), 0.2, h->GetMean(), 0.6*h->GetRMS());
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::gauss2c(double lo, double hi) {
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_gauss2c", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_gauss2c", fName.c_str()), iF_gauss2c, lo, hi, 5);
  f->SetParNames("const", "peak", "sigma", "fraction", "sigma2");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::gauss2c(TH1 *h) {
  if (0 == h) {
    cout << "empty histogram pointer" << endl;
    return 0;
  }
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_gauss2c", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_gauss2c", fName.c_str()), iF_gauss2c, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 5);
  f->SetParNames("const", "peak", "sigma", "fraction", "sigma2");
  f->SetParameters(h->GetMaximum(), h->GetMean(), 0.2*h->GetRMS(), 0.2, 0.6*h->GetRMS());
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::gauss3(double lo, double hi) {
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_gauss3", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_gauss3", fName.c_str()), iF_gauss3, lo, hi, 9);
  f->SetParNames("const", "peak", "sigma", "fraction2", "mean2", "sigma2", "fraction3", "mean3", "sigma3");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::gauss3(TH1 *h) {
  if (0 == h) {
    cout << "empty histogram pointer" << endl;
    return 0;
  }
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_gauss3", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_gauss3", fName.c_str()), iF_gauss3, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 9);
  f->SetParNames("const", "peak", "sigma", "fraction2", "mean2", "sigma2", "fraction3", "mean3", "sigma3");
  f->SetParameters(h->GetMaximum(), h->GetMean(), 0.2*h->GetRMS(),
		   0.3, h->GetMean(), 0.6*h->GetRMS(),
		   0.1, h->GetMean(), 1.0*h->GetRMS()
		   );
  return f;
}

// ----------------------------------------------------------------------
TF1* initFunc::pisat(double norm) {
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_pisat", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_pisat", fName.c_str()), iF_pisat, 0., 100., 1);
  f->SetParNames("norm");
  f->SetParameter(0, norm);
  return f;
}

// ----------------------------------------------------------------------
TF1* initFunc::kstarsat(double norm) {
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_kstarsat", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_kstarsat", fName.c_str()), iF_kstarsat, 0., 100., 1);
  f->SetParNames("satnorm");
  f->SetParameter(0, norm);
  return f;
}

// ----------------------------------------------------------------------
TF1* initFunc::pol1(TH1 *h) {
  if (0 == h) {
    cout << "empty histogram pointer" << endl;
    return 0;
  }
  TF1 *f = new TF1(Form("%s_pol1", fName.c_str()), iF_pol1, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 2);
  f->SetParNames("offset", "slope");
  //   cout << "Created f1 from " << h->GetBinLowEdge(1) << " to " << h->GetBinLowEdge(h->GetNbinsX()+1) << endl;

  double p1(0.), p0(0.);
  initPol1(p0, p1, h);
  f->SetParameters(p0, p1);

  //   cout << "  initialized to " << p0 << " and " << p1 << endl;
  //   cout << "  -> Integrals:  func = "
  //        << f->Integral(h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1))/h->GetBinWidth(1)
  //        << " histogram = "
  //        << h->GetSumOfWeights()
  //        << endl;
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::pol1(TH1 *h, double lo, double hi) {
  fLo = lo;
  fHi = hi;
  return pol1(h);
}


// ----------------------------------------------------------------------
TF1* initFunc::expo(TH1 *h) {
  if (0 == h) {
    cout << "empty histogram pointer" << endl;
    return 0;
  }
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_expo", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_expo", fName.c_str()), iF_expo, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 2);
  f->SetParNames("base", "exp");
  //  cout << "Created f1 from " << h->GetBinLowEdge(1) << " to " << h->GetBinLowEdge(h->GetNbinsX()+1) << endl;

  double p1(0.), p0(0.);
  initExpo(p0, p1, h);
  f->SetParameters(p0, p1);
  f->Update();
  if (1) {
    cout << Form("  expo initialized to %f and %f", p0,  p1)
	 << "  -> Integrals:  func = "
	 << f->Integral(h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1))/h->GetBinWidth(1)
	 << " histogram = "
	 << h->GetSumOfWeights()
	 << endl;
  }
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::expo(TH1 *h, double lo, double hi) {
  fLo = lo;
  fHi = hi;
  return expo(h);
}


// ----------------------------------------------------------------------
TF1* initFunc::argus(TH1 *h) {
  if (0 == h) {
    cout << "empty histogram pointer" << endl;
    return 0;
  }
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_argus", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_argus", fName.c_str()), iF_argus, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 3);
  f->SetParNames("norm.", "expo.", "endpoint");
  f->SetParameter(0, h->GetMaximum());
  f->SetParameter(1, h->GetMaximum());
  f->SetParameter(2, 5.28);

  return f;
}



// ----------------------------------------------------------------------
TF1* initFunc::argus(TH1 *h, double lo, double hi) {
  fLo = lo;
  fHi = hi;
  return argus(h);
}

// ----------------------------------------------------------------------
TF1* initFunc::pol1BsBlind(TH1 *h) {
  if (0 == h) {
    cout << "empty histogram pointer" << endl;
    return 0;
  }
  TF1 *f = new TF1(Form("%s_pol1bsblind", fName.c_str()), iF_pol1_BsBlind, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 2);
  f->SetParNames("offset", "slope");
  //   cout << "Created f1 from " << h->GetBinLowEdge(1) << " to " << h->GetBinLowEdge(h->GetNbinsX()+1) << endl;

  int lbin(1), hbin(h->GetNbinsX()+1), EDG(4), NB(EDG+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  //   cout << "lbin = " << lbin << " hbin = " << hbin << endl;

  double ylo = h->Integral(lbin, lbin+EDG)/NB;
  double yhi = h->Integral(hbin-EDG, hbin)/NB;
  double dx = h->GetBinLowEdge(hbin) - h->GetBinLowEdge(lbin);

  //   cout << "ylo = " << ylo << " yhi = " << yhi << " dx = " << dx << endl;

  double p1 = (yhi-ylo)/dx;
  double p0 = ylo - p1*h->GetBinLowEdge(1);
  f->SetParameters(p0, p1);

  double yintlo = h->Integral(h->FindBin(4.90), h->FindBin(5.20));
  double yinthi = h->Integral(h->FindBin(5.45), h->FindBin(5.90));
  if (yinthi < 1.e-03 ||yintlo < 1.e-03) {
    f->SetParameter(1, 0.);
    f->FixParameter(1, 0.);
    cout << "pol1BsBlind:  fixed slope because of zero entries in one sideband" << endl;
  }


  //   cout << "  initialized to " << p0 << " and " << p1 << endl;
  //   cout << "  -> Integrals:  func = "
  //        << f->Integral(h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1))/h->GetBinWidth(1)
  //        << " histogram = "
  //        << h->GetSumOfWeights()
  //        << endl;
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::expoBsBlind(TH1 *h) {
  if (0 == h) {
    cout << "empty histogram pointer" << endl;
    return 0;
  }
  TF1 *f = (TF1*)gROOT->FindObject(Form("%s_expobsblind", fName.c_str()));
  if (f) delete f;
  int npar(2);
  f = new TF1(Form("%s_expobsblind", fName.c_str()), iF_expo_BsBlind, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), npar);
  f->SetParNames("offset", "slope");

  int lbin(1), hbin(h->GetNbinsX()+1), EDG(4), NB(EDG+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  double p0(50.), p1(-0.5);

  double dx = h->GetBinLowEdge(hbin) - h->GetBinLowEdge(lbin);
  double ylo = h->Integral(lbin, lbin+EDG)/NB;
  double yhi = h->Integral(hbin-EDG, hbin)/NB;

  double yintlo = h->Integral(h->FindBin(4.90), h->FindBin(5.20));
  double yinthi = h->Integral(h->FindBin(5.45), h->FindBin(5.90));

  if (yhi > 0) {
    p1 = (TMath::Log(yhi) - TMath::Log(ylo))/dx;
    p0 = ylo/TMath::Exp(p1*fLo);
  } else {
    p0 = 50.;
    if (yhi > ylo) {
      p1 = 1.;
    } else {
      p1 = -1.*yinthi/yintlo;
    }
    cout << "initFunc::expoBsBlind: using default parameters for function initialization!" << endl;
  }
  if (fVerbose) {
    cout << "fLo: " << fLo << " fHi: " << fHi << " dx = " << dx << endl;
    cout << "ylo: " << ylo << " yhi: " << yhi << " log: " << TMath::Log(yhi) << " .. " << TMath::Log(ylo) << endl;
    cout << "p0:  " << p0  << " p1:  " << p1 << endl;
  }

  f->SetParameters(p0, p1);
  applyLimits(f, "expoBsBlind");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::pol0(TH1 *h) {
  if (0 == h) {
    cout << "empty histogram pointer" << endl;
    return 0;
  }
  TF1 *f = new TF1(Form("%s_pol0", fName.c_str()), iF_pol0, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 1);
  f->SetParNames("constant");
  int lbin(1), hbin(h->GetNbinsX()), EDG(4), NB(EDG+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }
  double ylo = h->Integral(lbin, lbin+EDG)/NB;
  double yhi = h->Integral(hbin-EDG, hbin)/NB;
  double p0 = 0.5*(yhi+ylo);

  f->SetParameter(0, p0);
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::pol0BsBlind(TH1 *h) {
  if (0 == h) {
    cout << "empty histogram pointer" << endl;
    return 0;
  }
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_pol0bsblind", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_pol0bsblind", fName.c_str()), iF_pol0_BsBlind, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 1);
  f->SetParNames("constant");
  int lbin(1), hbin(h->GetNbinsX()), EDG(4), NB(EDG+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }
  double ylo = h->Integral(lbin, lbin+EDG)/NB;
  double yhi = h->Integral(hbin-EDG, hbin)/NB;
  double p0 = 0.5*(yhi+ylo);

  f->SetParameter(0, p0);
  return f;
}

// ----------------------------------------------------------------------
TF1* initFunc::crystalBall(TH1 *h, double peak, double sigma, double alpha, double tailLength) {
  int npar(5);
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_xb", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_xb", fName.c_str()), iF_cb, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), npar);
  f->SetParNames("peak", "sigma", "crossover", "tail", "normalization");
  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  f->SetParameters(peak, sigma, alpha, tailLength, h->GetBinContent(h->FindBin(peak)));
  applyLimits(f, "crystalBall");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::pol1CrystalBall(TH1 *h, double peak, double sigma, double alpha, double tailLength) {
  int npar(7);
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_pol1xb", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_pol1xb", fName.c_str()), iF_pol1_cb, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), npar);
  f->SetParNames("peak", "sigma", "crossover", "tail", "norm.", "constant", "slope");

  double p0, p1;
  initPol1(p0, p1, h);
  double g0 = h->GetBinContent(h->FindBin(peak));
  cout << "peak: " << peak
       << " sigma: " << sigma
       << " alpha: " << alpha
       << " tail: " << tailLength
       << " norm: " << g0
       << " p0: " << p0
       << " p1: " << p1
       << endl;
  f->SetParameters(peak, sigma, alpha, tailLength, g0, p0, p1);
  applyLimits(f, "pol1CrystalBall");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::crystalBallGauss(TH1 *h, double peak, double sigma, double alpha, double tailLength) {
  int npar(8);
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_xbgauss", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_xbgauss", fName.c_str()), iF_gauss_cb, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), npar);
  f->SetParNames("fracNormG", "peakG", "fracsigmaG", "peak", "sigma", "crossover", "tail", "norm.");

  double g0 = h->GetBinContent(h->FindBin(peak));
  double fracNormG = 0.05;
  double fracPeakG = 1.;
  double fracSigmaG = 2.;
  cout << "fracNormG: " << fracNormG
       << " peakG: " << fracPeakG
       << " fracsigmaG: " << fracSigmaG
       << " peak: " << peak
       << " sigma: " << sigma
       << " alpha: " << alpha
       << " tail: " << tailLength
       << " norm: " << g0
       << endl;
  f->SetParameters(fracNormG, fracPeakG, fracSigmaG, peak, sigma, alpha, tailLength, g0);
  applyLimits(f, "crystalBallGauss");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::crystalBallBiGauss(TH1 *h, double peak, double sigma, double alpha, double tailLength) {
  int npar(9);
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_xbbigauss", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_xbbigauss", fName.c_str()), iF_bigauss_cb, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), npar);
  f->SetParNames("fracNormG", "peakG", "fracsigGl", "fracsigGh", "peak", "sigma", "crossover", "tail", "norm.");

  double g0 = h->GetBinContent(h->FindBin(peak));
  double fracNormG = 0.05;
  double fracPeakG = 1.;
  double fracSigmaGl = 2.;
  double fracSigmaGh = 2.;
  cout << "fracNormG: " << fracNormG
       << " peakG: " << fracPeakG
       << " fracsigmaGl: " << fracSigmaGl
       << " fracsigmaGh: " << fracSigmaGh
       << " peak: " << peak
       << " sigma: " << sigma
       << " alpha: " << alpha
       << " tail: " << tailLength
       << " norm: " << g0
       << endl;
  f->SetParameters(fracNormG, fracPeakG, fracSigmaGl, fracSigmaGh, peak, sigma, alpha, tailLength, g0);
  applyLimits(f, "crystalBallBiGauss");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::pol0gauss(TH1 *h, double peak, double sigma) {
  int npar(4);
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_pol0gauss", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_pol0gauss", fName.c_str()), iF_pol0_gauss, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), npar);
  f->SetParNames("normalization", "peak", "sigma", "constant");

  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  double p0;
  initPol0(p0, h);
  double A  = p0 * (fHi - fLo);
  double g0 = (h->Integral(lbin, hbin) - A/h->GetBinWidth(1))*h->GetBinWidth(1);

  f->SetParameters(g0, peak, sigma, p0);
  applyLimits(f, "pol0gauss");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::pol0Gauss(TH1 *h, double peak, double sigma) {
  int npar(4);
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_pol0Gauss", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_pol0Gauss", fName.c_str()), iF_pol0_Gauss, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), npar);
  f->SetParNames("area", "peak", "sigma", "constant", "slope");

  int lbin(1), hbin(h->GetNbinsX());
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  cout << "fLo: " << fLo << " fHi: " << fHi << " lbin: " << lbin << " hbin: " << hbin << endl;

  double p0;
  initPol0(p0, h);
  double A  = p0*(fHi - fLo);
  double g0 = (h->Integral(lbin, hbin)*h->GetBinWidth(1) - A);

  f->SetParameters(g0, peak, sigma, p0);
  applyLimits(f, "pol0Gauss");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::pol1gauss(TH1 *h, double peak, double sigma) {
  int npar(5);
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_pol1gauss", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_pol1gauss", fName.c_str()), iF_pol1_gauss, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), npar);
  f->SetParNames("gconst", "peak", "sigma", "p1const", "p1slope");

  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  double p0, p1;
  initPol1(p0, p1, h);
  double A   = 0.5*p1*(fHi*fHi - fLo*fLo) + p0*(fHi - fLo);
  double g0 = (h->Integral(lbin, hbin) - A/h->GetBinWidth(1))*h->GetBinWidth(1)/TMath::Sqrt(2*TMath::Pi())/sigma;

  cout << "fLo: " << fLo << " fHi: " << fHi << " lbin: " << lbin << " hbin: " << hbin << endl;
  cout << "g0: " << g0 << " peak: " << peak << " sigma: " << sigma << " p0: " << p0 << " p1: " << p1 << endl;
  f->SetParameters(g0, peak, sigma, p0, p1);
  applyLimits(f, "pol1gauss");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::pol1Gauss(TH1 *h, double peak, double sigma) {
  int npar(5);
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_pol1Gauss", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_pol1Gauss", fName.c_str()), iF_pol1_Gauss, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), npar);
  f->SetParNames("area", "peak", "sigma", "constant", "slope");

  int lbin(1), hbin(h->GetNbinsX());
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }


  double p0, p1;
  initPol1(p0, p1, h);
  double A   = 0.5*p1*(fHi*fHi - fLo*fLo) + p0*(fHi - fLo);
  double g0 = (h->Integral(lbin, hbin)*h->GetBinWidth(1) - A);
  cout << "fLo: " << fLo << " fHi: " << fHi << " lbin: " << lbin << " hbin: " << hbin << endl;
  cout << "g0: " << g0 << " peak: " << peak << " sigma: " << sigma << " p0: " << p0 << " p1: " << p1 << endl;
  f->SetParameters(g0, peak, sigma, p0, p1);
  applyLimits(f, "pol1Gauss");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::pol1Landau(TH1 *h, double peak, double sigma) {
  int npar(5);
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_pol1Landau", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_pol1Landau", fName.c_str()), iF_pol1_landau, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), npar);
  f->SetParNames("mpv", "sigma", "norm", "constant", "slope");

  int lbin(1), hbin(h->GetNbinsX());
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }


  double p0, p1;
  initPol1(p0, p1, h);
  double A   = 0.5*p1*(fHi*fHi - fLo*fLo) + p0*(fHi - fLo);
  double g0 = (h->Integral(lbin, hbin)*h->GetBinWidth(1) - A);
  cout << "fLo: " << fLo << " fHi: " << fHi << " lbin: " << lbin << " hbin: " << hbin << endl;
  cout << "g0: " << g0 << " peak: " << peak << " sigma: " << sigma << " p0: " << p0 << " p1: " << p1 << endl;
  f->SetParameters(peak, sigma, g0, p0, p1);
  applyLimits(f, "pol1Landau");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::pol2local(TH1 *h, double width) {
  TF1 *f(0);

  int lbin(1), hbin(h->GetNbinsX());
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  //  while ((f = (TF1*)gROOT->FindObject(Form("%s_pol2local", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_pol2local", fName.c_str()), iF_pol2local, h->GetBinLowEdge(lbin), h->GetBinLowEdge(hbin+1), 3);
  f->SetParNames("p0", "p1", "max");

  double maxVal   = h->GetMaximum();
  double maxPlace = h->GetBinCenter(h->GetMaximumBin());
  int    ytopbin  = h->FindBin(maxPlace + width);
  double ytop     = h->GetBinContent(ytopbin);
  double xtop     = h->GetBinCenter(ytopbin);
  double slope    = (ytop-maxVal)/(xtop-maxPlace)/(xtop-maxPlace);
  f->SetParameters(maxVal, slope, maxPlace);
  f->SetLineWidth(2);
  return f;
}



// ----------------------------------------------------------------------
TF1* initFunc::expoGauss(TH1 *h, double peak, double sigma) {

  int npar(5);
  TF1 *f = (TF1*)gROOT->FindObject(Form("%s_expo_Gauss", fName.c_str()));
  if (f) delete f;
  f = new TF1(Form("%s_expo_Gauss", fName.c_str()), iF_expo_Gauss, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), npar);
  f->SetParNames("area", "peak", "sigma", "base", "exp");
  //  f->SetLineColor(kBlue);
  f->SetLineWidth(2);

  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  double p0(0.), p1(0.);
  initExpo(p0, p1, h);
  double A   = p0*(TMath::Exp(p1*fHi) - TMath::Exp(p1*fLo));

  double g0 = (h->Integral(lbin, hbin)*h->GetBinWidth(1) - A);
  g0 = 10.;

  if (fVerbose) cout << "A: " << A << " g0: " << g0
		     << " peak = " << peak << " sigma = " << sigma
		     << " p0: " << p0 << " p1: " << p1
		     << endl;

  f->SetParameters(g0, peak, sigma, p0, p1);
  applyLimits(f, "expoGauss");
  return f;
}



// ----------------------------------------------------------------------
TF1* initFunc::expoErrGauss(TH1 *h, double peak, double sigma, double preco) {

  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_expo_err_Gauss", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_expo_err_Gauss", fName.c_str()), iF_expo_err_Gauss, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 9);
  f->SetParNames("area", "peak", "sigma", "base", "exp", "err0", "err1", "err2");
  //  f->SetLineColor(kBlue);
  f->SetLineWidth(2);

  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  double p0, p1;
  initExpo(p0, p1, h);
  double A   = p0*(TMath::Exp(p1*fHi) - TMath::Exp(p1*fLo));

  double g0 = (h->Integral(lbin, hbin)*h->GetBinWidth(1) - A);

  double e0(preco), e0Min(preco-0.001), e0Max(preco+0.001); // step
  double e1(0.075), e1Min(0.050),       e1Max(0.100);       // resolution
  //  double e2(1.15),  e2Min(1.05),        e2Max(1.25);        // level
  double e2(0.),     e2Min(0.),        e2Max(0.1);        // level


  if (fVerbose) cout << "A: " << A << " g0: " << g0
		     << " e0: " << e0 << " e1: " << e1 << " e2: " << e2
		     << " p0: " << p0 << " p1: " << p1
		     << endl;

  f->SetParameters(g0, peak, sigma, p0, p1, e0, e1, e2, 0.05*g0);

  // -- FIXME: remove hard-coded limits!
  f->ReleaseParameter(0);     f->SetParLimits(0, 0., 1.e6);
  f->ReleaseParameter(1);     f->SetParLimits(1, 5.2, 5.45);
  f->ReleaseParameter(2);     f->SetParLimits(2, 0.3*sigma, 1.3*sigma);
  f->ReleaseParameter(3);
  f->ReleaseParameter(4);
  f->ReleaseParameter(5);     f->SetParLimits(5, e0Min, e0Max);
  f->ReleaseParameter(6);     f->SetParLimits(6, e1Min, e1Max);
  f->ReleaseParameter(7);     f->SetParLimits(7, e2Min, e2Max);
  f->ReleaseParameter(8);     //f->SetParLimits(8, 0, 0.05*g0);

  return f;

}


// ----------------------------------------------------------------------
TF1* initFunc::pol1Err2Gauss(double lo, double hi) {
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_pol1_err2_Gauss", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_pol1_err2_Gauss", fName.c_str()), iF_pol1_err2_Gauss, lo, hi, 8);
  f->SetParNames("area", "peak", "sigma", "const", "slope", "err0", "err1", "err2");
  //  f->SetLineColor(kBlue);
  f->SetLineWidth(2);
  return f;
}

// ----------------------------------------------------------------------
TF1* initFunc::pol1Err2Gauss(TH1 *h, double peak, double sigma, double preco) {
  cout << "FIXME: missing implementation!" << endl;
  return 0;
}

// ----------------------------------------------------------------------
TF1* initFunc::pol1Err2gauss2c(double lo, double hi) {
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_pol1_err2_gauss2c", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_pol1_err2_gauss2c", fName.c_str()), iF_pol1_err2_gauss2c, lo, hi, 10);
  f->SetParNames("norm", "peak", "sigma", "fraction2", "sigma2", "const", "slope", "step", "res", "level");
  //  f->SetLineColor(kBlue);
  f->SetLineWidth(2);
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::expoErr2Gauss(double lo, double hi) {
  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_expo_err2_Gauss", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_expo_err2_Gauss", fName.c_str()), iF_expo_err2_Gauss, lo, hi, 8);
  f->SetParNames("area", "peak", "sigma", "base", "exp", "err0", "err1", "err2");
  //  f->SetLineColor(kBlue);
  f->SetLineWidth(2);
  return f;
}

// ----------------------------------------------------------------------
TF1* initFunc::expoErr2Gauss(TH1 *h, double peak, double sigma, double preco) {

  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_expo_err2_Gauss", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_expo_err2_Gauss", fName.c_str()), iF_expo_err2_Gauss, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 8);
  f->SetParNames("area", "peak", "sigma", "base", "exp", "err0", "err1", "err2");
  //  f->SetLineColor(kBlue);
  f->SetLineWidth(2);

  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  double p0, p1;
  initExpo(p0, p1, h);
  double A   = p0*(TMath::Exp(p1*fHi) - TMath::Exp(p1*fLo));
  double g0 = (h->Integral(lbin, hbin)*h->GetBinWidth(1) - A);

  double e0(preco), e0Min(preco-0.001), e0Max(preco+0.001); // step
  double e1(0.075), e1Min(0.050),       e1Max(0.100);       // resolution
  double e2(1.),     e2Min(0.99),        e2Max(1.01);       // high-range level = 0


  if (fVerbose) cout << "A: " << A << " g0: " << g0
		     << " e0: " << e0 << " e1: " << e1 << " e2: " << e2
		     << " p0: " << p0 << " p1: " << p1
		     << endl;

  f->SetParameters(g0, peak, sigma, p0, p1, e0, e1, e2, 0.05*g0);

  // -- FIXME: remove hard-coded limits!
  f->ReleaseParameter(0);     f->SetParLimits(0, 0., 1.e6);
  f->ReleaseParameter(1);     f->SetParLimits(1, 5.2, 5.45);
  f->ReleaseParameter(2);     f->SetParLimits(2, 0.3*sigma, 1.3*sigma);
  f->ReleaseParameter(3);
  f->ReleaseParameter(4);
  f->ReleaseParameter(5);     f->SetParLimits(5, e0Min, e0Max);
  f->ReleaseParameter(6);     f->SetParLimits(6, e1Min, e1Max);
  f->ReleaseParameter(7);     f->SetParLimits(7, e2Min, e2Max);
  f->ReleaseParameter(8);     //f->SetParLimits(8, 0, 0.05*g0);

  return f;

}


// ----------------------------------------------------------------------
TF1* initFunc::expoErrGaussLandau(TH1 *h, double peak, double sigma, double preco) {

  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_expo_err_Gauss_landau", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_expo_err_Gauss_landau", fName.c_str()), iF_expo_err_Gauss_landau, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 12);
  f->SetParNames("area", "peak", "sigma", "base", "exp", "err0", "err1", "err2", "err3");
  f->SetParName(9, "mpvl");
  f->SetParName(10, "sigl");
  //f->SetParName(11, "consl");
  //  f->SetLineColor(kBlue);
  f->SetLineWidth(2);

  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  double p0, p1;
  initExpo(p0, p1, h);
  if (p0 > 1.e7) p0 = 1.e7;
  double A   = p0*(TMath::Exp(p1*fHi) - TMath::Exp(p1*fLo));

  double g0 = (h->Integral(lbin, hbin)*h->GetBinWidth(1) - A);

  double e0(preco),  e0Min(preco-0.001), e0Max(preco+0.001);
  double e1(0.075),  e1Min(0.050), e1Max(0.100);
  double e2(1.15), e2Min(1.05),  e2Max(1.25);

  if (fVerbose) cout << "A: " << A << " g0: " << g0
		     << " e0: " << e0 << " e1: " << e1 << " e2: " << e2
		     << " p0: " << p0 << " p1: " << p1
		     << endl;

  f->SetParameters(g0, peak, sigma, p0, p1, e0, e1, e2, 0.05*g0);

  // -- FIXME: remove hard-coded limits!
  // for the landau
  // Fix the mpv and sigma to the fit values from the Bu2JpsiPi for the ENDCAP channel
  f->FixParameter(9, 5.349);   // landau mpv
  f->FixParameter(10, 0.06454);  // landau sigma
  //f->SetParameter(11, 1.);    // landau const (NOT A PARAMETER, FIXED BY THE BRANCHING RATIO)

  f->ReleaseParameter(0);     f->SetParLimits(0, 0., 1.e7);
  f->ReleaseParameter(1);     f->SetParLimits(1, 5.2, 5.45);
  f->ReleaseParameter(2);     f->SetParLimits(2, 0.3*sigma, 1.3*sigma);
  f->ReleaseParameter(3);
  f->ReleaseParameter(4);
  f->ReleaseParameter(5);     f->SetParLimits(5, e0Min, e0Max);
  f->ReleaseParameter(6);     f->SetParLimits(6, e1Min, e1Max);
  f->ReleaseParameter(7);     f->SetParLimits(7, e2Min, e2Max);
  f->ReleaseParameter(8);     //f->SetParLimits(8, 0, 0.05*g0);

  return f;

}



// ----------------------------------------------------------------------
TF1* initFunc::expoErrgauss2c(TH1 *h, double peak, double sigma1, double sigma2, double preco) {

  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_expo_err_gauss2c", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_expo_err_gauss2c", fName.c_str()), iF_expo_err_gauss2c, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 11);
  f->SetParNames("const", "peak", "sigma", "f2ndG", "s2ndG", "base", "exp", "err0", "err1", "err2", "err3");

  f->SetLineWidth(2);

  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  double p0, p1;
  initExpo(p0, p1, h);
  if (p0 > 1.e7) p0 = 1.e7;

  double A = p0*(TMath::Exp(p1*fHi) - TMath::Exp(p1*fLo));
  double H = h->Integral(lbin, hbin)*h->GetBinWidth(1);
  double g0 = (H - A);

  // -- new version with values from DK 2012/01/26: 5.146,0.055,1.00
  double e0(preco),  e0Min(preco-0.010), e0Max(preco+0.010);
  double e1(0.055),  e1Min(0.045), e1Max(0.065);
  double e2(1.00), e2Min(0.8),  e2Max(1.2);

  if (fVerbose) {
    cout << "fLo = " << fLo << " fHi = " << fHi << endl;
    cout << "A: " << A << " g0: " << g0 << " e0: " << e0 << " e1: " << e1 << " e2: " << e2 << " p0: " << p0 << " p1: " << p1
	 << " H: " << H << endl;
  }
  f->SetParameters(g0, peak, sigma1, 0.2, sigma2, p0, p1, e0, e1, e2, 0.05*g0);

  // -- FIXME: remove hard-coded limits!?
  if (fLimit[0]) {
    cout << "initFunc::expoErrgauss2c> limiting par 0 from " << fLimitLo[0] << " .. " << fLimitHi[0] << endl;
    f->SetParLimits(0, fLimitLo[0], fLimitHi[0]);
  } else {
    f->SetParLimits(0, 0., 1.e7);
  }

  if (fLimit[1]) {
    cout << "initFunc::expoErrgauss2c> limiting par 1 from " << fLimitLo[1] << " .. " << fLimitHi[1] << endl;
    f->SetParLimits(1, fLimitLo[1], fLimitHi[1]);
  } else {
    f->SetParLimits(1, 5.2, 5.45);
  }

  if (fLimit[2]) {
    cout << "initFunc::expoErrgauss2c> limiting par 2 from " << fLimitLo[2] << " .. " << fLimitHi[2] << endl;
    f->SetParLimits(2, fLimitLo[2], fLimitHi[2]);
  } else {
    f->SetParLimits(2, 0.2*sigma1, 1.5*sigma1);
  }

  if (fLimit[3]) {
    cout << "initFunc::expoErrgauss2c> limiting par 3 from " << fLimitLo[3] << " .. " << fLimitHi[3] << endl;
    f->SetParLimits(3, fLimitLo[3], fLimitHi[3]);
  } else {
    //    f->SetParLimits(2, 0.2*sigma1, 1.5*sigma1);
  }

  if (fLimit[4]) {
    cout << "initFunc::expoErrgauss2c> limiting par 4 from " << fLimitLo[4] << " .. " << fLimitHi[4] << endl;
    f->SetParLimits(4, fLimitLo[4], fLimitHi[4]);
  } else {
    f->SetParLimits(4, 0.5*sigma2, 2.0*sigma2);
  }

  if (fLimit[5]) {
    cout << "initFunc::expoErrgauss2c> limiting par 5 from " << fLimitLo[5] << " .. " << fLimitHi[5] << endl;
    f->SetParLimits(5, fLimitLo[5], fLimitHi[5]);
  } else {
    //    f->SetParLimits(5, 0.5*sigma2, 2.0*sigma2);
  }

  if (fLimit[6]) {
    cout << "initFunc::expoErrgauss2c> limiting par 6 from " << fLimitLo[6] << " .. " << fLimitHi[6] << endl;
    f->SetParLimits(6, fLimitLo[6], fLimitHi[6]);
  } else {
    //    f->SetParLimits(6, 0.5*sigma2, 2.0*sigma2);
  }

  if (fLimit[7]) {
    cout << "initFunc::expoErrgauss2c> limiting par 7 from " << fLimitLo[7] << " .. " << fLimitHi[7] << endl;
    f->SetParLimits(7, fLimitLo[7], fLimitHi[7]);
  } else {
    f->SetParLimits(7, e0Min, e0Max);
  }

  if (fLimit[8]) {
    cout << "initFunc::expoErrgauss2c> limiting par 8 from " << fLimitLo[8] << " .. " << fLimitHi[8] << endl;
    f->SetParLimits(8, fLimitLo[8], fLimitHi[8]);
  } else {
    f->SetParLimits(8, e1Min, e1Max);
  }

  if (fLimit[9]) {
    cout << "initFunc::expoErrgauss2c> limiting par 9 from " << fLimitLo[9] << " .. " << fLimitHi[9] << endl;
    f->SetParLimits(9, fLimitLo[9], fLimitHi[9]);
  } else {
    f->SetParLimits(9, e2Min, e2Max);
  }

  if (fLimit[10]) {
    cout << "initFunc::expoErrgauss2c> limiting par 10 from " << fLimitLo[10] << " .. " << fLimitHi[10] << endl;
    f->SetParLimits(10, fLimitLo[10], fLimitHi[10]);
  } else {
    //    f->SetParLimits(10, e2Min, e2Max);
  }
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::expoErrgauss2(TH1 *h, double peak1, double sigma1, double peak2, double sigma2, double preco) {

  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_expo_err_gauss2", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_expo_err_gauss2", fName.c_str()), iF_expo_err_gauss2, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 12);
  f->SetParNames("const", "peak", "sigma", "f2ndG", "p2ndG", "s2ndG", "base", "exp", "err0", "err1", "err2");
  f->SetParName(11, "err3");

  f->SetLineWidth(2);

  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  double p0, p1;
  initExpo(p0, p1, h);
  if (p0 > 1.e7) p0 = 1.e7;

  double A   = p0*(TMath::Exp(p1*fHi) - TMath::Exp(p1*fLo));
  double H   = h->Integral(lbin, hbin)*h->GetBinWidth(1);

  double g0 = H - A;

  double e0(preco),  e0Min(preco-0.001), e0Max(preco+0.001);
  double e1(0.055),  e1Min(0.8*e1), e1Max(1.2*e1);
  double e2(1.00), e2Min(0.8*e2),  e2Max(1.2*e2);


  if (fVerbose)  cout << "A: " << A << " g0: " << g0
		      << " e0: " << e0 << " e1: " << e1 << " e2: " << e2
		      << " p0: " << p0 << " p1: " << p1
		      << endl;

  f->SetParameters(g0, peak1, sigma1, 0.2, peak2, sigma2, p0, p1, e0, e1, e2);
  f->SetParameter(11,  0.05*g0);

  // -- FIXME: remove hardcoded par limits!
  f->ReleaseParameter(0);     f->SetParLimits(0, 0., 1.e7);
  f->ReleaseParameter(1);     f->SetParLimits(1, 5.2, 5.45);
  f->ReleaseParameter(2);     f->SetParLimits(2, 0.2*sigma1, 1.5*sigma1);
  f->ReleaseParameter(3);     f->SetParLimits(3, 0., 10000.);
  f->ReleaseParameter(4);     f->SetParLimits(4, 5.2, 5.45);
  f->ReleaseParameter(5);     f->SetParLimits(5, 0.5*sigma2, 2.0*sigma2);
  f->ReleaseParameter(6);
  f->ReleaseParameter(7);
  f->ReleaseParameter(8);     f->SetParLimits(8, e0Min, e0Max);
  f->ReleaseParameter(9);     f->SetParLimits(9, e1Min, e1Max);
  f->ReleaseParameter(10);     f->SetParLimits(10, e2Min, e2Max);
  f->ReleaseParameter(11);     //f->SetParLimits(8, 0, 0.05*g0);
  return f;

}

// ----------------------------------------------------------------------
// 2 gaussian, 2nd one has a fixed shape and the area related to the area of the 1st gaussian
// fraction:
// -1 fraction of the 2nd gaussian is left free to be fitted,
// >=0. the fraction is fixed to the passed value
// preco:
// -1 the step error function is not used in the fit
// >=0. the step is used with the edge equal to preco
//
TF1* initFunc::expoErrgauss2f(TH1 *h, double peak1, double sigma1, double peak2, double sigma2, double fraction, double preco) {
  bool free2ndGauss = true, freeErr=false;
  if(preco < 0.) {freeErr = false; preco=0.;} // disable the edge error function fitting
  else           {freeErr = true;} // use the edge fit at edge=preco
  if(fraction < 0.) {free2ndGauss = true;} // float the area of the 2nd gaus
  else              {free2ndGauss = false;} // fix the area to fraction

  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_expo_err_gauss2f", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_expo_err_gauss2f", fName.c_str()), iF_expo_err_gauss2f, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 12);
  f->SetParNames("const", "peak", "sigma", "fraction","p2ndG", "s2ndG", "base", "exp", "err0", "err1", "err2");
  f->SetParName(11, "err3");

  f->SetLineWidth(2);

  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  double p0(0.), p1(0.);
  initExpo(p0, p1, h);
  double A   = p0*(TMath::Exp(p1*fHi) - TMath::Exp(p1*fLo));
  double g0 = (h->Integral(lbin, hbin)*h->GetBinWidth(1) - A);
  g0 = 10.;

  double e0(preco), e0Min(preco-0.001), e0Max(preco+0.001);
  double e1(0.050), e1Min(0.8*e1),      e1Max(1.2*e1);
  double e2(1.00),  e2Min(0.8*e2),      e2Max(1.2*e2);

  if (fVerbose) cout << "A: " << A << " g0: " << g0
		     << " e0: " << e0 << " e1: " << e1 << " e2: " << e2
		     << " p0: " << p0 << " p1: " << p1
		     << endl;

  f->SetParameters(g0, peak1, sigma1, 0.1, peak2, sigma2, p0, p1, e0, e1, e2);
  f->SetParameter(11,  0.05*g0);

  if(free2ndGauss) {
    f->ReleaseParameter(3);     f->SetParLimits(3, 0., 1.);
  } else {
    f->FixParameter(3, fraction);   // fix the 2nd gaus fraction
  }

  // for the 2nd gaussian
  // Fix the 2nd gaussian to the Kstar shape extracted from B0ToJpsiX
  f->FixParameter(4, peak2);   // 2nd gaus mean
  f->FixParameter(5, sigma2);  // 2nd gaus sigma

  f->ReleaseParameter(6);     // exp
  f->ReleaseParameter(7);     // exp
  f->ReleaseParameter(8);     f->SetParLimits(8, e0Min, e0Max);  // err()
  f->ReleaseParameter(9);     f->SetParLimits(9, e1Min, e1Max);
  f->ReleaseParameter(10);     f->SetParLimits(10, e2Min, e2Max);
  if(freeErr) {
    f->ReleaseParameter(8);     f->SetParLimits(8, e0Min, e0Max);  // err()
    f->ReleaseParameter(9);     f->SetParLimits(9, e1Min, e1Max);
    f->ReleaseParameter(10);     f->SetParLimits(10, e2Min, e2Max);
    f->ReleaseParameter(11);   f->SetParLimits(11, 0., 1000.);
  } else {
    f->FixParameter(8, e0);
    f->FixParameter(9, e1);
    f->FixParameter(10, e2);
    f->FixParameter(11, 0.0);
  }  // fix at 0

  return f;

}


// ----------------------------------------------------------------------
// Add a landau to expErrGauss2
// The position (mpv) and sigma of the landau are fixed to the values obtained from the Bu2JpsiPi fit
TF1* initFunc::expoErrgauss2Landau(TH1 *h, double peak1, double sigma1, double peak2, double sigma2, double preco) {

  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_expo_err_gauss2_landau", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_expo_err_gauss2_landau", fName.c_str()), iF_expo_err_gauss2_landau, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 14);
  f->SetParNames("const", "peak", "sigma", "f2ndG", "p2ndG", "s2ndG", "base", "exp", "err0", "err1", "err2");
  f->SetParName(11, "err3");
  f->SetParName(12, "mpvl");
  f->SetParName(13, "sigl");
  //f->SetParName(14, "constl");

  f->SetLineWidth(2);

  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  double p0, p1;
  initExpo(p0, p1, h);
  if (p0 > 1.e7) p0 = 1.e7;

  double A   = p0*(TMath::Exp(p1*fHi) - TMath::Exp(p1*fLo));
  double H   = h->Integral(lbin, hbin)*h->GetBinWidth(1);

  double g0 = H - A;

  double e0(preco),  e0Min(preco-0.001), e0Max(preco+0.001);
  double e1(0.055),  e1Min(0.8*e1), e1Max(1.2*e1);
  double e2(1.00), e2Min(0.8*e2),  e2Max(1.2*e2);


  cout << "A: " << A << " g0: " << g0 << " e0: " << e0 << " e1: " << e1 << " e2: " << e2 << " p0: " << p0 << " p1: " << p1 << endl;

  f->SetParameters(g0, peak1, sigma1, 0.2, peak2, sigma2, p0, p1, e0, e1, e2);
  f->SetParameter(11,  0.05*g0);

  // for the landau
  // Fix the mpv and sigma to the fit values from the Bu2JpisPi for the BARREL channel
  f->FixParameter(12, 5.357);   // landau mpv
  f->FixParameter(13, 0.04885);  // landau sigma
  //f->SetParameter(14, 1.);    // landau const (NOT A PARAMETER, FIXED BY THE BRANCHING RATIO)


  f->ReleaseParameter(0);     f->SetParLimits(0, 0., 1.e7);
  f->ReleaseParameter(1);     f->SetParLimits(1, 5.2, 5.45);
  f->ReleaseParameter(2);     f->SetParLimits(2, 0.2*sigma1, 1.5*sigma1);
  f->ReleaseParameter(3);     f->SetParLimits(3, 0., 10000.);
  f->ReleaseParameter(4);     f->SetParLimits(4, 5.2, 5.45);
  f->ReleaseParameter(5);     f->SetParLimits(5, 0.5*sigma2, 2.0*sigma2);
  f->ReleaseParameter(6);
  f->ReleaseParameter(7);
  f->ReleaseParameter(8);     f->SetParLimits(8, e0Min, e0Max);
  f->ReleaseParameter(9);     f->SetParLimits(9, e1Min, e1Max);
  f->ReleaseParameter(10);     f->SetParLimits(10, e2Min, e2Max);
  f->ReleaseParameter(11);     //f->SetParLimits(8, 0, 0.05*g0);

  return f;

}

// ----------------------------------------------------------------------
TF1* initFunc::pol1ErrGauss(TH1 *h, double peak, double sigma, double preco) {

  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_pol1_err_Gauss", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_pol1_err_Gauss", fName.c_str()), iF_pol1_err_Gauss, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 9);
  f->SetParNames("area", "peak", "sigma", "const", "slope", "err0", "err1", "err2", "err3");
  //  f->SetLineColor(kBlue);
  f->SetLineWidth(2);

  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  double p0, p1;
  initPol1(p0, p1, h);
  double A   = 0.5*p1*(fHi*fHi - fLo*fLo) + p0*(fHi - fLo);
  double g0 = (h->Integral(lbin, hbin)*h->GetBinWidth(1) - A);

  double e0(preco),  e0Min(preco-0.01), e0Max(preco+0.01);
  double e1(0.07),  e1Min(0.04), e1Max(0.09);
  double e2(1.112), e2Min(0.900),  e2Max(1.200);

  cout << "A: " << A << " g0: " << g0 << " e0: " << e0 << " e1: " << e1 << " e2: " << e2 << " p0: " << p0 << " p1: " << p1 << endl;

  f->SetParameters(g0, peak, sigma, p0, p1, e0, e1, e2, 0.05*g0);

  f->ReleaseParameter(0);     f->SetParLimits(0, 0., 1.e7);
  f->ReleaseParameter(1);     f->SetParLimits(1, 5.2, 5.45);
  f->ReleaseParameter(2);     f->SetParLimits(2, 0.010, 0.080);
  f->ReleaseParameter(3);
  f->ReleaseParameter(4);
  f->ReleaseParameter(5);     f->SetParLimits(5, e0Min, e0Max);
  f->ReleaseParameter(6);     f->SetParLimits(6, e1Min, e1Max);
  f->ReleaseParameter(7);     f->SetParLimits(7, e2Min, e2Max);
  f->ReleaseParameter(8);     //f->SetParLimits(8, 0, 0.2*g0);

  return f;

}


// ----------------------------------------------------------------------
TF1* initFunc::pol1gauss2c(TH1 *h, double peak, double sigma) {

  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_pol1_gauss2c", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_pol1_gauss2c", fName.c_str()), iF_pol1_gauss2c, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 7);
  f->SetParNames("norm", "peak", "sigma", "fraction", "sigma2", "constant", "slope");
  //  f->SetLineColor(kBlue);
  f->SetLineWidth(2);

  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }
  double xlo(h->GetBinLowEdge(lbin));
  double xhi(h->GetBinLowEdge(hbin));

  double p0, p1;
  initPol1(p0, p1, h);
  if (fVerbose) cout << "p0: " << p0 << " p1: " << p1 << endl;
  double A   = 0.5*p1*(xhi*xhi - xlo*xlo) + p0*(xhi - xlo);

  double sqrt2pi = 2.506628275;
  double gInt    = h->Integral(lbin, hbin) - A;

  double gaussN  = gInt/(2.*sqrt2pi*sigma)*h->GetBinWidth(1);

  if (fVerbose) cout << "initFunc> gaussN = " << gaussN << " peak = " << peak << " sigma = " << sigma << " p0 = " << p0 << " p1 = " << p1 << endl;
  f->SetParameters(gaussN, peak, sigma, 0.2, 1.8*sigma, p0, p1);
  f->ReleaseParameter(0);     f->SetParLimits(0, 0., 1.e7);
  f->ReleaseParameter(1);     f->SetParLimits(1, peak-0.1, peak+0.1);
  f->ReleaseParameter(2);     f->SetParLimits(2, sigma*0.4, sigma*1.3);
  f->ReleaseParameter(3);     f->SetParLimits(3, 0.01, 2.0);
  f->ReleaseParameter(4);     f->SetParLimits(4, sigma*1.2, sigma*10.0);
  applyLimits(f, "pol1gauss2c");
  return f;


}

// ----------------------------------------------------------------------
TF1* initFunc::expogauss2(TH1 *h, double peak, double sigma, double deltaPeak, double deltaSigma) {

  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_expo_gauss2", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_expo_gauss2", fName.c_str()), iF_expo_gauss2, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 8);
  f->SetParNames("norm", "peak", "sigma", "fraction", "peak2", "sigma2", "base", "exp");
  f->SetLineWidth(2);

  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  double p0, p1;
  initExpo(p0, p1, h);
  cout << "p0: " << p0 << " p1: " << p1 << endl;

  double gaussN  = h->GetBinContent(h->FindBin(peak));

  f->SetParameters(gaussN, peak, sigma, 0.2, peak + deltaPeak, sigma+deltaSigma, p0, p1);
  f->ReleaseParameter(0);     f->SetParLimits(0, 0., 1.e7);
  f->ReleaseParameter(1);     f->SetParLimits(1, peak-0.1, peak+0.1);
  f->ReleaseParameter(2);     f->SetParLimits(2, 0.2*sigma, 1.6*sigma);
  f->ReleaseParameter(3);     f->SetParLimits(3, 0.01, 2.0);
  f->ReleaseParameter(4);     f->SetParLimits(4, peak+deltaPeak-0.2, peak+deltaPeak+0.2);
  f->ReleaseParameter(5);     f->SetParLimits(5, sigma*1.01, sigma*10.0);

  return f;
}



// ----------------------------------------------------------------------
TF1* initFunc::expogauss2c(TH1 *h, double peak, double sigma, double scaleSigma) {

  if (scaleSigma < 1.) {
    scaleSigma = 1.1;
    cout << "second gaussian sigma should be larger than for the first Gaussian, resetting scaleSigma = " << scaleSigma << endl;
  }

  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_expo_gauss2c", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_expo_gauss2c", fName.c_str()), iF_expo_gauss2c, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 7);
  f->SetParNames("norm", "peak", "sigma", "fraction", "sigma2", "base", "exp");
  f->SetLineWidth(2);

  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  double p0, p1;
  initExpo(p0, p1, h);
  cout << "p0: " << p0 << " p1: " << p1 << endl;

  double gaussN  = h->GetBinContent(h->FindBin(peak));

  f->SetParameters(gaussN, peak, sigma, 0.2, sigma*scaleSigma*1.02, p0, p1);
  f->ReleaseParameter(0);     f->SetParLimits(0, 0., 1.e7);
  f->ReleaseParameter(1);
  f->ReleaseParameter(2);     f->SetParLimits(2, 0.2*sigma, sigma*scaleSigma);
  f->ReleaseParameter(3);     f->SetParLimits(3, 0.01, 0.4);
  f->ReleaseParameter(4);     f->SetParLimits(4, sigma*scaleSigma*1.01, sigma*scaleSigma*1.4);
  f->ReleaseParameter(5);     f->SetParLimits(5, 0., 1.e9);
  f->ReleaseParameter(6);

  return f;


}


// ----------------------------------------------------------------------
TF1* initFunc::pol1gauss2(TH1 *h, double peak, double sigma, double deltaPeak, double deltaSigma) {

  TF1 *f(0);
  while ((f = (TF1*)gROOT->FindObject(Form("%s_pol1_gauss2", fName.c_str())))) if (f) delete f;
  f = new TF1(Form("%s_pol1_gauss2", fName.c_str()), iF_pol1_gauss2, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 8);
  f->SetParNames("1norm", "peak", "sigma", "fraction", "peak2", "sigma2", "constant", "slope");
  f->SetLineWidth(2);

  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  double p0, p1;
  initPol1(p0, p1, h);
  cout << "p0: " << p0 << " p1: " << p1 << endl;

  double gaussN  = h->GetBinContent(h->FindBin(peak));

  f->SetParameters(gaussN, peak, sigma, 0.2, peak + deltaPeak, sigma+deltaSigma, p0, p1);
  f->ReleaseParameter(0);     f->SetParLimits(0, 0., 1.e7);
  f->ReleaseParameter(1);     f->SetParLimits(1, peak-0.1, peak+0.1);
  f->ReleaseParameter(2);     f->SetParLimits(2, 0.2*sigma, 1.6*sigma);
  f->ReleaseParameter(3);     f->SetParLimits(3, 0.01, 2.0);
  f->ReleaseParameter(4);     f->SetParLimits(4, peak+deltaPeak-0.2, peak+deltaPeak+0.2);
  f->ReleaseParameter(5);     f->SetParLimits(5, sigma*1.01, sigma*10.0);

  return f;


}


// ----------------------------------------------------------------------
void initFunc::initPol0(double &p0, TH1 *h) {

  int EDG(4), NB(EDG+1);
  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  double ylo = h->Integral(lbin, lbin+EDG)/NB;
  double yhi = h->Integral(hbin-EDG, hbin)/NB;

  p0  = 0.5 * (ylo + yhi);
}

// ----------------------------------------------------------------------
void initFunc::initPol1(double &p0, double &p1, TH1 *h) {

  int EDG(4), NB(EDG+1);
  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }
  double xlo = h->GetBinLowEdge(lbin);

  double dx = h->GetBinLowEdge(hbin) - xlo;
  double ylo = h->Integral(lbin, lbin+EDG)/NB;
  double yhi = h->Integral(hbin-EDG-1, hbin-1)/NB;

  p1  = (yhi-ylo)/dx;
  p0  = ylo - p1*xlo;
  if (1 ||fVerbose) {
    cout << "lbin: " << lbin << " hbin: " << hbin
	 << " ylo: " << ylo << " yhi: " << yhi << " dx: " << dx
	 << " p0: " << p0 << " p1: " << p1 << endl;
  }
}


// ----------------------------------------------------------------------
void initFunc::initExpo(double &p0, double &p1, TH1 *h) {

  int EDG(4), NB(EDG+1);
  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  } else {
    fLo = h->GetBinLowEdge(lbin);
    fHi = h->GetBinLowEdge(hbin+1);
  }

  double dx = h->GetBinLowEdge(hbin+1) - h->GetBinLowEdge(lbin);
  double ylo = h->Integral(lbin, lbin+EDG)/NB;
  double yhi = h->Integral(hbin-EDG-1, hbin-1)/NB;

  if (ylo > 0 && yhi > 0) {
    p1 = (TMath::Log(yhi) - TMath::Log(ylo))/dx;
    p0 = ylo/TMath::Exp(p1*fLo);
  } else {
    if (yhi > ylo) {
      p1 = 1.;
    } else {
      p1 = -1.;
    }
    p0 = 50.;
  }

  // if (p0 > 1.e10) {
  //   p0 = 1.;
  //   p1 = -10.;
  // }

  if (p0 > 1.e9) p0 = 1.e9;

  if (fVerbose) {
    cout << "initFunc::initExpo dx: " << dx << endl;
    cout << "initFunc::initExpo fLo: " << fLo << " fHi: " << fHi << endl;
    cout << "initFunc::initExpo ylo: " << ylo << " yhi: " << yhi << endl;
    cout << "initFunc::initExpo  p0:  " << p0 <<  " p1:  " << p1 << endl;
    cout << "initFunc::initExpo integral: " <<  lbin << " .. " <<  lbin+EDG << endl;
    cout << "initFunc::initExpo integral: " <<	hbin-EDG-1 << " .. " << hbin-1 << endl;
  }
}


// ----------------------------------------------------------------------
void initFunc::initExpoHS(double &p0, double &p1, double &p2, TH1 *h) {

  int EDG(4), NB(EDG+1);
  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  } else {
    fLo = h->GetBinLowEdge(lbin);
    fHi = h->GetBinLowEdge(hbin+1);
  }
  p2   = fLo;

  double dx = h->GetBinLowEdge(hbin+1) - h->GetBinLowEdge(lbin);
  double ylo = h->Integral(lbin, lbin+EDG)/NB;
  double yhi = h->Integral(hbin-EDG-1, hbin-1)/NB;

  if (ylo > 0 && yhi > 0) {
    p1 = (TMath::Log(yhi) - TMath::Log(ylo))/dx;
    p0 = ylo;
  } else {
    if (yhi > ylo) {
      p1 = 1.;
    } else {
      p1 = -1.;
    }
    p0 = 50.;
  }

  if (fVerbose) {
    cout << "initFunc::initExpoHS lbin: " << lbin << " hbin: " << hbin << endl;
    cout << "initFunc::initExpoHS dx: " << dx << endl;
    cout << "initFunc::initExpoHS fLo: " << fLo << " fHi: " << fHi << endl;
    cout << "initFunc::initExpoHS ylo: " << ylo << " yhi: " << yhi << endl;
    cout << "initFunc::initExpoHS  p0:  " << p0 <<  " p1:  " << p1 <<  " p2:  " << p2 << endl;
    cout << "initFunc::initExpoHS integral: " <<  lbin << " .. " <<  lbin+EDG << endl;
    cout << "initFunc::initExpoHS integral: " <<	hbin-EDG-1 << " .. " << hbin-1 << endl;
    double par[] = {p0, p1, p2};
    double x[] = {fLo+0.001};
    cout << "initFunc::initExpoHS expoHS(" << x[0] << ") = " << iF_expo_HS(x, par) << endl;
  }
}

// ----------------------------------------------------------------------
// Uses the usual Landau from ROOT
TF1* initFunc::land(TH1 *h, double mpv, double sigma) {

  TF1 *f = new TF1(Form("%s_landau", fName.c_str()), iF_landau, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 3);
  //f->SetParNames("peak", "sigma", "constant");

  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  cout << "fLo: " << fLo << " fHi: " << fHi << " lbin: " << lbin << " hbin: " << hbin << endl;

  double p0 = mpv, p1 = sigma, p2=1.0;

  f->SetParameters(p0, p1, p2);
  f->ReleaseParameter(0);     f->SetParLimits(0, 5.0, 5.5);
  f->ReleaseParameter(1);     f->SetParLimits(1, 0., 0.1);
  f->ReleaseParameter(2);     f->SetParLimits(2, 0., 1e7);

  return f;

}


// ----------------------------------------------------------------------
// Simpler analytical "landau"
TF1* initFunc::landsimp(TH1 *h, double mpv, double sigma) {

  TF1 *f = new TF1(Form("%s_land", fName.c_str()), iF_landausimp, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), 3);
  //f->SetParNames("area", "peak", "sigma", "constant", "slope");

  int lbin(1), hbin(h->GetNbinsX()+1);
  if (fLo < fHi) {
    lbin = h->FindBin(fLo);
    hbin = h->FindBin(fHi);
  }

  cout << "fLo: " << fLo << " fHi: " << fHi << " lbin: " << lbin << " hbin: " << hbin << endl;

  double p0 = mpv, p1 = sigma, p2=1.0;

  f->SetParameters(p0, p1, p2);
  f->ReleaseParameter(0);     f->SetParLimits(0, 5.0, 5.5);
  f->ReleaseParameter(1);     f->SetParLimits(1, 0., 0.1);
  f->ReleaseParameter(2);     f->SetParLimits(2, 0., 1e7);
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::phiKK(TH1 *h) {
  int npar(8);
  TF1 *f = (TF1*)gROOT->FindObject(Form("%s_phiKK", fName.c_str()));
  if (f) delete f;
  f = new TF1(Form("%s_phiKK", fName.c_str()), iF_argus_gauss2, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), npar);
  f->SetParNames("const", "peak", "sigma1", "fraction", "sigma2", "norm.", "expo.", "endpoint");
  f->SetLineWidth(2);

  f->SetParameter(0, 11.23);
  f->SetParameter(1, 1.019);
  f->SetParameter(2, 0.003);
  f->SetParameter(3, 0.3);
  f->SetParameter(4, 0.010);
  f->SetParameter(5, 10.);
  f->SetParameter(6, 1.);
  fixPar(7, -2.*MKAON);

  applyLimits(f, "phiKK");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::bupsik(TH1 *h) {
  TVirtualFitter::SetMaxIterations(10000);
  TVirtualFitter::SetPrecision(1.e-4);
  int npar(10);
  TF1 *f = (TF1*)gROOT->FindObject(Form("%s_bupsik", fName.c_str()));
  if (f) delete f;
  f = new TF1(Form("%s_bupsik", fName.c_str()), iF_bupsik, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), npar);
  f->SetParNames("gnorm", "gpeak", "gsig1", "gfrac", "gsig2");

  f->SetParName(5, "exp0");
  f->SetParName(6, "exp1");

  f->SetParName(7, "errstep");
  f->SetParName(8, "errres");
  f->SetParName(9, "errnorm");

  f->SetLineWidth(2);

  f->SetParameter(0, h->GetMaximum());
  f->SetParameter(1, 5.28);
  f->SetParameter(2, 0.03); limitPar(2, 0.010, 0.040);
  f->SetParameter(3, 0.20); limitPar(3, 0.010, 0.400);
  f->SetParameter(4, 0.05); limitPar(4, 0.030, 0.100);

  double a(-1.), b(-1.);
  fLo = 5.5; fHi = 5.8;
  initExpo(a, b, h);
  f->SetParameter(5, a);
  f->SetParameter(6, b);

  //  fixPar(7, 5.142);
  f->SetParameter(7, 5.142); limitPar(7, 5.137, 5.147);
  f->SetParameter(8, 0.04);  limitPar(8, 0.030, 0.060);
  f->SetParameter(9, 0.4*h->Integral(1, 5)/5);
  applyLimits(f, "bupsik");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::bupsik1(TH1 *h) {
  TVirtualFitter::SetMaxIterations(10000);
  TVirtualFitter::SetPrecision(1.e-4);
  int npar(15);
  TF1 *f = (TF1*)gROOT->FindObject(Form("%s_bupsik1", fName.c_str()));
  if (f) delete f;
  f = new TF1(Form("%s_bupsik1", fName.c_str()), iF_bupsik1, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), npar);
  f->SetParNames("gnorm", "gpeak", "gsigL", "gsigH", "cbmean", "cbsigma", "cbalpha", "cbtail", "cbnorm", "satBg");
  f->SetParName(10, "exp0");
  f->SetParName(11, "exp1");

  f->SetParName(12, "errstep");
  f->SetParName(13, "errres");
  f->SetParName(14, "errnorm");

  f->SetLineWidth(2);

  f->SetParameter(0, 0.1); limitPar(0, 0.1,  0.4);
  f->SetParameter(1, 1.0); limitPar(1, 0.995, 1.005);
  f->SetParameter(2, 1.1); limitPar(2, 1.00, 2.00);
  f->SetParameter(3, 1.1); limitPar(3, 1.00, 2.00);

  f->SetParameter(4, 5.28);
  f->SetParameter(5, 0.03);
  fixPar(6, 4.); //  f->SetParameter(6, 1.6); //limitPar(6, 1.0, 5.0);
  fixPar(7, 2.); //f->SetParameter(7, 4.0); //limitPar(7, 2.0, 8.0);
  f->SetParameter(8, 0.7*h->GetMaximum());
  fixPar(9, 0.015);

  double a(-1.), b(-1.);
  fLo = 5.45; fHi = 5.9;
  initExpo(a, b, h);
  f->SetParameter(10, a);
  f->SetParameter(11, b);

  fixPar(12, 5.142);
  f->SetParameter(13, 0.05);  limitPar(13, 0.040, 0.060);
  //  cout << "set err norm: " << 0.5*h->Integral(1, 5)/5 << endl;
  f->SetParameter(14, 0.5*h->Integral(1, 5)/5);
  applyLimits(f, "bupsik1");
  return f;
}


// ----------------------------------------------------------------------
TF1* initFunc::bspsiphi(TH1 *h, double sigma) {
  TVirtualFitter::SetMaxIterations(10000);
  TVirtualFitter::SetPrecision(1.e-4);
  int npar(9);
  TF1 *f = (TF1*)gROOT->FindObject(Form("%s_bspsiphi", fName.c_str()));
  if (f) delete f;
  f = new TF1(Form("%s_bspsiphi", fName.c_str()), iF_bspsiphi, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1), npar);
  f->SetParNames("gnorm", "gpeak", "gsig1", "gfrac", "gsig2", "satnorm");

  f->SetParName(6, "exp0");
  f->SetParName(7, "exp1");
  f->SetParName(8, "cutoff");

  f->SetLineWidth(2);

  f->SetParameter(0, h->GetMaximum());
  f->SetParameter(1, 5.37);
  f->SetParameter(2, sigma); limitPar(2, 0.5*sigma, 1.4*sigma);
  f->SetParameter(3, 0.20); limitPar(3, 0.010, 0.600);
  f->SetParameter(4, 2.*sigma); limitPar(4, 1.5*sigma, 4.*sigma);
  f->SetParameter(5, 0.005*f->GetParameter(0)); limitPar(5, 0.001*f->GetParameter(0), 0.02*f->GetParameter(0));

  double a(-1.), b(-1.), c(-1.);
  fLo = 5.2; fHi = 5.8;
  initExpoHS(a, b, c, h);
  f->SetParameter(6, a);
  f->SetParameter(7, b);
  // limitPar(6, 0.0, 1.e5);
  // limitPar(7, -10., 10.);
  fixPar(8, c);
  applyLimits(f, "bspsiphi");
  return f;
}



// ----------------------------------------------------------------------=
double initFunc::iF_int_pisat(double norm) {
  double g3par[] = {norm, 5.36416e+00, 4.83167e-02, 3.39117e-01, 5.46405e+00,
		    9.83832e-02, 9.36584e-02, 5.52549e+00, 2.69386e-01};
  double sqrt2pi = 2.506628275;
  double satintegral = sqrt2pi*g3par[0]*(g3par[2] + g3par[3]*g3par[5] + g3par[6]*g3par[8]);
  return satintegral;
}
