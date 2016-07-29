#ifndef CUTS_h
#define CUTS_h

#include <string>

struct cuts {
  int index;
  std::string xmlFile;
  double mBdLo, mBdHi, mBsLo, mBsHi;
  double etaMin, etaMax, pt;
  double m1pt, m2pt, metaMin, metaMax;
  double iso, chi2dof, alpha, fls3d, flsxy, docatrk;
  double closetrk, maxdoca;
  double pvip, pvips;
  double pvlip, pvlips, pv2lip, pv2lips;
  void dump() {
    std::cout << "chan: " << index
	      << Form(", %3.2f < mu(fwd) < %3.2f", metaMin, metaMax)
	      << Form(", pt m1/m2 > %3.1f/%3.1f", m1pt, m2pt)
	      << Form(", fls3d/xy > %3.1f/%3.1f", fls3d, flsxy)
	      << Form(", chi2dof < %3.1f, alpha < %4.3f", chi2dof, alpha)
	      << std::endl;
  }
};


#endif
