#ifndef CUTS_h
#define CUTS_h

#include <string>

struct cuts {
  int index;
  std::string bdtXml;
  double bdtCut, bdtMuPt;
  double mBdLo, mBdHi, mBsLo, mBsHi, mBuLo, mBuHi;
  double etaMin, etaMax, pt, phiMin, phiMax;
  double m1pt, m2pt, metaMin, metaMax, muonbdt;
  double iso, m1iso, m2iso, chi2dof, alpha, fls3d, flxyLo, flxyHi, flsxy, docatrk;
  double closetrk, closetrks1, closetrks2, closetrks3, maxdoca;
  double pvip, pvips;
  double pvlip, pvlips, pv2lip, pv2lips;
  void dump() {
    std::cout << "chan: " << index
	      << Form(", %3.2f < mu(fwd) < %3.2f", metaMin, metaMax)
	      << Form(", pt m1/m2 > %3.1f/%3.1f", m1pt, m2pt)
	      << Form(", fls3d/xy > %3.1f/%3.1f", fls3d, flsxy)
	      << Form(", %3.1f < flxy < %3.1f", flxyLo, flxyHi)
	      << Form(", chi2dof < %3.1f, alpha < %4.3f", chi2dof, alpha)
	      << std::endl;
  }
};

#endif
