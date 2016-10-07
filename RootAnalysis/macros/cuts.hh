#ifndef CUTS_h
#define CUTS_h

#include <string>

struct cuts {
  int index;
  std::string xmlFile;
  std::vector<int> l1seeds;
  double mBdLo, mBdHi, mBsLo, mBsHi;
  double etaMin, etaMax, pt;
  double m1pt, m2pt, metaMin, metaMax;
  double iso, chi2dof, alpha, fls3d, flxyLo, flxyHi, flsxy, docatrk;
  double closetrk, maxdoca;
  double pvip, pvips;
  double pvlip, pvlips, pv2lip, pv2lips;
  void dump() {
    std::string sl1seeds("");
    for (unsigned int i = 0; i < l1seeds.size(); ++i) sl1seeds += Form("%d ", l1seeds[i]);
    std::cout << "chan: " << index
	      << Form(", %3.2f < mu(fwd) < %3.2f", metaMin, metaMax)
	      << ", l1seeds: " << sl1seeds
	      << Form(", pt m1/m2 > %3.1f/%3.1f", m1pt, m2pt)
	      << Form(", fls3d/xy > %3.1f/%3.1f", fls3d, flsxy)
	      << Form(", %3.1f < flxy < %3.1f", flxyLo, flxyHi)
	      << Form(", chi2dof < %3.1f, alpha < %4.3f", chi2dof, alpha)
	      << std::endl;
  }
};

#endif
