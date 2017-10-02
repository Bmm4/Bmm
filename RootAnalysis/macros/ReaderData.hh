#ifndef READERDATA_h
#define READERDATA_h

struct ReaderData {
  float pt, eta, m1eta, m2eta, m1pt, m2pt, m1phi, m2phi;
  float fls3d, fl3d, alpha, maxdoca, pvip, pvips, iso, docatrk, chi2dof;
  float closetrk;
  float closetrks1, closetrks2, closetrks3;
  float m1iso, m2iso, pvdchi2, othervtx;
  float pvlips, pvlip;
  float pv2lips, pv2lip;
  float m;
};

#endif
