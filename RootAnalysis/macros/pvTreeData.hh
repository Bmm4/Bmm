#ifndef PVTREEDATA
#define PVTREEDATA

struct pvTreeData {
  float pt, eta, phi, m, npv, m1pt, m2pt;
  int chan, idx1, idx2, idx3;
  bool hlt1;

  float lz1, lz2, mult1, mult2, prob1, prob2, chi1, chi2, dz12, dzmin;
  float gfl, gt, flsxy, flxy, fls3d, fl3d, fl1, fl2, fl3, t1, t2, t3;
  float gs, s1, s2, s3;

  float gx, gy, gz, p1x, p1y, p1z, p1d, p2x, p2y, p2z, p2d, p3x, p3y, p3z, sx, sy, sz, dsv;
  float d1, a1, d2, a2, d3, a3;
};

#endif
