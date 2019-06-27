#ifndef REDTREEDATA
#define REDTREEDATA

#define NTRGMAX 40

struct redTreeData {
  Long64_t run, evt;
  int ls, tm, pr, procid, pvn, rr, l1s, ps, chan;
  double corrW8;
  bool json, cb, tis, reftrg;
  bool hlt1, tos, l1t;
  int ntrgtos[NTRGMAX];
  bool dcand;
  double bdt, bdt2, pvw8, pvz;
  double rlumi;

  bool gmuid, gmugmid, gmutmid, gmumvaid, gmupt, gmueta, gtqual, gtpt, gteta;

  double pvlip, pvlips, pv2lip, pv2lips, pvip, pvips, pvip3d, pvips3d, dzmin, dz12;
  int pvntrk, pv2ntrk;

  int q, type;
  double pt, eta, phi, tau, taue, tauxy, tauxye, m, me, cm, m3, m4, cosa, alpha, iso;
  int isotrk, closetrk, closetrks1, closetrks2, closetrks3;
  double chi2, dof, chi2dof, pchi2dof, fls3d, fl3d, flxy, fl3dE, flsxy, docatrk, docatrkbdt, maxdoca, lip, lipE, tip, tipE;

  // -- opposite side
  double osiso, osreliso, osmpt, osmptrel, osmdr;

  int m1q, m2q;
  double m1pt, m1eta, m1phi, m1ip, m1chi2;
  double m2pt, m2eta, m2phi, m2ip, m2chi2;
  double muhelicity;
  double kpt, keta, kphi, kip, kchi2;
  double k1pt, k1eta, k1phi, k2pt, k2eta, k2phi;
  double pipt, pieta, piphi;

  bool m1id, m2id, m1gmid, m2gmid, m1mvaid, m2mvaid, m1rmvaid, m2rmvaid;
  double m1mvabdt, m2mvabdt;
  double m1rmvabdt, m2rmvabdt;
  double m1trigm, m2trigm;
  int m1gt, m2gt, k1gt, k2gt, kgt, pigt;
  int m1pix, m1bpix, m1bpixl1, m1pv, m1trk;
  int m2pix, m2bpix, m2bpixl1, m2pv, m2trk;
  int kpix, kbpix, ktrk;

  // -- START for tracking studies
  int m1tkqual, m1alg, m1valhits, m1layerswithhits;
  double m1dz, m1dzE, m1d0, m1d0E, m1dsz, m1dszE, m1dxy, m1dxyE, m1valhitfraction, m1ptE, m1etaE, m1phiE;
  int m2tkqual, m2alg, m2valhits, m2layerswithhits;
  double m2dz, m2dzE, m2d0, m2d0E, m2dsz, m2dszE, m2dxy, m2dxyE, m2valhitfraction, m2ptE, m2etaE, m2phiE;
  int ktkqual, kalg, kvalhits, klayerswithhits;
  double kdz, kdzE, kd0, kd0E, kdsz, kdszE, kdxy, kdxyE, kvalhitfraction, kptE, ketaE, kphiE;
  // -- END for tracking studies

  double mudist, mudeltar;
  double g1pt, g2pt, g3pt, g4pt, g1eta, g2eta, g3eta, g4eta, gtau, gfl3d;
  int    g1id, g2id;
  double t1pt, t1eta, t2pt, t2eta, t3pt, t3eta, t4pt, t4eta;

  // -- psi variables
  double mpsi;
  double psipt, psieta, psiphi;
  double psicosa, psiflsxy, psiprob, psimaxdoca;

  // -- other resonances
  double mkk, mkpi, mkpi1, mkpi2;
  double phipt, phieta, phiphi;
  double phidr, kstardr;
  bool   kstarfail;

  double md0, dm, ptd0;

  double hm1pt, hm1eta, hm1phi, hm2pt, hm2eta, hm2phi;

  double pvdchi2, m1iso, m1xpdist, m2iso, m2xpdist, othervtx;

};

#endif
