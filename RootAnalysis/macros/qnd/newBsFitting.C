void bla() {
  TFile *_file0 = TFile::Open("s01/plotReducedOverlays.2012s01.root");
  TH1D *h = (TH1D*)_file0->Get("ad0bdt_bspsiphiData_tauMassAo");
  fitPsYield a(0, 1);
 a.fit1_Bs2JpsiPhi(h, -1, ".");
}
