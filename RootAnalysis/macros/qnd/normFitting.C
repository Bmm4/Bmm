void normFitting() {

  int imode(0);

  TFile *f = TFile::Open("s01/plotResults.2012s01.root");
  string  name = Form("hW8Norm_bdt2012s01_bupsikData_chan0");
  bool ok = f->cd("bupsikData");

  fitPsYield fpy(name);
  fpy.fitBu2JpsiKp(5, "test", imode, 5.0, 5.8, 0.025);

}
