// ----------------------------------------------------------------------
void sbs(string sample = "ad0bdt_bupsikData", string selection = "Cu", string file = "s01/plotReducedOverlays.2011s01.root") {
  TFile *f = TFile::Open(file.c_str());
  string hname(Form("%s_bdt", sample.c_str()));

  string sbsControlPlotsFileName = Form("sbsctrltest");

  AnalysisDistribution a(hname.c_str());
  a.fVerbose = 1;
  a.fControlPlotsFileName = sbsControlPlotsFileName;
  a.fDirectory = string("./sbsctrl");


  TH1D *h(0);
  if (string::npos != sample.find("bupsik")) {
    h = a.sbsDistributionExpoErrGauss(hname.c_str(), selection.c_str());
  }
  if (string::npos != sample.find("bspsiphi")) {
    h = a.sbsDistributionBs2JpsiPhi(hname.c_str(), selection.c_str());
  }
  if (h) {
    h->Draw();
  }
}
