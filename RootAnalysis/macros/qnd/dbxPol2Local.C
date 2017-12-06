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
void dbxPol2Local() {
  TH1D *h1 = new TH1D("h1", "", 200, -0.5, 0.5);
  TF1 *f1 = new TF1("f1", iF_gauss, -0.5, 0.5, 3);
  f1->SetParameters(1., 0., 0.1);
  h1->FillRandom("f1", 10000);
  h1->Draw();

  initFunc a;
  a.fLo = -0.1;
  a.fHi =  0.1;
  TF1 *F1 = a.pol2local(h1, 0.05);
  F1->Draw("same");

  h1->Fit(F1, "R", "", -0.05, 0.05);

}
