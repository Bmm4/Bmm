// ----------------------------------------------------------------------
void replaceAll(string &sInput, const string &oldString, const string &newString) {
  string::size_type foundpos = sInput.find(oldString);
  while (foundpos != string::npos)  {
    sInput.replace(sInput.begin() + foundpos, sInput.begin() + foundpos + oldString.length(), newString);
    foundpos = sInput.find(oldString);
  }
}


// ----------------------------------------------------------------------
void plotMC() {

  TFile *fm = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-combined-BuToJpsiKp-2016BF-s01.root");

  TTree *tm = (TTree*)fm->Get("candAnaBu2JpsiK/events");
  string cuts = "hlt1 && tos && l1t && chan == 0 && bdt > 0.3 ";
  cuts += " && kpt > 4.0 && dcand ";
  cuts += " && (m1pt>4.) && (m2pt>4.) && (pvips < 5.) && (flsxy > 4.) && (maxdoca < 0.08) && (chi2/dof < 5.)";
  cuts += " && (alpha<0.2) && (fls3d>4.) && (TMath::Abs(pvips) < 4.) && (TMath::Abs(pvip) < 0.02)";
  cuts += " && ktrk>4";


  vector<string> vlist;
  vlist.push_back("ktrk");
  vlist.push_back("kvalhits");
  vlist.push_back("kpix");
  vlist.push_back("klayerswithhits");
  vlist.push_back("TMath::Log10(kdszE)");
  vlist.push_back("TMath::Log10(kd0E)");
  vlist.push_back("TMath::Log10(kdxyE)");
  vlist.push_back("kptE");
  vlist.push_back("kchi2");
  vlist.push_back("kalg");

  string name, hname;
  TH1D *h1(0);
  for (unsigned int iv = 0; iv < vlist.size(); ++iv) {
    name = vlist[iv];
    replaceAll(name, "TMath::Log10(", "");
    replaceAll(name, ")", "");
    if (vlist[iv] == "ktrk") h1 = new TH1D(Form("%s", vlist[iv].c_str()), "", 20, 0., 20.);
    if (vlist[iv] == "kvalhits") h1 = new TH1D(Form("%s", vlist[iv].c_str()), "", 40, 0., 40.);
    if (vlist[iv] == "kpix") h1 = new TH1D(Form("%s", vlist[iv].c_str()), "", 10, 0., 10.);
    if (vlist[iv] == "klayerswithhits") h1 = new TH1D(Form("%s", vlist[iv].c_str()), "", 25, 0., 25.);
    if (vlist[iv] == "kptE") h1 = new TH1D(Form("%s", vlist[iv].c_str()), "", 30, 0., 0.3);
    if (vlist[iv] == "kchi2") h1 = new TH1D(Form("%s", vlist[iv].c_str()), "", 40, 0., 40.);
    if (vlist[iv] == "kalg") h1 = new TH1D(Form("%s", vlist[iv].c_str()), "", 15, 0., 15.);

    if (name == "kdszE") h1 = new TH1D(Form("%s", name.c_str()), "", 40, -3., -1.);
    if (name == "kd0E") h1 = new TH1D(Form("%s", name.c_str()), "", 40, -3., -1.);
    if (name == "kdxyE") h1 = new TH1D(Form("%s", name.c_str()), "", 40, -3., -1.);

    hname = h1->GetName();
    tm->Draw(Form("%s >> %s", vlist[iv].c_str(), hname.c_str()), cuts.c_str());
    h1->GetXaxis()->SetTitle(vlist[iv].c_str());

    name += ".pdf";
    c0->SaveAs(name.c_str());
  }
}
