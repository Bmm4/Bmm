void redVsSmall(string redFile, string smallFile, string varRed = "bdt", string varSmall = "same", string cut = "nada") {
  TFile *rFile = TFile::Open(redFile.c_str());
  TTree *rt = (TTree*)rFile->Get("candAnaMuMu/events");

  TFile *sFile = TFile::Open(smallFile.c_str());
  TKey *key(0);
  TIter next(sFile->GetListOfKeys());
  string tname("");
  TTree *st;
  while ((key = (TKey*)next())) {
    if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TTree")) continue;
    tname = key->GetName();
    st = (TTree*)sFile->Get(tname.c_str());
  }
  if (0 == st) {
    cout << "no tree in small file found" << endl;
    return;
  }

  gStyle->SetOptTitle(false);
  gStyle->SetOptStat(false);

  if (varSmall == "same") varSmall = varRed;

  st->SetLineColor(kBlue);
  string cuts = "chan==0 && muid";
  if (cut != "nada") cuts += Form(" && %s", cut.c_str());
  int nst = st->Draw(varSmall.c_str(), cuts.c_str(), "hist");
  rt->SetLineColor(kRed);
  cuts = "chan==0 && hlt1 && tos && l1t && bdt>-1 && gmuid";
  if (cut != "nada") cuts += Form(" && %s", cut.c_str());
  int nrt = rt->Draw(varRed.c_str(), cuts.c_str(), "samehist");

  tl->SetTextColor(kBlack); tl->DrawLatexNDC(0.2, 0.91, tname.c_str());
  tl->SetTextColor(kBlue);  tl->DrawLatexNDC(0.56, 0.95, Form("small tree: %d", nst));
  tl->SetTextColor(kRed);   tl->DrawLatexNDC(0.56, 0.91, Form("red. tree:   %d", nrt));

}
