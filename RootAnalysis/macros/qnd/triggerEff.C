// ----------------------------------------------------------------------
std::pair<double, double> triggerEff(string var = "hlt1", string muonselection = "(m1rmvabdt>-2.5 && m2rmvabdt>-2.5) && (m1rmvabdt>0.58 && m2rmvabdt>0.58)"
		, string file = "/scratch/ursl/bmm4/s01/bmm-mc-combined-BdToMuMu-2016BF-s01.root") {
  TFile *f = TFile::Open(file.c_str());
  TTree *t = (TTree*)f->Get("candAnaMuMu/events");

  string basecuts = "0==chan && bdt > 0.3";
  string allcuts = basecuts + " && " + muonselection;
  string passcuts = allcuts + " && " + var;

  cout << "allcuts:  " << allcuts << endl;
  cout << "passcuts: " << passcuts << endl;

  c0->Clear();
  c0->Divide(1,2);
  c0->cd(1);
  double pass = t->Draw("m", passcuts.c_str());
  c0->cd(2);
  double all  = t->Draw("m", allcuts.c_str());
  double eff(-1.), err(-1.);
  if (all > 0) {
    eff = pass/all;
    err = TMath::Sqrt(((pass+1)*(all-pass+1))/((all+3)*(all+2)*(all+2)));
  }

  cout << "pass/all = " << pass << "/" << all << " = " << eff << " +/- " << err << endl;
  pair<double, double> a;
  a.first = eff;
  a.second = err;
  return a;
}


// ----------------------------------------------------------------------
void all() {
  string fname("/scratch/ursl/bmm4/s01/bmm-mc-combined-BdToMuMu-2016BF-s01.root");
  string globalmuon("(m1rmvabdt>-2.5 && m2rmvabdt>-2.5)");
  string ams("(m1rmvabdt>-2.5 && m2rmvabdt>-2.5) && (m1rmvabdt<0.58 && m2rmvabdt<0.58)");
  string ams1("(m1rmvabdt>-2.5 && m2rmvabdt>-2.5) && (m1rmvabdt<0.58 || m2rmvabdt<0.58)");
  string muid("(m1rmvabdt>-2.5 && m2rmvabdt>-2.5) && (m1rmvabdt>0.58 && m2rmvabdt>0.58)");

  string exa1("(m1rmvabdt>-2.5 && m2rmvabdt>-2.5) && ((m1rmvabdt<0.58 && m2rmvabdt>0.58) || (m1rmvabdt>0.58 && m2rmvabdt<0.58))");
  string hipt("(m1rmvabdt>-2.5 && m2rmvabdt>-2.5) && ((m2rmvabdt<0.58 && m1rmvabdt>0.58))");
  string lopt("(m1rmvabdt>-2.5 && m2rmvabdt>-2.5) && ((m2rmvabdt>0.58 && m1rmvabdt<0.58))");


  std::pair<double, double> sgglobalmuon = triggerEff("hlt1", globalmuon, fname);
  c0->SaveAs("triggerEffSg-gm.pdf");

  std::pair<double, double> sgams = triggerEff("hlt1", ams, fname);
  c0->SaveAs("triggerEffSg-ams.pdf");

  std::pair<double, double> sgams1 = triggerEff("hlt1", ams1, fname);
  c0->SaveAs("triggerEffSg-ams1.pdf");

  std::pair<double, double> sgmuid = triggerEff("hlt1", muid, fname);
  c0->SaveAs("triggerEffSg-muid.pdf");

  std::pair<double, double> sgexa1 = triggerEff("hlt1", exa1, fname);
  c0->SaveAs("triggerEffSg-exa1.pdf");

  std::pair<double, double> sghipt = triggerEff("hlt1", hipt, fname);
  c0->SaveAs("triggerEffSg-hipt.pdf");

  std::pair<double, double> sglopt = triggerEff("hlt1", lopt, fname);
  c0->SaveAs("triggerEffSg-lopt.pdf");

  fname = "/scratch/ursl/bmm4/s01/bmm-rare-rhlt1.root";
  std::pair<double, double> bgglobalmuon = triggerEff("rhlt1", globalmuon, fname);
  c0->SaveAs("triggerEffBg-gm.pdf");

  std::pair<double, double> bgams = triggerEff("rhlt1", ams, fname);
  c0->SaveAs("triggerEffBg-ams.pdf");

  std::pair<double, double> bgams1 = triggerEff("rhlt1", ams1, fname);
  c0->SaveAs("triggerEffBg-ams1.pdf");

  std::pair<double, double> bgmuid = triggerEff("rhlt1", muid, fname);
  c0->SaveAs("triggerEffBg-muid.pdf");

  std::pair<double, double> bgexa1 = triggerEff("rhlt1", exa1, fname);
  c0->SaveAs("triggerEffBg-exa1.pdf");

  std::pair<double, double> bghipt = triggerEff("rhlt1", hipt, fname);
  c0->SaveAs("triggerEffBg-hipt.pdf");

  std::pair<double, double> bglopt = triggerEff("rhlt1", lopt, fname);
  c0->SaveAs("triggerEffBg-lopt.pdf");

  cout << "Selection       signal  background" << endl;
  cout << Form("both muons fail MVA muon ID         &$%4.3f\\pm%4.3f$ &$%4.3f\\pm%4.3f$", sgams.first, sgams.second, bgams.first, bgams.second) << endl;
  cout << Form("at least one muon fails MVA muon ID &$%4.3f\\pm%4.3f$ &$%4.3f\\pm%4.3f$", sgams1.first, sgams1.second, bgams1.first, bgams1.second) << endl;
  cout << Form("exactly one muon fails MVA muon ID  &$%4.3f\\pm%4.3f$ &$%4.3f\\pm%4.3f$", sgexa1.first, sgexa1.second, bgexa1.first, bgexa1.second) << endl;
  // cout << Form("hi-pt muon pass MVA muon ID         &$%4.3f\\pm%4.3f$ &$%4.3f\\pm%4.3f$", sghipt.first, sghipt.second, bghipt.first, bghipt.second) << endl;
  // cout << Form("lo-pt muon pass MVA muon ID         &$%4.3f\\pm%4.3f$ &$%4.3f\\pm%4.3f$", sglopt.first, sglopt.second, bglopt.first, bglopt.second) << endl;
  cout << Form("both muons are global muons         &$%4.3f\\pm%4.3f$ &$%4.3f\\pm%4.3f$", sgglobalmuon.first, sgglobalmuon.second, bgglobalmuon.first, bgglobalmuon.second) << endl;
  cout << Form("both muons pass MVA muon ID         &$%4.3f\\pm%4.3f$ &$%4.3f\\pm%4.3f$", sgmuid.first, sgmuid.second, bgmuid.first, bgmuid.second) << endl;


}


// ----------------------------------------------------------------------
void dbx() {

  string fname("/scratch/ursl/bmm4/s01/bmm-rare-rhlt1.root");

  string ams1("(m1rmvabdt>-2.5 && m2rmvabdt>-2.5) && (m1rmvabdt<0.58 || m2rmvabdt<0.58)");
  string exa1("(m1rmvabdt>-2.5 && m2rmvabdt>-2.5) && ((m1rmvabdt<0.58 && m2rmvabdt>0.58) || (m1rmvabdt>0.58 && m2rmvabdt<0.58))");
  string hipt("(m1rmvabdt>-2.5 && m2rmvabdt>-2.5) && ((m2rmvabdt<0.58 && m1rmvabdt>0.58))");
  string lopt("(m1rmvabdt>-2.5 && m2rmvabdt>-2.5) && ((m2rmvabdt>0.58 && m1rmvabdt<0.58))");

  if (0) {
    std::pair<double, double> bghipt = triggerEff("rhlt1", lopt, fname);
    c0->SaveAs("triggerEffBg-lopt.pdf");
  }

  if (0) {
    std::pair<double, double> bghipt = triggerEff("rhlt1", exa1, fname);
    c0->SaveAs("triggerEffBg-exa1.pdf");
  }

  if (1) {
    std::pair<double, double> bgams1 = triggerEff("rhlt1", ams1, fname);
    c0->SaveAs("triggerEffBg-ams1.pdf");
  }
}
