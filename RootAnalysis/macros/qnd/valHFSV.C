void val(string);


// ----------------------------------------------------------------------
void valAll() {
  val("mm"); 
  val("ks"); 
  val("bu");
  val("bs"); 
  val("dstar"); 
}


// ----------------------------------------------------------------------
void val(string particle = "bu" ) {

  string dir("candAnaMuMu"); 
  
  TFile *fNew = TFile::Open("new.root"); 
  TFile *fOld = TFile::Open("old.root");

  vector<string> hname;
  hname.push_back(Form("%smass", particle.c_str())); 
  hname.push_back(Form("%schi2", particle.c_str())); 
  hname.push_back(Form("%sschi2", particle.c_str())); 
  hname.push_back(Form("%sflsxy", particle.c_str())); 
  hname.push_back(Form("%sfls3d", particle.c_str())); 
  hname.push_back(Form("%spvips", particle.c_str())); 
  hname.push_back(Form("%sdocatrk", particle.c_str())); 
  hname.push_back(Form("%salpha", particle.c_str())); 
  hname.push_back(Form("%stip", particle.c_str())); 
  hname.push_back(Form("%slip", particle.c_str())); 
  hname.push_back(Form("%spvw8", particle.c_str())); 
  hname.push_back(Form("%sdm", particle.c_str())); 

  c0->SetWindowSize(1000, 800); 
  c0->Clear();
  c0->Divide(3, 4);
  c0->cd(1);
  
  TH1D *hNew(0), *hOld(0), *hRatio(0); 
  for (unsigned int i = 0; i < hname.size(); ++i) {
    hNew = (TH1D*)fNew->Get(Form("%s/%s", dir.c_str(), hname[i].c_str())); 
    if (0 == hNew) continue;
    hOld = (TH1D*)fOld->Get(Form("%s/%s", dir.c_str(), hname[i].c_str())); 
    c0->cd(i+1);
    hNew->SetMinimum(0.); 
    hNew->Draw();
    hOld->SetMarkerStyle(24);
    hOld->Draw("samep");
  }
  c0->SaveAs(Form("valHFSV-overlay-%s.pdf", particle.c_str())); 

  c0->SetWindowSize(1000, 800); 
  c0->Clear();
  c0->Divide(3, 4);
  c0->cd(1);
  for (unsigned int i = 0; i < hname.size(); ++i) {
    hNew = (TH1D*)fNew->Get(Form("%s/%s", dir.c_str(), hname[i].c_str())); 
    if (0 == hNew) continue;
    hOld = (TH1D*)fOld->Get(Form("%s/%s", dir.c_str(), hname[i].c_str())); 
    hRatio = (TH1D*)hNew->Clone(Form("ratio_%s", hNew->GetName()));
    hRatio->Clear();
    hRatio->Divide(hNew, hOld); 
    hRatio->SetTitle(hNew->GetName()); 
    hRatio->SetMinimum(0.); 
    hRatio->SetMaximum(1.2); 
    c0->cd(i+1);
    hRatio->Draw();
  }
  c0->SaveAs(Form("valHFSV-ratio-%s.pdf", particle.c_str())); 
  

}
