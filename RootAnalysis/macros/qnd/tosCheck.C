// ----------------------------------------------------------------------
void tosCheck(string filename = "/scratch/ursl/bmm4/s01/bmm-data-bmmMuOnia2011-s01.root",
	      string cuts = "chan==0 && m2pt>4 & &m>4.5&&m<6 && hlt1",
	      int run1 = 163268, int run2 = 180252) {
  TFile *file = TFile::Open(filename.c_str());
  int NBINS(20);
  TH1D *heff = new TH1D("heff", "", NBINS, run1, run2);
  heff->Sumw2();
  TH1D *h0 = new TH1D("h0", "", NBINS, run1, run2);
  TH1D *h1 = new TH1D("h1", "", NBINS, run1, run2);
  //  TTree *events = (TTree*)file->Get("candAnaMuMu/events");
  TTree *events = (TTree*)file->Get("candAnaBu2JpsiK/events");

  string acuts = cuts + "&&tos";
  events->Draw("run>>h0", cuts.c_str());
  events->Draw("run>>h1", acuts.c_str());

  setTitles(heff, "run", Form("Eff((tos | %s)", cuts.c_str()), 0.03, 1.5, 1.7);

  heff->Divide(h1, h0, 1., 1., "b");
  shrinkPad(0.13, 0.15, 0.13);
  heff->Draw("e");
  gStyle->SetOptStat(0);

}


void plotAll() {
  tosCheck("/scratch/ursl/bmm4/s01/bmm-data-bmmMuOnia2011-s01.root", "chan==0 && m2pt>4 & &m>4.5&&m<6 && hlt1 &&l1t", 163268, 180253);
  c0->SaveAs("tosCheckNo-chan0-2011.pdf");
  tosCheck("/scratch/ursl/bmm4/s01/bmm-data-bmmMuOnia2011-s01.root", "chan==1 && m2pt>4 & &m>4.5&&m<6 && hlt1 &&l1t", 163268, 180253);
  c0->SaveAs("tosCheckNo-chan1-2011.pdf");

  tosCheck("/scratch/ursl/bmm4/s01/bmm-data-bmmMuOnia2012-s01.root", "chan==0 && m2pt>4 & &m>4.5&&m<6 && hlt1 &&l1t", 190456, 208687);
  c0->SaveAs("tosCheckNo-chan0-2012.pdf");
  tosCheck("/scratch/ursl/bmm4/s01/bmm-data-bmmMuOnia2012-s01.root", "chan==1 && m2pt>4 & &m>4.5&&m<6 && hlt1 &&l1t", 190456, 208687);
  c0->SaveAs("tosCheckNo-chan1-2012.pdf");

  tosCheck("/scratch/ursl/bmm4/s01/bmm-data-bmmCharmonium2016BF-s01.root", "chan==0 && m2pt>4 & &m>4.5&&m<6 && hlt1 &&l1t", 272007, 277000);
  c0->SaveAs("tosCheckNo-chan0-2016BF.pdf");
  tosCheck("/scratch/ursl/bmm4/s01/bmm-data-bmmCharmonium2016BF-s01.root", "chan==1 && m2pt>4 & &m>4.5&&m<6 && hlt1 &&l1t", 272007, 277000);
  c0->SaveAs("tosCheckNo-chan1-2016BF.pdf");

  tosCheck("/scratch/ursl/bmm4/s01/bmm-data-bmmCharmonium2016GH-s01.root", "chan==0 && m2pt>4 & &m>4.5&&m<6 && hlt1 &&l1t", 277000, 284045);
  c0->SaveAs("tosCheckNo-chan0-2016GH.pdf");
  tosCheck("/scratch/ursl/bmm4/s01/bmm-data-bmmCharmonium2016GH-s01.root", "chan==1 && m2pt>4 & &m>4.5&&m<6 && hlt1 &&l1t", 277000, 284045);
  c0->SaveAs("tosCheckNo-chan1-2016GH.pdf");

}



void compare(string fname = "bmm-data-bmmMuOnia2011-s01.root") {


  TFile *f1 = TFile::Open(Form("/scratch/ursl/bmm4/s01/%s", fname.c_str()));
  TFile *f2 = TFile::Open(Form("/scratch/ursl/bmm4/s01bis/%s", fname.c_str()));

  string cuts("hlt1&&tos&&l1t&&gmuid&&bdt>0.25");


  TH1D *h10 = new TH1D("h10", "", 40, 4.9, 5.9);
  TH1D *h20 = new TH1D("h20", "", 40, 4.9, 5.9);

  TH1D *h11 = new TH1D("h11", "", 40, 4.9, 5.9);
  TH1D *h21 = new TH1D("h21", "", 40, 4.9, 5.9);

  TH1D *h1  = new TH1D("h1", "", 40, 4.9, 5.9);
  TH1D *h2  = new TH1D("h2", "", 40, 4.9, 5.9);

  TTree *t1 = (TTree*)f1->Get("candAnaMuMu/events");
  TTree *t2 = (TTree*)f2->Get("candAnaMuMu/events");

  t1->Draw("m>>h10", Form("chan==0 && %s", cuts.c_str()));
  t1->Draw("m>>h11", Form("chan==1 && %s", cuts.c_str()));

  t2->Draw("m>>h20", Form("chan==0 && %s", cuts.c_str()));
  t2->Draw("m>>h21", Form("chan==1 && %s", cuts.c_str()));


  zone(2,2);

  h10->Draw("hist");
  h20->Draw("samee");

  c0->cd(2);
  h1->Add(h20, h10, 1., -1.);
  h1->Draw();

  c0->cd(3);
  h11->Draw("hist");
  h21->Draw("samee");

  c0->cd(4);
  h2->Add(h21, h11, 1., -1.);
  h2->Draw();

  string pdf(fname);
  replaceAll(pdf, ".root", "");
  c0->SaveAs(Form("qnd/compare-%s.pdf", pdf.c_str()));

}
