#include "common/util.hh"

void pv2011() {
  vector<TFile *> v;
  v.push_back(TFile::Open("/scratch/ursl/bmm4/v06/bmm-data-bmmMuOnia2011-v06.root"));
  v.push_back(TFile::Open("/scratch/ursl/bmm4/v06/bmm-mc-Summer17_private-BuToJpsiKp-v06.root"));

  vector<string> hname;
  hname.push_back("npv");
  hname.push_back("pvz");
  hname.push_back("pv0z");

  TH1D *h(0);
  gStyle->SetOptStat(0);
  for (unsigned int ih = 0; ih < hname.size(); ++ih) {
    TLegend *t = newLegend(0.2, 0.7, 0.88, 0.88);
    t->SetTextSize(0.03);
    for (unsigned int i = 0; i < v.size(); ++i) {
      h = (TH1D*)(v[i]->Get(hname[ih].c_str()));
      h->Scale(1./h->GetSumOfWeights());
      if (0 == i) {
	h->SetMaximum(1.4*h->GetMaximum());
	h->Draw("e");
	t->AddEntry(h, Form("2011 (mean/RMS):  %5.4f/%5.4f", h->GetMean(), h->GetRMS()), "p");
      } else {
	if (1 == i) {
	  h->SetLineColor(kRed);
	  h->Draw("histsame");
	  t->AddEntry(h, Form("MC:                         %5.4f/%5.4f", h->GetMean(), h->GetRMS()), "l");
	}
      }
    }
    t->Draw();
    c0->SaveAs(Form("qnd/pv2011-%s.pdf", hname[ih].c_str()));
  }

  TTree *tt(0);
  TH1D *hc0 = new TH1D("hc0", "", 100, -25., 25.);
  TH1D *hc1 = new TH1D("hc1", "", 100, -25., 25.);
  TH1D *hc2 = new TH1D("hc2", "", 100, -25., 25.);
  for (unsigned int i = 0; i < v.size(); ++i) {
    tt = (TTree*)(v[i]->Get("candAnaBu2JpsiK/events"));
    tt->Draw(Form("pvz>>hc%d", i));
  }
  TLegend *t = newLegend(0.2, 0.7, 0.88, 0.88);
  t->SetTextSize(0.03);
  hc0->Scale(1./hc0->GetSumOfWeights());
  hc0->SetMaximum(1.4*hc0->GetMaximum());
  hc0->Draw("e");
  t->AddEntry(hc0, Form("2011 (mean/RMS):  %5.4f/%5.4f", hc0->GetMean(), hc0->GetRMS()), "p");

  hc1->Scale(1./hc1->GetSumOfWeights());
  hc1->SetLineColor(kRed);
  hc1->Draw("samehist");
  t->AddEntry(hc1, Form("MC:                            %5.4f/%5.4f", hc1->GetMean(), hc1->GetRMS()), "l");

  t->Draw();
  c0->SaveAs(Form("qnd/pv2011-cand.pdf"));

}



void pv2012() {
  vector<TFile *> v;
  v.push_back(TFile::Open("/scratch/ursl/bmm4/v06/bmm-data-bmmMuOnia2012-v06.root"));
  v.push_back(TFile::Open("/scratch/ursl/bmm4/v06/bmm-mc-Winter17_private-BuToJpsiKp-v06.root"));

  vector<string> hname;
  hname.push_back("npv");
  hname.push_back("pvz");
  hname.push_back("pv0z");

  TH1D *h(0);
  gStyle->SetOptStat(0);
  for (unsigned int ih = 0; ih < hname.size(); ++ih) {
    TLegend *t = newLegend(0.2, 0.7, 0.88, 0.88);
    t->SetTextSize(0.03);
    for (unsigned int i = 0; i < v.size(); ++i) {
      h = (TH1D*)(v[i]->Get(hname[ih].c_str()));
      h->Scale(1./h->GetSumOfWeights());
      if (0 == i) {
	h->SetMaximum(1.4*h->GetMaximum());
	h->Draw("e");
	t->AddEntry(h, Form("2012 (mean/RMS):  %5.4f/%5.4f", h->GetMean(), h->GetRMS()), "p");
      } else {
	if (1 == i) {
	  h->SetLineColor(kRed);
	  h->Draw("histsame");
	  t->AddEntry(h, Form("MC:                         %5.4f/%5.4f", h->GetMean(), h->GetRMS()), "l");
	}
      }
    }
    t->Draw();
    c0->SaveAs(Form("qnd/pv2012-%s.pdf", hname[ih].c_str()));
  }

  TTree *tt(0);
  TH1D *hc0 = new TH1D("hc0", "", 100, -25., 25.);
  TH1D *hc1 = new TH1D("hc1", "", 100, -25., 25.);
  TH1D *hc2 = new TH1D("hc2", "", 100, -25., 25.);
  for (unsigned int i = 0; i < v.size(); ++i) {
    tt = (TTree*)(v[i]->Get("candAnaBu2JpsiK/events"));
    tt->Draw(Form("pvz>>hc%d", i));
  }
  TLegend *t = newLegend(0.2, 0.7, 0.88, 0.88);
  t->SetTextSize(0.03);
  hc0->Scale(1./hc0->GetSumOfWeights());
  hc0->SetMaximum(1.4*hc0->GetMaximum());
  hc0->Draw("e");
  t->AddEntry(hc0, Form("2012 (mean/RMS):  %5.4f/%5.4f", hc0->GetMean(), hc0->GetRMS()), "p");

  hc1->Scale(1./hc1->GetSumOfWeights());
  hc1->SetLineColor(kRed);
  hc1->Draw("samehist");
  t->AddEntry(hc1, Form("MC:                            %5.4f/%5.4f", hc1->GetMean(), hc1->GetRMS()), "l");

  t->Draw();
  c0->SaveAs(Form("qnd/pv2012-cand.pdf"));

}

void pv() {
  vector<TFile *> v;
  v.push_back(TFile::Open("cbmm-data-bmmCharmonium2016H-v06.bmmReader.2016.root"));
  v.push_back(TFile::Open("cbmm-mc-RunIISpring16DR80-BuToJpsiK-v06.bmmReader.mix-Bu2JpsiK.root"));
  v.push_back(TFile::Open("cbmm-mc-bmmMoriond2017-BuToJpsiK_BMuonFilter-v06.bmmReader.mix-Bu2JpsiK.root"));

  vector<string> hname;
  hname.push_back("npv");
  hname.push_back("pvz");
  hname.push_back("pv0z");

  TH1D *h(0);
  gStyle->SetOptStat(0);
  for (unsigned int ih = 0; ih < hname.size(); ++ih) {
    TLegend *t = newLegend(0.2, 0.7, 0.88, 0.88);
    t->SetTextSize(0.03);
    for (unsigned int i = 0; i < v.size(); ++i) {
      h = (TH1D*)(v[i]->Get(hname[ih].c_str()));
      h->Scale(1./h->GetSumOfWeights());
      if (0 == i) {
	h->SetMaximum(1.4*h->GetMaximum());
	h->Draw("e");
	t->AddEntry(h, Form("2016H (mean/RMS):  %5.4f/%5.4f", h->GetMean(), h->GetRMS()), "p");
      } else {
	if (1 == i) {
	  h->SetLineColor(kRed);
	  t->AddEntry(h, Form("RunIISpring16DR80: %5.4f/%5.4f", h->GetMean(), h->GetRMS()), "l");
	}
	if (2 == i) {
	  h->SetLineColor(kBlue);
	  t->AddEntry(h, Form("TrancheIV:                %5.4f/%5.4f", h->GetMean(), h->GetRMS()), "l");
	}
	h->Draw("histsame");
      }
    }
    t->Draw();
    c0->SaveAs(Form("qnd/pv-%s.pdf", hname[ih].c_str()));
  }

  TTree *tt(0);
  TH1D *hc0 = new TH1D("hc0", "", 100, -25., 25.);
  TH1D *hc1 = new TH1D("hc1", "", 100, -25., 25.);
  TH1D *hc2 = new TH1D("hc2", "", 100, -25., 25.);
  for (unsigned int i = 0; i < v.size(); ++i) {
    tt = (TTree*)(v[i]->Get("candAnaBu2JpsiK/events"));
    tt->Draw(Form("pvz>>hc%d", i));
  }
  TLegend *t = newLegend(0.2, 0.7, 0.88, 0.88);
  t->SetTextSize(0.03);
  hc0->Scale(1./hc0->GetSumOfWeights());
  hc0->SetMaximum(1.4*hc0->GetMaximum());
  hc0->Draw("e");
  t->AddEntry(hc0, Form("2016H (mean/RMS):  %5.4f/%5.4f", hc0->GetMean(), hc0->GetRMS()), "p");

  hc1->Scale(1./hc1->GetSumOfWeights());
  hc1->SetLineColor(kRed);
  hc1->Draw("samehist");
  t->AddEntry(hc1, Form("RunIISpring16DR80: %5.4f/%5.4f", hc1->GetMean(), hc1->GetRMS()), "l");

  hc2->Scale(1./hc2->GetSumOfWeights());
  hc2->SetLineColor(kBlue);
  hc2->Draw("samehist");
  t->AddEntry(hc2, Form("TrancheIV:                %5.4f/%5.4f", hc2->GetMean(), hc2->GetRMS()), "l");

  t->Draw();
  c0->SaveAs(Form("qnd/pv-cand.pdf"));


}
