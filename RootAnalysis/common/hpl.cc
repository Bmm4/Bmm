#include "hpl.hh"

#include <iostream>
#include <fstream>
#include <iomanip>


#include "TVirtualPad.h"  // access to gPad
#include "TRegexp.h"

using std::cout;
using std::endl;

// ----------------------------------------------------------------------
void hpl(const char *h, const char *options) { 

  TH1 *h1 = (TH1*)gDirectory->Get(h); 
  if (!h1) {
    cout << h << " not found" << endl;
  } else {
    hpl(h1, options);
  }
}


// ----------------------------------------------------------------------
void hpl(TH1 *h, const char *options) { 

  TString parse(options);

  if (!h) {
    cout << "histogram not found" << endl;
  }

  // -- "same" shortcut
  if (parse.Contains(",s", TString::kIgnoreCase)) {
    parse += TString("same");
  }

  // -- limits
  float xmin(99.), xmax(98.),  ymin(99.), ymax(98.);
  TRegexp re1("(.*:.*)");
  TRegexp re2("(.*:.*,.*:.*)");
  TSubString bla2 = parse(re2);
  sscanf(bla2.Data(), "(%f:%f,%f:%f)", &xmin, &xmax, &ymin, &ymax);
  parse(re2) = "";
  if ((xmin < ymax) && (ymin < ymax)) {
    cout << "restricting x to " << xmin << " .. " << xmax << " and y to " << ymin << " .. " << ymax << endl;
    h->SetAxisRange(xmin, xmax);
    h->SetMinimum(ymin);
    h->SetMaximum(ymax);
  }
  TSubString  bla1 = parse(re1);
  sscanf(bla1.Data(), "(%f:%f)", &xmin, &xmax);
  parse(re1) = "";
  if (xmin < xmax) {
    cout << "restricting to " << xmin << " .. " << xmax << endl;
    h->SetAxisRange(xmin, xmax);
  }


  // -- normalize histogram to histogram already plotted
  if (parse.Contains("norm", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("norm"), "");
    TList *tl = gPad->GetListOfPrimitives();
    if (tl) {
      TListIter next(tl);
      TObject *to;
      double norma(-99.); 
      TH1 *h1;
      while ((to=next())) {
	if (to->InheritsFrom(TH1::Class())) {
	  h1 = (TH1*)to;
	  norma = h1->GetSumOfWeights();
	  cout << h1->GetName() << "  " << h1->GetTitle() << " with " << norma << " events "  << endl;
	  break; // Take first one
	}
      }
      if (norma > -99.) {
	h->Scale(norma/h->GetSumOfWeights());
      }
    }
  }
  


  // -- grids, logx
  if (parse.Contains("logy", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("logy"), "");
    gPad->SetLogy(1);
  }
  if (parse.Contains("liny", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("liny"), "");
    gPad->SetLogy(0);
  }
  if (parse.Contains("logx", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("logx"), "");
    gPad->SetLogx(1);
  }
  if (parse.Contains("linx", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("linx"), "");
    gPad->SetLogx(0);
  }
  if (parse.Contains("logz", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("logz"), "");
    gPad->SetLogz(1);
  }
  if (parse.Contains("linz", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("linz"), "");
    gPad->SetLogz(0);
  }
  if (parse.Contains("gridx", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("gridx"), "");
    gPad->SetGridx(1);
  }
  if (parse.Contains("gridy", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("gridy"), "");
    gPad->SetGridy(1);
  }
  if (parse.Contains("ngridx", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("ngridx"), "");
    gPad->SetGridx(0);
  }
  if (parse.Contains("ngridy", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("ngridy"), "");
    gPad->SetGridy(0);
  }

  // -- Colors, fillstyle
  if (parse.Contains("red", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("red"), "");
    h->SetMarkerColor(kRed);
    h->SetLineColor(kRed);
    h->SetFillColor(kRed);
  }
  else if (parse.Contains("green", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("green"), "");
    h->SetMarkerColor(kGreen);
    h->SetLineColor(kGreen);
    h->SetFillColor(kGreen);
  }
  else if (parse.Contains("yellow", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("yellow"), "");
    h->SetMarkerColor(kYellow);
    h->SetLineColor(kYellow);
    h->SetFillColor(kYellow);
  }
  else if (parse.Contains("blue", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("blue"), "");
    h->SetMarkerColor(kBlue);
    h->SetLineColor(kBlue);
    h->SetFillColor(kBlue);
  }
  else if (parse.Contains("cyan", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("cyan"), "");
    h->SetMarkerColor(kCyan);
    h->SetLineColor(kCyan);
    h->SetFillColor(kCyan);
  }
  else {
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
  }

  if (parse.Contains("dot", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("dot"), "");
    h->SetLineStyle(3);
  }
  else if (parse.Contains("dash", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("dash"), "");
    h->SetLineStyle(2);
  }

  if (parse.Contains("fill", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("fill"), "");
    h->SetLineColor(kBlack);
    h->SetFillStyle(1000);
  }
  else if (parse.Contains("upr", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("upr"), "");
    h->SetFillStyle(3004);
  }
  else if (parse.Contains("upl", TString::kIgnoreCase)) {
    parse.ReplaceAll(TString("upl"), "");
    h->SetFillStyle(3005);
  }
  else {
    h->SetFillStyle(0);
  }

  //  cout << parse << endl;
  h->DrawCopy(parse);
}
