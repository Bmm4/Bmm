{
#include <string>
#include "../common/util.hh"

  using namespace std;
  using namespace RooFit;
  using namespace RooStats;

  string version = gSystem->Getenv("VERSION");

  gROOT->Macro("cms-tdr.C");
  gROOT->ForceStyle();
  gStyle->SetPalette(1);

  gStyle->SetOptStat(1112211);


  cout << "Hallo Loading rootlibs: lib/libAnaClasses.so lib/libUtil.so" << endl;
  gSystem->Load("lib/libRecoilClasses.so");
  gSystem->Load("lib/libUtil.so");
  gSystem->Load("libTMVA.so");


  TLatex *tl = new TLatex();
  tl->SetNDC(kTRUE);
  tl->SetTextFont(42);

  TLine *pl = new TLine();

  // --- Cleanup if this is not the first call to rootlogon.C
  TCanvas *c = 0;
  c = (TCanvas*)gROOT->FindObject("c0"); if (c) c->Delete(); c = 0;
  p = (TPad*)gROOT->FindObject("p0"); if (p) p->Delete(); p = 0;
  // --- Create a new canvas.
  TCanvas *c0 = new TCanvas("c0","--c0--",2303,0,656,700);
  c0->ToggleEventStatus();
}
