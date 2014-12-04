{
#include <string>
  
  using namespace std;
  using namespace RooFit;
  using namespace RooStats;

  string version = gSystem->Getenv("VERSION"); 
  
  cout << "Loading rootlibs" << endl;
  gSystem->Load("libPhysics.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libHtml.so");
  
  gSystem->Load("libRooStats.so");
  gSystem->Load("libRooFit.so");
  gSystem->Load("libRooFitCore.so");
  gSystem->Load("libTMVA.so");

  gSystem->Load("lib/libAnaClasses.so");
  
  gROOT->Macro("cms-tdr.C");
  gROOT->ForceStyle();
  gStyle->SetPalette(1);
  
  TLatex *tl = new TLatex();
  tl->SetNDC(kTRUE);
  tl->SetTextFont(42); 
  
  TLine *pl = new TLine();  
  
  // --- Cleanup if this is not the first call to rootlogon.C
  TCanvas *c = 0;
  c = (TCanvas*)gROOT->FindObject("c0"); if (c) c->Delete(); c = 0;
  p = (TPad*)gROOT->FindObject("p0"); if (p) p->Delete(); p = 0;
  // --- Create a new canvas.
  TCanvas c0("c0","--c0--",2303,0,656,700);
  c0->ToggleEventStatus();
}



