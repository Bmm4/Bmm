#include <string>
#include <stdio.h>

// ----------------------------------------------------------------------
void drawPoint(TVector3 v, int color, int style, double size) {
  TPolyMarker3D *pm3d1 = new TPolyMarker3D(1);
  pm3d1->SetPoint(0, v.X(), v.Y(), v.Z());
  pm3d1->SetMarkerColor(color);
  pm3d1->SetMarkerStyle(style);
  pm3d1->SetMarkerSize(size);
  pm3d1->Draw();
}

// ----------------------------------------------------------------------
void drawVector(TVector3 v, TVector3 p, int color, int style) {
  TPolyLine3D *pl3d1 = new TPolyLine3D(2);
  pl3d1->SetPoint(0, v.X(), v.Y(), v.Z());
  pl3d1->SetPoint(1, v.X() + p.X(), v.Y() + p.Y(), v.Z() + p.Z());
  pl3d1->SetLineColor(color);
  pl3d1->SetLineStyle(style);
  pl3d1->Draw();

  TPolyMarker3D *pm3d1 = new TPolyMarker3D(1);
  pm3d1->SetPoint(0, v.X() + p.X(), v.Y() + p.Y(), v.Z() + p.Z());
  pm3d1->SetMarkerColor(color);
  pm3d1->SetMarkerStyle(20);
  pm3d1->SetMarkerSize(0.3);
  pm3d1->Draw();
}

// ----------------------------------------------------------------------
void readVector(string label, TVector3 &reco, TVector3 &gen, string fileName) {
  ifstream INS;
  string sline;
  INS.open(fileName);
  float rx(99.), ry(99.), rz(99.), gx(99.), gy(99.), gz(99.);
  while (getline(INS, sline)) {
    if (string::npos != sline.find("#")) continue;
    if (string::npos != sline.find(label)) {
      sline.replace(0, label.length()+2, "");
      cout << label << ": " << sline << endl;
      sscanf(sline.c_str(), "reco(%f, %f, %f) gen(%f, %f, %f)", &rx, &ry, &rz, &gx, &gy, &gz);
    }
  }
  reco.SetXYZ(rx, ry, rz);
  gen.SetXYZ(gx, gy, gz);
}



// ----------------------------------------------------------------------
void bmtkDraw(string fileName) {

  TVector3 pvr, pvg;
  TVector3 svr, svg;
  TVector3 tvr, tvg;
  TVector3 h1pr, h1pg;
  TVector3 h2pr, h2pg;
  TVector3 h3pr, h3pg;
  TVector3 mpr, mpg;
  TVector3 kpr, kpg;

  readVector("PV", pvr, pvg, fileName);
  readVector("SV", svr, svg, fileName);
  readVector("TV", tvr, tvg, fileName);
  readVector("had1", h1pr, h1pg, fileName);
  readVector("had2", h2pr, h2pg, fileName);
  readVector("had3", h3pr, h3pg, fileName);
  readVector("muon", mpr, mpg, fileName);
  readVector("kaon", kpr, kpg, fileName);

  double xmin(-0.6), xmax(0.6), ymin(-0.6), ymax(0.6), zmin(-10.), zmax(10.);
  zmin = static_cast<int>(pvr.Z());
  zmax = zmin + 2.;

  xmin = tvr.X() - TMath::Abs(pvr.X()-tvr.X());
  xmax = tvr.X() + TMath::Abs(pvr.X()-tvr.X());

  ymin = tvr.Y() - TMath::Abs(pvr.Y()-tvr.Y());
  ymax = tvr.Y() + TMath::Abs(pvr.Y()-tvr.Y());

  // readVector("muon", tvr, tvg, fileName);
  // readVector("kaon", tvr, tvg, fileName);

  TCanvas *cV3D = new TCanvas("cV3D","PolyLine3D & PolyMarker3D Window",200,10,500,500);
  // Creating a view
  TView3D *view = (TView3D*) TView::CreateView(1);
  view->ShowAxis();
  view->SetRange(-0.4, -0.3, -2., 0.4, 0.3, 2.);
  //  view->RotateView(45., 45., gPad);
  view->FrontView(gPad);
  view->TopView(gPad);

  double rsize(1.0);
  drawPoint(pvr, kBlack, 30, 2.0);
  drawPoint(svr, kBlack, 20, rsize);
  drawPoint(tvr, kRed, 20, rsize);

  double gsize(1.5);
  drawPoint(svg, kBlack, 28, gsize);
  drawPoint(tvg, kRed, 28, gsize);


  drawVector(tvr, h1pr, kMagenta, kSolid);
  drawVector(tvr, h2pr, kMagenta, kSolid);
  drawVector(tvr, h3pr, kMagenta, kSolid);

  drawVector(svr, mpr, kRed, kSolid);
  drawVector(svr, kpr, kGreen+1, kSolid);

  drawVector(svr, 3.*(tvr-svr), kBlack, kDashed);

  drawVector(TVector3(pvr.X(), pvr.Y(), -20.), TVector3(0., 0., 40.), kBlack, kDotted);

  return;

}


// ----------------------------------------------------------------------
void bmtkDrawEvent(int ievt, string fileName= "bmtkDraw-events/evt") {
  bmtkDraw(Form("%s-%d.txt", fileName.c_str(), ievt));
}


// ----------------------------------------------------------------------
void dumpEvents(string fileName, string startLine = "PV: reco", string endLine = "had3: reco") {
  string sline, prevLine;
  ifstream INS;
  INS.open(fileName);
  int cnt(0);
  while (getline(INS, sline)) {
    if (string::npos != sline.find(startLine)) {
      ofstream OS(Form("bmtkDraw-events/evt-%d.txt", cnt));
      OS << "# " << prevLine << endl;
      OS << sline << endl;
      while (getline(INS, sline)) {
	OS << sline << endl;
	if (string::npos != sline.find(endLine)) break;
      }
      OS.close();
      ++cnt;

    }
    prevLine = sline;
  }
}
