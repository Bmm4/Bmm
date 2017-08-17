#include "candAnaDstar.hh"
#include <cmath>
#include <string>

#include "common/HFMasses.hh"
#include "danekUtils.h"

#define DO_TESTS
//#define MC_HISTOS
//#define OLD_OBJ_MARK  // old way of marking active hlt modules

using namespace std;

namespace {
  TVector3 DSVertex(0,0,0), DZVertex(0,0,0), PV(0,0,0);
  TVector3 DSMom(0,0,0), DZMom(0,0,0), PiSlowMom(0,0,0), PiMom(0,0,0), KMom(0,0,0);
  bool MYDEBUG=false;
}

// ----------------------------------------------------------------------
candAnaDstar::candAnaDstar(bmmReader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  cout << "==> candAnaDstar: name = " << name << ", reading cutsfile " << cutsFile << endl;

  readCuts(cutsFile, 1);
  if (fVerbose>10) MYDEBUG=true;

}


// ----------------------------------------------------------------------
candAnaDstar::~candAnaDstar() {
  cout << "==> candAnaDstar: destructor..." << endl;
}


//---------------------------------------------------------------------------------------------------------------
// Main candidate processing for Dstar
void candAnaDstar::candAnalysis() {
  static int count0=0, count1=0, count2=0, count3=0, count4=0, count5=0, count6=0;
  //static int tcount1=0, tcount2=0;
  if(MYDEBUG) cout<<" In candAnaDstar::candAnalysis "<<endl;
  fPreselection = false;  //  reset
  fmds=0.; // this is to singal for MC that the event did not pass presselection
  fmdz=0.;
  fpt=0.;

  if (0 == fpCand) return; // skip if no cand

  if(MYDEBUG) cout<<" In candAnaDstar::candAnalysis, call candAna::candAnalysis "<<endl;
  candAna::candAnalysis();  // call the main analysis

  // now do Dstar specific
  count0++;

  ((TH1D*)fHistDir->Get("Status"))->Fill(0.);

  //cout<<fhltType<<" "<<int(fhltType/1000)<<" "<<(fhltType%1000)<<endl;

  if (fVerbose>2) {
    cout<<" Dstar candidate "<<fpCand->fType<<" in event "<<fEvt<<" v "<<fVerbose;
    cout << " with mass = " << fpCand->fMass <<" cand num "<<count0<<endl;
    if (fVerbose>10) {
      cout << "DUMP HFDstarCandidate  " <<endl;
      dumpHFDstarCand(fpCand);
      //dumpHFTruthCand(fpCand);
      //doTest(fpCand,0); // testing all collections
    }
  }

  TAnaCand *pC(0);
  TAnaTrack *pK, *pPi, *pPis;
  double fls3d(0.), flsxy(0.), prob(0.), chi2(0.), alpha(0.), dr(0.);
  int piIndex=-1, KIndex=-1, piSlowIndex=-1;

  // -- D0
  if (fpCand->fDau1 < 0  || fpCand->fDau2 < 0) {
    if(fVerbose>1) {cout << "pCand->fDauX = -1!!! " << fpCand->fType << " skip event "<<endl; fpCand->dump();}
    return;  // skip if no daughters
  }
  ((TH1D*)fHistDir->Get("Status"))->Fill(1.);

  pC = fpEvt->getCand(fpCand->fDau1);  // D0 candidate
  pK = 0;
  pPi =0;
  piSlowIndex = fpEvt->getSigTrack(fpCand->fSig1)->fIndex; // slow pion index
  int pPisId = fpCand->fSig1; // slow pion index
  pPis = fpEvt->getSigTrack(pPisId); // slow pi from D*
  TVector3 piSlowMom = pPis->fPlab; // slow pi momentum vector

  if (fVerbose>8 ) {
    cout<<" found D0 "<<pC->fType<<endl;
    cout<<" found slow pi "<<pPis->fMCID<<" "
	<<pPis->fPlab.Perp()<<" "<<pPis->fPlab.Eta()<<" "<<pPis->fPlab.Phi()<<" "
	<<piSlowIndex<<endl;
  }

  // loop over D0 tracks
  for (int id = pC->fSig1; id <= pC->fSig2; ++id) {
    int index = fpEvt->getSigTrack(id)->fIndex;

    if (211 == fpEvt->getSigTrack(id)->fMCID) {  // pion
     //pPi = fpEvt->getRecTrack(id);
      pPi = fpEvt->getSigTrack(id);
      piIndex = index;
      if (fVerbose>8)
	cout<<" found pi "<<pPi->fMCID<<" "<<pPi->fPlab.Perp()<<" "<<pPi->fPlab.Eta()
	    <<" "<<pPi->fPlab.Phi() << " "<<piIndex<<endl;
    } else {  // kaon
      //pK = fpEvt->getRecTrack(id);
      pK = fpEvt->getSigTrack(id);
      KIndex = index;
      if (fVerbose>8 )
	cout<<" found K "<<pK->fMCID<<" "<<pK->fPlab.Perp()<<" "<<pK->fPlab.Eta()
	    <<" "<<pK->fPlab.Phi() << " "<<KIndex<<endl;
    }
  }

  if(pPi == 0 || pK==0) {
    if(fVerbose>0) {cout << " pi or K not found " << fpCand->fType << endl; fpCand->dump();}
    return;  // skip if no daughters
  }
  ((TH1D*)fHistDir->Get("Status"))->Fill(2.);


  int tm = 0;
  int mcOk = false;
  if(fIsMC) {
    mcOk = anaMC(false);  // check if dstar in MC

    if(mcOk==1) {
      //truthMatch return an interger, 0-no match, 1-correct match, -1-pi&K swapped
      tm = truthMatch(fpCand,fVerbose); // check truth matching
      if ( tm == 1 ) { // do only for matched MC
	if (fVerbose>10) cout << " MC matched -> " << fpCand->fType <<endl;
      } if (tm==-1) {
	cout<<" mixed pi-K"<<endl;
      } else {
	if (fVerbose>1) cout<<" Dstar cand not matched "<<endl;
	//anaMC(true); //dumpHFTruthCand(fpCand);
      }
            
    } else { // no clean Dstar in MC
      if(fVerbose>10) cout<<" Clean Dstar does not exist in MC "<<mcOk<<endl;
      //anaMC(true); //dumpHFTruthCand(fpCand);
    }
    
  } // end if MC

  // Pt, eta
  double pt    = fpCand->fPlab.Perp(); // D* cand pt
  double ptdz  = pC->fPlab.Perp();  // D0
  double ptPis = pPis->fPlab.Perp();
  double ptPi  = pPi->fPlab.Perp();
  double ptK   = pK->fPlab.Perp();

  double eta   = fpCand->fPlab.Eta();
  //double etaPis= pPis->fPlab.Eta();
  double etaPi = pPi->fPlab.Eta();
  double etaK  = pK->fPlab.Eta();

  // charge
  int qk   = pK->fQ;
  int qpi  = pPi->fQ;
  int qpis = pPis->fQ;

  // masses
  double mdstar = fpCand->fMass;
  double mdz = pC->fMass;
  double dm = mdstar - mdz;

  // D* vertex
  TAnaVertex sv = fpCand->fVtx;  // D* vertex
  TVector3 sv1 = sv.fPoint;  // D* decay vertex position
  int pvidx = (fpCand->fPvIdx > -1? fpCand->fPvIdx : 0);  // D* PV index
  TVector3 pv =  fpEvt->getPV(pvidx)->fPoint;  // Dstar vertex

  // D0 vertex
  TAnaVertex svD0 = pC->fVtx;  // D0 vertex
  TVector3 sv2 = svD0.fPoint; // D0 decay vertex position
  //int pvidx2 = (pC->fPvIdx > -1? pC->fPvIdx : 0);  // D0 PV index, DOES NOT HAVE A PV
  //TVector3 pv2 =  fpEvt->getPV(pvidx2)->fPoint;    // It is always -1

  //cout<< fpCand->fPvIdx <<" "<< pC->fPvIdx <<endl;

  if(fpCand->fPvIdx==-1 ) ((TH1D*)fHistDir->Get("Status"))->Fill(3.);


  //((TH1D*)fHistDir->Get("Status"))->Fill(4.);

  ((TH1D*)fHistDir->Get("pvidx"))->Fill(float(pvidx));

  // Calculate angles
  TVector3 t1(sv1-pv), t2(sv2-pv), t3(sv2-sv1);
  //cout<<t1.Z()<<" "<<t2.Z()<<" "<<t3.Z()<<endl;


  //fls3d = sv.fD3d/sv.fD3dE; // D*
  //flsxy = sv.fDxy/sv.fDxyE;
  fls3d = svD0.fD3d/svD0.fD3dE; //  use D0
  flsxy = svD0.fDxy/svD0.fDxyE;
  prob  = svD0.fProb;
  chi2  = svD0.fChi2;

  //alpha = t1.Angle(pCand->fPlab);  // D* pointing angle
  alpha  = t2.Angle(pC->fPlab); // D0 angle
  //falpha2 = t3.Angle(pC->fPlab); // D0 angle versus SV2-SV1
  dr = piSlowMom.Angle(fpCand->fPlab); // pislow openinig
  double dR = piSlowMom.Angle(pC->fPlab); // pislow openinig versus D0, this is what is used in HFDstar

  if(fVerbose>8 ) {
    cout<<" PVs "<<pvidx<<" "
	<<pv.X()<<" "<<pv.Y()<<" "<<pv.Z()<<" "
	<<sv1.X()<<" "<<sv1.Y()<<" "<<sv1.Z()<<" "
	<<sv2.X()<<" "<<sv2.Y()<<" "<<sv2.Z()<<" "
	<<t1.Mag()<<" "<<t2.Mag()<<" "<<t3.Mag()<<endl;
    cout<<"  "<<sv.fD3d<<" "<<sv.fDxy<<" "<<svD0.fD3d<<" "<<svD0.fDxy<<" "
	<<alpha<<" "<<t2.Angle(pC->fPlab)<<" "<<t3.Angle(pC->fPlab)<<" "<<t1.Angle(t3)<<endl;
  }


  ((TH1D*)fHistDir->Get("Status"))->Fill(10.);
  count1++;

  // Now the selection cuts cut
  if(fVerbose>8 ) cout<<"Check pre-cuts "<<endl;

  // Now histogram
  const bool doHisto = true;
  if(doHisto) {  //
    ((TH1D*)fHistDir->Get("mds"))->Fill(mdstar);
    ((TH1D*)fHistDir->Get("mdz"))->Fill(mdz);
    ((TH2D*)fHistDir->Get("h2d"))->Fill(mdz, dm);
    ((TH1D*)fHistDir->Get("dm"))->Fill(dm);
  }

  // Cuts always done
  // skip wrong sign decys
  int fVerbose0 = fVerbose;
  if( (qpi+qpis)==0 ) {if(fVerbose0>3) cout<<" failed qpi+qpis cut "<<qpi<<" "<<qpis<<endl; return;}
  ((TH1D*)fHistDir->Get("Status"))->Fill(11.);  //11
  // check sign
  if( (qk+qpi)!=0 ) {if(fVerbose0>3) cout<<" failed q cut "<<qpi<<" "<<qk<<endl; return;} // normal cuts 
  //if( (qk+qpi)==0 ) {if(fVerbose0>3) cout<<" failed q cut "<<qpi<<" "<<qk<<endl; return;} // reversed cuts 
  ((TH1D*)fHistDir->Get("Status"))->Fill(12.);
  // limit dm to +-160MeV
  if( dm<0.130 || dm>0.160 ) {if(fVerbose0>3) cout<<" failed dm cut "<<dm<<endl; return;}
  //if(fVerbose0>3 ) cout<<"Passed qcut "<<endl;
  ((TH1D*)fHistDir->Get("Status"))->Fill(13.);
  // Restablish mass cuts from MSSW
  if( mdz<1.76 || mdz>1.96 ) {if(fVerbose0>3) cout<<" failed mdz cut "<<mdz<<endl; return;}
  ((TH1D*)fHistDir->Get("Status"))->Fill(14.);
  if( mdstar<1.91 || mdstar>2.11 ) {if(fVerbose0>3) cout<<" failed mdstar cut "<<mdstar<<endl; return;}
  ((TH1D*)fHistDir->Get("Status"))->Fill(15.);

  if(fVerbose0>3 ) cout<<"Passed init-cuts "<<endl;

  // Now histogram
  if(doHisto) {  //

    ((TH1D*)fHistDir->Get("h1"))->Fill(t1.Mag());
    ((TH1D*)fHistDir->Get("h2"))->Fill(t2.Mag());
    ((TH1D*)fHistDir->Get("h3"))->Fill(t3.Mag());

    ((TH1D*)fHistDir->Get("h4"))->Fill( (pC->fPlab).Angle(fpCand->fPlab) );
    ((TH1D*)fHistDir->Get("h5"))->Fill(sv.fD3d);
    ((TH1D*)fHistDir->Get("h6"))->Fill(svD0.fD3d);

    ((TH1D*)fHistDir->Get("h7"))->Fill(t2.Angle(pC->fPlab));
    ((TH1D*)fHistDir->Get("h8"))->Fill(t3.Angle(pC->fPlab));
    ((TH1D*)fHistDir->Get("h9"))->Fill(t1.Angle(t3));
    ((TH1D*)fHistDir->Get("h10"))->Fill(t1.Angle(fpCand->fPlab));

    ((TH1D*)fHistDir->Get("fls3d"))->Fill(fls3d);
    ((TH1D*)fHistDir->Get("flsxy"))->Fill(flsxy);
    ((TH1D*)fHistDir->Get("prob"))->Fill(prob);
    ((TH1D*)fHistDir->Get("chi2"))->Fill(chi2);
    ((TH1D*)fHistDir->Get("alpha"))->Fill(alpha);
    ((TH1D*)fHistDir->Get("pt"))->Fill(pt);
    ((TH1D*)fHistDir->Get("dr"))->Fill(dr);
    ((TH1D*)fHistDir->Get("dR"))->Fill(dR);

    ((TH1D*)fHistDir->Get("ptdz"))->Fill(ptdz);
    ((TH1D*)fHistDir->Get("ptPis"))->Fill(ptPis);
    ((TH1D*)fHistDir->Get("ptPi"))->Fill(ptPi);
    ((TH1D*)fHistDir->Get("ptK"))->Fill(ptK);

    if(tm==1) ((TH1D*)fHistDir->Get("dm1"))->Fill(dm);
    if(mcOk) ((TH1D*)fHistDir->Get("dm2"))->Fill(dm);

    //if(tm==1) tcount1++;
  }  // if


  // Standard selection cuts
  if(1) { // skip it only for special testing
    //if (prob < 0.05) {if(fVerbose>3) cout<<" failed prob "<<prob<endl; return;}
    //((TH1D*)fHistDir->Get("Status"))->Fill(14.);

    if(ptPi<3.5 || ptK<3.5) {if(fVerbose0>3) cout<<" failed pi/k pt cut "<<ptPi<<" "<<ptK<<endl; return;}
    ((TH1D*)fHistDir->Get("Status"))->Fill(16.);

    if (ptPis < 0.4) {if(fVerbose0>3) cout<<" failed pt slow pt cut "<<ptPis<<endl; return;}
    ((TH1D*)fHistDir->Get("Status"))->Fill(17.);

    if (dr > 0.15) {if(fVerbose0>3) cout<<" failed dr cut "<<dr<<endl; return;}
    ((TH1D*)fHistDir->Get("Status"))->Fill(18.);

    if (chi2 > 2.0) {if(fVerbose0>3) cout<<" failed chis2 cut "<<chi2<<endl; return;}
    ((TH1D*)fHistDir->Get("Status"))->Fill(19.);

    if (pt < 5) {if(fVerbose0>3) cout<<" failed pt cut "<<pt<<endl; return;}
    ((TH1D*)fHistDir->Get("Status"))->Fill(20.);

    if (alpha > 0.3) {if(fVerbose0>3) cout<<" failed alpha cut "<<alpha<<endl; return;}
    ((TH1D*)fHistDir->Get("Status"))->Fill(21.);

    if (fls3d < 2) {if(fVerbose0>3) cout<<" failed fls3d cut "<<fls3d<<endl; return;}
    ((TH1D*)fHistDir->Get("Status"))->Fill(22.);

  } // skip for testing
  // no cuts after this step

  if(fVerbose0>3 ) cout<<"Passed pre-cuts "<<endl;
  count2++;

  // do some printing
  if(MYDEBUG) {
    cout<<" dstar: pi - "<< pPi->fIndex<<" "<<pPi->fMuID<<" "<<pPi->fMuIndex
	<<" K - "<<pK->fIndex<<" "<<pK->fMuID<<" "<<pK->fMuIndex
	<<endl;
    doTest(fpCand,40); // print trigger info
    doTest(fpCand,30); // print muon info
  }

  // check the trigger info
  hltInfo.clear(); // to save trigger info

  // check matching of trigger track to offline tracks
  double match1dr5 = 100., match2dr5 = 100.;
  bool passMatch = doTriggerMatchingForDs(match1dr5,match2dr5);
  ((TH1D*)fHistDir->Get("dr4"))->Fill(match1dr5);
  ((TH1D*)fHistDir->Get("dr4"))->Fill(match2dr5);
  // redefine the matching variable from CandAna
  fmatchTrigs = passMatch; // this will overwrite the value from CanDana
  fmatch1dr5 = match1dr5;
  fmatch2dr5 = match2dr5;

  // Look at muid
  int mid1 = 0, mid2= 0;
  if(MYDEBUG) cout<< " mu-index "<<pPi->fMuIndex <<" "<<pK->fMuIndex <<endl;
  if (pPi->fMuIndex > -1) mid1 = fpEvt->getMuon(pPi->fMuIndex)->fMuID;
  if (pK->fMuIndex > -1)  mid2 = fpEvt->getMuon(pK->fMuIndex)->fMuID;

  // all bits except ecal (CAL BITS ARE NEVER SET?)
  //bool muid11 = ( (mid1&0x7FFF) != 0);
  //bool muid12 = ( (mid2&0x7FFF) != 0);
  // global muons
  bool muid11 = ( (mid1&0x2) != 0);
  bool muid12 = ( (mid2&0x2) != 0);
  // global or tracker muons or standalone muons
  bool muid21 = ( (mid1&0x7) != 0);
  bool muid22 = ( (mid2&0x7) != 0);
  // global or tracker muons
  bool muid31 = ( (mid1&0x6) != 0);
  bool muid32 = ( (mid2&0x6) != 0);

  fb6=muid11; // global muons
  fb7=muid12;

  // tight muons
  bool muid1 = tightMuon(pPi);  // true for good/tight  muons
  bool muid2 = tightMuon(pK);
  // BDT muons
  fmuidmva1 = 0; // mvaMuon(pPi,fmva1);
  fmuidmva2 = 0; // mvaMuon(pK, fmva2);

  if(muid1) ((TH1D*)fHistDir->Get("Status"))->Fill(31.);
  if(muid2) ((TH1D*)fHistDir->Get("Status"))->Fill(32.);

  if(MYDEBUG) {
    cout<<" id 1: tig "<<muid1<<" mv "<<fmuidmva1<<" all "
	<<muid11<<" "<<muid21<<" "<<muid31<<" "<<mid1<<" "<<fmva1<<endl;
    cout<<" id 2: tig "<<muid2<<" mv "<<fmuidmva2<<" all "
	<<muid12<<" "<<muid22<<" "<<muid32<<" "<<mid2<<" "<<fmva2<<endl;
  }


  // use the matching function from candAna()
  //                               anyTrig muonOnly anyModule
  //fmatch1dr = doTriggerMatchingR(pPi,false,true,false);  // see if it matches HLT muon
  //fmatch2dr = doTriggerMatchingR(pK, false,true,false);  // see if it matches HLT muon

  // for testing only
  //                               anyTrig muonOnly
  //fmatch1dr1 = doTriggerMatchingR(pPi,false,false,false);  // see if it matches HLT muon
  //fmatch2dr1 = doTriggerMatchingR(pK, false,false,false);  // see if it matches HLT muon

  //fmatch1dr2 = doTriggerMatchingR(pPi,true,false,false);  // see if it matches HLT muon
  //fmatch2dr2 = doTriggerMatchingR(pK, true,false,false);  // see if it matches HLT muon

  //fmatch1dr3 = doTriggerMatchingR(pPi,true,true,false);  // see if it matches HLT muon
  //fmatch2dr3 = doTriggerMatchingR(pK, true,true,false);  // see if it matches HLT muon

  // use matching function from candAnaDstar() FOR TESTING
  //int idxt1=-1, idxt2=-1, idxt3=-1;
  //fmatch1dr4 = 0.; // doTriggerMatchingTest(idxt1,1,0);  // see which track  matches best HLT muon1
  //fmatch2dr4 = 0.; // doTriggerMatchingTest(idxt2,2,0);  // see which track  matches HLT muon2
  //double fmatch3dr4 = 0.; // doTriggerMatchingTest(idxt3,3,0);  // see if it matches HLT muon3, if exist
  //cout<<" match for trig muon3 "<<idxt3<<" "<<fmatch3dr4<<endl;

  //bool mumatch1 = (fmatch1dr<0.02); // see if it matches HLT muon
  //bool mumatch2 = (fmatch2dr<0.02); // see if it matches HLT muon
  //bool mumatch1 = (fmatch1dr<0.01); // see if it matches HLT muon
  //bool mumatch2 = (fmatch2dr<0.01); // see if it matches HLT muon
  //if(mumatch1) ((TH1D*)fHistDir->Get("Status"))->Fill(35.);
  //if(mumatch2) ((TH1D*)fHistDir->Get("Status"))->Fill(36.);

  //((TH1D*)fHistDir->Get("dr1"))->Fill(fmatch1dr);
  //((TH1D*)fHistDir->Get("dr1"))->Fill(fmatch2dr);
  //((TH1D*)fHistDir->Get("dr2"))->Fill(fmatch1dr1);
  //((TH1D*)fHistDir->Get("dr2"))->Fill(fmatch2dr1);
  //((TH1D*)fHistDir->Get("dr9"))->Fill(fmatch1dr2);
  //((TH1D*)fHistDir->Get("dr9"))->Fill(fmatch2dr2);
  //((TH1D*)fHistDir->Get("dr7"))->Fill(fmatch1dr3);
  //((TH1D*)fHistDir->Get("dr7"))->Fill(fmatch2dr3);
  //((TH1D*)fHistDir->Get("dr8"))->Fill(fmatch1dr4);
  //((TH1D*)fHistDir->Get("dr8"))->Fill(fmatch2dr4);


  // if(MYDEBUG) {
  //   //if(muid1 && (match1dr2 != match1dr4))
  //   cout<<" match 1 "<<fmatch1dr<<" "<<fmatch1dr1<<" "<<fmatch1dr2<<" "<<fmatch1dr3<<" "<<fmatch1dr4
  // 	<<" "<<mumatch1<<endl;
  //   //if(muid2&&(match2dr2 != match2dr4))
  //   cout<<" match 2 "<<fmatch2dr<<" "<<fmatch2dr1<<" "<<fmatch2dr2<<" "<<fmatch2dr3<<" "<<fmatch2dr4
  // 	<<" "<<mumatch2<<endl;
  // }

  // kink finder
  double chiPi = -99.; // pPi->fChi2;
  double chiK = -99.; // pK->fChi2;
  //int mid1 = 0, mid2= 0;
  if (pPi->fMuIndex > -1) {
    chiPi= fpEvt->getMuon(pPi->fMuIndex)->fMuonChi2;
    //mid1 = fpEvt->getMuon(pPi->fMuIndex)->fMuID;
  }
  if (pK->fMuIndex > -1)  {
    chiK = fpEvt->getMuon(pK->fMuIndex)->fMuonChi2;
    //mid2 = fpEvt->getMuon(pK->fMuIndex)->fMuID;
  }


  // Check the trigger matching
  // Make sure that there is a trigger which does not involve the pi & K
  const float drCut = 0.025; // was 0.02
  fveto=false;
  //                             muonsOnly matchPt allModules dR histoOffset
  bool fveto1 = doTriggerVeto(pPi, true,false,true,drCut,1); // pi
  bool fveto2 = doTriggerVeto(pK,  true,false,true,drCut,2); // K
  bool fveto3 = doTriggerVeto(pPis,true,false,true,drCut,3); // piSlow
  fveto = fveto1 || fveto2 || fveto3;

  // new veto
  fb4 = analyzeHltInfo(false);  // at least 1 clean trig
  fb5 = analyzeHltInfo(true);  // all triggers clean


#ifdef DO_TESTS
  int is1 = doTest(fpCand,42);  // number of triggers (right DS)
  int is21 = doTest(fpCand,51); // print trigger L1 obj info
  int is22 = doTest(fpCand,52); // print trigger L2 obj info
  int is23 = doTest(fpCand,53); // print trigger L3 obj info
  int is24 = doTest(fpCand,54); // print trigger L3 no-mu obj info
  int is31 = doTest(fpCand,61); // print trigger mu candidates L1
  int is32 = doTest(fpCand,62); // print trigger mu candidates L2
  int is33 = doTest(fpCand,63); // print trigger mu candidates L2

  int is4 = doTest(fpCand,71); // print hlt info vectors
  // is4 - number of triggers with no overlap with pi/K
  // hltInfo.size() number of all triggers
  // is4>0,  (hltInfo.size()
  int allTriggers = hltInfo.size();
  int badTriggers = allTriggers - is4;

  fitmp4 = allTriggers;
  fitmp2=is4; // number of good triggers
  fitmp3=-1;

  if(allTriggers>0) {
    fitmp3=badTriggers; // number of bad triggers, exclude 0 triggers case

    ((TH1D*)fHistDir->Get("htest20"))->Fill(float(is4));
    ((TH1D*)fHistDir->Get("htest21"))->Fill(float(badTriggers)/float(allTriggers));
    ((TH1D*)fHistDir->Get("htest16"))->Fill(float(allTriggers),float(is4));

    ((TH1D*)fHistDir->Get("htest0"))->Fill(float(badTriggers));
    ((TH1D*)fHistDir->Get("htest2"))->Fill(float(is1));
    ((TH1D*)fHistDir->Get("htest3"))->Fill(float(is21));
    ((TH1D*)fHistDir->Get("htest4"))->Fill(float(is22));
    ((TH1D*)fHistDir->Get("htest5"))->Fill(float(is23));
    ((TH1D*)fHistDir->Get("htest6"))->Fill(float(is24));
    ((TH1D*)fHistDir->Get("htest7"))->Fill(float(is33));
    ((TH1D*)fHistDir->Get("htest8"))->Fill(float(is32));
    ((TH1D*)fHistDir->Get("htest9"))->Fill(float(is31));

    ((TH1D*)fHistDir->Get("htest10"))->Fill(float(is1),float(is33));
    ((TH1D*)fHistDir->Get("htest11"))->Fill(float(is1),float(is32));
    ((TH1D*)fHistDir->Get("htest12"))->Fill(float(is1),float(is31));
    ((TH1D*)fHistDir->Get("htest13"))->Fill(float(is1),float(is23));
    ((TH1D*)fHistDir->Get("htest14"))->Fill(float(is1),float(is22));
    ((TH1D*)fHistDir->Get("htest15"))->Fill(float(is1),float(is21));

    if( fb4 != (fitmp2>0) ) cout<<" ERROR4 "<<fb4<<" "<<fitmp2<<endl;
    if( fb5 != (fitmp3==0) ) cout<<" ERROR5 "<<fb5<<" "<<fitmp3<<endl;

  }

  bool switchPrintout=false;
  if(1 &&  fitmp2>0 && fveto==1 ) {
    switchPrintout = true;
    cout<<" hlt info "<<hltInfo.size()<<" "<<is1
	<<" L1/2/3/noMu "<<is21<<"/"<<is22<<"/"<<is23<<"/"<<is24
	<<" L3/2/1 "<<is33<<"/"<<is32<<"/"<<is31<<" "<<is4<<endl;;
    cout<< fitmp4 << " "<< fitmp3 <<" "<< fitmp2 << endl;

    doTest(fpCand,40); // print trigger info
    doTest(fpCand,59); // print trigger obj info
    doTest(fpCand,60); // print trigger mu candidates

    cout<<" veto info "<<" "<<fveto<<" "<<fveto1<<" "<<fveto2<<" "<<fveto3<<endl;
    // loop over muons
    for(vector<int>::iterator iter=hltInfo.begin(); iter!=hltInfo.end(); iter++)
      cout<<hex<<*iter<<dec<<" "<<endl;

  }
#endif // DO_TESTS

  //bool fb11 = doTriggerVeto(pPi, false,true,true,drCut,0); //
  //bool fb12 = doTriggerVeto(pK,  false,true,true,drCut,0); //
  //bool fb13 = doTriggerVeto(pPis,false,true,true,drCut,0); //
  //fb1 = fb11 || fb12 || fb13;

  // include pt match
  bool fb21 = doTriggerVeto(pPi, true,true,true,drCut,0); //
  bool fb22 = doTriggerVeto(pK,  true,true,true,drCut,0); //
  bool fb23 = doTriggerVeto(pPis,true,true,true,drCut,0); //
  fb2 = fb21 || fb22 || fb23;

  // include non-muons
  bool fb31 = doTriggerVeto(pPi, false,false,true,drCut,0); //
  bool fb32 = doTriggerVeto(pK,  false,false,true,drCut,0); //
  bool fb33 = doTriggerVeto(pPis,false,false,true,drCut,0); //
  fb3 = fb31 || fb32 || fb33;

  // test drcut
  bool tm1 = doTriggerVeto(pPi, true,false,true,0.01,1); // pi
  bool tm2 = doTriggerVeto(pK,  true,false,true,0.01,2); // K
  bool tm3 = doTriggerVeto(pPis,true,false,true,0.01,0); // piSlow
  ftmp1 = tm1 || tm2 || tm3;

  tm1 = doTriggerVeto(pPi, true,false,true,0.02,1); // pi was 0.025
  tm2 = doTriggerVeto(pK,  true,false,true,0.02,2); // K
  tm3 = doTriggerVeto(pPis,true,false,true,0.02,0); // piSlow
  ftmp2 = tm1 || tm2 || tm3;

  tm1 = doTriggerVeto(pPi, true,false,true,0.025,1); // pi
  tm2 = doTriggerVeto(pK,  true,false,true,0.025,2); // K
  tm3 = doTriggerVeto(pPis,true,false,true,0.025,0); // piSlow
  ftmp3 = tm1 || tm2 || tm3;

  tm1 = doTriggerVeto(pPi, true,false,true,0.05,1); // pi
  tm2 = doTriggerVeto(pK,  true,false,true,0.05,2); // K
  tm3 = doTriggerVeto(pPis,true,false,true,0.05,0); // piSlow
  ftmp4 = tm1 || tm2 || tm3;

  tm1 = doTriggerVeto(pPi, true,false,true,0.075,1); // pi
  tm2 = doTriggerVeto(pK,  true,false,true,0.075,2); // K
  tm3 = doTriggerVeto(pPis,true,false,true,0.075,0); // piSlow
  ftmp4 = tm1 || tm2 || tm3;

  tm1 = doTriggerVeto(pPi, true,false,true,0.1,1); // pi
  tm2 = doTriggerVeto(pK,  true,false,true,0.1,2); // K
  tm3 = doTriggerVeto(pPis,true,false,true,0.1,0); // piSlow
  ftmp6 = tm1 || tm2 || tm3;

  tm1 = doTriggerVeto(pPi, true,false,true,0.15,1); // pi
  tm2 = doTriggerVeto(pK,  true,false,true,0.15,2); // K
  tm3 = doTriggerVeto(pPis,true,false,true,0.15,0); // piSlow
  ftmp7 = tm1 || tm2 || tm3;

  tm1 = doTriggerVeto(pPi, true,false,true,0.2,1); // pi
  tm2 = doTriggerVeto(pK,  true,false,true,0.2,2); // K
  tm3 = doTriggerVeto(pPis,true,false,true,0.2,0); // piSlow
  ftmp8 = tm1 || tm2 || tm3;

  tm1 = doTriggerVeto(pPi, true,false,true,0.3,1); // pi
  tm2 = doTriggerVeto(pK,  true,false,true,0.3,2); // K
  tm3 = doTriggerVeto(pPis,true,false,true,0.3,0); // piSlow
  ftmp9 = tm1 || tm2 || tm3;


  if(MYDEBUG)
    cout<<" veto "<<fveto<<" "<<fb2<<" "<<fb3<<endl;

#ifdef DO_TESTS
  // harder cuts, like final, for testing only
  if(0) { // enable for printout of very selected events

    // very tight, peak only
    if( dm<0.143 || dm>0.148 ) {if(fVerbose>3) cout<<" failed dm cut "<<dm<<endl; return;}  //
    // wide, also background
    //if( dm<0.135 || dm>0.155 ) {if(fVerbose>3) cout<<" failed dm cut "<<dm<<endl; return;}
    if( mdz<1.82 || mdz>1.91 ) {if(fVerbose>3) cout<<" failed mdz cut "<<mdz<<endl; return;}
    if( mdstar<1.97 || mdstar>2.05 ) {if(fVerbose>3) cout<<" failed mdstar cut "<<mdstar<<endl; return;}
    if(ptPi<4.0 || ptK<4.0) {if(fVerbose>3) cout<<" failed pi/k pt cut "<<ptPi<<" "<<ptK<<endl; return;}
    if (ptPis < 0.5) {if(fVerbose>3) cout<<" failed pt slow pt cut "<<ptPis<<endl; return;}
    if (dr > 0.08) {if(fVerbose>3) cout<<" failed dr cut "<<dr<<endl; return;}
    if (chi2 > 2.0) {if(fVerbose>3) cout<<" failed chis2 cut "<<chi2<<endl; return;}
    if (pt < 6.0) {if(fVerbose>3) cout<<" failed pt cut "<<pt<<endl; return;}
    if (alpha > 0.15) {if(fVerbose>3) cout<<" failed alpha cut "<<alpha<<endl; return;}
    if (fls3d < 2.0) {if(fVerbose>3) cout<<" failed fls3d cut "<<fls3d<<endl; return;}

    // json cut
    if(!fJSON) return; // check jason
    // triggered
    if(!fGoodHLT) return;
    // 1 tigger
    //if(fhltType>=1000) return;  // veto multiple triggers
    // veto
    if(fveto) return;  // veto trigger matched pi

    // muid
    if(!(muid1||muid2) ) return; // select misidentified pi->mu & K->mu

    switchPrintout=true;
  } // if test

  // get the trigger type
  //int itest = doTest(fpCand,41); // get the trigger type
  //if(itest<=0) cout<<" ERROR: not trigger "<<itest<<endl;
  //int dum=0;
  //fitmp1 = itest;

  //if( (itest==1) && ((fmatch1dr4>0.1) || (fmatch2dr4>0.1)) ) switchPrintout=true;
  //if( fGoodHLT && fmatchTrigs && (fhltType<1000) && ((fmatch1dr4>0.1) || (fmatch2dr4>0.1))) switchPrintout=true;

  if(switchPrintout) {  // if print candidates
    dumpAll();

    //static int ic=0;
    //ic++;
    cout<<"Dstar candidate "<<count0<<" dm "<<dm<<" muid1 "<<muid1<<" muid2 "<<muid2
	<<" index pi/K/piSlow "<<piIndex<<" "<<KIndex<<" "<<piSlowIndex<<endl;

    //cout<<fHLTPath<<" "<<fhltType<<" "<<fGoodHLT<<endl;

    //cout<<" json "<<fJSON<<" hlttype "<<fhltType<<" dm "<<dm<<" fb3 "<<fb3<<" muid1 "<<muid1
    //  <<" muid2 "<<muid2<<" count "<<ic<<endl;
    //cout<<" Dstar candidate "<<fpCand->fType<<" in event "<<fEvt<<" run "<<fRun;
    //cout << " with mass = " << fpCand->fMass <<" cand num "<<count0<<endl;
    //cout<<" veto "<<fveto<<" "<<fb1<<" "<<fb2<<" "<<fb3<<endl;
    // cout<<" matching "
    // 	<<fmatch1dr<<"/"<<fmatch2dr<<" "
    // 	<<fmatch1dr1<<"/"<<fmatch2dr1<<" "
    // 	<<fmatch1dr2<<"/"<<fmatch2dr2<<" "
    // 	<<fmatch1dr3<<"/"<<fmatch2dr3<<" "
    // 	<<mumatch1<<" "<<mumatch2<<endl;
    //cout<<" trig match to all tracks "<<fmatch1dr4<<"/"<<fmatch2dr4<<"/"<<fmatch3dr4<<" "
    //	<<idxt1<<"/"<<idxt2<<"/"<<idxt3<<endl;


    //cout << "DUMP HFDstarCandidate  " <<endl;
    //doTest(fpCand,0);    // print all
    //dumpHFDstarCand(fpCand);

    int fVerbose0 = fVerbose;
    fVerbose=100;
    cout<<" PI veto "<<endl;
    bool tmp  = doTriggerVeto(pPi, true,false,true,0.025,0); // pi
    cout<< " veto = "<<tmp<<endl;
    cout<<" K veto "<<endl;
    tmp = doTriggerVeto(pK,  true,false,true,0.025,0); // K
    cout<< " veto = "<<tmp<<endl;
    cout<<" PI-slow veto "<<endl;
    tmp = doTriggerVeto(pPis,true,false,true,0.025,0); // piSlow
    cout<< " veto = "<<tmp<<endl;


    ////float r1 = doTriggerMatchingR(pPi,false,true,false);  // see if it matches HLT muon
    ////float r2 = doTriggerMatchingR(pK, false,true,false);  // see if it matches HLT muon
    //float r1 = doTriggerMatchingR(pPi,true,true,false);  // see if it matches HLT muon
    //float r2 = doTriggerMatchingR(pK, true,true,false);  // see if it matches HLT muon
    //cout<<tmp<<" "<<r1<<" "<<r2<<endl;

    fVerbose=fVerbose0;

  } // if special printout

#endif // DO_TESTS


  //  match to offline muons
  double dr1 = matchToMuon(pPi,true); // skip same track muons
  double dr2 = matchToMuon(pK,true);
  double dr11 = matchToMuon(pPi,false); // do not skip same track muons
  double dr12 = matchToMuon(pK,false);

  if(MYDEBUG)
    cout<<" match to offline muons  "<<dr1<<"/"<<dr11<<" "<<dr2<<"/"<<dr12<<endl;
  if(muid1) ((TH1D*)fHistDir->Get("dr5"))->Fill(dr1);
  if(muid2) ((TH1D*)fHistDir->Get("dr5"))->Fill(dr2);
  ((TH1D*)fHistDir->Get("dr6"))->Fill(dr2);
  ((TH1D*)fHistDir->Get("dr6"))->Fill(dr1);
  //((TH1D*)fHistDir->Get("dr10"))->Fill(dr12);
  //((TH1D*)fHistDir->Get("dr10"))->Fill(dr11);

  if(dr1<0.02)  ((TH1D*)fHistDir->Get("Status"))->Fill(41.);
  if(dr2<0.02)  ((TH1D*)fHistDir->Get("Status"))->Fill(42.);
  //if(dr11<0.02)  ((TH1D*)fHistDir->Get("Status"))->Fill(43.);
  //if(dr12<0.02)  ((TH1D*)fHistDir->Get("Status"))->Fill(44.);

  // match offline muons to muon triggers, save in a vector
  int numHltMuon = doMuonTriggerMatching();
  if(MYDEBUG) cout<<" doMuonTriggerMatching "<<numHltMuon<<endl;

  // match tracks to trigger matched muons, IS IT STILL USEFULL?
  int idx1=-1, idx2=-1;
  // needs doMuonTriggerMatching to be run 1st
  bool matchToTrigMu1 = matchToTriggeredMuon(pPi, idx1);
  bool matchToTrigMu2 = matchToTriggeredMuon(pK, idx2);
  if(MYDEBUG) cout<<" matchToTriggeredMuon "<<matchToTrigMu1<<" "<<idx1<<" "
		  <<matchToTrigMu2<<" "<<idx2<<endl;
  if(matchToTrigMu1)  ((TH1D*)fHistDir->Get("Status"))->Fill(45.);
  if(matchToTrigMu2)  ((TH1D*)fHistDir->Get("Status"))->Fill(46.);

  // Match with a Jpsi in this event   IS IT USEFULL?
  //idx1=-1; idx2=-1;
  //bool foundJpsi = getJpsi(idx1, idx2); // needs doMuonTriggerMatching to be run 1st
  //bool RejectPion = ( (idx1==pPi->fIndex) || (idx2==pPi->fIndex) );
  // bool RejectKaon = ( (idx1==pK->fIndex)  || (idx2==pK->fIndex) );
  //if(MYDEBUG) cout<<" found jpis "<<foundJpsi<<" "<<RejectPion<<" "<<RejectKaon<<endl;
  //if(foundJpsi)  ((TH1D*)fHistDir->Get("Status"))->Fill(47.);
  //if(RejectPion)  ((TH1D*)fHistDir->Get("Status"))->Fill(48.);
  //if(RejectKaon)  ((TH1D*)fHistDir->Get("Status"))->Fill(49.);

  // Isolation (something has changed in nCloseTracks)
  //                       dcaCut(cm) ptCut(GeV)
  //int close1 = nCloseTracks(fpCand,0.03, 0.5); // around D*
  int close2 = 0; // nCloseTracks(pC,    0.03, 0.5); // around D0
  //                                      dca   R    Pt
  //double iso1 = isoClassicWithDOCA(fpCand, 0.05,0.7, 0.9); // D*
  double iso2 = isoClassicWithDOCA(pC,     0.05,0.7, 0.9); // D0

  //if(tm==1) tcount2++;

  // Select event for the redtree
  fPreselection = true;  // select this event for the standrad redtree
  if(MYDEBUG) cout<<" preselection "<<fPreselection<<" "<<fGoodHLT<<endl;
  //
  fCandTM = tm; // candidate matched to MC
  ftm= int(mcOk); // candidate exist in MC
  fmds=mdstar;
  fmdz=mdz;
  fchi2=chi2;
  falpha=alpha;
  ffls3d=fls3d;
  fqpis=qpis;
  fdr=dr;

  feta=eta;
  fetapi=etaPi;
  fetak=etaK;

  fpt=pt;
  fptdz=ptdz;
  fptpis=ptPis;
  fptpi=ptPi;
  fptk=ptK;

  fmuid1 = muid1;
  fmuid2 = muid2;
  fmumat1 = 0; // mumatch1;
  fmumat2 = 0; //mumatch2;
  fmudr1 = dr1;
  fmudr2 = dr2;

  fchipi = chiPi;
  fchik = chiK;
  fiso = iso2;
  fnclose = close2;

  //cout<<count0<<" "<<count1<<" "<<count2<<" "<<count3<<" "<<count4<<" "<<count5<<" "<<count6<<endl;
  if(fVerbose>0) {if(count0%10 == 0) cout<<count0<<" "<<count1<<" "<<count2<<" "<<count3<<" "<<count4<<" "<<count5<<" "<<count6<<endl;}
  else           {if(count0%100 == 0) cout<<count0<<" "<<count1<<" "<<count2<<" "<<count3<<" "<<count4<<" "<<count5<<" "<<count6<<endl;}

}

// ----------------------------------------------------------------------
// dump all information about the event
void candAnaDstar::candEvaluation() {

  if(MYDEBUG) cout<<" preselection "<<fPreselection<<" "<<fGoodHLT<<endl;
  // make some more preselection tests if needed

}
// ----------------------------------------------------------------------
// dump all information about the event
void candAnaDstar::dumpAll() {
    static int ic=0;
    ic++;
    cout<<"Dstar in event "<<fEvt<<" run "<<fRun<<" json "<<fJSON<<" count "<<ic<<endl;
    cout<<" HLT info "<<fHLTPath<<" "<<fhltType<<" "<<fGoodHLT<<" "<<fmatchTrigs<<endl;

    //cout<<" Dstar candidate "<<fpCand->fType<<" in event "<<fEvt<<" run "<<fRun;
    cout << " with mass = " << fpCand->fMass;
    cout<<" veto "<<fveto<<" "<<fb2<<" "<<fb3<<endl;
    //cout<<" muid "<<muid1<<"/"<<muid2<<" "
    //	<<fmatch1dr<<"/"<<fmatch2dr<<" "
    //	<<fmatch1dr1<<"/"<<fmatch2dr1<<" "
    //	<<fmatch1dr2<<"/"<<fmatch2dr2<<" "
    //	<<fmatch1dr3<<"/"<<fmatch2dr3<<" "
    //	<<mumatch1<<" "<<mumatch2<<endl;
    //cout<<" trig match to all tracks "<<fmatch1dr4<<"/"<<fmatch2dr4<<"/"<<fmatch3dr4<<" "
    //	<<idxt1<<"/"<<idxt2<<"/"<<idxt3<<endl;
    //cout << "DUMP HFDstarCandidate  " <<endl;
    doTest(fpCand,0);    // print all
    dumpHFDstarCand(fpCand);
}

// ----------------------------------------------------------------------
// Loop over all trigger confirmed muons, match to track
bool candAnaDstar::matchToTriggeredMuon(TAnaTrack *pt, int &idx) {
  const bool PRINT = false;
  double dRMin = 99.;

  //cout<<"matched muons "<<hltMatchedMuons.size()<<endl;
  if(hltMatchedMuons.size() <2) return 0;

  TVector3 track = pt->fPlab;  // test track momentum
  int it0 = pt->fIndex;
  if(PRINT) {
    cout<<" check track "<<it0
	<< "pt,eta,phi: " << pt->fPlab.Perp() << " " << pt->fPlab.Eta() << " " << pt->fPlab.Phi() << endl;
  }

  int select=-1;
  // loop over muons
  for(vector<int>::iterator iter=hltMatchedMuons.begin(); iter!=hltMatchedMuons.end(); iter++) {
    //cout<<*iter<<" "<<endl;
    TSimpleTrack *pTrack = fpEvt->getSimpleTrack(*iter);
    //cout<<*iter<<" "<<pTrack<<endl;
    if(pTrack == 0) continue;

    //pTrack->dump();

    TVector3 mom = pTrack->getP(); // momentum
    int index = pTrack->getIndex(); // same as itrk
    int q     = pTrack->getCharge();
    int qual  = pTrack->getHighPurity();
    int muonId= pTrack->getMuonID(); // muon id
    int pvidx = pTrack->getPvIndex();
    int genidx=pTrack->getGenIndex();
    int inds  = pTrack->getIndices();
    int bits  = pTrack->getBits();
    double pt = mom.Perp();

    if(PRINT)
      cout<<" track "<<*iter<<" idx "<<index<<" "<<q<<" "<<qual<<" muon "<<muonId<<" "<<pvidx<<" gen "<<genidx<<" "
	  <<hex<<inds<<" "<<bits<<dec<<" "<<pt<<" "<<mom.Phi()<<" "<<mom.Eta()<<endl;

    double dR = mom.DeltaR(track);
    if(PRINT) cout<<index<<" "<<pt<<" "<<dR<<endl;
    if(dR<dRMin) {dRMin=dR; select=index;} // select the best fit

  } // for loop

  if(PRINT && select==it0)
    cout<<" the matched muon has the same track id as the track under test "<<it0<<endl;

  bool match = (dRMin<0.02);
  if(match) idx = select;

  return match;
}
// ----------------------------------------------------------------------
// Loop over all trigger confirmed muons
int candAnaDstar::getJpsi(int &idx1, int &idx2) {
  const bool PRINT = false;

  //cout<<"matched muons "<<hltMatchedMuons.size()<<endl;
  if(hltMatchedMuons.size() <2) return 0;

  double m0 = 9999.;
  int num=0;
  for(vector<int>::iterator iter=hltMatchedMuons.begin(); iter!=hltMatchedMuons.end(); iter++) {
    //cout<<*iter<<" "<<endl;
    TSimpleTrack *pTrack = fpEvt->getSimpleTrack(*iter);
    //cout<<*iter<<" "<<pTrack<<endl;
    if(pTrack == 0) continue;

    //pTrack->dump();

    TVector3 mom = pTrack->getP(); // momentum
    int index = pTrack->getIndex(); // same as itrk
    int q     = pTrack->getCharge();
    int qual  = pTrack->getHighPurity();
    int muonId= pTrack->getMuonID(); // muon id
    int pvidx = pTrack->getPvIndex();
    int genidx=pTrack->getGenIndex();
    int inds  = pTrack->getIndices();
    int bits  = pTrack->getBits();
    double pt = mom.Perp();

    if(PRINT)
      cout<<" track "<<*iter<<" idx "<<index<<" "<<q<<" "<<qual<<" muon "<<muonId<<" "<<pvidx<<" gen "<<genidx<<" "
	  <<hex<<inds<<" "<<bits<<dec<<" "<<pt<<" "<<mom.Phi()<<" "<<mom.Eta()<<endl;


    for(vector<int>::iterator iter2=iter+1; iter2!=hltMatchedMuons.end(); iter2++) {
      //cout<<*iter2<<" "<<endl;
      TSimpleTrack *pTrack2 = fpEvt->getSimpleTrack(*iter2);
      //cout<<*iter2<<" "<<pTrack2<<endl;
      if(pTrack2 == 0) continue;

      TVector3 mom2 = pTrack2->getP(); // momentum
      //double pt2 = mom.Perp();

      double m = danekUtils::twoBodyDecayMass(mom, mom2, MMUON, MMUON);
      if( m>2.99 && m<3.19 ) {
	num++;
	if( abs(m-3.0969) < abs(m0-3.0969) ) {m0=m; idx1=*iter; idx2=*iter2;}
      }

      //cout<<*iter2<<" "<<pt2<<" "<<m<<" "<<m0<<endl;

    }
    // AnaTrack access does not work
    //TAnaTrack *rTrack = fpEvt->getRecTrack(*iter);
    //cout<<rTrack<<" "<<fpEvt->nRecTracks()<<endl;
    //if(rTrack == 0) continue;
    //mom = rTrack->fPlab; // momentum
    //pt = mom.Perp();
    //cout<<pt<<endl;

  } // for loop

  return num;
}
//---------------------------------------------------
// To analyze the MC event
// Looks for a clean Dstar->pi + D0->pi+K (clean means not other particles, e.g. gammas)
// Exits with True whenever at least one clean Dstar is found.
bool candAnaDstar::anaMC(bool print) {
  if(MYDEBUG) print = true;
  
#ifdef MC_HISTOS
  fmcmds=-1;
  fmcmdz=-1;
  fmcpt=-1;
  fmcptdz=-1;
  fmcptpis=-1;
  fmcptpi=-1;
  fmcptk=-1;
#endif

  //int numGenCands = fpEvt->nGenT();
  int numGenCands = fpEvt->nGenCands();
  if(print) cout<<" candAnaDstar::anaMC found gen cands "<<numGenCands<<endl;
  if(numGenCands<=0) return (false);

  TGenCand *pCand=0;
  bool foundDs = false, foundDz = false, foundPiSlow = false, foundPi=false, foundK=false, foundPV=false;
  int indexDs=-1, indexDz=-1, indexSlowPi=-1, indexPi=-1, indexK=-1;
  int qds =0, qk=0, qpi=0, qpis=0;
  int pC0 = 0;

  for (int it = 0; it < numGenCands; ++it) {
    //pCand = fpEvt->getGenT(it);
    pCand = fpEvt->getGenCand(it);
    //if(print) pCand->dump();   
    //if(print) cout <<it<< " "<<pCand->fID<<endl;

    //<<pCand->fQ<<" "<<pCand->fStatus<<" "
    //		   <<pCand->fMom1<<" "<<pCand->fMom2<<" "<<pCand->fDau1<<" "<<pCand->fDau2<<" "

    //if (TRUTHCAND == TMath::Abs(pC->fID)) {
    //for (int id = pB->fDau1; id <= pB->fDau2; ++id) {
    //pC = fpEvt->getGenTWithIndex(id);
    
    foundDs = false, foundDz = false, foundPiSlow = false, foundPi=false, foundK=false;

    if( !foundPV && ( abs(pCand->fID) == 5 || abs(pCand->fID) == 4  ) )
      {foundPV=true; PV=pCand->fV; if(print) cout<<" PV "<<PV.Z()<<endl;}   // get PV

    if( ( abs(pCand->fID) != 413) ) continue;        // look for Dstar,  skip others

    //if(print) cout <<it<< " "<<pCand->fID<<endl;
    if(print) cout <<" DS "<<it<< " " 
		   << pCand->fNumber << " "<<pCand->fID<<" "<<pCand->fQ<<" "<<pCand->fStatus<<" "
		   <<pCand->fMom1<<" "<<pCand->fMom2<<" "<<pCand->fDau1<<" "
		   <<pCand->fDau2<<" "<<pCand->fP.Perp()<<" "<<pCand->fV.Z()<<" "
		   <<pCand->fTag<<" "<<pCand->fMass<<" "<<pCand->fTime<<" "<<endl;

    qds = pCand->fQ;
    DSVertex = (pCand->fV);
    DSMom = (pCand->fP.Vect());
    foundDs = true;
    indexDs=it;
    pC0 = it;

    //int i1 = (pCand->fDau2)-(pCand->fDau1)+1;
    //if(i1!=2) {cout<<" number of daughters1 "<<i1<<endl;}
    //if(i1!=2) {continue;} // fpEvt->dumpGenBlock();}

    for(int id=(pCand->fDau1);id<=(pCand->fDau2);++id) { // check daughters
      //TGenCand *dau = fpEvt->getGenT(id);
      //TGenCand *dau = fpEvt->getGenTWithIndex(id);
      TGenCand *dau = fpEvt->getGenCand(id);
      if(print) cout<<" index "<<id<<endl;

      if( abs(dau->fID) == 421 ) { //  D0

	foundDz=true;
	indexDz=id;
	if(print) cout <<" D0 "<<dau->fNumber << " "<<dau->fID<<" "<<dau->fQ<<" "
		       <<dau->fMom1<<" "<<dau->fMom2<<" "<<dau->fDau1<<" "
		       <<dau->fDau2<<" "<<dau->fP.Perp()<<" "<<dau->fV.Z()<<endl;
	//TVector3 v1 = dau->fP.Vect();
	DZMom = (dau->fP.Vect());

	// D0 should have 2 daughters pi&K
	// Sometimes there is also a gamma in addition, reject it
	int i2 = (dau->fDau2)-(dau->fDau1)+1;
	if(i2!=2) 
	  {if(print)cout<<" number of D0 daughters not 2: "<<i2<<" skip"<<endl;continue;}

	for(int igd=(dau->fDau1);igd<=(dau->fDau2);++igd) { // check grand-daughters
	  //TGenCand *gdau = fpEvt->getGenTWithIndex(igd);
	  //TGenCand * gdau = fpEvt->getGenT(igd);
	  TGenCand *gdau = fpEvt->getGenCand(igd);
	  
	  if(print) cout<<" index "<<igd<<" "<<id<<endl;

	  //TVector3 v2 = gdau->fP.Vect();
	  if( abs(gdau->fID) == 321) {  // kaon
	    foundK = true;
	    indexK=igd;
	    id++;
	    if(print) cout<<" index "<<igd<<" "<<id<<endl;

	    KMom = (gdau->fP.Vect());
	    //float pt = v2.Perp();
	    //float eta = v2.Eta();
	    //float phi = v2.Phi();
	    if(print) cout <<" K "<<gdau->fNumber << " "<<gdau->fID<<" "<<gdau->fQ<<" "
		 <<gdau->fMom1<<" "<<gdau->fMom2<<" "
		 <<gdau->fP.Perp()<<" "<<gdau->fV.Z()<<endl;

	    DZVertex = (gdau->fV);
	    qk = gdau->fQ;

	  } else if( abs(gdau->fID) == 211) {  // pion
	    if(print) cout<<" index "<<igd<<" "<<id<<endl;
	    foundPi = true;
	    indexPi=igd;
	    id++;
	    PiMom = (gdau->fP.Vect());
	    qpi = gdau->fQ;
	    //float pt = v2.Perp();
	    //float eta = v2.Eta();
	    //float phi = v2.Phi();
	    if(print) cout <<" Pi "<<gdau->fNumber << " "<<gdau->fID<<" "<<gdau->fQ<<" "
			   <<gdau->fMom1<<" "<<gdau->fMom2<<" "
			   <<gdau->fP.Perp()<<" "<<gdau->fV.Z()<<endl;

	  } //

	  if( foundPi && foundK) break;
	} // end granddaughter loop

      } else if( (abs(dau->fID)) == 211) { // slow pion

	if(print) cout<<" index "<<id<<endl;
	foundPiSlow = true;
	indexSlowPi=id;
	if(print) cout <<" Slow pi "<<dau->fNumber << " "<<dau->fID<<" "<<dau->fQ<<" "
		       <<dau->fMom1<<" "<<dau->fMom2<<" "<<dau->fDau1<<" "
		       <<dau->fDau2<<" "<<dau->fP.Perp()<<" "<<dau->fV.Z()<<endl;
	PiSlowMom = dau->fP.Vect();
	qpis = dau->fQ;

      }  //

      if(foundPiSlow && foundDz) break;

    } // daugther loop

    // exit if we find the right candidate
    //cout<<" ok "<<foundPV<<" "<<foundDs<<" "<<foundDz<<" "<<foundPiSlow<<" "<<foundPi<<" "<<foundK<<" "<<fEvt<<endl;
    if(foundPV && foundDs && foundDz && foundPiSlow && foundPi && foundK) break;

  } // gen part loop

  bool ok = foundPV && foundDs && foundDz && foundPiSlow && foundPi && foundK;

  if(!ok) { 
    if(print) cout<<" not ok "<<foundPV<<" "<<foundDs<<" "<<foundDz<<" "<<foundPiSlow<<" "<<foundPi<<" "<<foundK<<endl;

  } else {

    if(print) cout<<" found Dstar "<<indexDs<<" "<<indexDz<<" "<<indexK<<" "<<indexPi<<" "<<indexSlowPi<<endl;

    if( (qds != qpis) || (qds != -qk) || ( qpis != -qk) || (qpis != qpi) ) {
      cout<<pC0<<" wrong charge in MC Ds,K,pi,pi_slow ";
      cout<<qds<<" "<<qk<<" "<<qpi<<" "<<qpis<<endl;
    }

#ifdef MC_HISTOS
    int tmp = qpi + qpis;
    ((TH1D*)fHistDir->Get("h300"))->Fill(float(tmp));

    //TVector3 t1(DSVertex-PV), t2(DZVertex-PV), t3(DZVertex-DSVertex);
    //double a1 = t1.Angle(DSMom);  // D* pointing angle
    //double a2 = t2.Angle(DZMom);  // D0 pointing angle
    //double a3 = t3.Angle(DZMom);  // D0 pointing angle with respect PV
    //double a4 = t1.Angle(t3);     // SV1 versus SV2
//     ((TH1D*)fHistDir->Get("h311"))->Fill(t1.Mag());
//     ((TH1D*)fHistDir->Get("h312"))->Fill(t2.Mag());
//     ((TH1D*)fHistDir->Get("h313"))->Fill(t3.Mag());
//     ((TH1D*)fHistDir->Get("h320"))->Fill(a1);
//     ((TH1D*)fHistDir->Get("h317"))->Fill(a2);
//     ((TH1D*)fHistDir->Get("h318"))->Fill(a3);
//     ((TH1D*)fHistDir->Get("h319"))->Fill(a4);
//     ((TH1D*)fHistDir->Get("h314"))->Fill(PiSlowMom.Perp());
//     ((TH1D*)fHistDir->Get("h315"))->Fill(PiSlowMom.Angle(DZMom));
//     ((TH1D*)fHistDir->Get("h316"))->Fill(PiSlowMom.Angle(DSMom));
//     ((TH1D*)fHistDir->Get("h321"))->Fill(DZMom.Angle(DSMom));
//     //((TH1D*)fHistDir->Get("h323"))->Fill(t3.Angle(DSMom));
//     ((TH1D*)fHistDir->Get("h325"))->Fill(t2.Angle(DSMom));
//     ((TH1D*)fHistDir->Get("h328"))->Fill(DSMom.Perp());
//     ((TH1D*)fHistDir->Get("h329"))->Fill(DZMom.Perp());
//     ((TH1D*)fHistDir->Get("h330"))->Fill(KMom.Perp());
//     ((TH1D*)fHistDir->Get("h330"))->Fill(PiMom.Perp());
    fmcpt=DSMom.Perp();
    fmcptdz=DZMom.Perp();
    fmcptpis=PiSlowMom.Perp();
    fmcptpi=PiMom.Perp();
    fmcptk=KMom.Perp();
    //double angle = danekUtils::twoBodyDecayAngle(KMom, PiMom);
    //double pt    = danekUtils::twoBodyDecayMomPerp(KMom, PiMom);
    double m1    = danekUtils::twoBodyDecayMass(KMom, PiMom, MKAON, MPION);
    double m2    = danekUtils::twoBodyDecayMass(KMom, PiMom, MPION, MKAON);
    fmcmdz=m1;
    //TVector3 t4(KMom+PiMom);
//     ((TH1D*)fHistDir->Get("h341"))->Fill(pt);
//     ((TH1D*)fHistDir->Get("h343"))->Fill(t4.Perp());
    ((TH1D*)fHistDir->Get("h342"))->Fill(m2); //broad
    ((TH1D*)fHistDir->Get("h344"))->Fill(m1); // narrow
//     tmp = DSMom.Perp();
//     ((TH2D*)fHistDir->Get("h351"))->Fill(tmp,m2);
//     ((TH2D*)fHistDir->Get("h352"))->Fill(angle,m2);
//     ((TH2D*)fHistDir->Get("h353"))->Fill(pt,m2);
//     tmp = KMom.Perp();
//     ((TH2D*)fHistDir->Get("h354"))->Fill(tmp,m2);
//     tmp = PiMom.Perp();
//     ((TH2D*)fHistDir->Get("h355"))->Fill(tmp,m2);
//     tmp = PiSlowMom.Perp();
//     ((TH2D*)fHistDir->Get("h356"))->Fill(tmp,m2);
//     tmp = PiSlowMom.Angle(DZMom);
//     ((TH2D*)fHistDir->Get("h357"))->Fill(tmp,m2);
     double m21    = danekUtils::twoBodyDecayMass(PiSlowMom, DZMom, MPION, m1);
     double m22    = danekUtils::twoBodyDecayMass(PiSlowMom, DZMom, MPION, m2);
     fmcmds=m21;
//     double angle2 = danekUtils::twoBodyDecayAngle(DZMom, PiSlowMom);
//     double pt2    = danekUtils::twoBodyDecayMomPerp(DZMom, PiSlowMom);
     ((TH1D*)fHistDir->Get("h345"))->Fill(m21); // narrrow
     ((TH1D*)fHistDir->Get("h346"))->Fill(m22); // wide 
//     ((TH1D*)fHistDir->Get("h347"))->Fill(angle2);
//     ((TH1D*)fHistDir->Get("h348"))->Fill(pt2);
//     ((TH1D*)fHistDir->Get("h349"))->Fill(m21-m1);
//     ((TH1D*)fHistDir->Get("h350"))->Fill(m22-m2);
#endif // MC_HISTOS
  }
  return ok;
}


// ----------------------------------------------------------------------
// Missused for various testing and debugging
void candAnaDstar::dumpHFTruthCand(TAnaCand *pC) {
  // misused 
  fpEvt->dumpGenBlock();

}

//----------------------------------------------------------------------
// Loop over all offline muons, select the ones which match trigger-muon
// put them in a seperate list : hltMatchedMuons
// return the number of trigger matched offline muons
// argument mode is unused for the moment
int candAnaDstar::doMuonTriggerMatching(void) {
  bool PRINT = false;
  //if(verbose>0) PRINT = true;

  // Check muons
  if(PRINT) cout<<"List muon tracks "<< fpEvt->nMuons()<<endl;
  hltMatchedMuons.clear(); // to save offline muons which match hlt muons

  const float dRMin = 0.02;
  int muHlt = 0, muTightHlt = 0;
  for (int it = 0; it< fpEvt->nMuons(); ++it) { // loop over muons
    TAnaMuon *muon = fpEvt->getMuon(it); // get muon

    // some direct muon methods
    //muon->dump();

    // some direct muon methods
    int muonId = muon->fMuID; // muon ID bits
    int muonIdx = muon->fMuIndex;  // index in the muon list = it
    int genidx = muon->fGenIndex; // gen index

    // get track index
    int itrk = muon->fIndex; // index of the SimpleTrack

    TVector3 muonMom = muon->fPlab;
    double ptMuon  = muonMom.Perp();

    ((TH1D*)fHistDir->Get("hmupt1"))->Fill(ptMuon);

    if(itrk<0) continue; // skip muons without inner tracker info
    ((TH1D*)fHistDir->Get("hmupt2"))->Fill(ptMuon);

    if(PRINT) cout<<" muon "<<it<<" index "<<itrk<<" "<<muon->fQ<<" "<<muon->fGlobalPlab.Perp()<<" "<<ptMuon<<" "<<muon->fValidHits
		  <<" "<<muon->fMuonChi2<<" m-id "<<hex<<muonId<<dec<<" "<<muonIdx<<" gen "<<genidx<<endl;

    bool muid  = tightMuon(muon);  // does it work on muons?
    if(muid) {
      if(PRINT) cout<<" A tight muon"<<endl;
      ((TH1D*)fHistDir->Get("hmupt3"))->Fill(ptMuon);
    }

    // Now do the rigger matching
    double dR = doTriggerMatchingR(muon,false,true,false); // see if it matches HLT muon
    ((TH1D*)fHistDir->Get("dr3"))->Fill(dR);
    //if(muid) ((TH1D*)fHistDir->Get("dr*"))->Fill(dR);

    if(dR<dRMin) {
      muHlt++;
      if(muid) {  // consider only tight muons
	muTightHlt++;
	if(PRINT) cout<<" Matched to HLT "<<dR<<" "<<it<<" "<<itrk<<" "<<muid<<" "<<muHlt<<" "<<muTightHlt<<endl;
	//cout<<" Matched to HLT "<<dR<<" "<<it<<" "<<itrk<<" "<<muid<<" "<<muHlt<<" "<<muTightHlt<<endl;
	hltMatchedMuons.push_back(itrk);
      }
    }
    // Going through the SimpleTrack gives the same info
    //TSimpleTrack *pTrack = fpEvt->getSimpleTrack(itrk);
    //cout<<it<<" "<<itrk<<" "<<pTrack->getP().Perp()<<endl; // same as direct muon

  } // END MUON LOOP

  if(PRINT) cout<<" HLT matched muons in this event "<<muHlt<<" tight "<<muTightHlt<<endl;
  ((TH1D*)fHistDir->Get("h11"))->Fill(float(muHlt));
  ((TH1D*)fHistDir->Get("h12"))->Fill(float(muTightHlt));

  return muHlt;

}
//----------------------------------------------------------------------
int candAnaDstar::doTest(TAnaCand *pC, int mode) {
  //bool PRINT = false;
  bool PRINT = true;

  bool printSimple = ( (mode==0) || (mode>=10&&mode<=19) );
  bool printSignal = ( (mode==0) || (mode>=20&&mode<=29) );
  bool printMuons  = ( (mode==0) || (mode>=30&&mode<=39) );
  bool printHLT    = ( (mode==0) || (mode>=40&&mode<=49) );
  bool printHLTObj = ( (mode==0) || (mode>=50&&mode<=59) );
  bool printHLTMuObj = ( (mode==0) || (mode>=60&&mode<=69) );
  bool printHLTInfo = ( (mode==0) || (mode>=70&&mode<=79) );

  // disable printing for 11,21,31,41,51
  if( ((mode%10)!=0) && ((mode%10)!=9) ) PRINT=false;
  int status = -1;

  if(printSimple) {

  // Check simple tracks
  if(PRINT) cout<<"List all tracks "<<fpEvt->nRecTracks()<<"  "<<fpEvt->nSimpleTracks()<<endl;
  // loop over simple tracks
  for(int itrk=0; itrk<(fpEvt->nSimpleTracks()); itrk++ ) {  // if the simple track exists
    TSimpleTrack *pTrack = fpEvt->getSimpleTrack(itrk);
    //cout<<itrk<<" "<<pTrack<<endl;
    if(pTrack == 0) continue;

    //pTrack->dump();

    TVector3 mom = pTrack->getP(); // momentum
    int index = pTrack->getIndex(); // same as itrk
    int q     = pTrack->getCharge();
    int qual  = pTrack->getHighPurity();
    int muonId= pTrack->getMuonID(); // muon id
    int pvidx = pTrack->getPvIndex();
    int genidx=pTrack->getGenIndex();
    int inds  = pTrack->getIndices();
    int bits  = pTrack->getBits();
    double pt = mom.Perp();
    status++; // count for return

    if(PRINT) cout<<" track "<<itrk<<" idx "<<index<<" "<<q<<" "<<qual
		  <<" muon "<<muonId<<" "<<pvidx<<" gen "<<genidx<<" "
		  <<hex<<inds<<" "<<bits<<dec<<" "
		  <<pt<<" "<<mom.Eta()<<" "<<mom.Phi()<<endl;

    ((TH1D*)fHistDir->Get("htrackpt1"))->Fill(pt);

    if(muonId>0) {
      ((TH1D*)fHistDir->Get("htrackpt2"))->Fill(pt);
      if(PRINT) cout<<" muon "<<muonId;
      bool muid  = tightMuon(pTrack);  // does this work on simpleTracks?
      if(muid) {
	if(PRINT) cout<<", tight muon";
	((TH1D*)fHistDir->Get("htrackpt3"))->Fill(pt);
      }
      if(PRINT) cout<<endl;
    } // muon

  } // end track loop

  } // printSimple

  if(printSignal) {
  // Check signal tracks
  if(PRINT) cout<<"List signal tracks "<<pC->fSig1<<" "<<pC->fSig2<<endl;
  for(int it = pC->fSig1; it<=pC->fSig2; it++) {  // loop over index of signal tracka
    TAnaTrack *pT = fpEvt->getSigTrack(it); // this gives TAnaTrack
    //pT->dump(); // so all normal TAnaTracks work
    int itrk = pT->fIndex; // index of the SimpleTrack
    double pt=pT->fPlab.Perp();
    int muonId = pT->fMuID; // muon id
    status++;

    if(PRINT)cout <<" signal "<<it<<" "<<itrk <<" "<<pT->fMCID<<" "<<pT->fGenIndex<<" "<<pt<<" "<<pT->fQ<<" "
		  <<pT->fChi2 <<" "<<pT->fDof<<" "<<pT->fTrackQuality<<" "<<pT->fAlgorithm<<" "<<muonId<<" "
		  <<pT->fMuIndex << " "<<pT->fPvIdx<<endl;

    ((TH1D*)fHistDir->Get("hsigpt1"))->Fill(pt);


    // Check tight muon only if muonId>0
    if(muonId>0) {
      if(PRINT) cout<<" muon "<<muonId;
      ((TH1D*)fHistDir->Get("hsigpt2"))->Fill(pt);

      bool muid  = tightMuon(pT);  //
      if(muid) {
	if(PRINT) cout<<", tight muon";
	((TH1D*)fHistDir->Get("hsigpt3"))->Fill(pt);
      }
      if(PRINT) cout<<endl;
    } // if muon
  }  // end signal loop

  } // printSignal


  // MUONS
  if(printMuons) {

  if(PRINT) cout<<"List muon tracks "<< fpEvt->nMuons()<<endl;

  const float dRMin = 0.02;
  int muHlt = 0, muTightHlt = 0;
  for (int it = 0; it< fpEvt->nMuons(); ++it) { // loop over muons
    TAnaMuon *muon = fpEvt->getMuon(it); // get muon

    // some direct muon methods
    //muon->dump();

    // some direct muon methods
    int muonId = muon->fMuID; // muon ID bits
    int muonIdx = muon->fMuIndex;  // index in the muon list = it
    int genidx = muon->fGenIndex; // gen index

    // get track index
    int itrk = muon->fIndex; // index of the SimpleTrack
    if(itrk<0) {
      cout<<" A muon without inner track info, skip "<<it<<" "<<muonId<<" "<<muonIdx<<endl;
      continue; // skip muons without inner tracker info
    }

    TVector3 muonMom = muon->fPlab;
    double ptMuon  = muonMom.Perp();
    if(PRINT)
      cout<<" muon "<<it<<" index "<<itrk<<" "
	  <<ptMuon<<" "<<muonMom.Eta()<<" "<<muonMom.Phi()<<" "
	//<<muon->fValidHits<<" "<<muon->fMuonChi2<<" m-id "
	//<<muon->fQ<<" "<<muon->fGlobalPlab.Perp()<<" "
	  <<hex<<muonId<<dec<<" "<<muonIdx<<" gen "<<genidx<<endl;

    bool muid  = tightMuon(muon);  // does it work on muons?
    if(muid && PRINT) cout<<"  a tight muon"<<endl;

    // Now do the rigger matching
    double dR = doTriggerMatchingR(muon,false,true,false); // see if it matches HLT muon
    if(dR<dRMin) {
      muHlt++; // just count
      if(muid) {  // consider only tight muons
	muTightHlt++;  // just count
	if(PRINT)
	  cout<<"  matched to HLT "<<dR<<" "<<it<<" "<<itrk<<" "<<muid<<" "<<muHlt<<" "<<muTightHlt<<endl;
      }
    }
  } // END MUON LOOP

  if(PRINT) cout<<" HLT matched muons in this event "<<muHlt<<" tight "<<muTightHlt<<endl;

  status = muTightHlt;
  } // printMuons

  // HLT 40-49
  if(printHLT) {
    int codedtrig=0;
    int count=0;
    ((TH1D*)fHistDir->Get("htest1"))->Fill(0.);
    int num=0;
    bool triggerId[12];
    for(int i=0;i<12;++i) {triggerId[i]=false;}
    if(PRINT) cout<<" List HLT "<<NHLT<<endl;
    for (int i = 0; i < NHLT; ++i) {
      TString a = fpEvt->fHLTNames[i];
      int    ps = fpEvt->fHLTPrescale[i];
      bool wasRun = fpEvt->fHLTWasRun[i];
      bool result = fpEvt->fHLTResult[i];
      bool error  = fpEvt->fHLTError[i];

      if (wasRun && result) {
	bool rightDS = true;
	if(DSNAME!= "")  // do only if DSNAME was define (using PDTRIGGER)
	  rightDS = fpReader->pdTrigger()->triggerInPd(DSNAME, a.Data());

	if(PRINT) cout << "Triggered "<<a << " "<<ps<<" "<<error<<" "<<rightDS<<" "<<i;
	if(rightDS) {
	  ++count;
	  ((TH1D*)fHistDir->Get("htest1"))->Fill(19.);
	  int status0=getHLTId(a);
	  if(status0<12) triggerId[status0]=true;
	  if(PRINT) cout << " right DS "<<fhltType<<" "<<status0;
	  num++; // there is space for 5 triggers
	  if(num<6) codedtrig = (codedtrig * 100) + status0; // x1x2y1y2z1z2....
	  //cout<<status0<<" "<<num<<" "<<status<<endl;
	} // DS
	if(PRINT) cout<<endl;
      } // ifRun
    } // for loop
    if(PRINT) cout <<" num trigs "<<num<<" stat "<<count<<"/"<<codedtrig<<endl;

    for(int i=0;i<12;++i) {
      if(triggerId[i]) ((TH1D*)fHistDir->Get("htest1"))->Fill(float(i));
    }

    if(mode==42) status=count;
    else         status=codedtrig;

  } // if

  //
  if(printHLTObj) {  // mode 50-59
    // Look at Trig Object v2
    map<int, int, less<int> > hltObjL1Map,hltObjL2Map,hltObjL3Map,hltObjL3NoMuMap;
    int count_final=0;

    TTrgObjv2 *pTO;
    // if(PRINT) cout<<" Dump TTrgObjv2 "<<fpEvt->nTrgObjv2()<<endl;
    // for (int i = 0; i < fpEvt->nTrgObjv2(); ++i) {
    //   pTO = fpEvt->getTrgObjv2(i);
    //   //pTO->dump();

    //   int hltIndex = pTO->fHltIndex;


  //cout<<" size "<<hltObjMap.size()<<endl;
  map<unsigned int, unsigned int, less<unsigned int> >::iterator  ix;
  for(int i=0; i!=fpEvt->nTrgObjv2(); i++) { // loop over all objects
    pTO = fpEvt->getTrgObjv2(i);
    int hltIndex = pTO->fHltIndex;

    bool lastModule= false;
    bool activeModule = false;
    ix=hltObjMap.find(i);
    //cout<<i<<" "<<hltIndex<<" ";
    //cout<<i<<" selected hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
    //	<<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber;
    if(ix!=hltObjMap.end()) {
      activeModule=true; // signal the teh module belongs a to selected/active HLT
      int num = ix->first;  // index of the module
      int hltN = (ix->second & 0x7FFFFFFF); // index of the HLT path the module belongs to
      lastModule = ( (ix->second & 0x80000000) != 0); // last module flag
      //cout<<" "<<num<<" "<<hltN<<" "<<lastN;
      if(num != i) cout<<" very very wrong1 "<<num<<" "<<i<<endl;
#ifdef OLD_OBJ_MARK
      if(hltN != (hltIndex%1000) ) cout<<" very very wrong2 "<<hltN<<" "<<hltIndex<<endl;
      if((hltIndex<1000) ) cout<<" very very wrong3 "<<hltIndex<<endl;
      if(lastModule != (hltIndex>1000000) ) cout<<" very very wrong4 "<<lastModule<<" "<<hltIndex<<endl;
#else
      if(hltN != hltIndex ) cout<<" very very wrong2 "<<hltN<<" "<<hltIndex<<endl;
#endif
    }
    //} else cout<<" nothing in map ";
    //cout<<endl;

#ifdef OLD_OBJ_MARK
    if( (hltIndex>1000) != activeModule)
      cout<<" very very wrong5 "<<activeModule<<" "<<hltIndex<<endl;
    if( (hltIndex>100000) != lastModule)
      cout<<" very very wrong6 "<<lastModule<<" "<<hltIndex<<endl;
#endif // OLD_OBJ_MARK

      // list only the selected guys
      //continue; // this object was selected, matches our trigger list
      if(lastModule) count_final++;

      if(PRINT) {
	if( lastModule || ( (mode==59) && activeModule ) )  {
	  cout<<i<<" hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
	      <<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber;
	  if(!activeModule) cout<<endl;
	  else cout<<" : selected filter "<<endl;
	}
      } // PRINT

      int level=0;
      if( (pTO->fHltPath).Contains("L2")||(pTO->fLabel).Contains("L2")||(pTO->fType).Contains("L2") ) {
	level=2;
      } else if( (pTO->fHltPath).Contains("L1")||(pTO->fLabel).Contains("L1")||(pTO->fType).Contains("L1") ) {
	level=1;
      } else {
	level=3;
      }

      vector<int> muonIndex = pTO->fIndex;
      vector<int> muonID = pTO->fID;
      vector<TLorentzVector> muonP = pTO->fP;
      int num = muonIndex.size();
      //if(num<2 || num>99) cout<<" number of trigger muons is "<<num<<" - "
      //		     <<i<<" hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
      //		     <<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber<<endl;

      for(int n=0;n<num;++n) {
	int index = muonIndex[n];
	int id = muonID[n];
	TLorentzVector p = muonP[n];
	// use index as unique is
	//if(lastModule)  {  // select only the last object in the path
	if(activeModule)  {  // select all object in the path
	  if(level==1) {
	    hltObjL1Map[index]++;  // save the index
	  } else if(level==2) {
	    hltObjL2Map[index]++;  // save the index
	  } else if(level==3) {
	    if(abs(id)==13) // muons
	      hltObjL3Map[index]++;  // save the index
	    else // track
	      hltObjL3NoMuMap[index]++;  // save the index
	  } else {
	    cout<<" elevel not set "<<endl;
	  }
	}

	if(PRINT)
	  if( lastModule || ( (mode==59) && activeModule ) )
	    cout<<" "<<n<<" index "<<index<<" id "<<id<<" pt/eta/phi "<<p.Pt()<<" "<<p.Eta()<<" "<<p.Phi()<<endl;

      }

    }

    if(PRINT) {
    cout<<hltObjL1Map.size()<<" "<<hltObjL2Map.size()<<" "<<hltObjL3Map.size()<<" "<<hltObjL3NoMuMap.size()<<endl;
    for (map<int, int, less<int> >::iterator ix = hltObjL1Map.begin(); ix != hltObjL1Map.end(); ++ix)
      cout<<" L1 "<<ix->first<<" "<<ix->second<<endl;
    for (map<int, int, less<int> >::iterator ix = hltObjL2Map.begin(); ix != hltObjL2Map.end(); ++ix)
      cout<<" L2 "<<ix->first<<" "<<ix->second<<endl;
    for (map<int, int, less<int> >::iterator ix = hltObjL3Map.begin(); ix != hltObjL3Map.end(); ++ix)
      cout<<" L3 "<<ix->first<<" "<<ix->second<<endl;
    for (map<int, int, less<int> >::iterator ix = hltObjL3NoMuMap.begin(); ix != hltObjL3NoMuMap.end(); ++ix)
      cout<<" L3-NoMu "<<ix->first<<" "<<ix->second<<endl;
    }

    //for(map<unsigned int, unsigned int, less<unsigned int>>::iterator
    //	  iter=hltObjMap.begin(); iter!=hltObjMap.end(); iter++) {
    //cout<<hex<<*iter<<" "<<dec<<endl;
    //}


    if(mode==51)      status = hltObjL1Map.size();
    else if(mode==52) status = hltObjL2Map.size();
    else if(mode==53) status = hltObjL3Map.size();
    else if(mode==54) status = hltObjL3NoMuMap.size();
    else              status = count_final;

  } // end printHLTMuObj

  // Print old HLT object
  // if(printHLTMuObj) {  // mode 60-69
  //   int countl1=0, countl2=0, countl3=0;
  //   TTrgObj *p;
  //   if(PRINT) cout<<" Dump TTrgObj "<<fpEvt->nTrgObj()<<endl;
  //   for (int i = 0; i < fpEvt->nTrgObj(); ++i) {
  //     p = fpEvt->getTrgObj(i);
  //     if(PRINT) p->dump();
  //     //cout<<i<<"  "<< p->fLabel << " number " << p->fNumber <<" ID = "
  //     //    << p->fID << " pT = " << p->fP.Perp()
  //     //    << " eta = " << p->fP.Eta()<< " phi = " << p->fP.Phi() << " "
  //     //    <<p->fID << endl;

  //     if( (p->fLabel).Contains("L3") )      {countl3++;}
  //     else if( (p->fLabel).Contains("L2") ) {countl2++;}
  //     else if( (p->fLabel).Contains("L1") ) {countl1++;}
  //   } // end for

  //   if(mode==61) status=countl1;
  //   else if(mode==62) status=countl2;
  //   else  status=countl3;
  //   // Check L1, not really used
  //   if(0) {
  //     cout << "--------------------  L1" << endl;
  //     for (int i = 0; i < NL1T; ++i) {
  // 	//result = wasRun = error = false;
  // 	auto a = fpEvt->fL1TNames[i];
  // 	auto ps = fpEvt->fL1TPrescale[i];
  // 	auto result = fpEvt->fL1TResult[i];
  // 	auto mask  = fpEvt->fL1TMask[i];
  // 	//if (a.Contains("Mu"))
  // 	if (result) {
  // 	  cout << a <<  " mask: " << mask << " result: " << result << " ps: " << ps << endl;
  // 	}
  //     }
  //   }
  // } // printHLTMuObj

  // print the hlt-info vector 70-79
  if(printHLTInfo) {
    if(PRINT) cout<<" hlt-info "<<hltInfo.size()<<endl;
    // loop over muons
    int c1=0, c2=0, c3=0;
    for(vector<int>::iterator iter=hltInfo.begin(); iter!=hltInfo.end(); iter++) {
      c1++;
      bool accept = ( ((*iter)&0x7030) == 0x0030 );
      bool accept2 = ( ((*iter)&0x7000) == 0x0000 );
      if(accept) c2++;
      if(accept2) c3++;
      if(PRINT) cout<<hex<<*iter<<" "<<((*iter)&0x7030)<<dec<<" "<<accept<<" "<<c1<<" "<<c2<<endl;
    }
    if(PRINT) cout<<" all/accepted "<<c1<<"/"<<c2<<endl;
    status = c2;
  }

  return status;
}

//----------------------------------------------------------------------
bool candAnaDstar::analyzeHltInfo(bool allTrigClean) {
  const bool PRINT = false;
 bool status=false;
 if(PRINT) cout<<" hlt-info "<<hltInfo.size()<<endl;
 // loop over muons
 int c1=0, c2=0, c3=0;
 for(vector<int>::iterator iter=hltInfo.begin(); iter!=hltInfo.end(); iter++) {
   c1++;
   bool accept = ( ((*iter)&0x7030) == 0x0030 ); // match track and no overlap with pi/K
   bool accept2 = ( ((*iter)&0x7000) == 0x0000 ); // no overlap with pi/K
   if(accept) c2++;
   if(accept2) c3++;
   if(PRINT) cout<<hex<<*iter<<" "<<((*iter)&0x7030)<<dec<<" "<<accept<<" "<<c1<<" "<<c2<<endl;
 }
 if(PRINT) cout<<" all/accepted "<<c1<<"/"<<c2<<endl;

 if(allTrigClean) {status = (c1==c2);} // all triggers clean
 else             {status = (c2>0);} // at least one clean trigger

 return status; // it is !veto, so 1-ok, 0-veto the event
}
// ----------------------------------------------------------------------
int candAnaDstar::getHLTId(TString a) {
  int status=0;

  if(fYear == 2015 ) {
    // For Charmonium 2015
    if (a.Contains("_Bs_")) {status=1; fb8=true;}
    else if(a.Contains("QuadMuon0_Dimuon0_Jpsi")) {status=2; }
    else if(a.Contains("_Jpsi_Muon")) {status=3; }
    // 4 empty
    else if(a.Contains("_Jpsi_Displaced"))    {status=5; fb9=true;}
    else if(a.Contains("_JpsiTrk_Displaced")) {status=6; }
    else if(a.Contains("_Jpsi_NoVertexing"))  {status=7; }
    else if(a.Contains("_Jpsi_NoOS_"))   {status=7; }
    else if(a.Contains("_Mu7p5_L2Mu2_")) {status=11; }
    else if(a.Contains("_Mu7p5_Track"))  {status=10; }
    else if(a.Contains("_Jpsi_"))     {status=8; }
    else if(a.Contains("_PsiPrime"))  {status=9; }

  } else if(fYear==2012) {
    // For MuOnia 2012
    if (a.Contains("_Bs_")) {status=1; fb8=true;}
    else if(a.Contains("Dimuon3p5_SameSign")) {status=2; }
    else if(a.Contains("Dimuon0_Jpsi_Muon"))  {status=3; }
    else if(a.Contains("Dimuon0_Upsilon_Muon")) {status=4; }
    else if(a.Contains("_Upsilon_"))          {status=4; }
    else if(a.Contains("_Jpsi_Displaced"))    {status=5; fb9=true;}
    else if(a.Contains("_JpsiTk_Displaced"))  {status=6; }
    else if(a.Contains("_Jpsi_"))     {status=7; }
    else if(a.Contains("_Mu5_Track"))  {status=8; }
    else if(a.Contains("_Mu7_Track"))  {status=8; }
    else if(a.Contains("_PsiPrime_")) {status=9; }
    else if(a.Contains("_Tau2Mu_ItTrack_")) {status=10; }
    else if(a.Contains("_Mu5_L2Mu3_")) {status=11; }
  }


  return status;
}
// ----------------------------------------------------------------------
void candAnaDstar::dumpHFDstarCand(TAnaCand *pC) {
  TAnaTrack *pT(0);


  // -- D0 daughters
  if (pC->fDau1 < 0) {
    cout << "XXXXXXXXX cannot get daughter cand of " << pC->fType << endl;
    return;
  }

  TAnaCand *pD = fpEvt->getCand(pC->fDau1);
  cout << "HFDstarCand: idx = " << pC->fIndex << " type = " << pC->fType
       << " m* = " << pC->fMass << " m0 = " << pD->fMass << " dm = " << pC->fMass-pD->fMass << endl;

  // Check D0 tracks
  cout<<" loop over D0 tracks "<<endl;
  for (int id = pD->fSig1; id <= pD->fSig2; ++id) {

    pT = fpEvt->getSigTrack(id); // this gives TAnaTrack
    //pT->dump(); // so all normal TAnaTracks work
    //cout <<id<<" "<<pT->fIndex <<" "<<pT->fMCID<<" "<<pT->fGenIndex<<" "<<pT->fPlab.Perp()<<" "<<pT->fQ<<" "
    // <<pT->fChi2 <<" "<<pT->fDof<<" "<<pT->fTrackQuality<<" "<<pT->fAlgorithm<<" "<<pT->fMuID<<" "
    // <<pT->fMuIndex << " "<<pT->fPvIdx<<endl;

    int index = pT->fIndex; // index of the SimpleTrack
    int genidx = pT->fGenIndex;
    if (211 == pT->fMCID) {  // pion
      cout<<" pion index "<<index<<" id "<< pT->fMCID<<" pt/eta/phi "
	  <<pT->fPlab.Perp()<<" "<<pT->fPlab.Eta()<<" "<<pT->fPlab.Phi()<<" "<<" gen " << genidx<<endl;
    } else {  // kaon
      cout<<" kaon index "<<index<<" id "<< pT->fMCID<<" pt/eta/phi "
	  <<pT->fPlab.Perp()<<" "<<pT->fPlab.Eta()<<" "<<pT->fPlab.Phi()<<" "<<" gen " << genidx<<endl;
    }
  }


  // -- slow pion
  //cout<<" simpler way "<<endl;
  int it = pC->fSig1;  // index of the signal track
  pT = fpEvt->getSigTrack(it); // this gives TAnaTrack
  cout << " slow pion index "<<pT->fIndex<<" id "<<pT->fMCID << " pt/eta/phi "
       <<pT->fPlab.Perp()<<" "<<pT->fPlab.Eta()<<" "<<pT->fPlab.Phi()<<" gen " << pT->fGenIndex<<endl;
  //pT->dump();

}

// ----------------------------------------------------------------------
// -- works ONLY for this dedicated decay mode with the seq vtx fitter!
int candAnaDstar::truthMatch(TAnaCand *pCand, int verbose) {
  const int verboseOffset=9;
  if(verbose>(1+verboseOffset)) cout<<" truthMatch "<<verbose<<endl;
  
  // -- check slow pion, the signal track for Dstar should be the slow pion 
  int index = pCand->fSig1;
  TAnaTrack *pT = fpEvt->getSigTrack(index);
  //pT->dump();
  int genIndex = pT->fGenIndex;
  
  // if not gen index exit
  if (genIndex < 0) {if(verbose > verboseOffset) cout << "pT->fGenIndex < 0" << endl;return 0;}
  
  if(verbose>(1+verboseOffset)) cout<<" pi slow "<<index<<" "<<genIndex<<endl;
  
  // Pi slow gen info
  //TGenCand *candGenSlowPi = fpEvt->getGenTWithIndex(genIndex); // gen cand slow
  TGenCand *candGenSlowPi = fpEvt->getGenCand(genIndex); // gen cand 
  
  // exit if empty 
  if (0 == candGenSlowPi) 
    {if(verbose > verboseOffset)cout << "0 == candGenSlowPi" << endl;return 0;}
  // check PID, that it is a pion in MC
  if(211 != TMath::Abs(candGenSlowPi->fID)) 
    {if(verbose>verboseOffset) cout << "211 != TMath::Abs(pG->fID)" << endl;return 0;}

  int moSlowPion = candGenSlowPi->fMom1; // save index of the mother of slow
  if(verbose > 1+verboseOffset) 
    cout << "slow pion " << pT->fIndex << " OK, fGenIndex = " << pT->fGenIndex << " "
	 <<moSlowPion<<" "<<candGenSlowPi->fID<<" "<<endl;
  
  // Look for mother (Dstar)
  //TGenCand *pDS = fpEvt->getGenTWithIndex(moSlowPion);  // Dstar
  TGenCand *pDS = fpEvt->getGenCand(moSlowPion);  // Dstar 
  if ((0 == pDS) || 413 != TMath::Abs(pDS->fID)) { // not a Dstar
    if (verbose > verboseOffset) 
      cout << "(0 == pDS) || 413 != pG->fID, pG->fID"<<" moSlowPion = " << moSlowPion << endl;
    return 0;
  }
  
  // -- Now look at the data Dstar candidate, it should have a daughter D0
  if (pCand->fDau1 < 0) 
    {if (verbose >verboseOffset) cout << "no pCand->fDau1" << endl;return 0;}
  
  // Look at DATA D0 tracks
  TAnaCand *pC = fpEvt->getCand(pCand->fDau1); // D0
  // Check D0 tracks
  if(verbose>1+verboseOffset) cout<<" loop over D0 tracks "<<endl;
  // -- check D0 daughters
  int moIdx(-1);
  int count=0, countPi=0, countK=0, missK=0, missPi=0;

  // Loop over the 2 signal tracks from D0
  for (int id = pC->fSig1; id <= pC->fSig2; ++id) { // should be pi+K
    
    pT = fpEvt->getSigTrack(id); // this gives TAnaTrack
    //pT->dump(); // so all normal TAnaTracks work
    //cout <<id<<" "<<pT->fIndex <<" "<<pT->fMCID<<" "<<pT->fGenIndex<<" "<<pT->fPlab.Perp()<<" "<<pT->fQ<<" "
    // <<pT->fChi2 <<" "<<pT->fDof<<" "<<pT->fTrackQuality<<" "<<pT->fAlgorithm<<" "<<pT->fMuID<<" "
    // <<pT->fMuIndex << " "<<pT->fPvIdx<<endl;
    
    int index = pT->fIndex; // index of the SimpleTrack
    int genidx = pT->fGenIndex;
    int type = pT->fMCID;
    if (genidx < 0) 
      {if(verbose>verboseOffset) cout << "no pT->fGenIndex" << endl;return 0;}
    
    //TGenCand *pG = fpEvt->getGenTWithIndex(genidx);  // get GenCand of pi, K
    TGenCand *pG = fpEvt->getGenCand(genidx);  // get GenCand of pi, K
    if(pG<=0) {if (verbose > verboseOffset) cout<<" pG invalid "<<endl;  continue;} //
    
    if (moIdx < 0) { // find mother
      moIdx = pG->fMom1;
    } else { // should be the same
      if (moIdx != pG->fMom1) 
	{if(verbose > verboseOffset)
	    cout<<"pi & K mothers not the same: moIdx != pG->fMom1"<<endl;return 0;}
    }

    if (verbose > 1+verboseOffset)
      cout << "dau cand sigtrack " << id
 	   << " with type = " << type
 	   << " and gen ID = " << pG->fID
 	   << " at gen idx = " << genidx
 	   << endl;
    
    // the MC ID and the assigned ID from reco do not agree
    if((verbose>(1+verboseOffset)) && (TMath::Abs(type)!=TMath::Abs(pG->fID)) )
      cout<<"should be the same:TMath::Abs(type) != TMath::Abs(pG->fID), type=" 
	  << type << " pG->fID = " << pG->fID<< " track " << pT->fIndex << endl;
    
    count++;
    // Look at the GEN PID
    if( TMath::Abs(pG->fID) == 321 ) { // Kaon
      countK++;
      if(verbose > 1+verboseOffset) 
	{cout<<" kaon index "<<index<<" id "<< type <<" pt "<<pT->fPlab.Perp()<<" gen " << genidx<<endl;}
      if( TMath::Abs(type) != 321 ) {  // no a kaon
	missK++;
	if(verbose>verboseOffset) 
	  {cout << " Kaon identified as pion, type = " << type << " pG->fID = " 
		<< pG->fID<< " track " << pT->fIndex<<endl;}
      }
      
    } else if( TMath::Abs(pG->fID) == 211 ) { // Pion
      countPi++;
      if(verbose > 1+verboseOffset) 
	{cout<<" pion index "<<index<<" id "<< type <<" pt "<<pT->fPlab.Perp()
	     <<" gen " << genidx<<endl;}
      if( TMath::Abs(type) != 211 ) {
	missPi++;
	//verbose=100;
	if(verbose>verboseOffset) 
	  cout << " Pion identified as kaon, type = " << type << " pG->fID = " 
	       << pG->fID<< " track " << pT->fIndex<<" "<<missPi<<endl;
	
      } // not a pion 
    } // in if id
  } // for D0 loop 

  // We have the D0  Gen
  //TGenCand *candGenD0 = fpEvt->getGenTWithIndex(moIdx);  // get GenCand of D0
  TGenCand *candGenD0 = fpEvt->getGenCand(moIdx);  // get GenCand of D0
  if(verbose > 1+verboseOffset) 
    cout<<" D0 "<<moIdx <<" "<<candGenD0->fMom1<<endl;
  // -- Get gen-level D0
  if (moIdx < 0) {
    if (verbose > verboseOffset) cout << "pG->fMom1 < 0" << endl;
    return 0;
  }

  //Check that D0 has exacty 2 daughters
  if (candGenD0->fDau2 - candGenD0->fDau1 > 1) { // D0 has more daughters
    if (verbose > verboseOffset) cout << "Do fDau2 - fDau1 > 1" << endl;
    return 0;
  }

  // Get the DStar Gen
  int indexDS = candGenD0->fMom1;
  //TGenCand *candGenDS = fpEvt->getGenTWithIndex(indexDS);  // get GenCand of DS
  TGenCand *candGenDS = fpEvt->getGenCand(indexDS);  // get GenCand of DS
  // Compare the mother of D0 and the slow pion
  if (verbose > 1+verboseOffset) 
    cout<<"DS "<<indexDS<<" "<<candGenDS->fNumber<<" "<<moSlowPion<<endl;
  if (indexDS != moSlowPion) {
    if (verbose >verboseOffset) cout << "pG->fMom1 != moSlowPion" << endl;
    return 0;
  }

  //if (missPi>0 || missK>0)   
  if (verbose > 1+verboseOffset)   
    cout << "===> truth matching OK"<<" Found kaons:"<< countK<<" Found pions "
	 <<countPi<<" Found "<<count<<", Miss K/Pi "<<missK<<"/"<<missPi<<endl;

  if     (missK==0 && missPi==0) {return  1;} // select righ combination
  else if(missK==1 && missPi==1) 
    {cout<<" mixed pi/K, return -1"<<endl; return -1;} // select wrong combination
  
  return 0;
  }
  

// ----------------------------------------------------------------------
void candAnaDstar::moreBasicCuts() {
  cout << "   candAnaDstar: more basic cuts" << endl;
}


// ----------------------------------------------------------------------
void candAnaDstar::bookHist() {

  cout << "==>candAnaDstar: bookHist" << endl;
  candAna::bookHist();

  //fHistDir->cd();

  TH1 *h = new TH1D("Status", "Status", 50, -0.5, 49.5);

  // Dstar histos
  h = new TH1D("mds", "m(dstar)",100, 1.5, 2.5);
  h = new TH1D("mdz", "m(d0)",   100, 1.5, 2.5);
  h = new TH1D("dm", "delta(m)", 80, 0.13, 0.17);
  h = new TH1D("ncand", "ncand", 200, 0., 200);
  h = new TH1D("fls3d", "fls3d", 200, 0., 20);
  h = new TH1D("flsxy", "flsxy", 200, 0., 20);
  h = new TH1D("prob", "prob(chi2/dof)", 100, 0., 1.);
  h = new TH1D("chi2", "chi2", 100, 0., 10.);
  h = new TH1D("alpha", "alpha", 50, 0., 1.0);
  h = new TH1D("dr","dr",100, 0., 1);  // pislow versus D*
  h = new TH1D("dR","dR",100, 0., 1);  // pislow versus D0
  h = new TH1D("pt", "pT", 50, 0., 25);
  h = new TH1D("ptdz", "pT", 50, 0., 25);
  h = new TH1D("ptK",   "pT", 50, 0., 10);
  h = new TH1D("ptPi",  "pT", 50, 0., 10);
  h = new TH1D("ptPis", "pT", 50, 0., 5);
  TH2F *h2 = new TH2F("h2d", "m(d0) vs dm", 60, 1.8, 1.92, 60, 0.13, 0.16);

  h = new TH1D("dm1", "delta(m)", 60, 0.13, 0.16);
  h = new TH1D("dm2", "delta(m)", 60, 0.13, 0.16);

  // overview histos
  h = new TH1D("htrackpt1","all track pT", 100, 0., 100);
  h = new TH1D("htrackpt2","muon track pT", 100, 0., 100);
  h = new TH1D("htrackpt3","tight muon track pT", 100, 0., 100);
  h = new TH1D("hsigpt1","all signal pT", 100, 0., 100);
  h = new TH1D("hsigpt2","muon signal pT", 100, 0., 100);
  h = new TH1D("hsigpt3","tight muon signal pT", 100, 0., 100);
  h = new TH1D("hmupt1","all muon pT", 100, 0., 100);
  h = new TH1D("hmupt2","tracker muon pT", 100, 0., 100);
  h = new TH1D("hmupt3","tight muon pT", 100, 0., 100);

  // test histos
  h = new TH1D("dr1", "dr1", 400, 0., 2.);
  h = new TH1D("dr2", "dr2", 400, 0., 2.);
  h = new TH1D("dr3", "dr3", 400, 0., 2.);
  h = new TH1D("dr4", "dr4", 400, 0., 2.);
  h = new TH1D("dr5", "dr5", 400, 0., 2.);
  h = new TH1D("dr6", "dr6", 400, 0., 2.);
  h = new TH1D("dr7", "dr7", 400, 0., 2.);
  h = new TH1D("dr8", "dr8", 400, 0., 1.);
  h = new TH1D("dr9", "dr9", 400, 0., 2.);
  h = new TH1D("dr10", "dr10", 400, 0., 1.);
  h = new TH1D("dr11", "dr11", 400, 0., 1.);
  h = new TH1D("dr12", "dr12", 400, 0., 1.);
  h = new TH1D("dr13", "dr13", 400, 0., 1.);
  h = new TH1D("dr14", "dr14", 400, 0., 1.);
  h = new TH1D("dr15", "dr15", 400, 0., 1.);
  h = new TH1D("dr16", "dr16", 400, 0., 1.);
  h = new TH1D("dr17", "dr17", 400, 0., 1.);
  h = new TH1D("dr18", "dr18", 400, 0., 1.);
  h = new TH1D("dr19", "dr19", 400, 0., 1.);
  h = new TH1D("dr20", "dr20", 400, 0., 1.);
  h = new TH1D("dr21", "dr21", 400, 0., 1.);
  h = new TH1D("dr30", "dr30", 400, 0., 1.);
  h = new TH1D("dr31", "dr31", 400, 0., 1.);
  h = new TH1D("dr32", "dr32", 400, 0., 1.);
  h = new TH1D("dr33", "dr33", 400, 0., 1.);
  h = new TH1D("dr40", "dr40", 400, 0., 1.);
  h = new TH1D("dr41", "dr41", 400, 0., 1.);
  h = new TH1D("dr42", "dr42", 400, 0., 1.);
  h = new TH1D("dr43", "dr43", 400, 0., 1.);
  h = new TH1D("dr44", "dr44", 400, 0., 1.);

  h = new TH1D("pvidx", "PVidx", 50, 0., 50.);
  h = new TH1D("h1", "h1", 100, 0., 1);
  h = new TH1D("h2", "h2", 100, 0., 1);
  h = new TH1D("h3", "h3", 100, 0., 1);
  h = new TH1D("h4", "h4", 350, 0., 3.5);
  h = new TH1D("h5", "h5", 100, 0., 1);
  h = new TH1D("h6", "h6", 100, 0., 1);
  h = new TH1D("h7", "h7", 350, 0., 3.5);
  h = new TH1D("h8", "h8", 350, 0., 3.5);
  h = new TH1D("h9", "h9", 350, 0., 3.5);
  h = new TH1D("h10","h10",100, 0., 1);
  h = new TH1D("h11", "h11", 10, 0., 10);
  h = new TH1D("h12", "h12", 10, 0., 10);

  h = new TH1D("htest0", "htest0", 20, 0., 20);
  h = new TH1D("htest1", "htest1", 20, 0., 20);
  h = new TH1D("htest2", "htest2", 20, 0., 20);
  h = new TH1D("htest3", "htest3", 20, 0., 20);
  h = new TH1D("htest4", "htest4", 20, 0., 20);
  h = new TH1D("htest5", "htest5", 20, 0., 20);
  h = new TH1D("htest6", "htest6", 20, 0., 20);
  h = new TH1D("htest7", "htest7", 20, 0., 20);
  h = new TH1D("htest8", "htest8", 20, 0., 20);
  h = new TH1D("htest9", "htest9", 20, 0., 20);
  h = new TH1D("htest20", "htest20", 20, 0., 20);
  h = new TH1D("htest21", "htest21", 20, 0., 2.);

  h2 = new TH2F("htest10", "htest10", 10,0.,10.,6,0.,6.);
  h2 = new TH2F("htest11", "htest11", 10,0.,10.,6,0.,6.);
  h2 = new TH2F("htest12", "htest12", 10,0.,10.,6,0.,6.);
  h2 = new TH2F("htest13", "htest13", 10,0.,10.,6,0.,6.);
  h2 = new TH2F("htest14", "htest14", 10,0.,10.,6,0.,6.);
  h2 = new TH2F("htest15", "htest15", 10,0.,10.,6,0.,6.);
  h2 = new TH2F("htest16", "htest16", 10,0.,10.,10,0.,10.);
  //h2 = new TH2F("htest17", "htest17", 10,0.,10.,10,0.,10.);
  //h2 = new TH2F("htest18", "htest18", 10,0.,10.,10,0.,10.);
  //h2 = new TH2F("htest19", "htest19", 10,0.,10.,10,0.,10.);

  // MC histos
#ifdef MC_HISTOS
  cout<<"  HERE HERE "<<endl;
  h = new TH1D("h300", "h300", 10, -5., 5.);
  h = new TH1D("h342", "h342", 100, 1.5, 2.5);
  h = new TH1D("h344", "h344", 100, 1.5, 2.5);
  h = new TH1D("h345", "h345", 100, 1.5, 2.5);
  h = new TH1D("h346", "h346", 100, 1.5, 2.5);
#endif

  // STD: Add variables to the standard redtrees
  fTree->Branch("ftm",&ftm,"ftm/I");
  fTree->Branch("fmuid1",&fmuid1,"fmuid1/O");
  fTree->Branch("fmuid2",&fmuid2,"fmuid2/O");
  fTree->Branch("fmumat1",&fmumat1,"fmumat1/O");
  fTree->Branch("fmumat2",&fmumat2,"fmumat2/O");
  fTree->Branch("fmds",&fmds,"fmds/F");
  fTree->Branch("fmdz",&fmdz,"fmdz/F");

  fTree->Branch("ffls3d",&ffls3d,"ffls3d/F");
  fTree->Branch("fchi2",&fchi2,"fchi2/F");
  fTree->Branch("falpha",&falpha,"falpha/F");
  fTree->Branch("fqpis",&fqpis,"fqpis/F");
  fTree->Branch("fdr",&fdr,"fdr/F");

  fTree->Branch("fpt",&fpt,"fpt/F");
  fTree->Branch("fptdz",&fptdz,"fptdz/F");
  fTree->Branch("fptpis",&fptpis,"fptpis/F");
  fTree->Branch("fptpi",&fptpi,"fptpi/F");
  fTree->Branch("fptk",&fptk,"fptk/F");

  fTree->Branch("feta",&feta,"feta/F");
  fTree->Branch("fetapi",&fetapi,"fetapi/F");
  fTree->Branch("fetak",&fetak,"fetak/F");

  fTree->Branch("fchipi",&fchipi,"fchipi/F");
  fTree->Branch("fchik",&fchik,"fchik/F");
  fTree->Branch("fiso",&fiso,"fiso/F");
  fTree->Branch("fnclose",&fnclose,"fnclose/I");

  fTree->Branch("fmudr1",&fmudr1,"fmudr1/F");
  fTree->Branch("fmudr2",&fmudr2,"fmudr2/F");

  fTree->Branch("fmatch1dr",    &fmatch1dr,      "fmatch1dr/F");
  fTree->Branch("fmatch2dr",    &fmatch2dr,      "fmatch2dr/F");

  // MVA
  fTree->Branch("fmuidmva1",    &fmuidmva1,      "fmuidmva1/O");
  fTree->Branch("fmuidmva2",    &fmuidmva2,      "fmuidmva2/O");

  fTree->Branch("fmva1",    &fmva1,      "fmva1/D");
  fTree->Branch("fmva2",    &fmva2,      "fmva2/D");

  fTree->Branch("fveto",&fveto,   "fveto/O");

  // tests
  fTree->Branch("fmatch1dr1",    &fmatch1dr1,      "fmatch1dr1/F");
  fTree->Branch("fmatch2dr1",    &fmatch2dr1,      "fmatch2dr1/F");
  fTree->Branch("fmatch1dr2",    &fmatch1dr2,      "fmatch1dr2/F");
  fTree->Branch("fmatch2dr2",    &fmatch2dr2,      "fmatch2dr2/F");
  fTree->Branch("fmatch1dr3",    &fmatch1dr3,      "fmatch1dr3/F");
  fTree->Branch("fmatch2dr3",    &fmatch2dr3,      "fmatch2dr3/F");
  fTree->Branch("fmatch1dr4",    &fmatch1dr4,      "fmatch1dr4/F");
  fTree->Branch("fmatch2dr4",    &fmatch2dr4,      "fmatch2dr4/F");
  fTree->Branch("fmatch1dr5",    &fmatch1dr5,      "fmatch1dr5/F");
  fTree->Branch("fmatch2dr5",    &fmatch2dr5,      "fmatch2dr5/F");

  fTree->Branch("fmatchTrigs",&fmatchTrigs,   "fmatchTrigs/O");
  fTree->Branch("fb2",&fb2,   "fb2/O");
  fTree->Branch("fb3",&fb3,   "fb3/O");
  fTree->Branch("fb4",&fb4,   "fb4/O");
  fTree->Branch("fb5",&fb5,   "fb5/O");
  fTree->Branch("fb6",&fb6,   "fb6/O");
  fTree->Branch("fb7",&fb7,   "fb7/O");
  fTree->Branch("fb8",&fb8,   "fb8/O");
  fTree->Branch("fb9",&fb9,   "fb9/O");

  fTree->Branch("fitmp1",&fitmp1,   "fitmp1/I");
  fTree->Branch("fitmp2",&fitmp2,   "fitmp2/I");
  fTree->Branch("fitmp3",&fitmp3,   "fitmp3/I");
  fTree->Branch("fitmp4",&fitmp4,   "fitmp4/I");
  fTree->Branch("fitmp5",&fitmp5,   "fitmp5/I");

  fTree->Branch("ftmp1",&ftmp1,   "ftmp1/O");
  fTree->Branch("ftmp2",&ftmp2,   "ftmp2/O");
  fTree->Branch("ftmp3",&ftmp3,   "ftmp3/O");
  fTree->Branch("ftmp4",&ftmp4,   "ftmp4/O");
  fTree->Branch("ftmp5",&ftmp5,   "ftmp5/O");
  fTree->Branch("ftmp6",&ftmp6,   "ftmp6/O");
  fTree->Branch("ftmp7",&ftmp7,   "ftmp7/O");
  fTree->Branch("ftmp8",&ftmp8,   "ftmp8/O");
  fTree->Branch("ftmp9",&ftmp9,   "ftmp9/O");

  fTree->Branch("hltt",    &fhltType,           "hltt/I");
  fTree->Branch("hlt",     &fGoodHLT,           "hlt/O");

  // MC
#ifdef MC_HISTOS
  cout<<"  HERE HERE HERE "<<endl;

  fTree->Branch("fmcmds",   &fmcmds,     "fmcmds/F");
  fTree->Branch("fmcmdz",   &fmcmdz,     "fmcmdz/F");
  fTree->Branch("fmcpt",    &fmcpt,      "fmcpt/F");
  fTree->Branch("fmcptdz",  &fmcptdz,    "fmcptdz/F");
  fTree->Branch("fmcptpis", &fmcptpis,   "fmcptpis/F");
  fTree->Branch("fmcptpi",  &fmcptpi,    "fmcptpi/F");
  fTree->Branch("fmcptk",   &fmcptk,     "fmcptk/F");
#endif
}

// ----------------------------------------------------------------------
void candAnaDstar::readCuts(string filename, int dump) {
  candAna::readCuts(filename, dump);

  fCutFile = filename;

  if (dump) cout << "==> candAnaDstar: Reading " << fCutFile << " for cut settings" << endl;
  vector<string> cutLines;
  readFile(fCutFile, cutLines);

  char CutName[100];
  float CutValue;
  //int ok(0);

  char  buffer[200];
  fHistDir->cd();
  TH1D *hcuts = (TH1D*)fHistDir->Get("hcuts");
  hcuts->GetXaxis()->SetBinLabel(200, fCutFile.c_str());
  string cstring = "B cand";

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str());

    //ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);


  }

}

//-----------------------------------------------------------------------------------
// Find best dr match between track p and any simple track
// returns index of the best track
// muons=1, match to muons tracks only
double candAnaDstar::findBestMatchingTrack(TLorentzVector p,int muons,int &idx) {

  double deltaRmin1(100);
  //double trigMatchDeltaPt1 = 99.;
  bool Print = false;
  int match = -1;
  idx=-1;

  // loop over simple tracks
  for(int itrk=0; itrk<(fpEvt->nSimpleTracks()); itrk++ ) {  // if the simple track exists
    TSimpleTrack *pTrack = fpEvt->getSimpleTrack(itrk);
    //cout<<itrk<<" "<<pTrack<<endl;
    if(pTrack == 0) continue;

    TVector3 mom = pTrack->getP(); // momentum
    int index = pTrack->getIndex(); // same as itrk
    //int q     = pTrack->getCharge();
    //int qual  = pTrack->getHighPurity();
    int muonId= pTrack->getMuonID(); // muon id
    //int pvidx = pTrack->getPvIndex();
    //int genidx=pTrack->getGenIndex();
    //int inds  = pTrack->getIndices();
    //int bits  = pTrack->getBits();
    //double pt = mom.Perp();

    if( (muons==1) &&  (muonId==0) ) continue; // skip non muons

    TLorentzVector tp;
    tp.SetPtEtaPhiM(mom.Perp(),mom.Eta(),mom.Phi(),MMUON); // assume a muon

    // check direction matching
    double deltaR1 = p.DeltaR(tp);

    if(Print) cout<<" track "<<itrk<<" idx "<<index<<" "<<mom.Perp()<<" "<<mom.Eta()<<" "<<mom.Phi()
		  <<" "<<deltaR1<<endl;

    if (deltaR1<deltaRmin1) {
      deltaRmin1=deltaR1;
      match = itrk;
      idx=index;

      if (Print) {cout << " selected "<< idx <<endl;}

      // check now the pt matching
      //double trigMatchDeltaPt=999.;
      //if (tp.Rho() > 0.) trigMatchDeltaPt = TMath::Abs(p.Rho()  - tp.Rho())/tp.Rho();
      //if( trigMatchDeltaPt < deltaPtMatch ) trigMatchDeltaPt1=trigMatchDeltaPt;

    } // if direction match

  } // end track loop

  if (Print) {cout << " selected "<< idx <<" "<<match<<" "<<deltaRmin1<<endl;}

  return deltaRmin1;
}

// ----------------------------------------------------------------------
void candAnaDstar::triggerHLT() {

  if ( HLTRANGE.begin()->first == "PDTRIGGER") { //

    const bool skipL12 = true;
    fGoodHLT = false;
    fhltType = -1;
    fHLTPath = "";
    hltObjMap.clear(); // reset the trig object map
    
    TString a;
    int ps(0);
    bool result(false), wasRun(false), error(false);
    
    bool pdTrigger = false;
    if (fVerbose > 2) cout << "PDTRIGGER requested... " << endl;
    fpReader->pdTrigger()->setHLTKey(fRun, fpReader->getFile());  // needed for the new way
    pdTrigger=true;
    
    // Check HLT
    // For every passed HLT look for a matching tigger from our list.
    // If if it confirmed by our list than match it with an object in the TrgObjv2 list
    // Mark the TrigObjv2 object my add a large number to the index.
    // Like this it can be recogised in the track match search.
    //  if ( (fVerbose>9) || (fVerbose==-32) ) cout<<" event "<<fEvt<<endl;
    bool isMuonTrigger=false; // just for diagnostics
    int foundNumHltObjects=0;
    int foundNumHlts=0;
    
    
    for (int i = 0; i < NHLT; ++i) {
      result = wasRun = error = false;
      a = fpEvt->fHLTNames[i];
      ps = fpEvt->fHLTPrescale[i];
      wasRun = fpEvt->fHLTWasRun[i];
      result = fpEvt->fHLTResult[i];
      error  = fpEvt->fHLTError[i];
      
      if (wasRun && result) { // passed
	
	if (fVerbose>1  || (-32 == fVerbose) ) cout << "passed: " << a << endl;
	if ((a == "digitisation_step")
	    || (a == "L1simulation_step")
	    || (a == "digi2raw_step")
	    || (a == "HLTriggerFinalPath")
	    || (a == "raw2digi_step")
	    || (a == "reconstruction_step")
	    ) {
	  //cout<<" does this ever happen? " <<a<<endl;
	  continue; // skip, go to the next triggr
	}

	// loop over our list of HLTs and look for matching
	// We assume that an event can be matched to only one trigger from our list,
	// that is they have to be exclusive.
	bool good = false;
	if (pdTrigger) { // PDTRIGGER mode - accept all triggers which fire and are on the DS list
	  
	  // check that this trigger belongs to our DS
	  bool rightDS = fpReader->pdTrigger()->triggerInPd(DSNAME, a.Data());
	  if (fVerbose>9) cout<<" check hlt-path "<<a.Data()<<" DS name "<<DSNAME<<" included? "<<rightDS<<endl;
	  if (rightDS) { // hlt_path in this DS
	    // check if this is a L1/L2 type trigger or L3, assume that L1/L2 is always in the definition
	    bool isL1L2 = a.Contains("L1") || a.Contains("L2");
	    if (skipL12 && isL1L2) {
	      if (fVerbose>1) cout<<" HIT-path os L1/L2 type, skip "<<a.Data()<<endl;
	      continue;
	    } else good=true; // accept, assume it is L3
	    
	  } else  { // hlt_path not in this DS
	    //if (fVerbose>1)
	    //cout<<" HIT-path not in this DS, skip it: "<<a.Data()<<" "<<rightDS<<" DS name: "<<DSNAME<<endl;
	    continue;
	  } // end if
	  
	} // if pdTrigger


	if (good) {  // for matched hlt paths select the trigger object
	  // this trigger matched one in our list
	  foundNumHlts++;
	  fGoodHLT = true;
	  fhltType = i;
	  fHLTPath = a.Data();
	  isMuonTrigger = a.Contains("Mu") || a.Contains("mu") || a.Contains("MU");
	  //if (ps!=1) cout<<"prescale not one "<<a.Data()<<" "<<ps<<endl;
	  
	  if(!pdTrigger) continue; // skip, so only when pdTrigger is selected
	  
	  bool foundHltObject = false;
	  int countModules=0, lastIndex=-1;
	  TTrgObjv2 *pTO;
	  if ( (fVerbose>9) || (fVerbose==-32))
	    cout<<" TTrgObjv2 objects, size= "<<fpEvt->nTrgObjv2()<<endl;
	  for (int i = 0; i < fpEvt->nTrgObjv2(); ++i) {  // loop over all saved hlt objects
	  pTO = fpEvt->getTrgObjv2(i);
	  
	  if (a == pTO->fHltPath) { // found the right one, matched the hlt-name
	    // this trig object matches a passed and selected trigger
	    foundHltObject=true;
	    foundNumHltObjects++;
	    countModules++;
            lastIndex=i;
	    int hltIndex = pTO->fHltIndex; // HLT path index
	    
            // mark all module as selected by adding a large number, so >1000.
            // do it either here (ALL MODULES)  or below (LAST MODULE)
	    hltObjMap[i]= (hltIndex & 0x7FFFFFFF); // make sure the uper bit is free
	    
	    vector<int> muonIndex = pTO->fIndex;
	    vector<int> muonID = pTO->fID;
	    vector<TLorentzVector> muonP = pTO->fP;
	    int num = muonIndex.size();
	    
	    if ( (fVerbose>9) || (fVerbose==-32)) {
	      cout<<" matched: "<<pTO->fHltPath<<" hlt-index: "<<pTO->fHltIndex<<" module label: "
		  <<pTO->fLabel<<" type: "<<pTO->fType<<" num of particles: "<<num<<endl;
	      pTO->dump();
	    }
	    
	  } // if matched
	  } // end for loop

	// only mark the last module for matching, do it by setting the highest bit
	  if (lastIndex>-1) {
	    pTO = fpEvt->getTrgObjv2(lastIndex);
	    hltObjMap[lastIndex]= (hltObjMap[lastIndex] | 0x80000000); // set top bit for final
	  }
	  
	  if (!foundHltObject && isMuonTrigger)
	    cout<<"Warning: canAna::triggerSelection: matching trigger module not found! "
		<<a<<endl;
	} // end if fGoodHLT
      } // if passed
    } // end for loop hlt

    // This is just to mark the trigger type and number of triggers for this event
    // the type only valid for the last type if there are more than 1
    //if (fhltType>999) {cout<<fhltType<<endl; fhltType -= 1000;}
    //fhltType = fhltType + (1000 * (foundNumHlts-1));
    
    // Diagnostics printout
    //  if ( (fVerbose>9) || (fVerbose==-32))
    //    cout<<" number of found matching hlt objects: "<<foundNumHltObjects<<endl;
    
    // Diagnostics printout
    // if (pdTrigger && isMuonTrigger && (foundNumHltObjects==0)) {
    //   cout<<" No matching hlt objects found:  "<<endl;
    //   //cout<< " Passed triggers: ";
    //   for (int i = 0; i < NHLT; ++i) {
    //     result = wasRun = false;
    //     a = fpEvt->fHLTNames[i];
    //     //ps = fpEvt->fHLTPrescale[i];
    //     wasRun = fpEvt->fHLTWasRun[i];
    //     result = fpEvt->fHLTResult[i];
    //     //error  = fpEvt->fHLTError[i];
    //     if (wasRun && result) cout << a << " ";
    //   }
    //   cout<<endl;
    
    //}
  
    // TESTS
    if (fVerbose>999) {  // Just testing
      cout << " ----------------------------------------------------------------------" << endl;
      
      TTrgObjv2 *pTO;
      cout<<" Dump TTrgObjv2 "<<fpEvt->nTrgObjv2()<<endl;
      for (int i = 0; i < fpEvt->nTrgObjv2(); ++i) {
	pTO = fpEvt->getTrgObjv2(i);
	vector<int> muonIndex = pTO->fIndex;
	vector<int> muonID = pTO->fID;
	vector<TLorentzVector> muonP = pTO->fP;
	int num = muonIndex.size();
	cout<<i<<" hlt "<<pTO->fHltPath<<" hlt-index "<<pTO->fHltIndex<<" module label "
	    <<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber<<" "<<num<<endl;
	
	for(int n=0;n<num;++n) {
	  int index = muonIndex[n];
	  int id = muonID[n];
	  TLorentzVector p = muonP[n];
	  cout<<n<<" index "<<index<<" id "<<id<<" pt/eta/phi "<<p.Pt()<<" "<<p.Eta()<<" "<<p.Phi()<<endl;
	}
	
      }
      
      // print the hlt object map
      for(map<unsigned int, unsigned int, less<unsigned int> >::iterator iter=hltObjMap.begin();
	  iter!=hltObjMap.end(); ++iter) {
	cout<<hex<<iter->first<<" "<<iter->second<<" "<<dec<<endl;
      }
      
    } // end testing
    
    if (false == fGoodHLT && fVerbose > 1) cout << "------->  event NOT triggered!" << endl;
  
    //-------------
  
  } else {
    
    
    fGoodHLT1 = false;
    fhltType = -1;
    fHLT1Path = "nada";
    
    TString a;
    string sa;
    int ps(0);
    bool result(false), wasRun(false), error(false);
    int verbose(fVerbose);
    
    // -- NOTRIGGER, just accept the event
    if (HLTRANGE.begin()->first == "NOTRIGGER") {
      if (verbose>2) cout << "NOTRIGGER requested... " << endl;
      fGoodHLT1 = true;
      fL1T      = true;
      fHLT1Path = "NOTRIGGER";
      fHltPrescale = 1;
      return;
    }
    if (fVerbose == -32) {
      for (int i = 0; i < NHLT; ++i) {
	result = wasRun = error = false;
	a = fpEvt->fHLTNames[i];
	ps = fpEvt->fHLTPrescale[i];
	wasRun = fpEvt->fHLTWasRun[i];
	result = fpEvt->fHLTResult[i];
	error  = fpEvt->fHLTError[i];
	
	if (result) { // passed
	  cout << "triggerHLT::result: " << a << " wasrun = " << wasRun << " ps = " << ps << " run = "
	       << fRun << " ls = " << fLS << " json = " << fJSON
	       << endl;
	}
      }
    }
    
    hltPathInfo hpi;
    string spath;
    int rmin, rmax;
    bool good(false);
    for (map<string, pair<int, int> >::iterator imap = HLTRANGE.begin(); imap != HLTRANGE.end(); ++imap) {
      spath = imap->first;
      rmin = imap->second.first;
      rmax = imap->second.second;
      if (fRun < rmin) continue;
      if (fRun > rmax) continue;
      if (fpReader->fHltPathInfo[spath].result) {
	sa = spath;
	ps = fpReader->fHltPathInfo[spath].prescale;
	good = true;
	break;
      }
      string sas = spath.substr(0, sa.rfind("_v")+2);
      if (fpReader->fHltPathInfo[sas].result) {
	good = true;
	break;
      }
    }
    if (good) {
      fHltPrescale = ps;
      fHLT1Path    = sa;
      fGoodHLT1    = true;
      return;
    }
    
  }
  
} // end tiggerHLT

//----------------------------------------------------------------------------
// Match trigger tracks to simple tracks
// Loop over HLT objects, use only the last module in the path
// muonNum - select the muon (track) number in the trigger object
// trigLevel - select trigger level L1/2/3
double candAnaDstar::doTriggerMatchingTest(int &idx, int muonNum, int trigLevel) {
  int indx1=-1;
  //const bool localPrint = true;
  const int verboseThr = 99;
  bool localPrint = (fVerbose > verboseThr);
  int mu1match(-1);
  string hlt1;
  double deltaRmin1(100);
  TTrgObjv2 *pTO;
  const bool anyTrig = false;
  const bool muonsOnly = false;
  int level=3;
  int count=0;
  fb8=false; // reset the Bs
  fb9=false; // and Jpsi flags

  if (localPrint) cout<<" doTriggerMatchingTest " <<trigLevel<<" "<<muonNum<< endl;

  //cout<<" size "<<hltObjMap.size()<<endl;
  map<unsigned int, unsigned int, less<unsigned int> >::iterator  ix;
  for(int i=0; i!=fpEvt->nTrgObjv2(); i++) { // loop over all objects
    pTO = fpEvt->getTrgObjv2(i);
    int hltIndex = pTO->fHltIndex;

    bool lastModule= false;
    bool activeModule = false;
    ix=hltObjMap.find(i);
    //cout<<i<<" "<<hltIndex<<" ";
    //cout<<i<<" selected hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
    //	<<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber;
    if(ix!=hltObjMap.end()) {
      activeModule=true; // signal the teh module belongs a to selected/active HLT
      int num = ix->first;  // index of the module
      int hltN = (ix->second & 0x7FFFFFFF); // index of the HLT path the module belongs to
      lastModule = ( (ix->second & 0x80000000) != 0); // last module flag
      //cout<<" "<<num<<" "<<hltN<<" "<<lastN;
      if(num != i) cout<<" very very wrong1 "<<num<<" "<<i<<endl;
#ifdef OLD_OBJ_MARK
      if(hltN != (hltIndex%1000) ) cout<<" very very wrong2 "<<hltN<<" "<<hltIndex<<endl;
      if((hltIndex<1000) ) cout<<" very very wrong3 "<<hltIndex<<endl;
      if(lastModule != (hltIndex>1000000) ) cout<<" very very wrong4 "<<lastModule<<" "<<hltIndex<<endl;
#else
      if(hltN != hltIndex ) cout<<" very very wrong2 "<<hltN<<" "<<hltIndex<<endl;
#endif

    }
    //} else cout<<" nothing in map ";
    //cout<<endl;
#ifdef OLD_OBJ_MARK
    if( (hltIndex>1000) != activeModule)
      cout<<" very very wrong5 "<<activeModule<<" "<<hltIndex<<endl;
    if( (hltIndex>100000) != lastModule)
      cout<<" very very wrong6 "<<lastModule<<" "<<hltIndex<<endl;
#endif

  // for(int i=0; i!=fpEvt->nTrgObjv2(); i++) { // loop over all objects
  //   pTO = fpEvt->getTrgObjv2(i);  //pTO->dump();
  //   int hltIndex = pTO->fHltIndex;
    if(localPrint)
	cout<<" hlt "<<i<<" "<<pTO->fHltPath<<" hlt-index "<<hltIndex
	    <<" module label "<<pTO->fLabel<<" "<<activeModule<<endl;

    //if(anyTrig || (hltIndex>1000) ) { // this object was selected, or use all
    if(anyTrig || lastModule ) { // this object was selected, use final last one
      if(localPrint)
	cout<<i<<" selected hlt "<<i<<" "<<pTO->fHltPath<<" hlt-index "<<hltIndex
	    <<" module label "<<pTO->fLabel<<" type "<<pTO->fType
	    <<" number "<<pTO->fNumber<<endl;
      TString a = pTO->fHltPath;
      level=3;
      if(a.Contains("L1"))      {if(localPrint) cout<<" L1 "<<a<<endl; level=1;}
      else if(a.Contains("L2")) {if(localPrint) cout<<" L2 "<<a<<endl; level=2;}
      if( (trigLevel>0) && (trigLevel!=level) ) continue; // skip
      count++;
      int hltId = getHLTId(a);

      vector<int> muonIndex = pTO->fIndex;
      vector<int> muonID = pTO->fID;
      vector<TLorentzVector> muonP = pTO->fP;
      int num = muonIndex.size();
      for(int n=0;n<num;++n) {  // loop over particles in this module, usually 2
	if(muonNum != (n+1) ) continue; // select on 1 track from the object (pair)
	int index = muonIndex[n];
	int id = muonID[n];
	TLorentzVector p = muonP[n];
	// skip checking non-muon objects
	if( (abs(id) != 13) && muonsOnly) continue;

	if(localPrint) {
	  cout<<" particle "<<n<<" index "<<index<<" id "<<id
	      <<" pt/eta/phi "<<p.Pt()<<" "<<p.Eta()<<" "<<p.Phi()<<endl;
	}

	// check direction matching
	int idx0=-1;
	double deltaR1 = findBestMatchingTrack(p,0,idx0); // look at all tracks
	double deltaR2 = findBestMatchingTrack(p,1,idx0); // look only at muons

	if(localPrint) cout <<" dr "<<deltaR1 <<" idx "<<idx0<<" "<<n<<" "<<hltId<<endl;

	// select best
	if (deltaR1<deltaRmin1) {
	  deltaRmin1=deltaR1;
	  mu1match = n;
	  hlt1 = pTO->fLabel;
	  indx1=i;
	  idx=idx0;
	  if(localPrint) cout <<" better dr "<<deltaR1 <<" idx "<<idx0<<" "<<n<<endl;
	} // if match

	// histo here
	//cout<<deltaR1<<" "<<deltaRmin1<<" "<<hltId<<endl;
	((TH1D*)fHistDir->Get("dr10"))->Fill(deltaR1);
	((TH1D*)fHistDir->Get("dr30"))->Fill(deltaR2);
	if(hltId==1)  {
	  ((TH1D*)fHistDir->Get("dr11"))->Fill(deltaR1);
	  if(fGoodHLT) { // HLT passed
	    if(fhltType<1000) { // single trigger
	      ((TH1D*)fHistDir->Get("dr31"))->Fill(deltaR1);
	      ((TH1D*)fHistDir->Get("dr32"))->Fill(deltaR2);

	      // if((deltaRmin1>0.1)) {
	      // 	cout << " bad match "
	      // 	     <<idx<<" "<<indx1<<" "<< deltaR1 << " "<<deltaR2<<" "<<mu1match
	      // 	     <<" "<<hlt1<<" "<<fpEvt->nTrgObjv2()<<" "<<fpEvt->nSimpleTracks()
	      // 	     <<" "<<level<<" "<<muonNum<<" "<<count
	      // 	     <<" "<<fGoodHLT<<" "<<fHLTmatch<<" "<<fhltType<<" "<<fHLTPath
	      // 	     <<endl;
	      // 	cout<<" Bs "<<deltaR1<<" "<<fGoodHLT<<" "<<fHLTmatch<<" "<<fhltType<<" "<<fHLTPath <<endl;
	      // 	for(int i=0; i!=fpEvt->nTrgObjv2(); i++) { // loop over all objects
	      // 	  pTO = fpEvt->getTrgObjv2(i);  //pTO->dump();
	      // 	  int hltIndex = pTO->fHltIndex;
	      // 	  if((hltIndex>1000000) ) cout<<pTO->fHltPath<<"/"<<pTO->fLabel<<" ";
	      // 	} // for i
	      // 	cout<<endl;
	      // 	dumpAll(); // PrintAll
	      //} // if deltaR1
	    }
	  }
	}
	else if(hltId==2) ((TH1D*)fHistDir->Get("dr12"))->Fill(deltaR1);
	else if(hltId==3) ((TH1D*)fHistDir->Get("dr13"))->Fill(deltaR1);
	else if(hltId==4) ((TH1D*)fHistDir->Get("dr14"))->Fill(deltaR1);
	else if(hltId==5) ((TH1D*)fHistDir->Get("dr15"))->Fill(deltaR1);
	else if(hltId==6) ((TH1D*)fHistDir->Get("dr16"))->Fill(deltaR1);
	else if(hltId==7) ((TH1D*)fHistDir->Get("dr17"))->Fill(deltaR1);
	else if(hltId==8) ((TH1D*)fHistDir->Get("dr18"))->Fill(deltaR1);
	else if(hltId==9) ((TH1D*)fHistDir->Get("dr19"))->Fill(deltaR1);
	else if(hltId==10) ((TH1D*)fHistDir->Get("dr20"))->Fill(deltaR1);
	else               ((TH1D*)fHistDir->Get("dr21"))->Fill(deltaR1);

      } // end for loop n

    } // end if valid module

  } // loop over all modules

  if (localPrint)
    cout << " best match "
	 <<idx<<" "<<indx1<<" "<< deltaRmin1 << " "<<mu1match<<" "<<hlt1<<" "<<count<<endl;

  if((count>0) && (deltaRmin1>1.0)) {
    cout << " no match "
	 <<idx<<" "<<indx1<<" "<< deltaRmin1 << " "<<mu1match<<" "<<hlt1
	 <<" "<<fpEvt->nTrgObjv2()<<" "<<fpEvt->nSimpleTracks()
	 <<" "<<level<<" "<<muonNum<<" "<<count
	 <<" "<<fGoodHLT<<" "<<fmatchTrigs<<" "<<fhltType<<" "<<fHLTPath
	 <<endl;
    // for(int i=0; i!=fpEvt->nTrgObjv2(); i++) { // loop over all objects
    //   pTO = fpEvt->getTrgObjv2(i);  //pTO->dump();
    //   int hltIndex = pTO->fHltIndex;

  //cout<<" size "<<hltObjMap.size()<<endl;
  map<unsigned int, unsigned int, less<unsigned int> >::iterator  ix;
  for(int i=0; i!=fpEvt->nTrgObjv2(); i++) { // loop over all objects
    pTO = fpEvt->getTrgObjv2(i);
    int hltIndex = pTO->fHltIndex;

    bool lastModule= false;
    //bool activeModule = false;
    ix=hltObjMap.find(i);
    //cout<<i<<" "<<hltIndex<<" ";
    //cout<<i<<" selected hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
    //	<<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber;
    if(ix!=hltObjMap.end()) {
      //activeModule=true; // signal the teh module belongs a to selected/active HLT
      int num = ix->first;  // index of the module
      int hltN = (ix->second & 0x7FFFFFFF); // index of the HLT path the module belongs to
      lastModule = ( (ix->second & 0x80000000) != 0); // last module flag
      //cout<<" "<<num<<" "<<hltN<<" "<<lastN;
      if(num != i) cout<<" very very wrong1 "<<num<<" "<<i<<endl;
#ifdef OLD_OBJ_MARK
      if(hltN != (hltIndex%1000) ) cout<<" very very wrong2 "<<hltN<<" "<<hltIndex<<endl;
      if((hltIndex<1000) ) cout<<" very very wrong3 "<<hltIndex<<endl;
      if(lastModule != (hltIndex>1000000) ) cout<<" very very wrong4 "<<lastModule<<" "<<hltIndex<<endl;
#else
      if(hltN != hltIndex ) cout<<" very very wrong2 "<<hltN<<" "<<hltIndex<<endl;
#endif

    }
    //} else cout<<" nothing in map ";
    //cout<<endl;

#ifdef OLD_OBJ_MARK
    if( (hltIndex>1000) != activeModule)
      cout<<" very very wrong5 "<<activeModule<<" "<<hltIndex<<endl;
    if( (hltIndex>100000) != lastModule)
      cout<<" very very wrong6 "<<lastModule<<" "<<hltIndex<<endl;
#endif

      if(lastModule) cout<<pTO->fHltPath<<"/"<<pTO->fLabel<<" ";
    }
    cout<<endl;
  }

  return deltaRmin1;
}


//----------------------------------------------------------------------------
// Loop over HLT objects, match the 2 HLT muons to offline tracks
// Customised for Dstar events
// Return true if at least 1 trigger matches offline tracks
bool candAnaDstar::doTriggerMatchingForDs(double &dr1, double &dr2) {

  //const bool localPrint = true;
  bool localPrint = (fVerbose > 99);
  bool passed=false;
  double deltaCut = 0.025; // This is the matching cut
  int mu1match=-1;
  int mu2match=-1;
  int indx1=-1;
  int indx2=-1;
  string hlt1;
  string hlt2;
  double deltaRmin1=100;
  double deltaRmin2=100;
  TTrgObjv2 *pTO;
  int level=3;
  int count0=0, count1=0;
  int idum=-1;

  if (localPrint) cout<<" doTriggerMatchingForDs " << endl;

  //cout<<" size "<<hltObjMap.size()<<endl;
  map<unsigned int, unsigned int, less<unsigned int> >::iterator  ix;
  for(int i=0; i!=fpEvt->nTrgObjv2(); i++) { // loop over all objects
    pTO = fpEvt->getTrgObjv2(i);
    int hltIndex = pTO->fHltIndex;

    bool lastModule= false;
    bool activeModule = false;
    ix=hltObjMap.find(i);
    //cout<<i<<" "<<hltIndex<<" ";
    //cout<<i<<" selected hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
    //	<<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber;
    if(ix!=hltObjMap.end()) {
      activeModule=true; // signal the teh module belongs a to selected/active HLT
      int num = ix->first;  // index of the module
      int hltN = (ix->second & 0x7FFFFFFF); // index of the HLT path the module belongs to
      lastModule = ( (ix->second & 0x80000000) != 0); // last module flag
      //cout<<" "<<num<<" "<<hltN<<" "<<lastN;
      if(num != i) cout<<" very very wrong1 "<<num<<" "<<i<<endl;
#ifdef OLD_OBJ_MARK
      if(hltN != (hltIndex%1000) ) cout<<" very very wrong2 "<<hltN<<" "<<hltIndex<<endl;
      if((hltIndex<1000) ) cout<<" very very wrong3 "<<hltIndex<<endl;
      if(lastModule != (hltIndex>1000000) ) cout<<" very very wrong4 "<<lastModule<<" "<<hltIndex<<endl;
#else
      if(hltN != hltIndex ) cout<<" very very wrong2 "<<hltN<<" "<<hltIndex<<endl;
#endif

    }
    //} else cout<<" nothing in map ";
    //cout<<endl;


#ifdef OLD_OBJ_MARK
    if( (hltIndex>1000) != activeModule)
      cout<<" very very wrong5 "<<activeModule<<" "<<hltIndex<<endl;
    if( (hltIndex>100000) != lastModule)
      cout<<" very very wrong6 "<<lastModule<<" "<<hltIndex<<endl;
#endif

    //  for(int i=0; i!=fpEvt->nTrgObjv2(); i++) { // loop over all objects
    double deltaR1=100, deltaR2=100;
    int idx1=-1, idx2=-1;

    //pTO = fpEvt->getTrgObjv2(i);  //pTO->dump();
    //int hltIndex = pTO->fHltIndex;
    //cout<<" selected hlt "<<i<<" "<<pTO->fHltPath<<" hlt-index "<<hltIndex
    //	<<" module label "<<pTO->fLabel<<endl;

    //if(activeModule ) { // this object was selected, or use all
    if(lastModule) { // this object was selected, use final last one
      if(localPrint)
	cout<<" selected hlt "<<i<<" "<<pTO->fHltPath<<" hlt-index "<<hltIndex
	    <<" module label "<<pTO->fLabel<<" type "<<pTO->fType
	    <<" number "<<pTO->fNumber<<" "<<activeModule<<endl;
      TString a = pTO->fHltPath;
      int hltId = getHLTId(a);
      int hltinfo= (hltId&0xF); // store in the lowest 4 bits

      // look only at level3 hlt
      level=3;
      if(a.Contains("L2"))      {if(localPrint) cout<<" L2 "<<a<<endl; level=2;}
      else if(a.Contains("L1")) {if(localPrint) cout<<" L1 "<<a<<endl; level=1;}
      if( level!=3 ) {
	cout<<" does not happen since i took L2Mu out "<<endl;
	hltinfo  = hltinfo|0x8000; // signal L1/2
	continue; // skip
      }
      count0++;


      vector<int> muonIndex = pTO->fIndex;
      vector<int> muonID = pTO->fID;
      vector<TLorentzVector> muonP = pTO->fP;
      int num = muonIndex.size();

      if(num<2) {
	cout<<" does it happen "<<endl;
	continue; // skip hlts with less than 2 muons
      }

      for(int n=0;n<num;++n) {  // loop over particles in this module, usually 2
	int index = muonIndex[n];
	int id = abs(muonID[n]);  // muon=13
	TLorentzVector p = muonP[n];

	// skip checking non-muon objects
	//if( (id != 13) ) continue; // NO some triggers are muon+track

	if(localPrint) {
	  cout<<" particle "<<n<<" index "<<index<<" id "<<id
	      <<" pt/eta/phi "<<p.Pt()<<" "<<p.Eta()<<" "<<p.Phi()<<endl;
	}

	// check direction matching
	int idx0=-1;
	double deltaR=100.;
	deltaR = findBestMatchingTrack(p,0,idx0); // look at all tracks, not just muons
	double tmp = findBestMatchingTrack(p,1,idum); // look at muons, just for testing
	if(id==13) { // for muons
	  ((TH1D*)fHistDir->Get("dr41"))->Fill(deltaR);
	  ((TH1D*)fHistDir->Get("dr42"))->Fill(tmp);
	} else {  // not muons, just tracks
	  ((TH1D*)fHistDir->Get("dr40"))->Fill(deltaR);
	}
	//if(localPrint) cout <<" dr "<<deltaR <<" idx "<<idx0<<" "<<i<<" "<<n<<endl;
	if(n>2) ((TH1D*)fHistDir->Get("dr43"))->Fill(deltaR);

	if(n==0) {
	  deltaR1=deltaR; idx1=idx0;
	  if(deltaR<deltaCut) hltinfo = hltinfo | 0x10; // matches a track
	  if(tmp<deltaCut)    hltinfo = hltinfo | 0x100; // matches a muon
	} else if(n==1) {
	  deltaR2=deltaR; idx2=idx0;
	  if(deltaR<deltaCut) hltinfo = hltinfo | 0x20;
	  if(tmp<deltaCut)    hltinfo = hltinfo | 0x200; // matches a muon
	} else {
	  if(deltaR<deltaCut) hltinfo = hltinfo | 0x40;
	  if(tmp<deltaCut)    hltinfo = hltinfo | 0x400; // matches a muon
	}


      } // end for loop n

      // save hltinfo for each last module in the hlt path
      hltInfo.push_back(hltinfo);  //

      // Now check if we match
      if( (deltaR1<deltaCut) && (deltaR2<deltaCut) ) {
	passed = true;
	count1++;
      } //else cout<<" matching failed,  dr1 "<<deltaR1 <<" idx "<<idx1
      //	 <<" dr2 "<<deltaR2 <<" idx "<<idx2<<" selected hlt "<<i<<" "
      //	 <<pTO->fHltPath<<" hlt-index "<<hltIndex
      //	 <<" module label "<<pTO->fLabel<<" type "<<pTO->fType
      //	 <<" number "<<pTO->fNumber<<count0<<" "<<count1<<endl;

      // select best in case several HLTs are present
      if (deltaR1<deltaRmin1) {
	deltaRmin1=deltaR1;
	hlt1 = pTO->fLabel;
	indx1=idx1;
	mu1match=i;
      } // if match
      if (deltaR2<deltaRmin2) {
	deltaRmin2=deltaR2;
	hlt2 = pTO->fLabel;
	indx2=idx2;
	mu2match=i;
      } // if match
      if(localPrint)
	cout <<" dr1 "<<deltaR1 <<" idx "<<idx1<<" dr2 "<<deltaR2 <<" idx "<<idx2<<" "<<i<<endl;

    } // end if valid module
  } // loop over all modules

  dr1=deltaRmin1;
  dr2=deltaRmin2;

  if (localPrint)
    cout << " best match "<<passed<<" "<<count0<<" "<<count1<<" "
	 <<indx1<<" "<< deltaRmin1 << " "<<mu1match<<" "<<hlt1<<" "
	 <<indx2<<" "<< deltaRmin2 << " "<<mu2match<<" "<<hlt2<<" "
	 <<endl;

  return passed;
}

//--------------------------------------------------------------
bool candAnaDstar::doTriggerVeto(TAnaTrack *fp, bool muonsOnly, bool matchPt,
			    bool anyModule, float deltaRthr, int histoOffset) {

  const double deltaPtMatch(0.30); // the pt matching cut 0.15
  const int verboseThr = 20;
  bool localPrint = (fVerbose==-32) || (fVerbose > verboseThr);
  //localPrint=true;
  //double trigMatchDeltaPtAll1 = 99.;
  int indx1=-1; //
  int mu1match(-1); //
  string hlt1; // hlt2;
  double deltaRmin1(99.);
  double trigMatchDeltaPt1 = 99.; //
  //int level=3; // select printout for L2 hlt-paths
  int modulesSelected=0, modulesSingleMatched=0;
  double drMinAll=99., drMin=99.; // best
  bool matchS=false;
  TTrgObjv2 *pTO;
  TLorentzVector tlvMu1; // tlvMu2;

  if (localPrint) {
    cout<<" doTriggerVeto "<<deltaRthr<<" "<<deltaPtMatch<<endl;
    cout << "pt,eta,phi: " << fp->fPlab.Perp() << " " << fp->fPlab.Eta() << " " << fp->fPlab.Phi()<< endl;
  }

  tlvMu1.SetPtEtaPhiM(fp->fPlab.Perp(),fp->fPlab.Eta(),fp->fPlab.Phi(),MMUON); // assume a muon

  //cout<<" part "<<histoOffset<<" hltInfo size "<<hltInfo.size()<<endl;
  unsigned int ninfo=0; // count HLTs

  // for(int i=0; i!=fpEvt->nTrgObjv2(); i++) { // loop over all modules
  //   // this includes modules from all valid (DS selecte) HLT paths, not split into HLT
  //   pTO = fpEvt->getTrgObjv2(i);
  //   int hltIndex = pTO->fHltIndex;


  //cout<<" size "<<hltObjMap.size()<<endl;
  map<unsigned int, unsigned int, less<unsigned int> >::iterator  ix;
  for(int i=0; i!=fpEvt->nTrgObjv2(); i++) { // loop over all objects
    pTO = fpEvt->getTrgObjv2(i);
    int hltIndex = pTO->fHltIndex;

    bool lastModule= false;
    bool activeModule = false;
    ix=hltObjMap.find(i);
    //cout<<i<<" "<<hltIndex<<" ";
    //cout<<i<<" selected hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
    //	<<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber;
    if(ix!=hltObjMap.end()) {
      activeModule=true; // signal the teh module belongs a to selected/active HLT
      int num = ix->first;  // index of the module
      int hltN = (ix->second & 0x7FFFFFFF); // index of the HLT path the module belongs to
      lastModule = ( (ix->second & 0x80000000) != 0); // last module flag
      //cout<<" "<<num<<" "<<hltN<<" "<<lastN;
      if(num != i) cout<<" very very wrong1 "<<num<<" "<<i<<endl;
#ifdef OLD_OBJ_MARK
     if(hltN != (hltIndex%1000) ) cout<<" very very wrong2 "<<hltN<<" "<<hltIndex<<endl;
      if((hltIndex<1000) ) cout<<" very very wrong3 "<<hltIndex<<endl;
      if(lastModule != (hltIndex>1000000) ) cout<<" very very wrong4 "<<lastModule<<" "<<hltIndex<<endl;
#else
     if(hltN != hltIndex ) cout<<" very very wrong2 "<<hltN<<" "<<hltIndex<<endl;
#endif
    }
    //} else cout<<" nothing in map ";
    //cout<<endl;

#ifdef OLD_OBJ_MARK
    if( (hltIndex>1000) != activeModule)
      cout<<" very very wrong5 "<<activeModule<<" "<<hltIndex<<endl;
    if( (hltIndex>100000) != lastModule)
      cout<<" very very wrong6 "<<lastModule<<" "<<hltIndex<<endl;
#endif

    //if(hltIndex>1000) { // this object was selected, matches our trigger list
    // anymodule==false should reproduce the old 2015 results
    if( lastModule || (anyModule&&activeModule) ) { // this object was selected,

      if( (lastModule) ) {
	if(ninfo>=hltInfo.size()) {cout<<" something wrong in hltInfo size "<<hltInfo.size()<<" "<<ninfo<<endl;}
	//else {cout<<ninfo<<" "<<hex<<hltInfo[ninfo]<<dec<<" "<<endl;}
      }

      if(localPrint) cout<<i<<" selected hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
			 <<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber<<endl;

      TString a = pTO->fHltPath;
      //double deltaRthr= deltaRthr0; // default cut
      //level=3;
      //if(a.Contains("L2") || a.Contains("L1") ) { // check if this is an L2 or L1 path
      //deltaRthr=deltaRthr1; // extend for L1&L2 objects
      //level=2;
      //if (localPrint)
      //  cout<<" L2/L1 object, extend dr cut to  "<<a<<" "<<deltaRthr<<endl;
      //}

      // reset the best resuts for each trigger module
      bool match1=false;// match2=false;
      int m1=-1; // m2=-1;
      deltaRmin1 = 99.; //deltaRmin2=99.;
      trigMatchDeltaPt1 = 99.; // trigMatchDeltaPt2 = 99.;

      vector<int> muonIndex = pTO->fIndex;
      vector<int> muonID = pTO->fID;
      vector<TLorentzVector> muonP = pTO->fP;
      int num = muonIndex.size();
      for(int n=0;n<num;++n) {  // loop over particles in this module, usually 2
	int index = muonIndex[n];
	int id = muonID[n];
	TLorentzVector p = muonP[n];

	if(localPrint)
	  cout<<"trg-track: pt/eta/phi "<<p.Pt()<<"/"<<p.Eta()<<"/"<<p.Phi()<<" i/n "<<i<<"/"<<n<<endl;

	// Do we do it? Can be a non-muon in the trigger
	if( muonsOnly && (abs(id) != 13) ) { // if not muon trigger skip if requested
	  if(fVerbose>1)
	    cout<<" a none hlt-muon found in a trigger object "
		<<n<<" id "<<id<<" "<<pTO->fHltPath<<" "<<pTO->fLabel<<" "
		<<pTO->fType<<" skip it "<<endl;
	  continue;  // skip checking non-muon objects
	}

	// check direction matching
	double deltaR1 = p.DeltaR(tlvMu1);
	if(localPrint) {
	  cout<<" particle "<<n<<" index "<<index<<" id "<<id
	      <<" pt/eta/phi "<<p.Pt()<<"/"<<p.Eta()<<"/"<<p.Phi()<<" i/n "<<i<<"/"<<n
	      <<" dr "<<deltaR1 <<endl;
	}

	if ( (histoOffset>0) && (histoOffset<4) ) {
	  if(histoOffset==1)      {
	    ((TH1D*)fHistDir->Get("test11"))->Fill(deltaR1);
	  } else if(histoOffset==2) {
	    ((TH1D*)fHistDir->Get("test12"))->Fill(deltaR1);
	  } else if(histoOffset==3) {
	    ((TH1D*)fHistDir->Get("test13"))->Fill(deltaR1);
	  }

	  // check if it is below threshold
	  //if( (deltaR1<deltaRthr) &&  (lastModule) ) { // use final hlt modules
	  if( (deltaR1<deltaRthr) &&  activeModule ) { // use all moduls in the path
	    //cout<<" matched part "<<histoOffset<<" "<<deltaR1<<endl;
	    int tmp=1;
	    tmp = tmp<<(histoOffset-1); // shift according to the particle
	    hltInfo[ninfo] = hltInfo[ninfo] | (tmp<<12);
	  }
	}

	if(deltaR1<deltaRmin1) {  // select best (smallest) dR for this particle
	  if(localPrint) {cout << " mu selected "<< deltaR1 <<endl;}

	  if (deltaR1<deltaRthr) {  // check if it is below threshold
	    // check now the p matching
	    double trigMatchDeltaPt=999.;
	    if (fp->fPlab.Mag() > 0.)
	      //trigMatchDeltaPt = TMath::Abs(p.Rho()  - fp1->fPlab.Mag())/fp1->fPlab.Mag();
	      trigMatchDeltaPt = (p.Pt()  - fp->fPlab.Perp())/fp->fPlab.Perp();
	    if( !matchPt || (TMath::Abs(trigMatchDeltaPt) < deltaPtMatch) ) {  // check if it is good enough
	      trigMatchDeltaPt1=(trigMatchDeltaPt);
	      deltaRmin1=deltaR1;  // best match until now
	      match1=true;
	      m1=n;
	    } else {if(localPrint) cout<<" pt1 match failed "<<trigMatchDeltaPt<<endl;} // if pt match
	  } else {if(localPrint) cout<<" dr1 too large "<<deltaR1<<endl; }// if delta
	} // if direction match

      } // end for loop n, tracks in a trig object

      if (localPrint)
	cout << " match for this module "
	     <<m1<<" "<< deltaRmin1 <<" "<<trigMatchDeltaPt1<<endl;  // best for thsi object

      //if(histoOffset==0)      {((TH1D*)fHistDir->Get("test14"))->Fill(deltaRmin1);}
      //else if(histoOffset==2) {((TH1D*)fHistDir->Get("test15"))->Fill(deltaRmin1);}
      //else if(histoOffset==3) {((TH1D*)fHistDir->Get("test16"))->Fill(deltaRmin1);}

      if(match1) {  // best for this module. compare with previous modules
	  // select the best batch
	  if(deltaRmin1<drMin) { // a better match, save it
	    drMin=deltaRmin1; // best among modules
	    mu1match = m1;
	    hlt1 = pTO->fLabel;
	    indx1=i;
	    //trigMatchDeltaPtAll1 = trigMatchDeltaPt1;
	    //if(histoOffset==0)       {((TH1D*)fHistDir->Get("test17"))->Fill(drMin);}
	    //else  if(histoOffset==2) {((TH1D*)fHistDir->Get("test18"))->Fill(drMin);}
	    //else  if(histoOffset==2) {((TH1D*)fHistDir->Get("test19"))->Fill(drMin);}
	  }  // if min
	  if(localPrint) cout<<" best matching for module "<<i<<" "<<drMin<<endl;
      }

      if( (lastModule) ) { // last module for thsi HLT
	// next module will be for the next HLT
	modulesSelected++;
	if(localPrint) cout<<" for this hlt "<<drMin<<" "<<indx1<<" "<<mu1match<<endl;
	if(drMin<deltaRthr) {  // check if it is below threshold
	  modulesSingleMatched++;
	  matchS=true;
	  if(drMin<drMinAll) drMinAll=drMin;
	  if(localPrint) cout<<" best dr "<<i<<" "<<drMinAll<<endl;
	}
	drMin=999.;
	++ninfo; // count triggers paths
      }

    } // end if a valid trigger module, i

  } // loop over all modules

  //if(histoOffset==0)       {((TH1D*)fHistDir->Get("test20"))->Fill(drMin);}
  //else  if(histoOffset==2) {((TH1D*)fHistDir->Get("test21"))->Fill(drMin);}
  //else  if(histoOffset==3) {((TH1D*)fHistDir->Get("test22"))->Fill(drMin);}

  bool veto = false;
  //if( matchS && (drMin<deltaRthr1) )  // check versus the worse case, L1/2
  if( matchS && (modulesSelected==modulesSingleMatched) )  // check versus the worse case, L1/2
    {veto=true; if(localPrint) cout<<" single veto "<<drMinAll<<endl;}
  //if( matchS )
  //{veto=true; if(localPrint) cout<<" single veto "<<drMinAll<<endl;}

  if (localPrint) {
    cout<<" veto = "<<veto<<" "<<matchS<< " best match "<< drMinAll;
    cout<<" modules "<<modulesSelected<<" "<<modulesSingleMatched<<endl;
  }

  // For testing only
  // if (modulesSelected <= 0) {
  //   if(TString(fHLTPath).Contains("Mu") || TString(fHLTPath).Contains("mu"))  {
  //     cout<<" error: no module found "<<fHLTPath<<endl;
  //     for(int i=0; i!=fpEvt->nTrgObjv2(); i++) { // loop over all objects
  // 	pTO = fpEvt->getTrgObjv2(i);
  // 	int hltIndex = pTO->fHltIndex;
  // 	//if( (lastModule) || (anyModule&&(activeModule)) ) { // this object was selected,
  // 	cout<<i<<" selected hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
  // 	    <<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber<<endl;
  //     }
  //   }
  // }


  // testing
  // if( (modulesSelected-1) != (fhltType/1000) )
  //   if(TString(fHLTPath).Contains("Mu") || TString(fHLTPath).Contains("mu"))
  //     cout<<" modules no right? "<<modulesSelected<<" "<<modulesMatched<<" "<<modulesSingleMatched<<" "<<veto<<" "<<fhltType<<endl;


  return veto;
}

//-------------------------------------------------------------------------------

// A trigger matcher based on deltaR (from Frank) + pt matching.
// check 2 muons, use only the selected hlt objects which correspond to triggers
// which passed and were on out trigger list.
// Only consider trig objects which match our trigger list.
// The main cuts are: deltaRthr for DR and deltaPtMatch for pt
// uses TTrgObjv2
// bool candAna::doTriggerMatching(TAnaTrack *fp1, TAnaTrack *fp2) { // call the normal version with (true)
//   int indx1=-1, indx2=-1;
//   const double deltaRthr(0.02); // final cut, Frank had 0.5, change 0.020
//   const double deltaPtMatch(0.15); // the pt matching cut
//   const int verboseThr = 30;
//   //const bool localPrint = false;
//   bool localPrint = (fVerbose==-32) || (fVerbose > verboseThr);
//   int mu1match(-1), mu2match(-1);
//   string hlt1, hlt2;
//   double deltaRmin1(100),deltaRmin2(100);
//   double trigMatchDeltaPt1 = 99., trigMatchDeltaPt2 = 99.;
//   bool match=false;
//   TTrgObjv2 *pTO;
//   TLorentzVector tlvMu1, tlvMu2;

//   if (localPrint) {
//     cout << "mu1: pt,eta,phi: " << fp1->fPlab.Perp() << " " << fp1->fPlab.Eta() << " " << fp1->fPlab.Phi()<< endl;
//     cout << "mu2: pt,eta,phi: " << fp2->fPlab.Perp() << " " << fp2->fPlab.Eta() << " " << fp2->fPlab.Phi()<< endl;
//   }

//   tlvMu1.SetPtEtaPhiM(fp1->fPlab.Perp(),fp1->fPlab.Eta(),fp1->fPlab.Phi(),MMUON); // assume a muon
//   tlvMu2.SetPtEtaPhiM(fp2->fPlab.Perp(),fp2->fPlab.Eta(),fp2->fPlab.Phi(),MMUON); // assume a muon

//   map<unsigned int, unsigned int, less<unsigned int> >::iterator  ix;
//   for(int i=0; i!=fpEvt->nTrgObjv2(); i++) { // loop over all objects
//     pTO = fpEvt->getTrgObjv2(i);
//     int hltIndex = pTO->fHltIndex;

//     bool lastModule= false;
//     bool activeModule = false;
//     ix=hltObjMap.find(i);
//     //cout<<i<<" selected hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
//     //	<<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber;
//     if (ix!=hltObjMap.end()) {
//       activeModule=true; // signal the teh module belongs a to selected/active HLT
//       lastModule = ( (ix->second & 0x80000000) != 0);
//     }

//     if (lastModule) { // this object was selected, use last module
//       if (localPrint) cout<<i<<" selected hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
// 			 <<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber<<" "<<activeModule<<endl;

//       bool match1=false, match2=false;
//       vector<int> muonIndex = pTO->fIndex;
//       vector<int> muonID = pTO->fID;
//       vector<TLorentzVector> muonP = pTO->fP;
//       int num = muonIndex.size();
//       for(int n=0;n<num;++n) {  // loop over particles in this module, usually 2
// 	int index = muonIndex[n];
// 	int id = muonID[n];
// 	TLorentzVector p = muonP[n];

// 	if ( abs(id) != 13 ) { // if not muon trigger skip
// 	  if (localPrint) cout<<" a none hlt-muon found in a trigger object, skip it, id= "
// 			     <<id<<" "<<pTO->fHltPath<<" "<<pTO->fLabel<<" "<<pTO->fType<<endl;
// 	  continue;  // skip checking non-muon objects
// 	}

// 	// check direction matching
// 	double deltaR1 = p.DeltaR(tlvMu1);
// 	double deltaR2 = p.DeltaR(tlvMu2);

// 	if (localPrint) {
// 	  cout<<" particle"<<n<<" index "<<index<<" id "<<id
// 	      <<" pt/eta/phi "<<p.Pt()<<"/"<<p.Eta()<<"/"<<p.Phi()<<" i/n "<<i<<"/"<<n
// 	      <<" dr "<<deltaR1 <<" "<<deltaR2<<endl;
// 	}

// 	// muon 1
// 	if (deltaR1<deltaRmin1) {
// 	  deltaRmin1=deltaR1;  // best match until now
// 	  if (fVerbose > verboseThr || localPrint) {cout << " mu1 selected "<< deltaR1 <<endl;}
// 	    // check now the pt matching
// 	  double trigMatchDeltaPt=999.;
// 	  if (fp1->fPlab.Mag() > 0.) trigMatchDeltaPt = TMath::Abs(p.Rho()  - fp1->fPlab.Mag())/fp1->fPlab.Mag();
// 	  if ( trigMatchDeltaPt < deltaPtMatch ) {  // check if it is good enough
// 	    if (deltaR1<deltaRthr) {
// 	      trigMatchDeltaPt1=trigMatchDeltaPt;
// 	      mu1match = n;
// 	      hlt1 = pTO->fLabel;
// 	      indx1=i;
// 	      match1=true;
// 	    } // if delta
// 	  } // if pt match
// 	} // if direction match

// 	// muon 2
// 	if (deltaR2<deltaRmin2) {
// 	  deltaRmin2=deltaR2;
// 	  if (localPrint) {cout << " mu2 selected "<< deltaR2 <<endl;}
// 	    // check now the pt matching
// 	  double trigMatchDeltaPt=999.;
// 	  if (fp2->fPlab.Mag() > 0.) trigMatchDeltaPt = TMath::Abs(p.Rho()  - fp2->fPlab.Mag())/fp2->fPlab.Mag();
// 	  if ( trigMatchDeltaPt < deltaPtMatch ) {
// 	    if (deltaR2<deltaRthr) {
// 	      trigMatchDeltaPt2=trigMatchDeltaPt;
// 	      mu2match = n;
// 	      hlt2 = pTO->fLabel;
// 	      indx2=i;
// 	      match2=true;
// 	    } // if delta
// 	  } // if pt match
// 	} // if direction match
//       } // end for loop n

//       // check that at least one module matched both
//       match = match || (match1&&match2);

//     } // end if valid module

//   } // loop over all modules

//   if (localPrint)
//     cout << " best match "
// 	 <<indx1<<" "<< deltaRmin1 << " "<<mu1match<<" "<<hlt1<<" "<<trigMatchDeltaPt1<<" "
// 	 <<indx2<<" "<< deltaRmin2 << " "<<mu2match<<" "<<hlt2<<" "<<trigMatchDeltaPt2<<endl;

//   ((TH1D*)fHistDir->Get("test8"))->Fill(trigMatchDeltaPt1);
//   ((TH1D*)fHistDir->Get("test8"))->Fill(trigMatchDeltaPt2);
//   ((TH1D*)fHistDir->Get("test2"))->Fill(deltaRmin1);
//   ((TH1D*)fHistDir->Get("test2"))->Fill(deltaRmin2);

//   if (mu1match>-1) {
//     double tmp=fMu1TrigM;
//     fMu1TrigM = deltaRmin1;
//     if (tmp!=fMu1TrigM) cout<<"Warning:  two methods inconsistent-mu1 "<<tmp<<" "<<fMu1TrigM<<endl;
//   }
//   if (mu2match>-1) {
//     double tmp=fMu2TrigM;
//     fMu2TrigM = deltaRmin2;
//     if (tmp!=fMu2TrigM) cout<<"Warning:  two methods inconsistent-mu2 "<<tmp<<" "<<fMu2TrigM<<endl;
//   }

//   bool HLTmatch = false;
//   if (match && mu1match>-1 && mu2match>-1) {
//     if (indx1!=indx2) { // should never happen since usually we only have one selected trigger
//       cout<<"Warning:  best match for the two muons is to two different modules "<<indx1<<" "<<indx2<<endl;
//       HLTmatch=true;
//     } else if (mu1match==mu2match) { // matched to 2 same tracks, what to do? skip it?
//       cout<<"Error:  two muons matched to same particle "<<indx1<<" "<<indx2<<" "<<mu1match<<endl;
//     } else { // ok
//       if (localPrint) cout<<" matching OK"<<indx1<<endl;
//       HLTmatch=true;
//     }
//   }

//   return HLTmatch;
// }
//-------------------------------------------------------------------------------
// match track to trigger, return DR
// return the best, smalles DR
// anyTrig = true - match to any triggered which fired & created a muon object (might be from another DS)
//         = false - match only to the trigger selected from our list
// muonsOnly = true - uae only muon trigger particles to match
//           = false - use all trigger particles
double candAnaDstar::doTriggerMatchingR(TAnaTrack *fp1, bool anyTrig, bool muonsOnly, bool anyModule) {
  int indx1=-1;
  //const double deltaRthr(0.02); // final cut, Frank had 0.5, change 0.020
  const double deltaPtMatch(0.15); // the pt matching cut
  const int verboseThr = 30;
  //const bool localPrint = false;
  bool localPrint = (fVerbose==-32) || (fVerbose > verboseThr);
  int mu1match(-1);
  string hlt1;
  double deltaRmin1(100);
  double trigMatchDeltaPt1 = 99.;

  TTrgObjv2 *pTO;
  TLorentzVector tlvMu1;

  if (localPrint) {
    cout << " doTriggerMatchingR pt,eta,phi: " << fp1->fPlab.Perp() << " " << fp1->fPlab.Eta() << " " << fp1->fPlab.Phi()<< endl;
  }

  tlvMu1.SetPtEtaPhiM(fp1->fPlab.Perp(),fp1->fPlab.Eta(),fp1->fPlab.Phi(),MMUON); // assume a muon

  //cout<<" size "<<hltObjMap.size()<<endl;
  map<unsigned int, unsigned int, less<unsigned int> >::iterator  ix;
  for(int i=0; i!=fpEvt->nTrgObjv2(); i++) { // loop over all objects
    pTO = fpEvt->getTrgObjv2(i);
    int hltIndex = pTO->fHltIndex;

    bool lastModule= false;
    bool activeModule = false;
    ix=hltObjMap.find(i);
    //cout<<i<<" selected hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
    //	<<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber;
    if (ix!=hltObjMap.end()) {
      activeModule=true; // signal the teh module belongs a to selected/active HLT
      int num = ix->first;  // index of the module
      int hltN = (ix->second & 0x7FFFFFFF); // index of the HLT path the module belongs to
      lastModule = ( (ix->second & 0x80000000) != 0); // last module flag
      //cout<<" "<<num<<" "<<hltN<<" "<<lastN;
      if (num != i) cout<<" very very wrong1 "<<num<<" "<<i<<endl;
      if (hltN != hltIndex ) cout<<" very very wrong2 "<<hltN<<" "<<hltIndex<<endl;

    }

    if (anyTrig || lastModule || (anyModule&&activeModule) ) { // this object was selected, or use all
      if (localPrint) cout<<i<<" selected hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
			 <<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber<<endl;

      vector<int> muonIndex = pTO->fIndex;
      vector<int> muonID = pTO->fID;
      vector<TLorentzVector> muonP = pTO->fP;
      int num = muonIndex.size();
      for(int n=0;n<num;++n) {  // loop over particles in this module, usually 2
	int index = muonIndex[n];
	int id = muonID[n];
	TLorentzVector p = muonP[n];

	if ( abs(id) != 13 ) { // if not muon trigger skip
	  //cout<<" match track "<<id<<" "<<pTO->fHltPath<<" "<<pTO->fLabel<<" "<<pTO->fType<<endl;
	  if (muonsOnly) continue;  // skip checking non-muon objects
	}

	// check direction matching
	double deltaR1 = p.DeltaR(tlvMu1);

	if (localPrint) {
	  cout<<" particle"<<n<<" index "<<index<<" id "<<id
	      <<" pt/eta/phi "<<p.Pt()<<" "<<p.Eta()<<" "<<p.Phi()<<endl;
	  cout <<" mu1 "<< i<<" "<<pTO->fLabel << " "<<n <<" "<<deltaR1 <<endl;
	}

	// muon 1
	if (deltaR1<deltaRmin1) {
	  deltaRmin1=deltaR1;
	  mu1match = n;
	  hlt1 = pTO->fLabel;
	  indx1=i;

	  if (fVerbose > verboseThr || localPrint) {cout << " mu1 selected "<< deltaR1 <<endl;}

	  // check now the pt matching
	  double trigMatchDeltaPt=999.;
	  if (fp1->fPlab.Mag() > 0.) trigMatchDeltaPt = TMath::Abs(p.Rho()  - fp1->fPlab.Mag())/fp1->fPlab.Mag();
	  if ( trigMatchDeltaPt < deltaPtMatch ) trigMatchDeltaPt1=trigMatchDeltaPt;

	} // if direction match

      } // end for loop n
    } // end if valid module

  } // loop over all modules


  if (localPrint)
    cout << " best match "
	 <<indx1<<" "<< deltaRmin1 << " "<<mu1match<<" "<<hlt1<<" "<<trigMatchDeltaPt1<<endl;

  fTrigMatchDeltaPt=trigMatchDeltaPt1;

  if (mu1match>-1) {
    if (localPrint) cout<<" matching OK "<<deltaRmin1<<" "<<hlt1<<endl;
    //double tmp=fMu1TrigM;
    //fMu1TrigM = deltaRmin1;
    //if (tmp!=fMu1TrigM) cout<<" two methods inconsistent "<<tmp<<" "<<fMu1TrigM<<endl;
  }

  return deltaRmin1;
}

// // ---------------------------------------------------------------------------------
// // To match a single track to a trigger object (selected or all)
// // pt - track
// // anyTrig - if true use all trigger objects
// // calls doTriggerMatchingR()
//   bool candAna::doTriggerMatching(TAnaTrack *pt, bool anyTrig, bool muonsOnly, bool anyModule) {

//   const double deltaRthrsh(0.02); // final cut, Frank had 0.5, change 0.020

//   double dR = doTriggerMatchingR(pt,anyTrig,muonsOnly,anyModule);
//   ((TH1D*)fHistDir->Get("test6"))->Fill(dR);

//   bool HLTmatch = (dR<deltaRthrsh );
//   return HLTmatch;

// }

// //--------------------------------------------------------------
// bool candAna::doTriggerVeto(TAnaTrack *fp, bool muonsOnly, bool matchPt,
// 			    bool anyModule, float deltaRthr, int histoOffset) {

//   cout<<" OBSOLETE, DO NOT USE "<<endl;

//   // The valid code is in CandAnaDstar.cc
//   return false;

// }




//-------------------------
void candAnaDstar::genMatch() {
  //cout<<" in candAnaDstar::genMatch "<<endl;
  return;
}
void candAnaDstar::recoMatch() {
  return;
}
void candAnaDstar::candMatch() {
  return;
}
void candAnaDstar::efficiencyCalculation() {
  return;
}
