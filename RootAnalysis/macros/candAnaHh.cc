#include "candAnaHh.hh"
#include <cmath>
#include <string>

#include "common/HFMasses.hh"
#include "danekUtils.h"

//#define MYCODE

using namespace std;

namespace {
  TVector3 BdVertexGen(0,0,0), PVGen(0,0,0);
  TVector3 BdMomGen(0,0,0), Pi1MomGen(0,0,0), Pi2MomGen(0,0,0);
  const bool MYDEBUG = false;
}

// ----------------------------------------------------------------------
candAnaHh::candAnaHh(bmmReader *pReader, std::string name, std::string cutsFile) :
  candAna(pReader, name, cutsFile) {
  cout << "==> candAnaHh: name = " << name << ", reading cutsfile " << cutsFile << endl;
  BLIND = 0;
  readCuts(cutsFile, 1);

}


// ----------------------------------------------------------------------
candAnaHh::~candAnaHh() {
  cout << "==> candAnaHh: destructor..." << endl;
  tree->Write();
}


// ----------------------------------------------------------------------
// To analyze the MC event
bool candAnaHh::anaMC(TAna01Event *evt) {
  const bool print = false;

  fpEvt = evt;
  bool foundPV=false, foundBd = false, foundPi1 = false, foundPi2=false;
  int pC0 = 0;


  int numGenCands = fpEvt->nGenCands();
  if(print) cout << "Found " << numGenCands << " gen cands in event" << endl;
  for (int it = 0; it < numGenCands; ++it) {  // loop over all gen candidates

    foundBd = false, foundPi1 = false, foundPi2=false;

    TGenCand * pCand = fpEvt->getGenCand(it);

    if( !foundPV && (abs(pCand->fID) == 5) ) {foundPV=true; PVGen=pCand->fV; if(print) cout<<" PV "<<PVGen.Z()<<endl;}   // get PV

    if( ( abs(pCand->fID) != 511) ) continue;        // skip others

    if(print) cout <<" Bd "<<it<< " " << pCand->fNumber << " "<<pCand->fID<<" "<<pCand->fQ<<" "<<pCand->fStatus<<" "
		   <<pCand->fMom1<<" "<<pCand->fMom2<<" "<<pCand->fDau1<<" "
		   <<pCand->fDau2<<" "<<pCand->fP.Perp()<<" "<<pCand->fV.Z()<<endl;

    BdMomGen = (pCand->fP.Vect());
    foundBd = true;
    pC0 = it;

    // Look at daugthers
    int i1 = (pCand->fDau2)-(pCand->fDau1)+1;
    int i2=0;
    //if(i1!=2) {continue;} // fpEvt->dumpGenBlock();}
    if(i1!=2) { if(print) cout<<" number of daughters wrong skip "<<i1<<" "<<pC0<<endl; continue;} // fpEvt->dumpGenBlock();}

    for(int id=(pCand->fDau1);id<=(pCand->fDau2);++id) {
      TGenCand * dau = fpEvt->getGenCand(id);  // check daughters

      if( abs(dau->fID) == 211 ) { //  pions
	if(print) cout <<" Pi "<<dau->fNumber << " "<<dau->fID<<" "<<dau->fQ<<" "
		       <<dau->fMom1<<" "<<dau->fMom2<<" "<<dau->fDau1<<" "
		       <<dau->fDau2<<" "<<dau->fP.Perp()<<" "<<dau->fV.Z()<<endl;

	BdVertexGen = (dau->fV);  // Bd decay vertex

	i2++;
	if(i2==1)      {foundPi1=true; Pi1MomGen = (dau->fP.Vect());}
	else if(i2==2) {foundPi2=true; Pi2MomGen = (dau->fP.Vect());}

      }

    } // daugther loop

      if(foundBd && foundPi1 && foundPi2) break;

  } // gen part loop


  bool ok = foundPV && foundBd && foundPi1 && foundPi2;
  if(ok) {

    TVector3 t1(BdVertexGen-PVGen);
    double a0 = t1.Angle(BdMomGen);  // Bd pointing angle

    ((TH1D*)fHistDir->Get("h11"))->Fill(t1.Mag());
    ((TH1D*)fHistDir->Get("h20"))->Fill(a0);

    ((TH1D*)fHistDir->Get("h28"))->Fill(BdMomGen.Perp());
    ((TH1D*)fHistDir->Get("h29"))->Fill(Pi1MomGen.Perp());
    ((TH1D*)fHistDir->Get("h30"))->Fill(Pi2MomGen.Perp());

    double dr    = danekUtils::twoBodyDecayAngle(Pi1MomGen, Pi2MomGen);  // openning angle
    //double pt    = danekUtils::twoBodyDecayMomPerp(Pi1MomGen, Pi2MomGen);
    //double m1    = danekUtils::twoBodyDecayMass(Pi1MomGen, Pi2MomGen, MPION, MPION);
    //double m2    = danekUtils::twoBodyDecayMass(Pi1MomGen, Pi2MomGen, MMUON, MMUON);

    ((TH1D*)fHistDir->Get("h47"))->Fill(dr);
    //((TH1D*)fHistDir->Get("h41"))->Fill(pt);
    //((TH1D*)fHistDir->Get("h42"))->Fill(m2);
    //((TH1D*)fHistDir->Get("h44"))->Fill(m1);

    double a1 = t1.Angle(Pi1MomGen);  // pi1 angle
    double a2 = t1.Angle(Pi2MomGen);  // pi1 angle
    double perp1 = Pi1MomGen.Mag() * sin(a1);
    double perp2 = Pi2MomGen.Mag() * sin(a2);
    double perpdiff = perp1 - perp2;
    double perpsum = perp1 + perp2;

    ((TH1D*)fHistDir->Get("h12"))->Fill(a1);
    ((TH1D*)fHistDir->Get("h12"))->Fill(a2);
    ((TH1D*)fHistDir->Get("h13"))->Fill(perp1);
    ((TH1D*)fHistDir->Get("h13"))->Fill(perp2);
    ((TH1D*)fHistDir->Get("h14"))->Fill(perpdiff);
    ((TH1D*)fHistDir->Get("h15"))->Fill(perpsum);
    ((TH2D*)fHistDir->Get("h16"))->Fill(a1,a2);
    ((TH2D*)fHistDir->Get("h17"))->Fill(perp1,perp2);


  } else {
    if(print) {cout<<" Bd->pipi not found in Gen"<<endl; fpEvt->dumpGenBlock(); }

  } // if OK

  return ok;
}


//------------------------------------------------------------------------------------------
void candAnaHh::evtAnalysis(TAna01Event *evt) {

  if(MYDEBUG) cout<<" candAnaHh::evtAnalysis() event =  "<<fEvt<<endl;
  candAna::evtAnalysis(evt);

  if(MYDEBUG) candAna::play2();

#ifdef MYCODE

  fpEvt = evt;
  fcands=0;
  hhAnalysis();

#endif

  return;
}
// ----------------------------------------------------------------------
// this runs for each candidate
void candAnaHh::candAnalysis() {

  //if(MYDEBUG) cout<<" candAnaHh::candAnalysis() "<<endl;
  if (0 == fpCand) return;

  //if(MYDEBUG) cout<<" Call candAna::candAnalysis: "<<fpEvt<<" "<<fpCand<<endl;

  candAna::candAnalysis();

  TAnaTrack *p1 = fpEvt->getSigTrack(fpCand->fSig1);
  TAnaTrack *p2 = fpEvt->getSigTrack(fpCand->fSig2);

  if(MYDEBUG) {
    cout<< "IN Hh =================================== "<<endl;
    cout<< " signal track pointers "<<p1 <<" "<<p2<<" "<<p1->fIndex<<" "<<p2->fIndex<<endl;
    p1->dump();
    p2->dump();
  }

  //TAnaTrack *pPi1 = fpEvt->getRecTrack(p1->fIndex);
  //TAnaTrack *pPi2 = fpEvt->getRecTrack(p2->fIndex);
  //if(MYDEBUG) cout<<" rec track pointers "<<pPi1<<" "<<pPi2<<endl;

  TSimpleTrack *ps1 = fpEvt->getSimpleTrack(p1->fIndex); // get simple track with index i
  TSimpleTrack *ps2 = fpEvt->getSimpleTrack(p2->fIndex); // get simple track with index i
  if(MYDEBUG) {
    cout<<" simple track pointers "<<ps1<<" "<<ps2
	<< " index "<< ps1->getIndex() << " " << ps2->getIndex() << endl;
    //<< " with ID = " << fpEvt->getSimpleTrackMCID(ps->getIndex())
    //      <<" "<<cIdx[i]<<" "<<(ps->getP()).Perp()
    //      <<  endl;

    ps1->dump();
    ps2->dump();
  }

  int mid1 = 0, mid2= 0;
  if(MYDEBUG) cout<< " mu-index "<<p1->fMuIndex <<" "<<p2->fMuIndex <<endl;

  if (p1->fMuIndex > -1) mid1 = fpEvt->getMuon(p1->fMuIndex)->fMuID;
  if (p2->fMuIndex > -1) mid2 = fpEvt->getMuon(p2->fMuIndex)->fMuID;

  ((TH1D*)fHistDir->Get("testhh0"))->Fill(1.);

  // veto global and tracker muons
  bool muonid1 = ((mid1&0x6)!=0); // 1-muon, 0-not muon
  bool muonid2 = ((mid2&0x6)!=0); //
  bool antimuon_veto = muonid1 || muonid2; // 1- means muon, ignore, 0 - no muon, pass
  if(MYDEBUG) cout<< " mu id "<< mid1 <<" "<< mid2 <<" anti mu veto "<<antimuon_veto<<endl;
  if(muonid1) ((TH1D*)fHistDir->Get("testhh0"))->Fill(3.);
  if(muonid2) ((TH1D*)fHistDir->Get("testhh0"))->Fill(4.);
  if(antimuon_veto) ((TH1D*)fHistDir->Get("testhh0"))->Fill(2.);

  //bool muon_veto =  true; // 1- means no muon, passed, 0 - muon, ignore
  // veto muons with any bit
  //bool muon_veto = (mid1==0) && (mid2==0); // 1- means no muon, passed, 0 - muon, ignore
  // veto all bits except ecal (CAL BITS ARE NEVER SET?)
  //bool muon_veto = ( (mid1&0x7FFF) == 0) && ( (mid2&0x7FFF) == 0); // 1- means no muon, passed, 0 - muon, ignore
  // veto global and tracker muons and standalone muons
  //bool muon_veto = ( (mid1&0x7) == 0) && ( (mid2&0x7) == 0); // 1- means no muon, passed, 0 - muon, ignore
  //bool muon_veto = !tightMuon(pPi1) && !tightMuon(pPi2); // 1- means no muon, passed, 0 - muon, ignore
  //bool muon_veto = !goodMuon(pPi1) && !goodMuon(pPi2); // 1- means no muon, passed, 0 - muon, ignore
//   if( ((mid1&0x8000) != 0) || ((mid2&0x8000) != 0) ) cout<<hex<<mid1<<" "<<mid2<<dec<<endl;
//    if( muon_veto1 != muon_veto12 ) {
//       cout<<fPreselection<<" "<<muon_veto1<<" "<<muon_veto11<<" "<<muon_veto12<<" "<<muon_veto2<<" "<<muon_veto3
//  	 << " muid "<<hex<<mid1 << " .. "<<mid2<<dec<<" ";
//       cout<<goodMuon(pPi1)<<" "<<tightMuon(pPi1)<<" "<<goodMuon(pPi2)<<" "<<tightMuon(pPi2)<<endl;
//    }

  // checkmatching
  //                         track allTrig useMuonOnly
  float dr1 = doTriggerMatchingR(p1, false, true);
  float dr2 = doTriggerMatchingR(p2, false, true);
  ((TH1D*)fHistDir->Get("testhh11"))->Fill(dr1);
  ((TH1D*)fHistDir->Get("testhh11"))->Fill(dr2);
  if(MYDEBUG) cout<< " DR !All/Muons "<< dr1 <<" "<< dr2 <<endl;

  dr1 = doTriggerMatchingR(p1, false, false);
  dr2 = doTriggerMatchingR(p2, false, false);
  ((TH1D*)fHistDir->Get("testhh12"))->Fill(dr1);
  ((TH1D*)fHistDir->Get("testhh12"))->Fill(dr2);
  if(MYDEBUG) cout<< " DR !All/!Muons "<< dr1 <<" "<< dr2 <<endl;

  dr1 = doTriggerMatchingR(p1, true, true);
  dr2 = doTriggerMatchingR(p2, true, true);
  ((TH1D*)fHistDir->Get("testhh13"))->Fill(dr1);
  ((TH1D*)fHistDir->Get("testhh13"))->Fill(dr2);
  if(MYDEBUG) cout<< " DR All/Muons "<< dr1 <<" "<< dr2 <<endl;

  dr1 = doTriggerMatchingR(p1, true, false);
  dr2 = doTriggerMatchingR(p2, true, false);
  ((TH1D*)fHistDir->Get("testhh14"))->Fill(dr1);
  ((TH1D*)fHistDir->Get("testhh14"))->Fill(dr2);
  if(MYDEBUG) cout<< " DR All/!Muons "<< dr1 <<" "<< dr2 <<endl;
  //                            singleMatch muonsOnly matchPt
  bool veto1 = doTriggerVeto(p1,p2,false,true,true); //
  bool veto2 = doTriggerVeto(p1,p2,true,true,true); // use this, 1 track in trigger vetos the event
  if(!veto1) ((TH1D*)fHistDir->Get("testhh0"))->Fill(5.); // count accepeted candidates
  if(!veto2) ((TH1D*)fHistDir->Get("testhh0"))->Fill(6.);

  if(MYDEBUG) cout<< " trigger veto (double) "<< veto1 <<" (single) "<<veto2<<endl;

  // veto global and tracker muons and standalone muons
  bool muon11 = ( (mid1&0x7) != 0);
  bool muon12 = ( (mid2&0x7) != 0);
  // Thight muon
  bool muon21 = tightMuon(p1);
  bool muon22 = tightMuon(p2);

  if(MYDEBUG) cout<< " muon-is "
		  << muon11 <<"/"<<muon12<<" "
		  << muon21 <<"/"<<muon22<<" "
		//<< muon31 <<"/"<<muon32<<" "
		  <<endl;

  if(muon11) ((TH1D*)fHistDir->Get("testhh0"))->Fill(10.); // count accepeted candidates
  if(muon12) ((TH1D*)fHistDir->Get("testhh0"))->Fill(20.); // count accepeted candidates
  if(muon21) ((TH1D*)fHistDir->Get("testhh0"))->Fill(11.); // count accepeted candidates
  if(muon22) ((TH1D*)fHistDir->Get("testhh0"))->Fill(21.); // count accepeted candidates
  //if(muon31) ((TH1D*)fHistDir->Get("testhh0"))->Fill(12.); // count accepeted candidates
  //if(muon31) ((TH1D*)fHistDir->Get("testhh0"))->Fill(22.); // count accepeted candidates
  if(!veto1){
  if(muon11) ((TH1D*)fHistDir->Get("testhh0"))->Fill(13.); // count accepeted candidates
  if(muon12) ((TH1D*)fHistDir->Get("testhh0"))->Fill(23.); // count accepeted candidates
  if(muon21) ((TH1D*)fHistDir->Get("testhh0"))->Fill(14.); // count accepeted candidates
  if(muon22) ((TH1D*)fHistDir->Get("testhh0"))->Fill(24.); // count accepeted candidates
  //if(muon31) ((TH1D*)fHistDir->Get("testhh0"))->Fill(15.); // count accepeted candidates
  //if(muon32) ((TH1D*)fHistDir->Get("testhh0"))->Fill(25.); // count accepeted candidates
  }
  if(!veto2){
  if(muon11) ((TH1D*)fHistDir->Get("testhh0"))->Fill(16.); // count accepeted candidates
  if(muon12) ((TH1D*)fHistDir->Get("testhh0"))->Fill(26.); // count accepeted candidates
  if(muon21) ((TH1D*)fHistDir->Get("testhh0"))->Fill(17.); // count accepeted candidates
  if(muon22) ((TH1D*)fHistDir->Get("testhh0"))->Fill(27.); // count accepeted candidates
  //if(muon31) ((TH1D*)fHistDir->Get("testhh0"))->Fill(18.); // count accepeted candidates
  //if(muon32) ((TH1D*)fHistDir->Get("testhh0"))->Fill(28.); // count accepeted candidates
  }



  if( fpCand->fMass<HH_MLO || fpCand->fMass>HH_MHI ) fPreselection = 0;
  //if( MUON_VETO==1 ) fPreselection = fPreselection && !antimuon_veto;

  return;

}
// ----------------------------------------------------------------------
// returns
//  1 - veto=true : the di-track pair is in the trigger, it triggered the event
//  0 - veto=false : there is a trigger which is not associated with the di-track
// bool candAnaHh::doTriggerVeto(TAnaTrack *fp1, TAnaTrack *fp2, bool singleMatch) {
//   const double deltaRthr(0.02); // final cut, Frank had 0.5, change 0.020
//   const double deltaPtMatch(0.15); // the pt matching cut
//   const int verboseThr = 20;
//   bool localPrint = (fVerbose==-32) || (fVerbose > verboseThr);
//   if(MYDEBUG) localPrint=true;

//   double deltaRminAll1(99.),deltaRminAll2(99.);
//   double trigMatchDeltaPtAll1 = 99., trigMatchDeltaPtAll2 = 99.;

//   int indx1=-1, indx2=-1;
//   int mu1match(-1), mu2match(-1);
//   string hlt1, hlt2;
//   double deltaRmin1(99.),deltaRmin2(99.);
//   double trigMatchDeltaPt1 = 99., trigMatchDeltaPt2 = 99.;

//   int modulesSelected=0, modulesMatched=0, modulesSingleMatched=0;
//   double drMin=99.;
//   bool match = false, matchS=false;
//   TTrgObjv2 *pTO;
//   TLorentzVector tlvMu1, tlvMu2;

//   if (localPrint) {
//     cout << "1: pt,eta,phi: " << fp1->fPlab.Perp() << " " << fp1->fPlab.Eta() << " " << fp1->fPlab.Phi()<< endl;
//     cout << "2: pt,eta,phi: " << fp2->fPlab.Perp() << " " << fp2->fPlab.Eta() << " " << fp2->fPlab.Phi()<< endl;
//   }

//   tlvMu1.SetPtEtaPhiM(fp1->fPlab.Perp(),fp1->fPlab.Eta(),fp1->fPlab.Phi(),MMUON); // assume a muon
//   tlvMu2.SetPtEtaPhiM(fp2->fPlab.Perp(),fp2->fPlab.Eta(),fp2->fPlab.Phi(),MMUON); // assume a muon

//   for(int i=0; i!=fpEvt->nTrgObjv2(); i++) { // loop over all objects
//     pTO = fpEvt->getTrgObjv2(i);
//     //pTO->dump();
//     int hltIndex = pTO->fHltIndex;
//     if(hltIndex>1000) { // this object was selected, matches our trigger list
//       if(localPrint) cout<<i<<" selected hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
// 			 <<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber<<endl;
//       modulesSelected++;

//       // reset the best resuts for each trigger module
//       bool match1=false, match2=false;
//       int m1=-1, m2=-1;
//       deltaRmin1 = 99.; deltaRmin2=99.;
//       trigMatchDeltaPt1 = 99.; trigMatchDeltaPt2 = 99.;

//       vector<int> muonIndex = pTO->fIndex;
//       vector<int> muonID = pTO->fID;
//       vector<TLorentzVector> muonP = pTO->fP;
//       int num = muonIndex.size();
//       for(int n=0;n<num;++n) {  // loop over particles in this module, usually 2
// 	int index = muonIndex[n];
// 	int id = muonID[n];
// 	TLorentzVector p = muonP[n];

// 	if(localPrint)
// 	  cout<<"trg-track: pt/eta/phi "<<p.Pt()<<"/"<<p.Eta()<<"/"<<p.Phi()<<" i/n "<<i<<"/"<<n<<endl;

// 	// Do we do it? Can be a non-muon in the trigger
// 	if( abs(id) != 13 ) { // if not muon trigger skip
// 	  if(fVerbose>1)
// 	    cout<<" a none hlt-muon found in a trigger object "
// 		<<n<<" id "<<id<<" "<<pTO->fHltPath<<" "<<pTO->fLabel<<" "
// 		<<pTO->fType<<" skip it "<<endl;
// 	  continue;  // skip checking non-muon objects
// 	}

// 	// check direction matching
// 	double deltaR1 = p.DeltaR(tlvMu1);
// 	double deltaR2 = p.DeltaR(tlvMu2);

// 	if(localPrint) {
// 	  cout<<" particle "<<n<<" index "<<index<<" id "<<id
// 	      <<" pt/eta/phi "<<p.Pt()<<"/"<<p.Eta()<<"/"<<p.Phi()<<" i/n "<<i<<"/"<<n
// 	      <<" dr "<<deltaR1 <<" "<<deltaR2<<endl;
// 	}

// 	// muon 1
// 	if(deltaR1<deltaRmin1) {
// 	  deltaRmin1=deltaR1;  // best match until now
// 	  if (fVerbose > verboseThr || localPrint) {cout << " mu1 selected "<< deltaR1 <<endl;}
// 	    // check now the pt matching
// 	  double trigMatchDeltaPt=999.;
// 	  if (fp1->fPlab.Mag() > 0.) trigMatchDeltaPt = TMath::Abs(p.Rho()  - fp1->fPlab.Mag())/fp1->fPlab.Mag();
// 	  if( trigMatchDeltaPt < deltaPtMatch ) {  // check if it is good enough
// 	    if (deltaR1<deltaRthr) {
// 	      trigMatchDeltaPt1=trigMatchDeltaPt;
// 	      match1=true;
// 	      m1=n;
// 	    } // if delta
// 	  } // if pt match
// 	} // if direction match

// 	// muon 2
// 	if(deltaR2<deltaRmin2) {
// 	  deltaRmin2=deltaR2;
// 	  if (localPrint) {cout << " mu2 selected "<< deltaR2 <<endl;}
// 	    // check now the pt matching
// 	  double trigMatchDeltaPt=999.;
// 	  if (fp2->fPlab.Mag() > 0.) trigMatchDeltaPt = TMath::Abs(p.Rho()  - fp2->fPlab.Mag())/fp2->fPlab.Mag();
// 	  if( trigMatchDeltaPt < deltaPtMatch ) {
// 	    if (deltaR2<deltaRthr) {
// 	      trigMatchDeltaPt2=trigMatchDeltaPt;
// 	      match2=true;
// 	      m2=n;
// 	    } // if delta
// 	  } // if pt match
// 	} // if direction match
//       } // end for loop n, tracks in a trig object

//       if (localPrint)
// 	cout << " match for this module "
// 	     <<m1<<" "<< deltaRmin1 <<" "<<trigMatchDeltaPt1<<" "
// 	     <<m2<<" "<< deltaRmin2 <<" "<<trigMatchDeltaPt2<<endl;

//       ((TH1D*)fHistDir->Get("testhh3"))->Fill(trigMatchDeltaPt1);
//       ((TH1D*)fHistDir->Get("testhh3"))->Fill(trigMatchDeltaPt2);
//       ((TH1D*)fHistDir->Get("testhh1"))->Fill(deltaRmin1);
//       ((TH1D*)fHistDir->Get("testhh1"))->Fill(deltaRmin2);

//       if(match1 || match2) {
// 	  modulesSingleMatched++;
// 	  matchS=true;
//       }

//       // check if this module matched
//       if( (match1 && match2) ) {
// 	if(m1==m2) {
// 	  cout<<"Error:  matched to same particle "<<endl;
// 	} else { // ok
// 	  match=true;
// 	  modulesMatched++;
// 	  if(localPrint) cout<<" matching good for module "<<i<<endl;
// 	  // select the best batch
// 	  double dr = deltaRmin1 + deltaRmin2; // maybe product is better
// 	  if(dr<drMin) { // a better match, save it
// 	    drMin=dr;
// 	    mu1match = m1;
// 	    mu2match = m2;
// 	    hlt1 = pTO->fLabel;
// 	    indx1=i;
// 	    hlt2 = pTO->fLabel;  // redundant
// 	    indx2=i; // redundant
// 	    deltaRminAll1 = deltaRmin1;
// 	    deltaRminAll2 = deltaRmin2;
// 	    trigMatchDeltaPtAll1 = trigMatchDeltaPt1;
// 	    trigMatchDeltaPtAll2 = trigMatchDeltaPt2;
// 	  }
// 	}
//       }  // if match1&&match2

//     } // end if a valid trigger module, i

//   } // loop over all modules

//   ((TH1D*)fHistDir->Get("testhh4"))->Fill(trigMatchDeltaPtAll1);
//   ((TH1D*)fHistDir->Get("testhh4"))->Fill(trigMatchDeltaPtAll2);
//   ((TH1D*)fHistDir->Get("testhh2"))->Fill(deltaRminAll1);
//   ((TH1D*)fHistDir->Get("testhh2"))->Fill(deltaRminAll2);

//   bool veto = false;
//   if(singleMatch) { // check single matched only
//     if( matchS && ((modulesSelected-modulesSingleMatched)<=0) )
//       {veto=true; if(localPrint) cout<<" single veto "<<endl;}
//   } else { // singleMatch=false, use double
//     if( match && ((modulesSelected-modulesMatched)<=0) )
//       {veto=true; if(localPrint) cout<<" double veto "<<endl;}
//   }

//   if (localPrint) {
//     cout<<" veto = "<<veto<<" "<<match<<" "<<matchS<<endl;
//     cout<<" modules "<<modulesSelected<<" "<<modulesMatched<<" "<<modulesSingleMatched<<endl;
//     cout << " best match "
// 	 <<indx1<<" "<< deltaRminAll1 << " "<<mu1match<<" "<<hlt1<<" "<<trigMatchDeltaPtAll1<<" "
// 	 <<indx2<<" "<< deltaRminAll2 << " "<<mu2match<<" "<<hlt2<<" "<<trigMatchDeltaPtAll2<<endl;
//   }


//   return veto;
// }

// ----------------------------------------------------------------------
void candAnaHh::hhAnalysis() {

  int count = 0;

  TAnaCand *pCand(0);
  TAnaTrack *pPi1, *pPi2;
  double fl3d(0), fls3d(0.), flsxy(0.), prob(0.), chi2(0.), alpha(0.), dr(0.);
  int tm=0;
  int ncand(0);

  if(fVerbose>10)
    cout << "Evt: " << fEvt << " ----------------------------------------------------------------------" << endl;


  // -- loop over all seq vtx fit candidates for D*
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {  // just count candidates
    pCand = fpEvt->getCand(iC);
    if (211211 == pCand->fType || pCand->fType == 91 || pCand->fType == -91) {

      if(pCand->fType == 211211) ++ncand;
      if (fIsMC>0 && fVerbose>0) {
	TVector3 s = pCand->fVtx.fPoint;
	cout << " -> " << iC <<" "<< pCand->fType<<" "<<pCand->fMom;
	cout << " sig: " << pCand->fSig1 << " .. " << pCand->fSig2;
	cout << ": dau " << pCand->fDau1 << " .. " << pCand->fDau2;
	cout << " mass: " << pCand->fMass << " " << s.Z() <<" "<<pCand->fPlab.Mag()<<endl;
	cout << "DUMP HFDHh with mass = " << pCand->fMass << endl;
	dumpHFTruthCand(pCand);
      }
    }
  }


  if(fVerbose>0) cout<<" num of cands "<<ncand<<" "<<fVerbose<<" "<<fIsMC<<endl;
  ((TH1D*)fHistDir->Get("all_cands"))->Fill(ncand);

  ((TH1D*)fHistDir->Get("status"))->Fill(0.);

  // Check MC Gen
  bool ok = false;
  if(fIsMC) {
      ok = anaMC(fpEvt);
      if(fVerbose>10) cout<<" Correct candidate = "<<ok<<endl;
      if(ok) ((TH1D*)fHistDir->Get("status"))->Fill(1.);
  }

  if(ncand>0) {
    ((TH1D*)fHistDir->Get("status"))->Fill(2.);
    ((TH1D*)fHistDir->Get("full_cands"))->Fill(ncand);
  }

  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);

    if (211211 != pCand->fType) continue;  //select HH
    ((TH1D*)fHistDir->Get("status"))->Fill(10.);

    if(fVerbose>0) {cout<<" Dump candidate "<<endl; dumpHFHhCand(pCand); pCand->dump();}

    double candMass = pCand->fMass;
    TVector3 pBdMom  = pCand->fPlab;
    double candPt  = pBdMom.Perp();
    double candEta = pBdMom.Eta();
    double doca = pCand->fMaxDoca; // doca between the 2 pions

    // PV
    int pvidx = (pCand->fPvIdx > -1? pCand->fPvIdx : 0);  // PV index
    TAnaVertex *pv =  fpEvt->getPV(pvidx);
    TVector3 pvPos =  pv->fPoint;
    double pvNtrk = pv->getNtracks();
    double pvNdof = pv->fNdof;
    double pvAveW = ((pvNdof+2.)/2.)/pvNtrk;
    if(fVerbose>10) cout<<" PV "<<pvPos.Z()<<" "<<pvNtrk<<" "<<pvNdof<<" "<<pvAveW<<" "<<pv->fChi2<<" "<<pv->fProb<<endl;

    // HH decay vertex
    TAnaVertex sv = pCand->fVtx;  //  HH vertex
    TVector3 svPos = sv.fPoint;

    TVector3 svpv(svPos-pvPos);

    fl3d = sv.fD3d;
    fls3d = sv.fD3d/sv.fD3dE;
    flsxy = sv.fDxy/sv.fDxyE;
    prob  = sv.fProb;
    chi2  = sv.fChi2;

    //cout<<chi2<<" "<<prob<<" "<<sv.fNdof<<endl;

    alpha = svpv.Angle(pCand->fPlab);  // pointing angle
    //dr = piSlowMom.Angle(pCand->fPlab); // pislow openinig

    // Check that only 2 tracks come from teh candidate
    int indx1 = (pCand->fSig1); // pion1
    int indx2 = (pCand->fSig2); // pion2

    if(fVerbose>10) cout<<" signal tracks "<<(indx2-indx1)<<" "<<indx1<<" "<<indx2<<endl;
    if(indx1<0 || indx2<0 || (indx2-indx1+1)!=2 ) {
      if(fVerbose>0) cout << " Wrong number of candidate isgnal tracks " << indx1<<" "<<indx2 << endl;
      continue;
    }

    ((TH1D*)fHistDir->Get("status"))->Fill(10.);

    // Get pion tracks
    int pi1Id = fpEvt->getSigTrack(indx1)->fIndex;
    int pi2Id = fpEvt->getSigTrack(indx2)->fIndex;
    pPi1 = fpEvt->getRecTrack( pi1Id ); // track 1
    pPi2 = fpEvt->getRecTrack( pi2Id ); // track 2
    //pPi1->dump();
    //pPi2->dump();

    TVector3 pi1Mom = pPi1->fPlab; // pi 1 momentum vector
    TVector3 pi2Mom = pPi2->fPlab; // pi 2 momentum vector
    double pt1 = pi1Mom.Perp();
    double pt2 = pi2Mom.Perp();

    // Get the di-pion opening angle
    dr = pi1Mom.Angle(pi2Mom);

    if(fVerbose>10) {
      cout<<" pion1 "<<pi1Id<<" "<< pt1<<" "<<pPi1->fQ<<" "<<fpEvt->getSigTrack(indx1)->fMCID <<endl;
      cout<<" pion2 "<<pi2Id<<" "<< pt2<<" "<<pPi2->fQ<<" "<<fpEvt->getSigTrack(indx2)->fMCID <<endl;
    }


    double m1    = danekUtils::twoBodyDecayMass(pi1Mom, pi2Mom, MPION, MPION);
    double m2    = danekUtils::twoBodyDecayMass(pi1Mom, pi2Mom, MKAON, MKAON);
    double m3    = danekUtils::twoBodyDecayMass(pi1Mom, pi2Mom, MPION, MKAON);
    double m4    = danekUtils::twoBodyDecayMass(pi1Mom, pi2Mom, MKAON, MPION);
    double m5    = danekUtils::twoBodyDecayMass(pi1Mom, pi2Mom, MPROTON, MPION);
    double m6    = danekUtils::twoBodyDecayMass(pi1Mom, pi2Mom, MPROTON, MKAON);


    // skip candidates with the same charge pions
    //if(pPi1->fQ == pPi2->fQ) {continue;}  // Not needed anymore, done in main selection
    //((TH1D*)fHistDir->Get("status"))->Fill(12.);

    // Oppening angle variables
    double a1 = svpv.Angle(pi1Mom);  // pi1 angle
    double a2 = svpv.Angle(pi2Mom);  // pi1 angle
    double perp1 = pi1Mom.Mag() * sin(a1);
    double perp2 = pi2Mom.Mag() * sin(a2);
    double perpdiff = perp1 - perp2;
    double perpsum = perp1 + perp2;



    // truthMatch return an interger, 0-no match, 1-correct match,
    tm = 0;
    if(fIsMC) {

      tm = truthMatch(pCand,fVerbose); // check truth matching
      if(fVerbose>10) cout<<" Truth matching = "<<tm<<endl;

      if ( (tm == 1) && fVerbose>20 ) {
	  cout << " Truth matched cand -> " << pCand->fType;
	  cout << " sig: " << pCand->fSig1 << " .. " << pCand->fSig2;
	  cout << ": dau " << pCand->fDau1 << " .. " << pCand->fDau2;
	  cout << " mass: " << pCand->fMass << " " << tm<<endl;
	  cout << "DUMP HFDstarCandidate with mass = " << pCand->fMass << endl;
	  dumpHFHhCand(pCand);
	  fpEvt->dumpGenBlock();
      } // end if
    } // if MC

    if(tm==1) {
      ((TH1D*)fHistDir->Get("status"))->Fill(20.);
      if(ok) ((TH1D*)fHistDir->Get("status"))->Fill(30.);
    }


    ((TH1D*)fHistDir->Get("full_m"))->Fill(candMass);

    // Skip muon, Look at muid
    //bool muid1 = goodMuon(pPi1);  // true for good  muons
    //bool muid2 = goodMuon(pPi2);
    //bool muid1 = tightMuon(pPi1);  // true for good/tight  muons
    //bool muid2 = tightMuon(pPi2);
    //if( muid1 || muid2 ) {continue;}

    int mid1 = 0, mid2= 0;
    if (pPi1->fMuIndex > -1) {
      //chiPi= fpEvt->getMuon(pPi->fMuIndex)->fMuonChi2;
      mid1 = fpEvt->getMuon(pPi1->fMuIndex)->fMuID;
    }
    if (pPi2->fMuIndex > -1)  {
      //chiK = fpEvt->getMuon(pK->fMuIndex)->fMuonChi2;
      mid2 = fpEvt->getMuon(pPi2->fMuIndex)->fMuID;
    }
    if(mid1!=0) ((TH1D*)fHistDir->Get("h40"))->Fill(mid1);
    if(mid2!=0) ((TH1D*)fHistDir->Get("h40"))->Fill(mid2);

    // Skip, there was a muon
    //if( mid1!=0 || mid2!=0 ) cout<<" muon "<<mid1<<" "<<pPi1->fMuID<<" "<<mid2<<" "<<pPi2->fMuID
    //			 <<" "<<pPi1->fMuIndex<<" "<<pPi2->fMuIndex<< endl;
    if( MUON_VETO==1 && (mid1!=0 || mid2!=0) ) {continue;}


    ((TH1D*)fHistDir->Get("status"))->Fill(11.);
    if(tm==1) ((TH1D*)fHistDir->Get("status"))->Fill(21.);

    bool t1 = highPurity(pPi1);
    bool t2 = highPurity(pPi2);
    //cout<<t1<<" "<<t2<<endl;
    if( !t1 || !t2) continue; // skip bad tracks

    ((TH1D*)fHistDir->Get("status"))->Fill(12.);
    if(tm==1) ((TH1D*)fHistDir->Get("status"))->Fill(22.);

    // Histogram
    ((TH1D*)fHistDir->Get("all_fl3d"))->Fill(fl3d);
    ((TH1D*)fHistDir->Get("all_fls3d"))->Fill(fls3d);
    ((TH1D*)fHistDir->Get("all_flsxy"))->Fill(flsxy);
    ((TH1D*)fHistDir->Get("all_prob"))->Fill(prob);
    ((TH1D*)fHistDir->Get("all_chi2"))->Fill(chi2);
    ((TH1D*)fHistDir->Get("all_alpha"))->Fill(alpha);
    ((TH1D*)fHistDir->Get("all_pt"))->Fill(candPt);
    ((TH1D*)fHistDir->Get("all_m"))->Fill(candMass);
    ((TH1D*)fHistDir->Get("all_eta"))->Fill(candEta);
    ((TH1D*)fHistDir->Get("all_dr"))->Fill(dr);
    ((TH1D*)fHistDir->Get("all_ptPi"))->Fill(pt1);
    ((TH1D*)fHistDir->Get("all_ptPi"))->Fill(pt2);
    ((TH1D*)fHistDir->Get("all_doca"))->Fill(doca);
    ((TH1D*)fHistDir->Get("all_pvW"))->Fill(pvAveW);
    ((TH1D*)fHistDir->Get("all_pvid"))->Fill(float(pvidx));
    ((TH2D*)fHistDir->Get("h2d"))->Fill(candPt,candMass);



    if(ok && tm==1) {  // Do here comparison between RECO and GEN quantities
      double c1 = (PVGen-pvPos).Mag();     // PV distance
      double c2 = (BdVertexGen-svPos).Mag(); // SV distance

      //cout<<" RECO-MV vertex distance "<<c1<<" "<<c2<<" "<<c3<<endl;
      ((TH1D*)fHistDir->Get("h31"))->Fill(c1);
      ((TH1D*)fHistDir->Get("h32"))->Fill(c2);

      double a11 = Pi1MomGen.Angle(pi1Mom);  // pi1 direction
      double a12 = Pi2MomGen.Angle(pi2Mom);  // pi2 direction
      double a13 = BdMomGen.Angle(pBdMom);   //  Bd direction

      ((TH1D*)fHistDir->Get("h35"))->Fill(a11);
      ((TH1D*)fHistDir->Get("h36"))->Fill(a12);
      ((TH1D*)fHistDir->Get("h37"))->Fill(a13);

    } // if OK

    // Do the selection cuts
    if (doca > 0.025) continue;
    ((TH1D*)fHistDir->Get("status"))->Fill(13.);
    if(tm==1) ((TH1D*)fHistDir->Get("status"))->Fill(23.);

    if (pvAveW < 0.70) continue;
    ((TH1D*)fHistDir->Get("status"))->Fill(14.);
    if(tm==1) ((TH1D*)fHistDir->Get("status"))->Fill(24.);

    if (fls3d < 1) continue;
    ((TH1D*)fHistDir->Get("status"))->Fill(15.);
    if(tm==1) ((TH1D*)fHistDir->Get("status"))->Fill(25.);

    if (candPt < 6.0) continue;
    ((TH1D*)fHistDir->Get("status"))->Fill(16.);
    if(tm==1) ((TH1D*)fHistDir->Get("status"))->Fill(26.);

    if (chi2 > 3.0) continue;
    ((TH1D*)fHistDir->Get("status"))->Fill(17.);
    if(tm==1) ((TH1D*)fHistDir->Get("status"))->Fill(27.);

    if (dr > 1.3) continue;
    ((TH1D*)fHistDir->Get("status"))->Fill(18.);
    if(tm==1) ((TH1D*)fHistDir->Get("status"))->Fill(28.);

    if (alpha > 0.7 && alpha < 2.4) continue;
    ((TH1D*)fHistDir->Get("status"))->Fill(19.);
    if(tm==1) ((TH1D*)fHistDir->Get("status"))->Fill(29.);

    //if(pt1<4. || pt2<4.) continue;  //  does not do aything for data, cand reco is already with 4

    if(candMass<HH_MLO || candMass>HH_MHI) continue;
    ((TH1D*)fHistDir->Get("status"))->Fill(20.);
    if(tm==1) ((TH1D*)fHistDir->Get("status"))->Fill(30.);

    pair<int, int> pclose;
    //                       dcaCut(cm) dcas ptCut(GeV)
    pclose = nCloseTracks(pCand,0.03, 2, 0.5); // around Bd
    int close = pclose.first;
    //                                      dca   R    Pt
    double iso = isoClassicWithDOCA(pCand, 0.05,0.7, 0.9); // arount Bd


    ((TH1D*)fHistDir->Get("h38"))->Fill(iso);
    ((TH1D*)fHistDir->Get("h39"))->Fill(close);


    count++; // count selected candidates

    // Final histos after cuts

    // Histogram
    ((TH1D*)fHistDir->Get("fl3d"))->Fill(fl3d);
    ((TH1D*)fHistDir->Get("fls3d"))->Fill(fls3d);
    ((TH1D*)fHistDir->Get("flsxy"))->Fill(flsxy);
    ((TH1D*)fHistDir->Get("prob"))->Fill(prob);
    ((TH1D*)fHistDir->Get("chi2"))->Fill(chi2);
    ((TH1D*)fHistDir->Get("alpha"))->Fill(alpha);
    ((TH1D*)fHistDir->Get("pt"))->Fill(candPt);
    ((TH1D*)fHistDir->Get("m"))->Fill(candMass);
    ((TH1D*)fHistDir->Get("eta"))->Fill(candEta);
    ((TH1D*)fHistDir->Get("dr"))->Fill(dr);
    ((TH1D*)fHistDir->Get("ptPi"))->Fill(pt1);
    ((TH1D*)fHistDir->Get("ptPi"))->Fill(pt2);
    ((TH1D*)fHistDir->Get("doca"))->Fill(doca);
    ((TH1D*)fHistDir->Get("pvW"))->Fill(pvAveW);
    ((TH1D*)fHistDir->Get("pvid"))->Fill(float(pvidx));

    ((TH1D*)fHistDir->Get("m1"))->Fill(m1);
    ((TH1D*)fHistDir->Get("m2"))->Fill(m2);
    ((TH1D*)fHistDir->Get("m3"))->Fill(m3);
    ((TH1D*)fHistDir->Get("m4"))->Fill(m4);


    if ( !fIsMC || tm ==1 ) {  // histogram truth matched candidates
      // Histogram
      ((TH1D*)fHistDir->Get("mc_fl3d"))->Fill(fl3d);
      ((TH1D*)fHistDir->Get("mc_fls3d"))->Fill(fls3d);
      ((TH1D*)fHistDir->Get("mc_flsxy"))->Fill(flsxy);
      ((TH1D*)fHistDir->Get("mc_prob"))->Fill(prob);
      ((TH1D*)fHistDir->Get("mc_chi2"))->Fill(chi2);
      ((TH1D*)fHistDir->Get("mc_alpha"))->Fill(alpha);
      ((TH1D*)fHistDir->Get("mc_pt"))->Fill(candPt);
      ((TH1D*)fHistDir->Get("mc_m"))->Fill(candMass);
      ((TH1D*)fHistDir->Get("mc_eta"))->Fill(candEta);
      ((TH1D*)fHistDir->Get("mc_dr"))->Fill(dr);
      ((TH1D*)fHistDir->Get("mc_ptPi"))->Fill(pt1);
      ((TH1D*)fHistDir->Get("mc_ptPi"))->Fill(pt2);
      ((TH1D*)fHistDir->Get("mc_doca"))->Fill(doca);
      ((TH1D*)fHistDir->Get("mc_pvW"))->Fill(pvAveW);
      ((TH1D*)fHistDir->Get("mc_pvid"))->Fill(float(pvidx));


      ((TH1D*)fHistDir->Get("h22"))->Fill(a1);
      ((TH1D*)fHistDir->Get("h22"))->Fill(a2);
      ((TH1D*)fHistDir->Get("h23"))->Fill(perp1);
      ((TH1D*)fHistDir->Get("h23"))->Fill(perp2);
      ((TH1D*)fHistDir->Get("h24"))->Fill(perpdiff);
      ((TH1D*)fHistDir->Get("h25"))->Fill(perpsum);
      ((TH2D*)fHistDir->Get("h26"))->Fill(a1,a2);
      ((TH2D*)fHistDir->Get("h27"))->Fill(perp1,perp2);

    }  // if

    // Save in a tree, save only masses between 130-160MeV
    if( (fcands<10) ) {
      ftm[fcands] = tm;

      fm[fcands]=candMass;
      fpt[fcands]=candPt;
      fchi2[fcands]=chi2;
      falpha[fcands]=alpha;
      ffls3d[fcands]=fls3d;
      fdr[fcands]=dr;
      fptpi1[fcands]=pt1;
      fptpi2[fcands]=pt2;
      fdoca[fcands]=doca;
      fweight[fcands]=pvAveW;

      fclose[fcands]=close;
      fiso[fcands]=iso;
      fperp1[fcands]=perp1;
      fperp2[fcands]=perp2;

      fm1[fcands]=m1;
      fm2[fcands]=m2;
      //fm3[fcands]=m3;
      //fm4[fcands]=m4;
      fm3[fcands]=m5;
      fm4[fcands]=m6;

      fcands++;
    } // if fcands

  }  // candidate loop

  if(count>0) {
    ((TH1D*)fHistDir->Get("cands"))->Fill(float(count));
    ((TH1D*)fHistDir->Get("status"))->Fill(3.);
    if(fcands>0) tree->Fill();
  }



}

// ----------------------------------------------------------------------
void candAnaHh::dumpHFTruthCand(TAnaCand *pC) {
  TAnaTrack *pT(0);
  if( pC->fSig1 == -1 && pC->fSig2==-1 ) return;
  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getRecTrack(fpEvt->getSigTrack(i)->fIndex);
    pT->dump();
  }
}


// ----------------------------------------------------------------------
void candAnaHh::dumpHFHhCand(TAnaCand *pC) {
  TAnaTrack *pT(0);

  cout << "HFHhCand: idx = " << pC->fIndex << " type = " << pC->fType<< " m = " << pC->fMass <<endl;
  cout<<"DOCA "<<pC->fMinDoca<<" "<<pC->fMaxDoca<<endl;


//   int nsize = pC->fNstTracks.size();
//   cout<<" tracks "<<nsize<<endl;
//   if (nsize>0) {
//     for(int i = 0; i<nsize; ++i) {
//       int trkId = pC->fNstTracks[i].first;
//       double doca = pC->fNstTracks[i].second.first;
//       cout<<i<<" "<<trkId<<" "<<doca<<endl;
//     }
//   }

  // -- D0 daughters
  if ( pC->fSig1<0  || pC->fSig2<0 ) {
    cout << "XXXXXXXXX cannot get signal cand of " << pC->fType << endl;
    return;
  }

  // -- pion 1
  int indx = (pC->fSig1);
  pT = fpEvt->getRecTrack( fpEvt->getSigTrack(indx)->fIndex );
  cout << fpEvt->getSigTrack(indx)->fMCID << " " ;
  pT->dump();

  // -- pion 2
  indx = (pC->fSig2);
  pT = fpEvt->getRecTrack( fpEvt->getSigTrack(indx)->fIndex );
  cout << fpEvt->getSigTrack(indx)->fMCID << " " ;
  pT->dump();



}


// ----------------------------------------------------------------------
// -- works ONLY for this dedicated decay mode with the seq vtx fitter!
int candAnaHh::truthMatch(TAnaCand *pCand, int verbose) {

  // -- check pion1
  TAnaTrack *pT = fpEvt->getRecTrack(fpEvt->getSigTrack(pCand->fSig1)->fIndex);
  if(fVerbose>0) cout<<fpEvt->getSigTrack(pCand->fSig1)->fIndex<<" "
		<<pT->fGenIndex<<" "
		<<(fpEvt->getRecTrack(pT->fIndex)->fGenIndex)<<" "<<verbose<<endl;  // same as above

  if (pT->fGenIndex < 0) {
    if (verbose > 0) cout << "pT->fGenIndex < 0" << endl;
    return 0;
  }
  TGenCand  *pG = fpEvt->getGenCand(fpEvt->getRecTrack(pT->fIndex)->fGenIndex);
  if (0 == pG) {
    if (verbose > 0) cout << "0 == pG" << endl;
    return 0;
  }
  if(fVerbose>0) cout<< pG->fID<<" "<<endl;  // gen id
  if (211 != TMath::Abs(pG->fID)) {
    if (verbose > 0) cout << "211 != TMath::Abs(pG->fID)" << endl;
    return 0;
  }

  // Check the mother
  int mom = pG->fMom1;
  pG = fpEvt->getGenCand(pG->fMom1);
  if ((0 == pG) || 511 != TMath::Abs(pG->fID)) {
    if (verbose > 0) cout << "(0 == pG) || 511 != pG->fID, pG->fID = " << pG->fID  << endl;
    return 0;
  }


  // -- check pion2
  pT = fpEvt->getRecTrack(fpEvt->getSigTrack(pCand->fSig2)->fIndex);
  if(fVerbose>0) cout<<fpEvt->getSigTrack(pCand->fSig2)->fIndex<<" "
		     <<pT->fGenIndex<<" "
		     <<(fpEvt->getRecTrack(pT->fIndex)->fGenIndex)<<" "<<verbose<<endl;  // same as above

  if (pT->fGenIndex < 0) {
    if (verbose > 0) cout << "pT->fGenIndex < 0" << endl;
    return 0;
  }

  pG = fpEvt->getGenCand(fpEvt->getRecTrack(pT->fIndex)->fGenIndex);
  if (0 == pG) {
    if (verbose > 0) cout << "0 == pG" << endl;
    return 0;
  }
  if(fVerbose>0) cout<< pG->fID<<" "<<endl;  // gen id
  if (211 != TMath::Abs(pG->fID)) {
    if (verbose > 0) cout << "211 != TMath::Abs(pG->fID)" << endl;
    return 0;
  }

  // Check the mother
  mom = pG->fMom1;
  pG = fpEvt->getGenCand(pG->fMom1);
  if ((0 == pG) || 511 != TMath::Abs(pG->fID)) {
    if (verbose > 0) cout << "(0 == pG) || 511 != pG->fID, pG->fID = " << pG->fID  << " "<<mom<<endl;
    return 0;
  }


  // -- check daughters
  int type(0), moIdx(-1);
  int count=0, countPi=0; // , countK=0, missK=0, missPi=0;
  //int daus = (pC->fSig2) - (pC->fSig1);
  //cout<<" Bd daugthers "<<(daus+1)<<endl;
  for (int id = pCand->fSig1; id <= pCand->fSig2; ++id) {
    pT = fpEvt->getSigTrack(id);
    type = pT->fMCID;
    pT = fpEvt->getRecTrack(pT->fIndex);
    if (pT->fGenIndex < 0) {
      if (verbose > 0) cout << "no pT->fGenIndex" << endl;
      return 0;
    }
    pG = fpEvt->getGenCand(pT->fGenIndex);
    if (moIdx < 0) {
      moIdx = pG->fMom1;
    } else {
      if (moIdx != pG->fMom1) {
	if (verbose > 0) cout << "moIdx != pG->fMom1" << endl;
	return 0;
      }
    }
    if (verbose > 0)
      cout << "dau cand sigtrack " << id
 	   << " with type = " << type
 	   << " and gen ID = " << pG->fID
 	   << " at gen idx = " << pT->fGenIndex
 	   << endl;


    //if (TMath::Abs(type) !=  TMath::Abs(pG->fID)) {
    //if (verbose > 0) cout << "TMath::Abs(type) != TMath::Abs(pG->fID), type = " << type << " pG->fID = " << pG->fID
    //			    << " track " << pT->fIndex << endl;
    //  return 0;
    // }

    count++;

    if( TMath::Abs(pG->fID) == 211 ) { // Pion
      countPi++;
      //cout <<count<<" "<<countPi;
      // The check below does not work because we assing the muon mass to all tracks
      //if( TMath::Abs(type) != 211 ) {
      //cout << " Pion identified as ?, type = " << type << " pG->fID = " << pG->fID
      //     << " track " << pT->fIndex;
      //missPi++;
      //return 0;
      //}
      //cout<<endl;
    }

  }

  if     (countPi==2) return  1; // select righ combination

  return 0;
}


// ----------------------------------------------------------------------
//void candAnaHh::moreBasicCuts() {
//cout << "   candAnaHh: more basic cuts" << endl;
//}


// ----------------------------------------------------------------------
void candAnaHh::bookHist() {
  cout << "==>candAnaHh: bookHist" << endl;
  TH1 *h=NULL;

  candAna::bookHist();

  fHistDir->cd();
  h = new TH1D("testhh0", "stat", 100, -1., 99.);

  //h = new TH1D("testhh1", "dr", 1000, 0., 1.);
  //h = new TH1D("testhh2", "dr", 1000, 0., 1.);
  //h = new TH1D("testhh3", "dpt", 200, -1., 1.);
  //h = new TH1D("testhh4", "dpt", 200, -1., 1.);

  h = new TH1D("testhh11", "dr", 1000, 0, 1);
  h = new TH1D("testhh12", "dr", 1000, 0, 1);
  h = new TH1D("testhh13", "dr", 1000, 0, 1);
  h = new TH1D("testhh14", "dr", 1000, 0, 1);

  //return;
  //#endif
  //fHistDir->cd();

#ifdef MYCODE

  h = new TH1D("status", "status", 100, -0.5, 99.5);
//   h = new TH1D("mdz", "m(d0)", 70, 1.8, 2.5);
//   h = new TH1D("dm", "delta(m)", 60, 0.13, 0.16);
//   h = new TH1D("ncand", "ncand", 200, 0., 200);
//   h = new TH1D("fls3d", "fls3d", 60, 0., 20);
//   h = new TH1D("flsxy", "flsxy", 60, 0., 20);
//   h = new TH1D("prob", "prob(chi2/dof)", 100, 0., 1.);
//   h = new TH1D("chi2", "chi2", 100, 0., 10.);
//   h = new TH1D("alpha", "alpha", 50, 0., 1.0);
//   h = new TH1D("dr","dr",100, 0., 1);
//   h = new TH1D("pt", "pT", 50, 0., 25);
//   h = new TH1D("ptdz", "pT", 50, 0., 25);
//   h = new TH1D("ptK",   "pT", 50, 0., 10);
//   h = new TH1D("ptPi",  "pT", 50, 0., 10);
//   h = new TH1D("ptPis", "pT", 50, 0., 5);
  TH2D *h2 = new TH2D("h2d", "m vs pt", 50, 0., 25, 120, 0., 12.);

  int mass_size = 100;  // 100
  float mass_min = HH_MLO; // 4.;
  float mass_max = HH_MHI; //7.;
  h = new TH1D("full_m", "cand mass", mass_size,mass_min, mass_max);
  h = new TH1D("all_m", "cand mass", mass_size,mass_min, mass_max);
  h = new TH1D("all_cands", "cands", 200, 0., 200);
  h = new TH1D("full_cands", "cands", 200, 0., 200);
  h = new TH1D("all_fl3d", "fl3d", 50, 0., 5.);
  h = new TH1D("all_fls3d", "fls3d", 60, 0., 20);
  h = new TH1D("all_flsxy", "flsxy", 60, 0., 20);
  h = new TH1D("all_prob", "prob(chi2/dof)", 100, 0., 1.);
  h = new TH1D("all_chi2", "chi2", 100, 0., 10.);
  h = new TH1D("all_alpha", "alpha",175, 0., 3.5);
  h = new TH1D("all_dr","dr",100, 0., 3.);
  h = new TH1D("all_pt", "cand pT", 100, 0., 50.);
  h = new TH1D("all_ptPi","pions pt", 100, 0., 20);
  h = new TH1D("all_doca","pions doca", 100, 0., 0.2);
  h = new TH1D("all_pvW", "PV Weight", 100, 0., 2);
  h = new TH1D("all_eta", "cand eta", 60, -3.0, 3.0);
  h = new TH1D("all_pvid", "cand PVidx", 100, 0., 100);

  h = new TH1D("mc_m", "cand mass", mass_size, mass_min, mass_max);
  h = new TH1D("mc_ncand", "ncand", 200, 0., 200);
  h = new TH1D("mc_fl3d", "fl3d", 50, 0., 5.);
  h = new TH1D("mc_fls3d", "fls3d", 60, 0., 20);
  h = new TH1D("mc_flsxy", "flsxy", 60, 0., 20);
  h = new TH1D("mc_prob", "prob(chi2/dof)", 100, 0., 1.);
  h = new TH1D("mc_chi2", "chi2", 100, 0., 10.);
  h = new TH1D("mc_alpha", "alpha", 175, 0.,3.5);
  h = new TH1D("mc_dr","dr",100, 0., 3.);
  h = new TH1D("mc_pt", "cand pT", 100, 0., 50.);
  h = new TH1D("mc_ptPi",  "pions pt", 100, 0., 20);
  h = new TH1D("mc_doca","pions doca", 100, 0., 0.2);
  h = new TH1D("mc_pvW", "PV Weight", 100, 0., 2);
  h = new TH1D("mc_eta", "cand eta", 60, -3.0, 3.0);
  h = new TH1D("mc_pvid", "cand PVidx", 100, 0., 100);

  h = new TH1D("m1", "cand mass", mass_size, mass_min, mass_max);
  h = new TH1D("m2", "cand mass", mass_size, mass_min, mass_max);
  h = new TH1D("m3", "cand mass", mass_size, mass_min, mass_max);
  h = new TH1D("m4", "cand mass", mass_size, mass_min, mass_max);

  h = new TH1D("m", "cand mass", mass_size, mass_min, mass_max);
  h = new TH1D("ncand", "ncand", 200, 0., 200);
  h = new TH1D("fl3d", "fl3d", 50, 0., 5.);
  h = new TH1D("fls3d", "fls3d", 60, 0., 20);
  h = new TH1D("flsxy", "flsxy", 60, 0., 20);
  h = new TH1D("prob", "prob(chi2/dof)", 100, 0., 1.);
  h = new TH1D("chi2", "chi2", 100, 0., 10.);
  h = new TH1D("alpha", "alpha", 375, 0., 3.5);
  h = new TH1D("dr","dr",100, 0., 3.);
  h = new TH1D("pt", "cand pT", 100, 0., 50.);
  h = new TH1D("ptPi",  "pions pt", 100, 0., 20);
  h = new TH1D("doca","pions doca", 100, 0., 0.2);
  h = new TH1D("pvW", "PV Weight", 100, 0., 2);
  h = new TH1D("eta", "cand eta", 60, -3.0, 3.0);
  h = new TH1D("pvid", "cand PVidx", 100, 0., 100);


  //h2 = new TH2D("all_h2d", "m(d0) vs dm", 60, 1.8, 1.92, 60, 0.13, 0.16);

  h = new TH1D("cands", "cands", 100, 0., 100);

//   h = new TH1D("h1", "h1", 100, 0., 1);
//   h = new TH1D("h2", "h2", 100, 0., 1);
//   h = new TH1D("h3", "h3", 100, 0., 1);
//   h = new TH1D("h4", "h4", 350, 0., 3.5);
//   h = new TH1D("h5", "h5", 100, 0., 1);
//   h = new TH1D("h6", "h6", 100, 0., 1);
//   h = new TH1D("h7", "h7", 350, 0., 3.5);
//   h = new TH1D("h8", "h8", 350, 0., 3.5);
//   h = new TH1D("h9", "h9", 350, 0., 3.5);
//   h = new TH1D("h10","h10",100, 0., 1);

  h = new TH1D("h11", "h11", 100, 0., 1);
  h = new TH1D("h12", "h12", 100, 0., 1.);
  h = new TH1D("h13", "h13", 100, 0., 4.);
  h = new TH1D("h14", "h14", 100,-1., 1);
  h = new TH1D("h15", "h15", 100, 0., 10);
  h2 = new TH2D("h16", "h16", 20, 0., 1. ,20,0.,1.);
  h2 = new TH2D("h17", "h17", 40, 0., 4.,40,0.,4.);
//   h = new TH1D("h18", "h18", 350, 0., 3.5);
//   h = new TH1D("h19", "h19", 350, 0., 3.5);
  h = new TH1D("h20", "h20", 100, 0.,0.2);
//   h = new TH1D("h21", "h21", 350, 0., 3.5);
  h = new TH1D("h22", "h22", 100, 0., 1.);
  h = new TH1D("h23", "h23", 100, 0., 5.);
  h = new TH1D("h24", "h24", 100,-2., 2.);
  h = new TH1D("h25", "h25", 100, 0., 10);
  h2 = new TH2D("h26", "h26", 20, 0., 1. ,20,0.,1.);
  h2 = new TH2D("h27", "h27", 40, 0., 5.,40,0.,5.);

  h = new TH1D("h28", "h28", 80, 0., 40);
  h = new TH1D("h29", "h29", 80, 0., 40);
  h = new TH1D("h30", "h30", 80, 0., 40);

  h = new TH1D("h31", "h31", 100, 0., 1);
  h = new TH1D("h32", "h32", 100, 0., 1);
//   h = new TH1D("h33", "h33", 100, 0., 1);
//   h = new TH1D("h34", "h34", 350, 0., 0.35);
  h = new TH1D("h35", "h35", 350, 0., 0.35);
  h = new TH1D("h36", "h36", 350, 0., 0.35);
  h = new TH1D("h37", "h37", 350, 0., 0.35);
  h = new TH1D("h38", "h38", 120, 0., 1.2);
  h = new TH1D("h39", "h39", 40, 0., 40);
  h = new TH1D("h40", "h40", 100, 0.,1000000);

  h = new TH1D("h41", "h41",80, 0.0, 40.);
  h = new TH1D("h42", "h42",100, 0.0, 20.);
//   h = new TH1D("h43", "h43",50, 0.0, 20.);
  h = new TH1D("h44", "h44",100, 0.0, 20.);

//   h = new TH1D("h45", "h45",100,  1.,2.5);
//   h = new TH1D("h46", "h46",100,  1.,2.5);
  h = new TH1D("h47", "h47",100,  0., 3.0);
//   h = new TH1D("h48", "h48",50,  0., 20.);
//   h = new TH1D("h49", "h49",200,  0.13, 0.18);
//   h = new TH1D("h50", "h50",200,  0.13, 0.18);

//   h2 = new TH2D("h51", "h51", 40, 0., 20., 80, 1.5, 2.3);
//   h2 = new TH2D("h52", "h52", 35, 0., 3.5, 80, 1.5, 2.3);
//   h2 = new TH2D("h53", "h53", 40, 0., 20., 80, 1.5, 2.3);
//   h2 = new TH2D("h54", "h53", 40, 0., 10., 80, 1.5, 2.3);
//   h2 = new TH2D("h55", "h55", 40, 0., 10., 80, 1.5, 2.3);
//   h2 = new TH2D("h56", "h56", 40, 0.,  2., 80, 1.5, 2.3);
//   h2 = new TH2D("h57", "h57", 40, 0.,0.4,  80, 1.5, 2.3);

//   h = new TH1D("h60", "h60", 5,-2.5,2.5);
//   h = new TH1D("h61", "h61", 5,-2.5,2.5);
//   h = new TH1D("h62", "h62", 5,-2.5,2.5);
//   //h = new TH1D("h63", "h63", 5,-2.5,2.5);
//   //h = new TH1D("h64", "h64", 5,-2.5,2.5);
//   //h = new TH1D("h65", "h65", 5,-2.5,2.5);
//   //h = new TH1D("h66", "h66", 5,-2.5,2.5);
//   //h = new TH1D("h67", "h67", 5,-2.5,2.5);

//   h = new TH1D("h71", "dm",40,0.135,0.155);
//   h = new TH1D("h72", "dm",40,0.135,0.155);
//   h = new TH1D("h73", "dm",40,0.135,0.155);
//   h = new TH1D("h74", "dm",40,0.135,0.155);
//   h = new TH1D("h75", "dm",40,0.135,0.155);
//   h = new TH1D("h76", "dm",40,0.135,0.155);
//   h = new TH1D("h77", "dm",40,0.135,0.155);

  // tree = new TTree("hh","hh");
  // tree->Branch("fcands",&fcands,"fcands/I");
  // tree->Branch("ftm",ftm,"ftm[fcands]/I");
  // tree->Branch("fm",fm,"fm[fcands]/F");
  // tree->Branch("ffls3d",ffls3d,"ffls3d[fcands]/F");
  // tree->Branch("fchi2",fchi2,"fchi2[fcands]/F");
  // tree->Branch("falpha",falpha,"falpha[fcands]/F");
  // tree->Branch("fdr",fdr,"fdr[fcands]/F");
  // tree->Branch("fpt",fpt,"fpt[fcands]/F");
  // tree->Branch("fptpi1",fptpi1,"fptpi1[fcands]/F");
  // tree->Branch("fptpi2",fptpi2,"fptpi2[fcands]/F");
  // tree->Branch("fdoca",fdoca,"fdoca[fcands]/F");
  // tree->Branch("fweight",fweight,"fweight[fcands]/F");
  // tree->Branch("fclose",fclose,"fclose[fcands]/I");
  // tree->Branch("fiso",fiso,"fiso[fcands]/F");
  // tree->Branch("fperp1",fperp1,"fperp1[fcands]/F");
  // tree->Branch("fperp2",fperp2,"fperp2[fcands]/F");

  // tree->Branch("fm1",fm1,"fm1[fcands]/F");
  // tree->Branch("fm2",fm2,"fm2[fcands]/F");
  // tree->Branch("fm3",fm3,"fm3[fcands]/F");
  // tree->Branch("fm4",fm4,"fm4[fcands]/F");

#endif

}

// ----------------------------------------------------------------------
void candAnaHh::readCuts(string filename, int dump) {
  candAna::readCuts(filename, dump);

  fCutFile = filename;

  if (dump) cout << "==> candAnaHh: Reading " << fCutFile << " for cut settings" << endl;
  vector<string> cutLines;
  readFile(fCutFile, cutLines);

  char CutName[100];
  float CutValue;

  char  buffer[200];
  fHistDir->cd();
  TH1D *hcuts = (TH1D*)fHistDir->Get("hcuts");
  hcuts->GetXaxis()->SetBinLabel(200, fCutFile.c_str());
  string cstring = "B cand";
  int ibin;

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str());

    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "MUON_VETO")) {
      MUON_VETO = int(CutValue);
      if (dump) cout << "MUON_VETO:      " << MUON_VETO << endl;
      ibin = 211;
      hcuts->SetBinContent(ibin, MUON_VETO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: muon-veto :: %d", CutName, MUON_VETO));
    }

    if (!strcmp(CutName, "HH_MLO")) {
      HH_MLO = CutValue;
      if (dump) cout << "HH_MLO:      " << HH_MLO << endl;
      ibin = 212;
      hcuts->SetBinContent(ibin, HH_MLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: hh mass-low :: %3.1f", CutName, HH_MLO));
    }

    if (!strcmp(CutName, "HH_MHI")) {
      HH_MHI = CutValue;
      if (dump) cout << "HH_MHI:      " << HH_MHI << endl;
      ibin = 213;
      hcuts->SetBinContent(ibin, HH_MHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: hh mass-high :: %3.1f", CutName, HH_MHI));
    }
  } // end for

}


// ----------------------------------------------------------------------
void candAnaHh::genMatch() {
  cout << "genMatch()  function" << endl;
}

// ----------------------------------------------------------------------
void candAnaHh::genMatchOld() {
  cout << "genMatchOld()  function" << endl;
}


// ----------------------------------------------------------------------
void candAnaHh::recoMatch() {
  cout << "recoMatch() function" << endl;
}


// ----------------------------------------------------------------------
void candAnaHh::candMatch() {
  cout << "candMatch()  function" << endl;
}


// ----------------------------------------------------------------------
void candAnaHh::efficiencyCalculation() {
  cout << "efficiencyCalculation() function" << endl;
}

// ----------------------------------------------------------------------
void candAnaHh::processType() {
  cout << "processType() function" << endl;
}
