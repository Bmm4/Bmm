#ifndef PRESELECTION
#define PRESELECTION

#include <string>
#include <map>
#include "TMath.h"
#include "TH1.h"
#include "redTreeData.hh"

// ----------------------------------------------------------------------
std::string preselection();
TH1D* getPreselectionNumbers();

// ----------------------------------------------------------------------
bool preselection(redTreeData &b);
void printRedTreeEvt(redTreeData &b);


// ----------------------------------------------------------------------
class presel {
 public:
  presel();
  presel(TH1*);
  void setupMap();

  bool preselection(redTreeData &b, int verbose = -1);
  bool preselection1(redTreeData &b);
  std::string preselection();
  void setCut(std::string name, double value);
  bool passCut(std::string name, double value);
  TH1D* getPreselectionNumbers();

  double fPTMIN,      fPTMAX;
  double fM1PTMIN,    fM1PTMAX;
  double fM2PTMIN,    fM2PTMAX;

  double fFL3DMAX;
  double fCHI2DOFMAX;
  double fPVIPMAX,    fPVIPSMAX;
  double fPVLIPMAX,   fPVLIPSMAX;
  double fMAXDOCAMAX;
  int    fCLOSETRKMAX;
  double fFLSXYMIN;

  double fFLS3DMIN,   fFLS3DMAX;
  double fDOCATRKMAX;
  double fISOMIN;
  double fALPHAMAX;
  double fMASSERRORMAX;

  std::map<std::string, double> fCuts; // I think this is quite useless

};

#endif
