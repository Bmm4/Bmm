#ifndef PRESELECTION
#define PRESELECTION

#include <string>
#include "TMath.h"
#include "TH1.h"
#include "redTreeData.hh"

// ----------------------------------------------------------------------
std::string preselection();
TH1D* getPreselectionNumbers();

// ----------------------------------------------------------------------
bool preselection(redTreeData &b, int channel);
void printRedTreeEvt(redTreeData &b);

#endif

