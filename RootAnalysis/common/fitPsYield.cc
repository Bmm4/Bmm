#include "fitPsYield.hh"
#include "util.hh"

#include <sstream>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

// ----------------------------------------------------------------------
fitPsYield::fitPsYield(string fname, int verbose): fVerbose(verbose), fBaseName(fname) {

}


// ----------------------------------------------------------------------
fitPsYield::~fitPsYield() {
}
