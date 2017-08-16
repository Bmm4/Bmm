#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TMVA/Reader.h"

#include "setupReader.hh"

using namespace std;

// ----------------------------------------------------------------------
TMVA::Reader* setupReader(string xmlFile, ReaderData &rd) {
  TMVA::Reader *reader = new TMVA::Reader( "!Color:Silent" );

  // -- read in variables from weight file
  vector<string> allLines;
  char  buffer[2000];
  cout << "setupReader, open file " << xmlFile ;
  ifstream is(xmlFile);
  while (is.getline(buffer, 2000, '\n')) allLines.push_back(string(buffer));
  int nvars(-1);
  string::size_type m1, m2;
  string stype;
  //  cout << "  read " << allLines.size() << " lines " << endl;
  for (unsigned int i = 0; i < allLines.size(); ++i) {
    // -- parse and add variables
    if (string::npos != allLines[i].find("Variables NVar")) {
      m1 = allLines[i].find("=\"");
      stype = allLines[i].substr(m1+2, allLines[i].size()-m1-2-2);
      //      cout << "  " << stype << " variables" << endl;
      nvars = atoi(stype.c_str());
      if (-1 == nvars) continue;
      for (unsigned int j = i+1; j < i+nvars+1; ++j) {
	m1 = allLines[j].find("Expression=\"")+10;
	m2 = allLines[j].find("\" Label=\"");
	stype = allLines[j].substr(m1+2, m2-m1-2);
	//	cout << "ivar " << j-i << " variable string: ->" << stype << "<-" << endl;
	if (stype == "m1pt") {
	  //	  cout << "  adding m1pt" << endl;
	  reader->AddVariable( "m1pt", &rd.m1pt);
	}
	if (stype == "m2pt") {
	  //	  cout << "  adding m2pt" << endl;
	  reader->AddVariable( "m2pt", &rd.m2pt);
	}
	if (stype == "m1eta") {
	  //	  cout << "  adding m1eta" << endl;
	  reader->AddVariable( "m1eta", &rd.m1eta);
	}
	if (stype == "m2eta") {
	  reader->AddVariable( "m2eta", &rd.m2eta);
	  //	  cout << "  adding m2eta" << endl;
	}
	if (stype == "pt") {
	  //	  cout << "  adding pt" << endl;
	  reader->AddVariable( "pt", &rd.pt);
	}
	if (stype == "eta") {
	  //	  cout << "  adding eta" << endl;
	  reader->AddVariable( "eta", &rd.eta);
	}
	if (stype == "fls3d") {
	  //	  cout << "  adding fls3d" << endl;
	  reader->AddVariable( "fls3d", &rd.fls3d);
	}
	if (stype == "alpha") {
	  //	  cout << "  adding alpha" << endl;
	  reader->AddVariable( "alpha", &rd.alpha);
	}
	if (stype == "maxdoca") {
	  //	  cout << "  adding maxdoca" << endl;
	  reader->AddVariable( "maxdoca", &rd.maxdoca);
	}
	if (stype == "pvip") {
	  ///	  cout << "  adding pvip" << endl;
	  reader->AddVariable( "pvip", &rd.pvip);
	}
	if (stype == "pvips") {
	  //	  cout << "  adding pvips" << endl;
	  reader->AddVariable( "pvips", &rd.pvips);
	}
	if (stype == "iso") {
	  //	  cout << "  adding iso" << endl;
	  reader->AddVariable( "iso", &rd.iso);
	}
	if (stype == "docatrk") {
	  //	  cout << "  adding docatrk" << endl;
	  reader->AddVariable( "docatrk", &rd.docatrk);
	}
	if (stype == "closetrk") {
	  //	  cout << "  adding closetrk" << endl;
	  reader->AddVariable( "closetrk", &rd.closetrk);
	}
	if (stype == "closetrks1") {
	  //	  cout << "  adding closetrks1" << endl;
	  reader->AddVariable( "closetrks1", &rd.closetrks1);
	}
	if (stype == "closetrks2") {
	  //	  cout << "  adding closetrks2" << endl;
	  reader->AddVariable( "closetrks2", &rd.closetrks2);
	}
	if (stype == "closetrks3") {
	  //	  cout << "  adding closetrks3" << endl;
	  reader->AddVariable( "closetrks3", &rd.closetrks3);
	}
	if (stype == "chi2dof") {
	  //	  cout << "  adding chi2dof" << endl;
	  reader->AddVariable( "chi2dof", &rd.chi2dof);
	}
	if (stype == "m1iso") {
	  //	  cout << "  adding m1iso" << endl;
	  reader->AddVariable( "m1iso", &rd.m1iso);
	}
	if (stype == "m2iso") {
	  //	  cout << "  adding m2iso" << endl;
	  reader->AddVariable( "m2iso", &rd.m2iso);
	}
	if (stype == "pvdchi2") {
	  //	  cout << "  adding pvdchi2" << endl;
	  reader->AddVariable( "pvdchi2", &rd.pvdchi2);
	}
	if (stype == "othervtx") {
	  //	  cout << "  adding othervtx" << endl;
	  reader->AddVariable( "othervtx", &rd.othervtx);
	}
	if (stype == "pvlips") {
	  //	  cout << "  adding pvlips" << endl;
	  reader->AddVariable( "pvlips", &rd.pvlips);
	}
	if (stype == "pvlip") {
	  //	  cout << "  adding pvlip" << endl;
	  reader->AddVariable( "pvlip", &rd.pvlip);
	}
	if (stype == "pv2lips") {
	  //	  cout << "  adding pv2lips" << endl;
	  reader->AddVariable( "pv2lips", &rd.pv2lips);
	}
	if (stype == "pv2lip") {
	  //	  cout << "  adding pv2lip" << endl;
	  reader->AddVariable( "pv2lip", &rd.pv2lip);
	}
      }
      break;
    }
  }

  nvars = -1;
  for (unsigned int i = 0; i < allLines.size(); ++i) {
    // -- parse and add spectators
    if (string::npos != allLines[i].find("Spectators NSpec")) {
      m1 = allLines[i].find("=\"");
      stype = allLines[i].substr(m1+2, allLines[i].size()-m1-2-2);
      //      cout << "==> " << stype << endl;
      nvars = atoi(stype.c_str());
      if (-1 == nvars) continue;
      for (unsigned int j = i+1; j < i+nvars+1; ++j) {
	m1 = allLines[j].find("Expression=\"")+10;
	m2 = allLines[j].find("\" Label=\"");
	stype = allLines[j].substr(m1+2, m2-m1-2);
	//	cout << "ivar " << j-i << " spectator string: ->" << stype << "<-" << endl;
	if (stype == "m") {
	  //	  cout << "  adding m as spectator" << endl;
	  reader->AddSpectator( "m", &rd.m);
	}
      }
      break;
    }
  }

  // --- Book the MVA methods
  reader->BookMVA("BDT", xmlFile);
  cout << ",  booked reader " << reader << endl;
  return reader;
}
