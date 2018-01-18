#ifndef SETUPREADER_h
#define SETUPREADER_h

#include "ReaderData.hh"
#include "preselection.hh"

TMVA::Reader* setupReader(std::string xmlFile, ReaderData &rd);
TMVA::Reader* setupReader(std::string xmlFile, ReaderData &rd, presel &);

#endif
