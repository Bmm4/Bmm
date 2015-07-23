#ifndef REDTREEDATA
#define REDTREEDATA

struct RedTreeData {
  Long64_t run, evt; 
  int      ls;

  bool     muid, hlt, hltmatch, json; 
  int      type;
  float    pt, eta, phi, ptrel; 

  float    jpt, jeta, jphi; 
  
};

#endif
