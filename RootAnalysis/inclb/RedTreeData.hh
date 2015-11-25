#ifndef REDTREEDATA
#define REDTREEDATA

struct RedTreeData {
  Long64_t run, evt; 
  int      ls, ps;
  int      type, procid;
  
  bool     muid, hlt, hltmatch, json; 
  float    pt, eta, phi, ptrel; 
  float    ip3d, ip2d, ip3dE, ip2dE; 
  
  float    jpt, jeta, jphi; 
  
};

#endif
