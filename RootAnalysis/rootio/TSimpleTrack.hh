#ifndef TSIMPLETRACK
#define TSIMPLETRACK


#include "TObject.h"
#include "TVector3.h"

class TSimpleTrack: public TObject {

public:
  TSimpleTrack(int index = -1);
  ~TSimpleTrack() { };
  void clear(); 
  void dump();

  // -- if the indices run into overflow (or are invalid) they are set to the maximum value
  enum IDX {IDXTRK = 0,  // 12 bits for the track index, maximum = 4095
	    IDXGEN = 12, // 12 bits for the gen-block index, maximum = 4095
	    IDXPV = 24   //  8 bits for the PV index, maximum = 255
  };
  
  enum BIT {BITCHARGE = 0,  // charge (1 = positive, 0 = negative)
	BITHP = 1,      // highPurity (1 = highPurity, 0 = not highPurity)
	BITMU = 2       // muon ID (just one bit if this track is in the muon list)
  };

  TVector3 getP()           {return TVector3(fPx, fPy, fPz);}
  int getIndex()            {return (fIndices&0xFFF);}
  int getCharge()           {return (fBits&(1<<BITCHARGE)?1:-1);}
  int getHighPurity()       {return ((fBits&(1<<BITHP)) >> BITHP);}
  int getMuonID()           {return ((fBits&(1<<BITMU)) >> BITMU);}
  int getPvIndex()          {int r = ((fIndices&(0xFF<<IDXPV)) >> IDXPV); return (r<0xFF?r:-1);}
  int getGenIndex()         {int r = ((fIndices&(0xFFF<<IDXGEN)) >> IDXGEN); return (r<0xFFF?r:-1);}
  int getIndices()          {return fIndices;}
  int getBits()             {return fBits;}
  
  void setP(double x, double y, double z) {fPx = x; fPy = y; fPz = z;}
  void setIndex(int idx)    {if (idx<(1<<12)) fIndices |= idx; else fIndices |= 0xFFF;}
  void setCharge(int q)     {if (q>0) fBits |= (1<<BITCHARGE); else fBits &= ~(1<<BITCHARGE);}
  void setHighPurity(int f) {if (f>0) fBits |= (1<<BITHP); else fBits &= ~(1<<BITHP);}
  void setMuonID(int id)    {if (id>0)fBits |= (1<<BITMU); else fBits &= ~(1<<BITMU);}
  void setPvIndex(int id)   {if (id<(1<<8)) fIndices |= (id<<IDXPV); else fIndices |= (0xFF<<IDXPV);}
  void setGenIndex(int idx) {if (idx > -1 && idx<(1<<12)) fIndices |= (idx<<IDXGEN); else fIndices |= (0xFFF<<IDXGEN);}

  
private:

  int   fBits, fIndices; 
  float fPx, fPy, fPz; 

//   Int_t    fBits, fIndices; 
//   Float_t  fPx, fPy, fPz; 

  ClassDef(TSimpleTrack,1)

};

#endif
