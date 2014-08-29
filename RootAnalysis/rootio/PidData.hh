#ifndef PIDDATA
#define PIDDATA

#include "TString.h"
#include "TObject.h"
#include <iostream>
#include <fstream>

class PidData: public TObject {

public:

	  
  PidData() {_pmax=_pmin=_tmin=_tmax=_fmin=_fmax=_e=_s=0.;_pass=_tot=0.;};
  // -- Force recalculation of eff and error by calling with  eff < -1.  

  PidData(double p, double pmax, double th, double tmax, double phi, double fmax, 
	  double e = 0., double s = 0., 
	  double pass = 0., double passError = 0., double tot = 0., double totError = 0., 
	  int mode = 0);
	  
  PidData(const PidData& orig);

  ~PidData() {;}

  PidData operator = (const PidData &cmp); 
  Bool_t operator >= (PidData &cmp) const;
  Bool_t operator < (const PidData &cmp) const;
  Bool_t operator == ( const PidData& cmp ) const;


  Bool_t isCell(double p, double t, double f);  
  void   calcEffAndErr();
  void   print();
  void   combine(double pass, double tot);
  void   combine(double pass, double passE, double tot, double totE);
  void   merge(PidData *other);
  void   divide(PidData *other);
  void   subtract(PidData *other);
  Bool_t isZero();
  Bool_t sameValue(PidData *, double epsilon = 1.e-3);
  
  double getPctr() const  {return (_pmax + _pmin)/2.;}
  double getTctr() const  {return (_tmax + _tmin)/2.;}
  double getFctr() const  {return (_fmax + _fmin)/2.;}

  double getPmax() const  {return _pmax;}
  double getTmax() const  {return _tmax;}
  double getFmax() const  {return _fmax;}

  double getPmin() const  {return _pmin;}
  double getTmin() const  {return _tmin;}
  double getFmin() const  {return _fmin;}

  double getE() const     {return _e;}
  double getEff() const   {return getE();}
  double getS() const     {return _s;}
  double getErr() const   {return getS();}
  double getPass() const  {return _pass;}
  double getPassE() const {return _passE;}
  double getTot() const   {return _tot;}
  double getTotE() const  {return _totE;}

  void   setPmax(double p)     { _pmax = p;}
  void   setTmax(double th)   { _tmax = th;}
  void   setFmax(double phi) { _fmax = phi;}

  void   setPmin(double p)   { _pmin = p;}
  void   setTmin(double th)  { _tmin = th;}
  void   setFmin(double phi) { _fmin = phi;}

  void   setE(double x)     { _e = x;}
  void   setS(double x)     { _s = x;}
  void   setPass(double x)  { _pass = x;}
  void   setPassE(double x) { _passE = x;}
  void   setTot(double x)   { _tot = x;}
  void   setTotE(double x)  { _totE = x;}

  int    getEffAndErrMode() const  {return _mode;}
  void   setEffAndErrMode(int i) { _mode = i;}

private:
  double _pmin, _pmax;
  double _tmin, _tmax;
  double _fmin, _fmax;
  double _e, _s;
  double _pass, _tot;
  double _passE, _totE;
  int _mode;

  ClassDef(PidData,1) //Testing PidData

};

std::ostream & operator << (std::ostream& , const PidData&);


#endif
