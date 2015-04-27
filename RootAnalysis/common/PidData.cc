#include "PidData.hh"

#include "TMath.h"

using namespace std;

ClassImp(PidData)


// ----------------------------------------------------------------------
PidData::PidData(double pmin, double pmax, double tmin, double tmax, double fmin, double fmax, 
		 double e, double s, double pass, double passError, 
		 double tot, double totError, 
		 int mode) {
  _pmin = pmin;
  _tmin = tmin;
  _fmin = fmin; 

  _pmax = pmax;
  _tmax = tmax;
  _fmax = fmax; 

  _e = e;
  _s = s;

  _pass = pass;
  _passE = passError;
  _tot  = tot;
  _totE = totError;
  _mode = mode;

  if (_e < -1.) calcEffAndErr();

}

// ----------------------------------------------------------------------
PidData::PidData(const PidData& orig) {
  _pmin = orig.getPmin();
  _tmin = orig.getTmin();
  _fmin = orig.getFmin();

  _pmax = orig.getPmax();
  _tmax = orig.getTmax();
  _fmax = orig.getFmax();

  _e = orig.getE();
  _s = orig.getS();

  _pass  = orig.getPass();
  _passE = orig.getPassE();

  _tot   = orig.getTot();
  _totE  = orig.getTotE();

  _mode  = orig.getEffAndErrMode();
}


// ----------------------------------------------------------------------
PidData PidData::operator = (const PidData &cmp) {
  setPmax(cmp.getPmax());
  setTmax(cmp.getTmax());
  setFmax(cmp.getFmax());
  
  setPmin(cmp.getPmin());
  setTmin(cmp.getTmin());
  setFmin(cmp.getFmin());
  
  setE(cmp.getE());
  setS(cmp.getS());
  
  setEffAndErrMode(cmp.getEffAndErrMode());

  setPass(cmp.getPass());
  setTot(cmp.getTot());

  setPassE(cmp.getPassE());
  setTotE(cmp.getTotE());
  
  return *this;
}


// ----------------------------------------------------------------------
Bool_t PidData::operator >= (PidData &cmp) const {
  Bool_t answer = (getPmax() >= cmp.getPmax() && getTmax() >= cmp.getTmax() 
		   && getFmax() >= cmp.getFmax());
  return answer;
}


// ----------------------------------------------------------------------
Bool_t PidData::operator < (const PidData &cmp) const {
  Bool_t answer(kFALSE);
  int anInt(0);
  if (getPmax() < cmp.getPmax()) {
    answer = kTRUE;
  } else if (getPmax() == cmp.getPmax()) {
    if (getTmax() < cmp.getTmax()) {
      answer = kTRUE;
    } else if (getTmax() == cmp.getTmax()) {
      if (getFmax() < cmp.getFmax()) {
	answer = kTRUE;
        } else {
          anInt += 1;
        }
    } else {
      anInt += 10;
    }
  } else {
    anInt += 100;
    }
  return answer;
}


// ----------------------------------------------------------------------
Bool_t PidData::operator == (const PidData &cmp) const {
  return ( getPmax() == cmp.getPmax() && getTmax() == cmp.getTmax() 
	   && getFmax() == cmp.getFmax() ); 
}

// ----------------------------------------------------------------------
void PidData::divide(PidData *other) {
  // this/other
  _pass = -1.;
  _tot  = -1.;
  double effo = other->getE();
  double erro = other->getS();
  double e(-99.), s(-99.);
  if (effo > 0.) {
    e = _e/effo;
    s = TMath::Sqrt(((_s*_s)/(effo*effo)) + ((_e*_e*erro*erro)/(effo*effo*effo*effo)));
  } else {
    e = 1.;
    s = 0.;
  }
  _e = e;
  _s = s;
  //  cout << "eff = " << _e << " err = " << _s << endl;
}

// ----------------------------------------------------------------------
void PidData::subtract(PidData *other) {
  // this - other
  _pass = -1.;
  _tot  = -1.;
  double effo = other->getE();
  double erro = other->getS();
  double e(-99.), s(-99.);
  if (effo > 0.) {
    e = _e - effo;
    s = TMath::Sqrt(_s*_s + erro*erro);
  } else {
    e = 1.;
    s = 0.;
  }
  _e = e;
  _s = s;
  //  cout << "eff = " << _e << " err = " << _s << endl;
}
  

// ----------------------------------------------------------------------
void PidData::combine(double pass, double total) {
  if (_pass < 0) {
    cout << "PidData::combine> Warning: combining tables with pass < 0!" << endl;
  }
  _pass += pass;
  _tot  += total;  
  calcEffAndErr();
}
 

// ----------------------------------------------------------------------
void PidData::combine(double pass, double passE, double tot, double totE) {
  if (_pass < 0) {
    cout << "PidData::combine> Warning: combining tables with pass < 0!" << endl;
  }
  _pass += pass;
  _tot  += tot;  
  _passE = TMath::Sqrt(_passE*_passE + passE*passE);
  _totE  = TMath::Sqrt(_totE*_totE + totE*totE);
  calcEffAndErr();
}


// ----------------------------------------------------------------------
void PidData::merge(PidData *other) {
  _pmin = (other->getPmin() < getPmin() ? other->getPmin() : getPmin());
  _pmax = (other->getPmax() > getPmax() ? other->getPmax() : getPmax());
  _fmin = (other->getFmin() < getFmin() ? other->getFmin() : getFmin());
  _fmax = (other->getFmax() > getFmax() ? other->getFmax() : getFmax());
  _tmin = (other->getTmin() < getTmin() ? other->getTmin() : getTmin());
  _tmax = (other->getTmax() > getTmax() ? other->getTmax() : getTmax());

  _pass += other->getPass();
  _passE = TMath::Sqrt(other->getPassE()*other->getPassE() + _passE*_passE);
  _tot += other->getTot();
  _totE  = TMath::Sqrt(other->getTotE()*other->getTotE() + _totE*_totE);

  double e(0.), s(0.);
  double effo = other->getEff();
  double erro = other->getErr();
  if ((_s < 1.e-7) && (_e < 1.e-7)) {
    e = effo;
    s = erro;
  } else if ((_s > 0.) && (erro > 0.)) {
    e = (_e/(_s*_s) + effo/(erro*erro)) / (1./(_s*_s) + 1./(erro*erro));
    s = 1./TMath::Sqrt(1./(_s*_s) + 1./(erro*erro));
    //    cout << "eff1 = " << _e << " + eff2 = " << effo << " -> " << e << endl;
  } else {
    //    cout << "_s or erro == 0" << endl;
  }
  _e = e;
  _s = s;

  calcEffAndErr();
}


// ----------------------------------------------------------------------
void PidData::calcEffAndErr() {
  
  // do nothing for _tot < epsilon. The values have been calculated in combine() and merge()
  if (_tot <= 0.) return; 
    
  if (_mode == 0) {        // d'Agostini
    _e = (_pass+1)/(_tot+2);
    _s = sqrt(((_pass+1)*(_tot-_pass+1))/((_tot+3)*(_tot+2)*(_tot+2)));
  } else if (_mode == 1){  // binomial
    _e = (_pass)/(_tot);
    if (_passE > _totE) {
      _passE = _totE;
    }
    _s = TMath::Sqrt( ((_tot*_tot) - (2*_pass*_tot))*_passE*_passE  + (_pass*_pass)*_totE*_totE ) / ( _tot*_tot );
  } else if (_mode == 99){ // testing 
    _e = (_pass)/(_tot);
    _s = 0.;
  }
  
}


// ----------------------------------------------------------------------
Bool_t PidData::isCell(double p, double t, double f) {
  return ((getPmin() <= p) && (p < getPmax()) 
	  && (getTmin() <= t) && (t < getTmax()) 
	  && (getFmin() <= f) && (f < getFmax()));
}


// ----------------------------------------------------------------------
Bool_t PidData::isZero() {
  return ((TMath::Abs(getE()) < 1.e-6)
	&& (TMath::Abs(getS()) < 1.e-6)
	&& (TMath::Abs(getPass()) < 1.e-6)
	&& (TMath::Abs(getTot()) < 1.e-6));
}


// ----------------------------------------------------------------------
Bool_t PidData::sameValue(PidData *cmp, double eps) {
  return ((TMath::Abs(getE()-cmp->getE()) < eps) 
	  && (TMath::Abs(getS()-cmp->getS()) < eps));
}


// ----------------------------------------------------------------------
void PidData::print() {
  char line[200];
  sprintf(line, "%7.3f %7.3f %5.2f %5.2f %5.2f %5.2f %9.6f %9.6f %11.0f %11.0f %11.0f %11.0f", 
	  getPmin(), getPmax(), 
	  getTmin(), getTmax(), 
	  getFmin(), getFmax(), 
	  getE(),    getS(),
	  getPass(), getPassE(),
	  getTot(),  getTotE());
  cout << line << endl;
  
}


// ===============================================================================
ostream & operator << (ostream& o, const PidData& cmp) {
  char line[200];
  sprintf(line, "%7.3f %7.3f %5.2f %5.2f %5.2f %5.2f %9.6f %9.6f %11.0f %11.0f %11.0f %11.0f", 
	  cmp.getPmin(), cmp.getPmax(), 
	  cmp.getTmin(), cmp.getTmax(), 
	  cmp.getFmin(), cmp.getFmax(), 
	  cmp.getE(), cmp.getS(),
	  cmp.getPass(), cmp.getPassE(),
	  cmp.getTot(), cmp.getTotE());
  return o << line;
  
}
