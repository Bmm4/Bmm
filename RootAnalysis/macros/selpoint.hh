#ifndef SELPOINT_H
#define SELPOINT_H

#include <vector>

class selpoint {
public:
  selpoint();
  ~selpoint();

  void eval(int cat = 0, double w8 = 1.); 

  double fCnt[10]; // 0: signal, 1: first background, 2: second background, ...

  std::vector<std::pair<double *, double> > fSmallerThan;  // <variable pointer, cut value>
  std::vector<std::pair<double *, double> > fLargerThan;  // <variable pointer, cut value>


};

#endif
