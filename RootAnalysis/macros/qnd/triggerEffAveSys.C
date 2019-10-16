#include "../common/util.hh"

void triggerEffAveSys() {
  std::vector<double> v, e;
  // -- channel 0
  v.push_back(0.043);
  v.push_back(0.023);
  v.push_back(0.026);
  v.push_back(0.024);
  v.push_back(0.010);
  v.push_back(0.017);
  v.push_back(0.020);
  v.push_back(0.014);

  e.push_back(0.001);
  e.push_back(0.001);
  e.push_back(0.001);
  e.push_back(0.001);
  e.push_back(0.001);
  e.push_back(0.001);
  e.push_back(0.001);
  e.push_back(0.001);

  double ave, error, chi2;
  average(ave, error, v, e, chi2);
  cout << "Average channel 0: " << ave << "+/-" << error << endl;

  // -- channel 1
  v.clear(); e.clear();
  v.push_back(0.034);
  v.push_back(0.056);
  v.push_back(0.006);
  v.push_back(0.039);
  v.push_back(0.010);
  v.push_back(0.011);
  v.push_back(0.024);
  v.push_back(0.001);

  e.push_back(0.001);
  e.push_back(0.001);
  e.push_back(0.001);
  e.push_back(0.001);
  e.push_back(0.001);
  e.push_back(0.001);
  e.push_back(0.001);
  e.push_back(0.001);

  average(ave, error, v, e, chi2);
  cout << "Average channel 1: " << ave << "+/-" << error << endl;

}
