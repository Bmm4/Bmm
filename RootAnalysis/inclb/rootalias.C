void printAxisLabels(TH1 *h) {
  for (int i = 1; i < h->GetNbinsX(); ++i) {
    cout << h->GetXaxis()->GetBinLabel(i) << endl;
  }
}
