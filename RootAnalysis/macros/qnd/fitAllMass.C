void fitAllMass(string sample = "bupsikData", string sel = "Presel") {
  // plotReducedOverlays a2016BF("results", "plotResults.2016BF.files", "baseCuts.2016.cuts", "BF");
  // plotReducedOverlays a2016GH("results", "plotResults.2016GH.files", "baseCuts.2016.cuts", "GH");

  // plotReducedOverlays a2012("results", "plotResults.2012.files", "baseCuts.2012.cuts");

  plotReducedOverlays a2011("results", "plotResults.2011.files", "baseCuts.2011.cuts");
  //  a2011.makeOverlay("bspsiphiData", "bspsiphiMcComb", "bdt");
  a2011.makeOverlay("bmmData", "bdmmMcComb", "bdt");

  // a2016BF.plotMass(sample, sel);
  // a2016GH.plotMass(sample, sel);

  // a2012.plotMass(sample, sel);
  // a2011.plotMass(sample, sel);

}
