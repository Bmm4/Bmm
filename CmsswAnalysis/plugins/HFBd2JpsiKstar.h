#ifndef _HFBD2JPSIKSTAR_h_
#define _HFBD2JPSIKSTAR_h_

#include "Bmm/CmsswAnalysis/plugins/HFVirtualDecay.h"

// ----------------------------------------------------------------------
class HFBd2JpsiKstar : public HFVirtualDecay {
	public:
		explicit HFBd2JpsiKstar(const edm::ParameterSet&);

	protected:
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void dumpConfiguration();
		
		int           fPsiMuons;
		double        fPsiWindow, fPhiWindow, fB0Window;


		// insert by jmonroy 8-7-2015
		// additional variables

		double fK0Window

		////////////////////////////////////////////////////
};

#endif
