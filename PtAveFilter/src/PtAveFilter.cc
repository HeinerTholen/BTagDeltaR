// -*- C++ -*-
//
// Package:    PtAveFilter
// Class:      PtAveFilter
// 
/**\class PtAveFilter PtAveFilter.cc BTagDeltaR/PtAveFilter/src/PtAveFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Heiner Josef Antonius Tholen,68/128,2997,
//         Created:  Wed May  7 16:38:51 CEST 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"

//
// class declaration
//

class PtAveFilter : public edm::EDFilter {
   public:
      explicit PtAveFilter(const edm::ParameterSet&);
      ~PtAveFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
      edm::InputTag src_;
      double cut_;
};

//
// constructors and destructor
//
PtAveFilter::PtAveFilter(const edm::ParameterSet& iConfig):
    src_(iConfig.getParameter<edm::InputTag>("src")),
    cut_(iConfig.getParameter<double>("cut"))
{
}


PtAveFilter::~PtAveFilter()
{
}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
PtAveFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<reco::CandidateView> src;
    iEvent.getByLabel(src_, src);

    if (src->size() > 1 && (src->at(0).pt()+src->at(1).pt()) > 2.*cut_) {
        return true;
    } else {
        return false;
    }
}

// ------------ method called once each job just before starting event loop  ------------
void 
PtAveFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PtAveFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
PtAveFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
PtAveFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
PtAveFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
PtAveFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PtAveFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(PtAveFilter);
