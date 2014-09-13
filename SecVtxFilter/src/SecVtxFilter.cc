// -*- C++ -*-
//
// Package:    SecVtxFilter
// Class:      SecVtxFilter
// 
/**\class SecVtxFilter SecVtxFilter.cc BTagDeltaR/SecVtxFilter/src/SecVtxFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Heiner Josef Antonius Tholen,68/128,2997,
//         Created:  Thu Jun 19 18:21:54 CEST 2014
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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
//
// class declaration
//

class SecVtxFilter : public edm::EDFilter {
   public:
      explicit SecVtxFilter(const edm::ParameterSet&);
      ~SecVtxFilter();

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
      edm::InputTag pv_src_;
      edm::InputTag sv_src_;
      double min_3d_significance_;
};

//
// constructors and destructor
//
SecVtxFilter::SecVtxFilter(const edm::ParameterSet& iConfig):
    pv_src_(iConfig.getParameter<edm::InputTag>("pv_src")),
    sv_src_(iConfig.getParameter<edm::InputTag>("sv_src")),
    min_3d_significance_(iConfig.getParameter<double>("min_3d_significance"))
{
    produces<std::vector<reco::Vertex> >();
}


SecVtxFilter::~SecVtxFilter()
{
}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
SecVtxFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<std::vector<reco::Vertex> > sv_src;
    iEvent.getByLabel(sv_src_, sv_src);

    edm::Handle<std::vector<reco::Vertex> > pv_src;
    iEvent.getByLabel(pv_src_, pv_src);
    const reco::Vertex &pv = pv_src->at(0);

    std::vector<reco::Vertex>* out = new std::vector<reco::Vertex>();
    for (unsigned i=0; i<sv_src->size(); ++i) {
        const reco::Vertex sv = sv_src->at(i);
    	GlobalVector fd = GlobalVector(
    	    sv.x() - pv.x(),
            sv.y() - pv.y(),
            sv.z() - pv.z()
        );
        reco::SecondaryVertex sv_cls = reco::SecondaryVertex(pv, sv, fd, true);
        if (sv_cls.dist3d().significance() > min_3d_significance_) {
            out->push_back(sv);
        }
    }

    std::auto_ptr<std::vector<reco::Vertex> > pOut(out);
    iEvent.put(pOut);

    return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
SecVtxFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SecVtxFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
SecVtxFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
SecVtxFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
SecVtxFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
SecVtxFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SecVtxFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(SecVtxFilter);
