// -*- C++ -*-
//
// Package:    SecVtxInfoProducer
// Class:      SecVtxInfoProducer
// 
/**\class SecVtxInfoProducer SecVtxInfoProducer.cc BTagDeltaR/SecVtxInfoProducer/src/SecVtxInfoProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Heiner Josef Antonius Tholen,68/128,2997,
//         Created:  Thu Jun 26 13:55:54 CEST 2014
// $Id$
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

//
// class declaration
//

class SecVtxInfoProducer : public edm::EDProducer {
   public:
      explicit SecVtxInfoProducer(const edm::ParameterSet&);
      ~SecVtxInfoProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
      edm::InputTag pv_src_;
      edm::InputTag sv_src_;
};

//
// constructors and destructor
//
SecVtxInfoProducer::SecVtxInfoProducer(const edm::ParameterSet& iConfig):
    pv_src_(iConfig.getParameter<edm::InputTag>("pv_src")),
    sv_src_(iConfig.getParameter<edm::InputTag>("sv_src"))
{
    produces<std::vector<double> >("DxyVal"); 
    produces<std::vector<double> >("DxyErr"); 
    produces<std::vector<double> >("DxyzVal"); 
    produces<std::vector<double> >("DxyzErr"); 
}


SecVtxInfoProducer::~SecVtxInfoProducer()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
SecVtxInfoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<std::vector<reco::Vertex> > sv_src;
    iEvent.getByLabel(sv_src_, sv_src);

    edm::Handle<std::vector<reco::Vertex> > pv_src;
    iEvent.getByLabel(pv_src_, pv_src);
    const reco::Vertex &pv = pv_src->at(0);

    std::vector<double>* out_Dxy = new std::vector<double>();
    std::vector<double>* out_Dxy_err = new std::vector<double>();
    std::vector<double>* out_Dxyz = new std::vector<double>();
    std::vector<double>* out_Dxyz_err = new std::vector<double>();
    for (unsigned i=0; i<sv_src->size(); ++i) {
        const reco::Vertex sv = sv_src->at(i);
        GlobalVector fd = GlobalVector(
            sv.x() - pv.x(),
            sv.y() - pv.y(),
            sv.z() - pv.z()
        );
        reco::SecondaryVertex sv_cls = reco::SecondaryVertex(pv, sv, fd, true);
        out_Dxy->push_back(sv_cls.dist2d().value());
        out_Dxy_err->push_back(sv_cls.dist2d().error());
        out_Dxyz->push_back(sv_cls.dist3d().value());
        out_Dxyz_err->push_back(sv_cls.dist3d().error());
    }

    iEvent.put(std::auto_ptr<std::vector<double> >(out_Dxy), "DxyVal");
    iEvent.put(std::auto_ptr<std::vector<double> >(out_Dxy_err), "DxyErr");
    iEvent.put(std::auto_ptr<std::vector<double> >(out_Dxyz), "DxyzVal");
    iEvent.put(std::auto_ptr<std::vector<double> >(out_Dxyz_err), "DxyzErr");
}

// ------------ method called once each job just before starting event loop  ------------
void 
SecVtxInfoProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SecVtxInfoProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
SecVtxInfoProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
SecVtxInfoProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
SecVtxInfoProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
SecVtxInfoProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SecVtxInfoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SecVtxInfoProducer);
