// -*- C++ -*-
//
// Package:    GenPartSelectors
// Class:      GenPartFromZDecay
// 
/**\class GenPartSelectors GenPartFromZDecay.cc BTagDeltaR/GenPartSelectors/src/GenPartFromZDecay.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Heiner Josef Antonius Tholen,68/128,2997,
//         Created:  Wed Mar 19 15:55:40 CET 2014
// $Id$
//
//


// system include files
#include <memory>
#include <vector>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//
// class declaration
//

class GenPartFromZDecay : public edm::EDProducer {
   public:
      explicit GenPartFromZDecay(const edm::ParameterSet&);
      ~GenPartFromZDecay();

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
      bool inverted_;
      edm::InputTag src_;
      std::map<const reco::GenParticle*, bool> hasZMomMap_;
};

//
// constructors and destructor
//
GenPartFromZDecay::GenPartFromZDecay(const edm::ParameterSet& iConfig):
    inverted_(iConfig.getParameter<bool>("inverted")),
    src_(iConfig.getParameter<edm::InputTag>("src"))
{
    produces<std::vector<reco::GenParticle> >();
}


GenPartFromZDecay::~GenPartFromZDecay()
{
}

// ------------ method called to produce the data  ------------
void
GenPartFromZDecay::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;

    Handle<vector<reco::GenParticle> > gens;
    iEvent.getByLabel(src_, gens);

    for (vector<reco::GenParticle>::const_iterator i = gens->begin(); i != gens->end(); ++i){
        if (!i->numberOfMothers()) {
            hasZMomMap_[&(*i)] = false;
        } else if (i->mother()->pdgId() == 23 && i->mother()->status() == 3) {
            hasZMomMap_[&(*i)] = true;
        } else {
            hasZMomMap_[&(*i)] = hasZMomMap_[(const reco::GenParticle*) &(*i->mother())];
        }
    }

    vector<reco::GenParticle>* out = new vector<reco::GenParticle>();
    for (map<const reco::GenParticle*, bool>::const_iterator i = hasZMomMap_.begin();
        i != hasZMomMap_.end();
        ++i) {
        if (!inverted_ && i->second) {
            out->push_back(*i->first);
        } else if (inverted_ && !i->second) {
            out->push_back(*i->first);
        }

    }
    hasZMomMap_.clear();

    std::auto_ptr<std::vector<reco::GenParticle> > pOut(out);
    iEvent.put(pOut);
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenPartFromZDecay::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenPartFromZDecay::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
GenPartFromZDecay::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
GenPartFromZDecay::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
GenPartFromZDecay::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
GenPartFromZDecay::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenPartFromZDecay::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenPartFromZDecay);
