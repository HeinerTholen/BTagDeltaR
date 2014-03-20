// -*- C++ -*-
//
// Package:    GenPartSelectors
// Class:      GenPartFinalB
// 
/**\class GenPartSelectors GenPartFinalB.cc BTagDeltaR/GenPartSelectors/src/GenPartFinalB.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Heiner Josef Antonius Tholen,68/128,2997,
//         Created:  Wed Mar 19 15:56:17 CET 2014
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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "HepPDT/ParticleID.hh"
//
// class declaration
//

class GenPartFinalB : public edm::EDProducer {
   public:
      explicit GenPartFinalB(const edm::ParameterSet&);
      ~GenPartFinalB();

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
      edm::InputTag src_;
};

//
// constructors and destructor
//
GenPartFinalB::GenPartFinalB(const edm::ParameterSet& iConfig):
    src_(iConfig.getParameter<edm::InputTag>("src"))
{
   produces<std::vector<reco::GenParticle> >();
}


GenPartFinalB::~GenPartFinalB()
{
}

// ------------ method called to produce the data  ------------
void
GenPartFinalB::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;

    Handle<vector<reco::GenParticle> > gens;
    iEvent.getByLabel(src_, gens);

    vector<reco::GenParticle>* out = new vector<reco::GenParticle>();

    for (vector<reco::GenParticle>::const_iterator it = gens->begin(); it != gens->end(); ++it){
        HepPDT::ParticleID pid(it->pdgId());
        if(pid.isHadron() && pid.hasBottom()) {
            bool is_final = true;
            for (unsigned j = 0; j < it->numberOfDaughters(); ++j) {
                HepPDT::ParticleID pidDau(it->daughter(j)->pdgId());
                if(pidDau.isHadron() && pidDau.hasBottom()) {
                    is_final = false;
                    break;
                }
            }
            if (is_final) {
                out->push_back(*it);
            }
        }
    }

    std::auto_ptr<std::vector<reco::GenParticle> > pOut(out);
    iEvent.put(pOut);
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenPartFinalB::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenPartFinalB::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
GenPartFinalB::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
GenPartFinalB::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
GenPartFinalB::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
GenPartFinalB::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenPartFinalB::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenPartFinalB);
