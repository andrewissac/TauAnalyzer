// -*- C++ -*-
//
// Package:    Tau/TauAnalyzer
// Class:      TauAnalyzer
//
/**\class TauAnalyzer TauAnalyzer.cc Tau/TauAnalyzer/plugins/TauAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrew Issac
//         Created:  Tue, 29 Sep 2020 13:49:42 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/Utilities/interface/InputTag.h"
 #include "DataFormats/TrackReco/interface/Track.h"
 #include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class TauAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TauAnalyzer(const edm::ParameterSet&);
      ~TauAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      //edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<std::vector<pat::Tau>> slimmedTausToken;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TauAnalyzer::TauAnalyzer(const edm::ParameterSet& iConfig)
 //:
  //tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))
{
   //now do what ever initialization is needed
   edm::InputTag slimmedTauTag("slimmedTaus");
   slimmedTausToken = consumes<pat::TauCollection>(slimmedTauTag);
}


TauAnalyzer::~TauAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   //using namespace edm;

    edm::Handle<pat::TauCollection> slimmedTausCollection;
    iEvent.getByToken(slimmedTausToken, slimmedTausCollection);

   int i = 0;
   for (auto&& tau : *(slimmedTausCollection.product())){
      float decayMode = (float)tau.decayMode();
      float dxy = (float)tau.dxy();
      //float dz = (float)tau.dz();
      //auto correctedP4 = new reco::Candidate::PolarLorentzVector(tau.correctedP4(tau.currentJECLevel()));
      //reco::LeafCandidate::LorentzVector correctedP4 = tau.correctedP4(tau.currentJECLevel());
      std::cout << correctedP4 << std::endl;
      // float pt -> correctedP4.pt();
      // float eta -> correctedP4.eta();
      // float phi -> correctedP4.phi();
      // float m -> correctedP4.mass();
      float ecalEnergy = (float)tau.ecalEnergy();
      float hcalEnergy = (float)tau.hcalEnergy();
      float ip3d = (float)tau.ip3d();

      if(i > 90000){
         std::cout << "decayMode: " << decayMode << std::endl;
         std::cout << "dxy: " << dxy << std::endl;
         //std::cout << "dz: " << dz << std::endl;
         // std::cout << "pt: " << pt << std::endl;
         // std::cout << "eta: " << eta << std::endl;
         // std::cout << "phi: " << phi << std::endl;
         // std::cout << "m: " << m << std::endl;
         std::cout << "ecalEnergy: " << ecalEnergy << std::endl;
         std::cout << "hcalEnergy: " << hcalEnergy << std::endl;
         std::cout << "ip3d: " << ip3d << std::endl;
      }
      i++;
   }


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
TauAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TauAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauAnalyzer);
