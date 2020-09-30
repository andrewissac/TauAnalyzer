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
#include "TTree.h"
#include "TFile.h"
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

      TTree *tree;

      // Taus
      const static int max_tau = 1000;
      UInt_t value_tau_n;
      float value_tau_pt[max_tau];
      float value_tau_eta[max_tau];
      float value_tau_phi[max_tau];
      float value_tau_mass[max_tau];
      float value_tau_dxy[max_tau];
      float value_tau_decayMode[max_tau];
      float value_tau_ecalEnergy[max_tau];
      float value_tau_hcalEnergy[max_tau];
      float value_tau_ip3d[max_tau];
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
{
   //now do what ever initialization is needed
   edm::InputTag slimmedTauTag("slimmedTaus");
   slimmedTausToken = consumes<pat::TauCollection>(slimmedTauTag);
   // create tree
   edm::Service<TFileService> fs;
   tree = fs->make<TTree>("Events", "Events");

   // Create tree branches
   tree->Branch("Tau_pt", value_tau_pt, "Tau_pt/F");
   tree->Branch("Tau_eta", value_tau_eta, "Tau_eta/F");
   tree->Branch("Tau_phi", value_tau_phi, "Tau_phi/F");
   tree->Branch("Tau_mass", value_tau_mass, "Tau_mass/F");
   tree->Branch("Tau_dxy", value_tau_dxy, "Tau_dxy/F");
   tree->Branch("Tau_decayMode", value_tau_decayMode, "Tau_decayMode/F");
   tree->Branch("Tau_ecalEnergy", value_tau_ecalEnergy, "Tau_ecalEnergy/F");
   tree->Branch("Tau_hcalEnergy", value_tau_hcalEnergy, "Tau_hcalEnergy/F");
   tree->Branch("Tau_ip3d", value_tau_ip3d, "Tau_ip3d/F");
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

   value_tau_n = 0;
   for (auto it = slimmedTausCollection->begin(); it != slimmedTausCollection->end(); it++) {
      // how toget p4 lorentzvector:
      // const auto p4 = it->p4();
      // float pt = (float)p4.pt();
      // float eta = (float)p4.eta();
      // float phi = (float)p4.phi();
      // float mass = (float)p4.mass();

      value_tau_pt[value_tau_n] = (float)it->pt();
      value_tau_eta[value_tau_n] = (float)it->eta();
      value_tau_phi[value_tau_n] = (float)it->phi();
      value_tau_mass[value_tau_n] = (float)it->mass();
      value_tau_dxy[value_tau_n] = (float)it->dxy();
      value_tau_decayMode[value_tau_n] = (float)it->decayMode();
      value_tau_ecalEnergy[value_tau_n] = (float)it->ecalEnergy();
      value_tau_hcalEnergy[value_tau_n] = (float)it->hcalEnergy();
      value_tau_ip3d[value_tau_n] = (float)it->ip3d();

      value_tau_n++;
   }

   // Fill event to tree
   tree->Fill();


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
