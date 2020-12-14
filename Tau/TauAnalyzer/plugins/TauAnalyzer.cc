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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
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
      UInt_t value_tau_n;
      Int_t value_run;
      UInt_t value_lumi_block;
      ULong64_t value_event;
      float value_tau_pt;
      float value_tau_eta;
      float value_tau_phi;
      float value_tau_mass;
      float value_tau_dxy;
      float value_tau_decayMode;
      float value_tau_ecalEnergy;
      float value_tau_hcalEnergy;
      float value_tau_ip3d;

      // Generator particles
      // const static int max_gen = 1000;
      // UInt_t value_gen_n;
      // float value_gen_pt[max_gen];
      // float value_gen_eta[max_gen];
      // float value_gen_phi[max_gen];
      // float value_gen_mass[max_gen];
      // float value_gen_pdgid[max_gen]; // typically int
      // float value_gen_status[max_gen]; // typically int
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
TauAnalyzer::TauAnalyzer(const edm::ParameterSet& iConfig) //: isData(iConfig.getParameter<bool>("isData"))
{
   //now do what ever initialization is needed
   edm::InputTag slimmedTauTag("slimmedTaus");
   slimmedTausToken = consumes<pat::TauCollection>(slimmedTauTag);

   // create tree
   edm::Service<TFileService> fs;
   tree = fs->make<TTree>("Events", "Events");

   // Event information
   tree->Branch("run", &value_run);
   tree->Branch("luminosityBlock", &value_lumi_block);
   tree->Branch("event", &value_event);
   // Create tree branches
   // Taus
   tree->Branch("Tau_pt", &value_tau_pt, "Tau_pt/F");
   tree->Branch("Tau_eta", &value_tau_eta, "Tau_eta/F");
   tree->Branch("Tau_phi", &value_tau_phi, "Tau_phi/F");
   tree->Branch("Tau_mass", &value_tau_mass, "Tau_mass/F");
   tree->Branch("Tau_dxy", &value_tau_dxy, "Tau_dxy/F");
   tree->Branch("Tau_decayMode", &value_tau_decayMode, "Tau_decayMode/F");
   tree->Branch("Tau_ecalEnergy", &value_tau_ecalEnergy, "Tau_ecalEnergy/F");
   tree->Branch("Tau_hcalEnergy", &value_tau_hcalEnergy, "Tau_hcalEnergy/F");
   tree->Branch("Tau_ip3d", &value_tau_ip3d, "Tau_ip3d/F");
   // Generator particles
//   if (!isData) {
//       tree->Branch("GenPart_pt", value_gen_pt, "GenPart_pt/F");
//       tree->Branch("GenPart_eta", value_gen_eta, "GenPart_eta/F");
//       tree->Branch("GenPart_phi", value_gen_phi, "GenPart_phi/F");
//       tree->Branch("GenPart_mass", value_gen_mass, "GenPart_mass/F");
//       tree->Branch("GenPart_pdgId", value_gen_pdgid, "GenPart_pdgId/F");
//       tree->Branch("GenPart_status", value_gen_status, "GenPart_status/F");
//   }
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

   // Get Taus
   edm::Handle<pat::TauCollection> slimmedTausCollection;
   //iEvent.getByLabel(edm::InputTag("slimmedTaus"), slimmedTausCollection); <- for some reason doesn't
   iEvent.getByToken(slimmedTausToken, slimmedTausCollection);

   // Event information
   value_run = iEvent.run();
   value_lumi_block = iEvent.luminosityBlock();
   value_event = iEvent.id().event();

   value_tau_n = 0;
   for (auto it = slimmedTausCollection->begin(); it != slimmedTausCollection->end(); it++) {
      // p4 lorentzvector:
      // const auto p4 = it->p4();

      const float dxy = (float)it->dxy();
      const float decayMode = (float)it->decayMode();
      // Selection rules
      if(dxy < -990 || decayMode == 5 || decayMode == 6) continue;

      value_tau_pt = (float)it->pt();
      value_tau_eta = (float)it->eta();
      value_tau_phi= (float)it->phi();
      value_tau_mass = (float)it->mass();
      value_tau_dxy = dxy;
      value_tau_decayMode = decayMode;
      value_tau_ecalEnergy = (float)it->ecalEnergy();
      value_tau_hcalEnergy = (float)it->hcalEnergy();
      value_tau_ip3d = (float)it->ip3d();

      // Fill event to tree
      tree->Fill();

      value_tau_n++;
   }

   // Generator Particles handling
   // if(!isData){
   //    edm::Handle<reco::GenParticleCollection> genParticleCollection;
   //    iEvent.getByLabel(edm::InputTag("genParticles"), genParticleCollection);

   //    value_gen_n = 0;
   //    std::vector<reco::GenParticle> interestingGenParticles;
   //    for(auto it = genParticleCollection->begin(); it != genParticleCollection->end(); it++){
   //       const auto status = it->status();
   //       const auto pdgId = std::abs(it->pdgId());
   //       if (status == 2 && pdgId == 15) { // tau
   //          interestingGenParticles.emplace_back(*it);
   //       }
   //    }

   //    // stole this snippet from Stefan, will be useful later on
   //    // Match taus with gen particles and jets
   //    for (auto p = selectedTaus.begin(); p != selectedTaus.end(); p++) {
   //       // Gen particle matching
   //       auto p4 = p->p4();
   //       auto idx = findBestVisibleMatch(interestingGenParticles, p4); // TODO: Subtract the invisible parts only once.
   //       if (idx != -1) {
   //          auto g = interestingGenParticles.begin() + idx;
   //          value_gen_pt[value_gen_n] = g->pt();
   //          value_gen_eta[value_gen_n] = g->eta();
   //          value_gen_phi[value_gen_n] = g->phi();
   //          value_gen_mass[value_gen_n] = g->mass();
   //          value_gen_pdgid[value_gen_n] = g->pdgId();
   //          value_gen_status[value_gen_n] = g->status();
   //          value_tau_genpartidx[p - selectedTaus.begin()] = value_gen_n;
   //          value_gen_n++;
   //       }

   //       // Jet matching
   //       value_tau_jetidx[p - selectedTaus.begin()] = findBestMatch(selectedJets, p4);
   //    }
   // }


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
