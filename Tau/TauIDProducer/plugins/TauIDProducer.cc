// -*- C++ -*-
//
// Package:    Tau/TauIDProducer
// Class:      TauIDProducer
// 
/**\class TauIDProducer TauIDProducer.cc Tau/TauIDProducer/plugins/TauIDProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrew Issac
//         Created:  Thu, 07 Jan 2021 16:44:45 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

namespace nn_inputs{
   constexpr int NumberOfOutputs = 2;
   constexpr int NumberOfInputs = 9;
   namespace TauInputs{
      enum vars{
         Tau_pt = 0,
         Tau_eta,
         Tau_phi,
         Tau_mass,
         Tau_dxy,
         Tau_decayMode,
         Tau_ecalEnergy,
         Tau_hcalEnergy,
         Tau_ip3d
      }
   }
}



//
// class declaration
//

class TauIDProducer : public edm::stream::EDProducer<> {
   public:
      explicit TauIDProducer(const edm::ParameterSet&);
      ~TauIDProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::InputTag src_;
      edm::EDGetTokenT<std::vector<pat::Tau>> slimmedTausToken;

      const std::string graphPath;
      tensorflow::GraphDef graph;
      tensorflow::Session* session;
      std::unique_ptr<tensorflow::Tensor> tauTensor;

      tensorflow::Tensor getPredictions(edm::Event& iEvent, edm::Handle<TauCollection> taus);
      void createTauInputs(const pat::Tau& tau,const size_t& tau_index);
      void createOutputs(edm::Event& iEvent, const tensorflow::Tensor& pred, edm::Handle<TauCollection> slimmedTausCollection)
};

//
// constructors and destructor
//
TauIDProducer::TauIDProducer(const edm::ParameterSet& iConfig)
{
   //now do what ever other initialization is needed
   src_  = iConfig.getParameter<edm::InputTag>( "src" );
   produces<std::vector<bool>>( "tauID" ).setBranchAlias( "tauID" );

   edm::InputTag slimmedTauTag("slimmedTaus");
   slimmedTausToken = consumes<pat::TauCollection>(slimmedTauTag);

   graphFilePath = "/work/aissac/Tau/ML/output_smallDataset_2020-12-14_17-12-20/model_2020-12-14_17-31-57/saved_model.pb"

   tauTensor = std::make_unique<tensorflow::Tensor>(
      tensorflow::DT_FLOAT, tensorflow::TensorShape{1, nn_inputs::NumberOfInputs});
   )

   // load graph and add it to session
   tensorflow::SessionOptions options;
   tensorflow::setThreading(options, 1);
   graphDef = tensorflow::loadGraphDef(graphFilePath);
   session = tensorflow::createSession(graphDef, options);
}


TauIDProducer::~TauIDProducer()
{
   tensorflow::closeSession(session.second);
}



//
// member functions
//

void createTauInputs(const pat::Tau& tau, const size_t& tau_index){
   tensorflow::Tensor& inputs = *tauTensor;
   inputs.flat<float>().setZero();

   const auto& get = [&](int var_index) -> float& { return inputs.matrix<float>()(0, var_index); };

   // TODO: select dxy < -999 || decayMode == 5 || decayMode == 6 here?
   namespace nn = nn_inputs::TauInputs;
   get(nn::Tau_pt) = tau.pt();
   get(nn::Tau_eta) = tau.eta();
   get(nn::Tau_phi) = tau.phi();
   get(nn::Tau_mass) = tau.mass();
   get(nn::Tau_dxy) = tau.dxy();
   get(nn::Tau_decayMode) = tau.decayMode();
   get(nn::Tau_ecalEnergy) = tau.ecalEnergy();
   get(nn::Tau_hcalEnergy) = tau.hcalEnergy();
   get(nn::Tau_ip3d) = tau.ip3d();
}

tensorflow::Tensor getPredictions(edm::Handle<pat::TauCollection> tausCollection){
   tensorflow::Tensor predictions(tensorflow::DT_FLOAT, {static_cast<int>(tausCollection->size()), nn_inputs::NumberOfOutputs});
   
   for (size_t tau_index = 0; tau_index != tausCollection->size(); ++tau_index) {

      std::vector<tensorflow::Tensor> pred_vector;
      createTauInputs(tausCollection->at(tau_index), tau_index);
      tensorflow::run(
            &(session),
            {{"input", *tauTensor}},
            {"predictions"},
            &pred_vector
         );
      for (int k = 0; k < nn_inputs::NumberOfOutputs; ++k) {
         const float pred = pred_vector[0].flat<float>()(k);
         if (!(pred >= 0 && pred <= 1))
         throw cms::Exception("TauIDProducer")
               << "invalid prediction = " << pred << " for tau_index = " << tau_index << ", pred_index = " << k;
         predictions.matrix<float>()(tau_index, k) = pred;
      }
   }
   return predictions;
}

// ------------ method called to produce the data  ------------
void
TauIDProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

      // Get Taus
   edm::Handle<pat::TauCollection> slimmedTausCollection;
   //iEvent.getByLabel(edm::InputTag("slimmedTaus"), slimmedTausCollection); <- for some reason doesn't
   iEvent.getByToken(slimmedTausToken, slimmedTausCollection);
 
   // TODO: use same selection as TauAnalyzer
   // std::vector<int> selectedTausIndices;
   // value_tau_n = 0;
   // for (auto it = slimmedTausCollection->begin(); it != slimmedTausCollection->end(); it++) {
   //    const float dxy = (float)it->dxy();
   //    const float decayMode = (float)it->decayMode();
   //    // Selection rules
   //    if(dxy < -990 || decayMode == 5 || decayMode == 6) continue;

   //    value_tau_n++;
   // }

   const tensorflow::Tensor& pred = getPredictions(iEvent, slimmedTausCollection);
   createOutput(iEvent, pred, slimmedTausCollection);
   
 
}

void TauIDProducer::createOutputs(edm::Event& iEvent, const tensorflow::Tensor& pred, edm::Handle<TauCollection> slimmedTausCollection){
   // TODO: ??? put predictions in event
   auto result = 
   event.put(std::move(result),)
}


// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
TauIDProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
TauIDProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
TauIDProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
TauIDProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TauIDProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TauIDProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauIDProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauIDProducer);
