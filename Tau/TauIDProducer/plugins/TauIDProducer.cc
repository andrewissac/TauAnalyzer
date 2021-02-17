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
#include "DataFormats/PatCandidates/interface/PATTauDiscriminator.h"

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
      };
   }
}



//
// class declaration
//

class TauIDProducer : public edm::stream::EDProducer<> {
   public:
      using TauDiscriminator = pat::PATTauDiscriminator;
      using TauCollection = std::vector<pat::Tau>;
      using TauRef = edm::Ref<TauCollection>;
      using TauRefProd = edm::RefProd<TauCollection>;
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
      edm::EDGetTokenT<TauCollection> slimmedTausToken;

      std::string graphFilePath;
      tensorflow::GraphDef* graphDef;
      tensorflow::Session* session;
      std::unique_ptr<tensorflow::Tensor> tauInputTensor;

      tensorflow::Tensor getPredictions(edm::Handle<TauCollection> taus);
      void createTauInputs(const pat::Tau& tau,const size_t& tau_index);
      void createOutputs(edm::Event& iEvent, const tensorflow::Tensor& pred, edm::Handle<TauCollection> slimmedTausCollection);
      std::unique_ptr<pat::PATTauDiscriminator> getOutputFromPredTensor(
         const edm::Handle<TauCollection>& taus,
         const tensorflow::Tensor& pred) const;
};

//
// constructors and destructor
//
TauIDProducer::TauIDProducer(const edm::ParameterSet& iConfig)
{
   //now do what ever other initialization is needed

   edm::InputTag slimmedTauTag("slimmedTausNewID");
   slimmedTausToken = consumes<pat::TauCollection>(slimmedTauTag);

   graphFilePath = "/work/aissac/Tau/ML/output_smallDataset_2020-12-14_17-12-20/model_2021-01-20_00-31-38/tf_model.pb";

   tauInputTensor = std::make_unique<tensorflow::Tensor>(tensorflow::DT_FLOAT, tensorflow::TensorShape{1, nn_inputs::NumberOfInputs});

   // load graph and add it to session
   tensorflow::SessionOptions options;
   tensorflow::setThreading(options, 1);
   graphDef = tensorflow::loadGraphDef(graphFilePath);
   session = tensorflow::createSession(graphDef, options);

   produces<TauDiscriminator>();
}


TauIDProducer::~TauIDProducer()
{
   tensorflow::closeSession(session);
}



//
// member functions
//

void TauIDProducer::createTauInputs(const pat::Tau& tau, const size_t& tau_index){
   tensorflow::Tensor& inputs = *tauInputTensor;
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

tensorflow::Tensor TauIDProducer::getPredictions(edm::Handle<TauCollection> tausCollection){
   tensorflow::Tensor predictions(tensorflow::DT_FLOAT, {static_cast<int>(tausCollection->size()), nn_inputs::NumberOfOutputs});
   for (size_t tau_index = 0; tau_index != tausCollection->size(); ++tau_index) {
      std::vector<tensorflow::Tensor> pred_vector;
      createTauInputs(tausCollection->at(tau_index), tau_index);
      tensorflow::run(
            session,
            // {{"input", *tauInputTensor}},
            // {"predictions"},
            {{"x:0", *tauInputTensor}}, // python train.py -> freezing the graph changed the input layer name to x:0 
            {"Identity:0"}, // python train.py -> freezing the graph changed the input layer name to Identity:0 
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

std::unique_ptr<pat::PATTauDiscriminator> TauIDProducer::getOutputFromPredTensor(const edm::Handle<TauCollection>& taus,
                                                                                       const tensorflow::Tensor& pred) const {
   std::unique_ptr<pat::PATTauDiscriminator> output = std::make_unique<pat::PATTauDiscriminator>(edm::RefProd<TauCollection>(taus));
   for (size_t tau_index = 0; tau_index < taus->size(); ++tau_index) {
      // TODO: check if (tau_index, 0) or (tau_index, 1)
      float x = pred.matrix<float>()(tau_index, 0);
      output->setValue(tau_index, x);
   }
   return output;
}

void TauIDProducer::createOutputs(edm::Event& iEvent, const tensorflow::Tensor& pred, edm::Handle<TauCollection> slimmedTausCollection){
   auto result = getOutputFromPredTensor(slimmedTausCollection, pred);
   iEvent.put(std::move(result));
}

// ------------ method called to produce the data  ------------
void
TauIDProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // Get Taus
   edm::Handle<TauCollection> slimmedTausCollection;
   iEvent.getByToken(slimmedTausToken, slimmedTausCollection);
 
   // TODO: use same selection as TauAnalyzer?
   // std::vector<int> selectedTausIndices;
   // value_tau_n = 0;
   // for (auto it = slimmedTausCollection->begin(); it != slimmedTausCollection->end(); it++) {
   //    const float dxy = (float)it->dxy();
   //    const float decayMode = (float)it->decayMode();
   //    // Selection rules
   //    if(dxy < -990 || decayMode == 5 || decayMode == 6) continue;

   //    value_tau_n++;
   // }

   const tensorflow::Tensor& pred = getPredictions(slimmedTausCollection);
   createOutputs(iEvent, pred, slimmedTausCollection);
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
