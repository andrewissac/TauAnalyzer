import FWCore.ParameterSet.Config as cms
from RecoTauTag.RecoTau.tools import runTauIdMVA

updatedTauName = "slimmedTausNewID"

process = cms.Process("TauID")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/6A62DE01-E641-E811-AFF1-008CFAC94038.root'
    )
)

tauIdEmbedder = runTauIdMVA.TauIDEmbedder(
    process, cms, updatedTauName = updatedTauName,
    toKeep = ["2017v2", "deepTau2017v2p1"]
    )

tauIdEmbedder.runTauID()
tauSrc_InputTag = cms.InputTag('slimmedTausNewID')# to be taken for any n-tuplizer

process.out = cms.OutputModule(
    "PoolOutputModule",fileName = cms.untracked.string("output_TauID.root"),
    outputCommands = cms.untracked.vstring("drop *", "keep *_" + updatedTauName + "_*_*")
    )

process.p = cms.Path(process.rerunMvaIsolationSequence*process.slimmedTausNewID)# put prior to n-tuplizer
process.ep = cms.EndPath(process.out)