import FWCore.ParameterSet.Config as cms
from RecoTauTag.RecoTau.tools import runTauIdMVA
import sys
from os import path

updatedTauName = "slimmedTausNewID"
tauIDProducerName = "TauIDProducer"

process = cms.Process("TauID")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

tauIdEmbedder = runTauIdMVA.TauIDEmbedder(
    process, cms, updatedTauName = updatedTauName,
    toKeep = ["2017v2", "deepTau2017v2p1"]
    )

tauIdEmbedder.runTauID()
tauSrc_InputTag = cms.InputTag('slimmedTausNewID')# to be taken for any n-tuplizer

# arguments_jobID = int(sys.argv[2]) 
# inputfile = str(sys.argv[3]) # this takes the filename of the root file passed in from job.sh from running cmsRun mypath/ConfFile_cfg.py myfile.root (in this case: sys.argv[3] = myfile.root )
# dataset = str(sys.argv[4])
# outputdir = str(sys.argv[5])

# process.source = cms.Source("PoolSource",
#     fileNames = cms.untracked.vstring(inputfile)
# )

# process.out = cms.OutputModule(
#     "PoolOutputModule",fileName = cms.untracked.string(path.join(outputdir,"TauIDProducer_output_{:04d}_{}.root".format(arguments_jobID, dataset))),
#     outputCommands = cms.untracked.vstring("drop *", "keep *_" + updatedTauName + "_*_*")
#     )

TauIDProducer = cms.EDProducer("TauIDProducer")
setattr(process, tauIDProducerName, TauIDProducer)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # genuine tau
    #    '/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/6A62DE01-E641-E811-AFF1-008CFAC94038.root'
        # tau fake
        "/store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/70000/782395CC-8C72-E811-87C7-7CD30AD0A7AA.root"
    )
)

process.out = cms.OutputModule(
    "PoolOutputModule",fileName = cms.untracked.string("TauIDProducer_outputtest.root"),
    outputCommands = cms.untracked.vstring("drop *", "keep *_" + updatedTauName + "_*_*", "keep *_" + tauIDProducerName + "_*_*")
    )

process.p = cms.Path(process.rerunMvaIsolationSequence * process.slimmedTausNewID * getattr(process, tauIDProducerName))# put prior to n-tuplizer
process.ep = cms.EndPath(process.out)