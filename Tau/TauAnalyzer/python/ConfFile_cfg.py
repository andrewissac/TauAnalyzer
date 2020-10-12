import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("TAU")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('TrackNumber')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# process.source = cms.Source("PoolSource",
#                                 # replace 'myfile.root' with the source file you want to use
#                                 fileNames = cms.untracked.vstring(
#             #'file:/ceph/sbrommer/tau_id/dy_testsample/FEB3954C-4942-E811-8A09-008CFAC91A38.root'
#             #'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/6A62DE01-E641-E811-AFF1-008CFAC94038.root'
#                 )
#                             )

arguments_jobID = int(sys.argv[2]) 
inputfile = str(sys.argv[3]) # this takes the filename of the root file passed in from job.sh from running cmsRun mypath/ConfFile_cfg.py myfile.root (in this case: sys.argv[3] = myfile.root )
dataset = str(sys.argv[4])

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(inputfile)
                            )

process.tauEDAnalyzer = cms.EDAnalyzer('TauAnalyzer')

process.TFileService = cms.Service( "TFileService", fileName=cms.string("output_{:04d}_{}.root".format(arguments_jobID, dataset)))
#process.TFileService = cms.Service( "TFileService", fileName=cms.string("test.root"))

process.p = cms.Path(process.tauEDAnalyzer)
