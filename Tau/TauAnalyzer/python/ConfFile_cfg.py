import FWCore.ParameterSet.Config as cms

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

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/ceph/sbrommer/tau_id/dy_testsample/FEB3954C-4942-E811-8A09-008CFAC91A38.root'
                )
                            )

process.tauEDAnalyzer = cms.EDAnalyzer('TauAnalyzer'
                              )

process.p = cms.Path(process.tauEDAnalyzer)
