import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/ceph/sbrommer/tau_id/dy_testsample/FEB3954C-4942-E811-8A09-008CFAC91A38.root'
                )
                            )

process.demo = cms.EDAnalyzer('TauAnalyzer'
                              )

process.p = cms.Path(process.demo)
