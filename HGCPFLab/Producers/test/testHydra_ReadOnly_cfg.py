import FWCore.ParameterSet.Config as cms

process = cms.Process('TESTHYDRA')

#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/l/lmastrol/CMSSW_8_0_X_2016-01-26-1100/src/HGCPFLab/Producers/test/step3_HydraOnly_Ele35.root'))
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/l/lmastrol/LatestHydra/CMSSW_8_0_X_2016-02-07-2300/src/HGCPFLab/Producers/test/step3_HydraOnly_Photon35.root'))
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/l/lmastrol/LatestHydra/CMSSW_8_0_X_2016-02-07-2300/src/HGCPFLab/Producers/test/step3_HydraOnly_Mu100.root'))
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.ExampleReader = cms.EDProducer("ExampleHydraPFProducer",HydraTag=cms.InputTag("Hydra"))

process.FakeClusterGen = cms.EDProducer("HydraFakeClusterBuilder",HydraTag=cms.InputTag("Hydra"),
                                        SplitRecHits=cms.bool(False),
                                        UseGenParticles=cms.bool(True),
                                        MinDebugEnergy=cms.untracked.double(30.)
                                       )

process.FakeClusterCaloFace = cms.EDProducer("HydraFakeClusterBuilder",HydraTag=cms.InputTag("Hydra"),
                                             SplitRecHits=cms.bool(False),
                                             UseGenParticles=cms.bool(False),
                                             MinDebugEnergy=cms.untracked.double(30.)
                                        )
process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands =  cms.untracked.vstring('keep *'),
    #fileName = cms.untracked.string('file:step3b_hydra_photons_1ev.root'),
	fileName = cms.untracked.string('file:testMu100.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("histoMu.root"),
      closeFileFast = cms.untracked.bool(True)
  )
  
  
#process.testSequence = cms.Sequence(process.ExampleReader+process.FakeClusterGen+process.FakeClusterCaloFace)
process.testSequence = cms.Sequence(process.FakeClusterCaloFace)
process.p  = cms.Path(process.testSequence)
process.ep = cms.EndPath(process.RECOSIMoutput)
