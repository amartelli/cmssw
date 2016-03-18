import FWCore.ParameterSet.Config as cms

process = cms.Process('TESTHYDRA')


process.load('Configuration.Geometry.GeometryExtended2023DevReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023Dev_cff')

#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/l/lmastrol/CMSSW_8_0_X_2016-01-26-1100/src/HGCPFLab/Producers/test/step3_HydraOnly_Ele35.root'))
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/l/lmastrol/LatestHydra/CMSSW_8_0_X_2016-02-07-2300/src/HGCPFLab/Producers/test/step3_HydraOnly_Photon35.root'))
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/l/lmastrol/LatestHydra/CMSSW_8_0_X_2016-02-07-2300/src/HGCPFLab/Producers/test/step3_HydraOnly_Mu100.root'))

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
#        'file:/tmp/amartell/HydraProduced_calib.root'
#        'file:/afs/cern.ch/work/l/lmastrol/public/xArabella/step3_HydraOnly_Ph35_HydraFix.root'
'file:HydraProduced_calib.root'
        ))



process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
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

process.HydraCaloCalibrator = cms.EDProducer("HydraCalibrator",
                                             HGCalUncalibRecHitCollection=cms.VInputTag('HGCalUncalibRecHit:HGCEEUncalibRecHits',
                                                                                        'HGCalUncalibRecHit:HGCHEFUncalibRecHits'),
                                             HydraTag=cms.InputTag("Hydra"),
                                             SplitRecHits=cms.bool(False),
                                             UseGenParticles=cms.bool(True),
                                             MinDebugEnergy=cms.untracked.double(30.)
                                        )


process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands =  cms.untracked.vstring('keep *'),
                                         fileName = cms.untracked.string('file:photons_Calib.root'),
#	fileName = cms.untracked.string('file:test_Calibrated.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
)

process.TFileService = cms.Service("TFileService",
#      fileName = cms.string("histoMu_calib.root"),
      fileName = cms.string("histoPho_calib.root"),
      closeFileFast = cms.untracked.bool(True)
  )
  
  
#process.testSequence = cms.Sequence(process.ExampleReader+process.FakeClusterGen+process.FakeClusterCaloFace)
#process.testSequence = cms.Sequence(process.FakeClusterCaloFace)
process.testSequence = cms.Sequence(process.HydraCaloCalibrator)
process.p  = cms.Path(process.testSequence)
process.ep = cms.EndPath(process.RECOSIMoutput)
