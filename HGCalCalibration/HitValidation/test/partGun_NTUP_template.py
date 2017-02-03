import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Demo")
process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


import FWCore.Utilities.FileUtils as FileUtils
readFiles = cms.untracked.vstring()
readFiles.extend(FileUtils.loadListFromFile ('INPUTFILELIST') )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
                            fileNames = readFiles,
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                            )



process.ana = cms.EDAnalyzer('HGCalHitCalibration',
                             detector = cms.string("all"),
                             rawRecHits = cms.bool(True)
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("file:OUTFILE.root")

                                   )

#process.hgcalLayerClusters.minClusters = cms.uint32(0)
#process.hgcalLayerClusters.realSpaceCone = cms.bool(True)

process.p = cms.Path(process.ana)
#process.p = cms.Path(process.hgcalLayerClusters+process.ana)
