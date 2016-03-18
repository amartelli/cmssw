# Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions PH2_1K_FB_V6::All -n 10 --eventcontent RECOSIM -s RAW2DIGI,L1Reco,RECO --datatier GEN-SIM-RECO --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023HGCalMuon --geometry Extended2023HGCalMuon,Extended2023HGCalMuonReco --magField 38T_PostLS1 --filein file:step2.root --fileout file:step3.root --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#process.load('Configuration.Geometry.GeometryExtended2023Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023DevReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023_cff')
process.load('Configuration.Geometry.GeometryExtended2023Dev_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
'file:testHGCalLocalReco.root'
        )
    
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('step3 nevts:5'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands =  cms.untracked.vstring(),
    fileName = cms.untracked.string('file:./HydraProduced_calib.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.Hydra = cms.EDProducer('HydraProducer',
                               HGCalUncalibRecHitCollection=cms.VInputTag('HGCalUncalibRecHit:HGCEEUncalibRecHits',
	                                                                      'HGCalUncalibRecHit:HGCHEFUncalibRecHits'),
                               HGCRecHitCollection=cms.VInputTag("particleFlowRecHitHGC"),
                               GenParticleCollection=cms.InputTag("genParticles"),
                               RecTrackCollection=cms.InputTag("generalTracks"),
                               SimTrackCollection=cms.InputTag("g4SimHits"),
                               SimVertexCollection=cms.InputTag("g4SimHits"),
                               SimHitCollection = cms.VInputTag('g4SimHits:HGCHitsEE',
                                                                'g4SimHits:HGCHitsHEfront')
                                                                #'g4SimHits:HGCHitsHEback')
                               )


process.reconstruction += process.Hydra

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.Hydra)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.reconstruction_step,process.RECOSIMoutput_step)

# customisation of the process.

process.RECOSIMoutput.outputCommands =  cms.untracked.vstring("drop *",
                                                              "keep *_Hydra_*_*",
                                                              "keep *_HGCalUncalibRecHit_*_*", 
                                                              "keep *_particleFlowRecHitHGC*__*",
                                                              "keep *GenParticle*_genParticles_*_*",
                                                              "keep *_g4SimHits_HGCHits*_*",
                                                              "keep *SimTrack*_g4SimHits_*_*",
                                                              "keep *SimVertex*_g4SimHits_*_*",
                                                              "keep *_generalTracks_*_*")
