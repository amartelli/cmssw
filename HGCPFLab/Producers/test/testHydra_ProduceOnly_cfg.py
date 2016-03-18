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
    #fileNames = cms.untracked.vstring('file:testHGCalLocalReco.root')
    fileNames = cms.untracked.vstring(
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_98.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_63.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_96.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_57.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_99.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_82.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_22.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_91.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_31.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_77.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_85.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_37.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_40.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_54.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_52.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_20.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_49.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_8.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_17.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_18.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_70.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_69.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_32.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_71.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_21.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_26.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_48.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_56.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_95.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_81.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_90.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_4.root',
        'root://cmsxrootd.fnal.gov//store/user/lmastrol/outputCrabTest/CRAB3_MC_ElectronGun_Feb16/160202_114208/0000/testHGCalLocalReco_28.root'
	)

)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('step3 nevts:-1'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands =  cms.untracked.vstring(),
    fileName = cms.untracked.string('file:step3_HydraOnly_Ele35.root'),
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
                                                              "keep *_particleFlowRecHitHGC*__*",
                                                              "keep *GenParticle*_genParticles_*_*",
                                                              "keep *_g4SimHits_HGCHits*_*",
                                                              "keep *SimTrack*_g4SimHits_*_*",
                                                              "keep *SimVertex*_g4SimHits_*_*",
                                                              "keep *_generalTracks_*_*")
