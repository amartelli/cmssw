# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: test_data_1025 --data -s NANO --data --eventcontent NANOAOD --datatier NANOAOD --filein /store/data/Run2018D/ParkingBPH1/MINIAOD/PromptReco-v2/000/321/833/00000/A8836E36-73AE-E811-AF6E-FA163E66D13C.root --no_exec -n 100 --conditions 102X_dataRun2_Prompt_v11 --era Run2_2018 --customise_commands=process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False))) --customise_commands=process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True) --customise_commands=process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")


runBToKstll = True
eleFinalState = True
lowPtEleFinalState = False
lowPtAndPfEleFinalState = True
kstarFinalState = False
useLostLeadMuonTracks = False
useLostSubLeadLepTracks = False ## cannot be true if lowPtEleFinalState == true
useLostChHadrTracks = True

import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('NANO',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('/store/data/Run2018D/ParkingBPH1/MINIAOD/PromptReco-v2/000/321/833/00000/A8836E36-73AE-E811-AF6E-FA163E66D13C.root'),
    fileNames = cms.untracked.vstring('/store/data/Run2018D/ParkingBPH5/MINIAOD/20Mar2019-v1/120000/0071842F-3A26-0D43-92F2-E2376273008E.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('test_data_10213 nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('test_data_10213_NANO.root'),
    outputCommands = process.NANOAODEventContent.outputCommands
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Prompt_v12', '')

# Path and EndPath definitions
process.nanoAOD_step = cms.Path(process.nanoSequence)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
process.schedule = cms.Schedule(process.nanoAOD_step,process.endjob_step,process.NANOAODoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeData 

#call to customisation function nanoAOD_customizeData imported from PhysicsTools.NanoAOD.nano_cff
process = nanoAOD_customizeData(process)

##the following should be enough                                                   
if runBToKstll:
    from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeBToKstll
    process = nanoAOD_customizeBToKstll(process)
if eleFinalState:
    from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeEleFinalState
    process = nanoAOD_customizeEleFinalState(process)
if lowPtEleFinalState:
    from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeLowPtEleFinalState
    process = nanoAOD_customizeLowPtEleFinalState(process)
if lowPtAndPfEleFinalState:
    from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeLowPtAndPfEleFinalState
    process = nanoAOD_customizeLowPtAndPfEleFinalState(process)
if kstarFinalState:
    from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeKstarFinalState
    process = nanoAOD_customizeKstarFinalState(process)
if useLostLeadMuonTracks:
    from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeLostLeadMuonTracks
    process = nanoAOD_customizeLostLeadMuonTracks(process)    
if useLostSubLeadLepTracks:
    from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeLostSubLeadLepTracks
    process = nanoAOD_customizeLostSubLeadLepTracks(process)
if useLostChHadrTracks:
    from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeLostChHadrTracks
    process = nanoAOD_customizeLostChHadrTracks(process)


# End of customisation functions

# Customisation from command line

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion 
