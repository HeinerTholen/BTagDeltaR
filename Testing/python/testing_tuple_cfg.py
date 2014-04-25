
###################################################################### init ###
import FWCore.ParameterSet.Config as cms

process = cms.Process("SecondaryVertex")

# load the full reconstraction configuration, to make sure we're getting all needed dependencies
#process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# for the conditions
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['startup']

# message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.categories += (["BLA"])
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))


########################################################### event selection ###
process.genPartFinalB = cms.EDProducer(
    "GenPartFinalB",
    src=cms.InputTag("genParticles"),
)

process.beePairs = cms.EDProducer(
    "CandidatePairProducer",
    src=cms.InputTag("genPartFinalB"),
)

cutstr_dr = "deltaR(daughterPtr(0).eta, daughterPtr(0).phi, daughterPtr(1).eta, daughterPtr(1).phi)"
process.closeBeePairs = cms.EDFilter(
    "CandViewSelector",
    src=cms.InputTag("beePairs"),
    cut=cms.string(cutstr_dr+" < 0.8"),
    filter=cms.bool(True),
)

################################################################# SV config ###

#process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("BTagDeltaR.My2ndVtxConfig.My2ndVtxConfig_cff")

# path
process.p = cms.Path(
    process.genPartFinalB *
    process.beePairs *
    process.closeBeePairs *
    process.my2ndVtxSequence
)


######################################################################## io ###
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM1.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM2.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM3.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM4.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM5.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM6.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM7.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM8.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM9.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM10.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM101.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM112.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM103.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM114.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM105.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM116.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM117.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM118.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM119.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM20.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM61.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM60.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM63.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM66.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM67.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM68.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM69.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM70.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM71.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM72.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM73.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM74.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM75.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM77.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM78.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM79.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM80.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM81.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM82.root",
      "file:/nfs/dust/cms/user/tholenhe/samples/Zbb_batch1/ZbbhadronicAODSIM83.root",
    )
)

process.out = cms.OutputModule("PoolOutputModule",
   outputCommands = cms.untracked.vstring(
       'keep *',
       'drop *_genPartFinalB_*_*',
       'drop *_beePairs_*_*',
       'drop *_closeBeePairs_*_*',
   ),
   fileName = cms.untracked.string(
       '/nfs/dust/cms/user/tholenhe/samples/recoPlusSV.root'
   ),
   SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
)

process.output = cms.EndPath(
   process.out
)
