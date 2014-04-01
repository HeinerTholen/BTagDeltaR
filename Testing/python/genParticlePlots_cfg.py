
import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(
#    "file:/afs/desy.de/user/t/tholenhe/xxl-af-cms/samples/syncExercise53.root"
    "file:/afs/desy.de/user/t/tholenhe/xxl-af-cms/samples/Zbb_batch1/ZbbhadronicAODSIM1.root",
    "file:/afs/desy.de/user/t/tholenhe/xxl-af-cms/samples/Zbb_batch1/ZbbhadronicAODSIM2.root"
    )
)

process.TFileService = cms.Service(
    "TFileService",
    fileName=cms.string('FILESERVICE.root')
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.categories += (["BLA"])
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))


process.genPartFromZ = cms.EDProducer(
    "GenPartFromZDecay",
    src=cms.InputTag("genParticles"),
    inverted=cms.bool(False),
)

process.genPartFinalB = cms.EDProducer(
    "GenPartFinalB",
    src=cms.InputTag("genPartFromZ"),
)

process.beePairs = cms.EDProducer(
    "CandidatePairProducer",
    src=cms.InputTag("genPartFinalB"),
)

process.finalBeeMult = cms.EDAnalyzer(
    "MultiplicityHisto",
    src=cms.InputTag("genPartFinalB"),
    nbins=cms.untracked.int32(10),
)

import MyUtility.PythonUtil.utility as utility
process.finalBeePt = utility.make_histo_analyzer(
    "genPartFinalB",
    "pt",
    50,
)

process.finalBeeEta = utility.make_histo_analyzer(
    "genPartFinalB",
    "eta",
    50,
)

process.finalBeePhi = utility.make_histo_analyzer(
    "genPartFinalB",
    "phi",
    50,
)

cutstr_dr = "deltaR(daughterPtr(0).eta, daughterPtr(0).phi, daughterPtr(1).eta, daughterPtr(1).phi)"
process.beePairDeltaR = utility.make_histo_analyzer(
    "beePairs",
    cutstr_dr,
    500,
)






process.zeeDaughters = cms.EDFilter(
    "CandViewSelector",
    src=cms.InputTag("genPartFromZ"),
    cut=cms.string("abs(pdgId) == 5 && mother.pdgId == 23"),
)

process.zeeDauPairs = cms.EDProducer(
    "CandidatePairProducer",
    src=cms.InputTag("zeeDaughters"),
)

process.boostedZeeDauPairs = cms.EDFilter(
    "CandViewSelector",
    src=cms.InputTag("zeeDauPairs"),
    cut=cms.string("pt > 100."),
)

process.zeeDauDeltaR = utility.make_histo_analyzer(
    "boostedZeeDauPairs",
    cutstr_dr,
    50,
)

process.zeeDauMass = utility.make_histo_analyzer(
    "boostedZeeDauPairs",
    "mass",
    50
)

process.zeeDauPt = utility.make_histo_analyzer(
    "boostedZeeDauPairs",
    "pt",
    50
)



process.p = cms.Path(
    process.genPartFromZ
    * process.genPartFinalB
    * process.beePairs
    * process.finalBeeMult
    * process.finalBeePt
    * process.finalBeeEta
    * process.finalBeePhi
    * process.beePairDeltaR

    * process.zeeDaughters
    * process.zeeDauPairs
    * process.boostedZeeDauPairs
    * process.zeeDauDeltaR
    * process.zeeDauMass
    * process.zeeDauPt
)