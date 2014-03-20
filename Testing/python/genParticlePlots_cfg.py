
import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(
    "file:/afs/desy.de/user/t/tholenhe/xxl-af-cms/samples/syncExercise53.root"
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


process.genPartFinalB = cms.EDProducer(
    "GenPartFinalB",
    src=cms.InputTag("genParticles"),
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
    50, 0., 500.
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

process.p = cms.Path(
    process.genPartFinalB
    * process.finalBeeMult
    * process.finalBeePt
    * process.finalBeeEta
    * process.finalBeePhi
)