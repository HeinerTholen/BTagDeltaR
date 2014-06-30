import FWCore.ParameterSet.Config as cms

process = cms.Process("MERGE")

process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(
        'file:/nfs/dust/cms/user/tholenhe/samples/'
        'DoubleVtxEffTTdilep/TTdilep_presel_100_1_304.root'
    )
)

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands=cms.untracked.vstring(
        'drop *',
        'keep *_genParticles_*_*',
        'keep *Vertex*_*_*_*',
        'keep *_*DistInfo_*_*',
        'keep *_*Weight*_*_*',
    ),
    fileName=cms.untracked.string('merging.root'),
    SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('p')),
)
process.outPath = cms.EndPath(process.out)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.categories += (["BLA"])
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.load("MyUtility.EvtWeightPU.evtWeightPU_cff")

process.bToCharmDecayVertexMergedDistInfo = cms.EDProducer(
    'SecVtxInfoProducer',
    pv_src=cms.InputTag('offlinePrimaryVertices'),
    sv_src=cms.InputTag('bToCharmDecayVertexMerged'),
)

process.twoVtxFilter = cms.EDFilter("VertexCountFilter",
    src=cms.InputTag('inclusiveMergedVerticesFiltered'),
    minNumber=cms.uint32(2),
    maxNumber=cms.uint32(999)
)

process.p = cms.Path(
    process.twoVtxFilter *
    process.bToCharmDecayVertexMergedDistInfo *
    process.puWeight
)

