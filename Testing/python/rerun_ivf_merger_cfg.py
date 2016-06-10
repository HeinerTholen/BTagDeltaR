import FWCore.ParameterSet.Config as cms
import glob



process = cms.Process('TEST')

process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(
        map(lambda s: 'file:%s' % s, glob.glob(
            "/nfs/dust/cms/user/tholenhe/samples/DoubleVtxEffTTdilep/"
        )[:10])
    )
)

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(
        'keep *',
    ),
    fileName = cms.untracked.string(
        '/nfs/dust/cms/user/tholenhe/samples/test_rerun_vtx_merger.root'
    ),
    #SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
)

from Configuration.StandardSequences.Reconstruction_cff import *
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *
inclusiveMergedVerticesFiltered.vertexFilter.maxDeltaRToVtxMomentum = cms.double(0.7)
inclusiveMergedVerticesFiltered.vertexFilter.maxDeltaRToJetAxis = cms.double(0.7)


process.p = cms.Path(

)
process.p += process.my2ndVtxSequence
