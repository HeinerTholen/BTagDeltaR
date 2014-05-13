
from Configuration.StandardSequences.Reconstruction_cff import *
from Configuration.StandardSequences.MagneticField_cff import *
from Configuration.Geometry.GeometryIdeal_cff import *

# ak5pf
MyAk5PFJetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
    j2tParametersVX,
    jets = cms.InputTag("ak5PFJets")
)

MyAk5PFImpactParameterPFTagInfos = impactParameterTagInfos.clone(
    jetTracks = "MyAk5PFJetTracksAssociatorAtVertex"
)

MyAk5PFSecondaryVertexTagInfos = secondaryVertexTagInfos.clone(
    trackIPTagInfos = cms.InputTag("MyAk5PFImpactParameterPFTagInfos")
)
MyAk5PFSecondaryVertexTagInfos.trackSelection.jetDeltaRMax = 0.5
MyAk5PFSecondaryVertexTagInfos.vertexCuts.maxDeltaRToJetAxis = 0.5

# ak7pf
MyAk7PFJetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
    j2tParametersVX,
    jets = cms.InputTag("ak7PFJets")
)

MyAk7PFImpactParameterPFTagInfos = impactParameterTagInfos.clone(
    jetTracks = "MyAk7PFJetTracksAssociatorAtVertex"
)

MyAk7PFSecondaryVertexTagInfos = secondaryVertexTagInfos.clone(
    trackIPTagInfos = cms.InputTag("MyAk7PFImpactParameterPFTagInfos")
)
MyAk7PFSecondaryVertexTagInfos.trackSelection.jetDeltaRMax = 0.7
MyAk7PFSecondaryVertexTagInfos.vertexCuts.maxDeltaRToJetAxis = 0.7

# ivf
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *
inclusiveMergedVerticesFiltered.vertexFilter.maxDeltaRToVtxMomentum = cms.double(0.7)
inclusiveMergedVerticesFiltered.vertexFilter.maxDeltaRToJetAxis = cms.double(0.7)


MyAk5PFInclusiveSecondaryVertexFinderTagInfos = inclusiveSecondaryVertexFinderTagInfos.clone(
    extSVDeltaRToJet = cms.double(0.5),
    trackIPTagInfos = cms.InputTag("MyAk5PFImpactParameterPFTagInfos")
)
MyAk5PFInclusiveSecondaryVertexFinderTagInfos.vertexCuts.maxDeltaRToJetAxis = 0.5
MyAk5PFInclusiveSecondaryVertexFinderTagInfos.trackSelection.jetDeltaRMax = 0.5

MyAk7PFInclusiveSecondaryVertexFinderTagInfos = inclusiveSecondaryVertexFinderTagInfos.clone(
    extSVDeltaRToJet = cms.double(0.7),
    trackIPTagInfos = cms.InputTag("MyAk7PFImpactParameterPFTagInfos")
)
MyAk7PFInclusiveSecondaryVertexFinderTagInfos.vertexCuts.maxDeltaRToJetAxis = 0.7
MyAk7PFInclusiveSecondaryVertexFinderTagInfos.trackSelection.jetDeltaRMax = 0.7


my2ndVtxSequence = cms.Sequence(
    inclusiveVertexing *                    # ivf
    inclusiveMergedVerticesFiltered *       # ivf
    bToCharmDecayVertexMerged *             # ivf
    MyAk5PFJetTracksAssociatorAtVertex *
    MyAk5PFImpactParameterPFTagInfos *
    MyAk5PFSecondaryVertexTagInfos *
    MyAk7PFJetTracksAssociatorAtVertex *
    MyAk7PFImpactParameterPFTagInfos *
    MyAk7PFSecondaryVertexTagInfos *
    MyAk5PFInclusiveSecondaryVertexFinderTagInfos *
    MyAk7PFInclusiveSecondaryVertexFinderTagInfos
)

my2ndVtxEventCont = cms.untracked.vstring(
       'keep *_MyAk*_*_*',
       'keep *_bToCharmDecayVertexMerged_*_*',
       'keep *_inclusiveMergedVertices*_*_*',
       'keep *_inclusiveMergedVertices*_*_*',
       'keep *_inclusiveVertexFinder_*_*',
       'keep *_trackVertexArbitrator_*_*',
       'keep *_vertexMerger_*_*',
)