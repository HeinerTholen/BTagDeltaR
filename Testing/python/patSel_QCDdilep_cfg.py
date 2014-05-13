import sys

import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('standard')
options.register('runOnMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "decide if run on MC or data")
options.register('outputFile', 'dijet.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, "name of output file")
#if( hasattr(sys, "argv") ):
#  options.parseArguments()


process = cms.Process( 'PAT' )


### ======================================================================== ###
###                                                                          ###
###                                 Constants                                ###
###                            (user job steering)                           ###
###                                                                          ###
### ======================================================================== ###


### Data or MC?
runOnMC = options.runOnMC

### Switch on/off selection steps

# Step 1
# (trigger selection: QuadJet50_Jet40)
useTrigger      = True
# Step 2
# (good vertex selection)
useGoodVertex   = True
# Step 3a
# (6 jets: pt > 30 GeV & |eta| < 2.4)
# Step 3b
# (4 jets: pt > 60 GeV, 5 jets: pt > 50 GeV, 6 jets: pt > 30 GeV, for all jets: |eta| < 2.4)
# (the cuts for the 4 and 5 leading jets are configurable via jetCutHard / jetCutMedium respectivly)

### Reference selection

from TopQuarkAnalysis.Configuration.patRefSel_refAllJets import *
# Muons
#muonCut = ''
# Electrons
#electronCut = ''
# Jets
#jetCut          = ''
#veryLooseJetCut = 'pt > 35.'
#looseJetCut     = 'pt > 50.'
#tightJetCut     = 'pt > 60.'

# Trigger and trigger object
#triggerSelectionData       = ''
#triggerObjectSelectionData = ''
#triggerSelectionMC       = ''
#triggerObjectSelectionMC = ''

### Particle flow

postfix = 'PF'

# subtract charged hadronic pile-up particles (from wrong PVs)
# effects also JECs
usePFnoPU       = False # before any top projection
usePfIsoLessCHS = True # switch to new PF isolation with L1Fastjet CHS

# other switches for PF top projections (default: all 'True')
useNoMuon     = True # before electron top projection
useNoElectron = True # before jet top projection
useNoJet      = True # before tau top projection
useNoTau      = True # before MET top projection

# cuts used in top projections
from TopQuarkAnalysis.Configuration.patRefSel_PF2PAT import *
# vertices
#pfVertices = 'goodOfflinePrimaryVertices'
#pfD0Cut   = 0.2
#pfDzCut   = 0.5
# muons
#pfMuonSelectionCut = 'pt > 5.'
useMuonCutBasePF = True # use minimal (veto) muon selection cut on top of 'pfMuonSelectionCut'
#muonCutPF = 'pt > 10. && abs(eta) < 2.5'
#pfMuonIsoConeR03 = False
#pfMuonCombIsoCut = 0.2
# electrons
#pfElectronSelectionCut  = 'pt > 5. && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfLostHits < 2'
useElectronCutBasePF  = True # use minimal (veto) electron selection cut on top of 'pfElectronSelectionCut'
#electronCutPF = 'pt > 20. && abs(eta) < 2.5'
#pfElectronIsoConeR03 = True
#pfElectronCombIsoCut  = 0.2

### JEC levels

# levels to be accessible from the jets
# jets are corrected to L3Absolute (MC), L2L3Residual (data) automatically, if enabled here
# and remain uncorrected, if none of these levels is enabled here
useL1FastJet    = True  # needs useL1Offset being off, error otherwise
useL1Offset     = False # needs useL1FastJet being off, error otherwise; not available from current GT!!!
useL2Relative   = True
useL3Absolute   = True
useL2L3Residual = True
useL5Flavor     = False
useL7Parton     = False

typeIMetCorrections = True

### Input

# list of input files
inputFiles = [
    "file:/nfs/dust/cms/user/tholenhe/samples/QCD_test.root"
]   # overwritten, if "useRelVals" is 'True'


# maximum number of events
maxEvents = options.maxEvents

### Conditions

# GlobalTags (w/o suffix '::All')
globalTagData = 'FT53_V21A_AN6::All'
globalTagMC   = 'START53_V27::All'

### Output

# output file
outputFile = options.outputFile

# event frequency of Fwk report
fwkReportEvery = 100

# switch for 'TrigReport'/'TimeReport' at job end
wantSummary = True


###                              End of constants                            ###
###                                                                          ###
### ======================================================================== ###


###
### Basic configuration
###

process.load( "TopQuarkAnalysis.Configuration.patRefSel_basics_cff" )
process.MessageLogger.cerr.FwkReport.reportEvery = fwkReportEvery
process.options.wantSummary = wantSummary
if runOnMC:
    process.GlobalTag.globaltag = globalTagMC
else:
    process.GlobalTag.globaltag = globalTagData


###
### Input configuration
###

process.load( "TopQuarkAnalysis.Configuration.patRefSel_inputModule_cfi" )
process.source.fileNames = inputFiles
process.maxEvents.input  = maxEvents


###
### Output configuration
###

process.load( "TopQuarkAnalysis.Configuration.patRefSel_outputModule_cff" )
# output file name
process.out.fileName = outputFile
# event content
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out.outputCommands += patEventContent
# clear event selection
process.out.SelectEvents.SelectEvents = []


###
### Cleaning and trigger selection configuration
###

### Trigger selection
triggerSelection = "HLT_DiPFJetAve80_v*"
from TopQuarkAnalysis.Configuration.patRefSel_triggerSelection_cff import triggerResults
process.step1 = triggerResults.clone(
    triggerConditions = [ triggerSelection ]
)

### Good vertex selection
process.load( "TopQuarkAnalysis.Configuration.patRefSel_goodVertex_cfi" )
process.step2 = process.goodOfflinePrimaryVertices.clone( filter = True )

### Event cleaning
process.load( 'TopQuarkAnalysis.Configuration.patRefSel_eventCleaning_cff' )
process.trackingFailureFilter.VertexSource = cms.InputTag( pfVertices )


###
### PAT/PF2PAT configuration
###

process.load( "PhysicsTools.PatAlgos.patSequences_cff" )

### Check JECs

# JEC set
jecSet = 'AK5PF'
if usePFnoPU:
    jecSet += 'chs'

# JEC levels
if useL1FastJet and useL1Offset:
    sys.exit( 'ERROR: switch off either "L1FastJet" or "L1Offset"' )
jecLevels = []
if useL1FastJet:
    jecLevels.append( 'L1FastJet' )
if useL1Offset:
    jecLevels.append( 'L1Offset' )
if useL2Relative:
    jecLevels.append( 'L2Relative' )
if useL3Absolute:
    jecLevels.append( 'L3Absolute' )
if useL2L3Residual and not runOnMC:
    jecLevels.append( 'L2L3Residual' )
if useL5Flavor:
    jecLevels.append( 'L5Flavor' )
if useL7Parton:
    jecLevels.append( 'L7Parton' )

### Switch configuration
from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
usePF2PAT( process
         , runPF2PAT           = True
         , runOnMC             = runOnMC
         , jetAlgo             = jetAlgo
         , postfix             = postfix
         , jetCorrections      = ( jecSet
                                 , jecLevels
                                 )
         , typeIMetCorrections = typeIMetCorrections
         , pvCollection        = cms.InputTag( pfVertices )
         )

if useMuonCutBasePF:
    pfMuonSelectionCut += ' && %s'%( muonCutPF )
if useElectronCutBasePF:
    from TopQuarkAnalysis.Configuration.patRefSel_pfIdentifiedElectrons_cfi import pfIdentifiedElectrons
    setattr( process, 'pfIdentifiedElectrons' + postfix, pfIdentifiedElectrons )
    getattr( process, 'pfIdentifiedElectrons' + postfix ).src = cms.InputTag( 'pfElectronsFromVertex' + postfix )
    getattr( process, 'pfSelectedElectrons'   + postfix ).src = cms.InputTag( 'pfIdentifiedElectrons' + postfix )
    getattr( process, 'patPF2PATSequence' + postfix ).replace( getattr( process, 'pfSelectedElectrons' + postfix )
                                                           , getattr( process, 'pfIdentifiedElectrons' + postfix ) + getattr( process, 'pfSelectedElectrons' + postfix )
                                                           )
    pfElectronSelectionCut += ' && %s'%( electronCutPF )

getattr( process, 'pfNoPileUp'   + postfix ).enable = usePFnoPU
getattr( process, 'pfNoMuon'     + postfix ).enable = useNoMuon
getattr( process, 'pfNoElectron' + postfix ).enable = useNoElectron
getattr( process, 'pfNoJet'      + postfix ).enable = useNoJet
getattr( process, 'pfNoTau'      + postfix ).enable = useNoTau

if useL1FastJet:
    getattr( process, 'pfPileUpIso' + postfix ).checkClosestZVertex = usePfIsoLessCHS

getattr( process, 'pfMuonsFromVertex' + postfix ).d0Cut    = pfD0Cut
getattr( process, 'pfMuonsFromVertex' + postfix ).dzCut    = pfDzCut
getattr( process, 'pfSelectedMuons'   + postfix ).cut = pfMuonSelectionCut
getattr( process, 'pfIsolatedMuons'   + postfix ).doDeltaBetaCorrection = True
getattr( process, 'pfIsolatedMuons'   + postfix ).deltaBetaFactor       = -0.5
getattr( process, 'pfIsolatedMuons'   + postfix ).isolationCut          = pfMuonCombIsoCut
if pfMuonIsoConeR03:
    getattr( process, 'pfIsolatedMuons' + postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'muPFIsoValueCharged03' + postfix )
                                                                                                )
    getattr( process, 'pfIsolatedMuons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'muPFIsoValuePU03' + postfix )
    getattr( process, 'pfIsolatedMuons' + postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'muPFIsoValueNeutral03' + postfix )
                                                                                                , cms.InputTag( 'muPFIsoValueGamma03' + postfix )
                                                                                                )
    getattr( process, 'pfMuons' + postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'muPFIsoValueCharged03' + postfix )
                                                                                        )
    getattr( process, 'pfMuons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'muPFIsoValuePU03' + postfix )
    getattr( process, 'pfMuons' + postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'muPFIsoValueNeutral03' + postfix )
                                                                                        , cms.InputTag( 'muPFIsoValueGamma03' + postfix )
                                                                                        )
    getattr( process, 'patMuons' + postfix ).isolationValues.pfNeutralHadrons   = cms.InputTag( 'muPFIsoValueNeutral03' + postfix )
    getattr( process, 'patMuons' + postfix ).isolationValues.pfChargedAll       = cms.InputTag( 'muPFIsoValueChargedAll03' + postfix )
    getattr( process, 'patMuons' + postfix ).isolationValues.pfPUChargedHadrons = cms.InputTag( 'muPFIsoValuePU03' + postfix )
    getattr( process, 'patMuons' + postfix ).isolationValues.pfPhotons          = cms.InputTag( 'muPFIsoValueGamma03' + postfix )
getattr( process, 'patMuons' + postfix ).isolationValues.pfChargedHadrons   = cms.InputTag( 'muPFIsoValueCharged03' + postfix )
getattr( process, 'pfElectronsFromVertex' + postfix ).d0Cut    = pfD0Cut
getattr( process, 'pfElectronsFromVertex' + postfix ).dzCut    = pfDzCut
getattr( process, 'pfSelectedElectrons'   + postfix ).cut = pfElectronSelectionCut
getattr( process, 'pfIsolatedElectrons'   + postfix ).doDeltaBetaCorrection = True # applies EA corrections here!
getattr( process, 'pfIsolatedElectrons'   + postfix ).deltaBetaFactor       = -1.0
getattr( process, 'pfIsolatedElectrons'   + postfix ).isolationCut          = pfElectronCombIsoCut
if pfElectronIsoConeR03:
    from EgammaAnalysis.ElectronTools.electronIsolatorFromEffectiveArea_cfi import elPFIsoValueEA03
    setattr( process, 'elPFIsoValueEA03' + postfix, elPFIsoValueEA03 )
    getattr( process, 'elPFIsoValueEA03' + postfix ).pfElectrons = cms.InputTag( 'pfSelectedElectrons' + postfix )
    getattr( process, 'patPF2PATSequence' + postfix ).replace( getattr( process, 'pfSelectedElectrons' + postfix )
                                                           , getattr( process, 'pfSelectedElectrons' + postfix ) + getattr( process, 'elPFIsoValueEA03' + postfix )
                                                           )
    getattr( process, 'pfIsolatedElectrons' + postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'elPFIsoValueCharged03PFId' + postfix )
                                                                                                     )
    #getattr( process, 'pfIsolatedElectrons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'elPFIsoValuePU03PFId' + postfix )
    getattr( process, 'pfIsolatedElectrons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'elPFIsoValueEA03' + postfix ) # EA corrections
    getattr( process, 'pfIsolatedElectrons' + postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'elPFIsoValueNeutral03PFId' + postfix )
                                                                                                    , cms.InputTag( 'elPFIsoValueGamma03PFId'   + postfix )
                                                                                                    )
    getattr( process, 'pfElectrons' + postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'elPFIsoValueCharged03PFId' + postfix )
                                                                                             )
    #getattr( process, 'pfElectrons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'elPFIsoValuePU03PFId' + postfix )
    getattr( process, 'pfElectrons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'elPFIsoValueEA03' + postfix ) # EA corrections
    getattr( process, 'pfElectrons' + postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'elPFIsoValueNeutral03PFId' + postfix )
                                                                                            , cms.InputTag( 'elPFIsoValueGamma03PFId'   + postfix )
                                                                                            )
    getattr( process, 'patElectrons' + postfix ).isolationValues.pfNeutralHadrons   = cms.InputTag( 'elPFIsoValueNeutral03PFId' + postfix )
    getattr( process, 'patElectrons' + postfix ).isolationValues.pfChargedAll       = cms.InputTag( 'elPFIsoValueChargedAll03PFId' + postfix )
    getattr( process, 'patElectrons' + postfix ).isolationValues.pfPUChargedHadrons = cms.InputTag( 'elPFIsoValuePU03PFId' + postfix )
    getattr( process, 'patElectrons' + postfix ).isolationValues.pfPhotons          = cms.InputTag( 'elPFIsoValueGamma03PFId' + postfix )
    getattr( process, 'patElectrons' + postfix ).isolationValues.pfChargedHadrons   = cms.InputTag( 'elPFIsoValueCharged03PFId' + postfix )
    getattr( process, 'patElectrons' + postfix ).isolationValues.user               = cms.VInputTag( cms.InputTag( "elPFIsoValueEA03%s"%( postfix ) ) )
else:
    from EgammaAnalysis.ElectronTools.electronIsolatorFromEffectiveArea_cfi import elPFIsoValueEA04
    setattr( process, 'elPFIsoValueEA04' + postfix, elPFIsoValueEA04 )
    getattr( process, 'elPFIsoValueEA04' + postfix ).pfElectrons = cms.InputTag( 'pfSelectedElectrons' + postfix )
    getattr( process, 'patPF2PATSequence' + postfix ).replace( getattr( process, 'pfSelectedElectrons' + postfix )
                                                           , getattr( process, 'pfSelectedElectrons' + postfix ) + getattr( process, 'elPFIsoValueEA04' + postfix )
                                                           )
    getattr( process, 'pfIsolatedElectrons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'elPFIsoValueEA04' + postfix ) # EA corrections
    getattr( process, 'pfElectrons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'elPFIsoValueEA04' + postfix ) # EA corrections
    getattr( process, 'patElectrons' + postfix ).isolationValues.user = cms.VInputTag( cms.InputTag( "elPFIsoValueEA04%s"%( postfix ) ) )


from PhysicsTools.PatAlgos.tools.coreTools import *

from TopQuarkAnalysis.Configuration.patRefSel_refAllJets_cfi import *

# remove MC matching, object cleaning, photons and taus

if not runOnMC:
    runOnData( process
           , names = [ 'PFAll' ]
           , postfix = postfix
           )
removeSpecificPATObjects( process
                        , names = [ 'Photons', 'Taus' ]
                        , postfix = postfix
                        ) # includes 'removeCleaning'

# additional event content has to be (re-)added _after_ the call to 'removeCleaning()':
process.out.outputCommands += [ 'keep edmTriggerResults_*_*_*'
                              , 'keep *_hltTriggerSummaryAOD_*_*'
                              # vertices and beam spot
                              , 'keep *_offlineBeamSpot_*_*'
                              , 'keep *_offlinePrimaryVertices*_*_*'
                              , 'keep *_goodOfflinePrimaryVertices*_*_*'
                              ]
if runOnMC:
    process.out.outputCommands += [ 'keep GenEventInfoProduct_*_*_*'
                                , 'keep recoGenParticles_*_*_*'
                                , 'keep *_addPileupInfo_*_*'
                                ]


###
### Scheduling
###

# MVA electron ID

process.load( "EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi" )
process.eidMVASequence = cms.Sequence(
  process.mvaTrigV0
+ process.mvaNonTrigV0
)

### Heiners preselection
process.jetPreselector = cms.EDFilter("CandViewSelector",
    src=cms.InputTag("ak5PFJets"),
    cut=cms.string('pt > 20.')
)
process.jetPreFilter = cms.EDFilter(
    "PATCandViewCountFilter",
    src=cms.InputTag('jetPreselector'),
    minNumber=cms.uint32(1),
    maxNumber=cms.uint32(2)
)
process.preselection = cms.Sequence(
    process.jetPreselector *
    process.jetPreFilter
)

### Heiners dijet selection
process.selectedPatJetsPF.cut = 'pt > 15.'

process.twoJetFilter = cms.EDFilter(
    "PATCandViewCountFilter",
    src=cms.InputTag('selectedPatJetsPF'),
    minNumber=cms.uint32(2),
    maxNumber=cms.uint32(2)
)

cutstr_avept = "(daughterPtr(0).pt+daughterPtr(1).pt)/2."
process.averageJetPt = cms.EDFilter("CandViewSelector",
    src=cms.InputTag("dijetPair"),
    cut=cms.string(cutstr_avept + ' > 105.'),
    filter=cms.bool(True)
)

process.dijetPair = cms.EDProducer(
    "CandidatePairProducer",
    src=cms.InputTag("selectedPatJetsPF")
)
cutstr_dr = "deltaR(" \
            "0., daughterPtr(0).phi, " \
            "0., daughterPtr(1).phi)" # only delta phi
process.dijetDeltaRFilter = cms.EDFilter("CandViewSelector",
    src=cms.InputTag("dijetPair"),
    cut=cms.string(cutstr_dr + ' > 2.7'),
    filter=cms.bool(True)
)

process.dijetSequence = cms.Sequence(
    process.twoJetFilter *
    process.dijetPair *
    process.averageJetPt *
    process.dijetDeltaRFilter
)

### Heiners additional event content
process.out.outputCommands += ['keep *_*_*_RECO']

### Heiners additional vertexing
import BTagDeltaR.My2ndVtxConfig.My2ndVtxConfig_cff as my_sv
process.extend(my_sv)
process.out.outputCommands += my_sv.my2ndVtxEventCont


# The paths

process.p = cms.Path()
process.p += process.step1
process.p += process.preselection
process.p += process.goodOfflinePrimaryVertices
process.p += process.step2
process.p += process.eventCleaning
if runOnMC:
    process.p += process.eventCleaningMC
else:
    process.p += process.eventCleaningData
process.p += process.eidMVASequence
process.p += getattr( process, 'patPF2PATSequence' + postfix )
process.p += process.dijetSequence
process.p += process.my2ndVtxSequence
process.out.SelectEvents.SelectEvents.append( 'p' )
