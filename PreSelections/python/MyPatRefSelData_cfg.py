import sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing


options = VarParsing ('standard')
options.register('runOnMC', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "decide if run on MC or data")
options.register('outputFile',
                 #'/nfs/dust/cms/user/tholenhe/samples/TTdilep_presel.root',
                 'TTdilep_presel.root',
                 VarParsing.multiplicity.singleton, VarParsing.varType.string, "name of output file")
#if( hasattr(sys, "argv") ):
    #options.parseArguments()


process = cms.Process( 'PAT' )


### ======================================================================== ###
###                                                                          ###
###                                 Constants                                ###
###                            (user job steering)                           ###
###                                                                          ###
### ======================================================================== ###


### Data or MC?
runOnMC = options.runOnMC

from TopQuarkAnalysis.Configuration.patRefSel_refMuJets import *

# Trigger and trigger object
#triggerSelectionData       = ''
#triggerObjectSelectionData = ''
#triggerSelectionMC       = ''
#triggerObjectSelectionMC = ''

### Particle flow

postfix = 'PF'

# subtract charged hadronic pile-up particles (from wrong PVs)
# effects also JECs
usePFnoPU       = True # before any top projection
usePfIsoLessCHS = True # switch to new PF isolation with L1Fastjet CHS

# other switches for PF top projections (default: all 'True')
useNoMuon     = True # before electron top projection
useNoElectron = True # before jet top projection
useNoJet      = True # before tau top projection
useNoTau      = True # before MET top projection

# cuts used in top projections
# vertices
#pfVertices  = 'goodOfflinePrimaryVertices'
#pfD0Cut     = 0.2
#pfDzCut     = 0.5
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
#pfElectronCombIsoCut = 0.15

### JEC levels

# levels to be accessible from the jets
# jets are corrected to L3Absolute (MC), L2L3Residual (data) automatically, if enabled here
# and remain uncorrected, if none of these levels is enabled here
useL1FastJet    = True  # needs useL1Offset being off, error otherwise
useL1Offset     = False # needs useL1FastJet being off, error otherwise; not available from current GT!!!
useL2Relative   = True
useL3Absolute   = True
useL2L3Residual = True  # takes effect only on data
useL5Flavor     = False
useL7Parton     = False

typeIMetCorrections = True

### Input

# list of input files
useRelVals = True # if 'Falsed', "inputFiles" is used
inputFiles = ["file:/nfs/dust/cms/user/tholenhe/samples/TTdilep.root"] # overwritten, if "useRelVals" is 'True'

# maximum number of events
maxEvents = options.maxEvents

### Conditions

# GlobalTags
globalTagData = 'FT53_V21A_AN6::All'
globalTagMC   = 'START53_V27::All'

### Output

# output file
outputFile = options.outputFile

# event frequency of Fwk report
fwkReportEvery = 10

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
triggerSelection = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v* OR HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*"
from TopQuarkAnalysis.Configuration.patRefSel_triggerSelection_cff import triggerResults
process.step0a = triggerResults.clone(
  triggerConditions = [ triggerSelection ]
)

### Good vertex selection
process.load( "TopQuarkAnalysis.Configuration.patRefSel_goodVertex_cfi" )
process.step0b = process.goodOfflinePrimaryVertices.clone( filter = True )

### Event cleaning
process.load( 'TopQuarkAnalysis.Configuration.patRefSel_eventCleaning_cff' )
process.trackingFailureFilter.VertexSource = cms.InputTag( pfVertices )
process.step0c = process.eventCleaning
if runOnMC:
  process.step0c += process.eventCleaningMC
else:
  process.step0c += process.eventCleaningData


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

getattr( process, 'pfMuonsFromVertex' + postfix ).d0Cut = pfD0Cut
getattr( process, 'pfMuonsFromVertex' + postfix ).dzCut = pfDzCut
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

getattr( process, 'pfElectronsFromVertex' + postfix ).d0Cut = pfD0Cut
getattr( process, 'pfElectronsFromVertex' + postfix ).dzCut = pfDzCut
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

from TopQuarkAnalysis.Configuration.patRefSel_refMuJets_cfi import *

# remove MC matching, object cleaning, objects etc.
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
                              , 'keep *_*_*_RECO'
                              ]
if runOnMC:
  process.out.outputCommands += [ 'keep GenEventInfoProduct_*_*_*'
                                , 'keep recoGenParticles_*_*_*'
                                , 'keep *_addPileupInfo_*_*'
                                ]




###
### Additional configuration
###

### Muons

intermediatePatMuons.src = cms.InputTag( 'selectedPatMuons' + postfix )
setattr( process, 'intermediatePatMuons' + postfix, intermediatePatMuons )

goodPatMuons.muonSource   = cms.InputTag( 'intermediatePatMuons' + postfix )
goodPatMuons.vertexSource = cms.InputTag( pfVertices )
setattr( process, 'goodPatMuons' + postfix, goodPatMuons )

step1.src = cms.InputTag( 'goodPatMuons' + postfix )
setattr( process, 'step1' + postfix, step1 )
step2.src = cms.InputTag( 'selectedPatMuons' + postfix )
setattr( process, 'step2' + postfix, step2 )


### Jets



### Electrons

step3 = cms.EDFilter(
  "PATCandViewCountFilter"
, src = cms.InputTag( 'selectedPatElectrons' )
, minNumber = cms.uint32( 1 )
, maxNumber = cms.uint32( 1 )
)

step3.src = cms.InputTag( 'selectedPatElectrons' + postfix )
setattr( process, 'step3' + postfix, step3 )

process.out.outputCommands.append( 'keep *_goodPatMuons*_*_*' )
process.out.outputCommands.append( 'keep *_selectedPatElectrons*_*_*' )

###
### Selection configuration
###

### Muons

getattr( process, 'patMuons' + postfix ).usePV      = muonsUsePV
getattr( process, 'patMuons' + postfix ).embedTrack = muonEmbedTrack

getattr( process, 'selectedPatMuons' + postfix ).cut = muonCut

getattr( process, 'intermediatePatMuons' + postfix ).cut = "pt > 20." #signalMuonCut

getattr( process, 'goodPatMuons' + postfix ).maxDZ = muonVertexMaxDZ

### Jets

### Electrons

getattr( process, 'patElectrons' + postfix ).electronIDSources = electronIDSources

getattr( process, 'selectedPatElectrons' + postfix ).cut = electronCut


###
### Scheduling
###

# MVA electron ID

process.load( "EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi" )
process.eidMVASequence = cms.Sequence(
  process.mvaTrigV0
+ process.mvaNonTrigV0
)

# The additional sequence

patAddOnSequence = cms.Sequence(
  getattr( process, 'intermediatePatMuons' + postfix )
* getattr( process, 'goodPatMuons'         + postfix )
)
setattr( process, 'patAddOnSequence' + postfix, patAddOnSequence )


# Heiners additional vertexing
import BTagDeltaR.My2ndVtxConfig.My2ndVtxConfig_cff as my_sv
process.extend(my_sv)
process.out.outputCommands += my_sv.my2ndVtxEventCont

# The paths

process.p = cms.Path()
process.p += process.step0a
process.p += process.goodOfflinePrimaryVertices
process.p += process.step0b
process.p += process.step0c
process.p += process.eidMVASequence
process.p += getattr( process, 'patPF2PATSequence' + postfix )
process.p += getattr( process, 'patAddOnSequence' + postfix )
process.p += getattr( process, 'step1' + postfix )
process.p += getattr( process, 'step2' + postfix )
process.p += getattr( process, 'step3' + postfix )
process.p += process.my2ndVtxSequence
process.out.SelectEvents.SelectEvents.append( 'p' )
