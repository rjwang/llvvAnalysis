import FWCore.ParameterSet.Config as cms

process = cms.Process("DMAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        )



process.mainAnalyzer = cms.EDAnalyzer('MainAnalyzer',
    dtag = cms.string('llvv'),
    isMC = cms.bool(True),
    verbose = cms.bool(False),
    isPythia8 = cms.bool(True),

    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    verticesTag = cms.InputTag("offlineSlimmedPrimaryVertices"),

    rhoAll = cms.InputTag("fixedGridRhoAll"),
    rhoFastjetAll = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoFastjetAllCalo = cms.InputTag("fixedGridRhoFastjetAllCalo"),
    rhoFastjetCentralCalo = cms.InputTag("fixedGridRhoFastjetCentralCalo"),
    rhoFastjetCentralChargedPileUp = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp"),
    rhoFastjetCentralNeutral = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),

    muonsTag = cms.InputTag("slimmedMuons"),

    electronsTag = cms.InputTag("slimmedElectrons"),
	## 25ns
    	electronVetoIdTag = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
	electronLooseIdTag = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
	electronMediumIdTag = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    	electronTightIdTag = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
	electronHEEPIdTag = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
	## 50ns
#    	electronVetoIdTag = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-CSA14-50ns-V1-standalone-veto"),
#	electronLooseIdTag = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-CSA14-50ns-V1-standalone-loose"),
#	electronMediumIdTag = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-CSA14-50ns-V1-standalone-medium"),
#    	electronTightIdTag = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-CSA14-50ns-V1-standalone-tight"),
#	electronHEEPIdTag = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51"),


    tausTag = cms.InputTag("slimmedTaus"),


    photonsTag = cms.InputTag("slimmedPhotons"),
    jetsTag = cms.InputTag("slimmedJets"),
    fatjetsTag = cms.InputTag("slimmedJetsAK8"),
    metsTag = cms.InputTag("slimmedMETs"),
    metsPuppiTag = cms.InputTag("slimmedMETsPuppi"),
    metFilterBitsTag = cms.InputTag("TriggerResults"),
    packedTag = cms.InputTag("packedGenParticles"),
    prunedTag = cms.InputTag("prunedGenParticles"),
    genJetsTag = cms.InputTag("slimmedGenJets"),


    ##trigger
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),

    	DoubleMuTrigs = cms.vstring("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
					"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
					"HLT_Mu27_TkMu8_v"),

    	DoubleEleTrigs = cms.vstring("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
					"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v",
					"HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v"),

    	SingleMuTrigs = cms.vstring("HLT_IsoMu20_v",
					"HLT_IsoTkMu20_v"),

    	SingleEleTrigs = cms.vstring("HLT_Ele27_eta2p1_WPLoose_Gsf_v"),

   	MuEGTrigs = cms.vstring("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v",
					"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"),

	#DoubleTauTrigs = cms.vstring("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v")
)
