import FWCore.ParameterSet.Config as cms

process = cms.Process("DMAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        )


##___________________________HCAL_Noise_Filter________________________________||
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False)
process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
   reverseDecision = cms.bool(False)
)

process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
   reverseDecision = cms.bool(False)
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
	#
	# ID decisions (common to all formats)
	#
	eleMediumIdMapTrig = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
	eleTightIdMapTrig  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
	#
	# ValueMaps with MVA results
	#
	mvaValuesMapTrig     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
	mvaCategoriesMapTrig = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),

	## 50ns
#    	electronVetoIdTag = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-CSA14-50ns-V1-standalone-veto"),
#	electronLooseIdTag = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-CSA14-50ns-V1-standalone-loose"),
#	electronMediumIdTag = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-CSA14-50ns-V1-standalone-medium"),
#    	electronTightIdTag = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-CSA14-50ns-V1-standalone-tight"),
#	electronHEEPIdTag = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51"),


    tausTag = cms.InputTag("slimmedTaus"),

    photonsTag = cms.InputTag("slimmedPhotons"),

    jetsTag = cms.InputTag("slimmedJets"),
    jetsPuppiTag = cms.InputTag("slimmedJetsPuppi"),
    fatjetsTag = cms.InputTag("slimmedJetsAK8"),

    metsTag = cms.InputTag("slimmedMETs"),
    metsNoHFTag = cms.InputTag("slimmedMETsNoHF"),
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
					"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"
				   ),

    	DoubleEleTrigs = cms.vstring("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
					),

    	SingleMuTrigs = cms.vstring("HLT_IsoMu20_v",
					"HLT_IsoTkMu20_v"),

    	SingleEleTrigs = cms.vstring(#"HLT_Ele27_eta2p1_WPLoose_Gsf_v",
					"HLT_Ele23_WPLoose_Gsf_v",
					#"HLT_Ele22_eta2p1_WPLoose_Gsf_v"
					),

   	MuEGTrigs = cms.vstring("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v",
					"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"),

	#DoubleTauTrigs = cms.vstring("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v")
)

