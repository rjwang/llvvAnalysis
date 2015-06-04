import FWCore.ParameterSet.Config as cms

#process = cms.Process("genAnalyzer")

#process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 5000
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
#                                        SkipEvent = cms.untracked.vstring('ProductNotFound')
#                                        )


genAnalyzer = cms.EDAnalyzer('GenAnalyzer',
	GenParticles = cms.untracked.InputTag("genParticles"),
	isPythia8 = cms.bool(True),
	GenJets = cms.untracked.InputTag("ak5GenJets")
#	GenJets = cms.untracked.InputTag("ak5GenJetsNoMuNoNu")

)

