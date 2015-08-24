import FWCore.ParameterSet.Config as cms

from llvvAnalysis.DMAnalysis.mainAnalyzer_cfi import *

process.mainAnalyzer.isMC = cms.bool(True)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


#load run conditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#------ Declare the correct global tag ------#

#global tag for CSA14 25ns 20 PU (asymptotic alignment and calibration) scenario
#process.GlobalTag.globaltag = 'PLS170_V7AN1::All'
#global tag for CSA14 50ns 40 PU (more pessimistic alignment and calibration) scenario
#process.GlobalTag.globaltag = 'PLS170_V6AN1::All'
#global tag for PHYS14 asymptotic 25ns scenario
#process.GlobalTag.globaltag = 'PHYS14_25_V3::All'
#process.GlobalTag.globaltag = 'MCRUN2_74_V9A::All'

process.GlobalTag.globaltag = 'MCRUN2_74_V9'



#
# Cut-based Electron ID
#

#
# Set up electron ID (VID framework)
#
## https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate
dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
## 25ns
#my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff']
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']
## 50ns
#my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_CSA14_50ns_V1_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff']



#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)




process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    )
)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("analysis.root")
                                  )






process.p = cms.Path(process.egmGsfElectronIDSequence * process.mainAnalyzer)
