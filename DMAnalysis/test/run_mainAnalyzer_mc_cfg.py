import FWCore.ParameterSet.Config as cms

from llvvAnalysis.DMAnalysis.mainAnalyzer_cfi import *

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
process.GlobalTag.globaltag = 'PHYS14_25_V1'


#
# Cut-based Electron ID
#
# Set up everything that is needed to compute electron IDs and
# add the ValueMaps with ID decisions into the event data stream
#

# Load tools and function definitions
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
# overwrite a default parameter: for miniAOD, the collection name is a slimmed one
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')

from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)

# Define which IDs we want to produce
# Each of these two example IDs contains all four standard
# cut-based ID working points (only two WP of the PU20bx25 are actually used here).
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V1_miniAOD_cff']
#Add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

# Do not forget to add the egmGsfElectronIDSequence to the path,
# as in the example below!

#
# Cut-based Electron ID
#



process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    )
)

process.mainAnalyzer.isPythia8 = cms.bool(True)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("analysis.root")
                                  )






process.p = cms.Path(process.egmGsfElectronIDSequence * process.mainAnalyzer)
