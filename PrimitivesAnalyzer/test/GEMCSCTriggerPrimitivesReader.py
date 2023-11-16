import FWCore.ParameterSet.Config as cms
#from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
from Configuration.Eras.Era_Run3_cff import Run3
# from Configuration.Eras.Era_Phase2_cff import Phase2

#process = cms.Process('analyzer',Phase2C9)
process = cms.Process('analyzer',Run3)
# process = cms.Process('analyzer', Phase2)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')

from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_v10', '')


process.MessageLogger.cerr.FwkReport.reportEvery = 5000

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register ('nEvents',
			-1, #Max number of events 
			VarParsing.multiplicity.singleton, 
			VarParsing.varType.int, 
			"Number of events")
options.parseArguments()

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(options.nEvents)
)
process.maxEvents.input = cms.untracked.int32(-1)


process.source = cms.Source("PoolSource", 
				fileNames = cms.untracked.vstring(options.inputFiles)
			# 	inputCommands = cms.untracked.vstring(
			# "keep *", 
			# # "drop TotemTimingDigiedmDetSetVector_totemTimingRawToDigi_TotemTiming_reRECO", 
			# # "drop TotemTimingRecHitedmDetSetVector_totemTimingRecHits__reRECO"
			# )
				)


process.source.fileNames.append("file:lcts2.root")

process.options = cms.untracked.PSet(
                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                        )
outfile = "out_GEMCSCTriggerPrimitivesReader.root"
process.TFileService = cms.Service("TFileService", fileName = cms.string(outfile))

process.GEMCSCTriggerPrimitivesReader = cms.EDAnalyzer('GEMCSCTriggerPrimitivesReader', 
	process.MuonServiceProxy,
	CSCLCTProducerData = cms.untracked.string("muonCSCDigis"),
	CSCLCTProducerEmul = cms.untracked.string("cscTriggerPrimitiveDigis"),
    debug = cms.bool(True),
)

process.p = cms.Path(process.GEMCSCTriggerPrimitivesReader)
