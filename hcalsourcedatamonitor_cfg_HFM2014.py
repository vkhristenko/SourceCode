import FWCore.ParameterSet.Config as cms

process = cms.Process("HCALSourceDataMonitor")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['com10'] ## == GR_R_53_V16::All in 5_3_7
from CondCore.DBCommon.CondDBSetup_cfi import *

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#	Take care of input RUn numbers
import sys

if len(sys.argv)<4:
	print "### ERROR: No Run Files were provided"
	print "### Exiting..."
	sys.exit(1)

vRuns = []
numRuns = len(sys.argv)-3
outstr = ''
for i in range(numRuns):
	iRun = i+3
	str = 'root://eoscms//eos/cms/store/group/comm_hcal/LS1/USC_%s.root' % (
			sys.argv[iRun])
	outstr += '_' + sys.argv[iRun]
	vRuns.append(str)
	print "### To be processed: %s" % str

# input files
process.source = cms.Source("HcalTBSource",
    fileNames = cms.untracked.vstring(
		vRuns
    )
)

# TB data unpacker
process.tbunpack = cms.EDProducer("HcalTBObjectUnpacker",
    HcalSlowDataFED = cms.untracked.int32(-1),
    HcalSourcePositionFED = cms.untracked.int32(12),
    HcalTriggerFED = cms.untracked.int32(1),
    fedRawDataCollectionTag = cms.InputTag('rawDataCollector')
)

# histo unpacker
process.load("EventFilter.HcalRawToDigi.HcalHistogramRawToDigi_cfi")
process.hcalhistos.HcalFirstFED = cms.untracked.int32(700)
process.hcalhistos.FEDs = cms.untracked.vint32(718, 719, 720, 721, 722, 723)

process.hcalSourceDataMon = cms.EDAnalyzer('HCALSourceDataMonitor',
    RootFileName = cms.untracked.string('ntuples_HFM2014/HFM2014_' + sys.argv[2] +
		outstr + '.root'),
    PrintRawHistograms = cms.untracked.bool(False),
    SelectDigiBasedOnTubeName = cms.untracked.bool(True),
	AnalysisType = cms.untracked.int32(0),
	Verbosity = cms.untracked.int32(0)
)

process.p = cms.Path(process.tbunpack
                     *process.hcalhistos
                     *process.hcalSourceDataMon
                    )

