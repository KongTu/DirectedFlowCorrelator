import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
#process.Timing = cms.Service("Timing")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options   = cms.untracked.PSet( wantSummary =
cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("Configuration.StandardSequences.Digi_cff")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.GlobalTag.globaltag = '75X_dataRun2_v12'

process.PAprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
    filter = cms.bool(True),   # otherwise it won't filter the events
)

#Reject beam scraping events standard pp configuration
process.NoScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAHighMultiplicity7/AOD/PromptReco-v1/000/285/480/00000/02BA31E5-08AF-E611-AAA3-FA163ED00180.root'
#'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/480/00000/32A34AA3-2CAF-E611-9C0D-FA163E8F093D.root'
#'root://cmsxrootd.fnal.gov//store/user/davidlw/PAHighMultiplicity7/RecoSkim2016_Pbp_V0Cascade_v1/170301_191341/0000/pPb_HM_1.root'
#'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/480/00000/4C0D189A-1BAF-E611-B5CA-FA163EAF1F45.root'
#'/store/hidata/HIRun2015/HIMinimumBias5/AOD/02May2016-v1/00000/002C1765-9B2E-E611-BA4B-F01FAFD5992D.root',
#'file:/afs/cern.ch/work/z/ztu/public/forZhenyu/00BF2A3E-7F93-E611-9644-F01FAFD9C9D0.root'
#'/store/hidata/HIRun2015/HIMinimumBias5/AOD/02May2016-v1/00000/00C836AF-0730-E611-A227-F01FAFD8F9BA.root'
#'file:/afs/cern.ch/work/z/ztu/CDDF/CMSSW_7_5_8_patch3/src/GenerateMassTree/DataTree/D0PbPb/pPb_HM_10.root'
#'file:PbPb_MB_101.root'
'file:PbPb_MB_1.root'
),

    secondaryFileNames = cms.untracked.vstring(
#'file:/afs/cern.ch/work/z/ztu/CDDF/CMSSW_7_5_8_patch3/src/GenerateMassTree/DataTree/D0PbPb/7A18BA88-7E91-E611-8582-0025901ACB5A.root'
#'file:6460B262-3DA7-E511-95FF-02163E014397.root'
'file:4EB647DE-5FA7-E511-9424-02163E013439.root',	
'file:00A973C7-5FA7-E511-9950-02163E014451.root'	
	)	
)

process.load("DirectedFlowCorrelator.DirectedFlowCorrelator.d0directedflowcorrelator_cfi")
#define the cuts
process.ana.useCentrality = True
process.ana.doEffCorrection = True
process.ana.doD0EffCorrection = False
process.ana.doHiReco = True
process.ana.trackName = "hiGeneralTracks"
process.ana.vertexName = "hiSelectedVertex"
process.ana.Nmin = 60
process.ana.Nmax = 160
process.ana.eff = 1

#loose cut
process.ana.D0VtxProbLow = cms.untracked.vdouble(0.0,0.0,0.04)
process.ana.D03DAngleHigh = cms.untracked.vdouble(0.16,0.12,0.16)
process.ana.D0DlosLow = cms.untracked.vdouble(3.0,3.5,3.5)

#tight cut
#process.ana.D0VtxProbLow = cms.untracked.vdouble(0.08,0.04,0.12)
#process.ana.D03DAngleHigh = cms.untracked.vdouble(0.08,0.06,0.08)
#process.ana.D0DlosLow = cms.untracked.vdouble(5.0,6.0,6.0)

process.TFileService = cms.Service("TFileService",fileName = cms.string("test_d0.root"))
process.p = cms.Path(  #process.hfCoincFilter3 *
                       #process.PAprimaryVertexFilter *
                       #process.NoScraping *
 		       process.ana)
