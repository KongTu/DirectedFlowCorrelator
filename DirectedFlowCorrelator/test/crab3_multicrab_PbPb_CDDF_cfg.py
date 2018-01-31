from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

import FWCore.ParameterSet.Config as cms
#load the cfi file and rewrite cross section parameter each time:
process = cms.Process('Demo')
process.load("DirectedFlowCorrelator.DirectedFlowCorrelator.directedflowcorrelator_cfi")

outputName = "multicrab_CDDF_v11"

config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.allowUndistributedCMSSW = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'directedflowcorrelator_cfg.py'
config.Data.allowNonValidInputDataset = True
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 30
config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/ChargedParticle_DirectedFlow/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = outputName

config.Site.storageSite = 'T2_CH_CERN'

if __name__ == '__main__':
   from CRABAPI.RawCommand import crabCommand
   from CRABClient.ClientExceptions import ClientException
   from httplib import HTTPException

   config.General.workArea = outputName

   def submit(config):
      try:
           crabCommand('submit', config = config)
      except HTTPException as hte:
           print "Failed submitting task: %s" % (hte.headers)
      except ClientException as cle:
          print "Failed submitting task: %s" % (cle)
   
   for num in range(0,1):

  	RequestName = outputName + "_MB5_Golden_" + str(num)
   	DataSetName = '/HIMinimumBias5/HIRun2015-02May2016-v1/AOD'
   	config.Data.lumiMask = 'Cert_262548-263757_PromptReco_HICollisions15_JSON_v2.txt'
   	config.General.requestName = RequestName
   	config.Data.inputDataset = DataSetName
   	submit(config)
	
	RequestName = outputName + "_MB6_Golden_" + str(num)
   	DataSetName = '/HIMinimumBias6/HIRun2015-02May2016-v1/AOD'
   	config.Data.lumiMask = 'Cert_262548-263757_PromptReco_HICollisions15_JSON_v2.txt'
   	config.General.requestName = RequestName
   	config.Data.inputDataset = DataSetName
   	submit(config)
        
	
	RequestName = outputName + "_MB7_Golden_" + str(num)
   	DataSetName = '/HIMinimumBias7/HIRun2015-02May2016-v1/AOD'
   	config.Data.lumiMask = 'Cert_262548-263757_PromptReco_HICollisions15_JSON_v2.txt'
   	config.General.requestName = RequestName
   	config.Data.inputDataset = DataSetName
   	submit(config)


	RequestName = outputName + "_MB5_TrackerOnly_" + str(num)
        DataSetName = '/HIMinimumBias5/HIRun2015-02May2016-v1/AOD'
        config.Data.lumiMask = 'Cert_263685-263757_PromptReco_HICollisions15_TrackerOnly_JSON.txt'
        config.General.requestName = RequestName
        config.Data.inputDataset = DataSetName
        submit(config)

        RequestName = outputName + "_MB6_TrackerOnly_" + str(num)
        DataSetName = '/HIMinimumBias6/HIRun2015-02May2016-v1/AOD'
        config.Data.lumiMask = 'Cert_263685-263757_PromptReco_HICollisions15_TrackerOnly_JSON.txt'
        config.General.requestName = RequestName
        config.Data.inputDataset = DataSetName
        submit(config)


        RequestName = outputName + "_MB7_TrackerOnly_" + str(num)
        DataSetName = '/HIMinimumBias7/HIRun2015-02May2016-v1/AOD'
        config.Data.lumiMask = 'Cert_263685-263757_PromptReco_HICollisions15_TrackerOnly_JSON.txt'
        config.General.requestName = RequestName
        config.Data.inputDataset = DataSetName
        submit(config)	
# python crab3_ppTrackingAnalyzer.py to execute 
# ./multicrab -c status -w crab_projects/ to check status 
