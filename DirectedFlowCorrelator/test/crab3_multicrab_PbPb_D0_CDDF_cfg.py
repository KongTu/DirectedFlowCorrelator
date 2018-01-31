from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

import FWCore.ParameterSet.Config as cms
#load the cfi file and rewrite cross section parameter each time:
process = cms.Process('Demo')
process.load("DirectedFlowCorrelator.DirectedFlowCorrelator.d0directedflowcorrelator_cfi")

outputName = "multicrab_D0_CDDF_v11_loose_HI_PeriDatasets"

config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.allowUndistributedCMSSW = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'd0directedflowcorrelator_cfg.py'
config.Data.allowNonValidInputDataset = True
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 7
config.Data.ignoreLocality = True
config.Data.useParent = True
config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/D0_DirectedFlow/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = outputName

config.Site.whitelist = ['T2_US_Vanderbilt']
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
   

   for num in range(0,4):
	"""
        RequestName = outputName + "Golden_MB_" + str(num+1)
        DataSetName = '/HIMinimumBias' + str(num+1) + '/davidlw-RecoSkim2015_D0_hireco_Golden_MB_tight_v2-51b0d5d038f3f98617a2426a9749f15c/USER'
        config.General.requestName = RequestName
        config.Data.inputDataset = DataSetName
        submit(config)
        """
     
   for num in range(0,11):
	
	"""
        RequestName = outputName + "TrackerOnly_MB_" + str(num+1)
        DataSetName = '/HIMinimumBias' + str(num+1) + '/davidlw-RecoSkim2015_D0_hireco_TrackerOnly_MB_v2-51b0d5d038f3f98617a2426a9749f15c/USER'
        config.General.requestName = RequestName
        config.Data.inputDataset = DataSetName
        submit(config)
	"""

   for num in range(0,1):

	"""
	RequestName = outputName + "TrackerOnly_DoubleMu_0"
        DataSetName = '/HIOniaL1DoubleMu0/davidlw-RecoSkim2015_D0_hireco_TrackerOnly_MB_v2-51b0d5d038f3f98617a2426a9749f15c/USER'
        config.General.requestName = RequestName
        config.Data.inputDataset = DataSetName
        submit(config)

	RequestName = outputName + "TrackerOnly_DoubleMu_0B"
        DataSetName = '/HIOniaL1DoubleMu0B/davidlw-RecoSkim2015_D0_hireco_TrackerOnly_MB_v2-51b0d5d038f3f98617a2426a9749f15c/USER'
        config.General.requestName = RequestName
        config.Data.inputDataset = DataSetName
        submit(config)

	RequestName = outputName + "TrackerOnly_DoubleMu_0C"
        DataSetName = '/HIOniaL1DoubleMu0C/davidlw-RecoSkim2015_D0_hireco_TrackerOnly_MB_v2-51b0d5d038f3f98617a2426a9749f15c/USER'
        config.General.requestName = RequestName
        config.Data.inputDataset = DataSetName
        submit(config)

	RequestName = outputName + "TrackerOnly_DoubleMu_0D"
        DataSetName = '/HIOniaL1DoubleMu0D/davidlw-RecoSkim2015_D0_hireco_TrackerOnly_MB_v2-51b0d5d038f3f98617a2426a9749f15c/USER'
        config.General.requestName = RequestName
        config.Data.inputDataset = DataSetName
        submit(config)

	#print 'double check the q2 cut is from %r to %r' % (q2Range[num], q2Range[num+1])
	#print 'double check the Ntrk cut is from %r to %r' % (ntrkRange[num], ntrkRange[num+1])  	
	"""

	RequestName = outputName + "Golden_MB5_" + str(num)
        DataSetName = '/HIMinimumBias5/davidlw-RecoSkim2015_D0_hireco_Golden_v1-0994449c32ecc197ef90c1037e0ac608/USER'
	config.General.requestName = RequestName
        config.Data.inputDataset = DataSetName
        submit(config)
	
 
        RequestName = outputName + "Golden_MB6_" + str(num)
        DataSetName = '/HIMinimumBias6/davidlw-RecoSkim2015_D0_hireco_Golden_v1-34d7cdadb6520433feaa67ecd8716835/USER'
	config.General.requestName = RequestName
        config.Data.inputDataset = DataSetName
        submit(config)
        
        RequestName = outputName + "Golden_MB7_" + str(num)
        DataSetName = '/HIMinimumBias7/davidlw-RecoSkim2015_D0_hireco_Golden_v1-34d7cdadb6520433feaa67ecd8716835/USER'
	config.General.requestName = RequestName
        config.Data.inputDataset = DataSetName
        submit(config)
     	
	RequestName = outputName + "Golden_MB3_" + str(num)
        DataSetName = '/HIMinimumBias3/davidlw-RecoSkim2015_D0_hireco_Golden_v1-34d7cdadb6520433feaa67ecd8716835/USER'
	config.General.requestName = RequestName
        config.Data.inputDataset = DataSetName
        submit(config)
        
        RequestName = outputName + "Golden_MB1_" + str(num)
        DataSetName = '/HIMinimumBias1/davidlw-RecoSkim2015_D0_hireco_Golden_v1-34d7cdadb6520433feaa67ecd8716835/USER'
        config.General.requestName = RequestName
        config.Data.inputDataset = DataSetName
        submit(config) 
        

	RequestName = outputName + "TrackerOnly_MB5_" + str(num)
        DataSetName = '/HIMinimumBias5/davidlw-RecoSkim2015_D0_hireco_TrackerOnly_v1-34d7cdadb6520433feaa67ecd8716835/USER'
	config.General.requestName = RequestName
        config.Data.inputDataset = DataSetName
        submit(config)

        
        RequestName = outputName + "TrackerOnly_MB6_" + str(num)
        DataSetName = '/HIMinimumBias6/davidlw-RecoSkim2015_D0_hireco_TrackerOnly_v1-34d7cdadb6520433feaa67ecd8716835/USER'
	config.General.requestName = RequestName
        config.Data.inputDataset = DataSetName
        submit(config)
        
        RequestName = outputName + "TrackerOnly_MB7_" + str(num)
        DataSetName = '/HIMinimumBias7/davidlw-RecoSkim2015_D0_hireco_TrackerOnly_v1-34d7cdadb6520433feaa67ecd8716835/USER'
	config.General.requestName = RequestName
        config.Data.inputDataset = DataSetName
        submit(config) 
	
	 
#pp reco
	"""
	#process.ana.Nmin = ntrkRange[num]
	#process.ana.Nmax = ntrkRange[num+1]
  	RequestName = outputName + "_MB5_" + str(num)
   	DataSetName = '/HIMinimumBias5/davidlw-RecoSkim2015_D0_pprereco_Golden_v1-45eb8d5205ff290372e262360cc383f2/USER'
   	config.General.requestName = RequestName
   	config.Data.inputDataset = DataSetName
   	submit(config)
        	
        #process.ana.Nmin = ntrkRange[num]
        #process.ana.Nmax = ntrkRange[num+1]
	RequestName = outputName + "_MB6_" + str(num)
   	DataSetName = '/HIMinimumBias6/davidlw-RecoSkim2015_D0_pprereco_Golden_v1-45eb8d5205ff290372e262360cc383f2/USER'
   	config.General.requestName = RequestName
   	config.Data.inputDataset = DataSetName
   	submit(config)
        
	#process.ana.Nmin = ntrkRange[num]
        #process.ana.Nmax = ntrkRange[num+1]
	RequestName = outputName + "_MB7_" + str(num)
   	DataSetName = '/HIMinimumBias7/davidlw-RecoSkim2015_D0_pprereco_Golden_v1-45eb8d5205ff290372e262360cc383f2/USER'
   	config.General.requestName = RequestName
   	config.Data.inputDataset = DataSetName
   	submit(config)
	"""
# python crab3_ppTrackingAnalyzer.py to execute 
# ./multicrab -c status -w crab_projects/ to check status 
