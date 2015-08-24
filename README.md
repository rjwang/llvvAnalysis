llvvAnalysis
==============

General:

	git clone git@github.com:rjwang/llvvAnalysis.git


Prepare:

	setenv SCRAM_ARCH slc6_amd64_gcc491
	scramv1 project CMSSW CMSSW_7_4_2
	cd CMSSW_7_4_2/src
	cmsenv
	git cms-merge-topic ikrav:egm_id_747_v2

MiniAOD PHYS14:

	source /afs/cern.ch/cms/cmsset_default.csh
	voms-proxy-init
	xrdcp root://xrootd.unl.edu//store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/06C61714-7E6C-E411-9205-002590DB92A8.root /tmp/`whoami`/
	xrdcp root://xrootd.unl.edu//store/mc/Phys14DR/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/484D51C6-2673-E411-8AB0-001E67398412.root /tmp/`whoami`/
