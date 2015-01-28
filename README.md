llvvAnalysis
==============

General:
	
	git clone git@github.com:rjwang/llvvAnalysis.git 


Prepare:

	setenv SCRAM_ARCH slc6_amd64_gcc481
	scramv1 project CMSSW CMSSW_7_2_2
	cd CMSSW_7_2_2/src
	cmsenv

MiniAOD test:

	source /afs/cern.ch/cms/cmsset_default.csh
	voms-proxy-init
	xrdcp root://xrootd.unl.edu//store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5EC2A65C-7A6C-E411-94D2-002590DB92A8.root /tmp/`whoami`/
