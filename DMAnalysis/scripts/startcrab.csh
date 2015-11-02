#!/bin/csh
source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.csh
eval `scramv1 runtime -csh`
source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.csh

cmsenv;

voms-proxy-init -voms cms -valid 999:59
