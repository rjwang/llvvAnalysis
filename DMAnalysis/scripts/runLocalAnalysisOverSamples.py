#!/usr/bin/env python
import os,sys
import json
import ROOT
import getopt
import commands

#print usage
def usage() :
    print ' '
    print 'runLocalAnalysisOverSamples.py [options]'
    print '  -s : submit to queue'
    print '  -e : executable name'
    print '  -j : json file containing the samples'
    print '  -d : input dir with the event summaries'
    print '  -o : output directory'
    print '  -c : templated configuration file to steer the job'
    print '  -l : luminosity (pb)'
    print '  -p : extra parameters configure'
    print '  -t : tag to match sample'
    print '  -g : log file from submitting to queue' #RJ
    print ' '
    exit(-1)

"""
Gets the value of a given item
(if not available a default value is returned)
"""
def getByLabel(desc,key,defaultVal=None) :
    try :
        return desc[key]
    except KeyError:
        return defaultVal


#parse the options
try:
     # retrive command line options
     shortopts  = "s:e:j:d:o:c:l:p:t:h:g:?" #RJ
     opts, args = getopt.getopt( sys.argv[1:], shortopts )
except getopt.GetoptError:
     # print help information and exit:
     print "ERROR: unknown options in argument %s" % sys.argv[1:]
     usage()
     sys.exit(1)

subtoBatch=False
requirementtoBatch='type==SLC6_64&&pool>30000'
samplesDB=''
theExecutable=''
inputdir=''
outdir=''
lumi=1
cfg_file=''
split=1
segment=0
params=''
onlytag='all'
queuelog=''

count=0
who = commands.getstatusoutput('whoami')[1]
SCRIPT = open('/tmp/'+who+'/SCRIPT_Submit2batch.sh',"w")
SCRIPT_L = open('/tmp/'+who+'/SCRIPT_Local.sh',"w")
SCRIPT_L.writelines('#!bin/sh \n\n')
SCRIPT_L.writelines('cd $CMSSW_BASE/src/llvvAnalysis/DMAnalysis/; \n\n')

for o,a in opts:
    if o in("-?", "-h"):
        usage()
        sys.exit(0)
    elif o in('-s'):
        subtoBatch=True
        queue=a
        if(queue=="True") : queue="2nd"
    elif o in('-j'): samplesDB = a
    elif o in('-e'): theExecutable = a
    elif o in('-d'): inputdir = a
    elif o in('-o'): outdir = a
    elif o in('-l'): lumi=float(a)
    elif o in('-c'): cfg_file = a
    elif o in('-p'): params = a
    elif o in('-t'): onlytag = a
    elif o in('-g'): queuelog = a #RJ

#open the file which describes the sample
jsonFile = open(samplesDB,'r')
procList=json.load(jsonFile,encoding='utf-8').items()

#run over sample
for proc in procList :

    #run over processes
    for desc in proc[1] :

        #run over items in process
        isdata=getByLabel(desc,'isdata',False)
        mctruthmode=getByLabel(desc,'mctruthmode',0)
	tag = getByLabel(desc,'tag','') #RJ
	print tag

	mytag = tag
	mytag = mytag.replace("#","")
	mytag = mytag.replace(" ","")
	mytag = mytag.replace("(","")
	mytag = mytag.replace(")","")
	mytag = mytag.replace("{","")
	mytag = mytag.replace("}","")
	mytag = mytag.replace("+","")
        ## split jobs by tag name
        SCRIPT_Temp = open('/tmp/'+who+'/SCRIPT_Local_'+mytag+'.sh',"w")
        SCRIPT_Temp.writelines('#!bin/sh \n\n')
        SCRIPT_Temp.writelines('cd $CMSSW_BASE/src/llvvAnalysis/DMAnalysis/; \n\n')

        data = desc['data']
        for d in data :
            origdtag = getByLabel(d,'dtag','')
            dtag = origdtag
            xsec = getByLabel(d,'xsec',-1)
            br = getByLabel(d,'br',[])
            suffix = str(getByLabel(d,'suffix' ,""))
            if(onlytag!='all') :
                if(dtag.find(onlytag)<0) : continue
            if(mctruthmode!=0) : dtag+='_filt'+str(mctruthmode)

            if(xsec>0 and not isdata) :
                for ibr in br :  xsec = xsec*ibr
            split=getByLabel(d,'split',1)


            mydtag = dtag
            mydtag = mydtag.replace("#","")
            mydtag = mydtag.replace(" ","")
            mydtag = mydtag.replace("(","")
            mydtag = mydtag.replace(")","")
            mydtag = mydtag.replace("{","")
            mydtag = mydtag.replace("}","")
            mydtag = mydtag.replace("+","")
            ## split jobs by dtag name
            SCRIPT_DTag = open('/tmp/'+who+'/SCRIPT_Local_'+mydtag+'.sh',"w")
            SCRIPT_DTag.writelines('#!bin/sh \n\n')
            SCRIPT_DTag.writelines('cd $CMSSW_BASE/src/llvvAnalysis/DMAnalysis/; \n\n')


	    for segment in range(0,split) :
                if(split==1):
	            	eventsFile=inputdir + '/' + origdtag + '.root'
                else:
                        eventsFile=inputdir + '/' + origdtag + '_' + str(segment) + '.root'

                if(eventsFile.find('/store/')==0)  : eventsFile = commands.getstatusoutput('cmsPfn ' + eventsFile)[1]
		if(eventsFile.find('?')>=0)  : eventsFile = eventsFile[:eventsFile.find('?')]

            	sedcmd = 'sed \"s%@input%' + eventsFile +'%;s%@outdir%' + outdir +'%;s%@isMC%' + str(not isdata) + '%;s%@mctruthmode%'+str(mctruthmode)+'%;s%@xsec%'+str(xsec)+'%;'
                sedcmd += 's%@cprime%'+str(getByLabel(d,'cprime',-1))+'%;'
                sedcmd += 's%@brnew%' +str(getByLabel(d,'brnew' ,-1))+'%;'
                sedcmd += 's%@suffix%' +suffix+'%;'
		sedcmd += 's%@tag%' +str(getByLabel(desc,'tag',-1))+'%;'#RJ
		sedcmd += 's%@wimpweights%' +str(getByLabel(d,'wimpweights',-1))+'%;'#RJ
            	if(params.find('@useMVA')<0) :          params = '@useMVA=False ' + params
                if(params.find('@weightsFile')<0) :     params = '@weightsFile= ' + params
                if(params.find('@evStart')<0) :         params = '@evStart=0 ' + params
                if(params.find('@evEnd')<0) :           params = '@evEnd=-1 ' + params
            	if(params.find('@saveSummaryTree')<0) : params = '@saveSummaryTree=False ' + params
            	if(params.find('@runSystematics')<0) :  params = '@runSystematics=False ' + params
            	if(len(params)>0) :
                    extracfgs = params.split(' ')
                    for icfg in extracfgs :
                        varopt=icfg.split('=')
                        if(len(varopt)<2) : continue
                        sedcmd += 's%' + varopt[0] + '%' + varopt[1] + '%;'
            	sedcmd += '\"'


		if(split==1):
	            	cfgfile=outdir +'/'+ dtag + suffix + '_cfg.py'
		else:
                        cfgfile=outdir +'/'+ dtag + suffix + '_' + str(segment) + '_cfg.py'
            	os.system('cat ' + cfg_file + ' | ' + sedcmd + ' > ' + cfgfile)
            	if(not subtoBatch) :
                	os.system(theExecutable + ' ' + cfgfile)
            	else :
			os.system('mkdir -p ' + queuelog)
			print('\033[33m submit2batch.sh -q'+queue+' -G'+queuelog+'/'+dtag+str(segment)+'.log'+' -R"' + requirementtoBatch + '" -J' + dtag + str(segment) + ' ${CMSSW_BASE}/bin/${SCRAM_ARCH}/wrapLocalAnalysisRun.sh ' + theExecutable + ' ' + cfgfile + '\033[0m')
			SCRIPT.writelines('submit2batch.sh -q'+queue+' -G'+queuelog+'/'+dtag+str(segment)+'.log'+' -R"' + requirementtoBatch + '" -J' + dtag + str(segment) + ' ${CMSSW_BASE}/bin/${SCRAM_ARCH}/wrapLocalAnalysisRun.sh ' + theExecutable + ' ' + cfgfile + '\n\n')
			SCRIPT_L.writelines(theExecutable + ' ' + cfgfile + ' >& '+queuelog+'/'+dtag+str(segment)+'.log'+' & \n\n')
			count = count + 1
			if count % 30 == 0: SCRIPT_L.writelines('sleep 25\n\n')

			SCRIPT_Temp.writelines(theExecutable + ' ' + cfgfile + ' >& '+queuelog+'/'+dtag+str(segment)+'.log'+' & \n\n')
			SCRIPT_DTag.writelines(theExecutable + ' ' + cfgfile + ' >& '+queuelog+'/'+dtag+str(segment)+'.log'+' & \n\n')
			#sys.exit(0)
			#os.system('submit2batch.sh -q'+queue+' -G'+queuelog+'/'+dtag+str(segment)+'.log'+' -R"' + requirementtoBatch + '" -J' + dtag + str(segment) + ' ${CMSSW_BASE}/bin/${SCRAM_ARCH}/wrapLocalAnalysisRun.sh ' + theExecutable + ' ' + cfgfile)

            SCRIPT_DTag.writelines('cd -;')
            SCRIPT_DTag.close()
            os.system('chmod u+x,g-r,o-r '+'/tmp/'+who+'/SCRIPT_Local_'+mydtag+'.sh ')
	    os.system('mkdir -p '+queuelog+'/split/')
            os.system('mv /tmp/'+who+'/SCRIPT_Local_'+mydtag+'.sh '+queuelog+'/split')

        SCRIPT_Temp.writelines('cd -;')
        SCRIPT_Temp.close()
	os.system('chmod u+x,g-r,o-r '+'/tmp/'+who+'/SCRIPT_Local_'+mytag+'.sh ')
	os.system('mkdir -p '+queuelog+'/combine/')
        os.system('mv /tmp/'+who+'/SCRIPT_Local_'+mytag+'.sh '+queuelog+'/combine/')

SCRIPT.close()
SCRIPT_L.writelines('cd -;')
SCRIPT_L.close()
os.system('chmod u+x,g-r,o-r '+'/tmp/'+who+'/SCRIPT_Submit2batch.sh '+'/tmp/'+who+'/SCRIPT_Local.sh ')
os.system('mkdir -p '+queuelog+'/all/')
os.system('mv /tmp/'+who+'/SCRIPT_Submit2batch.sh '+queuelog+'/all/')
os.system('mv /tmp/'+who+'/SCRIPT_Local.sh '+queuelog+'/all/')
os.system('cp $CMSSW_BASE/src/llvvAnalysis/DMAnalysis/scripts/splitlocaljobs.py '+queuelog+'/all/')
os.system('cp $CMSSW_BASE/src/llvvAnalysis/DMAnalysis/scripts/checkLocaljobs.py '+queuelog+'/all/')


##get good script for local run
#os.system('cp '+queuelog+'/combine/SCRIPT_Local_ZZ*.sh '+queuelog+'/')
#os.system('cp '+queuelog+'/combine/SCRIPT_Local_WZ*.sh '+queuelog+'/')
#os.system('cp '+queuelog+'/combine/SCRIPT_Local_Wjets*.sh '+queuelog+'/')
#os.system('cp '+queuelog+'/combine/SCRIPT_Local_WW*.sh '+queuelog+'/')
#os.system('cp '+queuelog+'/combine/SCRIPT_Local_Singletop.sh '+queuelog+'/')
#os.system('cp '+queuelog+'/combine/SCRIPT_Local_tbart.sh '+queuelog+'/')
#os.system('cp '+queuelog+'/combine/SCRIPT_Local_HZZd*.sh '+queuelog+'/')
#os.system('cp '+queuelog+'/split/SCRIPT_Local_Data8TeV_*.sh '+queuelog+'/')
#os.system('cp '+queuelog+'/split/SCRIPT_Local_MC8TeV_DY*.sh '+queuelog+'/')


#run plotter over results
#if(not subtoBatch) :
#    os.system('mkdir -p ' + outdir + '/plots')
#    #os.system('runPlotterOverSamples.py -j ' + samplesDB + ' -l ' + str(lumi) + ' -i ' + outdir + ' -o ' + outdir + '/plots -d localAnalysis')
#    os.system('runPlotter --json ' + samplesDB + ' --iLumi ' + str(lumi) + ' --inDir ' + outdir + ' --outDir ' + outdir + '/plots')
