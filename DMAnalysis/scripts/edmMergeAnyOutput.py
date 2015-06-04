#!/usr/bin/env python

import os,sys
import getopt
import commands
CopyRights  = '\033[92m'
CopyRights += '####################################\n'
CopyRights += '#     edmMergeOutput.py Script     #\n'
CopyRights += '#           Ren-Jie Wang           #\n'
CopyRights += '#       renjie.wang@cern.ch        #\n'
CopyRights += '#             Mar 2014             #\n'
CopyRights += '####################################\n'
CopyRights += '\033[0m'


dir='-1'
#nsplit='-1'
print CopyRights

def help() :
   print '  -c --> gives a directory'
   #print '  -n --> # of split files'

#parse the options
try:
   # retrive command line options
   shortopts  = "c:n:?"
   opts, args = getopt.getopt( sys.argv[1:], shortopts )
except getopt.GetoptError:
   # print help information and exit:
   print "ERROR: unknown options in argument %s" % sys.argv[1:]
   help()
   sys.exit(1)

for o,a in opts:
   if o in("-?", "-h"):
      help()
      sys.exit(1)
   elif o in('-c'): dir = a
   #elif o in('-n'): nsplit = a


if(dir == '-1'):# or nsplit=='-1'):
   help()
   sys.exit(1)


count=0
status, allFiles = commands.getstatusoutput('ls ' + dir + '/res/Events_*.root | grep root ')
for line in allFiles.split():
        count = count+1
	#print line.split('\n')[0]
print 'Total Files: '+str(count)#+' nsplit: '+nsplit

#sys.exit(1)

outputList=''

MERGYFILE = open('run_merge.sh',"w")
MERGYFILE.writelines('#!bin/sh \n')

jcount=0
fcount=0
for line in allFiles.split():
	jobFile=line.split('\n')
	#print jobFile[0]
	jcount += 1
	outputList = outputList+'file:'+jobFile[0]
	#if jcount!=count: outputList += ','

	if (jcount%100==0) or (jcount==count) :
			exe='edmCopyPickMerge inputFiles='+outputList+' outputFile=/tmp/rewang/'+dir+'_'+str(fcount)+'.root maxSize=2000000'
			#print exe
			#os.system(exe)
			print 'writing to: script_edmMerge_'+str(fcount)+'.sh'
			SCRIPT = open('script_edmMerge_'+str(fcount)+'.sh',"w")
			SCRIPT.writelines(exe+';\n')
			SCRIPT.close()
			MERGYFILE.writelines('sh script_edmMerge_'+str(fcount)+'.sh >& script_edmMerge_'+str(fcount)+'.log & \n')
			#print '\n'
			fcount += 1
			outputList=''
			continue

	if jcount!=count: outputList += ','

MERGYFILE.close()
