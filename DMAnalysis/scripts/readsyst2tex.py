#!/usr/bin/env python

import os,sys
import getopt
import commands
import math


DAT='-1'

def help() :
   print '\033[92m --------------------------------- \033[0m '
   print '\033[92m python readsyst2tex.py -f <.dat> \033[0m '
   print '\033[92m --------------------------------- \033[0m '

#parse the options
try:
   # retrive command line options
   shortopts  = "f:?"
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
   elif o in('-f'): DAT = a

if DAT == '-1':
   help()
   sys.exit(1)

print 'Input DAT: '+DAT


#### start to read dat file
f = open(DAT)
lines = f.readlines()
f.close()


obs = 0

### get procName, Yields and procTag
for m in lines:
	lline = m.split()
	if lline[0]=='process' and lline[1]!='0': procName = lline
	if lline[0]=='process' and lline[1]=='0': procTag = lline
	if lline[0]=='rate': yields = lline
	if lline[0]=='observation':
		print lline
		## sum over all obs
		for x in range(1, len(lline)):
			obs += float(lline[x])
	#if lline.count('lnN')==1:
	#	print lline

#print '\n'
print yields
print procName
print procTag


print '--------------------------------------'
for proc in range(len(procName)):

	if proc==0: continue
	#print procName[proc]
	#print yields[proc]
	#print '---------------------'

	rel_syst=0
	rel_stat=0

	for m in lines:
		lline = m.split()
		if lline.count('lnN')==1:
			#print lline
			if lline[proc+1] != '-':
			  if '_stat_' in lline[0]:
				#print 'stat:'+lline[proc+1]
				rel_stat = float(lline[proc+1])-1
			  else:
				#print lline[proc+1]
				rel_syst += (float(lline[proc+1])-1)*(float(lline[proc+1])-1)

			#sys.exit(1)

	rel_syst = math.sqrt(rel_syst)

	#print '---------------------'
	#print rel_stat
	#print rel_syst
	val_stat = rel_stat*float(yields[proc])
	val_syst = rel_syst*float(yields[proc])
	#print '---------------------'
	#print procName[proc]+': '+yields[proc]+' $\pm$ '+str(val_stat)+' (stat.)'+' $\pm$ ' + str(val_syst) + ' (syst.)'
	#sys.exit(1)


	pm='$\pm$'

	if( ('%0.2f'%float(yields[proc]) != '0.00') and ('%0.2f'%val_stat != '0.00') and ('%0.2f'%val_syst != '0.00')):
		print '%5s '%procName[proc] + '%0.2f'%float(yields[proc]) + ' %4s '%pm + '%0.2f'%val_stat + ' %4s '%pm + '%0.2f'%val_syst
	elif( ('%0.3f'%float(yields[proc]) != '0.000') and ('%0.3f'%val_stat != '0.000') and ('%0.3f'%val_syst != '0.000')):
			print '%5s '%procName[proc] + '%0.3f'%float(yields[proc]) + ' %4s '%pm + '%0.3f'%val_stat + ' %4s '%pm + '%0.3f'%val_syst
	else:
		print '%5s '%procName[proc] + '%0.4f'%float(yields[proc]) + ' %4s '%pm + '%0.4f'%val_stat + ' %4s '%pm + '%0.4f'%val_syst




##########################
## calculate stat., syst. for total bkg
##########################
totalbkg_syst = []
totalbkg_stat = []

for m in lines:
	lline = m.split()
	if lline.count('lnN')==1:
		syst_val = 0
		stat_val = 0
		#print lline
		for proc in range(len(procName)):
		  if proc>0 and float(procTag[proc])>0:
			#print procTag[proc]
			if lline[proc+1] != '-':
			  #print 'yields: ' + yields[proc] + ' --> ' + lline[proc+1]
			  if '_stat_' in lline[0]:
				stat_val += (float(lline[proc+1])-1.)*float(yields[proc])
			  else:
				syst_val += (float(lline[proc+1])-1.)*float(yields[proc])

		#print syst_val
		totalbkg_syst.append(syst_val)
		totalbkg_stat.append(stat_val)
		#print '\n\n'
#print totalbkg_syst
#print totalbkg_stat
#################
tot_bkg=0
tot_bkg_stat=0
tot_bkg_syst=0

for isyst in totalbkg_syst:
	#print isyst
	tot_bkg_syst += isyst*isyst

for istat in totalbkg_stat:
        #print istat
        tot_bkg_stat += istat*istat

for proc in range(len(procName)):
  	if proc>0 and float(procTag[proc])>0:
		tot_bkg += float(yields[proc])


tot_bkg_syst = math.sqrt(tot_bkg_syst)
tot_bkg_stat = math.sqrt(tot_bkg_stat)

#print tot_bkg_syst
#print tot_bkg_stat

print '--------------------------------------'
print '%5s '%'Totbkg' + '%0.2f'%tot_bkg + ' %4s '%pm + '%0.2f'%tot_bkg_stat + ' %4s '%pm + '%0.2f'%tot_bkg_syst
print '%5s '%'Data  ' + '%0.2f'%obs
print '--------------------------------------'


