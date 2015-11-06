#!/usr/bin/env python
import os,sys
import getopt
import commands

status, pwd = commands.getstatusoutput('echo ${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/html/')
fromfile = pwd+'index_template.html'
tofile = open(pwd+'index.html',"w")
status, date = commands.getstatusoutput('date')

with open(fromfile) as fp:

    for line in fp:
        if 'Created on DateToInsert by' in line:
                line = line.replace('DateToInsert',date)
	tofile.writelines(line)

tofile.close()

