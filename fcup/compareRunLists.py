#!/usr/bin/python

import glob, numpy

# run lists 
myList = open('/u/home/dmriser/mydoc/analysis/root_scripts/fcup/goodruns.txt','r')
nhList = open('/volatile/clas/clas12/dmriser/nathanFiles/AnalysisCode/programFiles/dataFiles.txt','r')

# all runs from mss
runs = glob.glob('/mss/clas/e1f/data/clas_0*')
possible = glob.glob('/volatile/clas/clas12/dmriser/analysis/ef1_analysis/root_files/clas_0*')
wpossible = glob.glob('/mss/clas/e1f/production/pass1/v1/data/run_0*')
unique_runs     = []
unique_possible = []
unique_wpossible = []

# get unique runs 
for run in runs:
    unique_runs.append(run.split('clas_0')[1].split('.A')[0])
    
unique_runs = numpy.unique(unique_runs)

for run in possible:
    unique_possible.append(run.split('clas_0')[1].split('.A')[0])

unique_possible = numpy.unique(unique_possible)

for run in wpossible:
    unique_wpossible.append(run.split('run_0')[1].split('.')[0])

unique_wpossible = numpy.unique(unique_wpossible)

# get my list 
mine = []

for run in myList: 
    mine.append(run.split('skim/')[1].split('.root')[0])

# get his list 
his = []

for run in nhList: 
    his.append(run.split('clas_0')[1].split('.A')[0])

his = numpy.unique(his)

# compare them 

print " RUN:   DR:   WG: "

bothY  = 0
bothN  = 0
heNmeY = 0
heYmeN = 0
heYmex = 0

for run in unique_runs: 
    iHasIt  = 'x'
    heHasIt = 'x'

    if (run in unique_possible):
        iHasIt = 'n'

    if (run in unique_wpossible):
        heHasIt = 'n'
 
    if (iHasIt is 'n') and (run in mine):
        iHasIt = 'y'
        
    if (heHasIt is 'n') and (run in his):
        heHasIt = 'y'

    if (heHasIt is 'y') and (iHasIt is 'y'):
        bothY+=1

    if (heHasIt is 'n') and (iHasIt is 'n'):
        bothN+=1

    if (heHasIt is 'n') and (iHasIt is 'y'):
        heNmeY+=1

    if (heHasIt is 'y') and (iHasIt is 'n'):
        heYmeN+=1

    if (heHasIt is 'y') and (iHasIt is 'x'):
        heYmex+=1

    print " %s   %s   %s " % (run, iHasIt, heHasIt)

print " ################################################# "
print " # runs from MSS         : ", len(unique_runs)
print " # of files I have before: ", len(unique_possible)
print " # of files I have before: ", len(unique_wpossible)
print " # of good files I have  : ", len(mine)
print " # of good files WG had  : ", len(his) 
print " # of both yes           : ", bothY
print " # of both no            : ", bothN
print " # WG: y  DR: n          : ", heYmeN
print " # WG: n  DR: y          : ", heNmeY
print " # WG: y  DR: x          : ", heYmex
print " ################################################# "
