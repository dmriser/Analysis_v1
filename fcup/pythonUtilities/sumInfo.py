#!/usr/bin/python 

import glob 

# ---------- params ------- 
testNumber = 3 
fPath = '/volatile/clas12/dmriser/analysis/e1f_analysis/fca/'

# -------------------------

runFile = glob.glob(fPath + str(testNumber) + '*')
allFiles = glob.glob(fPath + str(testNumber) + '*')

rFile = open(runFile[0])
total = int(rFile.read().split(' ')[0])/9624000.0

print " dQ total for run: ", total ," uC"
tot = 0

for f in allFiles:
    cFile = open(f)
    fca = cFile.read().split()[0] 
    tot += int(fca)
    print " dQ for %s : %s " % (f, fca)  

print " sum over files : ", int(tot)/9624000.0

