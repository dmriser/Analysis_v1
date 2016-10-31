#!/usr/bin/python 

import glob 

# ---------- params ------- 
testNumber = 38222 
fPath = '/volatile/clas/clas12/dmriser/analysis/e1f_analysis/fca/'

# -------------------------

runFile = glob.glob(fPath + str(testNumber) + '*')
allFiles = glob.glob(fPath + 'clas_0' + str(testNumber) + '*')

rFile = open(runFile[0])
total = int(rFile.read().split(' ')[0])/9624.00 *1e-3

print " dQ total for run: ", total ," uC"
tot = 0

for f in allFiles:
    cFile = open(f)
    fca = cFile.read() 
    tot += int(fca)
    print " dQ for %s : %s " % (f, fca)  

print " sum over files : ", int(tot)/9624.00 * 1e-3
