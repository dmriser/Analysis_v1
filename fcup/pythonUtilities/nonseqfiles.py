#!/usr/bin/python

import glob, numpy


runs = glob.glob('/volatile/clas12/dmriser/analysis/ef1_analysis/root_files/clas_*')
unique_runs = []

for run in runs:
    unique_runs.append(run.split('clas_')[1].split('.A')[0])
    
unique_runs = numpy.unique(unique_runs)

jump_errors = []

for run in unique_runs:
    path = str('/volatile/clas12/dmriser/analysis/ef1_analysis/root_files/clas_*' + run + '*')
    subfiles = glob.glob(path)
    
    stubs = []
    
    for stub in subfiles:
        stubs.append(stub.split('.A')[1].split('.')[0])

    stubs = numpy.unique(stubs)
        
    for x in range(1,len(stubs)):
        if (int(stubs[x])-int(stubs[x-1]) != 1):
            print " jump error in %d from %d to %d " % (int(run),int(stubs[x-1]),int(stubs[x]))
            jump_errors.append(int(run))

jump_errors = numpy.unique(jump_errors)

print jump_errors
