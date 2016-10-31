#!/usr/bin/python 

from subprocess import call

# job parameters 
n_files    = 11625
chunk_size = 100
n_jobs     = int(n_files/chunk_size) +1

print " > doing %s jobs \n" % (n_jobs)

# write jsubs and send 
for i in range(0,n_jobs): 
    
    jsub = open('temp.' + str(i) + '.jsub','w')
    jsub.write(
''' 
JOBNAME: e1f_rates
TRACK:   simulation
PROJECT: E1F

OTHER_FILES:/u/home/dmriser/mydoc/analysis/root_scripts/fcup/ratesFBF
/u/home/dmriser/mydoc/analysis/root_scripts/fcup/goodrunsROOT.txt \n
'''
)

    if ((n_files-i*chunk_size) > chunk_size): 
        jsub.write('COMMAND: ratesFBF ' + str(chunk_size) + ' ' + str(i*chunk_size) + '\n')
    else: 
        jsub.write('COMMAND: ratesFBF ' + str(n_files%chunk_size) + ' ' + str(i*chunk_size) + '\n')

    jsub.write('OUTPUT_DATA: ratesFBF.root \n')
    jsub.write('OUTPUT_TEMPLATE:/u/home/dmriser/mydoc/analysis/root_scripts/fcup/out/ratesFBF.' + str(i) + '.root \n')
    jsub.close()

    call(['jsub','temp.' + str(i) + '.jsub'])
    call(['rm','temp.' + str(i) + '.jsub'])
