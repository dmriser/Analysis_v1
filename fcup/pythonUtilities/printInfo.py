import numpy as np
import matplotlib.pyplot as plt

infoFile = '/u/home/dmriser/mydoc/analysis/root_scripts/Analysis_v2/lists/runsNew.info'

runs = []
charges = []

with open(infoFile,'rU') as file:
    for line in file:
        stubs = line.split(' ')
        charge = stubs[1]
        run = stubs[0]
        scale = stubs[3]
        
        runs.append(run)
        charges.append(float(charge)/float(scale))
 

charges_np = np.array(charges)
print(' Total Charge = ', charges_np.sum() )

