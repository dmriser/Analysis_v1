#FaradayCup.git

Source code for Faraday Cup studies. 

       David Riser 
University of Connecticut 
	4-5-16

Under Modification: Feb 6, 2017

(i) readFilesWithEntryNumber creates: 
    (a) accumulation for files, and runs by sutracting element by element
    (b) list of bad events, for which dQ = 0 give by STARTN STOPN in text files 
    (c) some quality control plots

(ii) ratesFBF:
     (a) loads accumulation from file, does PID, takes ratios 
     (b) gives plots
     (c) gives good run list based on guassian cut on dN/dQ for all files

NOTE: 
      This code is used in conjunction with skim.C, 
      found at https://github.com/dmriser/Skim.git


     


