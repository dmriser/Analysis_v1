#!/bin/tcsh 

set runs = `ls /volatile/clas/clas12/dmriser/analysis/ef1_analysis/root_files/clas_0*.root | cut -d'_' -f4 | cut -d'.' -f1 | sort | uniq`

foreach run ($runs)
	echo "doing $run"
	set targets = `ls /volatile/clas/clas12/dmriser/analysis/ef1_analysis/root_files/clas_$run.*`
	./skim $targets 
end
