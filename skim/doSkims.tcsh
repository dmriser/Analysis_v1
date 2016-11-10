#!/bin/tcsh 

set pathToFiles = "/volatile/clas/clas12/dmriser/analysis/e1f_analysis/debug_jobs2"

set runs = `ls $pathToFiles/clas_0*.root | cut -d'_' -f4 | cut -d'.' -f1 | sort | uniq`

foreach run ($runs)
	echo "doing $run"
	set targets = `ls $pathToFiles/clas_$run.*`
	./skim $targets 
end
