#!/bin/bash
# this script changes the names of the telocate and temp absence caller in the "final_results" folder for each strain
# USE: renmae_final_results.sh



#find ./*/final_results/  -name "*_temp_absence_nonredundant.bed" -exec echo {}  \;
#find ./*/final_results/  -name "*_telocate_nonredundant.bed" -exec  echo {} \;
#mv intermediate.txt all_nonredundant.txt
file=${1}

#rename --prepend=no_gff_

while read -r sample
do
	cd ${sample}/final_results/
	mv ${sample}_temp_absence_nonredundant.bed org_${sample}_temp_absence_nonredundant.bed_org
	mv ${sample}_telocate_nonredundant.bed org_${sample}_telocate_nonredundant.bed_org
	cd ../..

done < "${file}"
