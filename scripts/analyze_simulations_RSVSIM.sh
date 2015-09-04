#!/bin/bash
# this script calculates the averages of multiple simulation runs (such as average TPR and FDR)
# and genreates graphs/plots of the results with R
# NOTE: comment out "calculate_mean("retroseq")" and "calculate_mean("telocate")" in family_TFPN_average4.py and TFPN_average4.py as needed
# NOTE: family_average script doesn't always require using the renamed consensus seqs
# USE:analyze_simulations_RSVSIM.sh in appropriate folder

#recalculate_edit=RECALCULATE_EDIT.sh
#recalculate=RECALCULATE.sh
average=run_average_TFPNs_RSVSIM.sh
family_average=run_family_TFPNs_RSVSIM.sh
graph_TFPN=graph_TFPR_distances_RSVSIM.R
graph_family_TFPN=graph_TFPN_families.R
script_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts


consensus=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta
consensus_renamed=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/AB-PR/consensus_wTC8.fasta 
TE_lengths=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/LENGTHS/lengths_plusTC8.txt

#bash ${script_dir}/${recalculate} $length
#merge the distance calculations
mkdir averages
#calculate averages
bash ${script_dir}/${average} $consensus # consensus seq file doesn't matter here
#only sometimes run the below line
#bash ${script_dir}/${family_average} $consensus_renamed
#graph the TPR and FDR
cd averages
cat BEDCOMPARE_MEANS.txt | sed 's/_temp/_N2_temp/g' > tmp && mv tmp BEDCOMPARE_MEANS.txt
cat BEDCOMPARE_MEANS.txt | sed 's/_telocate/_N2_telocate/g' > tmp && mv tmp BEDCOMPARE_MEANS.txt

#Rscript ${script_dir}/${graph_TFPN}
#cd ..
#plot the TPR and FDR per family
#only sometimes run the below line
#cd families
#Rscript ${script_dir}/${graph_family_TFPN}
#cd ..

